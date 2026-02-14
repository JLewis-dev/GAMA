# CORE LOGIC ==================================================================

# Data richness score

.extract_assembly_level <- function(x) {
  if (is.list(x) && 'assembly_level' %in% names(x)) {
    lvl <- .flatten_to_char(x$assembly_level)
    if (!is.na(lvl) && lvl != '') return(lvl)
  }
  if (is.list(x) && 'assemblystatus' %in% names(x)) {
    lvl <- .flatten_to_char(x$assemblystatus)
    if (!is.na(lvl) && lvl != '') return(lvl)
  }
  if (is.atomic(x)) {
    txt <- tolower(paste(x, collapse = ' '))
    if (grepl('complete',   txt)) return('Complete Genome')
    if (grepl('chromosome', txt)) return('Chromosome')
    if (grepl('scaffold',   txt)) return('Scaffold')
    if (grepl('contig',     txt)) return('Contig')
  }
  NA_character_
}

.extract_n50 <- function(x) {
  if (is.list(x)) return(x$contign50 %||% x$scaffoldn50 %||% NA_real_)
  NA_real_
}

.ASSEMBLY_WEIGHTS <- c(
'Complete Genome'          = 10,
'Complete Genome (latest)' = 10,
'Chromosome'               = 8,
'Chromosome level'         = 8,
'Scaffold'                 = 5,
'Scaffold level'           = 5,
'Contig'                   = 2,
'Contig level'             = 2
)

.score_species <- function(assembly, sra, biosample) {
  if (is.null(assembly) || (assembly$count %||% 0L) == 0L) {
    best_score  <- 0
    total_score <- 0
    richness    <- 0
  } else {
    IDS <- assembly$ids %||% character()
    if (!length(IDS)) {
      best_score  <- 0
      total_score <- 0
      richness    <- 0
    } else {
      SUMS <- .fetch_esummary_batched('assembly', IDS)
      LEVELS <- sapply(SUMS, .extract_assembly_level)
      STRUCT <- .ASSEMBLY_WEIGHTS[LEVELS]
      STRUCT[is.na(STRUCT)] <- 0
      N50 <- sapply(SUMS, .extract_n50)
      max_struct <- max(STRUCT)
      tied_idx   <- which(STRUCT == max_struct)
      best_idx <- if (length(tied_idx) > 1) {
        tied_idx[which.max(N50[tied_idx])]
      } else tied_idx
      best_score  <- STRUCT[best_idx]
      total_score <- sum(STRUCT)
      richness <- best_score + log1p(total_score - best_score)
    }
  }
  sra_score       <- 2 * log1p(sra$count       %||% 0L)
  biosample_score <-     log1p(biosample$count %||% 0L)
  score <- richness + sra_score + biosample_score
  list(
  score           = score,
  richness        = richness,
  sra_score       = sra_score,
  biosample_score = biosample_score
  )
}

# Ontology-driven classification

.normalise_strategy <- function(x) {
  if (is.null(x) || is.na(x)) return(NA_character_)
  x <- tolower(x)
  x <- trimws(x)
  x <- gsub('[_\\-\\/\\.]+', ' ', x)
  x <- gsub('[^a-z0-9 ]', ' ', x)
  x <- gsub('\\s+', ' ', x)
  trimws(x)
}

.extract_xml_tag <- function(expxml, tag) {
  if (is.null(expxml) || is.na(expxml) || !nzchar(expxml)) return(NA_character_)
  wrapped <- paste0('<ROOT>', expxml, '</ROOT>')
  doc <- tryCatch(
  xml2::read_xml(wrapped, options = c('NOBLANKS', 'NOERROR', 'NOWARNING')),
  error = function(e) NULL
  )
  if (is.null(doc)) return(NA_character_)
  xml2::xml_ns_strip(doc)
  tag <- toupper(tag)
  paths <- c(
  paste0('.//LIBRARY_DESCRIPTOR/', tag),
  paste0('.//Library_descriptor/', tag),
  paste0('.//', tag),
  paste0('.//', tolower(tag))
  )
  for (p in paths) {
    node <- xml2::xml_find_first(doc, p)
    if (!inherits(node, 'xml_missing')) {
      val <- xml2::xml_text(node)
      if (nzchar(trimws(val))) return(val)
    }
  }
  nodes <- xml2::xml_find_all(doc, './/*')
  hits  <- nodes[tolower(xml2::xml_name(nodes)) == tolower(tag)]
  if (length(hits) > 0) {
    val <- xml2::xml_text(hits[[1]])
    if (nzchar(trimws(val))) return(val)
  }
  NA_character_
}

.is_unknown <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(TRUE)
  x <- .normalise_strategy(x)
  x %in% c('unknown', 'not applicable', 'na', 'none', 'unspecified', 'unclassified')
}

.is_other_strategy <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(FALSE)
  .normalise_strategy(x) %in% c('other')
}

.extract_geo_accessions <- function(expxml) {
  if (is.null(expxml) || is.na(expxml) || !nzchar(expxml)) {
    return(list(GSE = character(), GSM = character()))
  }
  txt <- as.character(expxml)
  gse <- unique(unlist(regmatches(txt, gregexpr('\\bGSE\\d+\\b', txt, perl = TRUE))))
  gsm <- unique(unlist(regmatches(txt, gregexpr('\\bGSM\\d+\\b', txt, perl = TRUE))))
  gse <- gse[nzchar(gse)]
  gsm <- gsm[nzchar(gsm)]
  list(GSE = gse, GSM = gsm)
}

.collapse_acc <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  paste(sort(unique(x)), collapse = ';')
}

.is_geo_linked <- function(GSE, GSM) {
  length(GSE) > 0L || length(GSM) > 0L
}

.classify_strategy <- function(term, ontology) {
  if (is.na(term) || !nzchar(term) || .is_unknown(term)) {
    return(list(class = 'unknown', subclass = 'Unknown'))
  }
  for (class in names(ontology)) {
    for (subclass in names(ontology[[class]])) {
      if (term %in% ontology[[class]][[subclass]]) {
        return(list(class = class, subclass = subclass))
      }
    }
  }
  list(class = 'other', subclass = 'Other')
}

.classify_strategy_fallback <- function(strategy_raw,
strategy_norm,
source_raw,
selection_raw,
title_raw = NA_character_) {
  if (!.is_unknown(strategy_raw) && !.is_unknown(strategy_norm) && !.is_other_strategy(strategy_norm)) {
    return(list(class = 'other', subclass = 'Other'))
  }
  src <- .normalise_strategy(source_raw)
  sel <- .normalise_strategy(selection_raw)
  ttl <- .normalise_strategy(title_raw)
  if (!is.na(sel) && grepl('\\brestriction\\b|\\brestriction digest\\b|\\brestriction enzyme\\b|\\bgenotyping\\b|\\bgbs\\b|\\brad\\b|\\bdd\\s*rad\\b|\\bapeki\\b|\\bpst\\s*i\\b|\\bmsp\\s*i\\b|\\beco\\s*ri\\b', sel)) {
    return(list(class = 'genomic', subclass = 'RAD-seq'))
  }
  if (!is.na(src) && grepl('\\bgenomic\\b', src)) {
    if (!is.na(sel) && grepl('\\brestriction\\b|\\brestriction digest\\b', sel)) {
      return(list(class = 'genomic', subclass = 'RAD-seq'))
    }
    return(list(class = 'genomic', subclass = 'WGS'))
  }
  if (!is.na(src) && grepl('\\btranscriptomic\\b', src)) {
    if (!is.na(ttl) && grepl('\\bsmall\\b|\\bmirna\\b|\\bsirna\\b|\\bpirna\\b|\\blncrna\\b', ttl)) {
      return(list(class = 'transcriptomic', subclass = 'small-RNA'))
    }
    if (!is.na(ttl) && grepl('\\biso\\s*seq\\b|\\bisoseq\\b|\\bdirect rna\\b', ttl)) {
      return(list(class = 'transcriptomic', subclass = 'Long-read'))
    }
    return(list(class = 'transcriptomic', subclass = 'RNA-seq'))
  }
  if (!is.na(sel) && grepl('\\bpoly a\\b|\\bpolya\\b|\\bcdna\\b|\\bmrna\\b|\\bribo\\s*deplet\\b', sel)) {
    return(list(class = 'transcriptomic', subclass = 'RNA-seq'))
  }
  if (!is.na(src) && grepl('\\bepigenom', src)) {
    if (!is.na(ttl) && grepl('\\bcut tag\\b|\\bcutandtag\\b', ttl)) return(list(class = 'epigenomic', subclass = 'CUT&Tag'))
    if (!is.na(ttl) && grepl('\\bcut run\\b|\\bcutandrun\\b', ttl)) return(list(class = 'epigenomic', subclass = 'CUT&RUN'))
    if (!is.na(ttl) && grepl('\\batac\\b', ttl)) return(list(class = 'epigenomic', subclass = 'ATAC-seq'))
    if (!is.na(ttl) && grepl('\\bdnase\\b', ttl)) return(list(class = 'epigenomic', subclass = 'DNase-seq'))
    if (!is.na(ttl) && grepl('\\bfaire\\b', ttl)) return(list(class = 'epigenomic', subclass = 'FAIRE-seq'))
    if (!is.na(ttl) && grepl('\\bmnase\\b', ttl)) return(list(class = 'epigenomic', subclass = 'MNase-seq'))
    if (!is.na(ttl) && grepl('\\bchip\\b', ttl)) return(list(class = 'epigenomic', subclass = 'ChIP-seq'))
    if (!is.na(ttl) && grepl('\\bbisulfite\\b|\\bwgbs\\b|\\brrbs\\b|\\bmethyl\\b', ttl)) return(list(class = 'epigenomic', subclass = 'Bisulfite-seq'))
    return(list(class = 'epigenomic', subclass = 'Other'))
  }
  if (!is.na(ttl) && grepl('\\batac\\b|\\bdnase\\b|\\bfaire\\b|\\bmnase\\b|\\bchip\\b|\\bcut tag\\b|\\bcutandtag\\b|\\bcut run\\b|\\bcutandrun\\b', ttl)) {
    if (grepl('\\bcut tag\\b|\\bcutandtag\\b', ttl)) return(list(class = 'epigenomic', subclass = 'CUT&Tag'))
    if (grepl('\\bcut run\\b|\\bcutandrun\\b', ttl)) return(list(class = 'epigenomic', subclass = 'CUT&RUN'))
    if (grepl('\\batac\\b', ttl)) return(list(class = 'epigenomic', subclass = 'ATAC-seq'))
    if (grepl('\\bdnase\\b', ttl)) return(list(class = 'epigenomic', subclass = 'DNase-seq'))
    if (grepl('\\bfaire\\b', ttl)) return(list(class = 'epigenomic', subclass = 'FAIRE-seq'))
    if (grepl('\\bmnase\\b', ttl)) return(list(class = 'epigenomic', subclass = 'MNase-seq'))
    if (grepl('\\bchip\\b', ttl)) return(list(class = 'epigenomic', subclass = 'ChIP-seq'))
  }
  if (!is.na(ttl) && grepl('\\bbisulfite\\b|\\bwgbs\\b|\\brrbs\\b|\\bmethyl\\b', ttl)) {
    return(list(class = 'epigenomic', subclass = 'Bisulfite-seq'))
  }
  if (!is.na(src) && grepl('\\bchromatin\\b|\\bchromosome conformation\\b|\\bconformation\\b', src)) {
    if (!is.na(ttl) && grepl('\\bhi c\\b|\\bhic\\b', ttl)) return(list(class = 'chromatin', subclass = 'Hi-C'))
    if (!is.na(ttl) && grepl('\\bchia pet\\b', ttl)) return(list(class = 'chromatin', subclass = 'ChIA-PET'))
    if (!is.na(ttl) && grepl('\\btcc\\b', ttl)) return(list(class = 'chromatin', subclass = 'TCC'))
    if (!is.na(ttl) && grepl('\\b3c\\b|\\b4c\\b|\\b5c\\b|\\bcapture c\\b|\\bpromoter capture\\b|\\bhichip\\b|\\bplac\\b', ttl)) return(list(class = 'chromatin', subclass = '3C-based'))
    return(list(class = 'chromatin', subclass = 'Other'))
  }
  if (!is.na(ttl) && grepl('\\bhi c\\b|\\bhic\\b|\\bchia pet\\b|\\btcc\\b|\\b3c\\b|\\b4c\\b|\\b5c\\b|\\bcapture c\\b|\\bhichip\\b|\\bplac\\b|\\bpromoter capture\\b', ttl)) {
    if (grepl('\\bhi c\\b|\\bhic\\b', ttl)) return(list(class = 'chromatin', subclass = 'Hi-C'))
    if (grepl('\\bchia pet\\b', ttl)) return(list(class = 'chromatin', subclass = 'ChIA-PET'))
    if (grepl('\\btcc\\b', ttl)) return(list(class = 'chromatin', subclass = 'TCC'))
    return(list(class = 'chromatin', subclass = '3C-based'))
  }
  if (.is_unknown(strategy_raw) || .is_unknown(strategy_norm) || .is_other_strategy(strategy_norm)) {
    return(list(class = 'unknown', subclass = 'Unknown'))
  }
  list(class = 'other', subclass = 'Other')
}

.ONTOLOGY <- list(
  genomic = list(
    `WGS` = c(
      'wgs',
      'wga',
      'wcs',
      'wxs',
      'finishing',
      'whole genome sequencing',
      'whole genome seq',
      'whole genome shotgun',
      'whole genome shotgun sequencing',
      'whole genome resequencing',
      'genome sequencing',
      'genomic sequencing'
    ),
    `Amplicon-seq` = c(
      'amplicon',
      'amplicon seq',
      'amplicon sequencing',
      'targeted amplicon',
      'targeted amplicon sequencing',
      'amplicon based sequencing',
      '16s amplicon',
      '16s sequencing',
      '16s rrna sequencing',
      '18s amplicon',
      'its amplicon',
      'its sequencing',
      'metabarcoding'
    ),
    `RAD-seq` = c(
      'rad seq',
      'rad sequencing',
      'radseq',
      'restriction site associated dna sequencing',
      'restriction site associated sequencing',
      'dd rad',
      'dd radseq',
      'double digest radseq',
      'gbs',
      'genotyping by sequencing',
      'genotyping by seq'
    ),
    `Targeted-Capture` = c(
      'targeted capture',
      'targeted sequencing',
      'targeted seq',
      'targeted resequencing',
      'exome sequencing',
      'exome seq',
      'whole exome sequencing',
      'wes',
      'capture sequencing',
      'capture seq',
      'hybrid capture',
      'target enrichment',
      'panel sequencing',
      'gene panel'
    ),
    `Clone-based` = c(
      'clone',
      'cloneend',
      'poolclone',
      'clone end',
      'clone ends',
      'clone end sequencing',
      'clone based',
      'clone based sequencing',
      'bac end sequencing',
      'bac end seq',
      'fosmid end sequencing',
      'fosmid end seq',
      'cosmid end sequencing',
      'cosmid end seq',
      'pool clone',
      'pool clone sequencing'
    )
  ),
  transcriptomic = list(
    `RNA-seq` = c(
      'rna seq',
      'ssrna seq',
      'ribo seq',
      'fl cdna',
      'est',
      'rnaseq',
      'rna sequencing',
      'transcriptome sequencing',
      'transcriptome seq',
      'transcriptome profiling',
      'whole transcriptome sequencing',
      'whole transcriptome seq',
      'total rna seq',
      'total rna sequencing',
      'stranded rna seq',
      'strand specific rna seq',
      'poly a rna seq',
      'polya rna seq',
      'cdna sequencing',
      'cdna seq',
      'full length cdna',
      'full length transcript sequencing',
      'single cell rna seq',
      'scrna seq',
      'single nucleus rna seq',
      'snrna seq',
      'tag based rna seq',
      'tag seq',
      'digital gene expression',
      'dge'
    ),
    `small-RNA` = c(
      'ncrna seq',
      'mirna seq',
      'small rna seq',
      'small rna sequencing',
      'smrna seq',
      'smrna sequencing',
      'microrna sequencing',
      'microrna seq',
      'mirna sequencing',
      'non coding rna seq',
      'noncoding rna seq',
      'lncrna seq',
      'lncrna sequencing',
      'small noncoding rna seq',
      'small non coding rna seq',
      'sirna seq',
      'pirna seq'
    ),
    `Long-read` = c(
      'iso seq',
      'isoseq',
      'direct rna seq',
      'direct rna sequencing',
      'nanopore cdna',
      'nanopore cdna seq',
      'full length transcriptome'
    )
  ),
  epigenomic = list(
    `Bisulfite-seq` = c(
      'bisulfite seq',
      'mbd seq',
      'medip seq',
      'mre seq',
      'bisulfite sequencing',
      'whole genome bisulfite sequencing',
      'wgbs',
      'wgbss',
      'reduced representation bisulfite sequencing',
      'rrbs',
      'rrbss',
      'methylation sequencing',
      'methyl seq',
      'methylation seq',
      'methylome sequencing',
      'mbd sequencing',
      'medip sequencing',
      'mre sequencing'
    ),
    `ChIP-seq` = c(
      'chip seq',
      'rip seq',
      'chip sequencing',
      'chipseq',
      'chip dna seq',
      'chip dna sequencing',
      'chip exo',
      'chip exo seq',
      'chip exo sequencing',
      'rna immunoprecipitation sequencing',
      'rna immunoprecipitation seq',
      'rip sequencing'
    ),
    `CUT&RUN` = c(
      'cut run',
      'cutandrun',
      'cut run sequencing',
      'cutandrun sequencing'
    ),
    `CUT&Tag` = c(
      'cut tag',
      'cutandtag',
      'cut tag sequencing',
      'cutandtag sequencing'
    ),
    `ATAC-seq` = c(
      'atac seq',
      'atac sequencing',
      'atacseq',
      'assay for transposase accessible chromatin sequencing',
      'assay for transposase accessible chromatin seq',
      'single cell atac seq',
      'scatac seq',
      'single nucleus atac seq',
      'snatac seq'
    ),
    `DNase-seq` = c(
      'dnase hypersensitivity',
      'dnase seq',
      'dnase sequencing',
      'dnase i hypersensitivity',
      'dnase i hypersensitivity sequencing',
      'dnase i seq',
      'dnase i sequencing'
    ),
    `FAIRE-seq` = c(
      'faire seq',
      'faire sequencing',
      'formaldehyde assisted isolation of regulatory elements sequencing',
      'formaldehyde assisted isolation of regulatory elements seq'
    ),
    `MNase-seq` = c(
      'mnase seq',
      'mnase sequencing',
      'micrococcal nuclease sequencing',
      'micrococcal nuclease seq',
      'nucleosome mapping',
      'nucleosome positioning sequencing'
    ),
    `SELEX` = c(
      'selex',
      'selex seq',
      'selex sequencing',
      'ht selex',
      'high throughput selex'
    )
  ),
  chromatin = list(
    `Hi-C` = c(
      'hi c',
      'hic',
      'hi c seq',
      'hi c sequencing',
      'chromosome conformation capture carbon copy'
    ),
    `3C-based` = c(
      '3c',
      '4c',
      '5c',
      'capture c',
      'promoter capture c',
      'hichip',
      'plac seq',
      'plac-seq'
    ),
    `ChIA-PET` = c(
      'chia pet',
      'chia pet seq',
      'chia pet sequencing',
      'chromatin interaction analysis by paired end tag sequencing',
      'chromatin interaction analysis by paired end tag seq'
    ),
    `TCC` = c(
      'tethered chromatin conformation capture',
      'tcc',
      'tcc seq',
      'tcc sequencing',
      'tethered chromatin conformation capture sequencing'
    )
  ),
  other = list(
    `Other` = c(
      'custom',
      'custom sequencing',
      'custom protocol',
      'metagenomic',
      'metagenome sequencing',
      'metagenomic sequencing',
      'metatranscriptomic',
      'metatranscriptome sequencing',
      'metatranscriptomic sequencing',
      'proteomic',
      'proteomics',
      'metabolomic',
      'metabolomics',
      'chip array',
      'microarray',
      'genotyping array',
      'snp array',
      'array',
      'imaging',
      'optical mapping',
      'optical map',
      'validation',
      'pcr',
      'qpcr',
      'rt pcr',
      'rt qpcr',
      'library construction',
      'sequence',
      'sequencing',
      'test',
      'pilot',
      'control',
      'spike in',
      'spike in control'
    )
  )
)

.sra_metadata_core <- function(results, species = NULL) {
  if (!is.null(species)) {
    species <- intersect(species, names(results))
    if (!length(species)) stop('No matching species found')
    results <- results[species]
  }
  SRA_IDS <- unique(unlist(purrr::map(results, ~ .x$sra$ids %||% character())))
  if (!length(SRA_IDS)) {
    return(tibble::tibble(
    species       = character(),
    sra_id        = character(),
    strategy_raw  = character(),
    strategy_norm = character(),
    class         = character(),
    subclass      = character(),
    geo_linked    = logical(),
    gse_ids       = character(),
    gsm_ids       = character()
    ))
  }
  ID_TO_SPECIES <- unlist(purrr::map(names(results), function(sp) {
    ids <- results[[sp]]$sra$ids %||% character()
    if (!length(ids)) return(character())
    stats::setNames(rep(sp, length(ids)), ids)
  }))
  batch_size <- 200
  BATCHES <- split(SRA_IDS, ceiling(seq_along(SRA_IDS) / batch_size))
  pb <- utils::txtProgressBar(min = 0, max = length(BATCHES), style = 3)
  OUT <- list()
  for (i in seq_along(BATCHES)) {
    utils::setTxtProgressBar(pb, i)
    ids  <- BATCHES[[i]]
    SUMS <- .fetch_esummary_batched('sra', ids)
    SUMS <- .normalise_esummary_list(SUMS, ids)
    for (acc in names(SUMS)) {
      x <- SUMS[[acc]]
      strategy_raw  <- .extract_xml_tag(x$expxml, 'LIBRARY_STRATEGY')
      strategy_norm <- .normalise_strategy(strategy_raw)
      cls           <- .classify_strategy(strategy_norm, .ONTOLOGY)
      if (cls$class %in% c('unknown', 'other')) {
        rescue_ok <- (.is_unknown(strategy_raw) || .is_unknown(strategy_norm) || .is_other_strategy(strategy_norm))
        if (rescue_ok) {
          source_raw    <- .extract_xml_tag(x$expxml, 'LIBRARY_SOURCE')
          selection_raw <- .extract_xml_tag(x$expxml, 'LIBRARY_SELECTION')
          title_raw     <- .extract_xml_tag(x$expxml, 'TITLE')
          rescue <- .classify_strategy_fallback(
          strategy_raw  = strategy_raw,
          strategy_norm = strategy_norm,
          source_raw    = source_raw,
          selection_raw = selection_raw,
          title_raw     = title_raw
          )
          if (!is.null(rescue$class) && rescue$class %in% c('genomic', 'transcriptomic', 'epigenomic', 'chromatin')) {
            cls <- rescue
          } else if (.is_unknown(strategy_raw) || .is_unknown(strategy_norm) || .is_other_strategy(strategy_norm)) {
            cls <- list(class = 'unknown', subclass = 'Unknown')
          }
        }
      }
      if (cls$class == 'other' && (.is_unknown(strategy_raw) || .is_unknown(strategy_norm))) {
        cls <- list(class = 'unknown', subclass = 'Unknown')
      }
      geo <- .extract_geo_accessions(x$expxml)
      geo_linked <- .is_geo_linked(geo$GSE, geo$GSM)
      gse_ids <- .collapse_acc(geo$GSE)
      gsm_ids <- .collapse_acc(geo$GSM)
      sp <- ID_TO_SPECIES[[acc]] %||% NA_character_
      OUT[[acc]] <- tibble::tibble(
      species       = sp,
      sra_id        = acc,
      strategy_raw  = strategy_raw,
      strategy_norm = strategy_norm,
      class         = cls$class,
      subclass      = cls$subclass,
      geo_linked    = geo_linked,
      gse_ids       = gse_ids,
      gsm_ids       = gsm_ids
      )
    }
  }
  close(pb)
  dplyr::bind_rows(OUT)
}

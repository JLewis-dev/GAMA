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

.ASSEMBLY_CLASSES <- c('complete', 'chromosome', 'scaffold', 'contig')

.ASSEMBLY_WEIGHTS <- c(
complete   = 10,
chromosome = 8,
scaffold   = 5,
contig     = 2
)

.assembly_level_class <- function(x) {
  if (is.list(x) && !inherits(x, 'data.frame')) {
    x <- vapply(x, .flatten_to_char, character(1))
  }
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  txt <- tolower(x)
  hit <- !is.na(txt) & grepl('complete', txt)
  out[hit] <- 'complete'
  hit <- !is.na(txt) & grepl('chromosome', txt)
  out[hit] <- 'chromosome'
  hit <- !is.na(txt) & grepl('scaffold', txt)
  out[hit] <- 'scaffold'
  hit <- !is.na(txt) & grepl('contig', txt)
  out[hit] <- 'contig'
  out
}

.assembly_level_weight <- function(x) {
  cls <- .assembly_level_class(x)
  out <- unname(.ASSEMBLY_WEIGHTS[cls])
  out[is.na(out)] <- 0
  out
}

.best_assembly_n50 <- function(level, n50) {
  weight <- .assembly_level_weight(level)
  if (!length(weight) || all(weight <= 0)) return(NA_real_)
  n50 <- suppressWarnings(as.numeric(n50))
  best_weight <- max(weight, na.rm = TRUE)
  best_n50 <- n50[weight == best_weight]
  if (!length(best_n50) || all(is.na(best_n50))) return(NA_real_)
  max(best_n50, na.rm = TRUE)
}

.assembly_empty_metadata <- function(sp) {
  tibble::tibble(
  species       = sp,
  entrez_uid    = NA_character_,
  level         = NA_character_,
  n50           = NA_real_,
  coverage      = NA_real_,
  biosample     = NA_character_,
  bioproject    = NA_character_,
  submitter     = NA_character_,
  release_date  = NA_character_,
  ftp_path      = NA_character_
  )
}

.assembly_metadata_core <- function(results, species = NULL) {
  if (!is.null(species)) {
    species <- intersect(species, names(results))
    if (!length(species)) .gama_stop('No matching species found.')
    results <- results[species]
  }
  n <- length(results)
  if (!n) return(.assembly_empty_metadata(character())[0, ])
  HAS_ASSEMBLY <- vapply(results, function(x) .search_count(x$assembly) > 0L, logical(1))
  if (!any(HAS_ASSEMBLY)) return(dplyr::bind_rows(lapply(names(results), .assembly_empty_metadata)))
  batch_total <- sum(vapply(results[HAS_ASSEMBLY], function(x) {
    .search_batch_count(x$assembly, batch_size = 100)
  }, integer(1)))
  pb <- .pb_init(max(1L, batch_total))
  OUT <- list()
  idx <- 0L
  tick <- 0L
  append_sums <- function(SUMS, sp, OUT, idx) {
    if (!length(SUMS)) return(list(OUT = OUT, idx = idx))
    for (j in seq_along(SUMS)) {
      x <- SUMS[[j]]
      acc <- names(SUMS)[j]
      if (is.na(acc) || !nzchar(acc)) acc <- .esummary_uid(x)
      idx <- idx + 1L
      OUT[[idx]] <- tibble::tibble(
      species       = sp,
      entrez_uid    = acc,
      level         = .extract_assembly_level(x),
      n50           = .extract_n50(x),
      coverage      = as.numeric(.flatten_to_char(x$coverage)),
      biosample     = .flatten_to_char(x$biosampleaccn),
      bioproject    = .flatten_to_char(x$gb_bioprojects$bioprojectaccn),
      submitter     = .flatten_to_char(x$submitterorganization),
      release_date  = .flatten_to_char(x$asmreleasedate_genbank),
      ftp_path      = .flatten_to_char(x$ftppath_genbank)
      )
    }
    list(OUT = OUT, idx = idx)
  }
  for (i in seq_along(results)) {
    sp <- names(results)[i]
    asm <- results[[i]]$assembly
    idx_before <- idx
    if (!is.null(asm) && .search_count(asm) > 0L) {
      ids <- .search_ids(asm)
      count <- .search_count(asm)
      if (length(ids) && count <= length(ids)) {
        BATCHES <- split(ids, ceiling(seq_along(ids) / 100L))
        for (b in BATCHES) {
          SUMS <- .safe_entrez_summary('assembly', id = b)
          SUMS <- .normalise_esummary_list(SUMS, b)
          tmp <- append_sums(SUMS, sp, OUT, idx)
          OUT <- tmp$OUT
          idx <- tmp$idx
          tick <- tick + 1L
          .pb_tick(pb, tick)
        }
      } else if (.search_has_history(asm)) {
        STARTS <- seq.int(0L, count - 1L, by = 100L)
        for (start in STARTS) {
          size <- min(100L, count - start)
          SUMS <- .safe_entrez_summary(
          db          = 'assembly',
          web_history = asm$web_history,
          retstart    = start,
          retmax      = size
          )
          SUMS <- .normalise_esummary_history_list(SUMS)
          tmp <- append_sums(SUMS, sp, OUT, idx)
          OUT <- tmp$OUT
          idx <- tmp$idx
          tick <- tick + 1L
          .pb_tick(pb, tick)
        }
      } else if (length(ids)) {
        BATCHES <- split(ids, ceiling(seq_along(ids) / 100L))
        for (b in BATCHES) {
          SUMS <- .safe_entrez_summary('assembly', id = b)
          SUMS <- .normalise_esummary_list(SUMS, b)
          tmp <- append_sums(SUMS, sp, OUT, idx)
          OUT <- tmp$OUT
          idx <- tmp$idx
          tick <- tick + 1L
          .pb_tick(pb, tick)
        }
      }
    }
    if (idx == idx_before) {
      idx <- idx + 1L
      OUT[[idx]] <- .assembly_empty_metadata(sp)
    }
  }
  .pb_close(pb)
  dplyr::bind_rows(OUT)
}

.score_species <- function(assembly, sra, biosample) {
  if (is.null(assembly) || .search_count(assembly) == 0L) {
    best_score  <- 0
    total_score <- 0
    richness    <- 0
  } else {
    SUMS <- .fetch_search_summaries('assembly', assembly)
    if (!length(SUMS)) {
      best_score  <- 0
      total_score <- 0
      richness    <- 0
    } else {
      LEVELS <- sapply(SUMS, .extract_assembly_level)
      STRUCT <- .assembly_level_weight(LEVELS)
      N50 <- suppressWarnings(as.numeric(sapply(SUMS, .extract_n50)))
      max_struct <- max(STRUCT)
      tied_idx   <- which(STRUCT == max_struct)
      best_idx <- if (length(tied_idx) > 1) {
        tied_n50 <- N50[tied_idx]
        if (all(is.na(tied_n50))) tied_idx[1L] else tied_idx[which.max(tied_n50)]
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

# SRA ontology-driven classification

.normalise_strategy <- function(x) {
  if (is.null(x) || is.na(x)) return(NA_character_)
  x <- tolower(x)
  x <- trimws(x)
  x <- gsub('[_\\-\\/\\.]+', ' ', x)
  x <- gsub('[^a-z0-9 ]', ' ', x)
  x <- gsub('\\s+', ' ', x)
  trimws(x)
}

.extract_sra_ids <- function(x) {
  if (is.null(x) || (is.atomic(x) && !nzchar(as.character(x)))) {
    return(list(
    biosample  = NA_character_,
    bioproject = NA_character_
    ))
  }
  txt <- if (is.list(x)) {
    paste(unlist(x, recursive = TRUE, use.names = FALSE), collapse = ' ')
  } else {
    as.character(x)
  }
  txt <- as.character(txt)
  pat_sam <- '\\b(?:SAMN|SAMEA|SAMD)\\d+\\b'
  pat_prj <- '\\b(?:PRJNA|PRJEB|PRJDB|PRJDA)\\d+\\b'
  sam_matches <- unique(unlist(regmatches(txt, gregexpr(pat_sam, txt, perl = TRUE))))
  prj_matches <- unique(unlist(regmatches(txt, gregexpr(pat_prj, txt, perl = TRUE))))
  sam_val <- if (length(sam_matches)) sam_matches[1] else NA_character_
  prj_val <- if (length(prj_matches)) prj_matches[1] else NA_character_
  list(biosample = sam_val, bioproject = prj_val)
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
    return(list(class = 'unknown', subclass = 'unknown'))
  }
  for (class in names(ontology)) {
    for (subclass in names(ontology[[class]])) {
      if (term %in% ontology[[class]][[subclass]]) {
        return(list(class = class, subclass = subclass))
      }
    }
  }
  list(class = 'other', subclass = 'other')
}

.classify_strategy_fallback <- function(strategy_raw,
strategy_norm,
source_raw,
selection_raw,
title_raw = NA_character_) {
  if (!.is_unknown(strategy_raw) && !.is_unknown(strategy_norm) && !.is_other_strategy(strategy_norm)) {
    return(list(class = 'other', subclass = 'other'))
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
      return(list(class = 'transcriptomic', subclass = 'long-read'))
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
    if (!is.na(ttl) && grepl('\\bbisulfite\\b|\\bwgbs\\b|\\brrbs\\b|\\bmethyl\\b', ttl)) return(list(class = 'epigenomic', subclass = 'bisulfite-seq'))
    return(list(class = 'epigenomic', subclass = 'other'))
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
    return(list(class = 'epigenomic', subclass = 'bisulfite-seq'))
  }
  if (!is.na(src) && grepl('\\bchromatin\\b|\\bchromosome conformation\\b|\\bconformation\\b', src)) {
    if (!is.na(ttl) && grepl('\\bhi c\\b|\\bhic\\b', ttl)) return(list(class = 'chromatin', subclass = 'Hi-C'))
    if (!is.na(ttl) && grepl('\\bchia pet\\b', ttl)) return(list(class = 'chromatin', subclass = 'ChIA-PET'))
    if (!is.na(ttl) && grepl('\\btcc\\b', ttl)) return(list(class = 'chromatin', subclass = 'TCC'))
    if (!is.na(ttl) && grepl('\\b3c\\b|\\b4c\\b|\\b5c\\b|\\bcapture c\\b|\\bpromoter capture\\b|\\bhichip\\b|\\bplac\\b', ttl)) return(list(class = 'chromatin', subclass = '3C-based'))
    return(list(class = 'chromatin', subclass = 'other'))
  }
  if (!is.na(ttl) && grepl('\\bhi c\\b|\\bhic\\b|\\bchia pet\\b|\\btcc\\b|\\b3c\\b|\\b4c\\b|\\b5c\\b|\\bcapture c\\b|\\bhichip\\b|\\bplac\\b|\\bpromoter capture\\b', ttl)) {
    if (grepl('\\bhi c\\b|\\bhic\\b', ttl)) return(list(class = 'chromatin', subclass = 'Hi-C'))
    if (grepl('\\bchia pet\\b', ttl)) return(list(class = 'chromatin', subclass = 'ChIA-PET'))
    if (grepl('\\btcc\\b', ttl)) return(list(class = 'chromatin', subclass = 'TCC'))
    return(list(class = 'chromatin', subclass = '3C-based'))
  }
  if (.is_unknown(strategy_raw) || .is_unknown(strategy_norm) || .is_other_strategy(strategy_norm)) {
    return(list(class = 'unknown', subclass = 'unknown'))
  }
  list(class = 'other', subclass = 'other')
}

.ONTOLOGY <- list(
  genomic = list(
    `amplicon-seq` = c(
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
    `clone-based` = c(
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
    `targeted-capture` = c(
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
    )
  ),
  transcriptomic = list(
    `long-read` = c(
      'iso seq',
      'isoseq',
      'direct rna seq',
      'direct rna sequencing',
      'nanopore cdna',
      'nanopore cdna seq',
      'full length transcriptome'
    ),
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
    )
  ),
  epigenomic = list(
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
    `bisulfite-seq` = c(
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
    `Hi-C` = c(
      'hi c',
      'hic',
      'hi c seq',
      'hi c sequencing',
      'chromosome conformation capture carbon copy'
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
    `other` = c(
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

.sra_ontology_parameter_map <- function(ontology, target = 'class') {
  target <- .gama_validate_parameters(target, 'target', c('class', 'subclass'), multiple = FALSE, allow_null = FALSE)
  keys <- character()
  values <- character()
  add_parameter <- function(parameter, value) {
    key <- .gama_parameter_key(parameter)
    keep <- !is.na(key) & nzchar(key)
    if (!any(keep)) return(invisible(NULL))
    keys <<- c(keys, key[keep])
    values <<- c(values, rep(value, sum(keep)))
    invisible(NULL)
  }
  for (class in names(ontology)) {
    for (subclass in names(ontology[[class]])) {
      value <- if (identical(target, 'class')) class else subclass
      add_parameter(subclass, value)
      variants <- ontology[[class]][[subclass]] %||% character()
      if (length(variants) > 0L) add_parameter(variants, value)
    }
  }
  if (!length(keys)) return(stats::setNames(character(), character()))
  grouped <- split(values, keys)
  out <- vapply(grouped, function(x) {
    x <- unique(x)
    if (length(x) == 1L) x else NA_character_
  }, character(1))
  out[!is.na(out)]
}

.sra_metadata_core <- function(results, species = NULL) {
  if (!is.null(species)) {
    species <- intersect(species, names(results))
    if (!length(species)) .gama_stop('No matching species found.')
    results <- results[species]
  }
  HAS_SRA <- vapply(results, function(x) .search_count(x$sra) > 0L, logical(1))
  if (!any(HAS_SRA)) {
    return(tibble::tibble(
    species       = character(),
    entrez_uid    = character(),
    biosample     = character(),
    bioproject    = character(),
    strategy_raw  = character(),
    strategy_norm = character(),
    class         = character(),
    subclass      = character(),
    geo_linked    = logical(),
    gse_ids       = character(),
    gsm_ids       = character()
    ))
  }
  batch_total <- sum(vapply(results[HAS_SRA], function(x) {
    sra <- x$sra
    count <- .search_count(sra)
    ids <- .search_ids(sra)
    if (count == 0L) return(0L)
    if (length(ids) && count <= length(ids)) return(as.integer(ceiling(length(ids) / 100L)))
    if (.search_has_history(sra)) return(as.integer(ceiling(count / 100L)))
    if (length(ids)) return(as.integer(ceiling(length(ids) / 100L)))
    0L
  }, integer(1)))
  pb <- .pb_init(max(1L, batch_total))
  OUT <- list()
  idx <- 0L
  tick <- 0L
  append_sums <- function(SUMS, sp, OUT, idx) {
    if (!length(SUMS)) return(list(OUT = OUT, idx = idx))
    for (j in seq_along(SUMS)) {
      x <- SUMS[[j]]
      acc <- names(SUMS)[j]
      if (is.na(acc) || !nzchar(acc)) acc <- .esummary_uid(x)
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
            cls <- list(class = 'unknown', subclass = 'unknown')
          }
        }
      }
      if (cls$class == 'other' && (.is_unknown(strategy_raw) || .is_unknown(strategy_norm))) {
        cls <- list(class = 'unknown', subclass = 'unknown')
      }
      sra_ids <- .extract_sra_ids(x)
      geo <- .extract_geo_accessions(x$expxml)
      geo_linked <- .is_geo_linked(geo$GSE, geo$GSM)
      gse_ids <- .collapse_acc(geo$GSE)
      gsm_ids <- .collapse_acc(geo$GSM)
      idx <- idx + 1L
      OUT[[idx]] <- tibble::tibble(
      species       = sp,
      entrez_uid    = acc,
      biosample     = sra_ids$biosample,
      bioproject    = sra_ids$bioproject,
      strategy_raw  = strategy_raw,
      strategy_norm = strategy_norm,
      class         = cls$class,
      subclass      = cls$subclass,
      geo_linked    = geo_linked,
      gse_ids       = gse_ids,
      gsm_ids       = gsm_ids
      )
    }
    list(OUT = OUT, idx = idx)
  }
  for (i in seq_along(results)) {
    sp <- names(results)[i]
    sra <- results[[i]]$sra
    if (is.null(sra) || .search_count(sra) == 0L) next
    ids <- .search_ids(sra)
    if (length(ids) && .search_count(sra) <= length(ids)) {
      BATCHES <- split(ids, ceiling(seq_along(ids) / 100L))
      for (b in BATCHES) {
        SUMS <- .safe_entrez_summary('sra', id = b)
        SUMS <- .normalise_esummary_list(SUMS, b)
        tmp <- append_sums(SUMS, sp, OUT, idx)
        OUT <- tmp$OUT
        idx <- tmp$idx
        tick <- tick + 1L
        .pb_tick(pb, tick)
      }
      next
    }
    if (.search_has_history(sra)) {
      STARTS <- seq.int(0L, .search_count(sra) - 1L, by = 100L)
      for (start in STARTS) {
        size <- min(100L, .search_count(sra) - start)
        SUMS <- .safe_entrez_summary(
        db          = 'sra',
        web_history = sra$web_history,
        retstart    = start,
        retmax      = size
        )
        SUMS <- .normalise_esummary_history_list(SUMS)
        tmp <- append_sums(SUMS, sp, OUT, idx)
        OUT <- tmp$OUT
        idx <- tmp$idx
        tick <- tick + 1L
        .pb_tick(pb, tick)
      }
      next
    }
    if (length(ids)) {
      BATCHES <- split(ids, ceiling(seq_along(ids) / 100L))
      for (b in BATCHES) {
        SUMS <- .safe_entrez_summary('sra', id = b)
        SUMS <- .normalise_esummary_list(SUMS, b)
        tmp <- append_sums(SUMS, sp, OUT, idx)
        OUT <- tmp$OUT
        idx <- tmp$idx
        tick <- tick + 1L
        .pb_tick(pb, tick)
      }
    }
  }
  .pb_close(pb)
  dplyr::bind_rows(OUT)
}

# BioSample ontology-driven classification

.biosample_uid_from_input <- function(input_id) {
  x <- as.character(input_id %||% NA_character_)
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  if (grepl('^[0-9]+$', x)) return(x)
  NA_character_
}

.biosample_accession_from_input <- function(input_id) {
  x <- as.character(input_id %||% NA_character_)
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  if (grepl('^(?:SAMN|SAMEA|SAMD)\\d+$', x, perl = TRUE)) return(x)
  NA_character_
}

.biosample_first_match <- function(x, pattern) {
  if (is.null(x) || is.na(x) || !nzchar(as.character(x))) return(NA_character_)
  m <- regmatches(x, regexpr(pattern, x, perl = TRUE))
  if (!length(m) || is.na(m) || !nzchar(m)) return(NA_character_)
  as.character(m)
}

.biosample_normalise_name <- function(x) {
  if (length(x) == 0L) return(character())
  out <- as.character(x)
  out[is.na(out)] <- NA_character_
  out <- tolower(out)
  out <- gsub('[_\\-\\/\\\\]+', ' ', out)
  out <- gsub('[^a-z0-9 ]', ' ', out)
  out <- gsub('\\s+', ' ', out)
  out <- trimws(out)
  out[out == ''] <- NA_character_
  out
}

.biosample_normalise_value <- function(x) {
  if (length(x) == 0L) return(character())
  out <- as.character(x)
  out[is.na(out)] <- NA_character_
  out <- tolower(out)
  out <- gsub('https?://\\S+', ' ', out, perl = TRUE)
  out <- gsub('[\t\r\n]+', ' ', out)
  out <- gsub('[^a-z0-9 \\-_/\\.;:,\\(\\)\\[\\]]', ' ', out)
  out <- gsub('\\s+', ' ', out)
  out <- trimws(out)
  out[out == ''] <- NA_character_
  out
}

.biosample_tissue_normalise <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x <- gsub('[\t\r\n]+', ' ', x)
  x <- gsub('[\\(\\)\\[\\]\\{\\}]', ' ', x)
  x <- gsub('[^a-z0-9 ;,/| -]', ' ', x)
  x <- gsub('[-_]+', ' ', x)
  x <- gsub('\\s+', ' ', x)
  x <- gsub('\\b(?:efo\\s*\\d+|efo\\d+)\\b', ' ', x, perl = TRUE)
  x <- gsub('\\b\\d+(?:\\.\\d+)?\\s*(?:d|day|days|wk|wks|week|weeks|mo|mos|month|months|yr|yrs|year|years)\\b\\s*(?:old)?\\b', ' ', x, perl = TRUE)
  x <- gsub('\\b(?:d|day|days|wk|wks|week|weeks|mo|mos|month|months|yr|yrs|year|years)[- ]old\\b', ' ', x, perl = TRUE)
  x <- gsub('\\bv\\s*\\d+\\b', ' ', x, perl = TRUE)
  x <- gsub('\\bstage\\b', ' ', x, perl = TRUE)
  x <- gsub('\\b\\d+(?:st|nd|rd|th)\\b', ' ', x, perl = TRUE)
  x <- gsub('\\bleaf\\s*\\d+\\b', 'leaf', x, perl = TRUE)
  x <- gsub('\\bleaf\\d+\\b', 'leaf', x, perl = TRUE)
  x <- gsub('\\b\\d+(?:dap|das|dpi|hpi)\\b', ' ', x, perl = TRUE)
  x <- gsub('\\b(root|endosperm|kernel|seed|ear|tassel|husk|cob|petiole|stem|tissue)\\d+\\b', '\\1', x, perl = TRUE)
  x <- gsub('\\b\\d+\\b', ' ', x, perl = TRUE)
  x <- gsub('\\breplicates?\\b', ' ', x, perl = TRUE)
  x <- gsub('\\b(?:fresh|fully expanded|fully grown|field grown)\\b', ' ', x, perl = TRUE)
  x <- trimws(gsub('\\s+', ' ', x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

.biosample_missing_like_values <- function() {
  c(
    'na', 'n/a', 'null', '-', '--', '.',
    'not collected', 'not collect', 'not provided', 'not recorded',
    'not specified', 'not available', 'none', 'no tissue', 'unknown',
    'unspecified', 'missing', 'no data', 'no applicable',
    'nor applicable', 'not applicable', 'unavailable', 'not submitted',
    'not determined', 'n/ap'
  )
}

.biosample_is_missing_like_value <- function(x) {
  if (length(x) == 0L) return(logical())
  x <- .biosample_tissue_normalise(x)
  miss <- .biosample_missing_like_values()
  is.na(x) | !nzchar(x) | x %in% miss
}

.biosample_is_tissue_attribute <- function(attribute_name_norm, attribute_name_harmonised = NULL) {
  n <- .biosample_normalise_name(attribute_name_norm)
  h <- .biosample_normalise_name(attribute_name_harmonised)
  accepted <- c('tissue','organism part','cell type','tissue type','organ')
  (n %in% accepted) | (h %in% accepted)
}

.biosample_anatomy_class_levels <- function() {
  c('aerial', 'ground', 'reproductive', 'whole', 'in_vitro', 'other')
}

.biosample_collapse_anatomy_profile <- function(x, level = 'anatomy_class') {
  vals <- unique(as.character(x))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  levels <- if (identical(level, 'anatomy_subclass')) {
    .biosample_anatomy_subclass_levels()
  } else {
    .biosample_anatomy_class_levels()
  }
  vals <- vals[vals %in% c(levels, 'unknown')]
  if (length(vals) > 1L) return('mixed')
  if (length(vals) == 1L) return(vals)
  'unknown'
}

.biosample_anatomy_profile_levels <- function() {
  c(.biosample_anatomy_class_levels(), 'mixed', 'unknown')
}

.biosample_anatomy_profile_labels <- function(parse = FALSE) {
  labs <- c(
    aerial = 'aerial',
    ground = 'ground',
    reproductive = 'reproductive',
    whole = 'whole',
    in_vitro = if (isTRUE(parse)) 'italic(\'in vitro\')' else 'in vitro',
    other = 'other',
    mixed = 'mixed',
    unknown = 'unknown'
  )
  labs[.biosample_anatomy_profile_levels()]
}

.biosample_modality_levels <- function() {
  c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'mixed', 'unknown')
}

.biosample_anatomy_subclass_map <- function(ontology = .ANATOMY_ONTOLOGY) {
  rows <- lapply(
    names(ontology),
    function(cls) {
      tibble::tibble(
        anatomy_class = cls,
        anatomy_subclass = names(ontology[[cls]])
      )
    }
  )
  dplyr::bind_rows(rows)
}

.biosample_anatomy_subclass_levels <- function(ontology = .ANATOMY_ONTOLOGY) {
  unique(.biosample_anatomy_subclass_map(ontology)$anatomy_subclass)
}

.ANATOMY_ONTOLOGY <- list(
  aerial = list(
    leaf = list(
      leaf = list(
        variants = c('leaf', 'leaves', 'cauline leaf', 'cauline leaves', 'mature leaf', 'mature leaves', 'old leaf', 'trifoliate leaves', 'trifoliate', 'true leaf', 'true leaves', 'foliage', 'seedling leaf', 'seedling leaves', 'leaf tissue', 'leaf sample', 'leave', 'leaves tissue', 'leave tissue', 'leave tissues', 'leaaves', 'leafs', 'leat', 'plant leave', 'foliar', 'cauline'),
        ontology = list(namespace = 'PO', id = 'PO:0025034', label = 'leaf'),
        rank = 'organ'
      ),
      rosette_leaf = list(
        variants = c('leaf rosette', 'rosette', 'rosette and cauline leaves', 'rosette leaf', 'rosette leaves', 'rossette leaf', 'rosetteleaves'),
        ontology = list(namespace = 'PO', id = 'PO:0000014', label = 'rosette leaf'),
        rank = 'organ'
      ),
      flag_leaf = list(
        variants = c('flag leaf', 'flag leaves', 'flagleaf', 'flag leafs'),
        ontology = list(namespace = 'PO', id = 'PO:0020103', label = 'flag leaf'),
        rank = 'organ'
      ),
      young_leaf = list(
        variants = c('young leaf', 'young leaves', 'young leave', 'youngleaf'),
        ontology = list(namespace = 'PO', id = 'PO:0025034', label = 'leaf'),
        rank = 'organ'
      ),
      cotyledon = list(
        variants = c('cotyledon', 'cotyledons', 'cotyleden', 'cotyldon'),
        ontology = list(namespace = 'PO', id = 'PO:0020030', label = 'cotyledon'),
        rank = 'organ'
      ),
      leaf_blade = list(
        variants = c('blade', 'lamina', 'leaf blade', 'leaf lamina'),
        ontology = list(namespace = 'PO', id = 'PO:0020039', label = 'leaf lamina'),
        rank = 'suborgan'
      ),
      leaf_sheath = list(
        variants = c('leaf sheath', 'sheath', 'sheaths', 'leafsheath'),
        ontology = list(namespace = 'PO', id = 'PO:0020104', label = 'leaf sheath'),
        rank = 'suborgan'
      ),
      petiole = list(
        variants = c('petiole', 'petioles'),
        ontology = list(namespace = 'PO', id = 'PO:0020038', label = 'petiole'),
        rank = 'suborgan'
      ),
      leaflet = list(
        variants = c('leaflet', 'leaflets', 'terminal leaflet'),
        ontology = list(namespace = 'PO', id = 'PO:0020049', label = 'leaflet'),
        rank = 'suborgan'
      ),
      ligule = list(
        variants = c('ligule', 'ligules'),
        ontology = list(namespace = 'PO', id = 'PO:0020105', label = 'ligule'),
        rank = 'suborgan'
      ),
      auricle = list(
        variants = c('auricle', 'auricles', 'leaf auricle region'),
        ontology = list(namespace = 'PO', id = 'PO:0020106', label = 'leaf sheath auricle'),
        rank = 'suborgan'
      ),
      leaf_vein = list(
        variants = c('leaf vein', 'vein'),
        ontology = list(namespace = 'PO', id = 'PO:0020138', label = 'leaf lamina vein'),
        rank = 'suborgan'
      ),
      midvein = list(
        variants = c('midrib', 'midvein'),
        ontology = list(namespace = 'PO', id = 'PO:0020139', label = 'leaf midvein'),
        rank = 'suborgan'
      ),
      leaf_epidermis = list(
        variants = c('leaf epidermis', 'lower leaf epidermis', 'upper leaf epidermis'),
        ontology = list(namespace = 'PO', id = 'PO:0000047', label = 'leaf lamina epidermis'),
        rank = 'suborgan'
      ),
      leaf_pooled_numbered = list(
        variants = c('leaf pooled numbered'),
        ontology = list(namespace = 'PO', id = 'PO:0025034', label = 'leaf'),
        rank = 'meta'
      ),
      leaf_primordium = list(
        variants = c('leaf primordium', 'leaf primordia', 'vascular leaf primordium', 'vascular leaf primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0000017', label = 'vascular leaf primordium'),
        rank = 'suborgan'
      ),
      cotyledon_primordium = list(
        variants = c('cotyledon primordium', 'cotyledon primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0000015', label = 'cotyledon primordium'),
        rank = 'suborgan'
      )
    ),
    shoot_meristem = list(
      shoot_apex = list(
        variants = c('shoot apex', 'shoot tip', 'apical part', 'branch tip', 'middle branch tip', 'top branch tip', 'bottom branch tip'),
        ontology = list(namespace = 'PO', id = 'PO:0000037', label = 'shoot axis apex'),
        rank = 'suborgan'
      ),
      shoot_apical_meristem = list(
        variants = c('shoot apical meristem', 'vegetative shoot apical meristem', 'sam', 'shoot apical meristem sam', 'apical shoot meristem', 'shoot meristem'),
        ontology = list(namespace = 'PO', id = 'PO:0020148', label = 'shoot apical meristem'),
        rank = 'suborgan'
      ),
      axillary_meristem = list(
        variants = c('axillary meristem', 'axillae', 'axillary meristems', 'tiller primordium', 'tiller primordia', 'tillering primordium', 'tillering site'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      shoot_axis_apex = list(
        variants = c('shoot axis apex'),
        ontology = list(namespace = 'PO', id = 'PO:0000037', label = 'shoot axis apex'),
        rank = 'suborgan'
      )
    ),
    stem = list(
      aerial_part = list(
        variants = c('above ground tissue', 'above ground tissues', 'aerial', 'aerial part', 'aerial parts', 'whole aerial tissue', 'aerial tissue', 'shoot tissue', 'above-ground', 'above-ground plant', 'aboveground', 'aboveground tissue', 'aboveground tissues', 'aboveground plant parts', 'plant arial part', 'arial part', 'areal part', 'tiller above ground', 'green tissue', 'green tissues', 'green part', 'green parts', 'total green tissue', 'total green tissues', 'aboveground green tissue', 'above-ground green tissue', 'photosynthetic tissue', 'photosynthetic tissues'),
        ontology = list(namespace = 'PO', id = 'PO:0009006', label = 'shoot system'),
        rank = 'organ'
      ),
      stem = list(
        variants = c('stalk', 'stem', 'stem base', 'stem differentiating xylem', 'stem tip', 'stems', 'stalks', 'twig', 'twigs', 'branch', 'branches', 'culm', 'culms', 'shoot axis', 'shoot axis tissue', 'stem tissue', 'aerial stem', 'aerial stems', 'tiller', 'tillers', 'tiller base'),
        ontology = list(namespace = 'PO', id = 'PO:0009047', label = 'stem'),
        rank = 'organ'
      ),
      shoot = list(
        variants = c('aerial shoot', 'aerial shoots', 'developing leaves and shoot apex', 'shoot', 'shoots', 'the whole shoot above cotelydon', 'whole shoot', 'whole shoots', 'plumule', 'plumulel'),
        ontology = list(namespace = 'PO', id = 'PO:0009006', label = 'shoot system'),
        rank = 'organ'
      ),
      stolon = list(
        variants = c('stolon', 'stolons'),
        ontology = list(namespace = 'PO', id = 'PO:0003024', label = 'stolon'),
        rank = 'organ'
      ),
      internode = list(
        variants = c('first internode', 'internode', 'internodes', 'shoot internode', 'shoot internodes', 'stem internode', 'stem internodes'),
        ontology = list(namespace = 'PO', id = 'PO:0005005', label = 'shoot axis internode'),
        rank = 'organ'
      ),
      node = list(
        variants = c('first node', 'last node', 'node', 'coleoptilar nodes', 'tiller bases', 'nodes', 'shoot node', 'shoot nodes', 'stem node', 'stem nodes'),
        ontology = list(namespace = 'PO', id = 'PO:0005004', label = 'shoot axis node'),
        rank = 'organ'
      ),
      pith = list(
        variants = c('pith'),
        ontology = list(namespace = 'PO', id = 'PO:0006109', label = 'pith'),
        rank = 'suborgan'
      ),
      bark = list(
        variants = c('bark'),
        ontology = list(namespace = 'PO', id = 'PO:0004518', label = 'bark'),
        rank = 'suborgan'
      ),
      wood = list(
        variants = c('wood'),
        ontology = list(namespace = 'PO', id = 'PO:0005848', label = 'secondary xylem'),
        rank = 'suborgan'
      ),
      bud = list(
        variants = c('bud', 'dormant bud', 'buds', 'axillary buds', 'axilary buds', 'tiller buds', 'l1 tiller buds'),
        ontology = list(namespace = 'PO', id = 'PO:0000055', label = 'bud'),
        rank = 'organ'
      ),
      coleoptile = list(
        variants = c('coleoptile', 'coleoptiles', 'coleoptyle', 'coleotile', 'coleotiles', 'coleoptile tip', 'coleotile tip'),
        ontology = list(namespace = 'PO', id = 'PO:0020033', label = 'coleoptile'),
        rank = 'organ'
      ),
      hypocotyl = list(
        variants = c('hypocotyl', 'hypocotyls', 'apical hook'),
        ontology = list(namespace = 'PO', id = 'PO:0020100', label = 'hypocotyl'),
        rank = 'organ'
      ),
      epicotyl = list(
        variants = c('epicotyl'),
        ontology = list(namespace = 'PO', id = 'PO:0020035', label = 'epicotyl'),
        rank = 'organ'
      ),
      mesocotyl = list(
        variants = c('mesocotyl'),
        ontology = list(namespace = 'PO', id = 'PO:0020037', label = 'mesocotyl'),
        rank = 'organ'
      ),
      vascular_cambium = list(
        variants = c('cambium', 'vascular cambium'),
        ontology = list(namespace = 'PO', id = 'PO:0005597', label = 'vascular cambium'),
        rank = 'suborgan'
      ),
      cork_cambium = list(
        variants = c('cork cambium', 'phellogen'),
        ontology = list(namespace = 'PO', id = 'PO:0005599', label = 'cork cambium'),
        rank = 'suborgan'
      ),
      intercalary_meristem = list(
        variants = c('intercalary meristem'),
        ontology = list(namespace = 'PO', id = 'PO:0006073', label = 'intercalary meristem'),
        rank = 'suborgan'
      )
    )
  ),
  ground = list(
    root = list(
      root = list(
        variants = c('root', 'roots', 'tap root', 'taproot', 'primary root', 'lateral roots', 'hairy roots', 'storage root', 'whole root', 'whole roots', 'seedling root', 'seedling roots', 'firoots'),
        ontology = list(namespace = 'PO', id = 'PO:0009005', label = 'root'),
        rank = 'organ'
      ),
      fine_root = list(
        variants = c('fine root'),
        ontology = list(namespace = 'PO', id = 'PO:0009005', label = 'root'),
        rank = 'suborgan'
      ),
      lateral_root = list(
        variants = c('lateral root'),
        ontology = list(namespace = 'PO', id = 'PO:0020121', label = 'lateral root'),
        rank = 'suborgan'
      ),
      adventitious_root = list(
        variants = c('adventitious root', 'adventitious roots'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      crown_root = list(
        variants = c('crown root', 'crown roots'),
        ontology = list(namespace = 'PO', id = 'PO:0000043', label = 'crown root'),
        rank = 'suborgan'
      ),
      root_hair = list(
        variants = c('root hair', 'root hairs', 'roothair'),
        ontology = list(namespace = 'PO', id = 'PO:0000256', label = 'root hair cell'),
        rank = 'cell'
      ),
      root_epidermis = list(
        variants = c('root epidermis'),
        ontology = list(namespace = 'PO', id = 'PO:0006036', label = 'root epidermis'),
        rank = 'suborgan'
      ),
      cortex = list(
        variants = c('cortex', 'root cortex'),
        ontology = list(namespace = 'PO', id = 'PO:0005708', label = 'cortex'),
        rank = 'suborgan'
      ),
      endodermis = list(
        variants = c('endodermis'),
        ontology = list(namespace = 'PO', id = 'PO:0000252', label = 'endodermis'),
        rank = 'suborgan'
      ),
      pericycle = list(
        variants = c('pericycle', 'root pericycle'),
        ontology = list(namespace = 'PO', id = 'PO:0006203', label = 'pericycle'),
        rank = 'suborgan'
      ),
      stele = list(
        variants = c('root stele', 'stele'),
        ontology = list(namespace = 'PO', id = 'PO:0025197', label = 'stele'),
        rank = 'suborgan'
      ),
      root_nodule = list(
        variants = c('root nodule', 'nodule', 'nodules'),
        ontology = list(namespace = 'PO', id = 'PO:0003023', label = 'root nodule'),
        rank = 'suborgan'
      ),
      radicle = list(
        variants = c('radicle'),
        ontology = list(namespace = 'PO', id = 'PO:0020031', label = 'radicle'),
        rank = 'organ'
      ),
      seminal_root = list(
        variants = c('seminal root', 'seminal roots'),
        ontology = list(namespace = 'PO', id = 'PO:0000046', label = 'seminal root'),
        rank = 'suborgan'
      ),
      shoot_borne_root = list(
        variants = c('shoot borne root', 'shoot borne roots'),
        ontology = list(namespace = 'PO', id = 'PO:0000042', label = 'shoot-borne root'),
        rank = 'suborgan'
      ),
      prop_root = list(
        variants = c('prop root', 'prop roots'),
        ontology = list(namespace = 'PO', id = 'PO:0000044', label = 'prop root'),
        rank = 'suborgan'
      )
    ),
    root_meristem = list(
      root_tip = list(
        variants = c('primary root tip', 'root tip', 'root tips', 'elongation zone', 'roottip'),
        ontology = list(namespace = 'PO', id = 'PO:0000025', label = 'root tip'),
        rank = 'suborgan'
      ),
      root_apex = list(
        variants = c('root apex'),
        ontology = list(namespace = 'PO', id = 'PO:0020147', label = 'root apical meristem'),
        rank = 'suborgan'
      ),
      root_apical_meristem = list(
        variants = c('root apical meristem', 'ram', 'root meristem', 'quiescent center', 'quiescent centre'),
        ontology = list(namespace = 'PO', id = 'PO:0020147', label = 'root apical meristem'),
        rank = 'suborgan'
      ),
      root_cap = list(
        variants = c('central root cap', 'lateral root cap', 'root cap', 'columella'),
        ontology = list(namespace = 'PO', id = 'PO:0020123', label = 'root cap'),
        rank = 'suborgan'
      ),
      coleorhiza = list(
        variants = c('coleorhiza', 'coleorhiza hair'),
        ontology = list(namespace = 'PO', id = 'PO:0020034', label = 'coleorhiza'),
        rank = 'suborgan'
      ),
      root_primordium = list(
        variants = c('root primordium', 'root primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0005029', label = 'root primordium'),
        rank = 'suborgan'
      ),
      lateral_root_primordium = list(
        variants = c('lateral root primordium', 'lateral root primordia', 'pericyclic lateral root primordium', 'non-pericyclic lateral root primordium'),
        ontology = list(namespace = 'PO', id = 'PO:0000016', label = 'lateral root primordium'),
        rank = 'suborgan'
      )
    ),
    storage = list(
      tuber = list(
        variants = c('mature tuber', 'tuber', 'tuber pith', 'tubers', 'young tuber', 'tuber flesh', 'tuber skin'),
        ontology = list(namespace = 'PO', id = 'PO:0025522', label = 'tuber'),
        rank = 'organ'
      ),
      bulb = list(
        variants = c('bulb'),
        ontology = list(namespace = 'PO', id = 'PO:0025356', label = 'bulb'),
        rank = 'organ'
      ),
      rhizome = list(
        variants = c('rhizome'),
        ontology = list(namespace = 'PO', id = 'PO:0004542', label = 'rhizome'),
        rank = 'organ'
      ),
      corm = list(
        variants = c('corm', 'corms'),
        ontology = list(namespace = 'PO', id = 'PO:0025355', label = 'corm'),
        rank = 'organ'
      )
    )
  ),
  reproductive = list(
    flower = list(
      flower = list(
        variants = c('floral tissue', 'flower', 'flowers', 'young flower', 'inflorescence flower', 'floral'),
        ontology = list(namespace = 'PO', id = 'PO:0009046', label = 'flower'),
        rank = 'organ'
      ),
      receptacle = list(
        variants = c('receptacle', 'mid receptacle'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      hypanthium = list(
        variants = c('hypanthium'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      inflorescence = list(
        variants = c('ear inflorescence', 'female inflorescence', 'inflorescence', 'inflorescence stem', 'inflorescences', 'immature inflorescence', 'infloresence', 'inflorescense', 'inflorescenses', 'inflorensences', 'immature infloresence', 'immature inflorescense', 'influorescence'),
        ontology = list(namespace = 'PO', id = 'PO:0009049', label = 'inflorescence'),
        rank = 'organ'
      ),
      panicle = list(
        variants = c('panicle', 'panicles', 'young panicle', 'panical', 'panilce'),
        ontology = list(namespace = 'PO', id = 'PO:0030123', label = 'panicle inflorescence'),
        rank = 'organ'
      ),
      spike = list(
        variants = c('spike', 'spikes', 'young spike', 'young spikes'),
        ontology = list(namespace = 'PO', id = 'PO:0030117', label = 'spike inflorescence'),
        rank = 'organ'
      ),
      rachis = list(
        variants = c('rachis'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      nectary = list(
        variants = c('nectary', 'nectaries'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      awn = list(
        variants = c('awn', 'awns'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      tassel = list(
        variants = c('immature tassel', 'tassel', 'tassel primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0020126', label = 'tassel inflorescence'),
        rank = 'organ'
      ),
      ear = list(
        variants = c('ear', 'immature ear', 'ears', 'immature ears', 'ears ea', 'cob', 'cobs'),
        ontology = list(namespace = 'PO', id = 'PO:0020136', label = 'ear inflorescence'),
        rank = 'organ'
      ),
      peduncle = list(
        variants = c('peduncle', 'shank', 'ear shank'),
        ontology = list(namespace = 'PO', id = 'PO:0009053', label = 'peduncle'),
        rank = 'suborgan'
      ),
      bract = list(
        variants = c('bract', 'husk', 'husks'),
        ontology = list(namespace = 'PO', id = 'PO:0009055', label = 'bract'),
        rank = 'suborgan'
      ),
      pedicel = list(
        variants = c('pedicel', 'pedicels'),
        ontology = list(namespace = 'PO', id = 'PO:0030112', label = 'pedicel'),
        rank = 'suborgan'
      ),
      sepal = list(
        variants = c('sepal', 'sepals'),
        ontology = list(namespace = 'PO', id = 'PO:0009031', label = 'sepal'),
        rank = 'organ'
      ),
      petal = list(
        variants = c('petal', 'petals'),
        ontology = list(namespace = 'PO', id = 'PO:0009032', label = 'petal'),
        rank = 'organ'
      ),
      stamen = list(
        variants = c('stamen', 'stamens'),
        ontology = list(namespace = 'PO', id = 'PO:0009029', label = 'stamen'),
        rank = 'organ'
      ),
      filament = list(
        variants = c('filament', 'stamen filament'),
        ontology = list(namespace = 'PO', id = 'PO:0009067', label = 'filament'),
        rank = 'suborgan'
      ),
      anther = list(
        variants = c('anther', 'anthers'),
        ontology = list(namespace = 'PO', id = 'PO:0009066', label = 'anther'),
        rank = 'suborgan'
      ),
      pistil = list(
        variants = c('gynoecium', 'pistil', 'pistils'),
        ontology = list(namespace = 'PO', id = 'PO:0009062', label = 'gynoecium'),
        rank = 'suborgan'
      ),
      carpel = list(
        variants = c('carpel', 'carpels'),
        ontology = list(namespace = 'PO', id = 'PO:0009030', label = 'carpel'),
        rank = 'organ'
      ),
      stigma = list(
        variants = c('stigma', 'stigmas'),
        ontology = list(namespace = 'PO', id = 'PO:0009073', label = 'stigma'),
        rank = 'suborgan'
      ),
      style = list(
        variants = c('style', 'styles'),
        ontology = list(namespace = 'PO', id = 'PO:0009074', label = 'style'),
        rank = 'suborgan'
      ),
      ovary = list(
        variants = c('ovary', 'ovaries'),
        ontology = list(namespace = 'PO', id = 'PO:0009072', label = 'plant ovary'),
        rank = 'suborgan'
      ),
      ovule = list(
        variants = c('ovule', 'ovules'),
        ontology = list(namespace = 'PO', id = 'PO:0020003', label = 'plant ovule'),
        rank = 'organ'
      ),
      nucellus = list(
        variants = c('nucellus', 'nucelli'),
        ontology = list(namespace = 'PO', id = 'PO:0020020', label = 'nucellus'),
        rank = 'suborgan'
      ),
      integument = list(
        variants = c('integument', 'integuments'),
        ontology = list(namespace = 'PO', id = 'PO:0020021', label = 'plant ovule integument'),
        rank = 'suborgan'
      ),
      silk = list(
        variants = c('silk', 'silks'),
        ontology = list(namespace = 'PO', id = 'PO:0006488', label = 'silk'),
        rank = 'suborgan'
      ),
      spikelet = list(
        variants = c('spikelet', 'young spikelet', 'spikelets', 'spiklet', 'spikerlet'),
        ontology = list(namespace = 'PO', id = 'PO:0009051', label = 'spikelet'),
        rank = 'organ'
      ),
      floret = list(
        variants = c('floret', 'florets'),
        ontology = list(namespace = 'PO', id = 'PO:0009082', label = 'spikelet floret'),
        rank = 'organ'
      ),
      flower_head = list(
        variants = c('flower head', 'flower heads', 'seed head', 'seed heads', 'capitulum', 'capitula'),
        ontology = NULL,
        rank = 'organ'
      ),
      lemma = list(
        variants = c('lemma', 'lemmas'),
        ontology = list(namespace = 'PO', id = 'PO:0009037', label = 'lemma'),
        rank = 'suborgan'
      ),
      palea = list(
        variants = c('palea', 'palea tissue', 'paleas'),
        ontology = list(namespace = 'PO', id = 'PO:0009038', label = 'palea'),
        rank = 'suborgan'
      ),
      glume = list(
        variants = c('glume', 'glumes'),
        ontology = list(namespace = 'PO', id = 'PO:0009039', label = 'glume'),
        rank = 'suborgan'
      ),
      lodicule = list(
        variants = c('lodicule', 'lodicules'),
        ontology = list(namespace = 'PO', id = 'PO:0009036', label = 'lodicule'),
        rank = 'suborgan'
      ),
      pollen = list(
        variants = c('male gametophyte', 'mature pollen', 'pollen', 'pollen grain', 'pollen tube', 'microgametophyte vegetative cell', 'vegetative cell'),
        ontology = list(namespace = 'PO', id = 'PO:0025281', label = 'pollen'),
        rank = 'cell'
      ),
      microspore = list(
        variants = c('microspore', 'microspores'),
        ontology = list(namespace = 'PO', id = 'PO:0020048', label = 'microspore'),
        rank = 'cell'
      ),
      meiocyte = list(
        variants = c('meiocyte', 'meiocytes', 'male meiocytes', 'isolated meiocytes', 'sporocyte', 'sporocytes'),
        ontology = list(namespace = 'PO', id = 'PO:0006204', label = 'sporocyte'),
        rank = 'cell'
      ),
      egg_cell = list(
        variants = c('egg cell', 'egg'),
        ontology = NULL,
        rank = 'cell'
      ),
      sperm_cell = list(
        variants = c('sperm cell', 'sperm cells', 'sperm'),
        ontology = NULL,
        rank = 'cell'
      ),
      zygote = list(
        variants = c('zygote'),
        ontology = NULL,
        rank = 'cell'
      ),
      central_cell = list(
        variants = c('central cell'),
        ontology = NULL,
        rank = 'cell'
      ),
      embryo_sac = list(
        variants = c('embryo sac', 'female gametophyte'),
        ontology = list(namespace = 'PO', id = 'PO:0025074', label = 'embryo sac'),
        rank = 'organ'
      ),
      floral_meristem = list(
        variants = c('floral meristem', 'flower meristem'),
        ontology = list(namespace = 'PO', id = 'PO:0000229', label = 'flower meristem'),
        rank = 'suborgan'
      ),
      inflorescence_meristem = list(
        variants = c('inflorescence apex', 'inflorescence meristem'),
        ontology = list(namespace = 'PO', id = 'PO:0000230', label = 'inflorescence meristem'),
        rank = 'suborgan'
      ),
      flower_bud = list(
        variants = c('flower bud', 'flower buds', 'floral bud', 'floral buds', 'axillary flower bud', 'terminal flower bud', 'closed buds', 'buds floral'),
        ontology = list(namespace = 'PO', id = 'PO:0000056', label = 'flower bud'),
        rank = 'organ'
      ),
      inflorescence_bud = list(
        variants = c('inflorescence bud', 'inflorescence buds', 'axillary inflorescence bud', 'terminal inflorescence bud'),
        ontology = list(namespace = 'PO', id = 'PO:0000057', label = 'inflorescence bud'),
        rank = 'organ'
      ),
      anther_wall = list(
        variants = c('anther wall', 'pollen sac wall'),
        ontology = list(namespace = 'PO', id = 'PO:0000002', label = 'anther wall'),
        rank = 'suborgan'
      ),
      pollen_sac = list(
        variants = c('pollen sac', 'pollen sacs', 'microsporangium'),
        ontology = list(namespace = 'PO', id = 'PO:0025277', label = 'pollen sac'),
        rank = 'suborgan'
      ),
      stamen_primordium = list(
        variants = c('stamen primordium', 'stamen primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0004705', label = 'stamen primordium'),
        rank = 'suborgan'
      ),
      carpel_primordium = list(
        variants = c('carpel primordium', 'carpel primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0004703', label = 'carpel primordium'),
        rank = 'suborgan'
      ),
      petal_primordium = list(
        variants = c('petal primordium', 'petal primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0000021', label = 'petal primordium'),
        rank = 'suborgan'
      ),
      sepal_primordium = list(
        variants = c('sepal primordium', 'sepal primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0004704', label = 'sepal primordium'),
        rank = 'suborgan'
      ),
      ovule_primordium = list(
        variants = c('ovule primordium', 'ovule primordia'),
        ontology = list(namespace = 'PO', id = 'PO:0000018', label = 'ovule primordium'),
        rank = 'suborgan'
      ),
      gynoecium_primordium = list(
        variants = c('gynoecium primordium', 'gynoecium primordia', 'pistil primordium'),
        ontology = list(namespace = 'PO', id = 'PO:0000019', label = 'gynoecium primordium'),
        rank = 'suborgan'
      )
    ),
    fruit = list(
      fruit = list(
        variants = c('fruit', 'fruits'),
        ontology = list(namespace = 'PO', id = 'PO:0009001', label = 'fruit'),
        rank = 'organ'
      ),
      capsule = list(
        variants = c('capsule', 'capsules'),
        ontology = NULL,
        rank = 'organ'
      ),
      pod = list(
        variants = c('pod', 'pods', 'bean pod', 'bean pods', 'pea pod', 'pea pods', 'legume pod', 'legume pods'),
        ontology = list(namespace = 'PO', id = 'PO:0030100', label = 'legume fruit'),
        rank = 'organ'
      ),
      silique = list(
        variants = c('silique', 'siliques'),
        ontology = list(namespace = 'PO', id = 'PO:0030106', label = 'silique fruit'),
        rank = 'organ'
      ),
      pericarp = list(
        variants = c('fruit pericarp', 'pericarp', 'pericarps'),
        ontology = list(namespace = 'PO', id = 'PO:0009084', label = 'pericarp'),
        rank = 'suborgan'
      ),
      exocarp = list(
        variants = c('exocarp', 'fruit exocarp'),
        ontology = list(namespace = 'PO', id = 'PO:0009085', label = 'exocarp'),
        rank = 'suborgan'
      ),
      mesocarp = list(
        variants = c('fruit mesocarp', 'mesocarp'),
        ontology = list(namespace = 'PO', id = 'PO:0009087', label = 'mesocarp'),
        rank = 'suborgan'
      ),
      endocarp = list(
        variants = c('endocarp', 'fruit endocarp'),
        ontology = list(namespace = 'PO', id = 'PO:0009086', label = 'endocarp'),
        rank = 'suborgan'
      ),
      fruit_flesh = list(
        variants = c('fruit flesh', 'pulp', 'fruit cortex', 'fruit pulp', 'flesh'),
        ontology = list(namespace = 'PO', id = 'PO:0025037', label = 'fruit storage parenchyma'),
        rank = 'suborgan'
      ),
      fruit_peel = list(
        variants = c('fruit peel', 'fruit skin', 'peel', 'skin'),
        ontology = list(namespace = 'PO', id = 'PO:0009085', label = 'exocarp'),
        rank = 'suborgan'
      ),
      locular_tissue = list(
        variants = c('locular tissue', 'locule', 'locules'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      fruit_placenta = list(
        variants = c('placenta', 'fruit placenta'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      fruit_septum = list(
        variants = c('septum', 'fruit septum'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      fruit_valve = list(
        variants = c('fruit valve', 'fruit valves', 'valve'),
        ontology = list(namespace = 'PO', id = 'PO:0000033', label = 'fruit valve'),
        rank = 'suborgan'
      )
    ),
    seed = list(
      seed = list(
        variants = c('seed', 'seeds', 'whole germinating seeds', 'whole seed', 'whole seeds', 'pea', 'peas', 'bean', 'beans', 'pulse', 'pulses', 'legume seed', 'legume seeds'),
        ontology = list(namespace = 'PO', id = 'PO:0009010', label = 'seed'),
        rank = 'organ'
      ),
      caryopsis = list(
        variants = c('caryopsis', 'caryposis'),
        ontology = list(namespace = 'PO', id = 'PO:0030104', label = 'caryopsis fruit'),
        rank = 'organ'
      ),
      grain = list(
        variants = c('grain', 'whole grain', 'grains', 'immature grains', 'developing grains'),
        ontology = list(namespace = 'PO', id = 'PO:0030104', label = 'caryopsis fruit'),
        rank = 'organ'
      ),
      kernel = list(
        variants = c('kernal', 'kernel', 'kernels', 'young kernel'),
        ontology = list(namespace = 'PO', id = 'PO:0030104', label = 'caryopsis fruit'),
        rank = 'organ'
      ),
      seed_coat = list(
        variants = c('immature seed coat', 'seed coat', 'seed coat epidermis', 'seedcoat', 'testa'),
        ontology = list(namespace = 'PO', id = 'PO:0009088', label = 'seed coat'),
        rank = 'suborgan'
      ),
      embryo = list(
        variants = c('embryo', 'embryo axis', 'embryonic axis', 'somatic embryo', 'embryos', 'mature embryos', 'immature embryos', 'germinating embryos', 'imamture embryos', 'immature embyro', 'immature zygotic embryos', 'culured immature embryos', 'cultured immature embryos'),
        ontology = list(namespace = 'PO', id = 'PO:0009009', label = 'plant embryo'),
        rank = 'organ'
      ),
      endosperm = list(
        variants = c('basal endosperm transfer layer', 'endopserm', 'endosperm', 'endosperm cell', 'maternal transfer zone', 'endospern', 'endosprm', 'endopsperm', 'endorsperm', 'endosprem'),
        ontology = list(namespace = 'PO', id = 'PO:0009089', label = 'endosperm'),
        rank = 'organ'
      ),
      aleurone_layer = list(
        variants = c('aleurone', 'aleuron layer', 'aleurone cell', 'aleurone layer'),
        ontology = list(namespace = 'PO', id = 'PO:0005360', label = 'aleurone layer'),
        rank = 'suborgan'
      ),
      scutellum = list(
        variants = c('apical scutellum', 'scutellum', 'scutella'),
        ontology = list(namespace = 'PO', id = 'PO:0020110', label = 'scutellum'),
        rank = 'organ'
      ),
      embryo_proper = list(
        variants = c('embryo proper', 'plant embryo proper'),
        ontology = list(namespace = 'PO', id = 'PO:0000001', label = 'plant embryo proper'),
        rank = 'suborgan'
      )
    )
  ),
  whole = list(
    whole = list(
      whole_plant = list(
        variants = c('whole', 'whole organism', 'whole plant', 'whole plants', 'total plant', 'total plants', 'plant', 'full plants', 'the entire plant', 'whole ground plant', 'whole body', 't0 plants'),
        match = 'exact',
        ontology = list(namespace = 'PO', id = 'PO:0000003', label = 'whole plant'),
        rank = 'organ'
      ),
      whole_seedling = list(
        variants = c('whole seedling', 'whole seedlings', 'whole seedlings including roots', 'total seedlings', 'total seedling', 'whole seeding', 'whole seedings', 'whole seedllings', 'whole seedlngs', 'whole seedlilng'),
        match = 'exact',
        ontology = list(namespace = 'PO', id = 'PO:0008037', label = 'seedling'),
        rank = 'organ'
      ),
      seedling = list(
        variants = c('seedling', 'seedlings', 'etiolated seedling', 'etiolated seedlings', 'young seedling', 'young seedlings', 'plantlet', 'plantlets', 'seeding', 'seedings', 'sprout', 'sprouts', 'seedligns', 'seedlngs', 'seedllings', 'seedlilng'),
        match = 'exact',
        ontology = list(namespace = 'PO', id = 'PO:0008037', label = 'seedling'),
        rank = 'organ'
      )
    )
  ),
  in_vitro = list(
    in_vitro = list(
      tissue_culture = list(
        variants = c('cell culture', 'tissue culture', 'culture', 'cultured cells', 'in vitro', 'gel culture'),
        ontology = list(namespace = 'PO', id = 'PO:0000004', label = 'in vitro plant structure'),
        rank = 'organ'
      ),
      callus = list(
        variants = c('callus', 'callus tissue', 'embryogenic callus', 'calli', 'dedifferentiated cellular calli', 'trangenic calli', 'transgenic calli', 'cllus', 'callusun'),
        ontology = list(namespace = 'PO', id = 'PO:0005052', label = 'plant callus'),
        rank = 'organ'
      ),
      cell_suspension = list(
        variants = c('cell suspension', 'cell suspension culture', 'suspension cells', 'suspension cell', 'suspension cell cultures'),
        ontology = list(namespace = 'PO', id = 'PO:0000005', label = 'cultured plant cell'),
        rank = 'cell'
      ),
      protoplast = list(
        variants = c('guard cell protoplast', 'leaf protoplast', 'mesophyll protoplast', 'mysophyll protoplast', 'protoplast', 'protoplasts'),
        ontology = list(namespace = 'PO', id = 'PO:0000006', label = 'plant protoplast'),
        rank = 'cell'
      )
    )
  ),
  other = list(
    other = list(
      generic_meristem = list(
        variants = c('apical meristem', 'meristem', 'apex', 'apices', 'apical', 'plant apex', 'apical region', 'dissected meristems', 'meristems', 'meristematic cells', 'meristematic zone'),
        ontology = list(namespace = 'PO', id = 'PO:0009013', label = 'generic meristem'),
        rank = 'suborgan'
      ),
      generic_axis = list(
        variants = c('axis', 'plant axis'),
        ontology = NULL,
        rank = 'meta'
      ),
      vascular = list(
        variants = c('vascular', 'vascular bundle', 'vascular bundles', 'vasculature', 'vascular tissue'),
        ontology = list(namespace = 'PO', id = 'PO:0005020', label = 'vascular bundle'),
        rank = 'suborgan'
      ),
      xylem = list(
        variants = c('developing xylem', 'xylem'),
        ontology = list(namespace = 'PO', id = 'PO:0005352', label = 'xylem'),
        rank = 'suborgan'
      ),
      phloem = list(
        variants = c('phloem'),
        ontology = list(namespace = 'PO', id = 'PO:0005417', label = 'phloem'),
        rank = 'suborgan'
      ),
      epidermis = list(
        variants = c('epidermis', 'epidermis including guard cells'),
        ontology = list(namespace = 'PO', id = 'PO:0005679', label = 'plant epidermis'),
        rank = 'suborgan'
      ),
      bundle_sheath = list(
        variants = c('bundle sheath', 'bundle sheath strand'),
        ontology = list(namespace = 'PO', id = 'PO:0006023', label = 'bundle sheath'),
        rank = 'suborgan'
      ),
      bundle_sheath_cell = list(
        variants = c('bundle sheath cell', 'bundle sheath cells'),
        ontology = list(namespace = 'PO', id = 'PO:0025541', label = 'bundle sheath cell'),
        rank = 'cell'
      ),
      guard_cell = list(
        variants = c('guard cell', 'guard cells'),
        ontology = list(namespace = 'PO', id = 'PO:0000293', label = 'guard cell'),
        rank = 'cell'
      ),
      mesophyll_cell = list(
        variants = c('leaf mesophyll', 'mesophyll', 'mesophyll cell', 'mesophyll cells'),
        ontology = list(namespace = 'PO', id = 'PO:0004006', label = 'mesophyll cell'),
        rank = 'cell'
      ),
      trichome = list(
        variants = c('glandular trichome', 'glandular trichomes', 'leaf trichome', 'leaf trichomes', 'stem trichome', 'stem trichomes', 'trichome', 'trichomes'),
        ontology = list(namespace = 'PO', id = 'PO:0000282', label = 'trichome'),
        rank = 'cell'
      ),
      somatic_cell = list(
        variants = c('somatic cell', 'plant cell', 'plant cells', 'somatic cells', 'somatic', 'cell', 'cells', 'all cell types', 'multiple cell-types'),
        ontology = list(namespace = 'PO', id = 'PO:0009002', label = 'plant cell'),
        rank = 'cell'
      ),
      chloroplast = list(
        variants = c('chloroplast'),
        ontology = list(namespace = 'GO', id = 'GO:0009507', label = 'chloroplast'),
        rank = 'organelle'
      ),
      nucleus = list(
        variants = c('nucleus', 'nuclei', 'nuclear'),
        ontology = list(namespace = 'GO', id = 'GO:0005634', label = 'nucleus'),
        rank = 'organelle'
      ),
      crown = list(
        variants = c('crown', 'collet'),
        ontology = NULL,
        rank = 'suborgan'
      ),
      multi_tissue = list(
        variants = c('leaf and shoot', 'multi tissue', 'leaf shoot', 'flower and leaves', 'multiple', 'mixed', 'various', 'mixed tissue', 'mixed tissues', 'multiple tissues', 'various tissues', 'whole tissue', 'total tissue', 'whole parts', 'leaf and stem', 'leaves and stems'),
        patterns = c('\\bmultiple tissues?\\b', '\\bmixed tissues?\\b', '\\bvarious tissues?\\b', '\\bmultiple plant tissues?\\b', '\\bseveral tissues?\\b'),
        ontology = NULL,
        rank = 'meta'
      ),
      bulk_sample = list(
        variants = c('bulk sample', 'bulk', 'bulk tissue', 'pooled sample', 'pooled tissue'),
        patterns = c('^bulk tissues?\\b', '^bulk samples?\\b', '^pooled tissues?\\b', '^pooled samples?\\b'),
        ontology = NULL,
        rank = 'meta'
      ),
      unspecified_tissue = list(
        variants = c('unspecified tissue', 'other', 'others', 'misc', 'miscellaneous', 'tissue', 'tissue sample', 'sample tissue', 'plant sample', 'all', 'total', 'vegetative', 'vegetative tissue', 'unspecified plant tissue', 'whole tisues'),
        ontology = NULL,
        rank = 'meta'
      ),
      environmental_sample = list(
        variants = c('rhizosphere', 'rhizosphere soil', 'soil', 'bulk soil', 'root associated soil', 'bulksoil'),
        ontology = NULL,
        rank = 'meta'
      ),
      plant_ontology_term = list(
        variants = c('plant ontology term'),
        ontology = NULL,
        rank = 'meta'
      ),
      mitochondrion = list(
        variants = c('mitochondrion', 'mitochondria', 'mitochondrial'),
        ontology = list(namespace = 'GO', id = 'GO:0005739', label = 'mitochondrion'),
        rank = 'organelle'
      )
    )
  ),
  unknown = list(
    unknown = list(
      unknown = list(
        variants = character(),
        ontology = NULL,
        rank = 'meta'
      )
    )
  )
)

.biosample_anatomy_label <- function(anatomy_term, node) {
  if (!is.null(node$ontology) && !is.null(node$ontology$label) && nzchar(node$ontology$label)) {
    return(as.character(node$ontology$label))
  }
  gsub('_', ' ', as.character(anatomy_term), fixed = TRUE)
}

.biosample_anatomy_ontology_parameter_map <- function(target = 'term', ontology = .ANATOMY_ONTOLOGY) {
  target <- .gama_validate_parameters(
    target,
    'target',
    c('class', 'subclass', 'term'),
    multiple = FALSE,
    allow_null = FALSE
  )
  keys <- character()
  values <- character()
  add_parameter <- function(parameter, value) {
    parameter <- as.character(parameter %||% character())
    if (!length(parameter)) return(invisible(NULL))
    key <- .gama_parameter_key(parameter)
    keep <- !is.na(parameter) & nzchar(trimws(parameter)) & !is.na(key) & nzchar(key)
    if (!any(keep)) return(invisible(NULL))
    keys <<- c(keys, key[keep])
    values <<- c(values, rep(as.character(value), sum(keep)))
    invisible(NULL)
  }
  for (anatomy_class in names(ontology)) {
    for (anatomy_subclass in names(ontology[[anatomy_class]])) {
      terms <- ontology[[anatomy_class]][[anatomy_subclass]]
      if (identical(target, 'class')) {
        add_parameter(anatomy_class, anatomy_class)
        add_parameter(anatomy_subclass, anatomy_class)
      } else if (identical(target, 'subclass')) {
        add_parameter(anatomy_subclass, anatomy_subclass)
      }
      for (anatomy_term in names(terms)) {
        node <- terms[[anatomy_term]]
        value <- switch(
          target,
          class = anatomy_class,
          subclass = anatomy_subclass,
          term = anatomy_term
        )
        add_parameter(anatomy_term, value)
        add_parameter(.biosample_anatomy_label(anatomy_term, node), value)
        if (!is.null(node$ontology) && !is.null(node$ontology$label)) {
          add_parameter(node$ontology$label, value)
        }
        add_parameter(node$variants %||% character(), value)
      }
    }
  }
  if (!length(keys)) return(stats::setNames(character(), character()))
  keep <- !is.na(keys) & nzchar(keys) & !is.na(values) & nzchar(values)
  keys <- keys[keep]
  values <- values[keep]
  if (!length(keys)) return(stats::setNames(character(), character()))
  out <- values
  names(out) <- keys
  out[!duplicated(paste(names(out), out, sep = '\r'))]
}

.compile_biosample_anatomy_variant_map <- function(ontology = .ANATOMY_ONTOLOGY) {
  variant_rows <- list()
  pattern_rows <- list()
  for (anatomy_class in names(ontology)) {
    for (anatomy_subclass in names(ontology[[anatomy_class]])) {
      terms <- ontology[[anatomy_class]][[anatomy_subclass]]
      for (anatomy_term in names(terms)) {
        node <- terms[[anatomy_term]]
        variants <- unique(.biosample_tissue_normalise(node$variants %||% character()))
        variants <- variants[!is.na(variants) & nzchar(variants)]
        rank <- as.character(node$rank %||% 'suborgan')
        match_mode <- as.character(node$match %||% if (as.character(anatomy_class) %in% c('whole', 'other')) 'exact' else 'phrase')
        if (length(variants)) {
          variant_rows[[length(variant_rows) + 1L]] <- data.frame(
            term_norm = variants,
            anatomy_term = as.character(anatomy_term),
            anatomy_class = as.character(anatomy_class),
            rank = rank,
            match_mode = match_mode,
            stringsAsFactors = FALSE
          )
        }
        patterns <- unique(as.character(node$patterns %||% character()))
        patterns <- patterns[!is.na(patterns) & nzchar(patterns)]
        if (length(patterns)) {
          pattern_rows[[length(pattern_rows) + 1L]] <- data.frame(
            pattern = patterns,
            anatomy_term = as.character(anatomy_term),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  map <- if (length(variant_rows)) do.call(rbind, variant_rows) else data.frame(
    term_norm = character(),
    anatomy_term = character(),
    anatomy_class = character(),
    rank = character(),
    match_mode = character(),
    stringsAsFactors = FALSE
  )
  pattern_map <- if (length(pattern_rows)) do.call(rbind, pattern_rows) else data.frame(
    pattern = character(),
    anatomy_term = character(),
    stringsAsFactors = FALSE
  )
  if (!nrow(map)) {
    attr(map, 'max_n') <- 1L
    attr(map, 'idx') <- stats::setNames(character(), character())
    attr(map, 'idx_phrase') <- stats::setNames(character(), character())
    attr(map, 'patterns') <- pattern_map
    return(map)
  }
  map$term_norm <- .biosample_tissue_normalise(map$term_norm)
  map <- map[!is.na(map$term_norm) & nzchar(map$term_norm), , drop = FALSE]
  map$match_mode[!map$match_mode %in% c('exact', 'phrase')] <- 'phrase'
  if (!nrow(map)) {
    attr(map, 'max_n') <- 1L
    attr(map, 'idx') <- stats::setNames(character(), character())
    attr(map, 'idx_phrase') <- stats::setNames(character(), character())
    attr(map, 'patterns') <- pattern_map
    return(map)
  }
  amb <- tapply(map$anatomy_term, map$term_norm, function(x) length(unique(x)) > 1L)
  if (any(amb)) {
    bad <- names(amb)[which(amb)[1]]
    .gama_stop('Anatomy ontology variant collision after normalisation: ', bad)
  }
  map <- map[order(map$anatomy_term, map$term_norm), , drop = FALSE]
  map <- map[!duplicated(map$term_norm), , drop = FALSE]
  map$word_n <- lengths(strsplit(map$term_norm, ' +', perl = TRUE))
  phrase <- map[map$match_mode == 'phrase', , drop = FALSE]
  attr(map, 'max_n') <- if (nrow(phrase)) max(phrase$word_n, na.rm = TRUE) else 1L
  attr(map, 'idx') <- stats::setNames(map$anatomy_term, map$term_norm)
  attr(map, 'idx_phrase') <- stats::setNames(phrase$anatomy_term, phrase$term_norm)
  attr(map, 'patterns') <- pattern_map
  map
}

.compile_biosample_anatomy_lexicon <- function(ontology = .ANATOMY_ONTOLOGY) {
  rows <- list()
  for (anatomy_class in names(ontology)) {
    for (anatomy_subclass in names(ontology[[anatomy_class]])) {
      terms <- ontology[[anatomy_class]][[anatomy_subclass]]
      for (anatomy_term in names(terms)) {
        node <- terms[[anatomy_term]]
        ontology_namespace <- if (!is.null(node$ontology)) as.character(node$ontology$namespace %||% NA_character_) else NA_character_
        ontology_id <- if (!is.null(node$ontology)) as.character(node$ontology$id %||% NA_character_) else NA_character_
        ontology_label <- if (!is.null(node$ontology)) as.character(node$ontology$label %||% NA_character_) else NA_character_
        rows[[length(rows) + 1L]] <- data.frame(
          anatomy_term = as.character(anatomy_term),
          label = .biosample_anatomy_label(anatomy_term, node),
          anatomy_class = as.character(anatomy_class),
          anatomy_subclass = as.character(anatomy_subclass),
          ontology_namespace = ontology_namespace,
          ontology_id = ontology_id,
          ontology_label = ontology_label,
          rank = as.character(node$rank %||% 'suborgan'),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  lex <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(lex)) return(lex)
  rank_order <- c(organ = 1, suborgan = 2, cell = 3, organelle = 4, meta = 5)
  lex$rank_score <- unname(rank_order[lex$rank])
  lex$rank_score[is.na(lex$rank_score)] <- 99
  lex <- lex[order(lex$anatomy_term, lex$rank_score), , drop = FALSE]
  lex <- lex[!duplicated(lex$anatomy_term), , drop = FALSE]
  lex$rank_score <- NULL
  lex
}

.match_biosample_anatomy_value <- function(x, variant_map) {
  x0 <- .biosample_tissue_normalise(x)
  if (is.na(x0) || !nzchar(x0) || is.null(variant_map) || !nrow(variant_map)) return(character())
  idx <- attr(variant_map, 'idx')
  if (is.null(idx) || !length(idx)) idx <- stats::setNames(variant_map$anatomy_term, variant_map$term_norm)
  idx_phrase <- attr(variant_map, 'idx_phrase')
  if (is.null(idx_phrase)) {
    phrase <- variant_map[variant_map$match_mode == 'phrase', , drop = FALSE]
    idx_phrase <- stats::setNames(phrase$anatomy_term, phrase$term_norm)
  }
  pattern_map <- attr(variant_map, 'patterns')
  if (is.null(pattern_map)) pattern_map <- data.frame(pattern = character(), anatomy_term = character(), stringsAsFactors = FALSE)
  max_n <- as.integer(attr(variant_map, 'max_n') %||% 1L)
  max_n <- max(1L, min(8L, max_n))
  hits <- character()
  h0 <- unname(idx[x0])
  if (!is.na(h0) && nzchar(h0)) hits <- c(hits, h0)
  if (nrow(pattern_map)) {
    use <- vapply(pattern_map$pattern, function(p) grepl(p, x0, perl = TRUE), logical(1))
    if (any(use)) hits <- c(hits, pattern_map$anatomy_term[use])
  }
  parts <- unlist(strsplit(x0, '\\s*(?:;|,|/|\\||\\band\\b)\\s*', perl = TRUE))
  parts <- unique(parts[!is.na(parts) & nzchar(parts)])
  if (!length(parts)) return(unique(hits[!is.na(hits) & nzchar(hits)]))
  h_parts <- unname(idx[parts])
  h_parts <- h_parts[!is.na(h_parts) & nzchar(h_parts)]
  if (length(h_parts)) hits <- c(hits, h_parts)
  if (!length(idx_phrase)) return(unique(hits[!is.na(hits) & nzchar(hits)]))
  for (p in parts) {
    toks <- unlist(strsplit(p, ' ', fixed = TRUE))
    toks <- toks[!is.na(toks) & nzchar(toks)]
    if (!length(toks)) next
    avail <- rep(TRUE, length(toks))
    for (n in seq.int(max_n, 1L, by = -1L)) {
      if (length(toks) < n) next
      for (i in seq_len(length(toks) - n + 1L)) {
        if (!all(avail[i:(i + n - 1L)])) next
        phrase <- paste(toks[i:(i + n - 1L)], collapse = ' ')
        h <- unname(idx_phrase[phrase])
        if (!is.na(h) && nzchar(h)) {
          hits <- c(hits, h)
          avail[i:(i + n - 1L)] <- FALSE
        }
      }
    }
  }
  unique(hits[!is.na(hits) & nzchar(hits)])
}

.BIOSAMPLE_ANATOMY_CACHE <- new.env(parent = emptyenv())

.biosample_anatomy_ref <- function() {
  if (!exists('variant_map', envir = .BIOSAMPLE_ANATOMY_CACHE, inherits = FALSE)) {
    assign('variant_map', .compile_biosample_anatomy_variant_map(), envir = .BIOSAMPLE_ANATOMY_CACHE)
  }
  if (!exists('lex', envir = .BIOSAMPLE_ANATOMY_CACHE, inherits = FALSE)) {
    assign('lex', .compile_biosample_anatomy_lexicon(), envir = .BIOSAMPLE_ANATOMY_CACHE)
  }
  list(
    variant_map = get('variant_map', envir = .BIOSAMPLE_ANATOMY_CACHE, inherits = FALSE),
    lex = get('lex', envir = .BIOSAMPLE_ANATOMY_CACHE, inherits = FALSE)
  )
}

.biosample_empty_attributes_unique <- function() {
  data.frame(
    source_channel = character(),
    node_name = character(),
    attribute_name_raw = character(),
    attribute_name_norm = character(),
    attribute_name_harmonised = character(),
    attribute_value_raw = character(),
    attribute_value_norm = character(),
    attribute_unit_raw = character(),
    attribute_index = integer(),
    stringsAsFactors = FALSE
  )
}

.biosample_empty_attributes <- function() {
  data.frame(
    species = character(),
    input_id = character(),
    entrez_uid = character(),
    biosample_id = character(),
    parse_status = character(),
    source_channel = character(),
    node_name = character(),
    attribute_name_raw = character(),
    attribute_name_norm = character(),
    attribute_name_harmonised = character(),
    attribute_value_raw = character(),
    attribute_value_norm = character(),
    attribute_unit_raw = character(),
    attribute_index = integer(),
    stringsAsFactors = FALSE
  )
}

.biosample_empty_records <- function(keep_xml = FALSE) {
  out <- data.frame(
    species = character(),
    input_id = character(),
    entrez_uid = character(),
    biosample_id = character(),
    accession_present = logical(),
    parse_status = character(),
    sampledata_present = logical(),
    sampledata_xml_ok = logical(),
    n_attributes = integer(),
    bioproject = character(),
    title_raw = character(),
    description_raw = character(),
    organism_raw = character(),
    package_raw = character(),
    stringsAsFactors = FALSE
  )
  if (isTRUE(keep_xml)) out$sampledata_xml <- character()
  out
}

.biosample_make_diagnostics <- function(requested_species, n_species, n_species_id_pairs, n_unique_input_ids, n_cached, n_fetched, parse_status_counts) {
  list(
    module = '.biosample_metadata_core',
    requested_species = as.character(requested_species %||% character()),
    n_species = as.integer(n_species %||% 0L),
    n_species_id_pairs = as.integer(n_species_id_pairs %||% 0L),
    n_unique_input_ids = as.integer(n_unique_input_ids %||% 0L),
    n_cached = as.integer(n_cached %||% 0L),
    n_fetched = as.integer(n_fetched %||% 0L),
    parse_status_counts = parse_status_counts %||% stats::setNames(integer(), character())
  )
}

.biosample_normalise_cache <- function(x) {
  if (is.null(x)) return(stats::setNames(vector('list', 0L), character()))
  if (!is.list(x)) .gama_stop('BioSample core: `esummary_cache` must be a list or NULL.')
  nm <- names(x) %||% character()
  if (any(!nzchar(nm))) .gama_stop('BioSample core: `esummary_cache` must be a *named* list (names = input_id).')
  x
}

.biosample_collect_id_map <- function(results) {
  if (is.null(results) || !is.list(results) || !length(results)) return(data.frame(species = character(), input_id = character(), stringsAsFactors = FALSE))
  spp <- names(results) %||% character()
  out <- list()
  for (sp in spp) {
    ids <- results[[sp]]$biosample$ids %||% character()
    ids <- as.character(ids)
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (!length(ids)) next
    out[[length(out) + 1L]] <- data.frame(species = rep(sp, length(unique(ids))), input_id = unique(ids), stringsAsFactors = FALSE)
  }
  if (!length(out)) return(data.frame(species = character(), input_id = character(), stringsAsFactors = FALSE))
  out <- do.call(rbind, out)
  out <- out[!duplicated(paste(out$species, out$input_id)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.biosample_extract_record_fields <- function(x, input_id) {
  uid_guess <- .flatten_to_char(x$uid %||% x$id %||% x$entrez_uid)
  uid_guess <- if (!is.na(uid_guess) && nzchar(uid_guess)) uid_guess else .biosample_uid_from_input(input_id)
  acc_guess <- .flatten_to_char(x$accession %||% x$accessionname %||% x$biosample %||% x$sample_accession)
  if (is.na(acc_guess) || !nzchar(acc_guess)) acc_guess <- .biosample_accession_from_input(input_id)
  if (is.na(acc_guess) || !nzchar(acc_guess)) {
    payload_txt <- tryCatch(paste(unlist(x, recursive = TRUE, use.names = FALSE), collapse = ' '), error = function(e) '')
    acc_guess <- .biosample_first_match(payload_txt, '\\b(?:SAMN|SAMEA|SAMD)\\d+\\b')
  }
  list(
    entrez_uid = uid_guess %||% NA_character_,
    biosample_id = acc_guess %||% NA_character_,
    title_raw = .flatten_to_char(x$title %||% x$sample_name %||% x$samplename),
    description_raw = .flatten_to_char(x$description %||% x$desc),
    organism_raw = .flatten_to_char(x$organism %||% x$organismname %||% x$taxonomy_name),
    package_raw = .flatten_to_char(x$package %||% x$package_name %||% x$model),
    sampledata_xml = .flatten_to_char(x$sampledata %||% x$sample_data %||% x$sampledataxml)
  )
}

.biosample_read_sampledata_xml <- function(sampledata_xml) {
  if (is.null(sampledata_xml) || is.na(sampledata_xml) || !nzchar(sampledata_xml)) return(NULL)
  sampledata_xml <- as.character(sampledata_xml)
  sampledata_xml <- gsub('^\\s*<\\?xml[^>]*\\?>\\s*', '', sampledata_xml, perl = TRUE)
  if (length(sampledata_xml) < 1L || is.na(sampledata_xml[1]) || !nzchar(trimws(sampledata_xml[1]))) return(NULL)
  wrapped <- paste0('<ROOT>', sampledata_xml, '</ROOT>')
  doc <- tryCatch(xml2::read_xml(wrapped, options = c('RECOVER','NONET','NOBLANKS','NOERROR','NOWARNING')), error = function(e) NULL)
  if (is.null(doc)) return(NULL)
  xml2::xml_ns_strip(doc)
  doc
}

.biosample_extract_bioproject <- function(doc) {
  if (is.null(doc)) return(NA_character_)
  links <- xml2::xml_find_all(doc, './/Link')
  vals <- character()
  if (length(links)) {
    target <- tolower(xml2::xml_attr(links, 'target'))
    label <- tolower(xml2::xml_attr(links, 'label'))
    type <- tolower(xml2::xml_attr(links, 'type'))
    target[is.na(target)] <- ''
    label[is.na(label)] <- ''
    type[is.na(type)] <- ''
    use <- grepl('bioproject', target) | grepl('bioproject', label) | grepl('bioproject', type)
    vals <- xml2::xml_text(links[use], trim = TRUE)
  }
  txt <- paste(c(vals, xml2::xml_text(doc, trim = TRUE)), collapse = ' ')
  hit <- regmatches(txt, gregexpr('PRJ[A-Z][A-Z0-9]+', txt, perl = TRUE))[[1]]
  hit <- unique(hit[!is.na(hit) & nzchar(hit)])
  if (!length(hit)) return(NA_character_)
  paste(sort(hit), collapse = '; ')
}

.biosample_extract_attribute_rows <- function(doc) {
  attrs <- xml2::xml_find_all(doc, './/Attribute')
  if (!length(attrs)) return(.biosample_empty_attributes_unique())
  attr_name <- xml2::xml_attr(attrs, 'attribute_name')
  harm_name <- xml2::xml_attr(attrs, 'harmonized_name')
  disp_name <- xml2::xml_attr(attrs, 'display_name')
  units_raw <- xml2::xml_attr(attrs, 'unit')
  units_raw2 <- xml2::xml_attr(attrs, 'units')
  vals <- xml2::xml_text(attrs, trim = TRUE)
  name_raw <- attr_name
  need_name <- is.na(name_raw) | !nzchar(trimws(name_raw))
  name_raw[need_name] <- harm_name[need_name]
  need_name <- is.na(name_raw) | !nzchar(trimws(name_raw))
  name_raw[need_name] <- disp_name[need_name]
  unit_final <- units_raw
  need_unit <- is.na(unit_final) | !nzchar(trimws(unit_final))
  unit_final[need_unit] <- units_raw2[need_unit]
  out <- data.frame(
    source_channel = rep('attribute', length(attrs)),
    node_name = rep('Attribute', length(attrs)),
    attribute_name_raw = as.character(name_raw),
    attribute_name_norm = .biosample_normalise_name(as.character(name_raw)),
    attribute_name_harmonised = .biosample_normalise_name(as.character(harm_name)),
    attribute_value_raw = as.character(vals),
    attribute_value_norm = .biosample_normalise_value(as.character(vals)),
    attribute_unit_raw = as.character(unit_final),
    attribute_index = seq_along(attrs),
    stringsAsFactors = FALSE
  )
  out$attribute_name_raw <- trimws(out$attribute_name_raw)
  out$attribute_value_raw <- trimws(out$attribute_value_raw)
  out$attribute_unit_raw <- trimws(out$attribute_unit_raw)
  out$attribute_name_raw[out$attribute_name_raw == ''] <- NA_character_
  out$attribute_value_raw[out$attribute_value_raw == ''] <- NA_character_
  out$attribute_unit_raw[out$attribute_unit_raw == ''] <- NA_character_
  keep <- !(is.na(out$attribute_name_raw) & is.na(out$attribute_value_raw))
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.biosample_parse_esummary_record <- function(x, input_id, keep_xml = FALSE, keep_bioproject = FALSE) {
  if (is.null(x)) {
    uid_in <- .biosample_uid_from_input(input_id)
    acc_in <- .biosample_accession_from_input(input_id)
    rec <- data.frame(
      input_id = as.character(input_id),
      entrez_uid = uid_in,
      biosample_id = acc_in,
      accession_present = !is.na(acc_in),
      parse_status = 'missing_esummary',
      sampledata_present = FALSE,
      sampledata_xml_ok = FALSE,
      n_attributes = 0L,
      bioproject = NA_character_,
      title_raw = NA_character_,
      description_raw = NA_character_,
      organism_raw = NA_character_,
      package_raw = NA_character_,
      stringsAsFactors = FALSE
    )
    if (isTRUE(keep_xml)) rec$sampledata_xml <- NA_character_
    return(list(record = rec, attributes = .biosample_empty_attributes_unique()))
  }
  meta <- .biosample_extract_record_fields(x = x, input_id = input_id)
  sampledata_xml <- meta$sampledata_xml
  if (is.na(sampledata_xml) || !nzchar(sampledata_xml)) {
    rec <- data.frame(
      input_id = as.character(input_id),
      entrez_uid = meta$entrez_uid,
      biosample_id = meta$biosample_id,
      accession_present = !is.na(meta$biosample_id) && nzchar(meta$biosample_id),
      parse_status = 'missing_sampledata',
      sampledata_present = FALSE,
      sampledata_xml_ok = FALSE,
      n_attributes = 0L,
      bioproject = NA_character_,
      title_raw = meta$title_raw,
      description_raw = meta$description_raw,
      organism_raw = meta$organism_raw,
      package_raw = meta$package_raw,
      stringsAsFactors = FALSE
    )
    if (isTRUE(keep_xml)) rec$sampledata_xml <- sampledata_xml
    return(list(record = rec, attributes = .biosample_empty_attributes_unique()))
  }
  doc <- .biosample_read_sampledata_xml(sampledata_xml)
  if (is.null(doc)) {
    rec <- data.frame(
      input_id = as.character(input_id),
      entrez_uid = meta$entrez_uid,
      biosample_id = meta$biosample_id,
      accession_present = !is.na(meta$biosample_id) && nzchar(meta$biosample_id),
      parse_status = 'xml_parse_failed',
      sampledata_present = TRUE,
      sampledata_xml_ok = FALSE,
      n_attributes = 0L,
      bioproject = NA_character_,
      title_raw = meta$title_raw,
      description_raw = meta$description_raw,
      organism_raw = meta$organism_raw,
      package_raw = meta$package_raw,
      stringsAsFactors = FALSE
    )
    if (isTRUE(keep_xml)) rec$sampledata_xml <- sampledata_xml
    return(list(record = rec, attributes = .biosample_empty_attributes_unique()))
  }
  attr_tbl <- .biosample_extract_attribute_rows(doc)
  if (nrow(attr_tbl)) {
    n_attr <- nrow(attr_tbl)
    attr_tbl <- data.frame(
      input_id = rep(as.character(input_id), n_attr),
      entrez_uid = rep(meta$entrez_uid, n_attr),
      biosample_id = rep(meta$biosample_id, n_attr),
      parse_status = rep('parsed', n_attr),
      attr_tbl,
      stringsAsFactors = FALSE
    )
    status <- 'parsed'
  } else {
    attr_tbl <- .biosample_empty_attributes_unique()
    status <- 'parsed_no_extractable_text'
    n_attr <- 0L
  }
  rec <- data.frame(
    input_id = as.character(input_id),
    entrez_uid = meta$entrez_uid,
    biosample_id = meta$biosample_id,
    accession_present = !is.na(meta$biosample_id) && nzchar(meta$biosample_id),
    parse_status = status,
    sampledata_present = TRUE,
    sampledata_xml_ok = TRUE,
    n_attributes = as.integer(n_attr),
    bioproject = if (isTRUE(keep_bioproject)) .biosample_extract_bioproject(doc) else NA_character_,
    title_raw = meta$title_raw,
    description_raw = meta$description_raw,
    organism_raw = meta$organism_raw,
    package_raw = meta$package_raw,
    stringsAsFactors = FALSE
  )
  if (isTRUE(keep_xml)) rec$sampledata_xml <- sampledata_xml
  list(record = rec, attributes = attr_tbl)
}

.biosample_progress_total <- function(results, esummary_cache = NULL, post_steps = 0L) {
  id_map <- .biosample_collect_id_map(results)
  if (!nrow(id_map)) return(0L)
  unique_ids <- unique(id_map$input_id)
  cache <- .biosample_normalise_cache(esummary_cache)
  cache_names <- names(cache) %||% character()
  missing_ids <- setdiff(unique_ids, cache_names)
  as.integer(length(missing_ids) + length(unique_ids) + as.integer(post_steps %||% 0L))
}

.biosample_metadata_core <- function(results, species = NULL, esummary_cache = NULL, batch_size = 100L, progress = interactive(), progress_state = NULL, include_empty_records = TRUE, keep_xml = FALSE, keep_bioproject = FALSE) {
  if (is.null(results) || !is.list(results) || !length(results)) .gama_stop('`.biosample_metadata_core()`: `results` must be a non-empty list from query_species().')
  if (is.null(names(results)) || any(!nzchar(names(results)))) .gama_stop('`.biosample_metadata_core()`: `results` must be a named list (species names expected).')
  if (!is.null(species)) {
    species <- unique(as.character(species))
    keep <- intersect(species, names(results))
    if (!length(keep)) .gama_stop('`.biosample_metadata_core()`: no matching species found in `results`.')
    results <- results[keep]
  }
  id_map <- .biosample_collect_id_map(results)
  if (!nrow(id_map)) {
    out <- .biosample_empty_attributes()
    rec <- .biosample_empty_records(keep_xml = keep_xml)
    diag <- .biosample_make_diagnostics(requested_species = names(results), n_species = length(results), n_species_id_pairs = 0L, n_unique_input_ids = 0L, n_cached = 0L, n_fetched = 0L, parse_status_counts = stats::setNames(integer(), character()))
    attr(out, 'biosample_records') <- rec
    attr(out, 'biosample_diagnostics') <- diag
    attr(out, 'biosample_esummary_cache') <- stats::setNames(vector('list', 0L), character())
    attr(out, 'query_info') <- attr(results, 'query_info')
    return(out)
  }
  close_progress <- is.null(progress_state)
  if (is.null(progress_state)) {
    progress_total <- .biosample_progress_total(results, esummary_cache = esummary_cache)
    progress_state <- .pb_state(progress_total, enabled = isTRUE(progress))
  }
  if (isTRUE(close_progress) && !is.null(progress_state)) on.exit(.pb_close_state(progress_state), add = TRUE)
  unique_ids <- unique(id_map$input_id)
  cache <- .biosample_normalise_cache(esummary_cache)
  cache_names <- names(cache) %||% character()
  cached_ids <- intersect(unique_ids, cache_names)
  missing_ids <- setdiff(unique_ids, cached_ids)
  if (length(missing_ids)) {
    fetched <- .fetch_esummary_batched('biosample', missing_ids, batch_size = as.integer(batch_size), progress_state = progress_state, progress_by = 'record')
    fetched <- .normalise_esummary_list(fetched, missing_ids)
    cache[missing_ids] <- fetched[missing_ids]
  }
  sums_by_input_id <- cache[unique_ids]
  names(sums_by_input_id) <- unique_ids
  record_rows <- vector('list', length(unique_ids))
  attr_rows_list <- vector('list', length(unique_ids))
  for (i in seq_along(unique_ids)) {
    input_id <- unique_ids[[i]]
    x <- sums_by_input_id[[input_id]] %||% NULL
    parsed <- .biosample_parse_esummary_record(x = x, input_id = input_id, keep_xml = keep_xml, keep_bioproject = keep_bioproject)
    record_rows[[i]] <- parsed$record
    attr_rows_list[[i]] <- parsed$attributes
    .pb_advance(progress_state)
  }
  records_unique <- do.call(rbind, record_rows)
  attrs_unique <- do.call(rbind, attr_rows_list)
  records <- merge(id_map, records_unique, by = 'input_id', all.x = TRUE, sort = FALSE)
  records <- records[, c('species','input_id', setdiff(names(records), c('species','input_id'))), drop = FALSE]
  if (nrow(attrs_unique)) {
    attrs <- merge(id_map, attrs_unique, by = 'input_id', all.y = TRUE, sort = FALSE)
    attrs$species <- as.character(attrs$species)
    attrs <- attrs[, c('species','input_id','entrez_uid','biosample_id','parse_status','source_channel','node_name','attribute_name_raw','attribute_name_norm','attribute_name_harmonised','attribute_value_raw','attribute_value_norm','attribute_unit_raw','attribute_index'), drop = FALSE]
  } else {
    attrs <- .biosample_empty_attributes()
  }
  if (isTRUE(include_empty_records)) {
    have <- unique(attrs$input_id)
    missing <- setdiff(records$input_id, have)
    if (length(missing)) {
      pad <- records[records$input_id %in% missing, c('species','input_id','entrez_uid','biosample_id','parse_status'), drop = FALSE]
      pad$source_channel <- NA_character_
      pad$node_name <- NA_character_
      pad$attribute_name_raw <- NA_character_
      pad$attribute_name_norm <- NA_character_
      pad$attribute_name_harmonised <- NA_character_
      pad$attribute_value_raw <- NA_character_
      pad$attribute_value_norm <- NA_character_
      pad$attribute_unit_raw <- NA_character_
      pad$attribute_index <- NA_integer_
      attrs <- rbind(attrs, pad)
    }
  }
  parse_counts <- table(records_unique$parse_status %||% character())
  diag <- .biosample_make_diagnostics(requested_species = names(results), n_species = length(results), n_species_id_pairs = nrow(id_map), n_unique_input_ids = length(unique_ids), n_cached = length(cached_ids), n_fetched = length(missing_ids), parse_status_counts = parse_counts)
  attr(attrs, 'biosample_records') <- records
  attr(attrs, 'biosample_diagnostics') <- diag
  attr(attrs, 'biosample_esummary_cache') <- sums_by_input_id
  attr(attrs, 'query_info') <- attr(results, 'query_info')
  attrs
}

.biosample_empty_tissue_classification <- function() {
  tibble::tibble(
    species = character(),
    input_id = character(),
    entrez_uid = character(),
    biosample_id = character(),
    tissue_raw = character(),
    tissue_norm = character(),
    anatomy_term = character(),
    anatomy_class = character(),
    anatomy_subclass = character(),
    rank = character(),
    ontology_namespace = character(),
    ontology_id = character(),
    ontology_label = character()
  )
}

.biosample_resolve_key <- function(biosample_id, input_id) {
  key <- as.character(biosample_id)
  src <- as.character(input_id)
  miss <- is.na(key) | !nzchar(key)
  key[miss] <- src[miss]
  key
}

.biosample_classify_tissue_attributes <- function(meta, species = NULL) {
  if (is.null(meta) || !is.data.frame(meta)) .gama_stop('`.biosample_classify_tissue_attributes()`: `meta` must be a data.frame.')
  req <- c('species', 'input_id', 'entrez_uid', 'biosample_id', 'parse_status', 'attribute_name_norm', 'attribute_name_harmonised', 'attribute_value_raw', 'attribute_value_norm')
  miss <- setdiff(req, names(meta))
  if (length(miss)) .gama_stop('`.biosample_classify_tissue_attributes()`: `meta` missing columns: ', paste(miss, collapse = ', '))
  tissue <- meta |>
    dplyr::filter(
      .data$parse_status == 'parsed',
      .biosample_is_tissue_attribute(.data$attribute_name_norm, .data$attribute_name_harmonised)
    ) |>
    dplyr::transmute(
      species = as.character(.data$species),
      input_id = as.character(.data$input_id),
      entrez_uid = as.character(.data$entrez_uid),
      biosample_id = .biosample_resolve_key(.data$biosample_id, .data$input_id),
      tissue_raw = as.character(.data$attribute_value_raw),
      tissue_norm = as.character(.data$attribute_value_norm)
    ) |>
    dplyr::filter(!is.na(.data$species), nzchar(.data$species), !is.na(.data$biosample_id), nzchar(.data$biosample_id))
  if (!is.null(species)) tissue <- tissue |> dplyr::filter(.data$species %in% .env$species)
  if (!nrow(tissue)) return(.biosample_empty_tissue_classification())
  ref <- .biosample_anatomy_ref()
  lex <- tibble::as_tibble(ref$lex)
  vals <- unique(as.character(tissue$tissue_norm))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  lookup <- if (length(vals)) {
    lapply(vals, function(val) {
      terms <- if (.biosample_is_missing_like_value(val)) character() else unique(.match_biosample_anatomy_value(val, ref$variant_map))
      if (!length(terms)) terms <- 'unknown'
      tibble::tibble(
        tissue_norm = val,
        anatomy_term = as.character(terms)
      )
    }) |>
      dplyr::bind_rows()
  } else {
    tibble::tibble(
      tissue_norm = character(),
      anatomy_term = character()
    )
  }
  tissue |>
    dplyr::left_join(lookup, by = 'tissue_norm', relationship = 'many-to-many') |>
    dplyr::mutate(
      anatomy_term = dplyr::if_else(
        .biosample_is_missing_like_value(.data$tissue_norm) | is.na(.data$anatomy_term) | !nzchar(dplyr::coalesce(.data$anatomy_term, '')),
        'unknown',
        .data$anatomy_term
      )
    ) |>
    dplyr::left_join(lex, by = 'anatomy_term') |>
    dplyr::mutate(
      anatomy_class = dplyr::if_else(.data$anatomy_term == 'unknown', 'unknown', as.character(.data$anatomy_class)),
      anatomy_subclass = dplyr::if_else(.data$anatomy_term == 'unknown', 'unknown', as.character(.data$anatomy_subclass)),
      rank = dplyr::if_else(.data$anatomy_term == 'unknown', 'meta', as.character(.data$rank)),
      ontology_namespace = dplyr::if_else(.data$anatomy_term == 'unknown', NA_character_, as.character(.data$ontology_namespace)),
      ontology_id = dplyr::if_else(.data$anatomy_term == 'unknown', NA_character_, as.character(.data$ontology_id)),
      ontology_label = dplyr::if_else(.data$anatomy_term == 'unknown', NA_character_, as.character(.data$ontology_label))
    ) |>
    dplyr::transmute(
      species = .data$species,
      input_id = .data$input_id,
      entrez_uid = .data$entrez_uid,
      biosample_id = .data$biosample_id,
      tissue_raw = .data$tissue_raw,
      tissue_norm = .data$tissue_norm,
      anatomy_term = .data$anatomy_term,
      anatomy_class = .data$anatomy_class,
      anatomy_subclass = .data$anatomy_subclass,
      rank = .data$rank,
      ontology_namespace = .data$ontology_namespace,
      ontology_id = .data$ontology_id,
      ontology_label = .data$ontology_label
    ) |>
    dplyr::distinct()
}

.biosample_anatomy_profile <- function(tissue_classification, species = NULL) {
  anatomy_levels <- .biosample_anatomy_profile_levels()
  subclass_levels <- c(.biosample_anatomy_subclass_levels(), 'mixed', 'unknown')
  out <- tibble::tibble(
    species = character(),
    biosample_id = character(),
    anatomy_class = character(),
    anatomy_subclass = character()
  )
  if (is.null(tissue_classification) || !is.data.frame(tissue_classification) || !nrow(tissue_classification)) return(out)
  req <- c('species', 'biosample_id', 'anatomy_class', 'anatomy_subclass')
  miss <- setdiff(req, names(tissue_classification))
  if (length(miss)) .gama_stop('`.biosample_anatomy_profile()`: `tissue_classification` missing columns: ', paste(miss, collapse = ', '))
  profiled <- tissue_classification |>
    dplyr::transmute(
      species = as.character(.data$species),
      biosample_id = as.character(.data$biosample_id),
      anatomy_class = as.character(.data$anatomy_class),
      anatomy_subclass = as.character(.data$anatomy_subclass)
    ) |>
    dplyr::filter(!is.na(.data$species), nzchar(.data$species), !is.na(.data$biosample_id), nzchar(.data$biosample_id))
  if (!is.null(species)) profiled <- profiled |> dplyr::filter(.data$species %in% .env$species)
  if (!nrow(profiled)) return(out)
  profiled |>
    dplyr::group_by(.data$species, .data$biosample_id) |>
    dplyr::summarise(
      anatomy_class = .biosample_collapse_anatomy_profile(.data$anatomy_class, level = 'anatomy_class'),
      anatomy_subclass = .biosample_collapse_anatomy_profile(.data$anatomy_subclass, level = 'anatomy_subclass'),
      .groups = 'drop'
    ) |>
    dplyr::filter(.data$anatomy_class %in% .env$anatomy_levels) |>
    dplyr::mutate(
      anatomy_subclass = dplyr::if_else(
        .data$anatomy_subclass %in% .env$subclass_levels,
        .data$anatomy_subclass,
        'unknown'
      )
    ) |>
    dplyr::distinct(.data$species, .data$biosample_id, .data$anatomy_class, .data$anatomy_subclass) |>
    dplyr::arrange(.data$species, .data$biosample_id, .data$anatomy_class, .data$anatomy_subclass)
}

.biosample_canonical_profile <- function(tissue_classification, species = NULL) {
  out <- tibble::tibble(
    species = character(),
    biosample_id = character(),
    anatomy_term = character(),
    anatomy_class = character(),
    anatomy_subclass = character(),
    rank = character(),
    ontology_namespace = character(),
    ontology_id = character(),
    ontology_label = character()
  )
  if (is.null(tissue_classification) || !is.data.frame(tissue_classification) || !nrow(tissue_classification)) return(out)
  req <- c('species', 'biosample_id', 'anatomy_term', 'anatomy_class', 'anatomy_subclass', 'rank', 'ontology_namespace', 'ontology_id', 'ontology_label')
  miss <- setdiff(req, names(tissue_classification))
  if (length(miss)) .gama_stop('`.biosample_canonical_profile()`: `tissue_classification` missing columns: ', paste(miss, collapse = ', '))
  out <- tissue_classification |>
    dplyr::transmute(
      species = as.character(.data$species),
      biosample_id = as.character(.data$biosample_id),
      anatomy_term = as.character(.data$anatomy_term),
      anatomy_class = as.character(.data$anatomy_class),
      anatomy_subclass = as.character(.data$anatomy_subclass),
      rank = as.character(.data$rank),
      ontology_namespace = as.character(.data$ontology_namespace),
      ontology_id = as.character(.data$ontology_id),
      ontology_label = as.character(.data$ontology_label)
    ) |>
    dplyr::filter(!is.na(.data$species), nzchar(.data$species), !is.na(.data$biosample_id), nzchar(.data$biosample_id), !is.na(.data$anatomy_term), nzchar(.data$anatomy_term))
  if (!is.null(species)) out <- out |> dplyr::filter(.data$species %in% .env$species)
  out |>
    dplyr::distinct(.data$species, .data$biosample_id, .data$anatomy_term, .keep_all = TRUE) |>
    dplyr::arrange(.data$species, .data$anatomy_class, .data$anatomy_term)
}

.biosample_modality_profile <- function(sra_profile) {
  modality_levels <- .biosample_modality_levels()
  out <- tibble::tibble(
    species = character(),
    biosample_id = character(),
    modality_class = character()
  )
  if (is.null(sra_profile) || !is.data.frame(sra_profile) || !nrow(sra_profile)) return(out)
  if (all(c('species', 'biosample_id', 'modality_class') %in% names(sra_profile))) {
    raw <- tibble::as_tibble(sra_profile) |>
      dplyr::transmute(
        species = as.character(.data$species),
        biosample_id = as.character(.data$biosample_id),
        modality_raw = as.character(.data$modality_class)
      )
  } else {
    req <- c('species', 'class')
    biosample_col <- if ('biosample_id' %in% names(sra_profile)) 'biosample_id' else if ('biosample' %in% names(sra_profile)) 'biosample' else NA_character_
    if (is.na(biosample_col) || any(!req %in% names(sra_profile))) {
      .gama_stop('`.biosample_modality_profile()`: `sra_profile` must contain `species`, a BioSample ID column (`biosample` or `biosample_id`), and `class`, or already-collapsed `modality_class`.')
    }
    raw <- tibble::as_tibble(sra_profile) |>
      dplyr::transmute(
        species = as.character(.data$species),
        biosample_id = as.character(.data[[biosample_col]]),
        modality_raw = as.character(.data$class)
      )
  }
  raw |>
    dplyr::mutate(
      species = trimws(.data$species),
      biosample_id = trimws(.data$biosample_id),
      modality_raw = trimws(.data$modality_raw)
    ) |>
    dplyr::filter(!is.na(.data$species), nzchar(.data$species), !is.na(.data$biosample_id), nzchar(.data$biosample_id)) |>
    dplyr::group_by(.data$species, .data$biosample_id) |>
    dplyr::summarise(
      modality_class = {
        vals <- unique(.data$modality_raw[!is.na(.data$modality_raw) & nzchar(.data$modality_raw)])
        vals <- vals[vals %in% modality_levels]
        if (!length(vals)) 'unknown' else if (length(vals) == 1L) vals else 'mixed'
      },
      .groups = 'drop'
    ) |>
    dplyr::arrange(.data$species, .data$biosample_id)
}

.biosample_term_heatmap_core <- function(BIO, SRA, species = NULL) {
  class_levels <- .biosample_anatomy_profile_levels()
  subclass_core <- .biosample_anatomy_subclass_map() |>
    dplyr::filter(
      !.data$anatomy_class %in% c('mixed', 'unknown'),
      !.data$anatomy_subclass %in% c('mixed', 'unknown')
    ) |>
    dplyr::distinct(.data$anatomy_subclass, .keep_all = TRUE)
  subclass_meta <- dplyr::bind_rows(
    subclass_core,
    tibble::tibble(
      anatomy_class = c('mixed', 'unknown'),
      anatomy_subclass = c('mixed', 'unknown')
    )
  ) |>
    dplyr::mutate(term_order = dplyr::row_number())
  if (anyDuplicated(subclass_meta$anatomy_subclass)) {
    .gama_stop('Duplicated BioSample anatomy subclass labels detected in interaction metadata.')
  }
  subclass_levels <- subclass_meta$anatomy_subclass
  modality_levels <- .biosample_modality_levels()
  BIO <- .gama_require_output(
    BIO,
    'summarise_biosample_availability',
    required_cols = c('species')
  )
  SRA <- .gama_require_output(SRA, 'summarise_sra_availability')
  canonical <- .gama_require_cache(
    BIO,
    attr_name = 'biosample_canonical_profile',
    required_cols = c('species', 'biosample_id', 'anatomy_subclass'),
    source = 'summarise_biosample_availability'
  )
  species_all <- unique(as.character(BIO$species))
  species_all <- species_all[!is.na(species_all) & nzchar(species_all)]
  species_use <- if (is.null(species)) {
    species_all
  } else {
    sp_in <- unique(as.character(species))
    sp_in <- sp_in[!is.na(sp_in) & nzchar(sp_in)]
    missing <- setdiff(sp_in, species_all)
    if (length(missing)) .gama_warn('Requested species not found in `BIO`: ', paste(missing, collapse = ', '), '. Dropping.')
    sp_in[sp_in %in% species_all]
  }
  if (!length(species_use)) {
    return(tibble::tibble(
      species = character(),
      modality_class = factor(levels = modality_levels),
      anatomy_class = factor(levels = class_levels),
      anatomy_subclass = character(),
      term_order = integer(),
      BioSample = integer()
    ))
  }
  sra_profile <- .gama_require_cache(
    SRA,
    attr_name = 'sra_profile',
    required_cols = c('species', 'class'),
    source = 'summarise_sra_availability'
  )
  if (!any(c('biosample', 'biosample_id') %in% names(sra_profile))) {
    .gama_input_error(
      'summarise_sra_availability',
      detected = .detect_gama_object(SRA),
      detail = "cache 'sra_profile' missing required BioSample ID column: biosample or biosample_id."
    )
  }
  modality <- .biosample_modality_profile(sra_profile)
  counts <- canonical |>
    dplyr::filter(.data$species %in% .env$species_use) |>
    dplyr::left_join(modality, by = c('species', 'biosample_id')) |>
    dplyr::filter(!is.na(.data$modality_class), nzchar(.data$modality_class)) |>
    dplyr::transmute(
      species = .data$species,
      biosample_id = .data$biosample_id,
      modality_class = .data$modality_class,
      anatomy_subclass = dplyr::coalesce(.data$anatomy_subclass, 'unknown')
    ) |>
    dplyr::mutate(
      anatomy_subclass = dplyr::if_else(
        .data$anatomy_subclass %in% .env$subclass_levels,
        .data$anatomy_subclass,
        'unknown'
      )
    ) |>
    dplyr::distinct(.data$species, .data$biosample_id, .data$modality_class, .data$anatomy_subclass) |>
    dplyr::group_by(.data$species, .data$biosample_id, .data$modality_class) |>
    dplyr::summarise(
      anatomy_subclass = .biosample_collapse_anatomy_profile(.data$anatomy_subclass, level = 'anatomy_subclass'),
      .groups = 'drop'
    ) |>
    dplyr::left_join(
      subclass_meta |>
        dplyr::select(anatomy_subclass, anatomy_class),
      by = 'anatomy_subclass'
    ) |>
    dplyr::mutate(anatomy_class = dplyr::coalesce(.data$anatomy_class, 'unknown')) |>
    dplyr::count(.data$species, .data$modality_class, .data$anatomy_class, .data$anatomy_subclass, name = 'BioSample')
  out_list <- vector('list', length(species_use))
  for (i in seq_along(species_use)) {
    sp <- species_use[[i]]
    counts_sp <- counts |>
      dplyr::filter(.data$species == .env$sp)
    grid_sp <- tidyr::expand_grid(
      modality_class = modality_levels,
      subclass_meta
    ) |>
      dplyr::mutate(species = sp) |>
      dplyr::left_join(
        counts_sp,
        by = c('species', 'modality_class', 'anatomy_class', 'anatomy_subclass')
      ) |>
      dplyr::transmute(
        species = .data$species,
        modality_class = factor(.data$modality_class, levels = modality_levels),
        anatomy_class = factor(as.character(.data$anatomy_class), levels = class_levels),
        anatomy_subclass = .data$anatomy_subclass,
        term_order = as.integer(.data$term_order),
        BioSample = as.integer(dplyr::coalesce(.data$BioSample, 0L))
      )
    out_list[[i]] <- grid_sp
  }
  dplyr::bind_rows(out_list)
}

.biosample_qc <- function(meta) {
  rec <- attr(meta, 'biosample_records')
  if (is.null(rec) || !is.data.frame(rec)) .gama_stop('`.biosample_qc()`: meta missing attr(meta, \'biosample_records\').')
  req <- c('species','input_id','parse_status','n_attributes','sampledata_present','sampledata_xml_ok')
  miss <- setdiff(req, names(rec))
  if (length(miss)) .gama_stop('`.biosample_qc()`: biosample_records missing columns: ', paste(miss, collapse = ', '))
  rec <- rec[!is.na(rec$input_id) & nzchar(rec$input_id), , drop = FALSE]
  spp <- unique(rec$species)
  out <- data.frame(species = spp, stringsAsFactors = FALSE)
  out$biosample_total <- as.integer(tapply(rec$input_id, rec$species, function(x) length(unique(x)))[spp])
  out$parsed_n <- as.integer(tapply(rec$parse_status == 'parsed', rec$species, sum, na.rm = TRUE)[spp])
  out$parsed_no_extractable_text_n <- as.integer(tapply(rec$parse_status == 'parsed_no_extractable_text', rec$species, sum, na.rm = TRUE)[spp])
  out$missing_esummary_n <- as.integer(tapply(rec$parse_status == 'missing_esummary', rec$species, sum, na.rm = TRUE)[spp])
  out$missing_sampledata_n <- as.integer(tapply(rec$parse_status == 'missing_sampledata', rec$species, sum, na.rm = TRUE)[spp])
  out$xml_parse_failed_n <- as.integer(tapply(rec$parse_status == 'xml_parse_failed', rec$species, sum, na.rm = TRUE)[spp])
  out$error_n <- out$missing_esummary_n + out$missing_sampledata_n + out$xml_parse_failed_n
  out$parsed_prop <- ifelse(out$biosample_total > 0, out$parsed_n / out$biosample_total, NA_real_)
  out$parsed_no_extractable_text_prop <- ifelse(out$biosample_total > 0, out$parsed_no_extractable_text_n / out$biosample_total, NA_real_)
  out$error_prop <- ifelse(out$biosample_total > 0, out$error_n / out$biosample_total, NA_real_)
  problem_ids <- rec[rec$parse_status != 'parsed', c('species','input_id','biosample_id','parse_status','sampledata_present','sampledata_xml_ok','n_attributes'), drop = FALSE]
  problem_ids <- problem_ids[order(problem_ids$species, problem_ids$parse_status, problem_ids$input_id), , drop = FALSE]
  list(qc = out, problem_ids = problem_ids)
}

biosample_qc <- function(meta) {
  .biosample_qc(meta)
}

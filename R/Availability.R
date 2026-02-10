# AVAILABILITY ================================================================
.query_assembly  <- function(sp) .ncbi_search('assembly',  sp)
.query_sra       <- function(sp) .ncbi_search('sra',       sp)
.query_biosample <- function(sp) .ncbi_search('biosample', sp)
#' Query NCBI databases for user-defined species sets using rentrez
#'
#' @param species Character vector of binomial species names.
#' @return A named list (one element per species) containing database search results, record IDs, accession counts, and query metadata.
#' @export
query_species <- function(species) {
  species <- unique(species)
  n <- length(species)
  pb <- .pb_start(n)
  RESULTS <- lapply(seq_along(species), function(i) {
    sp <- species[i]
    RES <- list(
      assembly  = .query_assembly(sp),
      sra       = .query_sra(sp),
      biosample = .query_biosample(sp)
    )
    .pb_tick(pb, i)
    RES
  })
  .pb_end(pb)
  names(RESULTS) <- species
  attr(RESULTS, 'query_info') <- .make_query_info(
    species = species,
    dbs     = c('assembly', 'sra', 'biosample'),
    retmax  = 999999
  )
  RESULTS
}
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
#' Summarise accession counts and compute data richness scores
#'
#' @param results List returned by [query_species()].
#' @return A tibble of per-species accession counts and composite data richness scores (class 'gdt_tbl').
#' @export
summarise_availability <- function(results) {
  SPECIES <- names(results)
  COMPONENTS <- lapply(results, function(x) .score_species(x$assembly, x$sra, x$biosample))
  OUT <- tibble::tibble(
    species   = SPECIES,
    Assembly  = purrr::map_int(results, ~ .x$assembly$count %||% 0L),
    SRA       = purrr::map_int(results, ~ .x$sra$count %||% 0L),
    BioSample = purrr::map_int(results, ~ .x$biosample$count %||% 0L),
    A         = purrr::map_dbl(COMPONENTS, 'richness'),
    S         = purrr::map_dbl(COMPONENTS, 'sra_score'),
    B         = purrr::map_dbl(COMPONENTS, 'biosample_score'),
    score     = purrr::map_dbl(COMPONENTS, 'score')
  )
  .as_gdt_table(OUT, results)
}

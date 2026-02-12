# AVAILABILITY ANALYSIS =======================================================

#' Query NCBI databases using a set of species
#'
#' @param species Character vector of binomial species names.
#' @return A named list (one element per species) containing database search results, record IDs, accession counts, and query metadata.
#' @export
query_species <- function(species) {
  species <- unique(species)
  n <- length(species)
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  RESULTS <- lapply(seq_along(species), function(i) {
    sp <- species[i]
    RES <- list(
    assembly  = .ncbi_search('assembly',  sp),
    sra       = .ncbi_search('sra',       sp),
    biosample = .ncbi_search('biosample', sp)
    )
    utils::setTxtProgressBar(pb, i)
    RES
  })
  close(pb)
  names(RESULTS) <- species
  attr(RESULTS, 'query_info') <- .make_query_info(
  species = species,
  dbs     = c('assembly', 'sra', 'biosample'),
  retmax  = 999999
  )
  RESULTS
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

#' Summarise SRA experimental modality composition
#'
#' @param results List returned by [query_species()].
#' @param species Optional character vector of species to include.
#' @param all Logical; if TRUE, include subclass-level counts.
#' @param include_geo Logical; if TRUE, append GEO linkage summaries.
#' @return A tibble of per-species SRA modality counts and totals (class 'gdt_tbl'); optionally includes subclass and GEO overlay columns.
#' @export
summarise_sra_availability <- function(results,
species     = NULL,
all         = FALSE,
include_geo = FALSE) {
  META <- .sra_metadata_core(results, species = species)
  CLASS <- META |>
  dplyr::count(species, class, name = 'count') |>
  tidyr::pivot_wider(
  names_from  = class,
  values_from = count,
  values_fill = 0
  )
  for (m in c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')) {
    if (!m %in% names(CLASS)) CLASS[[m]] <- 0L
  }
  SRA_TOTAL <- META |>
  dplyr::count(species, name = 'SRA')
  OUT <- CLASS |>
  dplyr::left_join(SRA_TOTAL, by = 'species') |>
  dplyr::select(species, SRA, genomic, transcriptomic, epigenomic, chromatin, other, unknown)
  if (all) {
    SUB <- META |>
    dplyr::count(species, subclass, name = 'count') |>
    tidyr::pivot_wider(
    names_from  = subclass,
    values_from = count,
    values_fill = 0
    )
    OUT <- dplyr::left_join(OUT, SUB, by = 'species')
  }
  if (!include_geo) return(.as_gdt_table(OUT, results))
  denom <- c('transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
  DENOM <- META |>
  dplyr::filter(.data$class %in% !!denom) |>
  dplyr::count(species, name = 'denom_total')
  GEO_DENOM <- META |>
  dplyr::filter(.data$class %in% !!denom) |>
  dplyr::filter(.data$geo_linked) |>
  dplyr::count(species, name = 'geo_linked_denom')
  GEO_CLASS <- META |>
  dplyr::filter(.data$geo_linked) |>
  dplyr::count(species, class, name = 'count') |>
  tidyr::pivot_wider(
  names_from  = class,
  values_from = count,
  values_fill = 0
  )
  for (m in c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')) {
    if (!m %in% names(GEO_CLASS)) GEO_CLASS[[m]] <- 0L
  }
  GEO_CLASS <- GEO_CLASS |>
  dplyr::transmute(
  species,
  genomic_geo        = as.integer(.data$genomic),
  transcriptomic_geo = as.integer(.data$transcriptomic),
  epigenomic_geo     = as.integer(.data$epigenomic),
  chromatin_geo      = as.integer(.data$chromatin),
  other_geo          = as.integer(.data$other),
  unknown_geo        = as.integer(.data$unknown)
  )
  OUT2 <- OUT |>
  dplyr::left_join(DENOM, by = 'species') |>
  dplyr::left_join(GEO_DENOM, by = 'species') |>
  dplyr::mutate(
  denom_total      = as.integer(.data$denom_total %||% 0L),
  geo_linked_denom = as.integer(.data$geo_linked_denom %||% 0L),
  geo_prop         = dplyr::if_else(.data$denom_total > 0,
  .data$geo_linked_denom / .data$denom_total,
  NA_real_)
  ) |>
  dplyr::left_join(GEO_CLASS, by = 'species')
  .as_gdt_table(OUT2, results)
}

#' Extract Assembly accession metadata
#'
#' @param results List returned by [query_species()].
#' @param species Optional character vector of species to include.
#' @param best Logical; if TRUE, return only the best assembly per species.
#' @return A tibble of assembly metadata (e.g., accession, assembly level, N50, coverage, BioSample/BioProject and release information).
#' @export
extract_assembly_metadata <- function(results, species = NULL, best = FALSE) {
  if (!is.null(species)) results <- results[species]
  n  <- length(results)
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  META <- lapply(seq_along(results), function(i) {
    sp  <- names(results)[i]
    res <- results[[i]]
    utils::setTxtProgressBar(pb, i)
    asm <- res$assembly
    if (is.null(asm) || (asm$count %||% 0L) == 0L) {
      return(tibble::tibble(
      species       = sp,
      accession     = NA_character_,
      level         = NA_character_,
      n50           = NA_real_,
      coverage      = NA_real_,
      biosample     = NA_character_,
      bioproject    = NA_character_,
      submitter     = NA_character_,
      release_date  = NA_character_,
      ftp_path      = NA_character_
      ))
    }
    IDS  <- asm$ids %||% character()
    SUMS <- .fetch_esummary_batched('assembly', IDS)
    SUMS <- .normalise_esummary_list(SUMS, IDS)
    do.call(dplyr::bind_rows, lapply(names(SUMS), function(acc) {
      x <- SUMS[[acc]]
      tibble::tibble(
      species       = sp,
      accession     = acc,
      level         = .extract_assembly_level(x),
      n50           = .extract_n50(x),
      coverage      = as.numeric(.flatten_to_char(x$coverage)),
      biosample     = .flatten_to_char(x$biosampleaccn),
      bioproject    = .flatten_to_char(x$gb_bioprojects$bioprojectaccn),
      submitter     = .flatten_to_char(x$submitterorganization),
      release_date  = .flatten_to_char(x$asmreleasedate_genbank),
      ftp_path      = .flatten_to_char(x$ftppath_genbank)
      )
    }))
  })
  close(pb)
  META <- dplyr::bind_rows(META)
  if (!best) return(.as_gdt_table(META, results))
  BEST <- lapply(seq_along(results), function(i) {
    sp  <- names(results)[i]
    res <- results[[i]]
    asm <- res$assembly
    if (is.null(asm) || (asm$count %||% 0L) == 0L) {
      return(tibble::tibble(species = sp, accession = NA_character_))
    }
    IDS  <- asm$ids %||% character()
    SUMS <- .fetch_esummary_batched('assembly', IDS)
    SUMS <- .normalise_esummary_list(SUMS, IDS)
    LEVELS <- sapply(SUMS, .extract_assembly_level)
    STRUCT <- .ASSEMBLY_WEIGHTS[LEVELS]
    STRUCT[is.na(STRUCT)] <- 0
    N50 <- sapply(SUMS, .extract_n50)
    max_struct <- max(STRUCT)
    tied_idx   <- which(STRUCT == max_struct)
    best_idx <- if (length(tied_idx) > 1) {
      tied_idx[which.max(N50[tied_idx])]
    } else tied_idx
    tibble::tibble(
    species   = sp,
    accession = IDS[best_idx]
    )
  })
  BEST <- dplyr::bind_rows(BEST)
  OUT  <- dplyr::inner_join(META, BEST, by = c('species', 'accession'))
  .as_gdt_table(OUT, results)
}

#' Extract SRA accession metadata and classify experimental modality
#'
#' @param results List returned by [query_species()].
#' @param species Optional character vector of species to include.
#' @param class Optional character vector of ontology classes to retain.
#' @param subclass Optional character vector of ontology subclasses to retain.
#' @param only_geo Logical; if TRUE, retain only GEO-linked experiments.
#' @return A tibble of experiment-level SRA metadata with strategy normalisation, ontology assignments and GEO linkage fields.
#' @export
extract_sra_metadata <- function(results,
species  = NULL,
class    = NULL,
subclass = NULL,
only_geo = FALSE) {
  META <- .sra_metadata_core(results, species = species)
  if (!is.null(class)) {
    META <- META |> dplyr::filter(.data$class %in% !!class)
  }
  if (!is.null(subclass)) {
    META <- META |> dplyr::filter(.data$subclass %in% !!subclass)
  }
  if (isTRUE(only_geo)) {
    META <- META |> dplyr::filter(.data$geo_linked)
  }
  .as_gdt_table(META, results)
}

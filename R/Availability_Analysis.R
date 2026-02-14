# AVAILABILITY ANALYSIS =======================================================

#' Query NCBI databases using species list
#'
#' Queries NCBI Assembly, SRA, and BioSample and returns per-species search
#' results in a named list. A provenance record is attached to the output as
#' the `query_info` attribute.
#'
#' @details
#' `query_species()` is the entry point for the NCBI query phase. Each species
#' is queried independently across the supported databases, and results are
#' stored in a per-species list with components `assembly`, `sra`, and
#' `biosample`.
#'
#' The returned object has an attribute `query_info` containing the tool
#' version, query timestamp (UTC), database names, and the search terms used.
#'
#' @param species Character vector of binomial species names (e.g.
#'   'Vigna angularis'). Duplicates are removed with [unique()].
#'
#' @return A named list with one element per species. Each element is a list
#'   with components `assembly`, `sra`, and `biosample` containing
#'   database-specific search results (including counts and record identifiers,
#'   depending on the internal search implementation). The output has a
#'   `query_info` attribute storing query provenance.
#'
#' @seealso [summarise_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' str(RESULTS, max.level = 2)
#' attr(RESULTS, 'query_info')
#' }
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

#' Summarise accession counts and compute data richness
#'
#' Collapses the output of [query_species()] into a species-level summary table
#' of accession counts (Assembly, SRA, BioSample) and composite data richness
#' scores.
#'
#' @details
#' For each species, accession counts are extracted from the underlying query
#' results and combined with a composite scoring system computed by the
#' internal scoring engine. The returned tibble is given class `gdt_tbl` and
#' retains query provenance via the `query_info` attribute, which is carried
#' over from the `results` object.
#'
#' @param results A list returned by [query_species()].
#'
#' @return A tibble with one row per species and the following columns:
#' \itemize{
#'   \item `species`: species name
#'   \item `Assembly`, `SRA`, `BioSample`: accession counts per database
#'   \item `A`, `S`, `B`: component scores used for the composite
#'   \item `score`: composite data richness score
#' }
#' The tibble has class `gdt_tbl` and carries a `query_info` attribute for
#' provenance.
#'
#' @seealso [query_species()], [plot_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SUMMARY <- summarise_availability(RESULTS)
#' print(SUMMARY)
#' attr(SUMMARY, 'query_info')
#' }
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

#' Summarise SRA modality composition
#'
#' Collapses experiment-level SRA metadata into species-level modality counts
#' using ontology-assigned classes and (optionally) subclasses.
#'
#' By default, the output includes class-level counts for the major modality
#' classes (genomic, transcriptomic, epigenomic, chromatin, other, unknown),
#' plus the total number of SRA experiments per species.
#'
#' GEO overlay:
#' When `include_geo = TRUE`, GEO-linked class counts are appended as
#' '<class>_geo' columns, alongside a GEO-compatible denominator summary and
#' proportion. The GEO-compatible denominator is defined as:
#' transcriptomic + epigenomic + chromatin + other + unknown.
#'
#' @param results A list returned by [query_species()].
#' @param species `NULL` (default) to include all species, or a character
#'   vector specifying which species to include.
#' @param all Logical; if `TRUE`, include subclass-level columns.
#' @param include_geo Logical; if `TRUE`, append GEO linkage summaries.
#'
#' @return A tibble with one row per species containing class-level counts and
#'   totals. When `all = TRUE`, subclass-level counts are included. When
#'   `include_geo = TRUE`, GEO summary columns are appended (e.g.
#'   `denom_total`, `geo_linked_denom`, `geo_prop`, and '<class>_geo').
#'   The tibble has class `gdt_tbl` and carries a `query_info` attribute.
#'
#' @seealso [extract_sra_metadata()], [plot_sra_availability()],
#'   [plot_sra_geo_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' print(SRA_SUMMARY)
#' }
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

#' Extract filtered Assembly metadata
#'
#' Retrieves and structures assembly metadata for one or more species using the
#' Assembly identifiers stored in the output of [query_species()]. Metadata are
#' returned in a tidy tibble and optionally reduced to a single 'best' assembly
#' per species.
#'
#' When `best = TRUE`, the best assembly is selected using
#' structural-weighting (assembly level) and N50 as a tie-breaker.
#'
#' @param results A list returned by [query_species()], containing Assembly IDs.
#' @param species `NULL` (default) to return assemblies for all species, or a
#'   character vector specifying which species to extract.
#' @param best Logical; if `TRUE`, return only the best assembly per species.
#'
#' @return A tibble with one row per assembly (or one row per species when
#'   `best = TRUE`). Fields include species, accession, assembly level, N50,
#'   coverage, BioSample/BioProject accessions, submitter, release date, and
#'   FTP path (where available). The tibble has class `gdt_tbl` and carries a
#'   `query_info` attribute for provenance.
#'
#' @seealso [query_species()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' ASM <- extract_assembly_metadata(RESULTS, best = TRUE)
#' print(ASM)
#' attr(ASM, 'query_info')
#' }
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

#' Extract filtered SRA metadata
#'
#' Retrieves experiment-level SRA metadata for one or more species using the
#' internal SRA metadata engine, then normalises sequencing strategy labels and
#' assigns curated modality classes and subclasses.
#'
#' For each experiment, the function:
#' \itemize{
#'   \item extracts the raw `LIBRARY_STRATEGY` value from the SRA XML
#'   \item normalises strategy strings
#'   \item assigns ontology-based `class` and `subclass`
#'   \item records GEO linkage fields (`geo_linked`, `gse_ids`, `gsm_ids`)
#' }
#'
#' Optional filters can restrict output to particular classes/subclasses, or
#' GEO-linked experiments only.
#'
#' @param results A list returned by [query_species()], containing SRA IDs.
#' @param species `NULL` (default) to include all species, or a character
#'   vector specifying which species to include.
#' @param class Optional character vector of modality classes to retain.
#' @param subclass Optional character vector of modality subclasses to retain.
#' @param only_geo Logical; if `TRUE`, retain only GEO-linked experiments.
#'
#' @return A tibble with one row per SRA experiment, including identifiers,
#'   strategy fields, ontology assignments, and GEO linkage columns. The tibble
#'   has class `gdt_tbl` and carries a `query_info` attribute for provenance.
#'
#' @seealso [query_species()], [summarise_sra_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA <- extract_sra_metadata(
#'   RESULTS,
#'   species = 'Vigna vexillata',
#'   class = 'genomic'
#' )
#' print(SRA)
#' attr(SRA, 'query_info')
#' }
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

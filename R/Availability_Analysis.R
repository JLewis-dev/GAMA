# AVAILABILITY ANALYSIS =======================================================

#' Query NCBI databases using a list of species names
#'
#' Queries NCBI Assembly, SRA, and BioSample and returns per-species search
#' results in a named list. A provenance record is attached to the output as
#' the `query_info` attribute.
#'
#' @details
#' `query_species()` is the entry point for the NCBI search phase. Each species
#' is queried independently across the supported databases, and results are
#' stored in a per-species list with components `assembly`, `sra`, and
#' `biosample`. When `synonyms` is supplied, results are collapsed under the
#' canonical species names using unique database record identifiers.
#'
#' The returned object has an attribute `query_info` containing the tool
#' version, timestamp (UTC), databases queried, search terms, and any
#' synonym groups used for canonical collapse.
#'
#' @param species Character vector of binomial species names (e.g.
#' `Vigna angularis`). Duplicates are removed with [unique()].
#' @param synonyms `NULL` (default) for no synonym collapse, or a named list
#' or named character vector mapping canonical species names to one or more
#' synonymous names. Query results for each canonical species are merged across
#' all supplied names using unique database record identifiers.
#'
#' @return A named list with one element per canonical species. Each element is
#' a list with components `assembly`, `sra`, and `biosample` containing
#' database-specific search results (including counts and record identifiers,
#' depending on the internal search implementation). When `synonyms` is used,
#' aliases are merged under the canonical species name without double counting
#' repeated record identifiers. The output has a `query_info` attribute storing
#' query provenance.
#'
#' @seealso [summarise_availability()], [summarise_assembly_availability()],
#' [extract_assembly_metadata()], [summarise_sra_availability()],
#' [extract_sra_metadata()], [summarise_biosample_availability()],
#' [extract_biosample_metadata()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(
#'   c('Vigna angularis', 'Vigna reflexo-pilosa', 'Vigna vexillata'),
#'   synonyms = list('Vigna reflexo-pilosa' = 'Vigna glabrescens')
#' )
#' print(RESULTS)
#' }
#' @export
query_species <- function(species, synonyms = NULL) {
  SPEC <- .prepare_synonym_map(species, synonyms = synonyms)
  QUERY_TERMS <- unique(unlist(SPEC$query_terms, use.names = FALSE))
  n <- length(QUERY_TERMS)
  pb <- .pb_init(3L * n)
  SEARCHES <- lapply(seq_along(QUERY_TERMS), function(i) {
    sp <- QUERY_TERMS[i]
    tick <- 3L * (i - 1L)
    assembly <- .ncbi_search('assembly', sp)
    .pb_tick(pb, tick + 1L)
    sra <- .ncbi_search('sra', sp)
    .pb_tick(pb, tick + 2L)
    biosample <- .ncbi_search('biosample', sp)
    .pb_tick(pb, tick + 3L)
    list(
    assembly  = assembly,
    sra       = sra,
    biosample = biosample
    )
  })
  .pb_close(pb)
  names(SEARCHES) <- QUERY_TERMS
  RESULTS <- lapply(SPEC$species, function(sp) {
    terms <- SPEC$query_terms[[sp]]
    list(
    assembly  = .collapse_searches(lapply(terms, function(term) SEARCHES[[term]]$assembly)),
    sra       = .collapse_searches(lapply(terms, function(term) SEARCHES[[term]]$sra)),
    biosample = .collapse_searches(lapply(terms, function(term) SEARCHES[[term]]$biosample))
    )
  })
  names(RESULTS) <- SPEC$species
  merged <- names(SPEC$synonyms)[lengths(SPEC$synonyms) > 0L]
  if (length(merged) > 0L) {
    if (length(merged) <= 10L) {
      lab <- vapply(merged, function(sp) {
        paste0(sp, ' <- ', paste(SPEC$synonyms[[sp]], collapse = ', '))
      }, character(1))
      .gama_msg('Collapsing synonym queries under canonical species names: ', paste(lab, collapse = '; '), '.')
    } else {
      .gama_msg('Collapsing synonym queries for ', length(merged), ' species.')
    }
  }
  no_hits <- SPEC$species[vapply(RESULTS, function(x) {
    a <- x$assembly$count %||% 0L
    s <- x$sra$count %||% 0L
    b <- x$biosample$count %||% 0L
    (a == 0L) && (s == 0L) && (b == 0L)
  }, logical(1))]
  if (length(no_hits) > 0L) {
    if (length(no_hits) <= 10L) {
      .gama_msg('No records found across Assembly/SRA/BioSample for: ', paste(no_hits, collapse = ', '), '.')
    } else {
      .gama_msg('No records found across Assembly/SRA/BioSample for ', length(no_hits), ' of ', length(SPEC$species), ' species.')
    }
  }
  attr(RESULTS, 'query_info') <- .make_query_info(
  species     = SPEC$species,
  dbs         = c('assembly', 'sra', 'biosample'),
  query_terms = SPEC$query_terms,
  synonyms    = SPEC$synonyms
  )
  RESULTS <- .set_gama_object(RESULTS, 'query_species')
  RESULTS
}

#' Summarise record counts and compute data richness
#'
#' Collapses the output of [query_species()] into a species-level summary table
#' of record counts (Assembly, SRA, BioSample) and data richness scores.
#'
#' @details
#' For each species, record counts are extracted from the underlying query
#' results and summarised into a single row. The returned tibble is given class
#' `gdt_tbl` and retains query provenance via the `query_info` attribute
#' carried over from the `results` object.
#'
#' The `score` column is computed as `A + S + B`, where `A`, `S`, and `B`
#' are transformed contributions from Assembly, SRA, and BioSample
#' availability. `A = best + ln(1 + total - best)`, with assembly weights of
#' Complete = 10, Chromosome = 8, Scaffold = 5, and Contig = 2; `best` is
#' the maximum-weight assembly, with ties broken by highest N50, and `total`
#' is the sum of all assembly weights. `S = 2 * ln(1 + SRA)`, and
#' `B = ln(1 + BioSample)`. This formulation prioritises high-quality
#' assemblies while applying diminishing returns to highly sampled taxa.
#'
#' @param results A list returned by [query_species()].
#'
#' @return A tibble with one row per species and the following columns:
#' \itemize{
#' \item `species`: species name
#' \item `Assembly`, `SRA`, `BioSample`: record counts per database
#' \item `A`, `S`, `B`: component scores used for the composite
#' \item `score`: composite data richness score
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
#' }
#' @export
summarise_availability <- function(results) {
  results <- .gama_require_output(results, 'query_species')
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
  .as_gdt_table(OUT, results, 'summarise_availability')
}

#' Summarise Assembly composition
#'
#' Collapses Assembly metadata into species-level assembly level counts.
#'
#' @details
#' `summarise_assembly_availability()` reports the total number of Assembly
#' records returned by [query_species()] and the number assigned to each
#' recognised assembly level: `complete`, `chromosome`, `scaffold`, and
#' `contig`. Assembly level counts are derived from Assembly esummary metadata;
#' the `Assembly` total is retained from the original query record count.
#'
#' The `best_n50` column is the highest N50 among assemblies at the highest
#' available recognised assembly level for each species. Assembly levels are
#' ranked as complete > chromosome > scaffold > contig. Where multiple
#' assemblies are present at the highest available level, the highest N50
#' among those assemblies is reported. If no recognised assembly level with an
#' available N50 is found, `best_n50` is returned as `NA_real_`.
#'
#' @param results A list returned by [query_species()].
#' @param species `NULL` (default) to include all species, or a character
#' vector specifying which species to include.
#'
#' @return A tibble with one row per species and the following columns:
#' \itemize{
#' \item `species`: species name
#' \item `Assembly`: total Assembly record count
#' \item `complete`, `chromosome`, `scaffold`, `contig`: assembly level counts
#' \item `best_n50`: highest N50 among assemblies at the highest available
#' recognised assembly level
#' }
#' The tibble has class `gdt_tbl` and carries a `query_info` attribute for
#' provenance.
#'
#' @seealso [query_species()], [plot_assembly_availability()],
#' [extract_assembly_metadata()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' ASM_SUMMARY <- summarise_assembly_availability(RESULTS)
#' print(ASM_SUMMARY)
#' }
#' @export
summarise_assembly_availability <- function(results, species = NULL) {
  results <- .gama_require_output(results, 'query_species')
  SPECIES_ALL <- names(results)
  SPECIES_USE <- if (is.null(species)) {
    SPECIES_ALL
  } else {
    sp_in <- as.character(species)
    sp_in <- sp_in[!is.na(sp_in) & nzchar(sp_in)]
    sp_in <- unique(sp_in)
    missing <- setdiff(sp_in, SPECIES_ALL)
    if (length(missing) > 0L) .gama_warn('Requested species not found in `results`: ', paste(missing, collapse = ', '), '. Dropping.')
    sp_in[sp_in %in% SPECIES_ALL]
  }
  LEVELS <- .ASSEMBLY_CLASSES
  if (!length(SPECIES_USE)) {
    OUT <- tibble::tibble(
    species   = character(),
    Assembly  = integer(),
    complete  = integer(),
    chromosome = integer(),
    scaffold  = integer(),
    contig    = integer(),
    best_n50  = numeric()
    )
    return(.as_gdt_table(OUT, results, 'summarise_assembly_availability'))
  }
  META <- .assembly_metadata_core(results, species = SPECIES_USE)
  META$level_class <- .assembly_level_class(META$level)
  CLASS <- META |>
  dplyr::filter(!is.na(.data$level_class)) |>
  dplyr::count(species, level_class, name = 'count') |>
  tidyr::pivot_wider(
  names_from  = level_class,
  values_from = count,
  values_fill = 0
  )
  for (m in LEVELS) if (!m %in% names(CLASS)) CLASS[[m]] <- 0L
  BEST <- META |>
  dplyr::group_by(species) |>
  dplyr::summarise(
  best_n50 = .best_assembly_n50(.data$level, .data$n50),
  .groups  = 'drop'
  )
  COUNTS <- tibble::tibble(
  species  = SPECIES_USE,
  Assembly = purrr::map_int(results[SPECIES_USE], ~ .x$assembly$count %||% 0L)
  )
  OUT <- COUNTS |>
  dplyr::left_join(CLASS, by = 'species') |>
  dplyr::left_join(BEST, by = 'species') |>
  dplyr::select(species, Assembly, dplyr::all_of(LEVELS), best_n50)
  OUT <- OUT |>
  dplyr::mutate(dplyr::across(dplyr::all_of(LEVELS), ~ {
    y <- as.integer(.x)
    y[is.na(y)] <- 0L
    y
  }))
  no_data <- OUT$species[OUT$Assembly == 0L]
  if (length(no_data) > 0L) {
    if (length(no_data) <= 10L) {
      .gama_msg('No Assembly records found for: ', paste(no_data, collapse = ', '), '. Returning zeros for counts.')
    } else {
      .gama_msg('No Assembly records found for ', length(no_data), ' of ', length(SPECIES_USE), ' species. Returning zeros for counts.')
    }
  }
  .as_gdt_table(OUT, results, 'summarise_assembly_availability')
}

#' Extract filtered Assembly metadata
#'
#' Retrieves and structures assembly metadata for one or more species using the
#' Assembly identifiers stored in the output of [query_species()]. Metadata are
#' returned in a tidy tibble and optionally reduced to a single best assembly
#' per species.
#'
#' When `best = TRUE`, the best assembly is selected using structural-weighting
#' (assembly level) and N50 as a tie-breaker.
#'
#' @param results A list returned by [query_species()], containing Assembly
#' IDs.
#' @param species `NULL` (default) to return assemblies for all species, or a
#' character vector specifying which species to extract.
#' @param best Logical; if `TRUE`, return only the best assembly per species.
#'
#' @return A tibble with one row per assembly (or one row per species when
#' `best = TRUE`). Fields include species, entrez_uid, assembly level, N50,
#' coverage, BioSample/BioProject accessions, submitter, release date, and FTP
#' path (where available). The tibble has class `gdt_tbl` and carries a
#' `query_info` attribute for provenance.
#'
#' @seealso [query_species()], [summarise_assembly_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' ASM <- extract_assembly_metadata(
#'   RESULTS,
#'   species = 'Vigna angularis',
#'   best = TRUE
#' )
#' print(ASM)
#' }
#' @export
extract_assembly_metadata <- function(results, species = NULL, best = FALSE) {
  results <- .gama_require_output(results, 'query_species')
  best <- .gama_validate_logical_parameter(best, 'best')
  RESULTS_USE <- results
  if (!is.null(species)) {
    species <- as.character(species)
    species <- species[!is.na(species) & nzchar(species)]
    missing <- setdiff(species, names(results))
    if (length(missing) > 0L) .gama_stop('Requested species not found in `results`: ', paste(missing, collapse = ', '), '.')
    RESULTS_USE <- results[species]
  }
  META <- .assembly_metadata_core(RESULTS_USE)
  if (!best) return(.as_gdt_table(META, results, 'extract_assembly_metadata'))
  BEST <- lapply(split(META, META$species), function(x) {
    STRUCT <- .assembly_level_weight(x$level)
    if (!length(STRUCT) || all(STRUCT <= 0)) return(x[1L, , drop = FALSE])
    N50 <- x$n50
    max_struct <- max(STRUCT)
    tied_idx   <- which(STRUCT == max_struct)
    best_idx <- if (length(tied_idx) > 1) {
      tied_n50 <- N50[tied_idx]
      if (all(is.na(tied_n50))) tied_idx[1L] else tied_idx[which.max(tied_n50)]
    } else tied_idx
    x[best_idx, , drop = FALSE]
  })
  OUT <- dplyr::bind_rows(BEST)
  .as_gdt_table(OUT, results, 'extract_assembly_metadata')
}

#' Summarise SRA modality composition
#'
#' Collapses record-level SRA metadata into species-level modality counts
#' using classes and subclasses assigned using the internal ontology.
#'
#' By default, the output includes class-level counts for the major modalities
#' (genomic, transcriptomic, epigenomic, chromatin, other, unknown), plus the
#' total number of SRA records per species.
#'
#' Profile cache (used by downstream summaries):
#' In addition to the summary table, this function attaches a cached, UID-level
#' profile as an attribute `sra_profile`. The profile contains (at minimum)
#' `species`, `entrez_uid`, `biosample`, `bioproject`, `class`, `subclass`, and
#' GEO linkage fields (`geo_linked`, `gse_ids`, `gsm_ids`). This cache is
#' intended to be re-used locally for downstream summaries that require
#' within-species structure (e.g. replication skew across BioProjects or
#' BioSample IDs) without re-querying NCBI. In particular,
#' [summarise_sra_skew()] and [summarise_interaction()] consume
#' `attr(x, 'sra_profile')` from the output of this function.
#'
#' GEO overlay:
#' GEO linkage fields are always cached in `attr(x, 'sra_profile')` regardless
#' of `include_geo`. When `include_geo = TRUE`, species-level GEO summary
#' columns are appended as `<class>_geo` columns, alongside a GEO-compatible
#' denominator summary (`denom_total`, `geo_linked_denom`) and proportion
#' (`geo_prop`). The GEO-compatible denominator is defined as: transcriptomic +
#' epigenomic + chromatin + other + unknown.
#'
#' @param results A list returned by [query_species()].
#' @param species `NULL` (default) to include all species, or a character
#' vector specifying which species to include.
#' @param all Logical; if `TRUE`, include subclass-level columns.
#' @param include_geo Logical; if `TRUE`, append GEO summary columns (GEO
#' fields are cached regardless).
#'
#' @return A tibble with one row per species containing class-level counts and
#' totals. When `all = TRUE`, subclass-level counts are included. When
#' `include_geo = TRUE`, GEO summary columns are appended (e.g. `denom_total`,
#' `geo_linked_denom`, `geo_prop`, and `<class>_geo`). GEO linkage fields are
#' cached regardless in `attr(x, 'sra_profile')`. The tibble has class
#' `gdt_tbl` and carries a `query_info` attribute. It also carries a cached
#' UID-level profile as attribute `sra_profile` (see `Profile cache`), plus
#' metadata in `sra_profile_info`.
#'
#' @seealso [query_species()], [plot_sra_availability()], [plot_sra_geo()],
#' [summarise_sra_skew()], [extract_sra_metadata()]
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
  results <- .gama_require_output(results, 'query_species')
  all <- .gama_validate_logical_parameter(all, 'all')
  include_geo <- .gama_validate_logical_parameter(include_geo, 'include_geo')
  SPECIES_ALL <- names(results)
  SPECIES_USE <- if (is.null(species)) {
    SPECIES_ALL
  } else {
    sp_in <- as.character(species)
    sp_in <- sp_in[!is.na(sp_in) & nzchar(sp_in)]
    sp_in <- unique(sp_in)
    missing <- setdiff(sp_in, SPECIES_ALL)
    if (length(missing) > 0L) .gama_warn('Requested species not found in `results`: ', paste(missing, collapse = ', '), '. Dropping.')
    sp_in[sp_in %in% SPECIES_ALL]
  }
  if (!length(SPECIES_USE)) {
    MODES <- c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
    OUT <- tibble::tibble(species = character(), SRA = integer())
    for (m in MODES) OUT[[m]] <- integer()
    if (isTRUE(include_geo)) {
      OUT$denom_total <- integer()
      OUT$geo_linked_denom <- integer()
      OUT$geo_prop <- numeric()
      for (m in MODES) OUT[[paste0(m, '_geo')]] <- integer()
    }
    OUT <- .as_gdt_table(OUT, results, 'summarise_sra_availability')
    attr(OUT, 'sra_profile') <- tibble::tibble(
    species    = character(),
    entrez_uid = character(),
    biosample  = character(),
    bioproject = character(),
    class      = character(),
    subclass   = character(),
    geo_linked = logical(),
    gse_ids    = character(),
    gsm_ids    = character()
    )
    attr(OUT, 'sra_profile_info') <- list(
    cached_at_utc    = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
    profile_time_utc = attr(OUT, 'query_info')$query_time_utc %||% NA_character_,
    id_col           = 'entrez_uid',
    fields           = c('species', 'entrez_uid', 'biosample', 'bioproject', 'class', 'subclass', 'geo_linked', 'gse_ids', 'gsm_ids')
    )
    return(OUT)
  }
  META <- .sra_metadata_core(results, species = SPECIES_USE)
  PROFILE <- META |>
  dplyr::select(species, entrez_uid, biosample, bioproject, class, subclass, geo_linked, gse_ids, gsm_ids)
  CLASS <- META |>
  dplyr::count(species, class, name = 'count') |>
  tidyr::pivot_wider(
  names_from  = class,
  values_from = count,
  values_fill = 0
  )
  MODES <- c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
  for (m in MODES) if (!m %in% names(CLASS)) CLASS[[m]] <- 0L
  SRA_TOTAL <- META |>
  dplyr::count(species, name = 'SRA')
  OUT <- CLASS |>
  dplyr::left_join(SRA_TOTAL, by = 'species') |>
  dplyr::select(species, SRA, dplyr::all_of(MODES))
  if (isTRUE(all)) {
    SUB <- META |>
    dplyr::count(species, subclass, name = 'count') |>
    tidyr::pivot_wider(
    names_from  = subclass,
    values_from = count,
    values_fill = 0
    )
    overlap <- intersect(setdiff(names(SUB), 'species'), setdiff(names(OUT), 'species'))
    if (length(overlap) > 0L) names(SUB)[match(overlap, names(SUB))] <- paste0(overlap, '_subclass')
    OUT <- dplyr::left_join(OUT, SUB, by = 'species')
  }
  if (isTRUE(include_geo)) {
    denom <- c('transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
    DENOM <- META |>
    dplyr::filter(.data$class %in% denom) |>
    dplyr::count(species, name = 'denom_total')
    GEO_DENOM <- META |>
    dplyr::filter(.data$class %in% denom, .data$geo_linked) |>
    dplyr::count(species, name = 'geo_linked_denom')
    GEO_CLASS <- META |>
    dplyr::filter(.data$geo_linked) |>
    dplyr::count(species, class, name = 'count') |>
    tidyr::pivot_wider(
    names_from  = class,
    values_from = count,
    values_fill = 0
    )
    for (m in MODES) if (!m %in% names(GEO_CLASS)) GEO_CLASS[[m]] <- 0L
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
    GEO_SUM <- tibble::tibble(species = SPECIES_USE) |>
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
    OUT <- OUT |>
    dplyr::left_join(GEO_SUM, by = 'species')
  }
  OUT <- tibble::tibble(species = SPECIES_USE) |>
  dplyr::left_join(OUT, by = 'species')
  if (nrow(OUT) > 0L) {
    count_cols <- setdiff(names(OUT), c('species', 'geo_prop'))
    OUT <- OUT |>
    dplyr::mutate(dplyr::across(dplyr::all_of(count_cols), ~ {
      y <- as.integer(.x)
      y[is.na(y)] <- 0L
      y
    }))
  }
  if ('SRA' %in% names(OUT) && length(SPECIES_USE) > 0L) {
    no_data <- OUT$species[OUT$SRA == 0L]
    if (length(no_data) > 0L) {
      if (length(no_data) <= 10L) {
        .gama_msg('No SRA records found for: ', paste(no_data, collapse = ', '), '. Returning zeros for counts.')
      } else {
        .gama_msg('No SRA records found for ', length(no_data), ' of ', length(SPECIES_USE), ' species. Returning zeros for counts.')
      }
    }
  }
  OUT <- .as_gdt_table(OUT, results, 'summarise_sra_availability')
  attr(OUT, 'sra_profile') <- PROFILE
  attr(OUT, 'sra_profile_info') <- list(
  cached_at_utc    = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
  profile_time_utc = attr(OUT, 'query_info')$query_time_utc %||% NA_character_,
  id_col           = 'entrez_uid',
  fields           = c('species', 'entrez_uid', 'biosample', 'bioproject', 'class', 'subclass', 'geo_linked', 'gse_ids', 'gsm_ids')
  )
  OUT
}

#' Summarise SRA replication skew across BioProjects / BioSample IDs
#'
#' Quantifies replication skew across independent units (BioProject or
#' BioSample ID) using the cached UID-level SRA profile produced by
#' [summarise_sra_availability()]. Optional filters restrict the analysis to a
#' single modality class.
#'
#' Profile cache (consumed by this function):
#' The input must carry a cached UID-level profile as attribute `sra_profile`
#' containing (at minimum) `species`, `entrez_uid`, `biosample`, `bioproject`,
#' and `class`. Each row in the profile corresponds to an Entrez UID.
#'
#' @details
#' The `eff` column is the *effective number of units* (Hill number of order
#' 2), computed as the inverse Simpson index: `eff = 1 / sum(p^2)`, where `p`
#' is the proportion of SRA records in each BioProject or BioSample ID.
#' Larger values indicate a more even spread; values near 1 indicate strong
#' concentration in few units.
#'
#' Records missing BioProject or BioSample IDs are excluded from the skew
#' calculation. A `skew_id_recovery` attribute reports the active denominator,
#' number of records included, number excluded, and recovery proportion.
#'
#' @param x A data.frame/tibble returned by summarise_sra_availability()
#' that has a cached profile attached as attribute `sra_profile`.
#' @param species Optional character vector of species names to filter the
#' output. If `NULL`, all species in x are returned.
#' @param unit Character scalar; either `bioproject` (default) or `biosample`.
#' @param class Optional character scalar specifying a single modality class.
#'
#' @return A tibble/data.frame with one row per species containing:
#' `species`, `BioProject`/`BioSample` (number of distinct units with records),
#' `class`, `min`, `q25`, `med`, `q75`, `max` (SRA records per unit), and
#' `eff` (effective number of units; inverse Simpson index).
#'
#' @seealso [summarise_sra_availability()], [plot_sra_skew()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' SRA_SKEW <- summarise_sra_skew(SRA_SUMMARY, class = 'transcriptomic')
#' print(SRA_SKEW)
#' }
#'
#' @export
summarise_sra_skew <- function(x, species = NULL, unit = 'bioproject', class = NULL) {
  x <- .gama_require_output(x, 'summarise_sra_availability')
  unit <- .gama_validate_parameters(
    unit,
    'unit',
    .gama_sra_skew_unit_parameters,
    multiple = FALSE,
    allow_null = FALSE
  )
  valid_class <- c(names(.ONTOLOGY), 'unknown')
  class_parameters <- .sra_ontology_parameter_map(.ONTOLOGY, 'class')
  class <- .gama_validate_parameters(
    class,
    'class',
    valid_class,
    variants = class_parameters,
    multiple = FALSE
  )
  prof <- .gama_require_cache(
    x,
    attr_name = 'sra_profile',
    required_cols = c('species', 'class', 'biosample', 'bioproject'),
    source = 'summarise_sra_availability'
  )
  if (is.null(species)) {
    if ('species' %in% names(x)) {
      species_all <- sort(unique(as.character(x$species)))
    } else {
      species_all <- sort(unique(as.character(prof$species)))
    }
  } else {
    species_all <- as.character(species)
  }
  species_all <- species_all[!is.na(species_all) & nzchar(species_all)]
  unit_col <- if (identical(unit, 'bioproject')) 'bioproject' else 'biosample'
  out_unit_label <- if (identical(unit, 'bioproject')) 'BioProject' else 'BioSample'
  out_cols <- c('species', out_unit_label, 'class', 'min', 'q25', 'med', 'q75', 'max', 'eff')
  if (!length(species_all)) {
    out <- as.data.frame(
      stats::setNames(
        replicate(length(out_cols), vector('list', 0), simplify = FALSE),
        out_cols
      )
    )
    for (nm in setdiff(out_cols, c('species', out_unit_label, 'class'))) out[[nm]] <- numeric(0)
    out[['species']] <- character(0)
    out[[out_unit_label]] <- integer(0)
    out[['class']] <- character(0)
    out <- .as_gdt_table(out, x, 'summarise_sra_skew')
    attr(out, 'sra_profile') <- prof[0, , drop = FALSE]
    attr(out, 'sra_profile_info') <- attr(x, 'sra_profile_info', exact = TRUE)
    attr(out, 'skew_id_recovery') <- .skew_id_recovery_table(
      species = character(),
      unit = character(),
      class = character(),
      records = integer(),
      included = integer()
    )
    return(out)
  }
  prof_use <- prof[prof$species %in% species_all, , drop = FALSE]
  if (!is.null(class)) prof_use <- prof_use[prof_use$class %in% class, , drop = FALSE]
  unit_values <- trimws(as.character(prof_use[[unit_col]]))
  prof_use[[unit_col]] <- unit_values
  has_unit <- !is.na(unit_values) & nzchar(unit_values)
  prof_skew <- prof_use[has_unit, , drop = FALSE]
  class_label <- if (is.null(class)) 'all' else as.character(class)
  recovery <- .skew_id_recovery_table(
    species = species_all,
    unit = out_unit_label,
    class = class_label,
    records = vapply(species_all, function(sp) sum(prof_use$species == sp), integer(1)),
    included = vapply(species_all, function(sp) sum(prof_skew$species == sp), integer(1))
  )
  n_records <- sum(recovery$records)
  n_excluded <- sum(recovery$excluded)
  if (n_excluded > 0L) {
    .gama_msg(
      'Excluding ', .gama_format_count(n_excluded), ' of ',
      .gama_format_count(n_records), ' SRA records with missing ',
      out_unit_label, ' IDs from skew calculation.'
    )
  }
  calc_one <- function(sp) {
    units <- prof_skew[[unit_col]][prof_skew$species == sp]
    counts <- as.numeric(table(units))
    counts <- counts[counts > 0]
    if (length(counts) == 0L) {
      return(data.frame(
        species = sp,
        ucount = 0L,
        min = NA_real_,
        q25 = NA_real_,
        med = NA_real_,
        q75 = NA_real_,
        max = NA_real_,
        eff = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    qs <- as.numeric(
      stats::quantile(counts, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 7)
    )
    p <- counts / sum(counts)
    eff <- 1 / sum(p^2)
    data.frame(
      species = sp,
      ucount = length(counts),
      min = min(counts),
      q25 = qs[1],
      med = qs[2],
      q75 = qs[3],
      max = max(counts),
      eff = eff,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, lapply(species_all, calc_one))
  names(out)[names(out) == 'ucount'] <- out_unit_label
  out$class <- if (is.null(class)) 'all' else as.character(class)
  if (is.null(species)) {
    out <- out[order(out$species), , drop = FALSE]
  } else {
    out$species <- factor(out$species, levels = species_all)
    out <- out[order(out$species), , drop = FALSE]
    out$species <- as.character(out$species)
  }
  out <- out[c('species', out_unit_label, 'class', 'min', 'q25', 'med', 'q75', 'max', 'eff')]
  out <- .as_gdt_table(out, x, 'summarise_sra_skew')
  attr(out, 'sra_profile') <- prof_skew
  attr(out, 'sra_profile_info') <- attr(x, 'sra_profile_info', exact = TRUE)
  attr(out, 'skew_id_recovery') <- recovery
  out
}

#' Extract filtered SRA metadata
#'
#' Retrieves record-level SRA metadata for one or more species, then
#' assigns curated modality classes and subclasses together with GEO linkage
#' fields.
#'
#' Optional filters restrict the output to particular modality classes,
#' subclasses, or GEO-linked SRA records only.
#'
#' @details
#' Accepted modality filters follow the internal GAMA classification scheme.
#' Top-level classes and their subclasses are:
#'
#' \describe{
#' \item{`genomic`}{`amplicon-seq`, `clone-based`, `RAD-seq`,
#' `targeted-capture`, `WGS`}
#' \item{`transcriptomic`}{`long-read`, `RNA-seq`, `small-RNA`}
#' \item{`epigenomic`}{`ATAC-seq`, `bisulfite-seq`, `CUT&RUN`,
#' `CUT&Tag`, `ChIP-seq`, `DNase-seq`, `FAIRE-seq`, `MNase-seq`, `SELEX`}
#' \item{`chromatin`}{`3C-based`, `ChIA-PET`, `Hi-C`, `TCC`}
#' \item{`other`}{`other`}
#' \item{`unknown`}{`unknown`}
#' }
#'
#' @param results A list returned by [query_species()], containing SRA IDs.
#' @param species `NULL` (default) to include all species, or a character
#' vector specifying which species to include.
#' @param class Optional character vector of modality classes to retain.
#' @param subclass Optional character vector of modality subclasses to retain.
#' @param only_geo Logical; if `TRUE`, retain only GEO-linked SRA records.
#'
#' @return A tibble with one row per SRA record, including identifiers,
#' strategy fields, ontology assignments, and GEO linkage columns. The tibble
#' has class `gdt_tbl` and carries a `query_info` attribute for provenance.
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
#' }
#' @export
extract_sra_metadata <- function(results,
species  = NULL,
class    = NULL,
subclass = NULL,
only_geo = FALSE) {
  results <- .gama_require_output(results, 'query_species')
  only_geo <- .gama_validate_logical_parameter(only_geo, 'only_geo')
  if (!is.null(species)) {
    species <- as.character(species)
    species <- species[!is.na(species) & nzchar(species)]
    missing <- setdiff(species, names(results))
    if (length(missing) > 0L) .gama_warn('Requested species not found in `results`: ', paste(missing, collapse = ', '), '. Dropping.')
    species <- species[species %in% names(results)]
  }
  valid_class <- c(names(.ONTOLOGY), 'unknown')
  valid_subclass <- c(
    unique(unlist(lapply(.ONTOLOGY, names), use.names = FALSE)),
    'unknown'
  )
  class_parameters <- .sra_ontology_parameter_map(.ONTOLOGY, 'class')
  subclass_parameters <- .sra_ontology_parameter_map(.ONTOLOGY, 'subclass')
  class <- .gama_validate_parameters(
    class,
    'class',
    valid_class,
    variants = class_parameters
  )
  subclass <- .gama_validate_parameters(
    subclass,
    'subclass',
    valid_subclass,
    variants = subclass_parameters
  )
  if (!is.null(species) && !length(species)) {
    META <- tibble::tibble(
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
    )
    return(.as_gdt_table(META, results, 'extract_sra_metadata'))
  }
  META <- .sra_metadata_core(results, species = species)
  if (!is.null(class)) {
    META <- META |> dplyr::filter(.data$class %in% .env$class)
  }
  if (!is.null(subclass)) {
    META <- META |> dplyr::filter(.data$subclass %in% .env$subclass)
  }
  if (isTRUE(only_geo)) {
    META <- META |> dplyr::filter(.data$geo_linked)
  }
  if (nrow(META) == 0L && (!is.null(species) || !is.null(class) || !is.null(subclass) || isTRUE(only_geo))) {
    .gama_msg('No SRA metadata records found for requested filters; returning empty table.')
  }
  .as_gdt_table(META, results, 'extract_sra_metadata')
}

#' Summarise BioSample anatomy composition
#'
#' Collapses BioSample sample-source metadata into species-level counts of
#' anatomy classes, subclasses, and canonical anatomy terms assigned using the
#' internal GAMA BioSample ontology.
#'
#' By default, the output includes class-level counts for the major anatomy
#' categories (`aerial`, `ground`, `reproductive`, `whole`, `in_vitro`,
#' `other`, `mixed`, and `unknown`), plus total and operable BioSample record
#' counts per species. Here, `whole` is reserved for whole-plant or
#' whole-seedling samples, not for whole sub-anatomy such as whole leaves,
#' roots, grains, or spikelets.
#'
#' A BioSample record is operable when it contains at least one accepted
#' sample-source attribute, regardless of whether the value can be assigned to
#' a curated anatomy term. Missing-like values and values with no ontological
#' match are retained as `unknown`. Broad, generic, or subcellular anatomy
#' values are retained as `other`.
#'
#' At class level, `mixed` is assigned when more than one anatomy class is
#' recovered from accepted fields within the same BioSample record.
#' Multi-subclass records within a single class remain assigned to that class
#' in this summary. Subclass-level mixed states are handled when
#' [summarise_interaction()] is called with `level = 'anatomy_subclass'`.
#'
#' Profile cache:
#' In addition to the summary table, this function attaches two cached
#' BioSample-level profiles: `biosample_anatomy_profile` and
#' `biosample_canonical_profile`. The anatomy profile stores one collapsed
#' anatomy class and subclass profile per operable BioSample record, together
#' with BioProject provenance where available. The canonical profile stores
#' row-level anatomy terms, subclasses, classes, ranks, and ontology fields.
#' These caches are reused by [summarise_biosample_skew()] and
#' [summarise_interaction()] without repeating the BioSample retrieval and
#' parsing stage.
#'
#' @param results A list returned by [query_species()].
#' @param species `NULL` (default) to include all species, or a character
#' vector specifying which species to include.
#' @param all Logical; if `TRUE`, include canonical anatomy-term columns.
#'
#' @return A tibble with one row per species containing total BioSample record
#' counts, operable BioSample record counts, and class-level anatomy counts.
#' When `all = TRUE`, canonical anatomy-term counts are also included. The
#' tibble has class `gdt_tbl` and carries a `query_info` attribute. It also
#' carries cached BioSample-level profiles as attributes
#' `biosample_anatomy_profile` and `biosample_canonical_profile`, plus metadata
#' in `biosample_anatomy_profile_info` and `biosample_canonical_profile_info`.
#'
#' @seealso [query_species()], [plot_biosample_availability()],
#' [summarise_biosample_skew()], [extract_biosample_metadata()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' BIO_SUMMARY <- summarise_biosample_availability(RESULTS)
#' print(BIO_SUMMARY)
#' }
#' @export
summarise_biosample_availability <- function(results, species = NULL, all = FALSE) {
  results <- .gama_require_output(results, 'query_species')
  all <- .gama_validate_logical_parameter(all, 'all')
  anatomy_levels <- .biosample_anatomy_profile_levels()
  canonical_fields <- c(
    'species',
    'biosample_id',
    'bioproject',
    'value_raw',
    'value_norm',
    'anatomy_term',
    'anatomy_class',
    'anatomy_subclass',
    'rank',
    'ontology_namespace',
    'ontology_id',
    'ontology_label'
  )
  SPECIES_ALL <- names(results)
  if (is.null(SPECIES_ALL) || any(!nzchar(SPECIES_ALL))) {
    .gama_input_error(
      'query_species',
      detected = .detect_gama_object(results),
      detail = 'missing species names.'
    )
  }
  SPECIES_USE <- if (is.null(species)) {
    SPECIES_ALL
  } else {
    sp_in <- unique(as.character(species))
    sp_in <- sp_in[!is.na(sp_in) & nzchar(sp_in)]
    missing <- setdiff(sp_in, SPECIES_ALL)
    if (length(missing) > 0L) .gama_warn('Requested species not found in `results`: ', paste(missing, collapse = ', '), '. Dropping.')
    sp_in[sp_in %in% SPECIES_ALL]
  }
  if (!length(SPECIES_USE)) {
    OUT <- tibble::tibble(
      species = character(),
      BioSample = integer(),
      operable = integer(),
      aerial = integer(),
      ground = integer(),
      reproductive = integer(),
      whole = integer(),
      in_vitro = integer(),
      other = integer(),
      mixed = integer(),
      unknown = integer()
    )
    OUT <- .as_gdt_table(OUT, results, 'summarise_biosample_availability')
    attr(OUT, 'biosample_anatomy_profile') <- tibble::tibble(
      species = character(),
      biosample_id = character(),
      bioproject = character(),
      anatomy_class = character(),
      anatomy_subclass = character()
    )
    attr(OUT, 'biosample_anatomy_profile_info') <- list(
      cached_at_utc = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
      id_col = 'biosample_id',
      fields = c('species', 'biosample_id', 'bioproject', 'anatomy_class', 'anatomy_subclass'),
      anatomy_levels = anatomy_levels
    )
    attr(OUT, 'biosample_canonical_profile') <- tibble::tibble(
      species = character(),
      biosample_id = character(),
      bioproject = character(),
      value_raw = character(),
      value_norm = character(),
      anatomy_term = character(),
      anatomy_class = character(),
      anatomy_subclass = character(),
      rank = character(),
      ontology_namespace = character(),
      ontology_id = character(),
      ontology_label = character()
    )
    attr(OUT, 'biosample_canonical_profile_info') <- list(
      cached_at_utc = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
      id_col = 'biosample_id',
      fields = canonical_fields
    )
    return(OUT)
  }
  OUT <- tibble::tibble(
    species = SPECIES_USE,
    BioSample = purrr::map_int(results[SPECIES_USE], ~ .x$biosample$count %||% 0L)
  )
  PROGRESS <- .pb_state(
    .biosample_progress_total(results[SPECIES_USE], post_steps = 4L),
    enabled = interactive()
  )
  if (!is.null(PROGRESS)) on.exit(.pb_close_state(PROGRESS), add = TRUE)
  META <- .biosample_metadata_core(
    results = results,
    species = SPECIES_USE,
    progress = FALSE,
    progress_state = PROGRESS,
    include_empty_records = FALSE,
    keep_bioproject = TRUE
  )
  CLASSIFIED <- .biosample_classify_tissue_attributes(META, species = SPECIES_USE)
  .pb_advance(PROGRESS)
  ANATOMY <- .biosample_anatomy_profile(CLASSIFIED, species = SPECIES_USE)
  .pb_advance(PROGRESS)
  CANONICAL <- .biosample_canonical_profile(CLASSIFIED, species = SPECIES_USE)
  .pb_advance(PROGRESS)
  operable_tbl <- ANATOMY |>
    dplyr::distinct(.data$species, .data$biosample_id) |>
    dplyr::count(.data$species, name = 'operable')
  bucket_wide <- ANATOMY |>
    dplyr::count(.data$species, .data$anatomy_class, name = 'n') |>
    tidyr::pivot_wider(
      names_from = anatomy_class,
      values_from = n,
      values_fill = 0
    )
  if (!nrow(bucket_wide)) bucket_wide <- tibble::tibble(species = SPECIES_USE)
  for (nm in anatomy_levels) if (!nm %in% names(bucket_wide)) bucket_wide[[nm]] <- 0L
  OUT <- OUT |>
    dplyr::left_join(operable_tbl, by = 'species') |>
    dplyr::left_join(bucket_wide, by = 'species') |>
    dplyr::mutate(
      operable = as.integer(dplyr::coalesce(.data$operable, 0L)),
      aerial = as.integer(dplyr::coalesce(.data$aerial, 0L)),
      ground = as.integer(dplyr::coalesce(.data$ground, 0L)),
      reproductive = as.integer(dplyr::coalesce(.data$reproductive, 0L)),
      whole = as.integer(dplyr::coalesce(.data$whole, 0L)),
      in_vitro = as.integer(dplyr::coalesce(.data$in_vitro, 0L)),
      other = as.integer(dplyr::coalesce(.data$other, 0L)),
      mixed = as.integer(dplyr::coalesce(.data$mixed, 0L)),
      unknown = as.integer(dplyr::coalesce(.data$unknown, 0L))
    ) |>
    dplyr::select(species, BioSample, operable, aerial, ground, reproductive, whole, in_vitro, other, mixed, unknown)
  if (isTRUE(all) && nrow(CANONICAL)) {
    CAN <- CANONICAL |>
      dplyr::filter(.data$anatomy_class != 'unknown', .data$anatomy_term != 'unknown') |>
      dplyr::count(.data$species, .data$anatomy_term, name = 'n') |>
      tidyr::pivot_wider(
        names_from = anatomy_term,
        values_from = n,
        values_fill = 0
      )
    if (nrow(CAN)) {
      OUT <- dplyr::left_join(OUT, CAN, by = 'species')
      major_cols <- c('species', 'BioSample', 'operable', .biosample_anatomy_profile_levels())
      extra_cols <- sort(setdiff(names(OUT), major_cols))
      OUT <- OUT[, c(major_cols, extra_cols), drop = FALSE]
    }
  }
  no_biosample <- OUT$species[OUT$BioSample == 0L]
  if (length(no_biosample) > 0L) {
    if (length(no_biosample) <= 10L) {
      .gama_msg('No BioSample records found for: ', paste(no_biosample, collapse = ', '), '. Returning zeros for counts.')
    } else {
      .gama_msg('No BioSample records found for ', length(no_biosample), ' of ', length(SPECIES_USE), ' species. Returning zeros for counts.')
    }
  }
  no_operable <- OUT$species[OUT$BioSample > 0L & OUT$operable == 0L]
  if (length(no_operable) > 0L) {
    if (length(no_operable) <= 10L) {
      .gama_msg('No operable BioSample records found for: ', paste(no_operable, collapse = ', '), '. Returning zeros for anatomy counts.')
    } else {
      .gama_msg('No operable BioSample records found for ', length(no_operable), ' of ', length(SPECIES_USE), ' species. Returning zeros for anatomy counts.')
    }
  }
  OUT <- .as_gdt_table(OUT, results, 'summarise_biosample_availability')
  attr(OUT, 'biosample_anatomy_profile') <- ANATOMY
  attr(OUT, 'biosample_anatomy_profile_info') <- list(
    cached_at_utc = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
    profile_time_utc = attr(OUT, 'query_info')$query_time_utc %||% NA_character_,
    id_col = 'biosample_id',
    fields = c('species', 'biosample_id', 'bioproject', 'anatomy_class', 'anatomy_subclass'),
    anatomy_levels = anatomy_levels
  )
  attr(OUT, 'biosample_canonical_profile') <- CANONICAL
  attr(OUT, 'biosample_canonical_profile_info') <- list(
    cached_at_utc = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
    profile_time_utc = attr(OUT, 'query_info')$query_time_utc %||% NA_character_,
    id_col = 'biosample_id',
    fields = canonical_fields
  )
  .pb_advance(PROGRESS)
  OUT
}

#' Summarise BioSample replication skew across BioProjects
#'
#' Quantifies replication skew across BioProjects using the cached
#' BioSample-level anatomy profile produced by
#' [summarise_biosample_availability()]. Optional filters restrict the analysis
#' to a single anatomy class.
#'
#' Profile cache (consumed by this function):
#' The input must carry a cached BioSample-level profile as attribute
#' `biosample_anatomy_profile` containing (at minimum) `species`,
#' `biosample_id`, `bioproject`, and `anatomy_class`. Each row in the profile
#' corresponds to one operable BioSample record.
#'
#' @details
#' The `eff` column is the *effective number of BioProjects* (Hill number of
#' order 2), computed as the inverse Simpson index: `eff = 1 / sum(p^2)`,
#' where `p` is the proportion of operable BioSample records in each
#' BioProject. Larger values indicate a more even spread; values near 1
#' indicate strong concentration in few BioProjects.
#'
#' Records missing BioProject IDs are excluded from the skew calculation. A
#' `skew_id_recovery` attribute reports the active denominator, number of
#' records included, number excluded, and recovery proportion.
#'
#' @param x A data.frame/tibble returned by
#' [summarise_biosample_availability()] that has a cached profile attached as
#' attribute `biosample_anatomy_profile`.
#' @param species Optional character vector of species names to filter the
#' output. If `NULL`, all species in x are returned.
#' @param anatomy_class Optional character scalar specifying a single anatomy
#' class.
#'
#' @return A tibble/data.frame with one row per species containing:
#' `species`, `BioProject` (number of distinct BioProjects with operable
#' BioSample records), `anatomy_class`, `min`, `q25`, `med`, `q75`, `max`
#' (operable BioSample records per BioProject), and `eff` (effective number of
#' BioProjects; inverse Simpson index).
#'
#' @seealso [summarise_biosample_availability()], [plot_biosample_skew()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' BIO_SUMMARY <- summarise_biosample_availability(RESULTS)
#' BIO_SKEW <- summarise_biosample_skew(BIO_SUMMARY, anatomy_class = 'aerial')
#' print(BIO_SKEW)
#' }
#'
#' @export
summarise_biosample_skew <- function(x, species = NULL, anatomy_class = NULL) {
  x <- .gama_require_output(x, 'summarise_biosample_availability')
  valid_anatomy_class <- .biosample_anatomy_profile_levels()
  anatomy_class_parameters <- .biosample_anatomy_ontology_parameter_map('class')
  anatomy_class <- .gama_validate_parameters(
    anatomy_class,
    'anatomy_class',
    valid_anatomy_class,
    variants = anatomy_class_parameters,
    multiple = FALSE
  )
  prof <- .gama_require_cache(
    x,
    attr_name = 'biosample_anatomy_profile',
    required_cols = c('species', 'biosample_id', 'bioproject', 'anatomy_class'),
    source = 'summarise_biosample_availability'
  )
  if (is.null(species)) {
    if ('species' %in% names(x)) {
      species_all <- sort(unique(as.character(x$species)))
    } else {
      species_all <- sort(unique(as.character(prof$species)))
    }
  } else {
    species_all <- as.character(species)
  }
  species_all <- species_all[!is.na(species_all) & nzchar(species_all)]
  out_cols <- c('species', 'BioProject', 'anatomy_class', 'min', 'q25', 'med', 'q75', 'max', 'eff')
  if (!length(species_all)) {
    out <- as.data.frame(
      stats::setNames(
        replicate(length(out_cols), vector('list', 0), simplify = FALSE),
        out_cols
      )
    )
    for (nm in setdiff(out_cols, c('species', 'BioProject', 'anatomy_class'))) out[[nm]] <- numeric(0)
    out[['species']] <- character(0)
    out[['BioProject']] <- integer(0)
    out[['anatomy_class']] <- character(0)
    out <- .as_gdt_table(out, x, 'summarise_biosample_skew')
    attr(out, 'biosample_anatomy_profile') <- prof[0, , drop = FALSE]
    attr(out, 'biosample_anatomy_profile_info') <- attr(x, 'biosample_anatomy_profile_info', exact = TRUE)
    attr(out, 'skew_id_recovery') <- .skew_id_recovery_table(
      species = character(),
      unit = character(),
      class = character(),
      records = integer(),
      included = integer()
    )
    return(out)
  }
  prof_use <- prof[prof$species %in% species_all, , drop = FALSE]
  if (!is.null(anatomy_class)) prof_use <- prof_use[prof_use$anatomy_class %in% anatomy_class, , drop = FALSE]
  unit_values <- trimws(as.character(prof_use$bioproject))
  prof_use$bioproject <- unit_values
  has_unit <- !is.na(unit_values) & nzchar(unit_values)
  prof_skew <- prof_use[has_unit, , drop = FALSE]
  anatomy_class_label <- if (is.null(anatomy_class)) 'all' else as.character(anatomy_class)
  recovery <- .skew_id_recovery_table(
    species = species_all,
    unit = 'BioProject',
    class = anatomy_class_label,
    records = vapply(species_all, function(sp) sum(prof_use$species == sp), integer(1)),
    included = vapply(species_all, function(sp) sum(prof_skew$species == sp), integer(1))
  )
  n_records <- sum(recovery$records)
  n_excluded <- sum(recovery$excluded)
  if (n_excluded > 0L) {
    .gama_msg(
      'Excluding ', .gama_format_count(n_excluded), ' of ',
      .gama_format_count(n_records),
      ' operable BioSample records with missing BioProject IDs from skew calculation.'
    )
  }
  calc_one <- function(sp) {
    units <- prof_skew$bioproject[prof_skew$species == sp]
    counts <- as.numeric(table(units))
    counts <- counts[counts > 0]
    if (length(counts) == 0L) {
      return(data.frame(
        species = sp,
        ucount = 0L,
        min = NA_real_,
        q25 = NA_real_,
        med = NA_real_,
        q75 = NA_real_,
        max = NA_real_,
        eff = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    qs <- as.numeric(
      stats::quantile(counts, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 7)
    )
    p <- counts / sum(counts)
    eff <- 1 / sum(p^2)
    data.frame(
      species = sp,
      ucount = length(counts),
      min = min(counts),
      q25 = qs[1],
      med = qs[2],
      q75 = qs[3],
      max = max(counts),
      eff = eff,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, lapply(species_all, calc_one))
  names(out)[names(out) == 'ucount'] <- 'BioProject'
  out$anatomy_class <- anatomy_class_label
  if (is.null(species)) {
    out <- out[order(out$species), , drop = FALSE]
  } else {
    out$species <- factor(out$species, levels = species_all)
    out <- out[order(out$species), , drop = FALSE]
    out$species <- as.character(out$species)
  }
  out <- out[c('species', 'BioProject', 'anatomy_class', 'min', 'q25', 'med', 'q75', 'max', 'eff')]
  out <- .as_gdt_table(out, x, 'summarise_biosample_skew')
  attr(out, 'biosample_anatomy_profile') <- prof_skew
  attr(out, 'biosample_anatomy_profile_info') <- attr(x, 'biosample_anatomy_profile_info', exact = TRUE)
  attr(out, 'skew_id_recovery') <- recovery
  out
}

#' Extract filtered BioSample metadata
#'
#' Retrieves BioSample sample-source metadata for one or more species, then
#' assigns curated anatomy classes, subclasses, and canonical anatomy terms
#' using the internal GAMA BioSample ontology.
#'
#' Optional filters restrict the output to selected anatomy classes,
#' subclasses, or canonical anatomy terms.
#'
#' @details
#' Accepted anatomy filters follow the internal GAMA BioSample anatomy scheme.
#' Top-level classes include `aerial`, `ground`, `reproductive`, `whole`,
#' `in_vitro`, `other`, and `unknown`; `whole` denotes whole-plant or
#' whole-seedling samples rather than whole sub-anatomy. Subclasses and terms
#' follow the curated ontology used by [summarise_biosample_availability()].
#'
#' Missing-like values and values with no ontological match are retained as
#' `unknown`. Broad, generic, or subcellular anatomy values are retained as
#' `other`.
#'
#' This function preserves term-level anatomy assignments. A single BioSample
#' record may therefore appear in multiple rows when more than one anatomy term
#' is recovered. The `anatomy_class_profile` and `anatomy_subclass_profile`
#' columns provide the corresponding collapsed BioSample-level class and
#' subclass states, including `mixed` where multiple classes or subclasses are
#' recovered from the same BioSample record.
#'
#' @param results A list returned by [query_species()], containing BioSample
#' IDs.
#' @param species `NULL` (default) to include all species, or a character
#' vector specifying which species to include.
#' @param anatomy_class Optional character vector of anatomy classes to retain.
#' @param anatomy_subclass Optional character vector of anatomy subclasses to
#' retain.
#' @param anatomy_term Optional character vector of canonical anatomy terms to
#' retain.
#'
#' @return A tibble with one row per recovered BioSample anatomy term. Columns
#' include species, Entrez UID, BioSample accession, BioProject accession, raw
#' and normalised sample-source values, row-level anatomy assignments, and
#' collapsed BioSample-level `anatomy_class_profile` and
#' `anatomy_subclass_profile` labels. The tibble has class `gdt_tbl` and
#' carries a `query_info` attribute for provenance.
#'
#' @seealso [query_species()], [summarise_biosample_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' BIO <- extract_biosample_metadata(
#'   RESULTS,
#'   species = 'Vigna vexillata',
#'   anatomy_subclass = 'leaf'
#' )
#' print(BIO)
#' }
#' @export
extract_biosample_metadata <- function(results,
species          = NULL,
anatomy_class    = NULL,
anatomy_subclass = NULL,
anatomy_term     = NULL) {
  results <- .gama_require_output(results, 'query_species')
  empty <- function() {
    tibble::tibble(
    species = character(),
    entrez_uid = character(),
    biosample = character(),
    bioproject = character(),
    value_raw = character(),
    value_norm = character(),
    anatomy_class = character(),
    anatomy_subclass = character(),
    anatomy_term = character(),
    anatomy_class_profile = character(),
    anatomy_subclass_profile = character()
    )
  }
  attach_biosample_attrs <- function(out, meta) {
    attr(out, 'biosample_records') <- attr(meta, 'biosample_records', exact = TRUE)
    attr(out, 'biosample_diagnostics') <- attr(meta, 'biosample_diagnostics', exact = TRUE)
    attr(out, 'biosample_esummary_cache') <- attr(meta, 'biosample_esummary_cache', exact = TRUE)
    out
  }
  if (!is.null(species)) {
    species <- as.character(species)
    species <- species[!is.na(species) & nzchar(species)]
    missing <- setdiff(species, names(results))
    if (length(missing) > 0L) .gama_warn('Requested species not found in `results`: ', paste(missing, collapse = ', '), '. Dropping.')
    species <- species[species %in% names(results)]
  }
  valid_anatomy_class <- c(.biosample_anatomy_class_levels(), 'unknown')
  valid_anatomy_subclass <- unique(c(.biosample_anatomy_subclass_levels(), 'unknown'))
  valid_anatomy_term <- unique(c(.biosample_anatomy_ref()$lex$anatomy_term, 'unknown'))
  anatomy_class_parameters <- .biosample_anatomy_ontology_parameter_map('class')
  anatomy_subclass_parameters <- .biosample_anatomy_ontology_parameter_map('subclass')
  anatomy_term_parameters <- .biosample_anatomy_ontology_parameter_map('term')
  anatomy_class <- .gama_validate_parameters(
    anatomy_class,
    'anatomy_class',
    valid_anatomy_class,
    variants = anatomy_class_parameters
  )
  anatomy_subclass <- .gama_validate_parameters(
    anatomy_subclass,
    'anatomy_subclass',
    valid_anatomy_subclass,
    variants = anatomy_subclass_parameters
  )
  anatomy_term <- .gama_validate_parameters(
    anatomy_term,
    'anatomy_term',
    valid_anatomy_term,
    variants = anatomy_term_parameters,
    show_values = FALSE,
    no_suggestion_hint = 'Try `anatomy_class` or `anatomy_subclass` for broader filtering.'
  )
  if (!is.null(species) && !length(species)) {
    OUT <- .as_gdt_table(empty(), results, 'extract_biosample_metadata')
    return(OUT)
  }
  SPECIES_USE <- if (is.null(species)) names(results) else species
  PROGRESS <- .pb_state(
    .biosample_progress_total(results[SPECIES_USE], post_steps = 2L),
    enabled = interactive()
  )
  if (!is.null(PROGRESS)) on.exit(.pb_close_state(PROGRESS), add = TRUE)
  META <- .biosample_metadata_core(
    results = results,
    species = species,
    progress = FALSE,
    progress_state = PROGRESS,
    keep_xml = FALSE,
    keep_bioproject = TRUE,
    include_empty_records = FALSE
  )
  HITS <- .biosample_classify_tissue_attributes(META, species = species) |>
    dplyr::transmute(
      species = .data$species,
      input_id = .data$input_id,
      entrez_uid = .data$entrez_uid,
      biosample = .data$biosample_id,
      bioproject = .data$bioproject,
      value_raw = .data$value_raw,
      value_norm = .data$value_norm,
      anatomy_class = .data$anatomy_class,
      anatomy_subclass = .data$anatomy_subclass,
      anatomy_term = .data$anatomy_term
    )
  .pb_advance(PROGRESS)
  if (!nrow(HITS)) {
    OUT <- .as_gdt_table(empty(), results, 'extract_biosample_metadata')
    OUT <- attach_biosample_attrs(OUT, META)
    return(OUT)
  }
  PROFILE <- HITS |>
    dplyr::group_by(.data$species, .data$biosample) |>
    dplyr::summarise(
      anatomy_class_profile = .biosample_collapse_anatomy_profile(.data$anatomy_class, level = 'anatomy_class'),
      anatomy_subclass_profile = .biosample_collapse_anatomy_profile(.data$anatomy_subclass, level = 'anatomy_subclass'),
      .groups = 'drop'
    )
  if (!is.null(anatomy_class)) {
    HITS <- HITS |> dplyr::filter(.data$anatomy_class %in% .env$anatomy_class)
  }
  if (!is.null(anatomy_subclass)) {
    HITS <- HITS |> dplyr::filter(.data$anatomy_subclass %in% .env$anatomy_subclass)
  }
  if (!is.null(anatomy_term)) {
    HITS <- HITS |> dplyr::filter(.data$anatomy_term %in% .env$anatomy_term)
  }
  OUT <- HITS |>
    dplyr::left_join(
      PROFILE,
      by = c('species', 'biosample')
    ) |>
    dplyr::transmute(
      species = .data$species,
      entrez_uid = .data$entrez_uid,
      biosample = .data$biosample,
      bioproject = .data$bioproject,
      value_raw = .data$value_raw,
      value_norm = .data$value_norm,
      anatomy_class = .data$anatomy_class,
      anatomy_subclass = .data$anatomy_subclass,
      anatomy_term = .data$anatomy_term,
      anatomy_class_profile = .data$anatomy_class_profile,
      anatomy_subclass_profile = .data$anatomy_subclass_profile
    ) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$species, .data$biosample, .data$anatomy_class, .data$anatomy_subclass, .data$anatomy_term)
  OUT <- .as_gdt_table(OUT, results, 'extract_biosample_metadata')
  OUT <- attach_biosample_attrs(OUT, META)
  .pb_advance(PROGRESS)
  OUT
}

#' Summarise interaction profiles
#'
#' Links cached SRA modality and BioSample anatomy profiles to summarise
#' modality-by-anatomy structure for each species.
#'
#' Profile caches:
#' The input `SRA` must carry the cached UID-level profile produced by
#' [summarise_sra_availability()]. The input `BIO` must carry the cached
#' BioSample-level anatomy profiles produced by
#' [summarise_biosample_availability()]. BioSample records are linked by shared
#' BioSample identifiers, then counted across SRA modality classes and either
#' BioSample anatomy classes or anatomy subclasses.
#'
#' @details
#' For each species, the `BioSample` column gives the observed number of linked
#' BioSample records in each modality-by-anatomy cell. Expected counts are
#' computed from the species-level marginal totals as
#' `expected = row_total * col_total / total`, where `row_total` is the total
#' number of linked records in a modality class, `col_total` is the total
#' number of linked records in an anatomy category, and `total` is the total
#' number of linked records for that species.
#'
#' The `residual` column is computed as
#' `(BioSample - expected) / sqrt(expected)`. Positive residuals indicate
#' modality-anatomy combinations that are over-represented relative to the
#' marginal distributions; negative residuals indicate under-represented
#' combinations.
#'
#' The meaning of `mixed` is tied to the requested anatomy resolution: with
#' `level = 'anatomy_class'`, `mixed` means more than one anatomy class was
#' recovered for a BioSample; with `level = 'anatomy_subclass'`, `mixed` means
#' more than one anatomy subclass was recovered.
#'
#' @param SRA A summary table returned by [summarise_sra_availability()].
#' @param BIO A summary table returned by [summarise_biosample_availability()].
#' @param level Anatomy resolution to summarise. One of `anatomy_class` or
#' `anatomy_subclass`.
#' @param species `NULL` (default) to include all species, or a character
#' vector specifying which species to include.
#'
#' @return A tibble with one row per species, modality class, and anatomy
#' category. Columns include `species`, `class`, `BioSample`, `expected`, and
#' `residual`, plus either `anatomy_class` or `anatomy_subclass`, depending on
#' `level`. The tibble has class `gdt_tbl`, carries a `query_info` attribute,
#' and includes interaction metadata in `interaction_info`.
#'
#' @seealso [summarise_sra_availability()],
#' [summarise_biosample_availability()], [plot_interaction()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' BIO_SUMMARY <- summarise_biosample_availability(RESULTS)
#' INTERACTION <- summarise_interaction(
#'   SRA_SUMMARY,
#'   BIO_SUMMARY,
#'   species = 'Vigna angularis'
#' )
#' print(INTERACTION)
#' }
#' @export
summarise_interaction <- function(SRA,
                                  BIO,
                                  level = 'anatomy_class',
                                  species = NULL) {
  level <- .gama_validate_parameters(level, 'level', c('anatomy_class', 'anatomy_subclass'), multiple = FALSE, allow_null = FALSE)
  BIO <- .gama_require_output(
    BIO,
    'summarise_biosample_availability',
    required_cols = c('species', 'BioSample')
  )
  SRA <- .gama_require_output(SRA, 'summarise_sra_availability')
  anatomy <- .gama_require_cache(
    BIO,
    attr_name = 'biosample_anatomy_profile',
    required_cols = c('species', 'biosample_id', 'anatomy_class'),
    source = 'summarise_biosample_availability'
  )
  sra_prof <- .gama_require_cache(
    SRA,
    attr_name = 'sra_profile',
    required_cols = c('species', 'class'),
    source = 'summarise_sra_availability'
  )
  if (!any(c('biosample', 'biosample_id') %in% names(sra_prof))) {
    .gama_input_error(
      'summarise_sra_availability',
      detected = .detect_gama_object(SRA),
      detail = 'cache \'sra_profile\' missing required BioSample ID column: biosample or biosample_id.'
    )
  }
  species_all <- unique(as.character(BIO$species))
  species_all <- species_all[!is.na(species_all) & nzchar(species_all)]
  if (!length(species_all)) .gama_stop('No species found in `BIO`.')
  if (is.null(species)) {
    species_use <- species_all
  } else {
    species_in <- unique(as.character(species))
    species_in <- species_in[!is.na(species_in) & nzchar(species_in)]
    missing_species <- setdiff(species_in, species_all)
    if (length(missing_species)) {
      .gama_warn('Requested species not found in input `BIO`: ', paste(missing_species, collapse = ', '), '. Dropping.')
    }
    species_use <- species_in[species_in %in% species_all]
    if (!length(species_use)) .gama_stop('No matching species found.')
  }
  add_residuals <- function(MAT, term_col) {
    MAT2 <- MAT |>
      dplyr::mutate(.interaction_term = as.character(.data[[term_col]]))
    total_by_species <- MAT2 |>
      dplyr::group_by(.data$species) |>
      dplyr::summarise(total_n = sum(.data$BioSample, na.rm = TRUE), .groups = 'drop')
    row_totals <- MAT2 |>
      dplyr::group_by(.data$species, .data$modality_class) |>
      dplyr::summarise(row_total = sum(.data$BioSample, na.rm = TRUE), .groups = 'drop')
    col_totals <- MAT2 |>
      dplyr::group_by(.data$species, .data$.interaction_term) |>
      dplyr::summarise(col_total = sum(.data$BioSample, na.rm = TRUE), .groups = 'drop')
    MAT2 |>
      dplyr::left_join(total_by_species, by = 'species') |>
      dplyr::left_join(row_totals, by = c('species', 'modality_class')) |>
      dplyr::left_join(col_totals, by = c('species', '.interaction_term')) |>
      dplyr::mutate(
        expected = dplyr::if_else(
          as.numeric(.data$total_n) > 0,
          as.numeric(.data$row_total) *
            as.numeric(.data$col_total) / as.numeric(.data$total_n),
          NA_real_
        ),
        residual = dplyr::case_when(
          .data$expected > 0 ~
            (as.numeric(.data$BioSample) - .data$expected) / sqrt(.data$expected),
          .data$BioSample == 0L & .data$expected == 0 ~ 0,
          TRUE ~ NA_real_
        )
      ) |>
      dplyr::select(-.interaction_term, -total_n, -row_total, -col_total)
  }
  match_report <- function(anatomy, modality, species_use) {
    matched <- anatomy |>
      dplyr::filter(.data$species %in% .env$species_use) |>
      dplyr::left_join(modality, by = c('species', 'biosample_id')) |>
      dplyr::transmute(
        species = .data$species,
        biosample_id = .data$biosample_id,
        matched_to_sra = !is.na(.data$modality_class) & nzchar(.data$modality_class)
      ) |>
      dplyr::distinct(.data$species, .data$biosample_id, .data$matched_to_sra)
    report <- matched |>
      dplyr::group_by(.data$species) |>
      dplyr::summarise(
        operable = dplyr::n(),
        matched_to_sra = sum(.data$matched_to_sra, na.rm = TRUE),
        .groups = 'drop'
      ) |>
      dplyr::mutate(
        unmatched_to_sra = .data$operable - .data$matched_to_sra,
        matched_prop = dplyr::if_else(
          .data$operable > 0L,
          .data$matched_to_sra / .data$operable,
          NA_real_
        )
      )
    tibble::tibble(species = species_use) |>
      dplyr::left_join(report, by = 'species') |>
      dplyr::mutate(
        operable = as.integer(dplyr::coalesce(.data$operable, 0L)),
        matched_to_sra = as.integer(dplyr::coalesce(.data$matched_to_sra, 0L)),
        unmatched_to_sra = as.integer(dplyr::coalesce(.data$unmatched_to_sra, 0L)),
        matched_prop = dplyr::if_else(
          .data$operable > 0L,
          .data$matched_to_sra / .data$operable,
          NA_real_
        )
      )
  }
  class_levels <- .biosample_anatomy_profile_levels()
  modality_levels <- .biosample_modality_levels()
  modality <- .biosample_modality_profile(sra_prof)
  report <- match_report(anatomy, modality, species_use)
  for (i in seq_len(nrow(report))) {
    .gama_msg(
      report$species[[i]],
      ': ',
      report$matched_to_sra[[i]],
      ' matched to SRA out of ',
      report$operable[[i]],
      ' operable BioSample records.'
    )
  }
  if (level == 'anatomy_class') {
    counts <- anatomy |>
      dplyr::filter(.data$species %in% .env$species_use) |>
      dplyr::left_join(modality, by = c('species', 'biosample_id')) |>
      dplyr::filter(!is.na(.data$modality_class), nzchar(.data$modality_class)) |>
      dplyr::transmute(
        species = .data$species,
        biosample_id = .data$biosample_id,
        modality_class = .data$modality_class,
        anatomy_class = dplyr::coalesce(.data$anatomy_class, 'unknown')
      ) |>
      dplyr::mutate(
        anatomy_class = dplyr::if_else(
          .data$anatomy_class %in% .env$class_levels,
          .data$anatomy_class,
          'unknown'
        )
      ) |>
      dplyr::distinct(.data$species, .data$biosample_id, .data$modality_class, .data$anatomy_class) |>
      dplyr::count(.data$species, .data$modality_class, .data$anatomy_class, name = 'BioSample')
    MAT <- tidyr::expand_grid(
      species = species_use,
      modality_class = modality_levels,
      anatomy_class = class_levels
    ) |>
      dplyr::left_join(
        counts,
        by = c('species', 'modality_class', 'anatomy_class')
      ) |>
      dplyr::transmute(
        species = .data$species,
        modality_class = as.character(.data$modality_class),
        anatomy_class = as.character(.data$anatomy_class),
        BioSample = as.integer(dplyr::coalesce(.data$BioSample, 0L))
      )
    term_col <- 'anatomy_class'
    term_meta <- tibble::tibble(
      term_order = seq_along(class_levels),
      x_term = class_levels,
      anatomy_class = class_levels
    )
  } else {
    SUB <- .biosample_term_heatmap_core(
      BIO = BIO,
      SRA = SRA,
      species = species_use
    )
    if (!nrow(SUB)) .gama_stop('No anatomy subclasses found for selected species.')
    term_meta <- SUB |>
      dplyr::distinct(.data$term_order, .data$anatomy_subclass, .data$anatomy_class) |>
      dplyr::arrange(.data$term_order) |>
      dplyr::transmute(
        term_order = .data$term_order,
        x_term = .data$anatomy_subclass,
        anatomy_class = as.character(.data$anatomy_class)
      )
    MAT <- SUB |>
      dplyr::transmute(
        species = .data$species,
        modality_class = as.character(.data$modality_class),
        anatomy_subclass = as.character(.data$anatomy_subclass),
        BioSample = as.integer(.data$BioSample)
      )
    term_col <- 'anatomy_subclass'
  }
  OUT <- MAT |>
    add_residuals(term_col = term_col) |>
    dplyr::arrange(
      match(.data$species, .env$species_use),
      match(.data$modality_class, .env$modality_levels),
      match(.data[[term_col]], .env$term_meta$x_term)
    ) |>
    dplyr::select(
      species,
      class = modality_class,
      dplyr::all_of(term_col),
      BioSample,
      expected,
      residual
    )
  qi <- attr(BIO, 'query_info', exact = TRUE)
  if (is.null(qi)) qi <- attr(SRA, 'query_info', exact = TRUE)
  if (is.null(qi)) .gama_warn('No query_info found on input objects.')
  attr(OUT, 'query_info') <- qi
  attr(OUT, 'interaction_info') <- list(
    level = level,
    term_col = term_col,
    term_levels = term_meta$x_term,
    modality_levels = modality_levels,
    term_meta = term_meta,
    match_report = report
  )
  OUT <- .set_gama_object(OUT, 'summarise_interaction')
  class(OUT) <- unique(c('gdt_tbl', class(OUT)))
  OUT
}

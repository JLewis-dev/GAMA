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
#' @seealso [summarise_availability()], [summarise_sra_availability()],
#' [extract_assembly_metadata()], [extract_sra_metadata()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(
#'   c('Vigna angularis', 'Vigna reflexo-pilosa', 'Vigna vexillata'),
#'   synonyms = list('Vigna reflexo-pilosa' = 'Vigna glabrescens')
#' )
#' str(RESULTS, max.level = 2)
#' attr(RESULTS, 'query_info')
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

#' Summarise accession counts and compute data richness
#'
#' Collapses the output of [query_species()] into a species-level summary table
#' of accession counts (Assembly, SRA, BioSample) and data richness scores.
#'
#' @details
#' For each species, accession counts are extracted from the underlying query
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
#' \item `Assembly`, `SRA`, `BioSample`: accession counts per database
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

#' Summarise SRA modality composition
#'
#' Collapses experiment-level SRA metadata into species-level modality counts
#' using ontology-assigned classes and (optionally) subclasses.
#'
#' By default, the output includes class-level counts for the major modalities
#' (genomic, transcriptomic, epigenomic, chromatin, other, unknown), plus the
#' total number of SRA experiments per species.
#'
#' Profile cache (used by [summarise_sra_skew()]):
#' In addition to the summary table, this function attaches a cached, UID-level
#' profile as an attribute `sra_profile`. The profile contains (at minimum)
#' `species`, `entrez_uid`, `biosample`, `bioproject`, `class`, `subclass`, and
#' GEO linkage fields (`geo_linked`, `gse_ids`, `gsm_ids`). This cache is
#' intended to be re-used locally for downstream summaries that require
#' within-species structure (e.g. replication skew across BioProjects or
#' BioSamples) without re-querying NCBI. In particular, [summarise_sra_skew()]
#' consumes `attr(x, 'sra_profile')` from the output of this function.
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

#' Extract filtered Assembly metadata
#'
#' Retrieves and structures assembly metadata for one or more species using the
#' Assembly identifiers stored in the output of [query_species()]. Metadata are
#' returned in a tidy tibble and optionally reduced to a single 'best' assembly
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
#' @seealso [query_species()], [summarise_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' ASM <- extract_assembly_metadata(RESULTS, best = TRUE)
#' print(ASM)
#' }
#' @export
extract_assembly_metadata <- function(results, species = NULL, best = FALSE) {
  results <- .gama_require_output(results, 'query_species')
  if (!is.null(species)) {
    species <- as.character(species)
    species <- species[!is.na(species) & nzchar(species)]
    missing <- setdiff(species, names(results))
    if (length(missing) > 0L) .gama_stop('Requested species not found in `results`: ', paste(missing, collapse = ', '), '.')
    results <- results[species]
  }
  n  <- length(results)
  pb <- .pb_init(n)
  META <- lapply(seq_along(results), function(i) {
    sp  <- names(results)[i]
    res <- results[[i]]
    .pb_tick(pb, i)
    asm <- res$assembly
    if (is.null(asm) || (asm$count %||% 0L) == 0L) {
      return(tibble::tibble(
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
      ))
    }
    SUMS <- .fetch_search_summaries('assembly', asm, batch_size = 100)
    if (!length(SUMS)) {
      return(tibble::tibble(
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
      ))
    }
    do.call(dplyr::bind_rows, lapply(seq_along(SUMS), function(j) {
      x <- SUMS[[j]]
      acc <- names(SUMS)[j]
      if (is.na(acc) || !nzchar(acc)) acc <- .esummary_uid(x)
      tibble::tibble(
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
    }))
  })
  .pb_close(pb)
  META <- dplyr::bind_rows(META)
  if (!best) return(.as_gdt_table(META, results, 'extract_assembly_metadata'))
  BEST <- lapply(split(META, META$species), function(x) {
    STRUCT <- .ASSEMBLY_WEIGHTS[x$level]
    STRUCT[is.na(STRUCT)] <- 0
    N50 <- x$n50
    max_struct <- max(STRUCT)
    tied_idx   <- which(STRUCT == max_struct)
    best_idx <- if (length(tied_idx) > 1) {
      tied_idx[which.max(N50[tied_idx])]
    } else tied_idx
    x[best_idx, , drop = FALSE]
  })
  OUT <- dplyr::bind_rows(BEST)
  .as_gdt_table(OUT, results, 'extract_assembly_metadata')
}

#' Extract filtered SRA metadata
#'
#' Retrieves experiment-level SRA metadata for one or more species, then
#' assigns curated modality classes and subclasses together with GEO linkage
#' fields.
#'
#' Optional filters restrict the output to particular modality classes,
#' subclasses, or GEO-linked experiments only.
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
#' `CUT&Tag`, `ChIP-seq`, `DNase-seq`, `FAIRE-seq`, `MNase-seq`,
#' `SELEX`}
#' \item{`chromatin`}{`3C-based`, `ChIA-PET`, `Hi-C`, `TCC`}
#' \item{`other`}{`other`}
#' }
#'
#' @param results A list returned by [query_species()], containing SRA IDs.
#' @param species `NULL` (default) to include all species, or a character
#' vector specifying which species to include.
#' @param class Optional character vector of modality classes to retain.
#' @param subclass Optional character vector of modality subclasses to retain.
#' @param only_geo Logical; if `TRUE`, retain only GEO-linked experiments.
#'
#' @return A tibble with one row per SRA experiment, including identifiers,
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
  if (!is.null(species)) {
    species <- as.character(species)
    species <- species[!is.na(species) & nzchar(species)]
    missing <- setdiff(species, names(results))
    if (length(missing) > 0L) .gama_warn('Requested species not found in `results`: ', paste(missing, collapse = ', '), '. Dropping.')
    species <- species[species %in% names(results)]
  }
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

#' Summarise SRA replication skew across BioProjects / BioSamples
#'
#' Quantifies replication skew across independent units (BioProject or
#' BioSample) using the cached UID-level SRA profile produced by
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
#' is the proportion of experiments in each BioProject/BioSample. Larger values
#' indicate a more even spread; values near 1 indicate strong concentration in
#' few units.
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
#' `class`, `min`, `q25`, `med`, `q75`, `max` (experiments per unit), and `eff`
#' (effective number of units; inverse Simpson index).
#'
#' @seealso [summarise_sra_availability()], [plot_sra_skew()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' SKEW <- summarise_sra_skew(SRA_SUMMARY, class = 'transcriptomic')
#' print(SKEW)
#' }
#'
#' @export
summarise_sra_skew <- function(x, species = NULL, unit = c('bioproject', 'biosample'), class = NULL) {
  x <- .gama_require_output(x, 'summarise_sra_availability')
  if (missing(unit)) unit <- 'bioproject'
  if (length(unit) != 1L) .gama_stop('`unit` must be a single value: \'bioproject\' or \'biosample\'.')
  unit <- match.arg(unit)
  if (!is.null(class) && length(class) != 1L) .gama_stop('`class` must be a single modality class (or NULL). Use one class per call.')
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
    return(out)
  }
  prof_use <- prof[prof$species %in% species_all, , drop = FALSE]
  if (!is.null(class)) prof_use <- prof_use[prof_use$class %in% class, , drop = FALSE]
  calc_one <- function(sp) {
    units <- prof_use[[unit_col]][prof_use$species == sp]
    units <- units[!is.na(units)]
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
  attr(out, 'sra_profile') <- prof_use
  attr(out, 'sra_profile_info') <- attr(x, 'sra_profile_info', exact = TRUE)
  out
}

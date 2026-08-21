# DIAGNOSTICS =================================================================

.diagnostics_rate <- function(numerator, denominator) {
  dplyr::if_else(
    denominator > 0L,
    as.numeric(numerator) / as.numeric(denominator),
    NA_real_
  )
}

.diagnostics_object <- function(x) {
  supported <- c(
    'summarise_sra_availability',
    'summarise_biosample_availability',
    'summarise_interaction',
    'summarise_sra_skew',
    'summarise_biosample_skew'
  )
  object <- .detect_gama_object(x)
  if (!object %in% supported) {
    .gama_stop(
      '`x` must be the output of one of: ',
      paste0(supported, '()', collapse = ', '),
      '; detected ', object, ' object.'
    )
  }
  object
}

.diagnostics_interaction_report <- function(x) {
  info <- .gama_require_cache(
    x,
    attr_name = 'interaction_info',
    required_cols = 'match_report',
    source = 'summarise_interaction'
  )
  report <- info$match_report
  required <- c('species', 'BioSample', 'operable', 'linked', 'unknown')
  if (!is.data.frame(report)) {
    .gama_input_error(
      'summarise_interaction',
      detected = .detect_gama_object(x),
      detail = 'cache \'interaction_info$match_report\' must be a data frame.'
    )
  }
  missing <- setdiff(required, names(report))
  if (length(missing)) {
    .gama_input_error(
      'summarise_interaction',
      detected = .detect_gama_object(x),
      detail = paste0(
        'cache \'interaction_info$match_report\' missing required columns: ',
        paste(missing, collapse = ', '),
        '.'
      )
    )
  }
  tibble::as_tibble(report)
}

.diagnostics_sra_skew_unit <- function(x) {
  units <- intersect(c('BioProject', 'BioSample'), names(x))
  if (length(units) != 1L) {
    .gama_input_error(
      'summarise_sra_skew',
      detected = .detect_gama_object(x),
      detail = paste0(
        'must contain exactly one active unit column: ',
        'BioProject or BioSample.'
      )
    )
  }
  if (identical(units, 'BioProject')) 'bioproject' else 'biosample'
}

.diagnostics_recoverability <- function(x, object) {
  switch(object,
    summarise_sra_availability = {
      data <- .gama_require_output(
        x,
        object,
        required_cols = c('species', 'SRA', 'unknown')
      )
      tibble::tibble(
        species = as.character(data$species),
        SRA = as.integer(data$SRA),
        unknown = as.integer(data$unknown),
        recoverability = .diagnostics_rate(
          data$SRA - data$unknown,
          data$SRA
        )
      )
    },
    summarise_biosample_availability = {
      data <- .gama_require_output(
        x,
        object,
        required_cols = c('species', 'BioSample', 'operable', 'unknown')
      )
      tibble::tibble(
        species = as.character(data$species),
        BioSample = as.integer(data$BioSample),
        operable = as.integer(data$operable),
        unknown = as.integer(data$unknown),
        recoverability = .diagnostics_rate(
          data$operable - data$unknown,
          data$BioSample
        )
      )
    },
    summarise_interaction = {
      data <- .diagnostics_interaction_report(x)
      tibble::tibble(
        species = as.character(data$species),
        BioSample = as.integer(data$BioSample),
        operable = as.integer(data$operable),
        linked = as.integer(data$linked),
        unknown = as.integer(data$unknown),
        recoverability = .diagnostics_rate(
          data$linked - data$unknown,
          data$BioSample
        )
      )
    },
    summarise_sra_skew = {
      recovery <- .gama_require_cache(
        x,
        attr_name = 'skew_id_recovery',
        required_cols = c('species', 'records', 'included'),
        source = object
      )
      tibble::tibble(
        species = as.character(recovery$species),
        SRA = as.integer(recovery$records),
        linked = as.integer(recovery$included),
        recoverability = .diagnostics_rate(
          recovery$included,
          recovery$records
        )
      )
    },
    summarise_biosample_skew = {
      recovery <- .gama_require_cache(
        x,
        attr_name = 'skew_id_recovery',
        required_cols = c('species', 'records', 'included'),
        source = object
      )
      tibble::tibble(
        species = as.character(recovery$species),
        BioSample = as.integer(recovery$records),
        linked = as.integer(recovery$included),
        recoverability = .diagnostics_rate(
          recovery$included,
          recovery$records
        )
      )
    }
  )
}

.diagnostics_records <- function(x, object) {
  switch(object,
    summarise_sra_availability = .gama_require_cache(
      x,
      attr_name = 'sra_profile',
      required_cols = c(
        'species',
        'entrez_uid',
        'biosample',
        'bioproject',
        'strategy_raw',
        'strategy_norm',
        'class'
      ),
      source = object
    ) |>
      dplyr::filter(.data$class == 'unknown') |>
      dplyr::transmute(
        species = as.character(.data$species),
        entrez_uid = as.character(.data$entrez_uid),
        biosample = as.character(.data$biosample),
        bioproject = as.character(.data$bioproject),
        strategy_raw = as.character(.data$strategy_raw),
        strategy_norm = as.character(.data$strategy_norm),
        class = as.character(.data$class)
      ),
    summarise_biosample_availability = .gama_require_cache(
      x,
      attr_name = 'biosample_canonical_profile',
      required_cols = c(
        'species',
        'entrez_uid',
        'biosample_id',
        'bioproject',
        'value_raw',
        'value_norm',
        'anatomy_class'
      ),
      source = object
    ) |>
      dplyr::filter(.data$anatomy_class == 'unknown') |>
      dplyr::transmute(
        species = as.character(.data$species),
        entrez_uid = as.character(.data$entrez_uid),
        biosample = as.character(.data$biosample_id),
        bioproject = as.character(.data$bioproject),
        value_raw = as.character(.data$value_raw),
        value_norm = as.character(.data$value_norm),
        anatomy_class = as.character(.data$anatomy_class)
      ) |>
      dplyr::distinct(),
    summarise_interaction = .gama_require_cache(
      x,
      attr_name = 'interaction_profile',
      required_cols = c(
        'species',
        'entrez_uid',
        'biosample',
        'bioproject',
        'strategy_raw',
        'strategy_norm',
        'class',
        'value_raw',
        'value_norm',
        'anatomy_class'
      ),
      source = object
    ) |>
      dplyr::filter(.data$class == 'unknown' | .data$anatomy_class == 'unknown') |>
      dplyr::transmute(
        species = as.character(.data$species),
        entrez_uid = as.character(.data$entrez_uid),
        biosample = as.character(.data$biosample),
        bioproject = as.character(.data$bioproject),
        strategy_raw = as.character(.data$strategy_raw),
        strategy_norm = as.character(.data$strategy_norm),
        class = as.character(.data$class),
        value_raw = as.character(.data$value_raw),
        value_norm = as.character(.data$value_norm),
        anatomy_class = as.character(.data$anatomy_class)
      ),
    summarise_sra_skew = {
      unit <- .diagnostics_sra_skew_unit(x)
      .gama_require_cache(
        x,
        attr_name = 'sra_profile',
        required_cols = c(
          'species',
          'entrez_uid',
          'biosample',
          'bioproject'
        ),
        source = object
      ) |>
        dplyr::filter(
          is.na(.data[[unit]]) |
            !nzchar(trimws(as.character(.data[[unit]])))
        ) |>
        dplyr::transmute(
          species = as.character(.data$species),
          entrez_uid = as.character(.data$entrez_uid),
          biosample = as.character(.data$biosample),
          bioproject = as.character(.data$bioproject)
        )
    },
    summarise_biosample_skew = .gama_require_cache(
      x,
      attr_name = 'biosample_anatomy_profile',
      required_cols = c(
        'species',
        'entrez_uid',
        'biosample_id',
        'bioproject'
      ),
      source = object
    ) |>
      dplyr::filter(
        is.na(.data$bioproject) |
          !nzchar(trimws(as.character(.data$bioproject)))
      ) |>
      dplyr::transmute(
        species = as.character(.data$species),
        entrez_uid = as.character(.data$entrez_uid),
        biosample = as.character(.data$biosample_id),
        bioproject = as.character(.data$bioproject)
      )
  )
}

.diagnostics_empty_records <- function(x, object) {
  text <- switch(object,
    summarise_sra_availability = 'No unknown SRA records found; returning empty table.',
    summarise_biosample_availability = paste0(
      'No operable BioSample records with unknown anatomy found; ',
      'returning empty table.'
    ),
    summarise_interaction = paste0(
      'No linked records with unknown SRA modality or BioSample anatomy found; ',
      'returning empty table.'
    ),
    summarise_sra_skew = paste0(
      'No SRA records with missing ',
      if (identical(.diagnostics_sra_skew_unit(x), 'bioproject')) {
        'BioProject'
      } else {
        'BioSample'
      },
      ' IDs found; returning empty table.'
    ),
    summarise_biosample_skew = paste0(
      'No operable BioSample records with missing BioProject IDs found; ',
      'returning empty table.'
    )
  )
  .gama_msg(text)
}

#' Report GAMA diagnostics
#'
#' Reports metadata recoverability and identifies records with unresolved
#' classifications or missing linkage identifiers in supported GAMA objects.
#'
#' @details
#' `report_diagnostics()` accepts output from [summarise_sra_availability()],
#' [summarise_sra_skew()], [summarise_biosample_availability()],
#' [summarise_biosample_skew()], or [summarise_interaction()]. Diagnostic data
#' are derived from summary fields, record-level profiles, and recovery data.
#'
#' `recoverability` is the proportion of the relevant record denominator for
#' which the required classification or linkage information was recovered.
#' Values range from 0 to 1 and are calculated from the returned count columns:
#' \itemize{
#' \item SRA availability: `(SRA - unknown) / SRA`
#' \item SRA skew: `linked / SRA`
#' \item BioSample availability: `(operable - unknown) / BioSample`
#' \item BioSample skew: `linked / BioSample`
#' \item interaction: `(linked - unknown) / BioSample`
#' }
#' `unknown` counts SRA records whose modality metadata remain missing or
#' uninformative after fallback parsing and BioSample records with missing-like
#' sample-source values or no anatomy ontology match. Non-missing, interpretable
#' SRA strategies outside the modality ontology are classified as `other`. In
#' interaction diagnostics, each linked BioSample is counted once when either
#' or both classifications are unknown. For skew diagnostics, `linked` counts
#' records carrying the active identifier. A zero denominator produces
#' `NA_real_`.
#'
#' With `view = 'recoverability'`, the function returns one row per species.
#' With `view = 'records'`, it returns unknown SRA modalities, operable
#' BioSample records with unknown anatomy, linked interaction records with
#' unknown modality or anatomy, or skew records missing the active linkage
#' identifier. If none are found, an empty table retaining the complete typed
#' column schema is returned and a message is emitted.
#'
#' @param x A data frame or tibble returned by [summarise_sra_availability()],
#' [summarise_sra_skew()], [summarise_biosample_availability()],
#' [summarise_biosample_skew()], or [summarise_interaction()].
#' @param view Diagnostic view: either `recoverability` (default) or `records`.
#'
#' @return A tibble with class `gdt_tbl` carrying the input `query_info`
#' attribute. The recoverability view contains one row per species, the
#' relevant denominator, diagnostic counts, and `recoverability`; the
#' records view contains object-specific identifiers and diagnostic fields.
#'
#' @seealso [summarise_sra_availability()], [summarise_sra_skew()],
#' [summarise_biosample_availability()], [summarise_biosample_skew()],
#' [summarise_interaction()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' BIO_SUMMARY <- summarise_biosample_availability(RESULTS)
#' DIAGNOSTICS <- report_diagnostics(BIO_SUMMARY)
#' print(DIAGNOSTICS)
#' }
#'
#' @export
report_diagnostics <- function(x, view = 'recoverability') {
  view <- .gama_validate_parameters(
    view,
    'view',
    c('recoverability', 'records'),
    multiple = FALSE,
    allow_null = FALSE
  )
  object <- .diagnostics_object(x)
  output <- switch(
    view,
    recoverability = .diagnostics_recoverability(x, object),
    records = .diagnostics_records(x, object)
  )
  if (identical(view, 'records') && !nrow(output)) {
    .diagnostics_empty_records(x, object)
  }
  .as_gdt_table(output, x)
}

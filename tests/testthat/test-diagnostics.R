diagnostics_cases <- list(
  sra_summary = list(
    input = 'SRA_SUMMARY_Arabidopsis_thaliana',
    recoverability =
      'SRA_SUMMARY_DIAGNOSTICS_RECOVERABILITY_Arabidopsis_thaliana',
    records = 'SRA_SUMMARY_DIAGNOSTICS_RECORDS_Arabidopsis_thaliana',
    recoverability_cols = c(
      'species', 'SRA', 'unknown', 'recoverability'
    ),
    records_cols = c(
      'species', 'entrez_uid', 'biosample', 'bioproject',
      'strategy_raw', 'strategy_norm', 'class'
    )
  ),
  sra_skew_bioproject = list(
    input = 'SRA_SKEW_BP_Arabidopsis_thaliana',
    recoverability = paste0(
      'SRA_SKEW_BP_DIAGNOSTICS_RECOVERABILITY_',
      'Arabidopsis_thaliana'
    ),
    records = paste0(
      'SRA_SKEW_BP_DIAGNOSTICS_RECORDS_',
      'Arabidopsis_thaliana'
    ),
    recoverability_cols = c(
      'species', 'SRA', 'linked', 'recoverability'
    ),
    records_cols = c(
      'species', 'entrez_uid', 'biosample', 'bioproject'
    )
  ),
  sra_skew_biosample = list(
    input = 'SRA_SKEW_BS_Arabidopsis_thaliana',
    recoverability = paste0(
      'SRA_SKEW_BS_DIAGNOSTICS_RECOVERABILITY_',
      'Arabidopsis_thaliana'
    ),
    records = paste0(
      'SRA_SKEW_BS_DIAGNOSTICS_RECORDS_',
      'Arabidopsis_thaliana'
    ),
    recoverability_cols = c(
      'species', 'SRA', 'linked', 'recoverability'
    ),
    records_cols = c(
      'species', 'entrez_uid', 'biosample', 'bioproject'
    )
  ),
  bio_summary = list(
    input = 'BIO_SUMMARY_Arabidopsis_thaliana',
    recoverability =
      'BIO_SUMMARY_DIAGNOSTICS_RECOVERABILITY_Arabidopsis_thaliana',
    records = 'BIO_SUMMARY_DIAGNOSTICS_RECORDS_Arabidopsis_thaliana',
    recoverability_cols = c(
      'species', 'BioSample', 'operable', 'unknown', 'recoverability'
    ),
    records_cols = c(
      'species', 'entrez_uid', 'biosample', 'bioproject',
      'value_raw', 'value_norm', 'anatomy_class'
    )
  ),
  bio_skew = list(
    input = 'BIO_SKEW_Arabidopsis_thaliana',
    recoverability =
      'BIO_SKEW_DIAGNOSTICS_RECOVERABILITY_Arabidopsis_thaliana',
    records = 'BIO_SKEW_DIAGNOSTICS_RECORDS_Arabidopsis_thaliana',
    recoverability_cols = c(
      'species', 'BioSample', 'linked', 'recoverability'
    ),
    records_cols = c(
      'species', 'entrez_uid', 'biosample', 'bioproject'
    )
  ),
  interaction_class = list(
    input = 'INTERACTION_CLASS_Arabidopsis_thaliana',
    recoverability =
      'INTERACTION_DIAGNOSTICS_RECOVERABILITY_Arabidopsis_thaliana',
    records = 'INTERACTION_DIAGNOSTICS_RECORDS_Arabidopsis_thaliana',
    recoverability_cols = c(
      'species', 'BioSample', 'operable', 'linked', 'unknown',
      'recoverability'
    ),
    records_cols = c(
      'species', 'entrez_uid', 'biosample', 'bioproject',
      'strategy_raw', 'strategy_norm', 'class', 'value_raw',
      'value_norm', 'anatomy_class'
    )
  ),
  interaction_subclass = list(
    input = 'INTERACTION_SUBCLASS_Arabidopsis_thaliana',
    recoverability =
      'INTERACTION_DIAGNOSTICS_RECOVERABILITY_Arabidopsis_thaliana',
    records = 'INTERACTION_DIAGNOSTICS_RECORDS_Arabidopsis_thaliana',
    recoverability_cols = c(
      'species', 'BioSample', 'operable', 'linked', 'unknown',
      'recoverability'
    ),
    records_cols = c(
      'species', 'entrez_uid', 'biosample', 'bioproject',
      'strategy_raw', 'strategy_norm', 'class', 'value_raw',
      'value_norm', 'anatomy_class'
    )
  )
)

expected_diagnostics_rate <- function(numerator, denominator) {
  out <- rep(NA_real_, length(denominator))
  has_denominator <- denominator > 0L
  out[has_denominator] <-
    as.numeric(numerator[has_denominator]) /
    as.numeric(denominator[has_denominator])
  out
}

expect_diagnostic_data <- function(x, expected) {
  testthat::expect_named(x, names(expected))
  testthat::expect_equal(nrow(x), nrow(expected))
  for (col in names(expected)) {
    testthat::expect_identical(x[[col]], expected[[col]])
  }
}

expect_empty_diagnostic_records <- function(x, regexp) {
  OUT <- NULL
  testthat::expect_message(
    OUT <- report_diagnostics(x, view = 'records'),
    regexp
  )
  expect_gdt_tbl(OUT)
  testthat::expect_equal(nrow(OUT), 0L)
}

test_that('Diagnostics fixtures are valid report_diagnostics outputs', {
  fixture_names <- unique(unlist(lapply(diagnostics_cases, function(case) {
    c(case$recoverability, case$records)
  }), use.names = FALSE))
  for (fixture in fixture_names) {
    expect_gdt_tbl(load_fixture(fixture))
  }
})

test_that('Diagnostics fixtures contain expected columns', {
  for (case in diagnostics_cases) {
    RECOVERABILITY <- load_fixture(case$recoverability)
    RECORDS <- load_fixture(case$records)
    expect_named(RECOVERABILITY, case$recoverability_cols)
    expect_named(RECORDS, case$records_cols)
  }
})

test_that('Diagnostics fixtures preserve query provenance', {
  for (case in diagnostics_cases) {
    INPUT <- load_fixture(case$input)
    for (view in c('recoverability', 'records')) {
      OUTPUT <- load_fixture(case[[view]])
      expect_identical(
        attr(OUTPUT, 'query_info', exact = TRUE),
        attr(INPUT, 'query_info', exact = TRUE)
      )
    }
  }
})

test_that('report_diagnostics returns expected recoverability summaries', {
  for (case_name in names(diagnostics_cases)) {
    case <- diagnostics_cases[[case_name]]
    INPUT <- load_fixture(case$input)
    EXPECTED <- load_fixture(case$recoverability)
    OBSERVED <- if (identical(case_name, 'sra_summary')) {
      report_diagnostics(INPUT)
    } else {
      report_diagnostics(INPUT, view = 'recoverability')
    }
    expect_equal(OBSERVED, EXPECTED)
  }
})

test_that('report_diagnostics returns expected diagnostic records', {
  for (case in diagnostics_cases) {
    INPUT <- load_fixture(case$input)
    EXPECTED <- load_fixture(case$records)
    OBSERVED <- suppressMessages(
      report_diagnostics(INPUT, view = 'records')
    )
    expect_equal(OBSERVED, EXPECTED)
  }
})

test_that('report_diagnostics applies correct recoverability formulas', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  SRA_DIAGNOSTICS <- report_diagnostics(SRA_SUMMARY)
  expect_equal(
    SRA_DIAGNOSTICS$recoverability,
    expected_diagnostics_rate(
      SRA_SUMMARY$SRA - SRA_SUMMARY$unknown,
      SRA_SUMMARY$SRA
    )
  )
  for (fixture in c(
    'SRA_SKEW_BP_Arabidopsis_thaliana',
    'SRA_SKEW_BS_Arabidopsis_thaliana'
  )) {
    SRA_SKEW <- load_fixture(fixture)
    RECOVERY <- attr(SRA_SKEW, 'skew_id_recovery', exact = TRUE)
    SRA_SKEW_DIAGNOSTICS <- report_diagnostics(SRA_SKEW)
    expect_equal(
      SRA_SKEW_DIAGNOSTICS$recoverability,
      expected_diagnostics_rate(RECOVERY$included, RECOVERY$records)
    )
  }
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  BIO_DIAGNOSTICS <- report_diagnostics(BIO_SUMMARY)
  expect_equal(
    BIO_DIAGNOSTICS$recoverability,
    expected_diagnostics_rate(
      BIO_SUMMARY$operable - BIO_SUMMARY$unknown,
      BIO_SUMMARY$BioSample
    )
  )
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  RECOVERY <- attr(BIO_SKEW, 'skew_id_recovery', exact = TRUE)
  BIO_SKEW_DIAGNOSTICS <- report_diagnostics(BIO_SKEW)
  expect_equal(
    BIO_SKEW_DIAGNOSTICS$recoverability,
    expected_diagnostics_rate(RECOVERY$included, RECOVERY$records)
  )
  INTERACTION <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  MATCH_REPORT <- attr(
    INTERACTION,
    'interaction_info',
    exact = TRUE
  )$match_report
  INTERACTION_DIAGNOSTICS <- report_diagnostics(INTERACTION)
  expect_equal(
    INTERACTION_DIAGNOSTICS$recoverability,
    expected_diagnostics_rate(
      MATCH_REPORT$linked - MATCH_REPORT$unknown,
      MATCH_REPORT$BioSample
    )
  )
})

test_that('report_diagnostics follows correct record-selection logic', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  SRA_EXPECTED <- attr(SRA_SUMMARY, 'sra_profile', exact = TRUE) |>
    dplyr::filter(.data$class == 'unknown') |>
    dplyr::transmute(
      species = as.character(.data$species),
      entrez_uid = as.character(.data$entrez_uid),
      biosample = as.character(.data$biosample),
      bioproject = as.character(.data$bioproject),
      strategy_raw = as.character(.data$strategy_raw),
      strategy_norm = as.character(.data$strategy_norm),
      class = as.character(.data$class)
    )
  SRA_OBSERVED <- report_diagnostics(SRA_SUMMARY, view = 'records')
  expect_diagnostic_data(SRA_OBSERVED, SRA_EXPECTED)
  for (fixture in c(
    'SRA_SKEW_BP_Arabidopsis_thaliana',
    'SRA_SKEW_BS_Arabidopsis_thaliana'
  )) {
    SRA_SKEW <- load_fixture(fixture)
    unit_col <- if ('BioProject' %in% names(SRA_SKEW)) {
      'bioproject'
    } else {
      'biosample'
    }
    SRA_SKEW_EXPECTED <- attr(
      SRA_SKEW,
      'sra_profile',
      exact = TRUE
    ) |>
      dplyr::filter(
        is.na(.data[[unit_col]]) |
          !nzchar(trimws(as.character(.data[[unit_col]])))
      ) |>
      dplyr::transmute(
        species = as.character(.data$species),
        entrez_uid = as.character(.data$entrez_uid),
        biosample = as.character(.data$biosample),
        bioproject = as.character(.data$bioproject)
      )
    SRA_SKEW_OBSERVED <- suppressMessages(
      report_diagnostics(SRA_SKEW, view = 'records')
    )
    expect_diagnostic_data(SRA_SKEW_OBSERVED, SRA_SKEW_EXPECTED)
  }
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  BIO_EXPECTED <- attr(
    BIO_SUMMARY,
    'biosample_canonical_profile',
    exact = TRUE
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
    dplyr::distinct()
  BIO_OBSERVED <- report_diagnostics(BIO_SUMMARY, view = 'records')
  expect_diagnostic_data(BIO_OBSERVED, BIO_EXPECTED)
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  BIO_SKEW_EXPECTED <- attr(
    BIO_SKEW,
    'biosample_anatomy_profile',
    exact = TRUE
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
  BIO_SKEW_OBSERVED <- report_diagnostics(BIO_SKEW, view = 'records')
  expect_diagnostic_data(BIO_SKEW_OBSERVED, BIO_SKEW_EXPECTED)
  for (fixture in c(
    'INTERACTION_CLASS_Arabidopsis_thaliana',
    'INTERACTION_SUBCLASS_Arabidopsis_thaliana'
  )) {
    INTERACTION <- load_fixture(fixture)
    INTERACTION_EXPECTED <- attr(
      INTERACTION,
      'interaction_profile',
      exact = TRUE
    ) |>
      dplyr::filter(
        is.na(.data$class) |
          .data$class == 'unknown' |
          .data$anatomy_class == 'unknown'
      ) |>
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
      )
    INTERACTION_OBSERVED <- report_diagnostics(
      INTERACTION,
      view = 'records'
    )
    expect_diagnostic_data(
      INTERACTION_OBSERVED,
      INTERACTION_EXPECTED
    )
  }
})

test_that('report_diagnostics handles zero recoverability denominators', {
  diagnostics_rate <- getFromNamespace('.diagnostics_rate', 'GAMA')
  expect_equal(
    diagnostics_rate(c(3L, 0L), c(4L, 0L)),
    c(0.75, NA_real_)
  )
})

test_that('report_diagnostics reports empty record views', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  PROFILE <- attr(SRA_SUMMARY, 'sra_profile', exact = TRUE)
  PROFILE$class <- rep('genomic', nrow(PROFILE))
  attr(SRA_SUMMARY, 'sra_profile') <- PROFILE
  expect_empty_diagnostic_records(
    SRA_SUMMARY,
    'No unknown SRA records found; returning empty table\\.'
  )
  SRA_SKEW <- load_fixture(
    'SRA_SKEW_BP_Arabidopsis_thaliana'
  )
  PROFILE <- attr(SRA_SKEW, 'sra_profile', exact = TRUE)
  ids <- trimws(as.character(PROFILE$bioproject))
  ids[is.na(ids) | !nzchar(ids)] <- 'PRJ_TEST'
  PROFILE$bioproject <- ids
  attr(SRA_SKEW, 'sra_profile') <- PROFILE
  expect_empty_diagnostic_records(
    SRA_SKEW,
    paste0(
      'No SRA records with missing BioProject IDs found; ',
      'returning empty table\\.'
    )
  )
  SRA_SKEW <- load_fixture(
    'SRA_SKEW_BS_Arabidopsis_thaliana'
  )
  PROFILE <- attr(SRA_SKEW, 'sra_profile', exact = TRUE)
  ids <- trimws(as.character(PROFILE$biosample))
  ids[is.na(ids) | !nzchar(ids)] <- 'SAMN_TEST'
  PROFILE$biosample <- ids
  attr(SRA_SKEW, 'sra_profile') <- PROFILE
  expect_empty_diagnostic_records(
    SRA_SKEW,
    paste0(
      'No SRA records with missing BioSample IDs found; ',
      'returning empty table\\.'
    )
  )
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  PROFILE <- attr(
    BIO_SUMMARY,
    'biosample_canonical_profile',
    exact = TRUE
  )
  PROFILE$anatomy_class <- rep('aerial', nrow(PROFILE))
  attr(BIO_SUMMARY, 'biosample_canonical_profile') <- PROFILE
  expect_empty_diagnostic_records(
    BIO_SUMMARY,
    paste0(
      'No operable BioSample records with unknown anatomy found; ',
      'returning empty table\\.'
    )
  )
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  PROFILE <- attr(
    BIO_SKEW,
    'biosample_anatomy_profile',
    exact = TRUE
  )
  ids <- trimws(as.character(PROFILE$bioproject))
  ids[is.na(ids) | !nzchar(ids)] <- 'PRJ_TEST'
  PROFILE$bioproject <- ids
  attr(BIO_SKEW, 'biosample_anatomy_profile') <- PROFILE
  expect_empty_diagnostic_records(
    BIO_SKEW,
    paste0(
      'No operable BioSample records with missing BioProject IDs found; ',
      'returning empty table\\.'
    )
  )
  INTERACTION <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  PROFILE <- attr(INTERACTION, 'interaction_profile', exact = TRUE)
  PROFILE$class <- rep('genomic', nrow(PROFILE))
  PROFILE$anatomy_class <- rep('aerial', nrow(PROFILE))
  attr(INTERACTION, 'interaction_profile') <- PROFILE
  expect_empty_diagnostic_records(
    INTERACTION,
    paste0(
      'No unlinked BioSample records or linked records with unknown ',
      'SRA modality or BioSample anatomy found; returning empty table\\.'
    )
  )
})

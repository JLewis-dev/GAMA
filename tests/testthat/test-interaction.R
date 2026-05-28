expected_interaction_values <- function(INTERACTION) {
  term_col <- if ('anatomy_class' %in% names(INTERACTION)) 'anatomy_class' else 'anatomy_subclass'
  INTERACTION |>
    tibble::as_tibble() |>
    dplyr::mutate(.interaction_term = as.character(.data[[term_col]])) |>
    dplyr::group_by(.data$species) |>
    dplyr::mutate(total_n = sum(.data$BioSample, na.rm = TRUE)) |>
    dplyr::group_by(.data$species, .data$class) |>
    dplyr::mutate(row_total = sum(.data$BioSample, na.rm = TRUE)) |>
    dplyr::group_by(.data$species, .data$.interaction_term) |>
    dplyr::mutate(col_total = sum(.data$BioSample, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      expected_calc = dplyr::if_else(
        as.numeric(.data$total_n) > 0,
        as.numeric(.data$row_total) * as.numeric(.data$col_total) / as.numeric(.data$total_n),
        NA_real_
      ),
      residual_calc = dplyr::case_when(
        .data$expected_calc > 0 ~ (as.numeric(.data$BioSample) - .data$expected_calc) / sqrt(.data$expected_calc),
        .data$BioSample == 0L & .data$expected_calc == 0 ~ 0,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::select(expected_calc, residual_calc)
}

synthetic_interaction_inputs <- function() {
  query_info <- list(
    tool_version = '0.3.0',
    query_time_utc = '2026-05-22T11:23:11Z',
    databases = c('assembly', 'sra', 'biosample'),
    terms = list('Synthetic species' = 'Synthetic species[Organism]'),
    synonyms = list('Synthetic species' = NULL)
  )
  SRA_PROFILE <- tibble::tibble(
    species = rep('Synthetic species', 6),
    entrez_uid = paste0('SRX', seq_len(6)),
    biosample = paste0('SAM', seq_len(6)),
    bioproject = rep('PRJ1', 6),
    class = c('genomic', 'genomic', 'genomic', 'transcriptomic', 'transcriptomic', 'transcriptomic'),
    subclass = c('WGS', 'WGS', 'WGS', 'RNA-seq', 'RNA-seq', 'RNA-seq'),
    geo_linked = rep(FALSE, 6),
    gse_ids = rep(NA_character_, 6),
    gsm_ids = rep(NA_character_, 6)
  )
  BIO_ANATOMY <- tibble::tibble(
    species = rep('Synthetic species', 6),
    biosample_id = paste0('SAM', seq_len(6)),
    anatomy_class = c('aerial', 'aerial', 'ground', 'aerial', 'ground', 'ground'),
    anatomy_subclass = c('leaf', 'leaf', 'root', 'leaf', 'root', 'root')
  )
  BIO_CANONICAL <- tibble::tibble(
    species = BIO_ANATOMY$species,
    biosample_id = BIO_ANATOMY$biosample_id,
    anatomy_term = BIO_ANATOMY$anatomy_subclass,
    anatomy_class = BIO_ANATOMY$anatomy_class,
    anatomy_subclass = BIO_ANATOMY$anatomy_subclass,
    rank = NA_character_,
    ontology_namespace = NA_character_,
    ontology_id = NA_character_,
    ontology_label = NA_character_
  )
  SRA_SUMMARY <- tibble::tibble(
    species = 'Synthetic species',
    SRA = 6L,
    genomic = 3L,
    transcriptomic = 3L,
    epigenomic = 0L,
    chromatin = 0L,
    other = 0L,
    unknown = 0L
  )
  BIO_SUMMARY <- tibble::tibble(
    species = 'Synthetic species',
    BioSample = 6L,
    operable = 6L,
    aerial = 3L,
    ground = 3L,
    reproductive = 0L,
    whole = 0L,
    in_vitro = 0L,
    other = 0L,
    mixed = 0L,
    unknown = 0L
  )
  attr(SRA_SUMMARY, 'query_info') <- query_info
  attr(SRA_SUMMARY, 'gama_object') <- 'summarise_sra_availability'
  attr(SRA_SUMMARY, 'sra_profile') <- SRA_PROFILE
  attr(SRA_SUMMARY, 'sra_profile_info') <- list(
    cached_at_utc = '2026-05-22T11:23:11Z',
    profile_time_utc = '2026-05-22T11:23:11Z',
    id_col = 'entrez_uid',
    fields = names(SRA_PROFILE)
  )
  class(SRA_SUMMARY) <- unique(c('gdt_tbl', class(SRA_SUMMARY)))
  attr(BIO_SUMMARY, 'query_info') <- query_info
  attr(BIO_SUMMARY, 'gama_object') <- 'summarise_biosample_availability'
  attr(BIO_SUMMARY, 'biosample_anatomy_profile') <- BIO_ANATOMY
  attr(BIO_SUMMARY, 'biosample_anatomy_profile_info') <- list(
    cached_at_utc = '2026-05-22T11:23:11Z',
    profile_time_utc = '2026-05-22T11:23:11Z',
    id_col = 'biosample_id',
    fields = names(BIO_ANATOMY)
  )
  attr(BIO_SUMMARY, 'biosample_canonical_profile') <- BIO_CANONICAL
  attr(BIO_SUMMARY, 'biosample_canonical_profile_info') <- list(
    cached_at_utc = '2026-05-22T11:23:11Z',
    profile_time_utc = '2026-05-22T11:23:11Z',
    id_col = 'biosample_id',
    fields = names(BIO_CANONICAL)
  )
  class(BIO_SUMMARY) <- unique(c('gdt_tbl', class(BIO_SUMMARY)))
  list(SRA = SRA_SUMMARY, BIO = BIO_SUMMARY)
}

test_that('INTERACTION_CLASS fixture is a valid summarise_interaction object', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  expect_gdt_tbl(INTERACTION_CLASS)
  expect_identical(attr(INTERACTION_CLASS, 'gama_object', exact = TRUE), 'summarise_interaction')
})

test_that('INTERACTION_SUBCLASS fixture is a valid summarise_interaction object', {
  INTERACTION_SUBCLASS <- load_fixture('INTERACTION_SUBCLASS_Arabidopsis_thaliana')
  expect_gdt_tbl(INTERACTION_SUBCLASS)
  expect_identical(attr(INTERACTION_SUBCLASS, 'gama_object', exact = TRUE), 'summarise_interaction')
})

test_that('INTERACTION fixtures contain expected columns', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  INTERACTION_SUBCLASS <- load_fixture('INTERACTION_SUBCLASS_Arabidopsis_thaliana')
  expect_named(INTERACTION_CLASS, c('species', 'class', 'anatomy_class', 'BioSample', 'expected', 'residual'))
  expect_named(INTERACTION_SUBCLASS, c('species', 'class', 'anatomy_subclass', 'BioSample', 'expected', 'residual'))
  expect_identical(unique(as.character(INTERACTION_CLASS$species)), 'Arabidopsis thaliana')
  expect_identical(unique(as.character(INTERACTION_SUBCLASS$species)), 'Arabidopsis thaliana')
})

test_that('INTERACTION fixtures preserve query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  INTERACTION_SUBCLASS <- load_fixture('INTERACTION_SUBCLASS_Arabidopsis_thaliana')
  expect_identical(
    attr(INTERACTION_CLASS, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
  expect_identical(
    attr(INTERACTION_SUBCLASS, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
})

test_that('INTERACTION fixtures carry interaction metadata', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  INTERACTION_SUBCLASS <- load_fixture('INTERACTION_SUBCLASS_Arabidopsis_thaliana')
  fixtures <- list(
    anatomy_class = INTERACTION_CLASS,
    anatomy_subclass = INTERACTION_SUBCLASS
  )
  required_info <- c('level', 'term_col', 'term_levels', 'modality_levels', 'term_meta', 'match_report')
  for (level in names(fixtures)) {
    INFO <- attr(fixtures[[level]], 'interaction_info', exact = TRUE)
    expect_false(is.null(INFO))
    expect_true(all(required_info %in% names(INFO)))
    expect_identical(INFO$level, level)
    expect_identical(INFO$term_col, level)
    expect_true(length(INFO$term_levels) > 0L)
    expect_true(length(INFO$modality_levels) > 0L)
    expect_gama_columns(INFO$term_meta, c('term_order', 'x_term', 'anatomy_class'))
    expect_gama_columns(INFO$match_report, c('species', 'operable', 'matched_to_sra', 'unmatched_to_sra', 'matched_prop'))
  }
})

test_that('INTERACTION_CLASS fixture summarises modality-by-anatomy-class counts', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  INFO <- attr(INTERACTION_CLASS, 'interaction_info', exact = TRUE)
  expect_true(all(INTERACTION_CLASS$BioSample >= 0L))
  expect_true(all(INTERACTION_CLASS$class %in% INFO$modality_levels))
  expect_true(all(INTERACTION_CLASS$anatomy_class %in% INFO$term_levels))
  expect_equal(sum(INTERACTION_CLASS$BioSample), sum(INFO$match_report$matched_to_sra))
  expect_gt(sum(INTERACTION_CLASS$BioSample), 0L)
})

test_that('INTERACTION_SUBCLASS fixture summarises modality-by-anatomy-subclass counts', {
  INTERACTION_SUBCLASS <- load_fixture('INTERACTION_SUBCLASS_Arabidopsis_thaliana')
  INFO <- attr(INTERACTION_SUBCLASS, 'interaction_info', exact = TRUE)
  expect_true(all(INTERACTION_SUBCLASS$BioSample >= 0L))
  expect_true(all(INTERACTION_SUBCLASS$class %in% INFO$modality_levels))
  expect_true(all(INTERACTION_SUBCLASS$anatomy_subclass %in% INFO$term_levels))
  expect_equal(sum(INTERACTION_SUBCLASS$BioSample), sum(INFO$match_report$matched_to_sra))
  expect_gt(sum(INTERACTION_SUBCLASS$BioSample), 0L)
})

test_that('summarise_interaction applies correct expected-count formula', {
  SYN <- synthetic_interaction_inputs()
  INTERACTION <- suppressMessages(summarise_interaction(SYN$SRA, SYN$BIO, level = 'anatomy_class'))
  genomic_aerial <- INTERACTION[INTERACTION$class == 'genomic' & INTERACTION$anatomy_class == 'aerial', ]
  genomic_ground <- INTERACTION[INTERACTION$class == 'genomic' & INTERACTION$anatomy_class == 'ground', ]
  transcriptomic_aerial <- INTERACTION[INTERACTION$class == 'transcriptomic' & INTERACTION$anatomy_class == 'aerial', ]
  transcriptomic_ground <- INTERACTION[INTERACTION$class == 'transcriptomic' & INTERACTION$anatomy_class == 'ground', ]
  expect_equal(unname(genomic_aerial$BioSample), 2L)
  expect_equal(unname(genomic_ground$BioSample), 1L)
  expect_equal(unname(transcriptomic_aerial$BioSample), 1L)
  expect_equal(unname(transcriptomic_ground$BioSample), 2L)
  expect_equal(unname(genomic_aerial$expected), 1.5)
  expect_equal(unname(genomic_ground$expected), 1.5)
  expect_equal(unname(transcriptomic_aerial$expected), 1.5)
  expect_equal(unname(transcriptomic_ground$expected), 1.5)
})

test_that('summarise_interaction applies correct Pearson residual formula', {
  SYN <- synthetic_interaction_inputs()
  INTERACTION <- suppressMessages(summarise_interaction(SYN$SRA, SYN$BIO, level = 'anatomy_class'))
  genomic_aerial <- INTERACTION[INTERACTION$class == 'genomic' & INTERACTION$anatomy_class == 'aerial', ]
  genomic_ground <- INTERACTION[INTERACTION$class == 'genomic' & INTERACTION$anatomy_class == 'ground', ]
  expected_positive <- (2 - 1.5) / sqrt(1.5)
  expected_negative <- (1 - 1.5) / sqrt(1.5)
  expect_equal(unname(genomic_aerial$residual), expected_positive)
  expect_equal(unname(genomic_ground$residual), expected_negative)
})

test_that('INTERACTION fixtures follow correct expected-count formula', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  INTERACTION_SUBCLASS <- load_fixture('INTERACTION_SUBCLASS_Arabidopsis_thaliana')
  CLASS_EXPECTED <- expected_interaction_values(INTERACTION_CLASS)
  SUBCLASS_EXPECTED <- expected_interaction_values(INTERACTION_SUBCLASS)
  expect_equal(unname(INTERACTION_CLASS$expected), unname(CLASS_EXPECTED$expected_calc))
  expect_equal(unname(INTERACTION_SUBCLASS$expected), unname(SUBCLASS_EXPECTED$expected_calc))
})

test_that('INTERACTION fixtures follow correct Pearson residual formula', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  INTERACTION_SUBCLASS <- load_fixture('INTERACTION_SUBCLASS_Arabidopsis_thaliana')
  CLASS_EXPECTED <- expected_interaction_values(INTERACTION_CLASS)
  SUBCLASS_EXPECTED <- expected_interaction_values(INTERACTION_SUBCLASS)
  expect_equal(unname(INTERACTION_CLASS$residual), unname(CLASS_EXPECTED$residual_calc))
  expect_equal(unname(INTERACTION_SUBCLASS$residual), unname(SUBCLASS_EXPECTED$residual_calc))
})

test_that('plot_interaction returns ggplot objects', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  INTERACTION_SUBCLASS <- load_fixture('INTERACTION_SUBCLASS_Arabidopsis_thaliana')
  p1 <- plot_interaction(INTERACTION_CLASS, value = 'count')
  p2 <- plot_interaction(INTERACTION_CLASS, value = 'residual')
  p3 <- plot_interaction(INTERACTION_SUBCLASS, value = 'count')
  p4 <- plot_interaction(INTERACTION_SUBCLASS, value = 'residual')
  expect_s3_class(p1, 'ggplot')
  expect_s3_class(p2, 'ggplot')
  expect_s3_class(p3, 'ggplot')
  expect_s3_class(p4, 'ggplot')
})

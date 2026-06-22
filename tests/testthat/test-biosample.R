expected_biosample_skew_eff <- function(BIO_SKEW, PROFILE) {
  vapply(seq_len(nrow(BIO_SKEW)), function(i) {
    sp <- as.character(BIO_SKEW$species[[i]])
    anatomy_class_i <- as.character(BIO_SKEW$anatomy_class[[i]])
    prof <- PROFILE[PROFILE$species == sp, , drop = FALSE]
    if (!identical(anatomy_class_i, 'all')) prof <- prof[prof$anatomy_class == anatomy_class_i, , drop = FALSE]
    units <- trimws(as.character(prof$bioproject))
    units <- units[!is.na(units) & nzchar(units)]
    counts <- as.numeric(table(units))
    counts <- counts[counts > 0]
    if (!length(counts)) return(NA_real_)
    p <- counts / sum(counts)
    1 / sum(p^2)
  }, numeric(1))
}

test_that('BIO_SUMMARY fixture is a valid summarise_biosample_availability object', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_gdt_tbl(BIO_SUMMARY)
  expect_identical(attr(BIO_SUMMARY, 'gama_object', exact = TRUE), 'summarise_biosample_availability')
})

test_that('BIO_SUMMARY fixture contains expected columns', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  class_levels <- getFromNamespace('.biosample_anatomy_profile_levels', 'GAMA')()
  expected_cols <- c('species', 'BioSample', 'operable', class_levels)
  expect_gama_columns(BIO_SUMMARY, expected_cols)
  expect_identical(BIO_SUMMARY$species, 'Arabidopsis thaliana')
  expect_equal(nrow(BIO_SUMMARY), 1L)
  expect_true(length(setdiff(names(BIO_SUMMARY), expected_cols)) > 0L)
})

test_that('BIO_SUMMARY fixture preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_identical(
    attr(BIO_SUMMARY, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
})

test_that('BIO_SUMMARY fixture carries BioSample profile caches', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  ANATOMY <- attr(BIO_SUMMARY, 'biosample_anatomy_profile', exact = TRUE)
  CANONICAL <- attr(BIO_SUMMARY, 'biosample_canonical_profile', exact = TRUE)
  anatomy_cols <- c('species', 'biosample_id', 'bioproject', 'anatomy_class', 'anatomy_subclass')
  canonical_cols <- c(
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
  expect_false(is.null(ANATOMY))
  expect_false(is.null(CANONICAL))
  expect_gama_columns(ANATOMY, anatomy_cols)
  expect_gama_columns(CANONICAL, canonical_cols)
  expect_gt(nrow(ANATOMY), 0L)
  expect_gt(nrow(CANONICAL), 0L)
  expect_true(all(ANATOMY$species == 'Arabidopsis thaliana'))
  expect_true(all(CANONICAL$species == 'Arabidopsis thaliana'))
  expect_true(any(!is.na(ANATOMY$biosample_id) & nzchar(ANATOMY$biosample_id)))
  expect_true(any(!is.na(ANATOMY$bioproject) & nzchar(ANATOMY$bioproject)))
  expect_true(any(!is.na(CANONICAL$biosample_id) & nzchar(CANONICAL$biosample_id)))
})

test_that('BioSample anatomy profile collapsing follows expected mixed-state logic', {
  collapse_profile <- getFromNamespace('.biosample_collapse_anatomy_profile', 'GAMA')
  expect_identical(collapse_profile(c('aerial', 'aerial')), 'aerial')
  expect_identical(collapse_profile(c('aerial', 'ground')), 'mixed')
  expect_identical(collapse_profile(c(NA_character_, '')), 'unknown')
  expect_identical(collapse_profile(character()), 'unknown')
  expect_identical(collapse_profile(c('leaf', 'leaf'), level = 'anatomy_subclass'), 'leaf')
  expect_identical(collapse_profile(c('leaf', 'root'), level = 'anatomy_subclass'), 'mixed')
})

test_that('BIO_SUMMARY fixture summarises anatomy-class counts', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  class_levels <- getFromNamespace('.biosample_anatomy_profile_levels', 'GAMA')()
  counts <- unname(as.integer(unlist(BIO_SUMMARY[class_levels], use.names = FALSE)))
  expect_true(all(counts >= 0L))
  expect_equal(sum(counts), unname(BIO_SUMMARY$operable))
})

test_that('BIO_SUMMARY fixture follows BioSample count structure', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_equal(
    unname(BIO_SUMMARY$BioSample),
    as.integer(RESULTS[['Arabidopsis thaliana']]$biosample$count)
  )
  expect_gte(unname(BIO_SUMMARY$BioSample), unname(BIO_SUMMARY$operable))
  expect_gt(unname(BIO_SUMMARY$BioSample), 0L)
  expect_gt(unname(BIO_SUMMARY$operable), 0L)
})

test_that('BIO_SUMMARY fixture contains recognised anatomy profile classes', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  ANATOMY <- attr(BIO_SUMMARY, 'biosample_anatomy_profile', exact = TRUE)
  class_levels <- getFromNamespace('.biosample_anatomy_profile_levels', 'GAMA')()
  class_values <- ANATOMY$anatomy_class[!is.na(ANATOMY$anatomy_class) & nzchar(ANATOMY$anatomy_class)]
  expect_true(any(class_values %in% class_levels))
  expect_true(all(class_values %in% class_levels))
})

test_that('BIO_SUMMARY fixture contains recognised anatomy profile subclasses', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  ANATOMY <- attr(BIO_SUMMARY, 'biosample_anatomy_profile', exact = TRUE)
  subclass_levels <- unique(c(getFromNamespace('.biosample_anatomy_subclass_levels', 'GAMA')(), 'mixed', 'unknown'))
  subclass_values <- ANATOMY$anatomy_subclass[!is.na(ANATOMY$anatomy_subclass) & nzchar(ANATOMY$anatomy_subclass)]
  expect_true(any(subclass_values %in% subclass_levels))
  expect_true(all(subclass_values %in% subclass_levels))
})

test_that('BIO_SUMMARY fixture contains recognised canonical anatomy terms', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  CANONICAL <- attr(BIO_SUMMARY, 'biosample_canonical_profile', exact = TRUE)
  anatomy_ref <- getFromNamespace('.biosample_anatomy_ref', 'GAMA')()
  valid_terms <- unique(c(anatomy_ref$lex$anatomy_term, 'unknown'))
  term_values <- CANONICAL$anatomy_term[!is.na(CANONICAL$anatomy_term) & nzchar(CANONICAL$anatomy_term)]
  expect_true(any(term_values %in% valid_terms))
  expect_true(all(term_values %in% valid_terms))
})

test_that('plot_biosample_availability returns a ggplot object', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  p <- plot_biosample_availability(BIO_SUMMARY)
  expect_s3_class(p, 'ggplot')
})

test_that('BIO_SKEW fixture is a valid summarise_biosample_skew object', {
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  expect_gdt_tbl(BIO_SKEW)
  expect_identical(attr(BIO_SKEW, 'gama_object', exact = TRUE), 'summarise_biosample_skew')
})

test_that('BIO_SKEW fixture contains expected columns', {
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  expect_named(BIO_SKEW, c('species', 'BioProject', 'anatomy_class', 'min', 'q25', 'med', 'q75', 'max', 'eff'))
  expect_identical(BIO_SKEW$species, 'Arabidopsis thaliana')
})

test_that('BIO_SKEW fixture preserves query provenance', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  expect_identical(
    attr(BIO_SKEW, 'query_info', exact = TRUE),
    attr(BIO_SUMMARY, 'query_info', exact = TRUE)
  )
})

test_that('BIO_SKEW fixture carries ID recovery diagnostics', {
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  expect_skew_id_recovery(BIO_SKEW, unit = 'BioProject')
})

test_that('summarise_biosample_skew applies correct inverse Simpson index formula', {
  PROFILE <- tibble::tibble(
    species = rep('Synthetic species', 10),
    biosample_id = paste0('SAM', seq_len(10)),
    bioproject = rep(c('PRJ1', 'PRJ2', 'PRJ3'), c(5, 3, 2)),
    anatomy_class = rep('aerial', 10),
    anatomy_subclass = rep('leaf', 10)
  )
  BIO_SUMMARY <- tibble::tibble(
    species = 'Synthetic species',
    BioSample = 10L,
    operable = 10L,
    aerial = 10L,
    ground = 0L,
    reproductive = 0L,
    whole = 0L,
    `in vitro` = 0L,
    other = 0L,
    mixed = 0L,
    unknown = 0L
  )
  attr(BIO_SUMMARY, 'query_info') <- list(
    tool_version = test_gama_version(),
    query_time_utc = '2026-05-31T07:30:00Z',
    databases = c('assembly', 'sra', 'biosample'),
    terms = list('Synthetic species' = 'Synthetic species[Organism]'),
    synonyms = list('Synthetic species' = NULL)
  )
  attr(BIO_SUMMARY, 'gama_object') <- 'summarise_biosample_availability'
  attr(BIO_SUMMARY, 'biosample_anatomy_profile') <- PROFILE
  attr(BIO_SUMMARY, 'biosample_anatomy_profile_info') <- list(
    cached_at_utc = '2026-05-31T07:30:00Z',
    profile_time_utc = '2026-05-31T07:30:00Z',
    id_col = 'biosample_id',
    fields = names(PROFILE)
  )
  class(BIO_SUMMARY) <- unique(c('gdt_tbl', class(BIO_SUMMARY)))
  BIO_SKEW <- summarise_biosample_skew(BIO_SUMMARY, anatomy_class = 'aerial')
  p <- c(5, 3, 2) / 10
  expected_eff <- 1 / sum(p^2)
  expect_equal(unname(BIO_SKEW$BioProject), 3L)
  expect_equal(unname(BIO_SKEW$eff), expected_eff)
})

test_that('BIO_SKEW fixture follows correct inverse Simpson index formula', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  PROFILE <- attr(BIO_SUMMARY, 'biosample_anatomy_profile', exact = TRUE)
  expect_equal(
    unname(BIO_SKEW$eff),
    unname(expected_biosample_skew_eff(BIO_SKEW, PROFILE))
  )
})

test_that('plot_biosample_skew returns a ggplot object', {
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  p <- plot_biosample_skew(BIO_SKEW)
  expect_s3_class(p, 'ggplot')
})

test_that('BIO fixture is a valid extract_biosample_metadata object', {
  BIO <- load_fixture('BIO_Arabidopsis_thaliana')
  expect_gdt_tbl(BIO)
  expect_identical(attr(BIO, 'gama_object', exact = TRUE), 'extract_biosample_metadata')
  expect_gt(nrow(BIO), 0L)
})

test_that('BIO fixture contains expected columns', {
  BIO <- load_fixture('BIO_Arabidopsis_thaliana')
  expected_cols <- c(
    'species',
    'entrez_uid',
    'biosample',
    'bioproject',
    'value_raw',
    'value_norm',
    'anatomy_class',
    'anatomy_subclass',
    'anatomy_term',
    'anatomy_class_profile',
    'anatomy_subclass_profile'
  )
  expect_named(BIO, expected_cols)
  expect_true(all(BIO$species == 'Arabidopsis thaliana'))
  expect_true(any(!is.na(BIO$entrez_uid) & nzchar(BIO$entrez_uid)))
  expect_true(any(!is.na(BIO$biosample) & nzchar(BIO$biosample)))
})

test_that('BIO fixture preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  BIO <- load_fixture('BIO_Arabidopsis_thaliana')
  expect_identical(
    attr(BIO, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
})

test_that('BIO fixture contains recognised anatomy classes', {
  BIO <- load_fixture('BIO_Arabidopsis_thaliana')
  valid_classes <- unique(c(getFromNamespace('.biosample_anatomy_class_levels', 'GAMA')(), 'unknown'))
  class_values <- BIO$anatomy_class[!is.na(BIO$anatomy_class) & nzchar(BIO$anatomy_class)]
  expect_true(any(class_values %in% valid_classes))
  expect_true(all(class_values %in% valid_classes))
  expect_false(any(class_values == 'mixed'))
})

test_that('BIO fixture contains recognised anatomy subclasses', {
  BIO <- load_fixture('BIO_Arabidopsis_thaliana')
  valid_subclasses <- unique(c(getFromNamespace('.biosample_anatomy_subclass_levels', 'GAMA')(), 'unknown'))
  subclass_values <- BIO$anatomy_subclass[!is.na(BIO$anatomy_subclass) & nzchar(BIO$anatomy_subclass)]
  expect_true(any(subclass_values %in% valid_subclasses))
  expect_true(all(subclass_values %in% valid_subclasses))
  expect_false(any(subclass_values == 'mixed'))
})

test_that('BIO fixture contains recognised anatomy terms', {
  BIO <- load_fixture('BIO_Arabidopsis_thaliana')
  anatomy_ref <- getFromNamespace('.biosample_anatomy_ref', 'GAMA')()
  valid_terms <- unique(c(anatomy_ref$lex$anatomy_term, 'unknown'))
  term_values <- BIO$anatomy_term[!is.na(BIO$anatomy_term) & nzchar(BIO$anatomy_term)]
  expect_true(any(term_values %in% valid_terms))
  expect_true(all(term_values %in% valid_terms))
})

test_that('BIO fixture preserves collapsed anatomy profile structure', {
  BIO <- load_fixture('BIO_Arabidopsis_thaliana')
  valid_class_profiles <- getFromNamespace('.biosample_anatomy_profile_levels', 'GAMA')()
  valid_subclass_profiles <- unique(c(getFromNamespace('.biosample_anatomy_subclass_levels', 'GAMA')(), 'mixed', 'unknown'))
  class_profiles <- BIO$anatomy_class_profile[!is.na(BIO$anatomy_class_profile) & nzchar(BIO$anatomy_class_profile)]
  subclass_profiles <- BIO$anatomy_subclass_profile[!is.na(BIO$anatomy_subclass_profile) & nzchar(BIO$anatomy_subclass_profile)]
  expect_true(any(class_profiles %in% valid_class_profiles))
  expect_true(all(class_profiles %in% valid_class_profiles))
  expect_true(any(subclass_profiles %in% valid_subclass_profiles))
  expect_true(all(subclass_profiles %in% valid_subclass_profiles))
})

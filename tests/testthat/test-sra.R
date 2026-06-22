expected_sra_skew_eff <- function(SRA_SKEW, PROFILE, unit_col) {
  vapply(seq_len(nrow(SRA_SKEW)), function(i) {
    sp <- as.character(SRA_SKEW$species[[i]])
    class_i <- as.character(SRA_SKEW$class[[i]])
    prof <- PROFILE[PROFILE$species == sp, , drop = FALSE]
    if (!identical(class_i, 'all')) prof <- prof[prof$class == class_i, , drop = FALSE]
    units <- prof[[unit_col]]
    units <- units[!is.na(units)]
    counts <- as.numeric(table(units))
    counts <- counts[counts > 0]
    if (!length(counts)) return(NA_real_)
    p <- counts / sum(counts)
    1 / sum(p^2)
  }, numeric(1))
}

test_that('SRA_SUMMARY fixture is a valid summarise_sra_availability object', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gdt_tbl(SRA_SUMMARY)
  expect_identical(attr(SRA_SUMMARY, 'gama_object', exact = TRUE), 'summarise_sra_availability')
})

test_that('SRA_SUMMARY fixture contains expected columns', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  classes <- c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
  geo_cols <- c('denom_total', 'geo_linked_denom', 'geo_prop', paste0(classes, '_geo'))
  expected_cols <- c('species', 'SRA', classes, geo_cols)
  expect_gama_columns(SRA_SUMMARY, expected_cols)
  expect_identical(SRA_SUMMARY$species, 'Arabidopsis thaliana')
  expect_equal(nrow(SRA_SUMMARY), 1L)
  expect_true(length(setdiff(names(SRA_SUMMARY), expected_cols)) > 0L)
})

test_that('SRA_SUMMARY fixture preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_identical(
    attr(SRA_SUMMARY, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
})

test_that('SRA_SUMMARY fixture carries SRA profile cache', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  PROFILE <- attr(SRA_SUMMARY, 'sra_profile', exact = TRUE)
  PROFILE_INFO <- attr(SRA_SUMMARY, 'sra_profile_info', exact = TRUE)
  expected_profile_cols <- c('species', 'entrez_uid', 'biosample', 'bioproject', 'class', 'subclass', 'geo_linked', 'gse_ids', 'gsm_ids')
  expected_info_fields <- c('cached_at_utc', 'profile_time_utc', 'id_col', 'fields')
  expect_false(is.null(PROFILE))
  expect_false(is.null(PROFILE_INFO))
  expect_gama_columns(PROFILE, expected_profile_cols)
  expect_true(all(expected_info_fields %in% names(PROFILE_INFO)))
  expect_gt(nrow(PROFILE), 0L)
  expect_identical(PROFILE_INFO$id_col, 'entrez_uid')
  expect_identical(PROFILE_INFO$fields, expected_profile_cols)
  expect_true(any(!is.na(PROFILE$species) & nzchar(PROFILE$species)))
  expect_true(any(!is.na(PROFILE$biosample) & nzchar(PROFILE$biosample)))
  expect_true(any(!is.na(PROFILE$class) & nzchar(PROFILE$class)))
  expect_true(any(!is.na(PROFILE$subclass) & nzchar(PROFILE$subclass)))
  expect_true(any(PROFILE$species == 'Arabidopsis thaliana' & !is.na(PROFILE$biosample) & nzchar(PROFILE$biosample)))
})

test_that('SRA_SUMMARY fixture summarises modality counts', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  PROFILE <- attr(SRA_SUMMARY, 'sra_profile', exact = TRUE)
  classes <- c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
  expected_counts <- vapply(classes, function(class) {
    as.integer(sum(PROFILE$class == class, na.rm = TRUE))
  }, integer(1))
  observed_counts <- unname(as.integer(unlist(SRA_SUMMARY[classes], use.names = FALSE)))
  expect_equal(unname(SRA_SUMMARY$SRA), nrow(PROFILE))
  expect_equal(observed_counts, unname(expected_counts))
  expect_equal(sum(observed_counts), nrow(PROFILE))
})

test_that('SRA_SUMMARY fixture follows GEO denominator logic', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  PROFILE <- attr(SRA_SUMMARY, 'sra_profile', exact = TRUE)
  classes <- c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
  denom_classes <- c('transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
  expected_denom <- as.integer(sum(PROFILE$class %in% denom_classes))
  expected_geo_denom <- as.integer(sum(PROFILE$class %in% denom_classes & PROFILE$geo_linked))
  expected_geo_prop <- if (expected_denom > 0L) expected_geo_denom / expected_denom else NA_real_
  expect_equal(unname(SRA_SUMMARY$denom_total), expected_denom)
  expect_equal(unname(SRA_SUMMARY$geo_linked_denom), expected_geo_denom)
  expect_equal(unname(SRA_SUMMARY$geo_prop), expected_geo_prop)
  for (class in classes) {
    geo_col <- paste0(class, '_geo')
    expect_equal(unname(SRA_SUMMARY[[geo_col]]), as.integer(sum(PROFILE$class == class & PROFILE$geo_linked)))
  }
})

test_that('plot_sra_availability returns a ggplot object', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  p <- plot_sra_availability(SRA_SUMMARY)
  expect_s3_class(p, 'ggplot')
})

test_that('plot_sra_geo returns a ggplot object', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  p <- plot_sra_geo(SRA_SUMMARY)
  expect_s3_class(p, 'ggplot')
})

test_that('SRA_SKEW_BIOPROJECT fixture is a valid summarise_sra_skew object', {
  SRA_SKEW_BIOPROJECT <- load_fixture('SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana')
  expect_gdt_tbl(SRA_SKEW_BIOPROJECT)
  expect_identical(attr(SRA_SKEW_BIOPROJECT, 'gama_object', exact = TRUE), 'summarise_sra_skew')
})

test_that('SRA_SKEW_BIOSAMPLE fixture is a valid summarise_sra_skew object', {
  SRA_SKEW_BIOSAMPLE <- load_fixture('SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana')
  expect_gdt_tbl(SRA_SKEW_BIOSAMPLE)
  expect_identical(attr(SRA_SKEW_BIOSAMPLE, 'gama_object', exact = TRUE), 'summarise_sra_skew')
})

test_that('SRA skew fixtures contain expected columns', {
  SRA_SKEW_BIOPROJECT <- load_fixture('SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana')
  SRA_SKEW_BIOSAMPLE <- load_fixture('SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana')
  expect_named(SRA_SKEW_BIOPROJECT, c('species', 'BioProject', 'class', 'min', 'q25', 'med', 'q75', 'max', 'eff'))
  expect_named(SRA_SKEW_BIOSAMPLE, c('species', 'BioSample', 'class', 'min', 'q25', 'med', 'q75', 'max', 'eff'))
  expect_identical(SRA_SKEW_BIOPROJECT$species, 'Arabidopsis thaliana')
  expect_identical(SRA_SKEW_BIOSAMPLE$species, 'Arabidopsis thaliana')
})

test_that('SRA skew fixtures preserve query provenance', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  SRA_SKEW_BIOPROJECT <- load_fixture('SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana')
  SRA_SKEW_BIOSAMPLE <- load_fixture('SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana')
  expect_identical(
    attr(SRA_SKEW_BIOPROJECT, 'query_info', exact = TRUE),
    attr(SRA_SUMMARY, 'query_info', exact = TRUE)
  )
  expect_identical(
    attr(SRA_SKEW_BIOSAMPLE, 'query_info', exact = TRUE),
    attr(SRA_SUMMARY, 'query_info', exact = TRUE)
  )
})

test_that('SRA skew fixtures carry ID recovery diagnostics', {
  SRA_SKEW_BIOPROJECT <- load_fixture('SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana')
  SRA_SKEW_BIOSAMPLE <- load_fixture('SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana')
  expect_skew_id_recovery(SRA_SKEW_BIOPROJECT, unit = 'BioProject')
  expect_skew_id_recovery(SRA_SKEW_BIOSAMPLE, unit = 'BioSample')
})

test_that('summarise_sra_skew applies correct inverse Simpson index formula', {
  PROFILE <- tibble::tibble(
    species = rep('Synthetic species', 10),
    entrez_uid = paste0('SRX', seq_len(10)),
    biosample = paste0('SAM', seq_len(10)),
    bioproject = rep(c('PRJ1', 'PRJ2', 'PRJ3'), c(5, 3, 2)),
    class = rep('genomic', 10),
    subclass = rep('WGS', 10),
    geo_linked = rep(FALSE, 10),
    gse_ids = rep(NA_character_, 10),
    gsm_ids = rep(NA_character_, 10)
  )
  SRA_SUMMARY <- tibble::tibble(
    species = 'Synthetic species',
    SRA = 10L,
    genomic = 10L,
    transcriptomic = 0L,
    epigenomic = 0L,
    chromatin = 0L,
    other = 0L,
    unknown = 0L
  )
  attr(SRA_SUMMARY, 'query_info') <- list(
    tool_version = test_gama_version(),
    query_time_utc = '2026-05-22T11:23:11Z',
    databases = c('assembly', 'sra', 'biosample'),
    terms = list('Synthetic species' = 'Synthetic species[Organism]'),
    synonyms = list('Synthetic species' = NULL)
  )
  attr(SRA_SUMMARY, 'gama_object') <- 'summarise_sra_availability'
  attr(SRA_SUMMARY, 'sra_profile') <- PROFILE
  attr(SRA_SUMMARY, 'sra_profile_info') <- list(
    cached_at_utc = '2026-05-22T11:23:11Z',
    profile_time_utc = '2026-05-22T11:23:11Z',
    id_col = 'entrez_uid',
    fields = names(PROFILE)
  )
  class(SRA_SUMMARY) <- unique(c('gdt_tbl', class(SRA_SUMMARY)))
  SRA_SKEW <- summarise_sra_skew(SRA_SUMMARY, unit = 'bioproject')
  p <- c(5, 3, 2) / 10
  expected_eff <- 1 / sum(p^2)
  expect_equal(unname(SRA_SKEW$BioProject), 3L)
  expect_equal(unname(SRA_SKEW$eff), expected_eff)
})

test_that('SRA skew fixtures follow correct inverse Simpson index formula', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  SRA_SKEW_BIOPROJECT <- load_fixture('SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana')
  SRA_SKEW_BIOSAMPLE <- load_fixture('SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana')
  PROFILE <- attr(SRA_SUMMARY, 'sra_profile', exact = TRUE)
  expect_equal(
    unname(SRA_SKEW_BIOPROJECT$eff),
    unname(expected_sra_skew_eff(SRA_SKEW_BIOPROJECT, PROFILE, 'bioproject'))
  )
  expect_equal(
    unname(SRA_SKEW_BIOSAMPLE$eff),
    unname(expected_sra_skew_eff(SRA_SKEW_BIOSAMPLE, PROFILE, 'biosample'))
  )
})

test_that('plot_sra_skew returns ggplot objects', {
  SRA_SKEW_BIOPROJECT <- load_fixture('SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana')
  SRA_SKEW_BIOSAMPLE <- load_fixture('SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana')
  p1 <- plot_sra_skew(SRA_SKEW_BIOPROJECT)
  p2 <- plot_sra_skew(SRA_SKEW_BIOSAMPLE)
  expect_s3_class(p1, 'ggplot')
  expect_s3_class(p2, 'ggplot')
})

test_that('SRA fixture is a valid extract_sra_metadata object', {
  SRA <- load_fixture('SRA_Arabidopsis_thaliana')
  expect_gdt_tbl(SRA)
  expect_identical(attr(SRA, 'gama_object', exact = TRUE), 'extract_sra_metadata')
  expect_gt(nrow(SRA), 0L)
})

test_that('SRA fixture contains expected columns', {
  SRA <- load_fixture('SRA_Arabidopsis_thaliana')
  expected_cols <- c('species', 'entrez_uid', 'biosample', 'bioproject', 'strategy_raw', 'strategy_norm', 'class', 'subclass', 'geo_linked', 'gse_ids', 'gsm_ids')
  expect_named(SRA, expected_cols)
  expect_true(all(SRA$species == 'Arabidopsis thaliana'))
  expect_true(any(!is.na(SRA$entrez_uid) & nzchar(SRA$entrez_uid)))
})

test_that('SRA fixture preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  SRA <- load_fixture('SRA_Arabidopsis_thaliana')
  expect_identical(
    attr(SRA, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
})

test_that('SRA fixture contains recognised modality classes', {
  SRA <- load_fixture('SRA_Arabidopsis_thaliana')
  valid_classes <- unique(c(names(getFromNamespace('.ONTOLOGY', 'GAMA')), 'other', 'unknown'))
  class_values <- SRA$class[!is.na(SRA$class) & nzchar(SRA$class)]
  expect_true(any(class_values %in% valid_classes))
  expect_true(all(class_values %in% valid_classes))
})

test_that('SRA fixture contains recognised modality subclasses', {
  SRA <- load_fixture('SRA_Arabidopsis_thaliana')
  ontology <- getFromNamespace('.ONTOLOGY', 'GAMA')
  valid_subclasses <- unique(c(unlist(lapply(ontology, names), use.names = FALSE), 'other', 'unknown'))
  subclass_values <- SRA$subclass[!is.na(SRA$subclass) & nzchar(SRA$subclass)]
  expect_true(any(subclass_values %in% valid_subclasses))
  expect_true(all(subclass_values %in% valid_subclasses))
})

test_that('SRA fixture preserves GEO linkage structure', {
  SRA <- load_fixture('SRA_Arabidopsis_thaliana')
  linked <- SRA$geo_linked
  has_gse <- !is.na(SRA$gse_ids) & nzchar(SRA$gse_ids)
  has_gsm <- !is.na(SRA$gsm_ids) & nzchar(SRA$gsm_ids)
  expect_type(linked, 'logical')
  expect_true(all(!linked | has_gse | has_gsm))
})

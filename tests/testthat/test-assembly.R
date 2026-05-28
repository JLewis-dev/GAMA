test_that('ASM_SUMMARY fixture is a valid summarise_assembly_availability object', {
  ASM_SUMMARY <- load_fixture('ASM_SUMMARY_Arabidopsis_thaliana')
  expect_gdt_tbl(ASM_SUMMARY)
  expect_identical(attr(ASM_SUMMARY, 'gama_object', exact = TRUE), 'summarise_assembly_availability')
})

test_that('ASM_SUMMARY fixture contains expected columns', {
  ASM_SUMMARY <- load_fixture('ASM_SUMMARY_Arabidopsis_thaliana')
  expect_named(ASM_SUMMARY, c('species', 'Assembly', 'complete', 'chromosome', 'scaffold', 'contig', 'best_n50'))
  expect_identical(ASM_SUMMARY$species, 'Arabidopsis thaliana')
  expect_equal(nrow(ASM_SUMMARY), 1L)
})

test_that('ASM_SUMMARY fixture preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  ASM_SUMMARY <- load_fixture('ASM_SUMMARY_Arabidopsis_thaliana')
  expect_identical(
    attr(ASM_SUMMARY, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
})

test_that('ASM_SUMMARY fixture summarises assembly level counts', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  ASM_SUMMARY <- load_fixture('ASM_SUMMARY_Arabidopsis_thaliana')
  ASM <- load_fixture('ASM_Arabidopsis_thaliana')
  assembly_level_class <- getFromNamespace('.assembly_level_class', 'GAMA')
  LEVELS <- c('complete', 'chromosome', 'scaffold', 'contig')
  level_class <- assembly_level_class(ASM$level)
  expected_counts <- vapply(LEVELS, function(level) {
    as.integer(sum(level_class == level, na.rm = TRUE))
  }, integer(1))
  observed_counts <- unname(as.integer(unlist(ASM_SUMMARY[LEVELS], use.names = FALSE)))
  expect_equal(unname(ASM_SUMMARY$Assembly), as.integer(RESULTS[['Arabidopsis thaliana']]$assembly$count))
  expect_equal(observed_counts, unname(expected_counts))
  expect_equal(sum(observed_counts), sum(!is.na(level_class)))
})

test_that('ASM_SUMMARY fixture follows best_n50 selection logic', {
  ASM_SUMMARY <- load_fixture('ASM_SUMMARY_Arabidopsis_thaliana')
  ASM <- load_fixture('ASM_Arabidopsis_thaliana')
  best_assembly_n50 <- getFromNamespace('.best_assembly_n50', 'GAMA')
  expected_best_n50 <- best_assembly_n50(ASM$level, ASM$n50)
  expect_equal(unname(ASM_SUMMARY$best_n50), unname(expected_best_n50))
})

test_that('plot_assembly_availability returns a ggplot object', {
  ASM_SUMMARY <- load_fixture('ASM_SUMMARY_Arabidopsis_thaliana')
  p <- plot_assembly_availability(ASM_SUMMARY)
  expect_s3_class(p, 'ggplot')
})

test_that('ASM fixture is a valid extract_assembly_metadata object', {
  ASM <- load_fixture('ASM_Arabidopsis_thaliana')
  expect_gdt_tbl(ASM)
  expect_identical(attr(ASM, 'gama_object', exact = TRUE), 'extract_assembly_metadata')
  expect_gt(nrow(ASM), 0L)
})

test_that('ASM fixture contains expected columns', {
  ASM <- load_fixture('ASM_Arabidopsis_thaliana')
  expect_named(ASM, c('species', 'entrez_uid', 'level', 'n50', 'coverage', 'biosample', 'bioproject', 'submitter', 'release_date', 'ftp_path'))
  expect_true(all(ASM$species == 'Arabidopsis thaliana'))
  expect_true(any(!is.na(ASM$entrez_uid) & nzchar(ASM$entrez_uid)))
})

test_that('ASM fixture preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  ASM <- load_fixture('ASM_Arabidopsis_thaliana')
  expect_identical(
    attr(ASM, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
})

test_that('ASM fixture contains no unrecognised assembly levels', {
  ASM <- load_fixture('ASM_Arabidopsis_thaliana')
  assembly_level_class <- getFromNamespace('.assembly_level_class', 'GAMA')
  non_missing <- !is.na(ASM$level) & nzchar(ASM$level)
  level_class <- assembly_level_class(ASM$level[non_missing])
  expect_true(any(non_missing))
  expect_false(any(is.na(level_class)))
})

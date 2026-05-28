test_that('SUMMARY fixture is a valid summarise_availability object', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gdt_tbl(SUMMARY)
  expect_identical(attr(SUMMARY, 'gama_object', exact = TRUE), 'summarise_availability')
})

test_that('SUMMARY fixture contains expected columns', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_named(SUMMARY, c('species', 'Assembly', 'SRA', 'BioSample', 'A', 'S', 'B', 'score'))
  expect_identical(SUMMARY$species, 'Arabidopsis thaliana')
  expect_equal(nrow(SUMMARY), 1L)
})

test_that('SUMMARY fixture preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_identical(
    attr(SUMMARY, 'query_info', exact = TRUE),
    attr(RESULTS, 'query_info', exact = TRUE)
  )
})

test_that('scoring uses expected assembly level weights', {
  assembly_level_weight <- getFromNamespace('.assembly_level_weight', 'GAMA')
  levels <- c(
    'Complete Genome',
    'Chromosome',
    'Scaffold',
    'Contig',
    'complete',
    'chromosome',
    'scaffold',
    'contig',
    'unrecognised',
    NA_character_
  )
  expect_equal(
    assembly_level_weight(levels),
    c(10, 8, 5, 2, 10, 8, 5, 2, 0, 0)
  )
})

test_that('Assembly fixture contains no unrecognised assembly levels', {
  ASM <- load_fixture('ASM_Arabidopsis_thaliana')
  assembly_level_weight <- getFromNamespace('.assembly_level_weight', 'GAMA')
  non_missing <- !is.na(ASM$level) & nzchar(ASM$level)
  expect_true(any(non_missing))
  expect_true(all(assembly_level_weight(ASM$level[non_missing]) > 0))
})

test_that('summarise_availability applies correct data richness formula', {
  fake_summaries <- list(
    a = list(assembly_level = 'Complete Genome', contign50 = 200),
    b = list(assembly_level = 'Scaffold', contign50 = 800),
    c = list(assembly_level = 'Contig', contign50 = 300)
  )
  RESULTS <- list(
    'Synthetic species' = list(
      assembly = list(count = 3L, ids = c('a', 'b', 'c')),
      sra = list(count = 9L),
      biosample = list(count = 4L)
    )
  )
  attr(RESULTS, 'query_info') <- list(
    tool_version = '0.3.0',
    query_time_utc = '2026-05-22T11:23:11Z',
    databases = c('assembly', 'sra', 'biosample'),
    terms = list('Synthetic species' = 'Synthetic species[Organism]'),
    synonyms = list('Synthetic species' = NULL)
  )
  attr(RESULTS, 'gama_object') <- 'query_species'
  testthat::local_mocked_bindings(
    .fetch_search_summaries = function(db, search, batch_size = 100) fake_summaries,
    .package = 'GAMA'
  )
  SUMMARY <- summarise_availability(RESULTS)
  expected_A <- 10 + log1p((10 + 5 + 2) - 10)
  expected_S <- 2 * log1p(9)
  expected_B <- log1p(4)
  expect_equal(unname(SUMMARY$A), expected_A)
  expect_equal(unname(SUMMARY$S), expected_S)
  expect_equal(unname(SUMMARY$B), expected_B)
  expect_equal(unname(SUMMARY$score), expected_A + expected_S + expected_B)
})

test_that('SUMMARY fixture follows correct data richness formula', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  ASM <- load_fixture('ASM_Arabidopsis_thaliana')
  assembly_level_weight <- getFromNamespace('.assembly_level_weight', 'GAMA')
  weights <- assembly_level_weight(ASM$level)
  expected_A <- if (length(weights) && any(weights > 0)) {
    best <- max(weights)
    best + log1p(sum(weights) - best)
  } else {
    0
  }
  counts <- unlist(SUMMARY[c('Assembly', 'SRA', 'BioSample')], use.names = FALSE)
  scores <- unlist(SUMMARY[c('A', 'S', 'B', 'score')], use.names = FALSE)
  expect_true(all(counts >= 0))
  expect_true(all(is.finite(scores)))
  expect_true(all(scores >= 0))
  expect_equal(unname(SUMMARY$A), expected_A)
  expect_equal(unname(SUMMARY$S), unname(2 * log1p(SUMMARY$SRA)))
  expect_equal(unname(SUMMARY$B), unname(log1p(SUMMARY$BioSample)))
  expect_equal(unname(SUMMARY$score), unname(SUMMARY$A + SUMMARY$S + SUMMARY$B))
})

test_that('plot_availability returns a ggplot object', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  p <- plot_availability(SUMMARY)
  expect_s3_class(p, 'ggplot')
})

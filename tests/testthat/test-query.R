test_that('RESULTS fixture is a valid query_species object', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_type(RESULTS, 'list')
  expect_named(RESULTS, 'Arabidopsis thaliana')
  expect_identical(attr(RESULTS, 'gama_object', exact = TRUE), 'query_species')
  expect_false(is.null(attr(RESULTS, 'query_info', exact = TRUE)))
})

test_that('RESULTS fixture stores one species with all queried databases', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  x <- RESULTS[['Arabidopsis thaliana']]
  expect_type(x, 'list')
  expect_named(x, c('assembly', 'sra', 'biosample'))
  for (db in c('assembly', 'sra', 'biosample')) {
    expect_type(x[[db]], 'list')
    expect_true('count' %in% names(x[[db]]))
    expect_true(is.numeric(x[[db]]$count))
    expect_length(x[[db]]$count, 1L)
    expect_gte(x[[db]]$count, 0)
  }
})

test_that('RESULTS fixture preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  query_info <- attr(RESULTS, 'query_info', exact = TRUE)
  expect_named(query_info, c('tool_version', 'query_time_utc', 'databases', 'terms', 'synonyms'))
  expect_match(query_info$query_time_utc, '^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$')
  expect_identical(query_info$databases, c('assembly', 'sra', 'biosample'))
  expect_named(query_info$terms, 'Arabidopsis thaliana')
  expect_identical(query_info$terms[['Arabidopsis thaliana']], 'Arabidopsis thaliana[Organism]')
  expect_named(query_info$synonyms, 'Arabidopsis thaliana')
  expect_null(query_info$synonyms[['Arabidopsis thaliana']])
})

test_that('fixture manifest preserves query provenance', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  MANIFEST <- load_fixture('Arabidopsis_thaliana_fixture_manifest')
  query_info <- attr(RESULTS, 'query_info', exact = TRUE)
  expect_true('query_time_utc' %in% names(MANIFEST))
  expect_true('tool_version' %in% names(MANIFEST))
  expect_true(all(MANIFEST$query_time_utc == query_info$query_time_utc))
  expect_true(all(MANIFEST$tool_version == query_info$tool_version))
})

test_that('fixture manifest records the complete fixture set', {
  MANIFEST <- load_fixture('Arabidopsis_thaliana_fixture_manifest')
  expected_objects <- c(
    'RESULTS_Arabidopsis_thaliana',
    'SUMMARY_Arabidopsis_thaliana',
    'ASM_SUMMARY_Arabidopsis_thaliana',
    'ASM_Arabidopsis_thaliana',
    'SRA_SUMMARY_Arabidopsis_thaliana',
    'SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana',
    'SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana',
    'SRA_Arabidopsis_thaliana',
    'BIO_SUMMARY_Arabidopsis_thaliana',
    'BIO_SKEW_Arabidopsis_thaliana',
    'BIO_Arabidopsis_thaliana',
    'INTERACTION_CLASS_Arabidopsis_thaliana',
    'INTERACTION_SUBCLASS_Arabidopsis_thaliana'
  )
  manifest_paths <- vapply(MANIFEST$file, fixture_path, character(1))
  expect_named(MANIFEST, c('object', 'file', 'size_bytes', 'query_time_utc', 'tool_version'))
  expect_identical(MANIFEST$object, expected_objects)
  expect_identical(MANIFEST$file, paste0(expected_objects, '.rds'))
  expect_true(all(file.exists(manifest_paths)))
  expect_equal(MANIFEST$size_bytes, unname(file.info(manifest_paths)$size))
  expect_true(all(MANIFEST$size_bytes > 0))
})

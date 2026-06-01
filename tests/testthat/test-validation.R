expect_gama_error <- function(expr, regexp = NULL) {
  expr <- substitute(expr)
  err <- tryCatch(eval(expr, parent.frame()), error = identity)
  testthat::expect_s3_class(err, 'error')
  if (!inherits(err, 'error')) return(invisible(err))
  msg <- conditionMessage(err)
  version <- getFromNamespace('.GAMA_VERSION', 'GAMA')
  testthat::expect_true(grepl(paste0('GAMA ', version, ' | '), msg, fixed = TRUE))
  if (!is.null(regexp)) testthat::expect_match(msg, regexp)
  invisible(err)
}

without_attr <- function(x, attr_name) {
  attr(x, attr_name) <- NULL
  x
}

test_that('GAMA validation errors use versioned messages', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_availability(SUMMARY, rank = 'not-a-rank'),
    'Invalid `rank` parameter'
  )
})

test_that('GAMA validation distinguishes GAMA and non-GAMA input objects', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_availability(tibble::tibble(x = 1)),
    'detected incompatible object'
  )
  expect_gama_error(
    summarise_availability(SRA_SUMMARY),
    'detected summarise_sra_availability object'
  )
})

test_that('query_species validation rejects invalid species input', {
  expect_gama_error(
    query_species(character()),
    '`species` must contain at least one valid species name'
  )
})

test_that('query_species validation rejects invalid synonym input', {
  expect_gama_error(
    query_species('Arabidopsis thaliana', synonyms = list('Arabidopsis lyrata')),
    '`synonyms` must be a named list or named character vector'
  )
})

test_that('summarise_availability validation rejects incompatible input objects', {
  expect_gama_error(
    summarise_availability(tibble::tibble(x = 1)),
    'output of query_species\\(\\).*detected incompatible object'
  )
})

test_that('summarise_availability validation reports detected GAMA object types', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_availability(SRA_SUMMARY),
    'output of query_species\\(\\).*detected summarise_sra_availability object'
  )
})

test_that('plot_availability validation rejects incompatible input objects', {
  expect_gama_error(
    plot_availability(tibble::tibble(x = 1)),
    'output of summarise_availability\\(\\).*detected incompatible object'
  )
})

test_that('plot_availability validation rejects invalid rank parameters', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_availability(SUMMARY, rank = 'higest'),
    'Invalid `rank` parameter.*Did you mean.*highest'
  )
})

test_that('plot_availability validation rejects invalid logical parameters', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_availability(SUMMARY, abbreviate = 'yes'),
    '`abbreviate` must be TRUE or FALSE'
  )
})

test_that('summarise_assembly_availability validation rejects incompatible input objects', {
  expect_gama_error(
    summarise_assembly_availability(tibble::tibble(x = 1)),
    'output of query_species\\(\\).*detected incompatible object'
  )
})

test_that('summarise_assembly_availability validation reports detected GAMA object types', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_assembly_availability(SUMMARY),
    'output of query_species\\(\\).*detected summarise_availability object'
  )
})

test_that('plot_assembly_availability validation rejects incompatible input objects', {
  expect_gama_error(
    plot_assembly_availability(tibble::tibble(x = 1)),
    'output of summarise_assembly_availability\\(\\).*detected incompatible object'
  )
})

test_that('plot_assembly_availability validation rejects invalid rank parameters', {
  ASM_SUMMARY <- load_fixture('ASM_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_assembly_availability(ASM_SUMMARY, rank = 'lowst'),
    'Invalid `rank` parameter.*Did you mean.*lowest'
  )
})

test_that('plot_assembly_availability validation rejects invalid logical parameters', {
  ASM_SUMMARY <- load_fixture('ASM_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_assembly_availability(ASM_SUMMARY, abbreviate = NA),
    '`abbreviate` must be TRUE or FALSE'
  )
})

test_that('extract_assembly_metadata validation rejects incompatible input objects', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    extract_assembly_metadata(SUMMARY),
    'output of query_species\\(\\).*detected summarise_availability object'
  )
})

test_that('extract_assembly_metadata validation rejects invalid logical parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_assembly_metadata(RESULTS, best = 'yes'),
    '`best` must be TRUE or FALSE'
  )
})

test_that('summarise_sra_availability validation rejects incompatible input objects', {
  expect_gama_error(
    summarise_sra_availability(tibble::tibble(x = 1)),
    'output of query_species\\(\\).*detected incompatible object'
  )
})

test_that('summarise_sra_availability validation reports detected GAMA object types', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_sra_availability(SUMMARY),
    'output of query_species\\(\\).*detected summarise_availability object'
  )
})

test_that('summarise_sra_availability validation rejects invalid logical parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_sra_availability(RESULTS, all = 'yes'),
    '`all` must be TRUE or FALSE'
  )
  expect_gama_error(
    summarise_sra_availability(RESULTS, include_geo = 'yes'),
    '`include_geo` must be TRUE or FALSE'
  )
})

test_that('plot_sra_availability validation rejects incompatible input objects', {
  expect_gama_error(
    plot_sra_availability(tibble::tibble(x = 1)),
    'output of summarise_sra_availability\\(\\).*detected incompatible object'
  )
})

test_that('plot_sra_availability validation rejects invalid rank parameters', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_sra_availability(SRA_SUMMARY, rank = 'inputt'),
    'Invalid `rank` parameter.*Did you mean.*input'
  )
})

test_that('plot_sra_availability validation rejects invalid logical parameters', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_sra_availability(SRA_SUMMARY, abbreviate = 'yes'),
    '`abbreviate` must be TRUE or FALSE'
  )
})

test_that('plot_sra_geo validation rejects incompatible input objects', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_sra_geo(SUMMARY),
    'output of summarise_sra_availability\\(\\).*detected summarise_availability object'
  )
})

test_that('plot_sra_geo validation rejects invalid rank parameters', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_sra_geo(SRA_SUMMARY, rank = 'higest'),
    'Invalid `rank` parameter.*Did you mean.*highest'
  )
})

test_that('summarise_sra_skew validation rejects incompatible input objects', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_sra_skew(SUMMARY),
    'output of summarise_sra_availability\\(\\).*detected summarise_availability object'
  )
})

test_that('summarise_sra_skew validation reports missing SRA profile caches', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  SRA_SUMMARY <- without_attr(SRA_SUMMARY, 'sra_profile')
  expect_gama_error(
    summarise_sra_skew(SRA_SUMMARY),
    "missing required cache 'sra_profile'"
  )
})

test_that('summarise_sra_skew validation rejects invalid unit parameters', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_sra_skew(SRA_SUMMARY, unit = 'project'),
    'Invalid `unit` parameter'
  )
})

test_that('summarise_sra_skew validation rejects invalid modality class parameters', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_sra_skew(SRA_SUMMARY, class = 'not-a-class'),
    'Invalid `class` parameter'
  )
})

test_that('plot_sra_skew validation rejects incompatible input objects', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_sra_skew(SRA_SUMMARY),
    'output of summarise_sra_skew\\(\\).*detected summarise_sra_availability object'
  )
})

test_that('plot_sra_skew validation rejects invalid rank parameters', {
  SRA_SKEW_BIOPROJECT <- load_fixture('SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana')
  expect_gama_error(
    plot_sra_skew(SRA_SKEW_BIOPROJECT, rank = 'lowst'),
    'Invalid `rank` parameter.*Did you mean.*lowest'
  )
})

test_that('extract_sra_metadata validation rejects incompatible input objects', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    extract_sra_metadata(SUMMARY),
    'output of query_species\\(\\).*detected summarise_availability object'
  )
})

test_that('extract_sra_metadata validation rejects invalid logical parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_sra_metadata(RESULTS, only_geo = 'yes'),
    '`only_geo` must be TRUE or FALSE'
  )
})

test_that('extract_sra_metadata validation rejects invalid modality class parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_sra_metadata(RESULTS, class = 'not-a-class'),
    'Invalid `class` parameter'
  )
})

test_that('extract_sra_metadata validation suggests recognised modality class parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_sra_metadata(RESULTS, class = 'genomik'),
    'Invalid `class` parameter.*Did you mean.*genomic'
  )
})

test_that('extract_sra_metadata validation rejects invalid modality subclass parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_sra_metadata(RESULTS, subclass = 'not-a-subclass'),
    'Invalid `subclass` parameter'
  )
})

test_that('extract_sra_metadata validation suggests recognised modality subclass parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_sra_metadata(RESULTS, subclass = 'rnaseq'),
    'Invalid `subclass` parameter.*Did you mean.*RNA-seq'
  )
})

test_that('extract_sra_metadata validation normalises recognised modality subclass variants', {
  ontology <- getFromNamespace('.ONTOLOGY', 'GAMA')
  valid_subclass <- c(unique(unlist(lapply(ontology, names), use.names = FALSE)), 'unknown')
  subclass_parameters <- getFromNamespace('.sra_ontology_parameter_map', 'GAMA')(ontology, 'subclass')
  suggest_parameter <- getFromNamespace('.gama_suggest_parameter', 'GAMA')
  expect_identical(
    suggest_parameter('GBS', valid_subclass, variants = subclass_parameters),
    'RAD-seq'
  )
})

test_that('extract_sra_metadata validation avoids unsafe modality subclass suggestions', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  err <- expect_gama_error(
    extract_sra_metadata(RESULTS, subclass = 'GBX'),
    'Invalid `subclass` parameter'
  )
  expect_false(grepl('Did you mean', conditionMessage(err), fixed = TRUE))
})

test_that('summarise_biosample_availability validation rejects incompatible input objects', {
  expect_gama_error(
    summarise_biosample_availability(tibble::tibble(x = 1)),
    'output of query_species\\(\\).*detected incompatible object'
  )
})

test_that('summarise_biosample_availability validation reports detected GAMA object types', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_biosample_availability(SUMMARY),
    'output of query_species\\(\\).*detected summarise_availability object'
  )
})

test_that('summarise_biosample_availability validation rejects invalid logical parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_biosample_availability(RESULTS, all = 'yes'),
    '`all` must be TRUE or FALSE'
  )
})

test_that('plot_biosample_availability validation rejects incompatible input objects', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_biosample_availability(SUMMARY),
    'output of summarise_biosample_availability\\(\\).*detected summarise_availability object'
  )
})

test_that('plot_biosample_availability validation rejects invalid rank parameters', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_biosample_availability(BIO_SUMMARY, rank = 'higest'),
    'Invalid `rank` parameter.*Did you mean.*highest'
  )
})

test_that('plot_biosample_availability validation rejects invalid logical parameters', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_biosample_availability(BIO_SUMMARY, abbreviate = 'yes'),
    '`abbreviate` must be TRUE or FALSE'
  )
})

test_that('summarise_biosample_skew validation rejects incompatible input objects', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_biosample_skew(SUMMARY),
    'output of summarise_biosample_availability\\(\\).*detected summarise_availability object'
  )
})

test_that('summarise_biosample_skew validation reports missing BioSample profile caches', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  BIO_SUMMARY <- without_attr(BIO_SUMMARY, 'biosample_anatomy_profile')
  expect_gama_error(
    summarise_biosample_skew(BIO_SUMMARY),
    "missing required cache 'biosample_anatomy_profile'"
  )
})

test_that('summarise_biosample_skew validation rejects invalid anatomy-class parameters', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_biosample_skew(BIO_SUMMARY, anatomy_class = 'not-a-class'),
    'Invalid `anatomy_class` parameter'
  )
})

test_that('plot_biosample_skew validation rejects incompatible input objects', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_biosample_skew(BIO_SUMMARY),
    'output of summarise_biosample_skew\\(\\).*detected summarise_biosample_availability object'
  )
})

test_that('plot_biosample_skew validation rejects invalid rank parameters', {
  BIO_SKEW <- load_fixture('BIO_SKEW_Arabidopsis_thaliana')
  expect_gama_error(
    plot_biosample_skew(BIO_SKEW, rank = 'lowst'),
    'Invalid `rank` parameter.*Did you mean.*lowest'
  )
})

test_that('extract_biosample_metadata validation rejects incompatible input objects', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    extract_biosample_metadata(SUMMARY),
    'output of query_species\\(\\).*detected summarise_availability object'
  )
})

test_that('extract_biosample_metadata validation rejects invalid anatomy-class parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_biosample_metadata(RESULTS, anatomy_class = 'not-a-class'),
    'Invalid `anatomy_class` parameter'
  )
})

test_that('extract_biosample_metadata validation suggests recognised anatomy-class parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_biosample_metadata(RESULTS, anatomy_class = 'aeriel'),
    'Invalid `anatomy_class` parameter.*Did you mean.*aerial'
  )
})

test_that('extract_biosample_metadata validation rejects invalid anatomy-subclass parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_biosample_metadata(RESULTS, anatomy_subclass = 'not-a-subclass'),
    'Invalid `anatomy_subclass` parameter'
  )
})

test_that('extract_biosample_metadata validation suggests recognised anatomy-subclass parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_biosample_metadata(RESULTS, anatomy_subclass = 'rooot'),
    'Invalid `anatomy_subclass` parameter.*Did you mean.*root'
  )
})

test_that('extract_biosample_metadata validation rejects invalid anatomy-term parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_biosample_metadata(RESULTS, anatomy_term = 'not-a-term'),
    'Invalid `anatomy_term` parameter'
  )
})

test_that('extract_biosample_metadata validation gives broad-filter hint for invalid anatomy-term parameters', {
  RESULTS <- load_fixture('RESULTS_Arabidopsis_thaliana')
  expect_gama_error(
    extract_biosample_metadata(RESULTS, anatomy_term = 'leef'),
    'Invalid `anatomy_term` parameter.*Try `anatomy_class` or `anatomy_subclass` for broader filtering'
  )
})

test_that('summarise_interaction validation rejects incompatible SRA input objects', {
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_interaction(BIO_SUMMARY, BIO_SUMMARY),
    'output of summarise_sra_availability\\(\\).*detected summarise_biosample_availability object'
  )
})

test_that('summarise_interaction validation rejects incompatible BioSample input objects', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_interaction(SRA_SUMMARY, SRA_SUMMARY),
    'output of summarise_biosample_availability\\(\\).*detected summarise_sra_availability object'
  )
})

test_that('summarise_interaction validation reports missing SRA profile caches', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  SRA_SUMMARY <- without_attr(SRA_SUMMARY, 'sra_profile')
  expect_gama_error(
    summarise_interaction(SRA_SUMMARY, BIO_SUMMARY),
    "missing required cache 'sra_profile'"
  )
})

test_that('summarise_interaction validation reports missing BioSample profile caches', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  BIO_SUMMARY <- without_attr(BIO_SUMMARY, 'biosample_anatomy_profile')
  expect_gama_error(
    summarise_interaction(SRA_SUMMARY, BIO_SUMMARY),
    "missing required cache 'biosample_anatomy_profile'"
  )
})

test_that('summarise_interaction validation rejects invalid level parameters', {
  SRA_SUMMARY <- load_fixture('SRA_SUMMARY_Arabidopsis_thaliana')
  BIO_SUMMARY <- load_fixture('BIO_SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    summarise_interaction(SRA_SUMMARY, BIO_SUMMARY, level = 'class'),
    'Invalid `level` parameter'
  )
})

test_that('plot_interaction validation rejects incompatible input objects', {
  SUMMARY <- load_fixture('SUMMARY_Arabidopsis_thaliana')
  expect_gama_error(
    plot_interaction(SUMMARY),
    'output of summarise_interaction\\(\\).*detected summarise_availability object'
  )
})

test_that('plot_interaction validation rejects invalid value parameters', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  expect_gama_error(
    plot_interaction(INTERACTION_CLASS, value = 'counts'),
    'Invalid `value` parameter.*Did you mean.*count'
  )
})

test_that('plot_interaction validation rejects invalid rank parameters', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  expect_gama_error(
    plot_interaction(INTERACTION_CLASS, rank = 'lowst'),
    'Invalid `rank` parameter.*Did you mean.*lowest'
  )
})

test_that('plot_interaction validation rejects invalid logical parameters', {
  INTERACTION_CLASS <- load_fixture('INTERACTION_CLASS_Arabidopsis_thaliana')
  expect_gama_error(
    plot_interaction(INTERACTION_CLASS, show_values = 'yes'),
    '`show_values` must be TRUE or FALSE'
  )
})

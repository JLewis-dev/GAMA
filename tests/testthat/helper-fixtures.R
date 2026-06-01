fixture_path <- function(file) {
  if (!grepl('\\.rds$', file)) file <- paste0(file, '.rds')
  testthat::test_path('object-fixtures', file)
}

load_fixture <- function(file) {
  path <- fixture_path(file)
  if (!file.exists(path)) stop('Fixture not found: ', path, call. = FALSE)
  readRDS(path)
}

expect_gdt_tbl <- function(x) {
  testthat::expect_s3_class(x, 'gdt_tbl')
  testthat::expect_false(is.null(attr(x, 'query_info')))
}

expect_gama_columns <- function(x, cols) {
  testthat::expect_true(all(cols %in% names(x)))
}

expect_skew_id_recovery <- function(x, unit = NULL, class_col = 'class') {
  RECOVERY <- attr(x, 'skew_id_recovery', exact = TRUE)
  expected_cols <- c('species', 'unit', class_col, 'records', 'included', 'excluded', 'id_recovery_prop')
  testthat::expect_false(is.null(RECOVERY))
  expect_gama_columns(RECOVERY, expected_cols)
  testthat::expect_false(any(is.na(RECOVERY$records)))
  testthat::expect_false(any(is.na(RECOVERY$included)))
  testthat::expect_false(any(is.na(RECOVERY$excluded)))
  testthat::expect_true(all(RECOVERY$records >= RECOVERY$included))
  testthat::expect_equal(RECOVERY$excluded, RECOVERY$records - RECOVERY$included)
  testthat::expect_true(all(is.na(RECOVERY$id_recovery_prop) | (RECOVERY$id_recovery_prop >= 0 & RECOVERY$id_recovery_prop <= 1)))
  if (!is.null(unit)) testthat::expect_true(all(RECOVERY$unit == unit))
}

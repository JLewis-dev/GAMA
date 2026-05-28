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

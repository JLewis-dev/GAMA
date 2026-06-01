devtools::load_all()

fixture_dir <- file.path('tests', 'testthat', 'object-fixtures')
rows_per_group <- 100000L

dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

save_fixture <- function(x, name) {
  saveRDS(x, file.path(fixture_dir, paste0(name, '.rds')), compress = 'xz')
  invisible(x)
}

compact_biosample_profile_fixture <- function(x, n = rows_per_group) {
  if (is.null(x)) return(NULL)
  group_cols <- intersect(
    c('species', 'anatomy_class', 'anatomy_subclass', 'anatomy_term'),
    names(x)
  )
  if (!length(group_cols)) return(dplyr::slice_head(x, n = n))
  x |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::slice_head(n = n) |>
    dplyr::ungroup()
}

compact_biosample_summary_fixture <- function(x, n = rows_per_group) {
  attr(x, 'biosample_anatomy_profile') <- compact_biosample_profile_fixture(attr(x, 'biosample_anatomy_profile', exact = TRUE), n = n)
  attr(x, 'biosample_canonical_profile') <- compact_biosample_profile_fixture(attr(x, 'biosample_canonical_profile', exact = TRUE), n = n)
  attr(x, 'biosample_anatomy_profile_info') <- NULL
  attr(x, 'biosample_canonical_profile_info') <- NULL
  x
}

compact_biosample_skew_fixture <- function(x, n = rows_per_group) {
  attr(x, 'biosample_anatomy_profile') <- compact_biosample_profile_fixture(attr(x, 'biosample_anatomy_profile', exact = TRUE), n = n)
  x
}

compact_biosample_extract_fixture <- function(x, n = rows_per_group) {
  OUT <- x |>
    dplyr::group_by(
      .data$anatomy_class_profile,
      .data$anatomy_subclass_profile,
      .data$anatomy_class,
      .data$anatomy_subclass
    ) |>
    dplyr::slice_head(n = n) |>
    dplyr::ungroup()
  attr(OUT, 'query_info') <- attr(x, 'query_info', exact = TRUE)
  attr(OUT, 'gama_object') <- attr(x, 'gama_object', exact = TRUE)
  class(OUT) <- class(x)
  OUT
}

RESULTS_Arabidopsis_thaliana <- query_species('Arabidopsis thaliana')

SUMMARY_Arabidopsis_thaliana <- summarise_availability(RESULTS_Arabidopsis_thaliana)

ASM_SUMMARY_Arabidopsis_thaliana <- summarise_assembly_availability(RESULTS_Arabidopsis_thaliana)
ASM_Arabidopsis_thaliana <- extract_assembly_metadata(RESULTS_Arabidopsis_thaliana)

SRA_SUMMARY_Arabidopsis_thaliana <- summarise_sra_availability(RESULTS_Arabidopsis_thaliana, all = TRUE, include_geo = TRUE)
SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana <- summarise_sra_skew(SRA_SUMMARY_Arabidopsis_thaliana, unit = 'bioproject')
SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana <- summarise_sra_skew(SRA_SUMMARY_Arabidopsis_thaliana, unit = 'biosample')
SRA_Arabidopsis_thaliana <- extract_sra_metadata(RESULTS_Arabidopsis_thaliana)

BIO_SUMMARY_Arabidopsis_thaliana <- summarise_biosample_availability(RESULTS_Arabidopsis_thaliana, all = TRUE)
BIO_SKEW_Arabidopsis_thaliana <- summarise_biosample_skew(BIO_SUMMARY_Arabidopsis_thaliana)
BIO_Arabidopsis_thaliana <- extract_biosample_metadata(RESULTS_Arabidopsis_thaliana)

INTERACTION_CLASS_Arabidopsis_thaliana <- summarise_interaction(SRA_SUMMARY_Arabidopsis_thaliana, BIO_SUMMARY_Arabidopsis_thaliana, level = 'anatomy_class')
INTERACTION_SUBCLASS_Arabidopsis_thaliana <- summarise_interaction(SRA_SUMMARY_Arabidopsis_thaliana, BIO_SUMMARY_Arabidopsis_thaliana, level = 'anatomy_subclass')

BIO_SUMMARY_Arabidopsis_thaliana <- compact_biosample_summary_fixture(BIO_SUMMARY_Arabidopsis_thaliana, n = rows_per_group)
BIO_SKEW_Arabidopsis_thaliana <- compact_biosample_skew_fixture(BIO_SKEW_Arabidopsis_thaliana, n = rows_per_group)
BIO_Arabidopsis_thaliana <- compact_biosample_extract_fixture(BIO_Arabidopsis_thaliana, n = rows_per_group)

save_fixture(RESULTS_Arabidopsis_thaliana, 'RESULTS_Arabidopsis_thaliana')
save_fixture(SUMMARY_Arabidopsis_thaliana, 'SUMMARY_Arabidopsis_thaliana')
save_fixture(ASM_SUMMARY_Arabidopsis_thaliana, 'ASM_SUMMARY_Arabidopsis_thaliana')
save_fixture(ASM_Arabidopsis_thaliana, 'ASM_Arabidopsis_thaliana')
save_fixture(SRA_SUMMARY_Arabidopsis_thaliana, 'SRA_SUMMARY_Arabidopsis_thaliana')
save_fixture(SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana, 'SRA_SKEW_BIOPROJECT_Arabidopsis_thaliana')
save_fixture(SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana, 'SRA_SKEW_BIOSAMPLE_Arabidopsis_thaliana')
save_fixture(SRA_Arabidopsis_thaliana, 'SRA_Arabidopsis_thaliana')
save_fixture(BIO_SUMMARY_Arabidopsis_thaliana, 'BIO_SUMMARY_Arabidopsis_thaliana')
save_fixture(BIO_SKEW_Arabidopsis_thaliana, 'BIO_SKEW_Arabidopsis_thaliana')
save_fixture(BIO_Arabidopsis_thaliana, 'BIO_Arabidopsis_thaliana')
save_fixture(INTERACTION_CLASS_Arabidopsis_thaliana, 'INTERACTION_CLASS_Arabidopsis_thaliana')
save_fixture(INTERACTION_SUBCLASS_Arabidopsis_thaliana, 'INTERACTION_SUBCLASS_Arabidopsis_thaliana')

query_info <- attr(RESULTS_Arabidopsis_thaliana, 'query_info', exact = TRUE)
if (is.null(query_info) || is.null(query_info$query_time_utc) || is.null(query_info$tool_version)) {
  stop('Missing query provenance on RESULTS_Arabidopsis_thaliana.', call. = FALSE)
}

Arabidopsis_thaliana_fixture_manifest <- tibble::tibble(
  object = c(
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
)

Arabidopsis_thaliana_fixture_manifest$file <- paste0(Arabidopsis_thaliana_fixture_manifest$object, '.rds')
Arabidopsis_thaliana_fixture_manifest$size_bytes <- file.info(file.path(fixture_dir, Arabidopsis_thaliana_fixture_manifest$file))$size
Arabidopsis_thaliana_fixture_manifest$query_time_utc <- query_info$query_time_utc
Arabidopsis_thaliana_fixture_manifest$tool_version <- query_info$tool_version

save_fixture(Arabidopsis_thaliana_fixture_manifest, 'Arabidopsis_thaliana_fixture_manifest')

print(Arabidopsis_thaliana_fixture_manifest)

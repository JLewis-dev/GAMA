devtools::load_all()

example_dir <- file.path('data-raw', 'examples')
figure_dir <- file.path('man', 'figures')

dir.create(example_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

save_example <- function(x, name) {
  saveRDS(x, file.path(example_dir, paste0(name, '.rds')), compress = 'xz')
  invisible(x)
}

save_example_plot <- function(x, name) {
  ggplot2::ggsave(
    filename = file.path(figure_dir, paste0(name, '.png')),
    plot = x,
    width = 10,
    height = 5,
    units = 'in',
    dpi = 200,
    bg = 'white',
    limitsize = FALSE
  )
  invisible(x)
}

SPECIES_example <- c(
  'Arabidopsis thaliana',
  'Glycine max',
  'Phaseolus vulgaris',
  'Vigna radiata'
)

RESULTS_example <- query_species(SPECIES_example)

SUMMARY_example <- summarise_availability(RESULTS_example)

SRA_SUMMARY_example <- summarise_sra_availability(RESULTS_example)

SKEW_example <- summarise_sra_skew(SRA_SUMMARY_example, class = 'transcriptomic')

BIO_SUMMARY_example <- summarise_biosample_availability(RESULTS_example)

INTERACTION_example <- summarise_interaction(
  SRA_SUMMARY_example,
  BIO_SUMMARY_example,
  level = 'anatomy_subclass',
  species = 'Arabidopsis thaliana'
)

save_example(SUMMARY_example, 'SUMMARY_example')
save_example(SRA_SUMMARY_example, 'SRA_SUMMARY_example')
save_example(SKEW_example, 'SKEW_example')
save_example(INTERACTION_example, 'INTERACTION_example')

Data_richness_plot <- plot_availability(SUMMARY_example)
Modality_plot <- plot_sra_availability(SRA_SUMMARY_example)
Skew_plot <- plot_sra_skew(SKEW_example)
Interaction_plot <- plot_interaction(INTERACTION_example)

if (!inherits(Interaction_plot, 'ggplot') && is.list(Interaction_plot)) {
  Interaction_plot <- Interaction_plot[[1L]]
}

save_example_plot(Data_richness_plot, 'Data_richness_plot')
save_example_plot(Modality_plot, 'Modality_plot')
save_example_plot(Skew_plot, 'Skew_plot')
save_example_plot(Interaction_plot, 'Interaction_plot')

query_info <- attr(RESULTS_example, 'query_info', exact = TRUE)
if (is.null(query_info) || is.null(query_info$query_time_utc) || is.null(query_info$tool_version)) {
  stop('Missing query provenance on RESULTS_example.', call. = FALSE)
}

example_manifest <- tibble::tibble(
  object = c(
    'SUMMARY_example',
    'SRA_SUMMARY_example',
    'SKEW_example',
    'INTERACTION_example'
  )
)

example_manifest$file <- paste0(example_manifest$object, '.rds')
example_manifest$size_bytes <- file.info(file.path(example_dir, example_manifest$file))$size
example_manifest$query_time_utc <- query_info$query_time_utc
example_manifest$tool_version <- query_info$tool_version

save_example(example_manifest, 'example_manifest')

print(example_manifest)

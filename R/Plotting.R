# PLOTTING ====================================================================

#' Plot data richness
#'
#' Visualises species-level data richness using the composite scores generated
#' during the query phase. Assembly, SRA, and BioSample contributions are shown
#' as stacked bar segments for each species.
#'
#' @param SUMMARY A tibble returned by [summarise_availability()], containing
#' species-level composite scores and component values (`A`, `S`, `B`).
#' @param rank Ordering of species on the x-axis. One of `highest` (default),
#' `lowest`, `A-Z`, `Z-A`, or `input`.
#' @param abbreviate Logical; if `TRUE` (default), abbreviate species names to
#' `G. species` style labels.
#' @param theme_fn A ggplot2 theme function (e.g. [ggplot2::theme_minimal()]).
#' @param colours Named character vector of fill colours for the `Assembly`,
#' `SRA`, and `BioSample` segments.
#'
#' @return A ggplot object showing stacked bar segments for each species.
#'
#' @seealso [summarise_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SUMMARY <- summarise_availability(RESULTS)
#' print(SUMMARY)
#' plot_availability(SUMMARY)
#' }
#' @export
plot_availability <- function(
SUMMARY,
rank        = 'highest',
abbreviate  = TRUE,
theme_fn    = ggplot2::theme_minimal,
colours     = c(
Assembly  = '#E3F9F2',
SRA       = '#A8E6CF',
BioSample = '#56C5A8'
)
) {
  rank <- .gama_validate_parameters(rank, 'rank', .gama_rank_parameters, multiple = FALSE, allow_null = FALSE)
  abbreviate <- .gama_validate_logical_parameter(abbreviate, 'abbreviate')
  req_cols <- c('species', 'score', 'A', 'S', 'B')
  SUMMARY <- .gama_require_output(SUMMARY, 'summarise_availability', required_cols = req_cols)
  SUMMARY$species_label <- if (abbreviate) .shorten_species(SUMMARY$species) else SUMMARY$species
  SUMMARY$species_label <- factor(
  SUMMARY$species_label,
  levels = {
    if (rank == 'input') {
      SUMMARY$species_label
    } else if (rank == 'highest') {
      SUMMARY$species_label[order(SUMMARY$score, decreasing = TRUE)]
    } else if (rank == 'lowest') {
      SUMMARY$species_label[order(SUMMARY$score, decreasing = FALSE)]
    } else if (rank == 'A-Z') {
      sort(SUMMARY$species_label, decreasing = FALSE)
    } else {
      sort(SUMMARY$species_label, decreasing = TRUE)
    }
  }
  )
  PLOT <- SUMMARY |>
  dplyr::select(species_label, score, A, S, B) |>
  tidyr::pivot_longer(
  cols      = c(A, S, B),
  names_to  = 'database',
  values_to = 'score_component'
  ) |>
  dplyr::mutate(
  database = dplyr::case_when(
  database == 'A' ~ 'Assembly',
  database == 'S' ~ 'SRA',
  database == 'B' ~ 'BioSample'
  ),
  database = factor(database, levels = c('Assembly', 'SRA', 'BioSample')),
  segment  = score_component
  )
  ggplot2::ggplot(
  data = PLOT,
  ggplot2::aes(x = species_label, y = segment, fill = database)
  ) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::scale_fill_manual(values = colours) +
  ggplot2::labs(
  x    = NULL,
  y    = 'Data richness score',
  fill = 'Database:'
  ) +
  theme_fn(base_size = 13) +
  ggplot2::theme(
  plot.margin     = ggplot2::margin(10, 40, 10, 10), legend.position = 'top',
  axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, face  = 'italic'),
  axis.title.y    = ggplot2::element_text(margin = ggplot2::margin(r = 10))
  )
}

#' Plot Assembly composition
#'
#' Visualises species-level Assembly composition using stacked horizontal bars.
#' Each bar shows the proportional contribution of recognised assembly levels
#' (`complete`, `chromosome`, `scaffold`, `contig`), with total Assembly
#' accession counts labelled.
#'
#' Operates on the wide-format summary returned by
#' [summarise_assembly_availability()].
#'
#' @param ASSEMBLY A wide-format Assembly summary table returned by
#' [summarise_assembly_availability()].
#' @param species `NULL` (default) to plot all species, or a character vector
#' of species to include.
#' @param rank Ordering of species. One of `highest` (default), `lowest`,
#' `A-Z`, `Z-A`, or `input`.
#' @param abbreviate Logical; if `TRUE` (default), abbreviate species names.
#' @param theme_fn A ggplot2 theme function.
#' @param colours Named character vector of assembly level colours.
#'
#' @return A ggplot object showing proportional assembly level profiles across
#' species.
#'
#' @seealso [summarise_assembly_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' ASM_SUMMARY <- summarise_assembly_availability(RESULTS)
#' plot_assembly_availability(ASM_SUMMARY)
#' }
#' @export
plot_assembly_availability <- function(
    ASSEMBLY,
    species    = NULL,
    rank       = 'highest',
    abbreviate = TRUE,
    theme_fn   = ggplot2::theme_minimal,
    colours    = c(
      complete   = '#E3F9F2',
      chromosome = '#A8E6CF',
      scaffold   = '#56C5A8',
      contig     = '#2FA083'
    )
) {
  rank <- .gama_validate_parameters(rank, 'rank', .gama_rank_parameters, multiple = FALSE, allow_null = FALSE)
  abbreviate <- .gama_validate_logical_parameter(abbreviate, 'abbreviate')
  LEVELS <- .ASSEMBLY_CLASSES
  req_cols <- c('species', 'Assembly', LEVELS, 'best_n50')
  ASSEMBLY <- .gama_require_output(ASSEMBLY, 'summarise_assembly_availability', required_cols = req_cols)
  CORE <- ASSEMBLY |>
    dplyr::select(species, Assembly, dplyr::all_of(LEVELS))
  if (!is.null(species)) {
    species <- as.character(species)
    species <- species[!is.na(species) & nzchar(species)]
    species <- unique(species)
    missing_sp <- setdiff(species, CORE$species)
    if (length(missing_sp)) .gama_warn(sprintf('Requested species not found in input `ASSEMBLY`: %s. Dropping.', paste(missing_sp, collapse = ', ')))
    CORE <- CORE |> dplyr::filter(.data$species %in% .env$species)
    if (!nrow(CORE)) .gama_stop('No matching species found.')
  }
  TOTAL <- CORE[, c('species', 'Assembly')]
  if (rank == 'highest') {
    TOTAL <- TOTAL[order(TOTAL$Assembly, decreasing = TRUE), ]
  } else if (rank == 'lowest') {
    TOTAL <- TOTAL[order(TOTAL$Assembly, decreasing = FALSE), ]
  } else if (rank == 'A-Z') {
    TOTAL <- TOTAL[order(TOTAL$species, decreasing = FALSE), ]
  } else if (rank == 'Z-A') {
    TOTAL <- TOTAL[order(TOTAL$species, decreasing = TRUE), ]
  } else if (rank == 'input') {
    TOTAL <- TOTAL
  }
  species_order <- TOTAL$species
  LONG <- CORE |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(LEVELS),
      names_to  = 'level',
      values_to = 'count'
    ) |>
    dplyr::mutate(
      prop = dplyr::if_else(.data$Assembly > 0, .data$count / .data$Assembly, 0)
    )
  LONG$level <- factor(LONG$level, levels = LEVELS)
  LONG$species_label <- if (abbreviate) {
    .shorten_species(LONG$species)
  } else {
    LONG$species
  }
  levels_in_plot <- if (abbreviate) {
    .shorten_species(species_order)
  } else {
    species_order
  }
  LONG$species_label <- factor(LONG$species_label, levels = rev(levels_in_plot))
  p <- ggplot2::ggplot(
    LONG,
    ggplot2::aes(x = species_label, y = prop, fill = level)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::coord_flip(clip = 'off') +
    ggplot2::labs(
      x    = NULL,
      y    = 'Proportion of Assembly accessions',
      fill = 'Level:'
    ) +
    theme_fn(base_size = 13) +
    ggplot2::theme(
      plot.margin        = ggplot2::margin(10, 40, 10, 10),
      legend.position    = 'top',
      axis.text.y        = ggplot2::element_text(size = 11, face = 'italic'),
      axis.text.x        = ggplot2::element_text(size = 11),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_blank()
    )
  totals_df <- TOTAL
  totals_df$species_label <- if (abbreviate) {
    .shorten_species(totals_df$species)
  } else {
    totals_df$species
  }
  totals_df$species_label <- factor(
    totals_df$species_label,
    levels = levels(LONG$species_label)
  )
  p +
    ggplot2::geom_text(
      data = totals_df,
      ggplot2::aes(x = species_label, y = 1.005, label = Assembly),
      hjust = 0,
      size  = 3.5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
}

#' Plot SRA modality composition
#'
#' Visualises species-level SRA modality composition using stacked horizontal
#' bars. Each bar shows the proportional contribution of major experimental
#' classes (`genomic`, `transcriptomic`, `epigenomic`, `chromatin`, `other`,
#' `unknown`), with total SRA counts labelled.
#'
#' Operates on the wide-format summary returned by
#' [summarise_sra_availability()].
#'
#' @param SRA A wide-format SRA summary table returned by
#' [summarise_sra_availability()].
#' @param species `NULL` (default) to plot all species, or a character vector
#' of species to include.
#' @param rank Ordering of species. One of `highest` (default), `lowest`,
#' `A-Z`, `Z-A`, or `input`.
#' @param abbreviate Logical; if `TRUE` (default), abbreviate species names.
#' @param theme_fn A ggplot2 theme function.
#' @param colours Named character vector of class colours.
#'
#' @return A ggplot object showing proportional SRA modality profiles across
#' species.
#'
#' @seealso [summarise_sra_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' plot_sra_availability(SRA_SUMMARY)
#' }
#' @export
plot_sra_availability <- function(
    SRA,
    species    = NULL,
    rank       = 'highest',
    abbreviate = TRUE,
    theme_fn   = ggplot2::theme_minimal,
    colours    = c(
      genomic        = '#E3F9F2',
      transcriptomic = '#A8E6CF',
      epigenomic     = '#56C5A8',
      chromatin      = '#2FA083',
      other          = '#166A55',
      unknown        = '#BDBDBD'
    )
) {
  rank <- .gama_validate_parameters(rank, 'rank', .gama_rank_parameters, multiple = FALSE, allow_null = FALSE)
  abbreviate <- .gama_validate_logical_parameter(abbreviate, 'abbreviate')
  req_cols <- c('species', 'SRA', 'genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
  SRA <- .gama_require_output(SRA, 'summarise_sra_availability', required_cols = req_cols)
  CORE <- SRA |>
    dplyr::select(
      species, SRA, genomic, transcriptomic, epigenomic, chromatin, other,
      unknown
    )
  if (!is.null(species)) {
    missing_sp <- setdiff(species, CORE$species)
    if (length(missing_sp)) .gama_warn(sprintf('Requested species not found in input `SRA`: %s. Dropping.', paste(missing_sp, collapse = ', ')))
    CORE <- CORE |> dplyr::filter(.data$species %in% .env$species)
    if (!nrow(CORE)) .gama_stop('No matching species found.')
  }
  TOTAL <- CORE[, c('species', 'SRA')]
  if (rank == 'highest') {
    TOTAL <- TOTAL[order(TOTAL$SRA, decreasing = TRUE), ]
  } else if (rank == 'lowest') {
    TOTAL <- TOTAL[order(TOTAL$SRA, decreasing = FALSE), ]
  } else if (rank == 'A-Z') {
    TOTAL <- TOTAL[order(TOTAL$species, decreasing = FALSE), ]
  } else if (rank == 'Z-A') {
    TOTAL <- TOTAL[order(TOTAL$species, decreasing = TRUE), ]
  } else if (rank == 'input') {
    TOTAL <- TOTAL
  }
  species_order <- TOTAL$species
  LONG <- CORE |>
    tidyr::pivot_longer(
      cols      = c(
        'genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other',
        'unknown'
      ),
      names_to  = 'class',
      values_to = 'count'
    ) |>
    dplyr::group_by(species) |>
    dplyr::mutate(prop = count / sum(count)) |>
    dplyr::ungroup()
  LONG$class <- factor(
    LONG$class,
    levels = c(
      'genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other',
      'unknown'
    )
  )
  LONG$species_label <- if (abbreviate) {
    .shorten_species(LONG$species)
  } else {
    LONG$species
  }
  levels_in_plot <- if (abbreviate) {
    .shorten_species(species_order)
  } else {
    species_order
  }
  LONG$species_label <- factor(LONG$species_label, levels = rev(levels_in_plot))
  p <- ggplot2::ggplot(
    LONG,
    ggplot2::aes(x = species_label, y = prop, fill = class)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::coord_flip(clip = 'off') +
    ggplot2::labs(
      x    = NULL,
      y    = 'Proportion of experiments',
      fill = 'Class:'
    ) +
    theme_fn(base_size = 13) +
    ggplot2::theme(
      plot.margin        = ggplot2::margin(10, 40, 10, 10),
      legend.position    = 'top',
      axis.text.y        = ggplot2::element_text(size = 11, face = 'italic'),
      axis.text.x        = ggplot2::element_text(size = 11),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_blank()
    )
  totals_df <- TOTAL
  totals_df$species_label <- if (abbreviate) {
    .shorten_species(totals_df$species)
  } else {
    totals_df$species
  }
  totals_df$species_label <- factor(
    totals_df$species_label,
    levels = levels(LONG$species_label)
  )
  p +
    ggplot2::geom_text(
      data = totals_df,
      ggplot2::aes(x = species_label, y = 1.005, label = SRA),
      hjust = 0,
      size  = 3.5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
}

#' Plot SRA modality GEO linkage overlay
#'
#' Visualises GEO linkage across SRA modality classes, excluding `genomic`.
#' Each modality is displayed as a 100% bar with a translucent background
#' fill, while the coloured segment represents the GEO-linked fraction.
#' Labels show `GEO-linked / Total` for each modality.
#'
#' Operates on the direct output of [summarise_sra_availability()].
#' GEO-linked counts are derived from the cached `attr(SRA, 'sra_profile')`,
#' which stores per-experiment `geo_linked` values regardless of
#' `include_geo`.
#'
#' If `species` is `NULL`, plots are generated for all species in the table. A
#' single species returns a ggplot object; multiple species return a named
#' list of ggplot objects.
#'
#' @param SRA A wide-format SRA summary table returned by
#' [summarise_sra_availability()].
#' @param species `NULL` (default) for all species, or a character vector of
#' species to plot.
#' @param rank Ordering of species (when `species = NULL`), or ordering applied
#' to the requested species vector. One of `highest`, `lowest`, `A-Z`, `Z-A`,
#' or `input`.
#' @param theme_fn A ggplot2 theme function.
#' @param colours Named character vector of class colours.
#' @param alpha_vals Named numeric vector giving alpha values for GEO-linked
#' vs not GEO-linked segments.
#'
#' @return A ggplot object (single species) or a named list of ggplot objects
#' (multiple species).
#'
#' @seealso [summarise_sra_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' plot_sra_geo(SRA_SUMMARY, species = 'Vigna vexillata')
#' }
#' @export
plot_sra_geo <- function(
    SRA,
    species    = NULL,
    rank       = 'highest',
    theme_fn   = ggplot2::theme_minimal,
    colours    = c(
      transcriptomic = '#A8E6CF',
      epigenomic     = '#56C5A8',
      chromatin      = '#2FA083',
      other          = '#166A55',
      unknown        = '#BDBDBD'
    ),
    alpha_vals = c(`Not GEO-linked` = 0.25, `GEO-linked` = 1)
) {
  rank <- .gama_validate_parameters(rank, 'rank', .gama_rank_parameters, multiple = FALSE, allow_null = FALSE)
  GEO_CLASSES <- c(
    'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown'
  )
  FIXED_ORDER <- rev(GEO_CLASSES)
  SRA <- .gama_require_output(
    SRA,
    'summarise_sra_availability',
    required_cols = c('species', GEO_CLASSES)
  )
  prof <- .gama_require_cache(
    SRA,
    attr_name = 'sra_profile',
    required_cols = c('species', 'class', 'geo_linked'),
    source = 'summarise_sra_availability'
  )
  GEO_COUNTS <- prof |>
    dplyr::filter(.data$class %in% GEO_CLASSES, .data$geo_linked) |>
    dplyr::count(species, class, name = 'linked') |>
    tidyr::complete(
      species = unique(SRA$species),
      class   = GEO_CLASSES,
      fill    = list(linked = 0L)
    )
  if (is.null(species)) {
    sp <- SRA$species
    if (rank != 'input') {
      ord <- SRA[, c('species', 'SRA')]
      if (rank == 'highest') ord <- ord[order(ord$SRA, decreasing = TRUE), ]
      if (rank == 'lowest') ord <- ord[order(ord$SRA, decreasing = FALSE), ]
      if (rank == 'A-Z') ord <- ord[order(ord$species, decreasing = FALSE), ]
      if (rank == 'Z-A') ord <- ord[order(ord$species, decreasing = TRUE), ]
      sp <- ord$species
    }
  } else {
    sp <- intersect(species, SRA$species)
    missing_sp <- setdiff(species, SRA$species)
    if (length(missing_sp)) {
      .gama_warn(sprintf(
        'Requested species not found in input `SRA`: %s. Dropping.',
        paste(missing_sp, collapse = ', ')
      ))
    }
    if (!length(sp)) .gama_stop('No matching species found.')
    if (rank == 'A-Z') sp <- sort(sp, decreasing = FALSE)
    if (rank == 'Z-A') sp <- sort(sp, decreasing = TRUE)
  }
  .plot_one <- function(one_species) {
    row <- SRA |> dplyr::filter(.data$species == one_species)
    if (!nrow(row)) .gama_stop('No matching species found.')
    total <- as.integer(row[1, GEO_CLASSES, drop = TRUE])
    linked <- GEO_COUNTS |>
      dplyr::filter(.data$species == one_species) |>
      dplyr::arrange(match(.data$class, GEO_CLASSES)) |>
      dplyr::pull(linked)
    total[is.na(total)] <- 0L
    linked[is.na(linked)] <- 0L
    prop_linked <- ifelse(total > 0, linked / total, 0)
    prop_not <- ifelse(total > 0, 1 - prop_linked, 0)
    SUM <- tibble::tibble(
      class       = factor(GEO_CLASSES, levels = FIXED_ORDER),
      total       = total,
      linked      = linked,
      label       = paste0(linked, '/', total),
      prop_linked = prop_linked,
      prop_not    = prop_not
    )
    LONG <- dplyr::bind_rows(
      SUM |> dplyr::transmute(
        class,
        status = 'Not GEO-linked',
        prop   = prop_not
      ),
      SUM |> dplyr::transmute(
        class,
        status = 'GEO-linked',
        prop   = prop_linked
      )
    )
    LONG$class <- factor(LONG$class, levels = FIXED_ORDER)
    LONG$status <- factor(
      LONG$status,
      levels = c('Not GEO-linked', 'GEO-linked')
    )
    ggplot2::ggplot(
      LONG,
      ggplot2::aes(x = class, y = prop, fill = class, alpha = status)
    ) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::geom_text(
        data = SUM,
        ggplot2::aes(x = class, y = 1.005, label = label),
        size = 3.5,
        vjust = 0,
        hjust = 0,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
      ggplot2::scale_alpha_manual(values = alpha_vals) +
      ggplot2::guides(
        fill  = ggplot2::guide_legend(reverse = TRUE),
        alpha = ggplot2::guide_legend(reverse = FALSE)
      ) +
      ggplot2::coord_flip(clip = 'off') +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 1.03),
        breaks = seq(0, 1, 0.25),
        labels = function(x) sprintf('%.2f', x)
      ) +
      ggplot2::labs(
        title = one_species,
        x     = NULL,
        y     = 'Proportion GEO-linked',
        fill  = 'Class:',
        alpha = NULL
      ) +
      theme_fn(base_size = 13) +
      ggplot2::theme(
        plot.margin     = ggplot2::margin(10, 60, 10, 10),
        legend.position = 'top',
        plot.title      = ggplot2::element_text(face = 'italic')
      )
  }
  if (length(sp) == 1) return(.plot_one(sp))
  OUT <- lapply(sp, .plot_one)
  names(OUT) <- sp
  OUT
}

#' Plot SRA replication skew
#'
#' Visualises replication skew using the summary output of
#' [summarise_sra_skew()]. Boxplots are drawn from pre-computed summary
#' statistics on a log10 y-axis.
#'
#' Optional labels show `eff=<x> (n=<y>)` for each species. When a cached
#' UID-level profile is available as `attr(SKEW, 'sra_profile')`, per-unit
#' points (one point per BioProject/BioSample) can be overlaid with
#' horizontal jitter.
#'
#' @param SKEW A tibble returned by [summarise_sra_skew()].
#' @param species `NULL` (default) to include all species, or a character
#'   vector of species to plot.
#' @param rank Ordering of species. One of `highest` (default), `lowest`,
#'   `A-Z`, `Z-A`, or `input`.
#' @param abbreviate Logical; if `TRUE` (default), abbreviate species labels.
#' @param theme_fn A ggplot2 theme function.
#' @param colours Named character vector for box fill and line colours.
#' @param show_points Logical; if `TRUE` (default), overlay per-unit points
#'   when cached profiles are available.
#' @param point_colour Character scalar giving the colour of overlaid points.
#' @param point_alpha Numeric alpha value for overlaid points.
#' @param show_labels Logical; if `TRUE` (default), label each box with
#'   `eff` and `n`.
#' @param label_digits Integer; decimal places for `eff` labels.
#'
#' @return A ggplot object.
#'
#' @seealso [summarise_sra_skew()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' SKEW <- summarise_sra_skew(SRA_SUMMARY, class = 'transcriptomic')
#' plot_sra_skew(SKEW)
#' }
#'
#' @export
plot_sra_skew <- function(
    SKEW,
    species      = NULL,
    rank         = 'highest',
    abbreviate   = TRUE,
    theme_fn     = ggplot2::theme_minimal,
    colours      = c(box = '#E3F9F2', line = 'black'),
    show_points  = TRUE,
    point_colour = '#56C5A8',
    point_alpha  = 0.25,
    show_labels  = TRUE,
    label_digits = 1L
) {
  rank <- .gama_validate_parameters(rank, 'rank', .gama_rank_parameters, multiple = FALSE, allow_null = FALSE)
  abbreviate <- .gama_validate_logical_parameter(abbreviate, 'abbreviate')
  show_points <- .gama_validate_logical_parameter(show_points, 'show_points')
  show_labels <- .gama_validate_logical_parameter(show_labels, 'show_labels')
  SKEW <- .gama_require_output(
    SKEW,
    'summarise_sra_skew',
    required_cols = c('species', 'class', 'min', 'q25', 'med', 'q75', 'max', 'eff')
  )
  unit_col <- if ('BioProject' %in% names(SKEW)) {
    'BioProject'
  } else if ('BioSample' %in% names(SKEW)) {
    'BioSample'
  } else {
    .gama_input_error(
      'summarise_sra_skew',
      detected = 'summarise_sra_skew',
      detail = 'missing `BioProject` or `BioSample` column.'
    )
  }
  CORE0 <- SKEW |> dplyr::select(species, dplyr::all_of(unit_col), class, min, q25, med, q75, max, eff)
  if (!is.null(species)) {
    missing_sp <- setdiff(species, CORE0$species)
    if (length(missing_sp)) .gama_warn(sprintf('Requested species not found in input `SKEW`: %s. Dropping.', paste(missing_sp, collapse = ', ')))
    CORE0 <- CORE0 |> dplyr::filter(.data$species %in% .env$species)
    if (!nrow(CORE0)) .gama_stop('No matching species found in `SKEW` for the requested filter.')
  }
  class_vals <- unique(CORE0$class)
  if (length(class_vals) != 1L) .gama_stop('plot_sra_skew() expects a single class per plot.')
  class_val <- class_vals[[1]]
  dropped <- unique(CORE0$species[is.na(CORE0$max)])
  CORE <- CORE0 |> dplyr::filter(!is.na(.data$max))
  if (!nrow(CORE)) {
    .gama_stop(sprintf('No data to plot for class \'%s\': all selected species have zero records for this class.', class_val))
  }
  if (length(dropped) > 0L) {
    .gama_msg(sprintf(
      'Dropping %d species with no \'%s\' data: %s.',
      length(dropped),
      class_val,
      paste(dropped, collapse = ', ')
    ))
  }
  if (rank == 'highest') CORE <- CORE[order(CORE$eff, decreasing = TRUE), ]
  if (rank == 'lowest')  CORE <- CORE[order(CORE$eff, decreasing = FALSE), ]
  if (rank == 'A-Z')     CORE <- CORE[order(CORE$species, decreasing = FALSE), ]
  if (rank == 'Z-A')     CORE <- CORE[order(CORE$species, decreasing = TRUE), ]
  CORE$species_label <- if (abbreviate) .shorten_species(CORE$species) else CORE$species
  CORE$species_label <- factor(CORE$species_label, levels = CORE$species_label)
  x_labs <- levels(CORE$species_label)
  CORE$x <- as.numeric(CORE$species_label)
  class_title <- if (identical(tolower(class_val), 'all')) 'Experiments' else paste0(toupper(substr(class_val, 1, 1)), substr(class_val, 2, nchar(class_val)))
  y_prefix <- if (identical(tolower(class_val), 'all')) 'Experiments' else paste0(class_title, ' experiments')
  y_lab <- paste0(y_prefix, ' per ', unit_col, ' (log10)')
  box_w <- 0.7
  cap_w <- 0.22
  BOX <- CORE |> dplyr::transmute(x, xmin = x - box_w / 2, xmax = x + box_w / 2, ymin = q25, ymax = q75, med = med, min = min, max = max)
  p <- ggplot2::ggplot(CORE, ggplot2::aes(x = x)) +
    ggplot2::geom_rect(data = BOX, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = colours[['box']], colour = NA, inherit.aes = FALSE) +
    ggplot2::geom_segment(data = BOX, ggplot2::aes(x = x, xend = x, y = ymax, yend = max), linewidth = 0.35, colour = colours[['line']], inherit.aes = FALSE) +
    ggplot2::geom_segment(data = BOX, ggplot2::aes(x = x, xend = x, y = ymin, yend = min), linewidth = 0.35, colour = colours[['line']], inherit.aes = FALSE) +
    ggplot2::geom_segment(data = BOX, ggplot2::aes(x = x - cap_w / 2, xend = x + cap_w / 2, y = max, yend = max), linewidth = 0.35, colour = colours[['line']], inherit.aes = FALSE) +
    ggplot2::geom_segment(data = BOX, ggplot2::aes(x = x - cap_w / 2, xend = x + cap_w / 2, y = min, yend = min), linewidth = 0.35, colour = colours[['line']], inherit.aes = FALSE) +
    ggplot2::geom_segment(data = BOX, ggplot2::aes(x = x - box_w / 2, xend = x + box_w / 2, y = med, yend = med), linewidth = 0.8, colour = colours[['line']], inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq_along(x_labs), labels = x_labs) +
    ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
    ggplot2::labs(x = NULL, y = y_lab) +
    theme_fn(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, face = 'italic', size = 11),
      axis.text.y = ggplot2::element_text(size = 11),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 40, 10, 10)
    )
  if (isTRUE(show_points)) {
    prof <- attr(SKEW, 'sra_profile', exact = TRUE)
    if (is.null(prof)) {
      .gama_msg('`show_points` is TRUE but no `sra_profile` found on `SKEW`; plotting boxplots only.')
    } else {
      unit_col_lower <- tolower(unit_col)
      counts_df <- prof |>
        dplyr::filter(.data$species %in% CORE$species) |>
        dplyr::count(species, .data[[unit_col_lower]], name = 'count') |>
        dplyr::mutate(species_label = if (abbreviate) .shorten_species(species) else species) |>
        dplyr::mutate(species_label = factor(species_label, levels = levels(CORE$species_label))) |>
        dplyr::mutate(x = as.numeric(species_label))
      p <- p + ggplot2::geom_jitter(
        data = counts_df, ggplot2::aes(x = x, y = count),
        width = 0.08, height = 0,
        size = 1.6, shape = 19,
        colour = point_colour, alpha = point_alpha,
        inherit.aes = FALSE
      )
    }
  }
  if (isTRUE(show_labels)) {
    lab_df <- CORE
    lab_df$label_y <- lab_df$max * 1.2
    lab_df$lab <- sprintf(paste0('eff = %0.', label_digits, 'f (n = %d)'), lab_df$eff, lab_df[[unit_col]])
    p <- p +
      ggplot2::geom_text(data = lab_df, ggplot2::aes(x = x, y = label_y, label = lab), vjust = 0, size = 3.5, inherit.aes = FALSE) +
      ggplot2::coord_cartesian(clip = 'off')
  }
  p
}

#' Plot BioSample anatomy composition
#'
#' Visualises species-level BioSample anatomy composition using stacked
#' horizontal bars. Each bar shows the proportional contribution of major
#' anatomy classes (`aerial`, `ground`, `reproductive`, `whole`, `in_vitro`,
#' `other`, `mixed`, `unknown`), with operable BioSample counts labelled.
#'
#' Operates on the wide-format summary returned by
#' [summarise_biosample_availability()]. Proportions are calculated over
#' operable BioSamples, where operable means that an accepted sample-source
#' attribute was present.
#'
#' @param BIO A wide-format BioSample summary table returned by
#' [summarise_biosample_availability()].
#' @param species `NULL` (default) to plot all species, or a character vector
#' of species to include.
#' @param rank Ordering of species. One of `highest` (default), `lowest`,
#' `A-Z`, `Z-A`, or `input`.
#' @param abbreviate Logical; if `TRUE` (default), abbreviate species names.
#' @param theme_fn A ggplot2 theme function.
#' @param colours Named character vector of class colours.
#'
#' @return A ggplot object showing proportional BioSample anatomy profiles
#' across species.
#'
#' @seealso [summarise_biosample_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' BIO_SUMMARY <- summarise_biosample_availability(RESULTS)
#' plot_biosample_availability(BIO_SUMMARY)
#' }
#' @export
plot_biosample_availability <- function(
    BIO,
    species    = NULL,
    rank       = 'highest',
    abbreviate = TRUE,
    theme_fn   = ggplot2::theme_minimal,
    colours    = c(
      aerial       = '#E3F9F2',
      ground       = '#A8E6CF',
      reproductive = '#63C7A8',
      whole        = '#329C87',
      in_vitro     = '#1F7A73',
      other        = '#145A52',
      mixed        = '#062A22',
      unknown      = '#BDBDBD'
    )
) {
  rank <- .gama_validate_parameters(rank, 'rank', .gama_rank_parameters, multiple = FALSE, allow_null = FALSE)
  abbreviate <- .gama_validate_logical_parameter(abbreviate, 'abbreviate')
  class_levels <- .biosample_anatomy_profile_levels()
  req_cols <- c('species', 'operable', class_levels)
  BIO <- .gama_require_output(
    BIO,
    'summarise_biosample_availability',
    required_cols = req_cols
  )
  CORE <- BIO |>
    dplyr::select(species, operable, dplyr::all_of(class_levels))
  if (!is.null(species)) {
    missing_sp <- setdiff(species, CORE$species)
    if (length(missing_sp)) .gama_warn(sprintf('Requested species not found in input `BIO`: %s. Dropping.', paste(missing_sp, collapse = ', ')))
    CORE <- CORE |> dplyr::filter(.data$species %in% .env$species)
    if (!nrow(CORE)) .gama_stop('No matching species found.')
  }
  dropped <- unique(CORE$species[is.na(CORE$operable) | CORE$operable == 0L])
  CORE <- CORE |> dplyr::filter(!is.na(.data$operable), .data$operable > 0L)
  if (!nrow(CORE)) .gama_stop('No data to plot: all selected species have zero operable BioSamples.')
  if (length(dropped) > 0L) {
    .gama_msg(sprintf(
      'Dropping %d species with zero operable BioSamples: %s.',
      length(dropped),
      paste(dropped, collapse = ', ')
    ))
  }
  TOTAL <- CORE[, c('species', 'operable')]
  if (rank == 'highest') {
    TOTAL <- TOTAL[order(TOTAL$operable, decreasing = TRUE), ]
  } else if (rank == 'lowest') {
    TOTAL <- TOTAL[order(TOTAL$operable, decreasing = FALSE), ]
  } else if (rank == 'A-Z') {
    TOTAL <- TOTAL[order(TOTAL$species, decreasing = FALSE), ]
  } else if (rank == 'Z-A') {
    TOTAL <- TOTAL[order(TOTAL$species, decreasing = TRUE), ]
  } else if (rank == 'input') {
    TOTAL <- TOTAL
  }
  species_order <- TOTAL$species
  LONG <- CORE |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(class_levels),
      names_to  = 'class',
      values_to = 'count'
    ) |>
    dplyr::mutate(prop = .data$count / .data$operable)
  LONG$class <- factor(LONG$class, levels = class_levels)
  LONG$species_label <- if (abbreviate) .shorten_species(LONG$species) else LONG$species
  levels_in_plot <- if (abbreviate) .shorten_species(species_order) else species_order
  LONG$species_label <- factor(LONG$species_label, levels = rev(levels_in_plot))
  legend_labels <- parse(text = unname(.biosample_anatomy_profile_labels(parse = TRUE)))
  p <- ggplot2::ggplot(LONG, ggplot2::aes(x = species_label, y = prop, fill = class)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = colours, breaks = class_levels, labels = legend_labels) +
    ggplot2::coord_flip(clip = 'off') +
    ggplot2::labs(x = NULL, y = 'Proportion of operable BioSamples', fill = 'Class:') +
    theme_fn(base_size = 13) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 40, 10, 10),
      legend.position = 'top',
      axis.text.y = ggplot2::element_text(size = 11, face = 'italic'),
      axis.text.x = ggplot2::element_text(size = 11),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank()
    )
  totals_df <- TOTAL
  totals_df$species_label <- if (abbreviate) .shorten_species(totals_df$species) else totals_df$species
  totals_df$species_label <- factor(totals_df$species_label, levels = levels(LONG$species_label))
  p +
    ggplot2::geom_text(
      data = totals_df,
      ggplot2::aes(x = species_label, y = 1.005, label = operable),
      hjust = 0,
      size = 3.5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
}

#' Plot SRA-BioSample interaction
#'
#' Visualises SRA-BioSample interaction summaries as modality-by-anatomy
#' heatmaps for one or more species.
#'
#' `plot_interaction()` operates on the summary returned by
#' [summarise_interaction()]. The plotted anatomy resolution is inherited from
#' the `INTERACTION` object: if the summary was generated with
#' `level = 'anatomy_class'`, columns represent anatomy classes; if it was
#' generated with `level = 'anatomy_subclass'`, columns represent anatomy
#' subclasses.
#'
#' @details
#' With `value = 'count'`, the heatmap shows the observed number of linked
#' BioSamples in each SRA modality-by-anatomy cell.
#'
#' With `value = 'residual'`, the heatmap shows Pearson residuals from the
#' species-level marginal expectation calculated by [summarise_interaction()].
#' Positive residuals indicate modality-anatomy combinations that are
#' over-represented relative to the marginal distributions; negative residuals
#' indicate under-represented combinations.
#'
#' When multiple species are plotted, species can be ordered with `rank`.
#' The ranking is applied to the filtered interaction summary before plots are
#' generated.
#'
#' @param INTERACTION A summary table returned by [summarise_interaction()].
#' @param species `NULL` (default) to plot all species, or a character vector
#' specifying which species to include.
#' @param value Heatmap value to display. One of `count` or `residual`.
#' @param rank Ordering of species. One of `highest` (default), `lowest`,
#' `A-Z`, `Z-A`, or `input`.
#' @param theme_fn A ggplot2 theme function.
#' @param low_fill Low count fill colour.
#' @param high_fill High count fill colour.
#' @param positive_low_fill Low positive residual fill colour.
#' @param positive_high_fill High positive residual fill colour.
#' @param negative_low_fill Low negative residual fill colour.
#' @param negative_high_fill High negative residual fill colour.
#' @param zero_fill Zero residual fill colour.
#' @param na_fill Missing value fill colour.
#' @param show_values Logical; if `TRUE` (default), show cell values.
#' @param value_size Numeric text size for cell values.
#'
#' @return A ggplot object for a single species, or a named list of ggplot
#' objects when multiple species are plotted.
#'
#' @seealso [summarise_interaction()], [summarise_sra_availability()],
#' [summarise_biosample_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' BIO_SUMMARY <- summarise_biosample_availability(RESULTS)
#' INTERACTION <- summarise_interaction(
#'   SRA_SUMMARY,
#'   BIO_SUMMARY,
#'   species = 'Vigna angularis'
#' )
#' plot_interaction(
#'   INTERACTION,
#'   value = 'count'
#' )
#' plot_interaction(
#'   INTERACTION,
#'   value = 'residual'
#' )
#' }
#' @export
plot_interaction <- function(INTERACTION,
                             species = NULL,
                             value = 'count',
                             rank = 'highest',
                             theme_fn = ggplot2::theme_minimal,
                             low_fill  = '#E3F9F2',
                             high_fill = '#56C5A8',
                             positive_low_fill  = '#FADDD4',
                             positive_high_fill = '#DC7D67',
                             negative_low_fill  = '#EDE5F6',
                             negative_high_fill = '#8D65B0',
                             zero_fill = '#F2F2F2',
                             na_fill   = '#BDBDBD',
                             show_values = TRUE,
                             value_size = 2.4) {
  value <- .gama_validate_parameters(value, 'value', c('count', 'residual'), multiple = FALSE, allow_null = FALSE)
  rank <- .gama_validate_parameters(rank, 'rank', .gama_rank_parameters, multiple = FALSE, allow_null = FALSE)
  show_values <- .gama_validate_logical_parameter(show_values, 'show_values')
  INTERACTION <- .gama_require_output(
    INTERACTION,
    'summarise_interaction',
    required_cols = c('species', 'BioSample', 'expected', 'residual')
  )
  if (!'class' %in% names(INTERACTION)) {
    if ('modality_class' %in% names(INTERACTION)) {
      INTERACTION <- INTERACTION |>
        dplyr::rename(class = modality_class)
    } else {
      .gama_input_error(
        'summarise_interaction',
        detected = .detect_gama_object(INTERACTION),
        detail = 'missing modality class column: class.'
      )
    }
  }
  info <- attr(INTERACTION, 'interaction_info', exact = TRUE)
  if (is.null(info)) info <- list()
  term_col <- if (!is.null(info$term_col)) info$term_col else if ('anatomy_class' %in% names(INTERACTION)) 'anatomy_class' else 'anatomy_subclass'
  if (!term_col %in% names(INTERACTION)) {
    .gama_input_error(
      'summarise_interaction',
      detected = .detect_gama_object(INTERACTION),
      detail = 'missing anatomy term column.'
    )
  }
  level <- if (!is.null(info$level)) info$level else term_col
  modality_levels <- info$modality_levels
  if (is.null(modality_levels)) modality_levels <- unique(as.character(INTERACTION$class))
  term_meta <- info$term_meta
  if (is.null(term_meta)) {
    term_levels <- unique(as.character(INTERACTION[[term_col]]))
    term_meta <- tibble::tibble(
      term_order = seq_along(term_levels),
      x_term = term_levels,
      anatomy_class = term_levels
    )
  }
  term_levels <- term_meta$x_term
  species_order <- function(tbl, species = NULL, rank = 'highest') {
    species_all <- unique(as.character(tbl$species))
    species_all <- species_all[!is.na(species_all) & nzchar(species_all)]
    if (!length(species_all)) .gama_stop('No species found in `INTERACTION`.')
    if (is.null(species)) {
      sp <- species_all
    } else {
      species_in <- unique(as.character(species))
      species_in <- species_in[!is.na(species_in) & nzchar(species_in)]
      missing_species <- setdiff(species_in, species_all)
      if (length(missing_species)) {
        .gama_warn('Requested species not found in input `INTERACTION`: ', paste(missing_species, collapse = ', '), '. Dropping.')
      }
      sp <- species_in[species_in %in% species_all]
      if (!length(sp)) .gama_stop('No matching species found.')
    }
    if (rank %in% c('highest', 'lowest')) {
      ord <- tbl |>
        dplyr::filter(.data$species %in% .env$sp) |>
        dplyr::group_by(.data$species) |>
        dplyr::summarise(metric_value = sum(.data$BioSample, na.rm = TRUE), .groups = 'drop')
      if (rank == 'highest') ord <- ord[order(ord$metric_value, decreasing = TRUE), , drop = FALSE]
      if (rank == 'lowest') ord <- ord[order(ord$metric_value, decreasing = FALSE), , drop = FALSE]
      sp <- ord$species
    }
    if (rank == 'A-Z') sp <- sort(sp, decreasing = FALSE)
    if (rank == 'Z-A') sp <- sort(sp, decreasing = TRUE)
    sp
  }
  transform_fill <- function(x) {
    x <- as.numeric(x)
    if (value == 'residual') return(sign(x) * log10(abs(x) + 1))
    log10(x + 1)
  }
  sp <- species_order(INTERACTION, species = species, rank = rank)
  class_core_levels <- .biosample_anatomy_class_levels()
  plot_one <- function(one_species) {
    MAT <- INTERACTION |>
      dplyr::filter(.data$species == .env$one_species) |>
      dplyr::mutate(
        x_term = factor(.data[[term_col]], levels = term_levels),
        class = factor(as.character(.data$class), levels = modality_levels),
        label = if (.env$value == 'count') {
          dplyr::if_else(.data$BioSample > 0L, as.character(.data$BioSample), '')
        } else {
          dplyr::if_else(
            !is.na(.data$residual) & abs(.data$residual) > 0,
            sprintf('%.2f', .data$residual),
            ''
          )
        },
        fill_raw = if (.env$value == 'count') {
          as.numeric(.data$BioSample)
        } else {
          as.numeric(.data$residual)
        },
        fill_value = transform_fill(.data$fill_raw)
      )
    if (!nrow(MAT)) .gama_stop('No interaction data found for ', one_species, '.')
    finite_raw <- MAT$fill_raw[is.finite(MAT$fill_raw)]
    if (value == 'count') {
      max_raw <- if (length(finite_raw)) max(finite_raw, na.rm = TRUE) else NA_real_
      legend_raw <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
      if (is.finite(max_raw) && max_raw > 0) legend_raw <- legend_raw[legend_raw <= max_raw]
      if (!length(legend_raw)) legend_raw <- 1L
      legend_breaks <- transform_fill(legend_raw)
      fill_scale <- ggplot2::scale_fill_gradient(
        name = 'Count',
        low = low_fill,
        high = high_fill,
        na.value = na_fill,
        breaks = legend_breaks,
        labels = legend_raw,
        guide = ggplot2::guide_colourbar(
          title.position = 'top',
          title.hjust = 0.5
        )
      )
    } else {
      max_abs <- if (length(finite_raw)) max(abs(finite_raw), na.rm = TRUE) else NA_real_
      has_residual <- is.finite(max_abs) && max_abs > 0
      residual_limit <- if (has_residual) max_abs else 1
      legend_pos <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
      legend_pos <- legend_pos[legend_pos <= residual_limit]
      if (length(legend_pos)) {
        legend_raw <- c(-rev(legend_pos), legend_pos)
      } else if (has_residual) {
        legend_raw <- c(-residual_limit, residual_limit)
      } else {
        legend_raw <- 0
      }
      legend_breaks <- transform_fill(legend_raw)
      residual_limits <- transform_fill(c(-residual_limit, residual_limit))
      if (residual_limit > 1) {
        residual_anchor_raw <- c(-residual_limit, -1, 0, 1, residual_limit)
        residual_colours <- c(
          negative_high_fill, negative_low_fill, zero_fill,
          positive_low_fill, positive_high_fill
        )
      } else {
        residual_anchor_raw <- c(-residual_limit, 0, residual_limit)
        residual_colours <- c(negative_low_fill, zero_fill, positive_low_fill)
      }
      residual_anchor_values <- transform_fill(residual_anchor_raw)
      residual_anchor_values <- (residual_anchor_values - residual_limits[1]) / diff(residual_limits)
      fill_scale <- ggplot2::scale_fill_gradientn(
        name = 'Residual',
        colours = residual_colours,
        values = residual_anchor_values,
        limits = residual_limits,
        na.value = na_fill,
        breaks = legend_breaks,
        labels = legend_raw,
        guide = ggplot2::guide_colourbar(
          title.position = 'top',
          title.hjust = 0.5
        )
      )
    }
    x_labels <- if (level == 'anatomy_class') {
      label_expr <- unname(.biosample_anatomy_profile_labels(parse = TRUE)[term_levels])
      label_expr[term_levels == 'in_vitro'] <- "italic('in vitro')"
      parse(text = label_expr)
    } else {
      label_text <- gsub('_', ' ', term_levels)
      label_expr <- ifelse(
        label_text == 'in vitro',
        "italic('in vitro')",
        paste0("'", label_text, "'")
      )
      parse(text = label_expr)
    }
    boundaries <- numeric()
    if (level == 'anatomy_subclass') {
      boundary_meta <- term_meta |>
        dplyr::arrange(.data$term_order) |>
        dplyr::mutate(anatomy_class = as.character(.data$anatomy_class))
      boundaries <- which(
        utils::head(boundary_meta$anatomy_class, -1L) !=
          utils::tail(boundary_meta$anatomy_class, -1L)
      ) + 0.5
    }
    p <- ggplot2::ggplot(
      MAT,
      ggplot2::aes(
        x = .data$x_term,
        y = .data$class,
        fill = .data$fill_value
      )
    ) +
      ggplot2::geom_tile(colour = 'white', linewidth = 0.35) +
      ggplot2::scale_x_discrete(
        limits = term_levels,
        labels = x_labels,
        drop = FALSE
      ) +
      ggplot2::scale_y_discrete(
        limits = rev(modality_levels),
        drop = FALSE
      ) +
      fill_scale +
      ggplot2::labs(
        title = one_species,
        x = 'Sample tissue',
        y = 'Class'
      ) +
      ggplot2::coord_fixed() +
      theme_fn(base_size = 13) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(10, 10, 10, 10),
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(face = 'italic'),
        legend.position = 'right',
        legend.title = ggplot2::element_text(size = 11),
        legend.text = ggplot2::element_text(size = 10)
      )
    if (length(boundaries)) {
      p <- p + ggplot2::geom_vline(
        xintercept = boundaries,
        colour = '#BDBDBD',
        linewidth = 0.35
      )
    }
    if (isTRUE(show_values)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        size = value_size
      )
    }
    attr(p, 'gama_interaction_data') <- MAT
    p
  }
  if (length(sp) == 1L) return(plot_one(sp))
  OUT <- lapply(sp, plot_one)
  names(OUT) <- sp
  OUT
}

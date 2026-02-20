# PLOTTING ====================================================================

#' Plot data richness
#'
#' Visualises species-level data richness using the composite scores generated
#' during the query phase. Assembly, SRA, and BioSample contributions are shown
#' as stacked bar segments for each species.
#'
#' @param SUMMARY A tibble returned by [summarise_availability()], containing
#'   species-level composite scores and component values (`A`, `S`, `B`).
#' @param rank Ordering of species on the x-axis. One of 'highest' (default),
#'   'lowest', 'A-Z', 'Z-A', or 'input'.
#' @param abbreviate Logical; if `TRUE` (default), abbreviate species names to
#'   'G. species' style labels.
#' @param theme_fn A ggplot2 theme function (e.g. [ggplot2::theme_minimal()]).
#' @param colours Named character vector of fill colours for the 'Assembly',
#'   'SRA', and 'BioSample' segments.
#'
#' @return A ggplot object showing stacked bar segments for each species.
#'
#' @seealso [query_species()], [summarise_availability()]
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
rank        = c('highest', 'lowest', 'A-Z', 'Z-A', 'input'),
abbreviate  = TRUE,
theme_fn    = ggplot2::theme_minimal,
colours     = c(
Assembly  = '#E3F9F2',
SRA       = '#A8E6CF',
BioSample = '#56C5A8'
)
) {
  rank <- match.arg(rank)
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

#' Plot SRA modality composition
#'
#' Visualises species-level SRA modality composition using stacked horizontal
#' bars. Each bar shows the proportional contribution of major SRA classes
#' (genomic, transcriptomic, epigenomic, chromatin, other, unknown), with total
#' SRA counts labelled.
#'
#' Operates on the wide-format summary returned by
#' [summarise_sra_availability()].
#'
#' @param SRA A wide-format SRA summary table returned by
#'   [summarise_sra_availability()].
#' @param species `NULL` (default) to plot all species, or a character vector
#'   of species to include.
#' @param rank Ordering of species. One of 'highest' (default), 'lowest',
#'   'A-Z', 'Z-A', or 'input'.
#' @param abbreviate Logical; if `TRUE` (default), abbreviate species names.
#' @param theme_fn A ggplot2 theme function.
#' @param colours Named character vector of class colours.
#'
#' @return A ggplot object showing proportional SRA modality profiles across
#'   species.
#'
#' @seealso [summarise_sra_availability()], [plot_sra_geo_availability()]
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
    rank       = c('highest', 'lowest', 'A-Z', 'Z-A', 'input'),
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
  rank <- match.arg(rank)
  CORE <- SRA |>
    dplyr::select(
      species, SRA, genomic, transcriptomic, epigenomic, chromatin, other,
      unknown
    )
  if (!is.null(species)) {
    CORE <- CORE |> dplyr::filter(.data$species %in% species)
    if (!nrow(CORE)) stop('No matching species found')
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
#' Visualises GEO linkage as an overlay on SRA modality classes. Each modality
#' is displayed as a 100% bar in a neutral background, with the modality colour
#' representing the GEO-linked fraction. Labels show `GEO-linked / Total` for
#' each modality.
#'
#' Operates on the wide-format output of [summarise_sra_availability()] when
#' `include_geo = TRUE`.
#'
#' If `species` is `NULL`, plots are generated for all species in the table.
#' A single species returns a ggplot object; multiple species returns a named
#' list of ggplot objects.
#'
#' @param SRA A wide-format SRA summary table returned by
#'   [summarise_sra_availability()] when `include_geo = TRUE`.
#' @param species `NULL` (default) for all species, or a character vector of
#'   species to plot.
#' @param classes Character vector of modality classes to display. Values are
#'   intersected with the fixed plotting order used internally.
#' @param rank Ordering of species (when `species = NULL`), or ordering applied
#'   to the requested species vector. One of 'highest', 'lowest', 'A-Z', 'Z-A',
#'   or 'input'.
#' @param theme_fn A ggplot2 theme function.
#' @param colours Named character vector of class colours.
#' @param alpha_vals Named numeric vector giving alpha values for GEO-linked
#'   vs not GEO-linked segments.
#'
#' @return A ggplot object (single species) or a named list of ggplot objects
#'   (multiple species).
#'
#' @seealso [summarise_sra_availability()], [plot_sra_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SRA_GEO <- summarise_sra_availability(RESULTS, include_geo = TRUE)
#' plot_sra_geo_availability(SRA_GEO, species = 'Vigna vexillata')
#' }
#' @export
plot_sra_geo_availability <- function(
SRA,
species    = NULL,
classes    = c('transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown'),
rank       = c('highest', 'lowest', 'A-Z', 'Z-A', 'input'),
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
  rank <- match.arg(rank)
  FIXED_ORDER <- c('unknown', 'other', 'chromatin', 'epigenomic', 'transcriptomic')
  classes <- intersect(FIXED_ORDER, classes)
  geo_cols <- paste0(classes, '_geo')
  if (!all(geo_cols %in% names(SRA))) {
    stop('GEO summary columns missing. Run summarise_sra_availability(..., include_geo = TRUE).')
  }
  if (is.null(species)) {
    sp <- SRA$species
    if (rank != 'input') {
      ord <- SRA[, c('species', 'SRA')]
      if (rank == 'highest') {
        ord <- ord[order(ord$SRA, decreasing = TRUE), ]
      } else if (rank == 'lowest') {
        ord <- ord[order(ord$SRA, decreasing = FALSE), ]
      } else if (rank == 'A-Z') {
        ord <- ord[order(ord$species, decreasing = FALSE), ]
      } else if (rank == 'Z-A') {
        ord <- ord[order(ord$species, decreasing = TRUE), ]
      }
      sp <- ord$species
    }
  } else {
    sp <- intersect(species, SRA$species)
    if (!length(sp)) stop('No matching species found')
    if (rank == 'A-Z') {
      sp <- sort(sp, decreasing = FALSE)
    } else if (rank == 'Z-A') {
      sp <- sort(sp, decreasing = TRUE)
    }
  }
  .plot_one <- function(one_species) {
    row <- SRA |> dplyr::filter(.data$species == !!one_species)
    if (!nrow(row)) stop('No matching species found')
    total  <- as.integer(row[1, classes,  drop = TRUE])
    linked <- as.integer(row[1, geo_cols, drop = TRUE])
    total[is.na(total)]   <- 0L
    linked[is.na(linked)] <- 0L
    prop_linked <- ifelse(total > 0, linked / total, 0)
    prop_not    <- ifelse(total > 0, 1 - prop_linked, 0)
    SUM <- tibble::tibble(
    class       = factor(classes, levels = FIXED_ORDER),
    total       = total,
    linked      = linked,
    label       = paste0(linked, '/', total),
    prop_linked = prop_linked,
    prop_not    = prop_not
    )
    LONG <- dplyr::bind_rows(
    SUM |> dplyr::transmute(class, status = 'Not GEO-linked', prop = prop_not),
    SUM |> dplyr::transmute(class, status = 'GEO-linked',     prop = prop_linked)
    )
    LONG$class  <- factor(LONG$class, levels = FIXED_ORDER)
    LONG$status <- factor(LONG$status, levels = c('Not GEO-linked', 'GEO-linked'))
    ggplot2::ggplot(
    LONG,
    ggplot2::aes(x = class, y = prop, fill = class, alpha = status)
    ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
    data = SUM,
    ggplot2::aes(x = class, y = 1.02, label = label),
    size = 3.5,
    vjust = 0,
    inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = colours, drop = FALSE) +
    ggplot2::scale_alpha_manual(values = alpha_vals) +
    ggplot2::guides(
    fill  = ggplot2::guide_legend(reverse = TRUE),
    alpha = ggplot2::guide_legend(reverse = FALSE)
    ) +
    ggplot2::coord_flip(clip = 'off') +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.12)) +
    ggplot2::labs(
    title = one_species,
    x     = NULL,
    y     = 'Proportion GEO-linked',
    fill  = 'Class:',
    alpha = NULL
    ) +
    theme_fn(base_size = 13) +
    ggplot2::theme(
    plot.margin     = ggplot2::margin(10, 40, 10, 10),
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
#' Visualises replication skew using the summary output of summarise_sra_skew().
#' Boxplots are drawn from pre-computed five-number summaries on a log10 y-axis
#' and optionally labelled as 'eff=<x> (n=<y>)'. Optionally overlays per-unit
#' points (one point per BioProject/BioSample) jittered horizontally when a
#' cached UID-level profile is available as attr(SKEW, 'sra_profile').
#'
#' @param SKEW A tibble returned by summarise_sra_skew().
#' @param species NULL (default) to include all species, or a character vector.
#' @param rank One of 'highest' (default), 'lowest', 'A-Z', 'Z-A', or 'input'.
#' @param abbreviate Logical; if TRUE (default), abbreviate species labels.
#' @param show_points Logical; if TRUE (default), overlay per-unit points.
#' @param point_colour Colour for the overlaid data points (character).
#' @param theme_fn A ggplot2 theme function.
#' @param colours Named character vector for box fill and line colours.
#' @param show_labels Logical; if TRUE (default), label each box with eff and n.
#' @param label_digits Integer; decimal places for eff labels.
#' @param point_alpha Numeric alpha (transparency) for overlaid points.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' SRA_SUMMARY <- summarise_sra_availability(RESULTS)
#' SKEW  <- summarise_sra_skew(SRA_SUMMARY)
#' plot_sra_skew(SKEW)
#' }
#'
#' @export
plot_sra_skew <- function(
    SKEW,
    species      = NULL,
    rank         = c('highest', 'lowest', 'A-Z', 'Z-A', 'input'),
    abbreviate   = TRUE,
    show_points  = TRUE,
    point_colour = '#56C5A8',
    theme_fn     = ggplot2::theme_minimal,
    colours      = c(box = '#E3F9F2', line = 'black'),
    show_labels  = TRUE,
    label_digits = 1L,
    point_alpha  = 0.25
) {
  rank <- match.arg(rank)
  if (!all(c('species', 'class', 'min', 'q25', 'med', 'q75', 'max', 'eff') %in% names(SKEW))) stop('Input must be the output of summarise_sra_skew().', call. = FALSE)
  unit_col <- if ('BioProject' %in% names(SKEW)) 'BioProject' else if ('BioSample' %in% names(SKEW)) 'BioSample' else stop('Expected a BioProject or BioSample column in `SKEW`.', call. = FALSE)
  CORE0 <- SKEW |> dplyr::select(species, dplyr::all_of(unit_col), class, min, q25, med, q75, max, eff)
  if (!is.null(species)) {
    CORE0 <- CORE0 |> dplyr::filter(.data$species %in% species)
    if (!nrow(CORE0)) stop('No matching species found in `SKEW` for the requested filter.', call. = FALSE)
  }
  class_vals <- unique(CORE0$class)
  if (length(class_vals) != 1L) stop('plot_sra_skew() expects a single class per plot.', call. = FALSE)
  class_val <- class_vals[[1]]
  dropped <- unique(CORE0$species[is.na(CORE0$max)])
  CORE <- CORE0 |> dplyr::filter(!is.na(.data$max))
  if (!nrow(CORE)) {
    stop(sprintf('No data to plot for class \'%s\': all selected species have zero records for this class.', class_val), call. = FALSE)
  }
  if (length(dropped) > 0L) {
    message(sprintf(
      'Dropping %d species with no \'%s\' data: %s',
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
      message('show_points = TRUE but no sra_profile found on SKEW; plotting boxplots only.')
    } else {
      unit_col_lower <- tolower(unit_col)
      counts_df <- prof |>
        dplyr::filter(.data$species %in% CORE$species) |>
        dplyr::count(species, !!rlang::sym(unit_col_lower), name = 'count') |>
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

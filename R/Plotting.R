# PLOTTING ====================================================================
#' Plot data richness
#'
#' @param SUMMARY Tibble returned by [summarise_availability()].
#' @param rank Character; species ordering method.
#' @param abbreviate Logical; abbreviate species names.
#' @param theme_fn ggplot2 theme function.
#' @param colours Named character vector of fill colours.
#' @return A ggplot object showing stacked bar charts of data richness components.
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
  shorten_species <- function(x) {
    vapply(
      strsplit(x, ' '),
      function(parts) {
        if (length(parts) < 2) return(x)
        paste0(substr(parts[1], 1, 1), '. ', parts[2])
      },
      character(1)
    )
  }
  SUMMARY$species_label <- if (abbreviate) shorten_species(SUMMARY$species) else SUMMARY$species
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
#' Plot proportional SRA experimental modality composition
#'
#' @param SRA Tibble returned by [summarise_sra_availability()].
#' @param species Optional character vector of species to include.
#' @param rank Character; species ordering method.
#' @param abbreviate Logical; abbreviate species names.
#' @param theme_fn ggplot2 theme function.
#' @param colours Named character vector of class colours.
#' @return A ggplot object showing proportional SRA modality composition across species.
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
    dplyr::select(species, SRA, genomic, transcriptomic, epigenomic, chromatin, other, unknown)
  if (!is.null(species)) {
    CORE <- CORE |> dplyr::filter(.data$species %in% species)
    if (!nrow(CORE)) stop('No matching species found')
  }
  TOTAL <- CORE[, c('species', 'SRA')]
  if (rank == 'highest') {
    TOTAL <- TOTAL[order(TOTAL$SRA, decreasing = FALSE), ]
  } else if (rank == 'lowest') {
    TOTAL <- TOTAL[order(TOTAL$SRA, decreasing = TRUE), ]
  } else if (rank == 'A-Z') {
    TOTAL <- TOTAL[order(TOTAL$species, decreasing = TRUE), ]
  } else if (rank == 'Z-A') {
    TOTAL <- TOTAL[order(TOTAL$species, decreasing = FALSE), ]
  }
  species_order <- TOTAL$species
  LONG <- CORE |>
    tidyr::pivot_longer(
      cols      = c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown'),
      names_to  = 'class',
      values_to = 'count'
    ) |>
    dplyr::group_by(species) |>
    dplyr::mutate(prop = count / sum(count)) |>
    dplyr::ungroup()
  LONG$class <- factor(
    LONG$class,
    levels = c('genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown')
  )
  shorten_species <- function(x) {
    vapply(
      strsplit(x, ' '),
      function(parts) {
        if (length(parts) < 2) return(x)
        paste0(substr(parts[1], 1, 1), '. ', parts[2])
      },
      character(1)
    )
  }
  LONG$species_label <- if (abbreviate) shorten_species(LONG$species) else LONG$species
  LONG$species_label <- factor(
    LONG$species_label,
    levels = if (abbreviate) shorten_species(species_order) else species_order
  )
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
  totals_df$species_label <- factor(
    if (abbreviate) shorten_species(totals_df$species) else totals_df$species,
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
#' Plot GEO linkage overlay for SRA experimental modalities
#'
#' @param x Tibble returned by [summarise_sra_availability()] with include_geo = TRUE
#' @param species Optional character vector of species to plot.
#' @param classes Character vector of modality classes to display.
#' @param rank Character; species ordering method.
#' @param theme_fn ggplot2 theme function.
#' @param colours Named character vector of class colours.
#' @param alpha_vals Named numeric vector controlling GEO-linked transparency.
#' @return A ggplot object (single species) or named list of ggplot objects showing per-modality GEO linkage.
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

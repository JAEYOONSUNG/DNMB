.dnmb_plot_padloc_module <- function(genbank_table, output_dir) {
  add_panel_header <- function(plot_obj, label, title, x = 0.0005, y = 0.968, size = 12, plot_y = 0.045, plot_h = 0.88) {
    cowplot::ggdraw() +
      cowplot::draw_plot(plot_obj, x = 0, y = plot_y, width = 1, height = plot_h) +
      cowplot::draw_label(paste0(label, "  ", title), x = x, y = y, hjust = 0, vjust = 1, fontface = "bold", size = size)
  }

  tbl <- .dnmb_contig_ordered_table(genbank_table)
  req <- c("PADLOC_system", "PADLOC_system_number")
  if (!nrow(tbl) || !all(req %in% names(tbl))) {
    return(NULL)
  }

  hits <- tbl[!is.na(tbl$PADLOC_system) & nzchar(tbl$PADLOC_system), , drop = FALSE]
  if (!nrow(hits)) {
    return(NULL)
  }

  hits$PADLOC_system_number <- suppressWarnings(as.integer(hits$PADLOC_system_number))
  hits$PADLOC_target_coverage <- suppressWarnings(as.numeric(hits$PADLOC_target_coverage))
  hits$PADLOC_hmm_coverage <- suppressWarnings(as.numeric(hits$PADLOC_hmm_coverage))

  system_summary <- hits |>
    dplyr::group_by(.data$PADLOC_system) |>
    dplyr::summarise(
      n_proteins = dplyr::n(),
      n_loci = dplyr::n_distinct(.data$PADLOC_system_number),
      best_target_cov = max(.data$PADLOC_target_coverage, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_loci), dplyr::desc(.data$n_proteins), .data$PADLOC_system)
  system_summary$PADLOC_system <- factor(system_summary$PADLOC_system, levels = rev(system_summary$PADLOC_system))
  palette <- stats::setNames(
    grDevices::colorRampPalette(c(
      "#F2A93B",
      "#F08F20",
      "#EB7C12",
      "#E16A0D",
      "#D25709",
      "#BE4306"
    ))(length(unique(system_summary$PADLOC_system))),
    unique(as.character(system_summary$PADLOC_system))
  )

  p_inventory <- ggplot2::ggplot(system_summary, ggplot2::aes(x = .data$n_loci, y = .data$PADLOC_system)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$PADLOC_system),
      width = 0.62,
      color = "grey35",
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0.05, label = .data$PADLOC_system),
      hjust = 0,
      size = 3.0,
      color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$n_loci + 0.08,
        label = paste0(.data$n_proteins, " proteins | best cov ", sprintf("%.2f", .data$best_target_cov))
      ),
      hjust = 0,
      size = 2.8,
      color = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.34))) +
    ggplot2::labs(title = NULL, x = "Loci detected", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  loci <- hits |>
    dplyr::group_by(.data$contig, .data$PADLOC_system_number, .data$PADLOC_system) |>
    dplyr::summarise(
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      genes = dplyr::n(),
      best_cov = max(.data$PADLOC_target_coverage, na.rm = TRUE),
      .groups = "drop"
    )
  contigs <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  # Only show contigs that have PADLOC loci
  hit_contigs <- unique(loci$contig)
  contigs <- contigs[contigs$contig %in% hit_contigs, , drop = FALSE]
  loci$label <- ifelse(
    duplicated(loci$PADLOC_system) | duplicated(loci$PADLOC_system, fromLast = TRUE),
    paste0(loci$PADLOC_system, " (", loci$PADLOC_system_number, ")"),
    as.character(loci$PADLOC_system)
  )
  # Facet labels: contig name with size
  contig_labels <- stats::setNames(
    paste0(contigs$sector_label, " (", scales::label_comma()(contigs$length_bp), " bp)"),
    contigs$contig
  )
  loci$contig_facet <- factor(contig_labels[loci$contig], levels = unname(contig_labels))
  contigs$contig_facet <- factor(contig_labels[contigs$contig], levels = unname(contig_labels))

  loci$midpoint <- (loci$start + loci$end) / 2

  p_layout <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contigs,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = 0.5, yend = 0.5),
      linewidth = 0.6, color = "grey78"
    ) +
    ggplot2::geom_rect(
      data = loci,
      ggplot2::aes(
        xmin = .data$start, xmax = .data$end,
        ymin = 0.48, ymax = 0.52,
        fill = .data$PADLOC_system
      ),
      color = "grey35", linewidth = 0.2, alpha = 0.95
    ) +
    ggrepel::geom_text_repel(
      data = loci,
      ggplot2::aes(x = .data$midpoint, y = 0.52, label = .data$label),
      size = 1.9, lineheight = 0.85,
      direction = "y", nudge_y = 0.03,
      segment.size = 0.2, segment.color = "grey55",
      max.overlaps = Inf, seed = 42,
      min.segment.length = 0,
      force = 1.5, force_pull = 0.3,
      box.padding = 0.12,
      ylim = c(0.54, 0.72)
    ) +
    ggplot2::facet_wrap(~contig_facet, ncol = 1, scales = "free_x",
                        strip.position = "top") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(labels = scales::label_comma()) +
    ggplot2::scale_y_continuous(limits = c(0.44, 0.74), expand = c(0, 0)) +
    ggplot2::labs(title = "PADLOC genome layout", x = "Genome coordinate (bp)", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 8.7, margin = ggplot2::margin(t = 1, b = 1)),
      strip.background = ggplot2::element_rect(fill = "grey97", colour = "grey85", linewidth = 0.3),
      panel.spacing.y = grid::unit(0.10, "lines")
    )

  top_hits <- hits |>
    dplyr::mutate(
      padloc_label = paste0(.data$locus_tag, " | ", .data$product),
      padloc_rank = dplyr::dense_rank(dplyr::desc(.data$PADLOC_target_coverage))
    ) |>
    dplyr::arrange(dplyr::desc(.data$PADLOC_target_coverage), dplyr::desc(.data$PADLOC_hmm_coverage), .data$locus_tag) |>
    dplyr::slice_head(n = 15)
  top_hits$padloc_label <- factor(top_hits$padloc_label, levels = rev(top_hits$padloc_label))

  p_hits <- ggplot2::ggplot(top_hits, ggplot2::aes(x = .data$PADLOC_target_coverage, y = .data$padloc_label)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$PADLOC_system),
      width = 0.62,
      color = "grey35",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$PADLOC_target_coverage + 0.02,
        label = paste0(.data$PADLOC_system, " | hmm ", sprintf("%.2f", .data$PADLOC_hmm_coverage))
      ),
      hjust = 0,
      size = 2.7,
      color = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(limits = c(0, max(1.05, max(top_hits$PADLOC_target_coverage, na.rm = TRUE) + 0.22))) +
    ggplot2::labs(title = NULL, x = "Target coverage", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "none"
    )

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "PADLOC_overview.pdf")
  n_contigs <- nrow(contigs)
  legend_plot <- p_layout +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      legend.text = ggplot2::element_text(size = 9.2),
      legend.key.width = grid::unit(1.3, "lines"),
      legend.key.height = grid::unit(0.9, "lines"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))
  common_legend <- cowplot::get_legend(legend_plot)
  legend_panel <- cowplot::ggdraw() + cowplot::draw_grob(common_legend, x = 0.5, y = 0.5, width = 0.98, height = 0.95, hjust = 0.5, vjust = 0.5)

  layout_height <- max(0.65, 0.40 * n_contigs + 0.25)
  composite <- cowplot::plot_grid(
    add_panel_header(p_inventory, "A", "PADLOC inventory", plot_y = 0.06, plot_h = 0.86),
    add_panel_header(p_layout, "B", "PADLOC genome layout", plot_y = 0.035, plot_h = 0.90),
    add_panel_header(p_hits, "C", "Top PADLOC protein calls", plot_y = 0.06, plot_h = 0.86),
    legend_panel,
    ncol = 1,
    rel_heights = c(0.86, layout_height, 1.10, 0.30)
  )
  total_h <- max(10.2, 1.9 + layout_height * 2.15 + 3.1)
  .dnmb_module_plot_save(composite, pdf_path, width = 10, height = total_h)
  list(pdf = pdf_path)
}

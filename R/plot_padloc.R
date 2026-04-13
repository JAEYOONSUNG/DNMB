.dnmb_plot_padloc_module <- function(genbank_table, output_dir) {
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
    grDevices::hcl.colors(length(unique(system_summary$PADLOC_system)), palette = "Dark 3"),
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
      color = "white",
      fontface = "bold"
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
    ggplot2::labs(title = "PADLOC inventory", x = "Loci detected", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(4, 4, 4, 18)
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
  contigs$track <- rev(seq_len(nrow(contigs)))
  loci$track <- contigs$track[match(loci$contig, contigs$contig)]
  loci$label <- ifelse(
    duplicated(loci$PADLOC_system) | duplicated(loci$PADLOC_system, fromLast = TRUE),
    paste0(loci$PADLOC_system, " (", loci$PADLOC_system_number, ")"),
    as.character(loci$PADLOC_system)
  )

  p_layout <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contigs,
      ggplot2::aes(x = 1, xend = .data$length_bp, y = .data$track, yend = .data$track),
      linewidth = 1.2,
      color = "grey78"
    ) +
    ggplot2::geom_rect(
      data = loci,
      ggplot2::aes(
        xmin = .data$start,
        xmax = .data$end,
        ymin = .data$track - 0.24,
        ymax = .data$track + 0.24,
        fill = .data$PADLOC_system
      ),
      color = "grey35",
      linewidth = 0.25,
      alpha = 0.95
    ) +
    ggplot2::geom_text(
      data = loci,
      ggplot2::aes(
        x = (.data$start + .data$end) / 2,
        y = .data$track + 0.36,
        label = .data$label
      ),
      size = 2.7,
      vjust = 0
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(
      breaks = contigs$track,
      labels = NULL,
      expand = ggplot2::expansion(mult = c(0.10, 0.16))
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::labs(title = "PADLOC genome layout", x = "Genome position (bp)", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(
      data = contigs,
      ggplot2::aes(x = -Inf, y = .data$track, label = .data$sector_label),
      hjust = 1.1, size = 2.8, color = "grey30", fontface = "bold"
    ) +
    ggplot2::coord_cartesian(clip = "off")

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
    ggplot2::labs(title = "Top PADLOC protein calls", x = "Target coverage", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "PADLOC_overview.pdf")
  composite <- cowplot::plot_grid(
    p_inventory,
    p_layout,
    p_hits,
    labels = c("A", "B", "C"),
    label_size = 14,
    label_fontface = "bold",
    label_x = 0,
    label_y = c(1.02, 1.02, 1.02),
    hjust = 0,
    ncol = 1,
    rel_heights = c(0.80, 1.35, 1.25)
  )
  .dnmb_module_plot_save(composite, pdf_path, width = 10, height = 11)
  list(pdf = pdf_path)
}

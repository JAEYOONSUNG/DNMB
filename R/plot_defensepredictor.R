.dnmb_plot_defensepredictor_module <- function(genbank_table,
                                               output_dir,
                                               threshold = .dnmb_defensepredictor_default_threshold()) {
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  req <- c("DefensePredictor_mean_log_odds", "DefensePredictor_score_band")
  if (!nrow(tbl) || !all(req %in% names(tbl))) {
    return(NULL)
  }

  tbl <- tbl[!is.na(tbl$DefensePredictor_mean_log_odds), , drop = FALSE]
  if (!nrow(tbl)) {
    return(NULL)
  }

  tbl$DefensePredictor_mean_log_odds <- suppressWarnings(as.numeric(tbl$DefensePredictor_mean_log_odds))
  tbl$DefensePredictor_sd_log_odds <- suppressWarnings(as.numeric(tbl$DefensePredictor_sd_log_odds))
  tbl$is_hit <- !is.na(tbl$DefensePredictor_mean_log_odds) & tbl$DefensePredictor_mean_log_odds >= as.numeric(threshold)[1]
  tbl$score_band_plot <- ifelse(tbl$is_hit, tbl$DefensePredictor_score_band, "DP_lt4")
  band_levels <- c("DP_10plus", "DP_8to10", "DP_6to8", "DP_4to6", "DP_lt4")
  band_palette <- c(
    DP_10plus = "#A50026",
    DP_8to10 = "#D73027",
    DP_6to8 = "#F46D43",
    DP_4to6 = "#FDAE61",
    DP_lt4 = "#BDBDBD"
  )
  tbl$score_band_plot <- factor(tbl$score_band_plot, levels = band_levels)

  p_hist <- ggplot2::ggplot(tbl, ggplot2::aes(x = .data$DefensePredictor_mean_log_odds)) +
    ggplot2::geom_histogram(
      bins = 45,
      fill = "#D81B60",
      color = "white",
      linewidth = 0.2,
      alpha = 0.9
    ) +
    ggplot2::geom_vline(
      xintercept = as.numeric(threshold)[1],
      linetype = "dashed",
      linewidth = 0.7,
      color = "#333333"
    ) +
    ggplot2::annotate(
      "text",
      x = as.numeric(threshold)[1],
      y = Inf,
      label = paste0("threshold = ", as.numeric(threshold)[1]),
      vjust = 1.6,
      hjust = -0.02,
      size = 3.1,
      color = "#333333"
    ) +
    ggplot2::labs(x = "Mean log-odds", y = "Proteins") +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "DefensePredictor score distribution",
                      hjust = -0.02, vjust = 1.5, size = 4, fontface = "bold") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 6, 5, 6)
    )

  top_hits <- tbl |>
    dplyr::arrange(dplyr::desc(.data$DefensePredictor_mean_log_odds), .data$locus_tag) |>
    dplyr::slice_head(n = 18) |>
    dplyr::mutate(
      dp_label = paste0(.data$locus_tag, " | ", .data$product),
      dp_label = factor(.data$dp_label, levels = rev(.data$dp_label))
    )
  p_top <- ggplot2::ggplot(top_hits, ggplot2::aes(x = .data$DefensePredictor_mean_log_odds, y = .data$dp_label)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$score_band_plot),
      width = 0.62,
      color = "grey30",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$DefensePredictor_mean_log_odds + 0.18,
        label = paste0("sd ", sprintf("%.2f", .data$DefensePredictor_sd_log_odds))
      ),
      hjust = 0,
      size = 2.7,
      color = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = band_palette, drop = FALSE) +
    ggplot2::labs(x = "Mean log-odds", y = NULL) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "Top DefensePredictor candidates",
                      hjust = -0.02, vjust = 1.5, size = 4, fontface = "bold") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 6, 5, 6)
    )

  contigs <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  contigs$track <- rev(seq_len(nrow(contigs)))
  tbl$track <- contigs$track[match(tbl$contig, contigs$contig)]
  tbl$midpoint <- (tbl$start + tbl$end) / 2

  p_layout <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contigs,
      ggplot2::aes(x = 1, xend = .data$length_bp, y = .data$track, yend = .data$track),
      linewidth = 1.1,
      color = "grey80"
    ) +
    ggplot2::geom_point(
      data = tbl,
      ggplot2::aes(
        x = .data$midpoint,
        y = .data$track,
        fill = .data$score_band_plot,
        size = pmax(.data$DefensePredictor_mean_log_odds, 0)
      ),
      shape = 21,
      color = "grey20",
      stroke = 0.22,
      alpha = 0.88
    ) +
    ggplot2::scale_fill_manual(values = band_palette, drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(1.1, 6.2), guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = contigs$track,
      labels = NULL,
      expand = ggplot2::expansion(mult = c(0.10, 0.16))
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0.005, 0.005))
    ) +
    ggplot2::labs(x = "Genome position (bp)", y = NULL, fill = "Score band") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 6, 5, 6)
    ) +
    ggplot2::geom_text(
      data = contigs,
      ggplot2::aes(x = 1, y = .data$track + 0.35, label = .data$sector_label),
      hjust = 0, size = 2.6, color = "grey40", fontface = "bold"
    ) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "DefensePredictor genome layout",
                      hjust = -0.02, vjust = 1.5, size = 4, fontface = "bold")

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "DefensePredictor_overview.pdf")
  # Layout order matches DefenseFinder: A=inventory, B=genome layout, C=candidates
  composite <- cowplot::plot_grid(
    p_hist,
    p_layout,
    p_top,
    labels = c("A", "B", "C"),
    label_size = 14,
    label_fontface = "bold",
    label_x = 0,
    label_y = c(1.02, 1.02, 1.02),
    hjust = 0,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(0.80, 1.35, 1.25)
  )
  .dnmb_module_plot_save(composite, pdf_path, width = 10, height = 11)
  list(pdf = pdf_path)
}

.dnmb_antidefense_module_prefix <- function(module_name) {
  module_name <- trimws(as.character(module_name)[1])
  switch(
    module_name,
    dbAPIS = "dbAPIS",
    AcrFinder = "AcrFinder",
    stop("Unsupported anti-defense module: ", module_name, call. = FALSE)
  )
}

.dnmb_antidefense_module_file_stub <- function(module_name) {
  module_name <- trimws(as.character(module_name)[1])
  switch(
    module_name,
    dbAPIS = "dbAPIS",
    AcrFinder = "AcrFinder",
    stop("Unsupported anti-defense module: ", module_name, call. = FALSE)
  )
}

.dnmb_antidefense_module_dir <- function(output_dir, module_name) {
  stub <- switch(
    trimws(as.character(module_name)[1]),
    dbAPIS = "dbapis",
    AcrFinder = "acrfinder",
    stop("Unsupported anti-defense module: ", module_name, call. = FALSE)
  )
  file.path(output_dir, paste0("dnmb_module_", stub))
}

.dnmb_antidefense_panel_header <- function(plot_obj,
                                           label,
                                           title,
                                           x = 0.0005,
                                           y = 0.968,
                                           size = 12,
                                           plot_y = 0.045,
                                           plot_h = 0.88) {
  cowplot::ggdraw() +
    cowplot::draw_plot(plot_obj, x = 0, y = plot_y, width = 1, height = plot_h) +
    cowplot::draw_label(
      paste0(label, "  ", title),
      x = x,
      y = y,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = size
    )
}

.dnmb_antidefense_display_palette <- function(values,
                                              colors = c(
                                                "#F2A93B",
                                                "#F08F20",
                                                "#EB7C12",
                                                "#E16A0D",
                                                "#D25709",
                                                "#BE4306"
                                              )) {
  values <- unique(stats::na.omit(as.character(values)))
  if (!length(values)) {
    return(stats::setNames(character(), character()))
  }
  stats::setNames(grDevices::colorRampPalette(colors)(length(values)), values)
}

.dnmb_antidefense_target_label <- function(x, max_len = 24L) {
  x <- .dnmb_integrated_defense_clean_text(x)
  x <- gsub("NAD\\+ reconstitution pathway \\(NARP\\)", "NARP", x, ignore.case = TRUE)
  x <- gsub("pyrimidine cyclase system for antiphage resistance \\(Pycsar\\)", "Pycsar", x, ignore.case = TRUE)
  x <- gsub("CRISPR[-–]Cas evasion by DNA repair", "CRISPR repair evasion", x, ignore.case = TRUE)
  .dnmb_integrated_defense_short_label(x, max_len = max_len)
}

.dnmb_antidefense_module_has_outputs <- function(output_dir, module_name) {
  module_dir <- .dnmb_antidefense_module_dir(output_dir, module_name)
  if (!dir.exists(module_dir)) {
    return(FALSE)
  }
  any(file.exists(list.files(module_dir, recursive = TRUE, full.names = TRUE)))
}

.dnmb_antidefense_module_hit_table <- function(genbank_table, module_name) {
  prefix <- .dnmb_antidefense_module_prefix(module_name)
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  family_col <- paste0(prefix, "_family_id")
  if (!family_col %in% names(tbl)) {
    return(data.frame())
  }

  hit_label_col <- paste0(prefix, "_hit_label")
  support_col <- paste0(prefix, "_support")
  evidence_col <- paste0(prefix, "_evidence_mode")
  role_col <- paste0(prefix, "_enzyme_role")
  subtype_col <- paste0(prefix, "_defense_type")
  clan_col <- paste0(prefix, "_clan_defense_type")
  pident_col <- paste0(prefix, "_pident")
  evalue_col <- paste0(prefix, "_i_evalue")
  score_col <- paste0(prefix, "_hit_score")
  hmm_cov_col <- paste0(prefix, "_hmm_coverage")
  query_cov_col <- paste0(prefix, "_query_coverage")

  tbl$module_family_id <- as.character(tbl[[family_col]])
  tbl$module_hit_label <- if (hit_label_col %in% names(tbl)) as.character(tbl[[hit_label_col]]) else NA_character_
  tbl$module_support <- if (support_col %in% names(tbl)) as.character(tbl[[support_col]]) else NA_character_
  tbl$module_evidence <- if (evidence_col %in% names(tbl)) as.character(tbl[[evidence_col]]) else NA_character_
  tbl$module_role <- if (role_col %in% names(tbl)) as.character(tbl[[role_col]]) else NA_character_
  tbl$module_subtype <- if (subtype_col %in% names(tbl)) as.character(tbl[[subtype_col]]) else NA_character_
  tbl$module_clan <- if (clan_col %in% names(tbl)) as.character(tbl[[clan_col]]) else NA_character_
  tbl$module_pident <- if (pident_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[pident_col]])) else NA_real_
  tbl$module_evalue <- if (evalue_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[evalue_col]])) else NA_real_
  tbl$module_hit_score <- if (score_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[score_col]])) else NA_real_
  tbl$module_hmm_coverage <- if (hmm_cov_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[hmm_cov_col]])) else NA_real_
  tbl$module_query_coverage <- if (query_cov_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[query_cov_col]])) else NA_real_

  keep <- !is.na(tbl$module_family_id) & nzchar(tbl$module_family_id)
  tbl <- tbl[keep, , drop = FALSE]
  if (!nrow(tbl)) {
    return(tbl)
  }

  tbl$module_display <- dplyr::coalesce(
    ifelse(!is.na(tbl$module_hit_label) & nzchar(tbl$module_hit_label), tbl$module_hit_label, NA_character_),
    tbl$module_family_id
  )
  tbl$module_group <- dplyr::coalesce(
    ifelse(!is.na(tbl$module_subtype) & nzchar(tbl$module_subtype), tbl$module_subtype, NA_character_),
    ifelse(!is.na(tbl$module_clan) & nzchar(tbl$module_clan), tbl$module_clan, NA_character_),
    ifelse(!is.na(tbl$module_role) & nzchar(tbl$module_role), tbl$module_role, NA_character_),
    tbl$module_display
  )
  tbl
}

.dnmb_antidefense_status_plot <- function(module_name, output_dir, status_text) {
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(.dnmb_antidefense_module_file_stub(module_name), "_overview.pdf"))
  p <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 0.5, y = 0.68, label = paste0(module_name, " anti-defense overview"), fontface = "bold", size = 6) +
    ggplot2::annotate("text", x = 0.5, y = 0.42, label = status_text, size = 4.8, color = "grey30") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")
  .dnmb_module_plot_save(p, pdf_path, width = 10, height = 4.8)
  list(pdf = pdf_path)
}

.dnmb_plot_antidefense_module <- function(genbank_table, output_dir, module_name) {
  tbl <- .dnmb_antidefense_module_hit_table(genbank_table, module_name)
  has_outputs <- .dnmb_antidefense_module_has_outputs(output_dir, module_name)
  if (!nrow(tbl)) {
    if (!has_outputs) {
      return(NULL)
    }
    return(.dnmb_antidefense_status_plot(
      module_name = module_name,
      output_dir = output_dir,
      status_text = "Module completed, but no anti-defense hits were detected for this genome."
    ))
  }

  palette <- .dnmb_defensefinder_palette(tbl$module_group)
  summary_tbl <- tbl |>
    dplyr::group_by(.data$module_group, .data$module_display) |>
    dplyr::summarise(
      n_hits = dplyr::n(),
      contigs = dplyr::n_distinct(.data$contig),
      best_pident = suppressWarnings(max(.data$module_pident, na.rm = TRUE)),
      best_evalue = suppressWarnings(min(.data$module_evalue, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_hits), .data$module_group, .data$module_display)

  summary_tbl$best_pident[!is.finite(summary_tbl$best_pident)] <- NA_real_
  summary_tbl$best_evalue[!is.finite(summary_tbl$best_evalue)] <- NA_real_
  summary_tbl$module_display <- factor(summary_tbl$module_display, levels = rev(unique(summary_tbl$module_display)))

  p_inventory <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = .data$n_hits, y = .data$module_display)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$module_group),
      width = 0.62, color = "grey35", linewidth = 0.2, alpha = 0.92
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$n_hits + 0.05,
        label = ifelse(
          is.na(.data$best_pident),
          paste0(.data$contigs, " contigs"),
          paste0(.data$contigs, " contigs | best pident ", format(round(.data$best_pident, 1), nsmall = 1))
        )
      ),
      hjust = 0, size = 2.7, color = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.35))) +
    ggplot2::labs(
      title = paste0("A   ", module_name, " inventory"),
      x = "Hit count",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.title.position = "plot"
    )

  layout_tbl <- tbl
  contig_levels <- unique(as.character(layout_tbl$contig))
  layout_tbl$contig <- factor(layout_tbl$contig, levels = rev(contig_levels))
  layout_tbl$start <- suppressWarnings(as.numeric(layout_tbl$start))
  layout_tbl$end <- suppressWarnings(as.numeric(layout_tbl$end))
  layout_tbl$mid <- (layout_tbl$start + layout_tbl$end) / 2
  layout_tbl$label_short <- .dnmb_integrated_defense_short_label(layout_tbl$module_display, max_len = 16L)

  p_layout <- ggplot2::ggplot(layout_tbl, ggplot2::aes(y = .data$contig)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$start, xend = .data$end, yend = .data$contig, color = .data$module_group),
      linewidth = 2.4, lineend = "round", alpha = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$mid, color = .data$module_group),
      size = 1.8, alpha = 0.95
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(x = .data$mid, label = .data$label_short, color = .data$module_group),
      size = 2.6,
      min.segment.length = 0,
      box.padding = 0.18,
      point.padding = 0.12,
      max.overlaps = 100
    ) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::labs(
      title = paste0("B   ", module_name, " genome layout"),
      x = "Genomic position (bp)",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.title.position = "plot",
      legend.position = "bottom"
    )

  top_tbl <- summary_tbl |>
    dplyr::mutate(
      best_evalue_label = ifelse(is.na(.data$best_evalue), "NA", format(signif(.data$best_evalue, 3), scientific = TRUE)),
      label = paste0(
        .data$module_display, " [", .data$module_group, "]",
        "\nHits: ", .data$n_hits,
        " | Contigs: ", .data$contigs,
        ifelse(is.na(.data$best_pident), "", paste0(" | Best pident: ", format(round(.data$best_pident, 1), nsmall = 1))),
        "\nBest e-value: ", .data$best_evalue_label
      )
    ) |>
    dplyr::slice_head(n = 12L)

  top_tbl$rank <- factor(seq_len(nrow(top_tbl)), levels = rev(seq_len(nrow(top_tbl))))
  p_hits <- ggplot2::ggplot(top_tbl, ggplot2::aes(y = .data$rank)) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0, label = .data$label),
      hjust = 0, vjust = 0.5, size = 3.0, lineheight = 1.08
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(title = paste0("C   Top ", module_name, " calls")) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(.dnmb_antidefense_module_file_stub(module_name), "_overview.pdf"))
  composite <- cowplot::plot_grid(
    p_inventory,
    p_layout,
    p_hits,
    ncol = 1,
    rel_heights = c(0.95, 1.45, 1.15),
    align = "v"
  )
  .dnmb_module_plot_save(composite, pdf_path, width = 10, height = 13.2)
  list(pdf = pdf_path)
}

.dnmb_plot_dbapis_module <- function(genbank_table, output_dir) {
  tbl <- .dnmb_antidefense_module_hit_table(genbank_table, "dbAPIS")
  has_outputs <- .dnmb_antidefense_module_has_outputs(output_dir, "dbAPIS")
  if (!nrow(tbl)) {
    if (!has_outputs) {
      return(NULL)
    }
    return(.dnmb_antidefense_status_plot(
      module_name = "dbAPIS",
      output_dir = output_dir,
      status_text = "Module completed, but no anti-defense hits were detected for this genome."
    ))
  }

  tbl <- tbl[!is.na(tbl$start) & !is.na(tbl$end), , drop = FALSE]
  if (!nrow(tbl)) {
    return(.dnmb_antidefense_status_plot(
      module_name = "dbAPIS",
      output_dir = output_dir,
      status_text = "dbAPIS hits were detected, but no genomic coordinates were available for plotting."
    ))
  }

  tbl$module_display <- as.character(tbl$module_display)
  tbl$module_group <- as.character(tbl$module_group)
  tbl$module_query_coverage <- suppressWarnings(as.numeric(tbl$module_query_coverage))
  tbl$module_hmm_coverage <- suppressWarnings(as.numeric(tbl$module_hmm_coverage))
  tbl$module_evalue <- suppressWarnings(as.numeric(tbl$module_evalue))
  tbl$module_hit_score <- suppressWarnings(as.numeric(tbl$module_hit_score))
  tbl$start <- suppressWarnings(as.numeric(tbl$start))
  tbl$end <- suppressWarnings(as.numeric(tbl$end))

  summary_tbl <- tbl |>
    dplyr::group_by(.data$module_display, .data$module_group) |>
    dplyr::summarise(
      n_proteins = dplyr::n(),
      n_contigs = dplyr::n_distinct(.data$contig),
      best_query_cov = max(.data$module_query_coverage, na.rm = TRUE),
      best_hmm_cov = max(.data$module_hmm_coverage, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_proteins), dplyr::desc(.data$best_query_cov), .data$module_display)

  summary_tbl$best_query_cov[!is.finite(summary_tbl$best_query_cov)] <- NA_real_
  summary_tbl$best_hmm_cov[!is.finite(summary_tbl$best_hmm_cov)] <- NA_real_
  summary_tbl$target_label <- .dnmb_antidefense_target_label(summary_tbl$module_group, max_len = 24L)
  summary_tbl$module_display <- factor(summary_tbl$module_display, levels = rev(unique(summary_tbl$module_display)))
  palette <- .dnmb_antidefense_display_palette(
    levels(summary_tbl$module_display),
    colors = c("#DDF6D2", "#BFEAAB", "#9CDB84", "#78CB60", "#52B846", "#2F9D31")
  )
  legend_labels <- stats::setNames(
    paste0(as.character(summary_tbl$module_display), " | ", summary_tbl$target_label),
    as.character(summary_tbl$module_display)
  )

  protein_label <- ifelse(summary_tbl$n_proteins == 1L, " protein", " proteins")
  inventory_detail <- ifelse(
    is.na(summary_tbl$best_query_cov),
    paste0(summary_tbl$n_proteins, protein_label, " | target ", summary_tbl$target_label),
    paste0(summary_tbl$n_proteins, protein_label, " | target ", summary_tbl$target_label, " | best cov ", sprintf("%.2f", summary_tbl$best_query_cov))
  )

  p_inventory <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = .data$n_proteins, y = .data$module_display)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$module_display),
      width = 0.62,
      color = "grey35",
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0.05, label = .data$module_display),
      hjust = 0,
      size = 3.0,
      color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$n_proteins + 0.08, label = inventory_detail),
      hjust = 0,
      size = 2.8,
      color = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.34))) +
    ggplot2::labs(title = NULL, x = "Hits detected", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  layout_tbl <- tbl |>
    dplyr::arrange(.data$contig, .data$start, .data$end) |>
    dplyr::group_by(.data$module_display) |>
    dplyr::mutate(
      display_rank = dplyr::row_number(),
      display_n = dplyr::n(),
      label = ifelse(.data$display_n > 1L, paste0(.data$module_display, " (", .data$display_rank, ")"), .data$module_display)
    ) |>
    dplyr::ungroup()

  contigs <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  hit_contigs <- unique(as.character(layout_tbl$contig))
  contigs <- contigs[contigs$contig %in% hit_contigs, , drop = FALSE]

  contig_labels <- stats::setNames(
    paste0(contigs$sector_label, " (", scales::label_comma()(contigs$length_bp), " bp)"),
    contigs$contig
  )
  layout_tbl$contig_facet <- factor(contig_labels[layout_tbl$contig], levels = unname(contig_labels))
  contigs$contig_facet <- factor(contig_labels[contigs$contig], levels = unname(contig_labels))
  layout_tbl$midpoint <- (layout_tbl$start + layout_tbl$end) / 2

  p_layout <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contigs,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = 0.5, yend = 0.5),
      linewidth = 0.6,
      color = "grey78"
    ) +
    ggplot2::geom_rect(
      data = layout_tbl,
      ggplot2::aes(
        xmin = .data$start, xmax = .data$end,
        ymin = 0.48, ymax = 0.52,
        fill = .data$module_display
      ),
      color = "grey35",
      linewidth = 0.2,
      alpha = 0.95
    ) +
    ggrepel::geom_text_repel(
      data = layout_tbl,
      ggplot2::aes(x = .data$midpoint, y = 0.52, label = .data$label),
      size = 1.9,
      lineheight = 0.85,
      direction = "y",
      nudge_y = 0.03,
      segment.size = 0.2,
      segment.color = "grey55",
      max.overlaps = Inf,
      seed = 42,
      min.segment.length = 0,
      force = 1.5,
      force_pull = 0.3,
      box.padding = 0.12,
      ylim = c(0.54, 0.72)
    ) +
    ggplot2::facet_wrap(~contig_facet, ncol = 1, scales = "free_x", strip.position = "top") +
    ggplot2::scale_fill_manual(values = palette, breaks = names(legend_labels), labels = unname(legend_labels[names(legend_labels)])) +
    ggplot2::scale_x_continuous(labels = scales::label_comma()) +
    ggplot2::scale_y_continuous(limits = c(0.44, 0.74), expand = c(0, 0)) +
    ggplot2::labs(title = "dbAPIS genome layout", x = "Genome coordinate (bp)", y = NULL, fill = NULL) +
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

  top_hits <- tbl |>
    dplyr::mutate(
      dbapis_label = paste0(.data$locus_tag, " | ", .data$product)
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$module_query_coverage),
      dplyr::desc(.data$module_hmm_coverage),
      .data$module_evalue,
      .data$locus_tag
    ) |>
    dplyr::slice_head(n = 15L)
  top_hits$dbapis_label <- factor(top_hits$dbapis_label, levels = rev(top_hits$dbapis_label))

  right_label <- ifelse(
    is.na(top_hits$module_hmm_coverage),
    .dnmb_antidefense_target_label(top_hits$module_group, max_len = 24L),
    paste0(.dnmb_antidefense_target_label(top_hits$module_group, max_len = 24L), " | hmm ", sprintf("%.2f", top_hits$module_hmm_coverage))
  )

  max_query_cov <- suppressWarnings(max(top_hits$module_query_coverage, na.rm = TRUE))
  if (!is.finite(max_query_cov)) {
    max_query_cov <- 1
  }

  p_hits <- ggplot2::ggplot(top_hits, ggplot2::aes(x = .data$module_query_coverage, y = .data$dbapis_label)) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$module_display),
      width = 0.62,
      color = "grey35",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$module_query_coverage + 0.02, label = right_label),
      hjust = 0,
      size = 2.7,
      color = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(limits = c(0, max(1.05, max_query_cov + 0.22))) +
    ggplot2::labs(title = NULL, x = "Query coverage", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "none"
    )

  legend_rows <- if (length(palette) > 8L) 2L else 1L
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
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, nrow = legend_rows, byrow = TRUE))
  common_legend <- cowplot::get_legend(legend_plot)
  legend_panel <- cowplot::ggdraw() +
    cowplot::draw_grob(common_legend, x = 0.5, y = 0.5, width = 0.98, height = 0.95, hjust = 0.5, vjust = 0.5)

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "dbAPIS_overview.pdf")
  n_contigs <- nrow(contigs)
  layout_height <- max(0.65, 0.40 * n_contigs + 0.25)
  composite <- cowplot::plot_grid(
    .dnmb_antidefense_panel_header(p_inventory, "A", "dbAPIS inventory", plot_y = 0.06, plot_h = 0.86),
    .dnmb_antidefense_panel_header(p_layout, "B", "dbAPIS genome layout", plot_y = 0.035, plot_h = 0.90),
    .dnmb_antidefense_panel_header(p_hits, "C", "Top dbAPIS protein calls", plot_y = 0.06, plot_h = 0.86),
    legend_panel,
    ncol = 1,
    rel_heights = c(0.86, layout_height, 1.10, 0.30)
  )
  total_h <- max(10.2, 1.9 + layout_height * 2.15 + 3.1)
  .dnmb_module_plot_save(composite, pdf_path, width = 10, height = total_h)
  list(pdf = pdf_path)
}

.dnmb_plot_acrfinder_module <- function(genbank_table, output_dir) {
  .dnmb_plot_antidefense_module(genbank_table, output_dir, module_name = "AcrFinder")
}

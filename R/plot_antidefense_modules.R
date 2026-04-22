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

  tbl$module_family_id <- as.character(tbl[[family_col]])
  tbl$module_hit_label <- if (hit_label_col %in% names(tbl)) as.character(tbl[[hit_label_col]]) else NA_character_
  tbl$module_support <- if (support_col %in% names(tbl)) as.character(tbl[[support_col]]) else NA_character_
  tbl$module_evidence <- if (evidence_col %in% names(tbl)) as.character(tbl[[evidence_col]]) else NA_character_
  tbl$module_role <- if (role_col %in% names(tbl)) as.character(tbl[[role_col]]) else NA_character_
  tbl$module_subtype <- if (subtype_col %in% names(tbl)) as.character(tbl[[subtype_col]]) else NA_character_
  tbl$module_clan <- if (clan_col %in% names(tbl)) as.character(tbl[[clan_col]]) else NA_character_
  tbl$module_pident <- if (pident_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[pident_col]])) else NA_real_
  tbl$module_evalue <- if (evalue_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[evalue_col]])) else NA_real_

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
  .dnmb_plot_antidefense_module(genbank_table, output_dir, module_name = "dbAPIS")
}

.dnmb_plot_acrfinder_module <- function(genbank_table, output_dir) {
  .dnmb_plot_antidefense_module(genbank_table, output_dir, module_name = "AcrFinder")
}

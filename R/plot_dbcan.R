.dnmb_plot_dbcan_module <- function(genbank_table, output_dir) {
  tbl <- as.data.frame(genbank_table, stringsAsFactors = FALSE)

  cgc_id_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_id", "dbCAN_cgc_id", "dbcan_cgc_id"))
  cgc_type_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_gene_type", "dbCAN_cgc_gene_type", "dbcan_cgc_gene_type"))

  if (is.null(cgc_id_col) || is.null(cgc_type_col)) {
    cgc_path <- file.path(output_dir, "dnmb_module_dbcan", "run_dbcan", "cgc_standard_out.tsv")
    if (!file.exists(cgc_path) || file.info(cgc_path)$size < 100) {
      stop("dbCAN plot: no CGC data found — neither in genbank_table columns ",
           "nor in dnmb_module_dbcan/run_dbcan/cgc_standard_out.tsv.",
           call. = FALSE)
    }
    cgc_raw <- utils::read.delim(cgc_path, sep = "\t", header = TRUE,
                                  stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(cgc_raw)) {
      stop("dbCAN plot: cgc_standard_out.tsv exists but has no data rows.",
           call. = FALSE)
    }
    col_map <- c("Protein ID" = "locus_tag", "CGC#" = "cgc_num",
                 "Contig ID" = "contig", "Gene Type" = "dbcan_cgc_gene_type",
                 "Gene Start" = "start", "Gene Stop" = "end",
                 "Gene Strand" = "direction", "Gene Annotation" = "dbcan_cgc_protein_family")
    for (old in names(col_map)) {
      if (old %in% names(cgc_raw)) names(cgc_raw)[names(cgc_raw) == old] <- col_map[old]
    }
    cgc_raw$start <- suppressWarnings(as.numeric(cgc_raw$start))
    cgc_raw$end <- suppressWarnings(as.numeric(cgc_raw$end))
    cgc_raw <- cgc_raw[!is.na(cgc_raw$start) & !is.na(cgc_raw$end), , drop = FALSE]
    if ("contig" %in% names(cgc_raw) && "cgc_num" %in% names(cgc_raw)) {
      cgc_raw$dbcan_cgc_id <- paste(cgc_raw$contig, cgc_raw$cgc_num, sep = "|")
    }
    family_col <- .dnmb_pick_column(tbl, c("dbCAN_family_id", "family_id"))
    if (!is.null(family_col) && "locus_tag" %in% names(cgc_raw)) {
      fam_map <- stats::setNames(as.character(tbl[[family_col]]), tbl$locus_tag)
      cgc_raw$dbCAN_family_id <- unname(fam_map[cgc_raw$locus_tag])
    } else {
      cgc_raw$dbCAN_family_id <- NA_character_
    }
    cgc_raw$dbcan_pul_substrate <- NA_character_
    tbl <- cgc_raw[!is.na(cgc_raw$dbcan_cgc_id) & nzchar(cgc_raw$dbcan_cgc_id), , drop = FALSE]
  } else {
    cgc_family_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_protein_family", "dbCAN_cgc_protein_family", "dbcan_cgc_protein_family"))
    substrate_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_pul_substrate", "dbCAN_pul_substrate", "dbcan_pul_substrate"))
    family_col <- .dnmb_pick_column(tbl, c("dbCAN_family_id", "family_id"))
    tbl$start <- suppressWarnings(as.numeric(tbl$start))
    tbl$end <- suppressWarnings(as.numeric(tbl$end))
    tbl <- tbl[!is.na(tbl$start) & !is.na(tbl$end), , drop = FALSE]
    tbl$dbcan_cgc_id <- as.character(tbl[[cgc_id_col]])
    tbl$dbcan_cgc_gene_type <- as.character(tbl[[cgc_type_col]])
    tbl$dbcan_cgc_protein_family <- if (!is.null(cgc_family_col)) as.character(tbl[[cgc_family_col]]) else NA_character_
    tbl$dbcan_pul_substrate <- if (!is.null(substrate_col)) as.character(tbl[[substrate_col]]) else NA_character_
    tbl$dbCAN_family_id <- if (!is.null(family_col)) as.character(tbl[[family_col]]) else NA_character_
    tbl <- tbl[!is.na(tbl$dbcan_cgc_id) & nzchar(tbl$dbcan_cgc_id), , drop = FALSE]
  }

  if (!nrow(tbl)) {
    stop("dbCAN plot: no CGC data to plot.", call. = FALSE)
  }

  known_types <- c("CAZyme", "TC", "TF", "STP")
  tbl$gene_type <- ifelse(
    is.na(tbl$dbcan_cgc_gene_type) | !nzchar(tbl$dbcan_cgc_gene_type) |
      !tbl$dbcan_cgc_gene_type %in% known_types,
    "other", tbl$dbcan_cgc_gene_type
  )

  type_palette <- c(
    CAZyme = "#0F766E",
    TC     = "#D97706",
    TF     = "#7C3AED",
    STP    = "#DC2626",
    other  = "#9CA3AF"
  )

  # --- Build CGC summary ---
  cgc_summary <- tbl |>
    dplyr::group_by(.data$dbcan_cgc_id) |>
    dplyr::summarise(
      contig = dplyr::first(.data$contig),
      n_genes = dplyr::n(),
      n_cazyme = sum(.data$gene_type == "CAZyme", na.rm = TRUE),
      cazyme_families = paste(
        sort(unique(stats::na.omit(.data$dbCAN_family_id[.data$gene_type == "CAZyme"]))),
        collapse = ", "
      ),
      substrate = dplyr::first(stats::na.omit(.data$dbcan_pul_substrate)),
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      gene_types = paste(sort(unique(.data$gene_type[.data$gene_type != "other" & .data$gene_type != "null"])), collapse = "+"),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_cazyme), dplyr::desc(.data$n_genes))

  cgc_summary$short_id <- sub("^.*\\|", "", cgc_summary$dbcan_cgc_id)

  # --- Panel A: Genome layout (full width, top) ---
  contigs <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  hit_contigs <- unique(cgc_summary$contig)
  contigs <- contigs[contigs$contig %in% hit_contigs, , drop = FALSE]

  cgc_summary$cgc_label <- ifelse(
    !is.na(cgc_summary$substrate) & nzchar(cgc_summary$substrate),
    paste0(cgc_summary$short_id, "\n", cgc_summary$substrate),
    cgc_summary$short_id
  )
  cgc_summary$midpoint <- (cgc_summary$start + cgc_summary$end) / 2

  contig_labels <- stats::setNames(
    paste0(contigs$sector_label, " (", scales::label_comma()(contigs$length_bp), " bp)"),
    contigs$contig
  )
  cgc_summary$contig_facet <- factor(contig_labels[cgc_summary$contig], levels = unname(contig_labels))
  contigs$contig_facet <- factor(contig_labels[contigs$contig], levels = unname(contig_labels))

  p_layout <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contigs,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = 0.5, yend = 0.5),
      linewidth = 0.8, color = "grey78"
    ) +
    ggplot2::geom_rect(
      data = cgc_summary,
      ggplot2::aes(
        xmin = .data$start, xmax = .data$end,
        ymin = 0.35, ymax = 0.65,
        fill = .data$n_cazyme
      ),
      color = "grey35", linewidth = 0.25, alpha = 0.95
    ) +
    ggrepel::geom_text_repel(
      data = cgc_summary,
      ggplot2::aes(x = .data$midpoint, y = 0.74, label = .data$cgc_label),
      size = 2.3, vjust = 0, lineheight = 0.85,
      direction = "x", nudge_y = 0.15,
      segment.size = 0.2, segment.color = "grey60",
      max.overlaps = 50, seed = 42, min.segment.length = 0.1
    ) +
    ggplot2::facet_wrap(~contig_facet, ncol = 1, scales = "free_x", strip.position = "top") +
    ggplot2::scale_fill_gradient(low = "#81C784", high = "#0F766E", name = "CAZymes") +
    ggplot2::scale_x_continuous(labels = scales::label_comma()) +
    ggplot2::labs(title = "dbCAN CGC genome layout", x = "Genome coordinate (bp)", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot",
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 9)
    )

  # --- Panel B: Combined inventory + gene type composition (compact) ---
  gene_type_counts <- tbl |>
    dplyr::group_by(.data$dbcan_cgc_id, .data$gene_type) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  gene_type_counts$short_id <- sub("^.*\\|", "", gene_type_counts$dbcan_cgc_id)
  cgc_order <- cgc_summary$short_id
  gene_type_counts$short_id <- factor(gene_type_counts$short_id, levels = rev(cgc_order))

  annot_map <- stats::setNames(
    vapply(seq_len(nrow(cgc_summary)), function(i) {
      parts <- c(
        if (nzchar(cgc_summary$cazyme_families[i])) cgc_summary$cazyme_families[i] else NULL,
        if (!is.na(cgc_summary$substrate[i]) && nzchar(cgc_summary$substrate[i]))
          cgc_summary$substrate[i] else NULL
      )
      if (length(parts)) paste(parts, collapse = " | ") else ""
    }, character(1)),
    cgc_summary$short_id
  )
  annot_df <- cgc_summary |>
    dplyr::mutate(
      short_id = factor(.data$short_id, levels = rev(cgc_order)),
      annot_text = unname(annot_map[as.character(.data$short_id)])
    )

  p_combined <- ggplot2::ggplot(gene_type_counts, ggplot2::aes(x = .data$n, y = .data$short_id, fill = .data$gene_type)) +
    ggplot2::geom_col(width = 0.7, color = "grey35", linewidth = 0.12) +
    ggplot2::geom_text(
      data = annot_df,
      ggplot2::aes(x = .data$n_genes + 0.3, y = .data$short_id, label = .data$annot_text),
      inherit.aes = FALSE,
      hjust = 0, size = 2.0, color = "grey30"
    ) +
    ggplot2::scale_fill_manual(values = type_palette, drop = FALSE, name = "Gene type") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.35))) +
    ggplot2::labs(
      title = paste0("CGC gene composition (", nrow(cgc_summary), " clusters)"),
      x = "Number of genes", y = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.title.position = "plot",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.3, "cm"),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  # --- Panel C: Substrate-CAZyme family network ---
  p_network <- .dnmb_dbcan_substrate_cazyme_plot(cgc_summary, tbl, type_palette)

  # --- Composite: A top, B+C bottom side by side ---
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "dbcan_cgc_overview.pdf")
  n_contigs <- nrow(contigs)
  n_cgc <- nrow(cgc_summary)

  bottom_row <- cowplot::plot_grid(
    p_combined, p_network,
    labels = c("B", "C"),
    label_size = 14, label_fontface = "bold",
    label_x = 0, label_y = c(1.02, 1.02), hjust = 0,
    ncol = 2, rel_widths = c(1.1, 0.9)
  )

  composite <- cowplot::plot_grid(
    p_layout,
    bottom_row,
    labels = c("A", ""),
    label_size = 14, label_fontface = "bold",
    label_x = 0, label_y = c(1.02, 1), hjust = 0,
    ncol = 1,
    rel_heights = c(0.30, 0.70)
  )

  page_h <- max(11, min(16, 3 + 0.22 * n_cgc))
  .dnmb_module_plot_save(composite, pdf_path, width = 14, height = page_h)
  list(pdf = pdf_path)
}


# Substrate → CAZyme family alluvial / tile plot
.dnmb_dbcan_substrate_cazyme_plot <- function(cgc_summary, tbl, type_palette) {
  caz_genes <- tbl[tbl$gene_type == "CAZyme" & !is.na(tbl$dbCAN_family_id) & nzchar(tbl$dbCAN_family_id), , drop = FALSE]
  if (!nrow(caz_genes)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No CAZyme families detected",
                               size = 3.5, color = "grey50"))
  }

  caz_genes$short_id <- sub("^.*\\|", "", caz_genes$dbcan_cgc_id)
  sub_map <- stats::setNames(cgc_summary$substrate, cgc_summary$short_id)
  caz_genes$substrate <- sub_map[caz_genes$short_id]
  caz_genes$substrate <- ifelse(
    is.na(caz_genes$substrate) | !nzchar(caz_genes$substrate),
    "unknown", caz_genes$substrate
  )

  links <- caz_genes |>
    dplyr::group_by(.data$substrate, .data$dbCAN_family_id) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$n))

  all_substrates <- unique(links$substrate)
  all_families <- unique(links$dbCAN_family_id)

  sub_pal <- stats::setNames(
    c(grDevices::hcl.colors(max(length(all_substrates), 3), palette = "Set 2"))[seq_along(all_substrates)],
    all_substrates
  )
  if ("unknown" %in% names(sub_pal)) sub_pal["unknown"] <- "#BDBDBD"

  # Tile heatmap: substrate (x) vs CAZyme family (y)
  links$substrate <- factor(links$substrate, levels = c(setdiff(all_substrates, "unknown"), "unknown"))
  fam_order <- links |>
    dplyr::group_by(.data$dbCAN_family_id) |>
    dplyr::summarise(total = sum(.data$n), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$total))
  links$dbCAN_family_id <- factor(links$dbCAN_family_id, levels = rev(fam_order$dbCAN_family_id))

  ggplot2::ggplot(links, ggplot2::aes(x = .data$substrate, y = .data$dbCAN_family_id)) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$n, fill = .data$substrate),
      shape = 21, color = "grey30", stroke = 0.3, alpha = 0.85
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$n),
      size = 2.3, color = "grey20"
    ) +
    ggplot2::scale_size_continuous(range = c(3, 12), name = "Genes", guide = "none") +
    ggplot2::scale_fill_manual(values = sub_pal, name = "Substrate", guide = "none") +
    ggplot2::labs(
      title = "Substrate \u2013 CAZyme family",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.title.position = "plot",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(4, 8, 4, 4)
    )
}

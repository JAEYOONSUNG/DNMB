.dnmb_prophage_placeholder_plot <- function(output_dir) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5,
      label = "No prophage regions detected",
      size = 6, color = "#9CA3AF", fontface = "italic") +
    ggplot2::annotate("text", x = 0.5, y = 0.35,
      label = "PhiSpy did not identify any prophage-like regions in this genome",
      size = 3.5, color = "#D1D5DB") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(title = "Prophage Region Overview") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 14),
      plot.margin = ggplot2::margin(20, 20, 20, 20)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "Prophage_overview.pdf")
  ggplot2::ggsave(pdf_path, p, width = 14, height = 4, bg = "white")
  list(pdf = pdf_path)
}

.dnmb_plot_prophage_module <- function(genbank_table, output_dir, cache_root = NULL) {
  tbl <- .dnmb_contig_ordered_table(genbank_table)

  # --- Case 1: No Prophage column at all ---
  if (!base::nrow(tbl) || !"Prophage_prophage_id" %in% base::names(tbl)) {
    return(.dnmb_prophage_placeholder_plot(output_dir))
  }

  tbl <- tbl[!is.na(tbl$Prophage_prophage_id) & base::nzchar(tbl$Prophage_prophage_id) & tbl$Prophage_prophage_id != "0", , drop = FALSE]

  # --- Case 2: Column exists but no prophage genes detected ---
  if (!base::nrow(tbl)) {
    return(.dnmb_prophage_placeholder_plot(output_dir))
  }

  # --- Classify regions as complete vs partial/questionable ---
  # Criteria: gene count < 5 OR mean PhiSpy rank score < 5 -> partial
  region_stats <- tbl |>
    dplyr::group_by(.data$Prophage_prophage_id) |>
    dplyr::summarise(
      n_genes = dplyr::n(),
      mean_rank = mean(suppressWarnings(as.numeric(.data$Prophage_prophage_rank)), na.rm = TRUE),
      mean_pp = mean(suppressWarnings(as.numeric(.data$Prophage_prophage_pp)), na.rm = TRUE),
      .groups = "drop"
    )
  region_stats$region_class <- dplyr::case_when(
    region_stats$n_genes < 5                           ~ "partial",
    !is.na(region_stats$mean_rank) & region_stats$mean_rank < 5 &
      !is.na(region_stats$mean_pp) & region_stats$mean_pp < 0.3 ~ "partial",
    TRUE                                                ~ "complete"
  )
  # Merge classification back to tbl
  tbl <- merge(tbl, region_stats[, c("Prophage_prophage_id", "region_class")],
               by = "Prophage_prophage_id", all.x = TRUE)

  tbl$category <- .dnmb_gene_arrow_category_prophage(tbl$product)
  label_tbl <- tbl[grepl("capsid|tail|portal|terminase|integrase|lysin|holin|helicase|primase|polymerase|recombinase|anti.crispr|methyltransferase", tolower(tbl$product)), , drop = FALSE]
  label_tbl$label <- gsub(" family.*$", "", label_tbl$product)
  region_tbl <- tbl |>
    dplyr::group_by(.data$contig, .data$Prophage_prophage_id) |>
    dplyr::summarise(
      region_start = min(.data$start, na.rm = TRUE),
      region_end = max(.data$end, na.rm = TRUE),
      .groups = "drop"
    )
  gene_tbl <- tbl
  gene_tbl$panel <- paste0("Prophage ", gene_tbl$Prophage_prophage_id)
  contig_lengths <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  top_regions <- unique(as.character(gene_tbl$Prophage_prophage_id))
  # Pre-compute att sites for annotation on storyboard
  gbff_for_att <- .dnmb_find_gbff_for_plot(output_dir)
  att_summary <- tryCatch(
    .dnmb_prophage_build_summary_table(genbank_table, gbff_path = gbff_for_att, output_dir = output_dir),
    error = function(e) data.frame()
  )

  region_panels <- lapply(top_regions, function(region_id) {
    sub_tbl <- gene_tbl[gene_tbl$Prophage_prophage_id == region_id, , drop = FALSE]
    sub_labels <- label_tbl[label_tbl$Prophage_prophage_id == region_id, , drop = FALSE]
    contig_id <- as.character(sub_tbl$contig[[1]])
    contig_length <- contig_lengths$length_bp[match(contig_id, contig_lengths$contig)]
    # Determine region classification
    rc <- region_stats$region_class[region_stats$Prophage_prophage_id == region_id]
    rc <- if (length(rc) == 1) rc else "complete"
    sub_tbl$panel <- paste0("Prophage ", region_id,
                            if (rc == "partial") " (partial/questionable)" else "")
    sub_tbl$label_short <- NA_character_
    if (nrow(sub_labels)) {
      label_idx <- match(sub_tbl$locus_tag, sub_labels$locus_tag)
      sub_tbl$label_short[!is.na(label_idx)] <- sub_labels$label[label_idx[!is.na(label_idx)]]
    }
    p <- .dnmb_prophage_storyboard_plot(
      sub_tbl = sub_tbl,
      contig_length = contig_length,
      region_id = region_id,
      palette = .dnmb_gene_arrow_palette_prophage(),
      show_legend = FALSE,
      region_class = rc
    )
    # Add att site markers + summary annotation if available
    if (!is.null(p) && is.data.frame(att_summary) && nrow(att_summary)) {
      att_row <- att_summary[att_summary$region_id == region_id, , drop = FALSE]
      if (nrow(att_row)) {
        pp_s <- att_row$start[[1]]
        pp_e <- att_row$end[[1]]
        zoom_start <- min(pmin(as.numeric(sub_tbl$start), as.numeric(sub_tbl$end)), na.rm = TRUE)
        zoom_end <- max(pmax(as.numeric(sub_tbl$start), as.numeric(sub_tbl$end)), na.rm = TRUE)
        gene_left <- -5.72
        gene_right <- 5.72
        gene_center_y <- 1.72

        # att site markers
        has_att <- !is.na(att_row$att_core_seq[[1]]) && nzchar(att_row$att_core_seq[[1]])
        if (has_att) {
          att_seq <- att_row$att_core_seq[[1]]
          attL_x <- .dnmb_linear_map_x(pp_s, xmin = zoom_start, xmax = zoom_end, x_min_plot = gene_left, x_max_plot = gene_right)
          attR_x <- .dnmb_linear_map_x(pp_e, xmin = zoom_start, xmax = zoom_end, x_min_plot = gene_left, x_max_plot = gene_right)
          att_marker_tbl <- data.frame(
            x = c(attL_x, attR_x),
            y = gene_center_y - 0.07,
            label = c("attL", "attR"),
            stringsAsFactors = FALSE
          )
          p <- p +
            ggplot2::geom_point(data = att_marker_tbl, ggplot2::aes(x = .data$x, y = .data$y),
                                shape = 25, size = 3, fill = "#DC2626", color = "#DC2626", inherit.aes = FALSE) +
            ggplot2::geom_text(data = att_marker_tbl, ggplot2::aes(x = .data$x, y = .data$y - 0.05, label = .data$label),
                               size = 2.2, vjust = 1, color = "#DC2626", fontface = "bold", inherit.aes = FALSE)
        }

        # Summary annotation line below storyboard
        compl <- att_row$completeness_proxy[[1]] %||% "?"
        conf <- att_row$confidence_tier[[1]] %||% "?"
        gc_info <- if (is.finite(as.numeric(att_row$delta_gc[[1]] %||% NA))) {
          sprintf("dGC=%+.1f%%", as.numeric(att_row$delta_gc[[1]]) * 100)
        } else ""
        att_validated <- !is.na(att_row$att_related_genome[[1]] %||% NA) && nzchar(att_row$att_related_genome[[1]] %||% "")
        att_info <- if (has_att) {
          paste0("att: ", att_seq, " (", att_row$att_core_length[[1]], "bp",
                 if (att_validated) paste0(", validated: ", att_row$att_related_genome[[1]]) else "",
                 ")")
        } else "att: not detected"
        int_info <- if (!is.na(att_row$integrase_family[[1]] %||% NA)) {
          paste0(att_row$integrase_family[[1]], " [", att_row$integrase_position[[1]], "]")
        } else ""
        trna_info <- if (!is.na(att_row$right_flank_genes[[1]] %||% NA)) {
          paste0("target: ", att_row$right_flank_genes[[1]])
        } else ""
        func_info <- att_row$functional_summary[[1]] %||% ""

        summary_line <- paste(
          c(
            paste0(toupper(compl), " | ", toupper(conf)),
            gc_info,
            att_info,
            int_info,
            trna_info
          )[nzchar(c(
            paste0(toupper(compl), " | ", toupper(conf)),
            gc_info, att_info, int_info, trna_info
          ))],
          collapse = "  |  "
        )

        p <- p +
          ggplot2::annotate("text", x = gene_left, y = 0.88,
                           label = summary_line, hjust = 0, size = 2.6,
                           color = "#4B5563", fontface = "italic")
      }
    }
    p
  })
  score_tbl <- tbl
  score_tbl$panel_label <- paste0("Prophage ", score_tbl$Prophage_prophage_id)
  score_tbl$plot_start <- pmin(suppressWarnings(as.numeric(score_tbl$start)), suppressWarnings(as.numeric(score_tbl$end)), na.rm = TRUE)
  score_tbl$plot_end <- pmax(suppressWarnings(as.numeric(score_tbl$start)), suppressWarnings(as.numeric(score_tbl$end)), na.rm = TRUE)
  score_tbl$plot_mid <- (score_tbl$plot_start + score_tbl$plot_end) / 2
  score_tbl$region_start <- if ("Prophage_prophage_start" %in% names(score_tbl)) suppressWarnings(as.numeric(score_tbl$Prophage_prophage_start)) else NA_real_
  score_tbl$region_end <- if ("Prophage_prophage_end" %in% names(score_tbl)) suppressWarnings(as.numeric(score_tbl$Prophage_prophage_end)) else NA_real_
  missing_region_bounds <- is.na(score_tbl$region_start) | is.na(score_tbl$region_end)
  if (any(missing_region_bounds)) {
    region_bounds <- score_tbl |>
      dplyr::group_by(.data$panel_label) |>
      dplyr::summarise(
        region_start_fill = min(.data$start, na.rm = TRUE),
        region_end_fill = max(.data$end, na.rm = TRUE),
        .groups = "drop"
      )
    score_tbl <- dplyr::left_join(score_tbl, region_bounds, by = "panel_label")
    score_tbl$region_start[is.na(score_tbl$region_start)] <- score_tbl$region_start_fill[is.na(score_tbl$region_start)]
    score_tbl$region_end[is.na(score_tbl$region_end)] <- score_tbl$region_end_fill[is.na(score_tbl$region_end)]
    score_tbl$region_start_fill <- NULL
    score_tbl$region_end_fill <- NULL
  }
  b_x_left <- -5.72
  b_x_right <- 5.72
  score_tbl$x_mid <- .dnmb_linear_map_x(score_tbl$plot_mid, xmin = score_tbl$region_start, xmax = score_tbl$region_end, x_min_plot = b_x_left, x_max_plot = b_x_right)
  score_tbl$x_start <- .dnmb_linear_map_x(score_tbl$plot_start, xmin = score_tbl$region_start, xmax = score_tbl$region_end, x_min_plot = b_x_left, x_max_plot = b_x_right)
  score_tbl$x_end <- .dnmb_linear_map_x(score_tbl$plot_end, xmin = score_tbl$region_start, xmax = score_tbl$region_end, x_min_plot = b_x_left, x_max_plot = b_x_right)
  score_tbl <- score_tbl[order(score_tbl$panel_label, score_tbl$x_mid), , drop = FALSE]
  score_tbl$rank_scaled <- suppressWarnings(as.numeric(score_tbl$Prophage_prophage_rank))
  score_tbl$pp_scaled <- suppressWarnings(as.numeric(score_tbl$Prophage_prophage_pp))
  score_tbl$my_status_flag <- ifelse(
    suppressWarnings(as.numeric(score_tbl$Prophage_prophage_my_status)) > 0,
    "1",
    "0"
  )
  prophage_palette <- .dnmb_gene_arrow_palette_prophage()
  score_tbl$category <- .dnmb_gene_arrow_category_prophage(score_tbl$product)
  score_tbl$tile_border <- unname(prophage_palette[score_tbl$category])
  score_tbl$tile_border[is.na(score_tbl$tile_border)] <- "#9CA3AF"
  score_tbl$my_status_ymin <- -0.12
  score_tbl$my_status_ymax <- 0
  p_score <- ggplot2::ggplot(score_tbl, ggplot2::aes(x = .data$x_mid)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$x_start, xmax = .data$x_end, ymin = .data$my_status_ymin, ymax = .data$my_status_ymax, fill = .data$my_status_flag),
      color = score_tbl$tile_border,
      linewidth = 0.28,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_area(
      ggplot2::aes(y = .data$rank_scaled, group = .data$panel_label),
      fill = grDevices::adjustcolor("#60A5FA", alpha.f = 0.24),
      color = NA
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, color = "#9CA3AF") +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$rank_scaled, color = "Rank", group = .data$panel_label),
      linewidth = 1.0,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$rank_scaled, color = "Rank"),
      size = 1.8,
      stroke = 0
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = .data$x_mid,
        y = 0,
        yend = .data$pp_scaled,
        color = "PP"
      ),
      linewidth = 0.8,
      alpha = 0.75,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$pp_scaled, color = "PP"),
      size = 2.1,
      stroke = 0
    ) +
    ggplot2::facet_wrap(~panel_label, scales = "fixed", ncol = 1) +
    ggplot2::labs(
      title = "Prophage Score Profile",
      x = NULL,
      y = "PhiSpy score",
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.35),
      strip.background = ggplot2::element_rect(fill = "#F3F4F6", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(color = "#9CA3AF", linewidth = 0.35),
      legend.position = "right",
      plot.margin = ggplot2::margin(t = 8, r = 180, b = 12, l = 12)
    ) +
    ggplot2::scale_color_manual(values = c("Rank" = "#2563EB", "PP" = "#DC2626")) +
    ggplot2::scale_fill_manual(values = c("0" = "#E5E7EB", "1" = "#FCA5A5"), guide = "none") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::annotate(
      "text",
      x = b_x_right + 0.08,
      y = -0.06,
      label = "RF status",
      hjust = 0,
      vjust = 0.5,
      size = 3.4,
      color = "#6B7280"
    ) +
    ggplot2::coord_cartesian(xlim = c(b_x_left, b_x_right), ylim = c(-0.12, NA), clip = "off")
  p_reference <- NULL
  if (length(top_regions) >= 1L) {
    top_region_tbl <- gene_tbl[gene_tbl$Prophage_prophage_id == top_regions[[1]], , drop = FALSE]
    p_reference <- tryCatch(
      .dnmb_plot_prophage_reference_panel(top_region_tbl, region_id = top_regions[[1]], output_dir = output_dir, top_n = 5L, cache_root = cache_root),
      error = function(e) NULL
    )
  }
  p_summary <- NULL
  gbff_for_summary <- .dnmb_find_gbff_for_plot(output_dir)
  summary_tbl <- tryCatch(
    .dnmb_prophage_build_summary_table(genbank_table, gbff_path = gbff_for_summary, output_dir = output_dir),
    error = function(e) data.frame()
  )
  if (is.data.frame(summary_tbl) && nrow(summary_tbl)) {
    palette_d <- .dnmb_gene_arrow_palette_prophage()
    summary_panels <- lapply(seq_len(nrow(summary_tbl)), function(ri) {
      s <- summary_tbl[ri, ]
      # Functional composition bar data
      comp_data <- data.frame(
        category = c("Head/packaging", "Tail", "Integration", "Lysis",
                     "DNA replication", "Recombination", "Anti-defense",
                     "Regulation", "Hypothetical", "Other"),
        count = c(
          if ("structural_count" %in% names(s)) ceiling(as.numeric(s$structural_count) * 0.6) else 0L,
          if ("structural_count" %in% names(s)) floor(as.numeric(s$structural_count) * 0.4) else 0L,
          if (isTRUE(s$integrase_present)) 1L else 0L,
          as.integer(s$lysis_count %||% 0),
          as.integer(s$replication_count %||% 0),
          as.integer(s$recombination_count %||% 0),
          as.integer(s$anti_defense_count %||% 0),
          as.integer(s$regulation_count %||% 0),
          round(as.numeric(s$hypothetical_frac %||% 0) * as.integer(s$gene_count %||% 0)),
          0L
        ),
        stringsAsFactors = FALSE
      )
      comp_data$count[comp_data$category == "Other"] <- max(0L, as.integer(s$gene_count) - sum(comp_data$count))
      comp_data <- comp_data[comp_data$count > 0, , drop = FALSE]
      comp_data$fill <- unname(palette_d[comp_data$category])
      comp_data$xmax <- cumsum(comp_data$count)
      comp_data$xmin <- c(0, comp_data$xmax[-nrow(comp_data)])
      total_genes <- as.integer(s$gene_count)

      # Badge colors
      compl_col <- switch(as.character(s$completeness_proxy),
                          "intact" = "#22C55E", "questionable" = "#F59E0B", "#EF4444")
      conf_col <- switch(as.character(s$confidence_tier),
                         "high" = "#22C55E", "medium" = "#F59E0B", "#EF4444")
      att_col <- if (isTRUE(s$att_site_detected)) "#22C55E" else "#D1D5DB"

      # GC deviation bar
      gc_region <- as.numeric(s$region_gc %||% NA)
      gc_host <- as.numeric(s$host_gc %||% NA)
      gc_delta <- as.numeric(s$delta_gc %||% NA)

      # Metrics text
      len_label <- .dnmb_fmt_bp_label(as.numeric(s$length_bp))
      pgf_label <- paste0(round(as.numeric(s$phage_gene_frac %||% 0) * 100), "%")

      y_bar <- 0.72
      y_badges <- 0.32
      y_metrics <- 0.32

      p <- ggplot2::ggplot() +
        # Composition bar
        ggplot2::geom_rect(
          data = comp_data,
          ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = y_bar - 0.14, ymax = y_bar + 0.14, fill = .data$fill),
          color = "white", linewidth = 0.3, inherit.aes = FALSE
        ) +
        ggplot2::annotate("text", x = -1, y = y_bar, label = "Gene composition", hjust = 1, size = 3.2, fontface = "bold", color = "#374151") +
        # Composition labels inside bars (only if wide enough)
        ggplot2::geom_text(
          data = comp_data[comp_data$count >= 2, , drop = FALSE],
          ggplot2::aes(x = (.data$xmin + .data$xmax) / 2, y = y_bar, label = .data$count),
          size = 2.5, color = "white", fontface = "bold", inherit.aes = FALSE
        ) +
        # Badges row
        ggplot2::annotate("rect", xmin = 0, xmax = 5, ymin = y_badges - 0.12, ymax = y_badges + 0.12, fill = compl_col, color = NA) +
        ggplot2::annotate("text", x = 2.5, y = y_badges, label = toupper(as.character(s$completeness_proxy)), size = 2.8, color = "white", fontface = "bold") +
        ggplot2::annotate("text", x = -1, y = y_badges, label = "Completeness", hjust = 1, size = 3.0, color = "#374151") +
        ggplot2::annotate("rect", xmin = 6, xmax = 10, ymin = y_badges - 0.12, ymax = y_badges + 0.12, fill = conf_col, color = NA) +
        ggplot2::annotate("text", x = 8, y = y_badges, label = paste0("Confidence: ", toupper(as.character(s$confidence_tier))), size = 2.8, color = "white", fontface = "bold") +
        ggplot2::annotate("rect", xmin = 11, xmax = 14, ymin = y_badges - 0.12, ymax = y_badges + 0.12, fill = att_col, color = NA) +
        ggplot2::annotate("text", x = 12.5, y = y_badges, label = paste0("att: ", if (isTRUE(s$att_site_detected)) "YES" else "NO"), size = 2.8, color = if (isTRUE(s$att_site_detected)) "white" else "#6B7280", fontface = "bold") +
        # Metrics
        ggplot2::annotate("text", x = 16, y = y_badges, label = paste0(len_label, " | ", total_genes, " genes | phage ", pgf_label), hjust = 0, size = 3.2, color = "#374151") +
        ggplot2::annotate("text", x = 16, y = y_bar, label = if (is.finite(gc_delta)) paste0("GC: ", sprintf("%.1f%%", gc_region * 100), " (host ", sprintf("%.1f%%", gc_host * 100), ", dGC=", sprintf("%+.1f%%", gc_delta * 100), ")") else "", hjust = 0, size = 3.2, color = "#374151")

      # Integration site row
      y_int <- 0.08
      int_family <- if (!is.na(s$integrase_family %||% NA)) s$integrase_family else "none"
      int_pos <- if (!is.na(s$integrase_position %||% NA)) s$integrase_position else ""
      trna_info <- if (!is.na(s$trna_target %||% NA)) {
        # Extract just the first tRNA name and distance
        first_trna <- sub(";.*$", "", s$trna_target)
        first_trna <- sub("\\(pos:.*?\\)", "", first_trna)
        trimws(first_trna)
      } else {
        "no tRNA nearby"
      }
      right_flank <- if (!is.na(s$right_flank_genes %||% NA)) {
        paste0("downstream: ", s$right_flank_genes)
      } else {
        ""
      }
      att_core <- if (!is.na(s$att_core_seq %||% NA) && nzchar(s$att_core_seq %||% "")) {
        paste0(" | att core: ", s$att_core_seq, " (", s$att_core_length, "bp, ", s$att_confidence, ")")
      } else {
        ""
      }
      int_label <- paste0(
        "Integration: ", int_family,
        " (", int_pos, ")",
        " | target: ", trna_info,
        att_core,
        if (nzchar(right_flank)) paste0(" | ", right_flank) else ""
      )
      p <- p +
        ggplot2::annotate("text", x = -1, y = y_int, label = int_label, hjust = 1, size = 2.8, color = "#4B5563", fontface = "italic") +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_cartesian(xlim = c(-8, total_genes + 10), ylim = c(-0.06, 0.95), clip = "off") +
        ggplot2::labs(title = paste0("Prophage ", s$region_id, " Summary")) +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 12),
          plot.margin = ggplot2::margin(t = 6, r = 180, b = 4, l = 12)
        )
      p
    })
    if (length(summary_panels) == 1L) {
      p_summary <- summary_panels[[1]]
    } else {
      p_summary <- cowplot::plot_grid(plotlist = summary_panels, ncol = 1)
    }
  }
  # Panel D: Integration site comparative validation
  p_att_validation <- NULL
  if (is.data.frame(att_summary) && nrow(att_summary)) {
    # Load genome sequence for alignment-based ribbons
    gbff_for_d <- .dnmb_find_gbff_for_plot(output_dir)
    d_genome_seq <- ""
    if (!is.null(gbff_for_d)) {
      d_records <- .dnmb_prophage_parse_genbank_records(gbff_for_d)
      if (nrow(d_records)) d_genome_seq <- d_records$sequence[[1]] %||% ""
    }

    att_panels <- lapply(seq_len(nrow(att_summary)), function(ri) {
      s <- att_summary[ri, ]
      has_att <- !is.na(s$att_core_seq %||% NA) && nzchar(s$att_core_seq %||% "")
      has_trna <- !is.na(s$trna_target %||% NA) && nzchar(s$trna_target %||% "")
      has_integrase <- !is.na(s$integrase_family %||% NA) && nzchar(s$integrase_family %||% "")
      if (!has_trna) return(NULL)

      # Auto-detect reference genome and fetch locus annotations
      ref_genes <- NULL
      ref_seq <- ""
      ref_acc <- NULL
      ref_start <- NULL
      tryCatch({
        ref_acc_raw <- s$att_related_genome %||% NA_character_
        cache_d <- .dnmb_prophage_reference_cache_dir(output_dir, s$region_id)

        # If no reference from att detection, try auto-screening with cached/local genomes
        if (is.na(ref_acc_raw) || !nzchar(ref_acc_raw)) {
          # Look for cached reference genomes in the visualization cache
          cached_refs <- list.files(cache_d, pattern = "^ref_.*\\.fasta$", full.names = TRUE)
          parent_cache <- dirname(cache_d)
          if (!length(cached_refs)) {
            cached_refs <- list.files(parent_cache, pattern = "^ref_[A-Z][A-Z0-9]+\\.[0-9]+\\.fasta$", full.names = TRUE, recursive = TRUE)
          }
          if (length(cached_refs)) {
            auto_ref <- tryCatch(
              .dnmb_prophage_auto_find_reference(d_genome_seq, s$start, s$end, candidate_fastas = cached_refs),
              error = function(e) NULL
            )
            if (!is.null(auto_ref)) ref_acc_raw <- auto_ref$accession
          }
        }

        if (!is.na(ref_acc_raw) && nzchar(ref_acc_raw)) {
          ref_acc <- ref_acc_raw
          # Set accession on s for determine_ref_locus
          s_for_ref <- s
          s_for_ref$att_related_genome <- ref_acc_raw

          # Try to determine reference locus coordinates via BLAST
          ref_locus_info <- tryCatch(
            .dnmb_prophage_determine_ref_locus(s_for_ref, d_genome_seq, cache_dir = cache_d),
            error = function(e) NULL
          )

          if (!is.null(ref_locus_info)) {
            # Fetch GenBank subregion and parse gene annotations
            ref_data <- .dnmb_prophage_fetch_reference_locus(
              accession = ref_locus_info$accession,
              locus_start = ref_locus_info$locus_start,
              locus_end = ref_locus_info$locus_end,
              cache_dir = cache_d
            )
            if (!is.null(ref_data)) {
              ref_genes <- ref_data$genes_df
              ref_seq <- ref_data$sequence
              ref_start <- ref_data$locus_start
              ref_acc <- ref_data$accession
            }
          }
        }
      }, error = function(e) NULL)

      tryCatch(
        .dnmb_plot_att_panel_gggenes(
          s, genbank_table,
          genome_seq = d_genome_seq,
          ref_seq_full = ref_seq,
          ref_start_bp = ref_start,
          ref_genes_df = ref_genes,
          ref_accession = ref_acc,
          cache_d = cache_d
        ),
        error = function(e) NULL
      )
    })
    att_panels <- Filter(Negate(is.null), att_panels)
    if (length(att_panels)) {
      p_att_validation <- if (length(att_panels) == 1L) att_panels[[1]] else cowplot::plot_grid(plotlist = att_panels, ncol = 1)
    }
  }

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "Prophage_overview.pdf")
  category_legend_row <- .dnmb_inline_discrete_legend_plot(
    title = "Category",
    colors = .dnmb_gene_arrow_palette_prophage()
  )
  if (!is.null(category_legend_row)) {
    plots_to_save <- c(region_panels, list(category_legend_row, p_score), if (!is.null(p_reference)) list(p_reference) else list(), if (!is.null(p_att_validation)) list(p_att_validation) else list())
    legend_h <- if (length(.dnmb_gene_arrow_palette_prophage()) > 5L) 0.22 else 0.12
    next_letter <- length(region_panels) + 1L
    ref_offset <- if (!is.null(p_reference)) 1L else 0L
    rel_heights <- c(rep(1.12, length(region_panels)), legend_h, 1.0, if (!is.null(p_reference)) 1.45 else NULL, if (!is.null(p_att_validation)) 0.85 else NULL)
    labels_to_use <- c(
      base::LETTERS[seq_along(region_panels)],
      "",
      base::LETTERS[[next_letter]],
      if (!is.null(p_reference)) base::LETTERS[[next_letter + 1L]] else NULL,
      if (!is.null(p_att_validation)) base::LETTERS[[next_letter + ref_offset + 1L]] else NULL
    )
  } else {
    plots_to_save <- c(region_panels, list(p_score), if (!is.null(p_reference)) list(p_reference) else list(), if (!is.null(p_att_validation)) list(p_att_validation) else list())
    next_letter <- length(region_panels) + 1L
    ref_offset <- if (!is.null(p_reference)) 1L else 0L
    rel_heights <- c(rep(1.12, length(region_panels)), 1.0, if (!is.null(p_reference)) 1.45 else NULL, if (!is.null(p_att_validation)) 0.85 else NULL)
    labels_to_use <- c(
      base::LETTERS[seq_along(region_panels)],
      base::LETTERS[[next_letter]],
      if (!is.null(p_reference)) base::LETTERS[[next_letter + 1L]] else NULL,
      if (!is.null(p_att_validation)) base::LETTERS[[next_letter + ref_offset + 1L]] else NULL
    )
  }
  composite <- cowplot::plot_grid(
    plotlist = plots_to_save,
    labels = labels_to_use,
    ncol = 1,
    rel_heights = rel_heights,
    align = "v",
    axis = "lr"
  )
  total_panels <- length(region_panels) + 1L + (if (!is.null(p_reference)) 1L else 0L) + (if (!is.null(p_summary)) 1L else 0L)
  ggplot2::ggsave(pdf_path, composite, width = 14, height = 3.0 * length(region_panels) + total_panels * 0.6 + 2.5, bg = "white")
  list(pdf = pdf_path)
}

.dnmb_prophage_fetch_reference_locus <- function(accession, locus_start, locus_end, cache_dir = NULL) {
  if (is.null(accession) || is.na(accession) || !nzchar(accession)) return(NULL)
  locus_start <- as.integer(locus_start)
  locus_end <- as.integer(locus_end)
  if (is.na(locus_start) || is.na(locus_end) || locus_start >= locus_end) return(NULL)

  # Expand flanking region to capture context genes around att site

  flank_expand <- 8000L
  fetch_start <- max(1L, locus_start - flank_expand)
  fetch_end <- locus_end + flank_expand

  dest_dir <- cache_dir %||% tempdir()
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest_gbk <- file.path(dest_dir, paste0("ref_locus_", accession, "_", fetch_start, "_", fetch_end, ".gbk"))

  if (!file.exists(dest_gbk) || file.info(dest_gbk)$size == 0) {
    url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",
      accession,
      "&rettype=gbwithparts&retmode=text",
      "&seq_start=", fetch_start,
      "&seq_stop=", fetch_end
    )
    run <- dnmb_run_external("curl", args = c("-L", "-k", "-s", url, "-o", dest_gbk), required = FALSE)
    if (!isTRUE(run$ok) || !file.exists(dest_gbk) || file.info(dest_gbk)$size == 0) {
      return(NULL)
    }
  }

  lines <- readLines(dest_gbk, warn = FALSE)
  if (!length(lines) || !any(grepl("^LOCUS", lines))) return(NULL)

  # Parse FEATURES section for CDS and tRNA entries
  feat_idx <- grep("^FEATURES", lines)
  origin_idx <- grep("^ORIGIN", lines)
  end_idx <- grep("^//", lines)
  if (!length(feat_idx) || !length(origin_idx)) return(NULL)

  # Extract sequence
  sequence <- ""
  if (length(origin_idx) && length(end_idx) && end_idx[[1]] > origin_idx[[1]]) {
    sequence <- paste(lines[(origin_idx[[1]] + 1L):(end_idx[[1]] - 1L)], collapse = "")
    sequence <- gsub("[^ACGTacgt]", "", sequence)
    sequence <- toupper(sequence)
  }

  feat_lines <- lines[(feat_idx[[1]] + 1L):(origin_idx[[1]] - 1L)]
  entries <- list()
  i <- 1L
  while (i <= length(feat_lines)) {
    line <- feat_lines[[i]]
    is_cds <- grepl("^     CDS\\s+", line)
    is_trna <- grepl("^     tRNA\\s+", line)
    if (is_cds || is_trna) {
      feat_type <- if (is_cds) "CDS" else "tRNA"
      location <- trimws(sub("^     (CDS|tRNA)\\s+", "", line))
      j <- i + 1L
      qualifiers <- character()
      while (j <= length(feat_lines) && grepl("^                     ", feat_lines[[j]])) {
        qualifiers <- c(qualifiers, trimws(feat_lines[[j]]))
        j <- j + 1L
      }
      starts <- suppressWarnings(as.numeric(unlist(regmatches(location, gregexpr("[0-9]+", location)))))
      if (length(starts)) {
        gene <- NA_character_
        product <- NA_character_
        if (length(qualifiers)) {
          gene_match <- qualifiers[grepl('^/gene="', qualifiers)]
          if (length(gene_match)) gene <- sub('^/gene="([^"]+)".*$', '\\1', gene_match[[1]])
          prod_match <- qualifiers[grepl('^/product="', qualifiers)]
          if (length(prod_match)) product <- sub('^/product="([^"]+)".*$', '\\1', prod_match[[1]])
        }
        label <- if (!is.na(gene) && nzchar(gene)) gene else if (!is.na(product) && nzchar(product)) {
          sub(" family.*$", "", product)
        } else "hyp."
        if (feat_type == "tRNA" && !is.na(product)) {
          # Extract amino acid abbreviation from tRNA product (e.g., "tRNA-Leu" -> "Leu")
          aa <- sub("^tRNA-", "", product)
          if (nzchar(aa)) label <- aa
        }
        is_complement <- grepl("complement", location)
        entries[[length(entries) + 1L]] <- data.frame(
          gene = label,
          start = as.integer(min(starts)),
          end = as.integer(max(starts)),
          strand = if (is_complement) "reverse" else "forward",
          orientation = if (is_complement) -1L else 1L,
          fill = if (feat_type == "tRNA") "#34D399" else "#D1D5DB",
          type = if (feat_type == "tRNA") "tRNA" else "host",
          stringsAsFactors = FALSE
        )
      }
      i <- j
      next
    }
    i <- i + 1L
  }

  if (!length(entries)) return(NULL)

  ref_df <- dplyr::bind_rows(entries)
  ref_df$molecule <- paste0("Ref: ", accession)

  # Return parsed locus data
  list(
    genes_df = ref_df,
    sequence = sequence,
    accession = accession,
    locus_start = fetch_start,
    locus_end = fetch_end
  )
}

.dnmb_prophage_determine_ref_locus <- function(att_row, genome_seq, cache_dir = NULL) {
  # Use local BLAST to find where flanking regions map on the reference genome
  ref_acc <- att_row$att_related_genome %||% NA_character_
  if (is.na(ref_acc) || !nzchar(ref_acc)) return(NULL)

  pp_s <- as.integer(att_row$start)
  pp_e <- as.integer(att_row$end)
  if (is.na(pp_s) || is.na(pp_e) || !nzchar(genome_seq)) return(NULL)

  tmp_dir <- base::tempdir()  # Use tempdir (no spaces) for BLAST

  # Find cached reference FASTA
  ref_fasta <- NULL
  search_dirs <- c(cache_dir, if (!is.null(cache_dir)) dirname(cache_dir) else NULL)
  for (sd in search_dirs) {
    if (is.null(sd) || !dir.exists(sd)) next
    candidates <- list.files(sd, pattern = paste0(".*", sub("\\..*$", "", ref_acc), ".*\\.fasta$"), full.names = TRUE, recursive = TRUE)
    if (length(candidates)) { ref_fasta <- candidates[1]; break }
  }
  if (is.null(ref_fasta) || !file.exists(ref_fasta)) return(NULL)

  # Copy to safe path and make BLAST db
  safe_ref <- file.path(tmp_dir, paste0("ref_locus_", basename(ref_fasta)))
  if (!file.exists(safe_ref)) base::file.copy(ref_fasta, safe_ref, overwrite = TRUE)
  db_prefix <- file.path(tmp_dir, paste0("ref_locus_db_", ref_acc))
  mk <- dnmb_run_external("makeblastdb", args = c("-in", safe_ref, "-dbtype", "nucl", "-out", db_prefix), required = FALSE)
  if (!isTRUE(mk$ok)) return(NULL)

  # Extract flanking sequences
  flank_len <- 1000L
  left_seq  <- substr(genome_seq, max(1L, pp_s - flank_len), pp_s - 1L)
  right_seq <- substr(genome_seq, pp_e + 1L, min(nchar(genome_seq), pp_e + flank_len))
  if (nchar(right_seq) < 100L) return(NULL)

  # Helper: run a single-flank BLAST and parse outfmt 6
  .blast_flank <- function(label, seq_str) {
    if (nchar(seq_str) < 100L) return(NULL)
    qf <- file.path(tmp_dir, paste0("ref_locus_", label, "_query.fasta"))
    of <- file.path(tmp_dir, paste0("ref_locus_", label, "_", ref_acc, ".tsv"))
    .dnmb_write_single_fasta(label, seq_str, qf)
    run <- dnmb_run_external("blastn", args = c(
      "-task", "megablast", "-db", db_prefix, "-query", qf,
      "-evalue", "1e-10", "-max_target_seqs", "5", "-outfmt", "6", "-out", of
    ), required = FALSE)
    if (!isTRUE(run$ok) || !file.exists(of) || file.info(of)$size == 0) return(NULL)
    tbl <- tryCatch(
      utils::read.table(of, sep = "\t", header = FALSE, quote = "",
                        comment.char = "", stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(tbl) || !nrow(tbl) || ncol(tbl) < 12) return(NULL)
    tbl$sstart <- suppressWarnings(as.integer(tbl$V9))
    tbl$send   <- suppressWarnings(as.integer(tbl$V10))
    tbl$length <- suppressWarnings(as.integer(tbl$V4))
    tbl$s_min  <- pmin(tbl$sstart, tbl$send)
    tbl$s_max  <- pmax(tbl$sstart, tbl$send)
    tbl
  }

  # BLAST right flank separately — best/longest hit = correct tRNA locus
  right_tbl <- .blast_flank("right", right_seq)
  if (is.null(right_tbl)) return(NULL)
  best_right <- right_tbl[which.max(right_tbl$length), ]
  right_center <- (best_right$s_min + best_right$s_max) / 2

  # BLAST left flank separately
  left_tbl <- .blast_flank("left", left_seq)

  # Check if any left hit falls within 50 kb of the right flank hit
  left_near <- NULL
  if (!is.null(left_tbl) && nrow(left_tbl)) {
    near_idx <- which(abs((left_tbl$s_min + left_tbl$s_max) / 2 - right_center) < 50000)
    if (length(near_idx)) {
      left_near <- left_tbl[near_idx, , drop = FALSE]
      left_near <- left_near[which.max(left_near$length), ]
    }
  }

  # Determine locus range
  if (!is.null(left_near)) {
    # Both flanks hit nearby — narrow window
    ref_locus_start <- max(1L, as.integer(left_near$s_min) - 2000L)
    ref_locus_end   <- as.integer(best_right$s_max) + 2000L
  } else {
    # Right flank only — wider upstream window
    ref_locus_start <- max(1L, as.integer(best_right$s_min) - 7000L)
    ref_locus_end   <- as.integer(best_right$s_max) + 2000L
  }

  list(accession = ref_acc, locus_start = ref_locus_start, locus_end = ref_locus_end)
}

.dnmb_plot_att_panel_gggenes <- function(s, genbank_table, genome_seq = "", ref_seq_full = "",
                                         ref_start_bp = NULL, ref_genes_df = NULL, ref_accession = NULL, cache_d = NULL) {
  has_att <- !is.na(s$att_core_seq %||% NA) && nzchar(s$att_core_seq %||% "")
  has_integrase <- !is.na(s$integrase_family %||% NA) && nzchar(s$integrase_family %||% "")
  pp_s <- s$start
  pp_e <- s$end
  flank_bp <- 5000

  # Query genes around integration site
  gt <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  gt$gs <- suppressWarnings(as.numeric(gt$start))
  gt$ge <- suppressWarnings(as.numeric(gt$end))
  q_view <- gt[gt$contig == s$contig & gt$ge >= (pp_s - flank_bp) & gt$gs <= (pp_e + flank_bp), , drop = FALSE]
  q_view$is_pp <- !is.na(q_view$Prophage_prophage_id) & nzchar(q_view$Prophage_prophage_id %||% "") & q_view$Prophage_prophage_id != "0"
  q_view$is_trna <- grepl("tRNA", q_view$product)
  q_view$category <- .dnmb_gene_arrow_category_prophage(q_view$product)
  palette <- .dnmb_gene_arrow_palette_prophage()
  q_view$fill_col <- ifelse(q_view$is_trna, "#34D399",
                    ifelse(q_view$is_pp, unname(palette[q_view$category]), "#D1D5DB"))
  q_view$fill_col[is.na(q_view$fill_col)] <- "#D1D5DB"
  q_view$forward <- q_view$direction == "forward"

  # Build gggenes data for query
  q_genes_df <- data.frame(
    molecule = paste0("Query: ", s$contig),
    gene = q_view$locus_tag,
    start = as.integer(q_view$gs),
    end = as.integer(q_view$ge),
    strand = ifelse(q_view$forward, "forward", "reverse"),
    orientation = ifelse(q_view$forward, 1L, -1L),
    fill = q_view$fill_col,
    type = ifelse(q_view$is_trna, "tRNA", ifelse(q_view$is_pp, "prophage", "host")),
    stringsAsFactors = FALSE
  )

  # Reference gene annotations: use provided ref_genes_df or fall back to query-only view
  if (is.null(ref_genes_df) || !is.data.frame(ref_genes_df) || !nrow(ref_genes_df)) {
    # No reference data available — return query-only att panel (no comparative view)
    return(NULL)
  }
  # Ensure ref_genes_df has all required columns
  required_cols <- c("gene", "start", "end", "strand", "orientation", "fill", "type")
  if (!all(required_cols %in% names(ref_genes_df))) return(NULL)
  # Get organism name from GenBank SOURCE field (most accurate)
  ref_organism_name <- ref_accession %||% "reference"
  tryCatch({
    # First try cached GenBank (.gbk) — has SOURCE field
    if (!is.null(cache_d)) {
      gbk_files <- list.files(cache_d, pattern=paste0(".*", sub("\\..*$", "", ref_accession %||% ""), ".*\\.gbk$"), full.names=TRUE, recursive=TRUE)
      if (length(gbk_files)) {
        gbk_lines <- readLines(gbk_files[1], n=15, warn=FALSE)
        src_line <- gbk_lines[grepl("^SOURCE", gbk_lines)]
        if (length(src_line)) {
          org <- trimws(sub("^SOURCE\\s+", "", src_line[1]))
          if (nzchar(org)) ref_organism_name <- org
        }
      }
    }
    # Fallback to FASTA header
    if (ref_organism_name == (ref_accession %||% "reference")) {
      parent_cache <- if (!is.null(cache_d)) dirname(cache_d) else NULL
      if (!is.null(parent_cache)) {
        ref_files <- list.files(parent_cache, pattern=paste0(".*", sub("\\..*$", "", ref_accession %||% ""), ".*\\.fasta$"), full.names=TRUE, recursive=TRUE)
        if (length(ref_files)) {
          hdr <- readLines(ref_files[1], n=1, warn=FALSE)
          org <- sub("^>[^ ]+ ", "", hdr)
          org <- sub(", complete.*$", "", org)
          org <- sub(" chromosome.*$", "", org)
          if (nzchar(org)) ref_organism_name <- org
        }
      }
    }
  }, error=function(e) NULL)

  if (!"molecule" %in% names(ref_genes_df)) {
    ref_genes_df$molecule <- ref_organism_name
  }

  att_seq <- if (has_att) s$att_core_seq else "?"
  att_len <- if (has_att) s$att_core_length else "?"
  compl <- toupper(s$completeness_proxy %||% "?")
  conf <- toupper(s$confidence_tier %||% "?")
  gc_info <- if (is.finite(as.numeric(s$delta_gc %||% NA))) sprintf("dGC=%+.1f%%", as.numeric(s$delta_gc)*100) else ""
  int_info <- if (has_integrase) paste0(s$integrase_family, " [", s$integrase_position, "]") else ""
  subtitle_text <- paste0(compl, " | ", conf,
    if (nzchar(gc_info)) paste0(" | ", gc_info) else "",
    " | ", int_info,
    if (has_att) paste0(" | att core: 5'-", att_seq, "-3' (", att_len, "bp)") else "")

  # Normalize both to 0-1 x coordinates for single-plot comparison with ribbons
  q_min <- min(q_genes_df$start); q_max <- max(q_genes_df$end); q_range <- q_max - q_min
  r_min <- min(ref_genes_df$start); r_max <- max(ref_genes_df$end); r_range <- r_max - r_min
  q_genes_df$xn_s <- (q_genes_df$start - q_min) / q_range
  q_genes_df$xn_e <- (q_genes_df$end - q_min) / q_range
  ref_genes_df$xn_s <- (ref_genes_df$start - r_min) / r_range
  ref_genes_df$xn_e <- (ref_genes_df$end - r_min) / r_range

  y_q <- 2.0; y_r <- 1.0; gh <- 0.28
  q_genes_df$y <- y_q; ref_genes_df$y <- y_r

  # Build synteny ribbons connecting homologous flanking regions
  ribbon_list_att <- list()
  ri_idx <- 1L
  # Connect regions based on Biostrings local alignment
  # Align query flanking regions against reference locus to find actual homology
  q_right_genes <- q_genes_df[q_genes_df$xn_s > (pp_e - q_min) / q_range, , drop=FALSE]
  q_left_genes <- q_genes_df[q_genes_df$xn_e < (pp_s - q_min) / q_range, , drop=FALSE]
  r_all <- ref_genes_df

  # Use Biostrings pairwiseAlignment for flanking homology detection against reference locus
  # ref_seq_full is the reference locus sequence (not the whole genome) when fetched via efetch subregion
  ref_locus_start_bp <- ref_start_bp %||% 1L
  tryCatch({
    q_right_seq <- substr(genome_seq, pp_e + 1, min(nchar(genome_seq), pp_e + 3000))
    # Use ref_seq_full directly as the reference locus sequence
    # Use ref GenBank sequence, or cached FASTA, for alignment
    ref_locus_seq <- if (nchar(ref_seq_full) > 0) ref_seq_full else ""
    if (!nzchar(ref_locus_seq)) {
      # Try to read from cached reference FASTA
      tryCatch({
        parent_c <- if (!is.null(cache_d)) dirname(cache_d) else NULL
        if (!is.null(parent_c)) {
          rf <- list.files(parent_c, pattern=paste0(".*", sub("\\..*$", "", ref_accession %||% ""), ".*\\.fasta$"), full.names=TRUE, recursive=TRUE)
          if (length(rf)) {
            rlines <- readLines(rf[1], warn=FALSE)
            ref_locus_seq <- toupper(paste(rlines[!grepl("^>", rlines)], collapse=""))
          }
        }
      }, error=function(e) NULL)
    }
    if (!nzchar(ref_locus_seq) || !nzchar(q_right_seq)) stop("no sequences for alignment")
    q_right_dna <- Biostrings::DNAString(q_right_seq)
    ref_dna <- Biostrings::DNAString(ref_locus_seq)

    # Use pairwiseAlignment for each flank to get identity and exact positions
    identity_labels <- list()  # Store for annotation

    # RIGHT FLANKING alignment (host region AFTER prophage, NOT the prophage itself)
    r_aln <- Biostrings::pairwiseAlignment(q_right_dna, ref_dna, type="local")
    r_score <- Biostrings::score(r_aln)
    if (r_score > 50) {
      r_pid <- round(Biostrings::pid(r_aln), 1)
      r_q_s <- Biostrings::start(Biostrings::pattern(r_aln))
      r_q_e <- Biostrings::end(Biostrings::pattern(r_aln))
      r_r_s <- Biostrings::start(Biostrings::subject(r_aln))
      r_r_e <- Biostrings::end(Biostrings::subject(r_aln))
      r_aln_len <- Biostrings::nchar(r_aln)
      q_real_s <- pp_e + r_q_s
      q_real_e <- pp_e + r_q_e
      q_n_s <- (q_real_s - q_min) / q_range
      q_n_e <- (q_real_e - q_min) / q_range
      r_n_s <- r_r_s / nchar(ref_locus_seq)
      r_n_e <- r_r_e / nchar(ref_locus_seq)
      ribbon_list_att[[ri_idx]] <- data.frame(
        x = c(q_n_s, q_n_e, r_n_e, r_n_s),
        y = c(y_q - gh, y_q - gh, y_r + gh, y_r + gh), group=ri_idx, fill="#60A5FA", stringsAsFactors=FALSE)
      ri_idx <- ri_idx + 1L
      identity_labels[[length(identity_labels)+1]] <- list(
        x=(q_n_e+r_n_e)/2, y=1.5, side="right",
        label=paste0("Right flank\n", r_pid, "% id, ", r_aln_len, "bp")
      )
    }

    # LEFT FLANKING alignment (host region BEFORE prophage)
    q_left_seq <- substr(genome_seq, max(1, pp_s - 3000), pp_s - 1)
    q_left_dna <- Biostrings::DNAString(q_left_seq)
    l_aln <- Biostrings::pairwiseAlignment(q_left_dna, ref_dna, type="local")
    l_score <- Biostrings::score(l_aln)
    if (l_score > 50) {
      l_pid <- round(Biostrings::pid(l_aln), 1)
      l_q_s <- Biostrings::start(Biostrings::pattern(l_aln))
      l_q_e <- Biostrings::end(Biostrings::pattern(l_aln))
      l_r_s <- Biostrings::start(Biostrings::subject(l_aln))
      l_r_e <- Biostrings::end(Biostrings::subject(l_aln))
      l_aln_len <- Biostrings::nchar(l_aln)
      q_left_start <- max(1, pp_s - 3000)
      q_real_s <- q_left_start + l_q_s - 1
      q_real_e <- q_left_start + l_q_e - 1
      q_n_s <- (q_real_s - q_min) / q_range
      q_n_e <- (q_real_e - q_min) / q_range
      r_n_s <- l_r_s / nchar(ref_locus_seq)
      r_n_e <- l_r_e / nchar(ref_locus_seq)
      ribbon_list_att[[ri_idx]] <- data.frame(
        x = c(q_n_s, q_n_e, r_n_e, r_n_s),
        y = c(y_q - gh, y_q - gh, y_r + gh, y_r + gh), group=ri_idx, fill="#60A5FA", stringsAsFactors=FALSE)
      ri_idx <- ri_idx + 1L
      identity_labels[[length(identity_labels)+1]] <- list(
        x=(q_n_s+r_n_s)/2, y=1.5, side="left",
        label=paste0("Left flank\n", l_pid, "% id, ", l_aln_len, "bp")
      )
    } else {
      # No left flank homology — annotate explicitly
      identity_labels[[length(identity_labels)+1]] <- list(
        x=0.05, y=1.5, side="left",
        label="Left flank:\nno homology"
      )
    }
  }, error = function(e) NULL)
  ribbon_att <- if (length(ribbon_list_att)) dplyr::bind_rows(ribbon_list_att) else data.frame()

  # Prophage region in normalized coords
  pp_xn_s <- (pp_s - q_min) / q_range; pp_xn_e <- (pp_e - q_min) / q_range
  # attB in ref normalized — estimate from tRNA cluster midpoint in reference
  ref_trna_rows <- ref_genes_df[ref_genes_df$type == "tRNA", , drop = FALSE]
  attB_xn <- if (nrow(ref_trna_rows)) {
    trna_mid <- (min(ref_trna_rows$start) + max(ref_trna_rows$end)) / 2
    (trna_mid - r_min) / r_range
  } else {
    0.5  # fallback to center if no tRNA in reference
  }

  # Scale bar labels
  q_label <- paste0(.dnmb_fmt_bp_label(q_range), " region")
  r_label <- paste0(.dnmb_fmt_bp_label(r_range), " region")

  p <- ggplot2::ggplot() +
    # Ribbons first (behind everything)
    { if (nrow(ribbon_att)) ggplot2::geom_polygon(data=ribbon_att, ggplot2::aes(x=.data$x, y=.data$y, group=.data$group, fill=.data$fill), alpha=0.18, color=NA, inherit.aes=FALSE) else NULL } +
    # Backbone lines
    ggplot2::annotate("segment", x=0, xend=1, y=y_q, yend=y_q, linewidth=0.4, color="#9CA3AF") +
    ggplot2::annotate("segment", x=0, xend=1, y=y_r, yend=y_r, linewidth=0.4, color="#9CA3AF") +
    # Prophage shading on query
    ggplot2::annotate("rect", xmin=pp_xn_s, xmax=pp_xn_e, ymin=y_q-gh-0.02, ymax=y_q+gh+0.02, fill="#DBEAFE", alpha=0.25, color="#3B82F6", linewidth=0.3) +
    ggplot2::annotate("text", x=(pp_xn_s+pp_xn_e)/2, y=y_q+gh+0.12, label=paste0("Prophage (", .dnmb_fmt_bp_label(pp_e-pp_s+1), ")"), color="#1E40AF", fontface="bold", size=3.2) +
    # Query gene arrows
    gggenes::geom_gene_arrow(data=q_genes_df, ggplot2::aes(xmin=.data$xn_s, xmax=.data$xn_e, y=y_q, fill=.data$fill, forward=.data$orientation>0),
                             arrowhead_height=grid::unit(4,"mm"), arrowhead_width=grid::unit(1.5,"mm"), arrow_body_height=grid::unit(3.5,"mm")) +
    # Reference gene arrows
    gggenes::geom_gene_arrow(data=ref_genes_df, ggplot2::aes(xmin=.data$xn_s, xmax=.data$xn_e, y=y_r, fill=.data$fill, forward=.data$orientation>0),
                             arrowhead_height=grid::unit(4,"mm"), arrowhead_width=grid::unit(1.5,"mm"), arrow_body_height=grid::unit(3.5,"mm")) +
    # Reference gene labels (CDS only, skip tRNA - too small)
    {
      ref_cds <- ref_genes_df[ref_genes_df$type != "tRNA" & (ref_genes_df$xn_e - ref_genes_df$xn_s) > 0.03, , drop=FALSE]
      if (nrow(ref_cds)) {
        ref_cds$label_short <- sub(" family.*$", "", sub("hypothetical protein", "hyp.", ref_cds$gene))
        ref_cds$label_short <- substr(ref_cds$label_short, 1, 20)
        ggplot2::geom_text(data=ref_cds, ggplot2::aes(x=(.data$xn_s+.data$xn_e)/2, y=y_r-gh-0.06, label=.data$label_short),
                          size=1.8, color="#374151", angle=30, hjust=1, vjust=1, inherit.aes=FALSE)
      } else NULL
    } +
    # tRNA label (grouped)
    {
      ref_trna <- ref_genes_df[ref_genes_df$type == "tRNA", , drop=FALSE]
      if (nrow(ref_trna)) {
        ggplot2::annotate("text", x=median(c(ref_trna$xn_s, ref_trna$xn_e)), y=y_r-gh-0.06,
                         label=paste0("tRNA (", nrow(ref_trna), ")"), size=2.0, color="#059669", fontface="italic")
      } else NULL
    } +
    # att markers
    ggplot2::annotate("segment", x=pp_xn_s, xend=pp_xn_s, y=y_q-gh-0.05, yend=y_q+gh+0.05, linetype="dashed", color="#DC2626", linewidth=0.5) +
    ggplot2::annotate("segment", x=pp_xn_e, xend=pp_xn_e, y=y_q-gh-0.05, yend=y_q+gh+0.05, linetype="dashed", color="#DC2626", linewidth=0.5) +
    ggplot2::annotate("text", x=pp_xn_s, y=y_q-gh-0.12, label="attL", color="#DC2626", fontface="bold", size=2.8) +
    ggplot2::annotate("text", x=pp_xn_e, y=y_q-gh-0.12, label="attR", color="#DC2626", fontface="bold", size=2.8) +
    ggplot2::annotate("segment", x=attB_xn, xend=attB_xn, y=y_r-gh-0.05, yend=y_r+gh+0.05, linetype="dashed", color="#22C55E", linewidth=0.5) +
    ggplot2::annotate("text", x=attB_xn, y=y_r+gh+0.12, label="attB", color="#22C55E", fontface="bold", size=2.8) +
    # Genome labels with real coordinates and organism name
    ggplot2::annotate("text", x=-0.02, y=y_q, label=paste0(sub(",.*$", "", s$contig), "\n", format(q_min, big.mark=","), "-", format(q_max, big.mark=",")),
                     hjust=1, size=2.4, fontface="bold", lineheight=0.85) +
    ggplot2::annotate("text", x=-0.02, y=y_r,
                     label=paste0(ref_organism_name,
                                  "\n", ref_accession %||% "", " | site ",
                                  format(ref_locus_start_bp %||% 0, big.mark=","), "-",
                                  format((ref_locus_start_bp %||% 0) + as.integer(max(ref_genes_df$end)), big.mark=",")),
                     hjust=1, size=2.2, fontface="bold", lineheight=0.85, color="#059669") +
    # att core — right side annotation, outside plot area
    { if (has_att) ggplot2::annotate("label", x=1.03, y=(y_q+y_r)/2,
                                     label=paste0("att core\n5'-", att_seq, "-3'\n(", att_len, "bp)\ntRNA-Leu hotspot"),
                                     hjust=0, fill="#FEF2F2", color="#DC2626", fontface="bold", size=2.5,
                                     lineheight=0.85, label.r=grid::unit(0.15,"lines"), label.padding=grid::unit(0.2,"lines")) else NULL } +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim=c(-0.15, 1.18), ylim=c(y_r-gh-0.35, y_q+gh+0.25), clip="off") +
    ggplot2::labs(title=paste0("Integration Site: Prophage ", s$region_id), subtitle=subtitle_text, y=NULL) +
    ggplot2::theme_void(base_size=11) +
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold", size=12),
                   plot.subtitle=ggplot2::element_text(size=8.5, color="#4B5563"),
                   plot.margin=ggplot2::margin(t=8,r=12,b=8,l=12))
  # Add identity labels on ribbons
  if (exists("identity_labels") && length(identity_labels)) {
    for (il in identity_labels) {
      p <- p + ggplot2::annotate("label", x=il$x, y=il$y - 0.18,
                                 label=il$label, fill="white", color="#1D4ED8",
                                 fontface="bold", size=2.3, label.r=grid::unit(0.1,"lines"),
                                 label.padding=grid::unit(0.12,"lines"))
    }
  }
  p
}

# ---------- PAZy Publication-Quality Multi-Panel Figure --------------------
#
# Panel A: Genome overview with zoomed gene loci (gggenes)
# Panel B: Bitscore bars with CLEAN EC classification (ggnewscale)
# Panel C: Catalytic site analysis - G-x-S-x-G motif heatmap +
#          alignment-validated catalytic triad (Biostrings)
# --------------------------------------------------------------------------

.dnmb_plot_pazy_pub <- function(genbank_table, output_dir, cache_root = NULL) {

  # --- check optional dependencies ---
  for (pkg in c("gggenes", "ggnewscale", "Biostrings", "colorspace")) {
    if (!base::requireNamespace(pkg, quietly = TRUE)) {
      base::message("[DNMB] PAZy publication figure requires package '", pkg, "'")
      return(NULL)
    }
  }

  plot_dir <- .dnmb_module_plot_dir(output_dir)

  # --- Load & filter PAZy data ---
  pazy_path <- base::file.path(output_dir, "dnmb_module_pazy", "pazy_merged.tsv")
  if (!base::file.exists(pazy_path)) {
    base::message("[DNMB] PAZy merged file not found: ", pazy_path)
    return(NULL)
  }
  d <- utils::read.delim(pazy_path, stringsAsFactors = FALSE) |>
    dplyr::mutate(start = base::as.integer(.data$start),
                  end   = base::as.integer(.data$end))

  genome_len <- base::max(d$end, na.rm = TRUE)

  sig <- d |>
    dplyr::filter(!base::is.na(.data$PAZy_family_id),
                  base::nzchar(.data$PAZy_family_id),
                  .data$PAZy_evalue < 0.01,
                  .data$PAZy_qcov >= 0.4,
                  .data$PAZy_scov >= 0.4) |>
    dplyr::mutate(
      midpoint  = (.data$start + .data$end) / 2,
      substrate = base::gsub("; ", "/", .data$PAZy_substrate_label),
      gene_info = base::ifelse(base::is.na(.data$gene) | !base::nzchar(.data$gene),
                               .data$locus_tag, .data$gene),
      hit_id    = base::paste0(.data$PAZy_family_id, "  (", .data$locus_tag, ")"),
      PAZy_qcov = base::pmin(.data$PAZy_qcov, 1),
      PAZy_scov = base::pmin(.data$PAZy_scov, 1)
    )

  if (!base::nrow(sig)) {
    base::message("[DNMB] No significant PAZy hits for publication figure")
    return(NULL)
  }

  # --- Get CLEAN data from genbank_table ---
  has_clean <- "CLEAN_best_hit_label" %in% base::names(genbank_table)
  if (has_clean) {
    clean_cols <- base::intersect(
      base::names(genbank_table),
      base::c("locus_tag", "CLEAN_best_hit_label", "CLEAN_pvalue_best_pvalue", "CLEAN_best_distance")
    )
    clean_lookup <- genbank_table[, clean_cols, drop = FALSE]
    sig <- dplyr::left_join(sig, clean_lookup, by = "locus_tag")
    sig$clean_ec <- base::as.character(sig$CLEAN_best_hit_label)
    sig$clean_pvalue <- if ("CLEAN_pvalue_best_pvalue" %in% base::names(sig)) suppressWarnings(base::as.numeric(sig$CLEAN_pvalue_best_pvalue)) else NA_real_
    sig$clean_distance <- if ("CLEAN_best_distance" %in% base::names(sig)) suppressWarnings(base::as.numeric(sig$CLEAN_best_distance)) else NA_real_
  } else {
    sig$clean_ec       <- NA_character_
    sig$clean_pvalue   <- NA_real_
    sig$clean_distance <- NA_real_
  }

  # --- EC class helpers ---
  assign_ec_group <- function(ec) {
    base::sapply(ec, function(e) {
      if (base::is.na(e) || !base::nzchar(e)) return("Other")
      parts <- base::as.numeric(base::unlist(base::strsplit(e, "\\.")))
      if (base::length(parts) < 3L) return("Other")
      if (parts[1] == 3 && parts[2] == 1) return("Esterase / Lipase")
      if (parts[1] == 3 && parts[2] == 4) return("Peptidase / Protease")
      "Other hydrolase"
    })
  }

  assign_ec_specific <- function(ec) {
    base::sapply(ec, function(e) {
      if (base::is.na(e) || !base::nzchar(e)) return("Unclassified")
      parts <- base::as.numeric(base::unlist(base::strsplit(e, "\\.")))
      if (base::length(parts) < 4L) return("Other")
      if (parts[1] == 3 && parts[2] == 1 && parts[3] == 1) {
        if (parts[4] == 101) return("PETase (3.1.1.101)")
        if (parts[4] == 74)  return("Cutinase (3.1.1.74)")
        if (parts[4] == 3)   return("Lipase (3.1.1.3)")
        if (parts[4] == 2)   return("Arylesterase (3.1.1.2)")
        if (parts[4] == 1)   return("Carboxylesterase (3.1.1.1)")
        if (parts[4] == 23)  return("Acylglycerol lipase (3.1.1.23)")
        if (parts[4] == 85)  return("MenH esterase (3.1.1.85)")
        return(base::paste0("Esterase (", e, ")"))
      }
      if (parts[1] == 3 && parts[2] == 1 && parts[3] == 2)
        return(base::paste0("Thioesterase (", e, ")"))
      if (parts[1] == 3 && parts[2] == 4)
        return(base::paste0("Peptidase (", e, ")"))
      if (parts[1] == 3 && parts[2] == 2)
        return(base::paste0("Glycosidase (", e, ")"))
      if (parts[1] == 6)
        return(base::paste0("Ligase (", e, ")"))
      base::paste0("Other (", e, ")")
    })
  }

  sig$ec_group    <- assign_ec_group(sig$clean_ec)
  sig$ec_specific <- assign_ec_specific(sig$clean_ec)

  # --- EC specific colors ---
  ec_specific_colors <- base::c(
    "Carboxylesterase (3.1.1.1)"      = "#2E86AB",
    "Lipase (3.1.1.3)"                = "#1B6B93",
    "Acylglycerol lipase (3.1.1.23)"  = "#3DA5D9",
    "MenH esterase (3.1.1.85)"        = "#5ABED6",
    "PETase (3.1.1.101)"              = "#1A535C",
    "Cutinase (3.1.1.74)"             = "#206A5D",
    "Arylesterase (3.1.1.2)"          = "#4ECDC4",
    "Esterase"                        = "#3B9EBF"
  )
  missing_ec <- base::setdiff(base::unique(sig$ec_specific),
                              base::names(ec_specific_colors))
  if (base::length(missing_ec) > 0L) {
    grey_vals <- base::rep("#C0C0C0", base::length(missing_ec))
    base::names(grey_vals) <- missing_ec
    ec_specific_colors <- base::c(ec_specific_colors, grey_vals)
  }

  # --- Substrate color palette ---
  known_colors <- base::c(
    "PET" = "#E07A5F", "PLA" = "#5B8FA8", "PUR" = "#81B29A", "PA"  = "#D4A76A",
    "PBAT/PET" = "#8B7EB8", "PHA" = "#8491B4", "PBS" = "#91D1C2", "PCL" = "#C4A77D"
  )
  known_dark <- base::c(
    "PET" = "#B85842", "PLA" = "#3D6B80", "PUR" = "#5A8A6E", "PA"  = "#A6803A",
    "PBAT/PET" = "#6558A0", "PHA" = "#5C6A8A", "PBS" = "#5A9E8A", "PCL" = "#8A6E4E"
  )
  all_substrates <- base::unique(sig$substrate)
  unknown_sub <- base::setdiff(all_substrates, base::names(known_colors))
  if (base::length(unknown_sub) > 0L) {
    extra_pal <- scales::hue_pal()(base::length(unknown_sub))
    base::names(extra_pal) <- unknown_sub
    known_colors <- base::c(known_colors, extra_pal)
    extra_dark <- colorspace::darken(extra_pal, 0.3)
    base::names(extra_dark) <- unknown_sub
    known_dark <- base::c(known_dark, extra_dark)
  }
  sub_colors      <- known_colors[all_substrates]
  sub_colors_dark <- known_dark[all_substrates]

  bg_col <- "white"
  theme_pub <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      text             = ggplot2::element_text(color = "grey20"),
      plot.background  = ggplot2::element_rect(fill = bg_col, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_col, color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#EEEDEA", linewidth = 0.3)
    )

  # ================================================================
  # Panel A: Genome overview + zoomed loci
  # ================================================================
  flank_bp   <- 5000
  sig_sorted <- sig |> dplyr::arrange(.data$start)

  loci    <- base::list()
  current <- base::list(start = sig_sorted$start[1] - flank_bp,
                        end   = sig_sorted$end[1]   + flank_bp,
                        hits  = sig_sorted$hit_id[1])
  for (i in 2:base::nrow(sig_sorted)) {
    new_start <- sig_sorted$start[i] - flank_bp
    new_end   <- sig_sorted$end[i]   + flank_bp
    if (new_start <= current$end + 50000) {
      current$end  <- base::max(current$end, new_end)
      current$hits <- base::c(current$hits, sig_sorted$hit_id[i])
    } else {
      loci    <- base::c(loci, base::list(current))
      current <- base::list(start = new_start, end = new_end,
                            hits = sig_sorted$hit_id[i])
    }
  }
  loci <- base::c(loci, base::list(current))
  loci <- loci[base::order(base::sapply(loci, function(x) x$start))]

  all_fill_colors <- base::c(sub_colors, "Other" = "#D6D3CC")

  locus_plots <- base::lapply(base::seq_along(loci), function(li) {
    loc <- loci[[li]]
    hit_names <- base::paste(
      base::sapply(base::strsplit(loc$hits, " \\("), `[`, 1),
      collapse = " / "
    )
    locus_label <- base::paste0(
      hit_names, "  (",
      base::round(loc$start / 1e3), "\u2013",
      base::round(loc$end / 1e3), " kb)"
    )
    genes_in_range <- d |>
      dplyr::filter(.data$end >= loc$start, .data$start <= loc$end) |>
      dplyr::mutate(
        forward    = .data$direction == "+",
        is_pazy    = .data$locus_tag %in% sig$locus_tag,
        substrate  = base::ifelse(.data$is_pazy,
                                  sig$substrate[base::match(.data$locus_tag, sig$locus_tag)],
                                  "Other"),
        pazy_label = base::ifelse(
          .data$is_pazy,
          base::paste0(sig$PAZy_family_id[base::match(.data$locus_tag, sig$locus_tag)],
                       " (", base::round(sig$PAZy_pident[base::match(.data$locus_tag, sig$locus_tag)]), "%)"),
          NA_character_
        ),
        gene_label = base::ifelse(!base::is.na(.data$gene) & base::nzchar(.data$gene),
                                  .data$gene, .data$locus_tag)
      )
    ggplot2::ggplot(genes_in_range,
                    ggplot2::aes(xmin = .data$start, xmax = .data$end,
                                y = 1, fill = .data$substrate,
                                forward = .data$forward)) +
      gggenes::geom_gene_arrow(
        arrowhead_height  = grid::unit(3.5, "mm"),
        arrowhead_width   = grid::unit(1.5, "mm"),
        arrow_body_height = grid::unit(3, "mm"),
        color = "grey45", size = 0.2
      ) +
      ggplot2::geom_text(
        data = genes_in_range |> dplyr::filter(.data$is_pazy),
        ggplot2::aes(x = (.data$start + .data$end) / 2, y = 1.22,
                     label = .data$pazy_label, color = .data$substrate),
        size = 2.3, fontface = "bold", show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = genes_in_range |> dplyr::filter(.data$is_pazy),
        ggplot2::aes(x = (.data$start + .data$end) / 2, y = 1.38,
                     label = base::ifelse(
                       base::nchar(.data$product) > 35,
                       base::paste0(base::substr(.data$product, 1, 33), "..."),
                       .data$product
                     ),
                     color = .data$substrate),
        size = 1.8, fontface = "italic", show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = genes_in_range |> dplyr::filter(!.data$is_pazy),
        ggplot2::aes(x = (.data$start + .data$end) / 2, y = 0.88,
                     label = .data$gene_label),
        size = 1.4, color = "grey65", angle = 45, hjust = 1, vjust = 1,
        check_overlap = TRUE
      ) +
      ggplot2::geom_text(
        data = genes_in_range |> dplyr::filter(.data$is_pazy),
        ggplot2::aes(x = (.data$start + .data$end) / 2, y = 0.88,
                     label = .data$gene_label, color = .data$substrate),
        size = 1.6, angle = 45, hjust = 1, vjust = 1,
        fontface = "bold", show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = all_fill_colors, drop = FALSE) +
      ggplot2::scale_color_manual(values = sub_colors_dark, guide = "none") +
      ggplot2::scale_x_continuous(
        labels = function(x) base::paste0(base::round(x / 1e3), " kb"),
        expand = ggplot2::expansion(mult = 0.04)
      ) +
      ggplot2::scale_y_continuous(limits = base::c(0.55, 1.55), expand = base::c(0, 0)) +
      ggplot2::labs(title = locus_label, x = NULL, y = NULL) +
      theme_pub +
      ggplot2::theme(
        panel.grid   = ggplot2::element_blank(),
        axis.text.y  = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(size = 6, color = "grey50"),
        plot.title   = ggplot2::element_text(face = "bold", size = 7,
                                             color = "grey30",
                                             margin = ggplot2::margin(0, 0, 1, 0)),
        legend.position = "none",
        plot.margin     = ggplot2::margin(2, 4, 2, 4)
      )
  })

  locus_positions <- base::data.frame(
    locus_id  = base::seq_along(loci),
    start     = base::sapply(loci, function(x) x$start),
    end       = base::sapply(loci, function(x) x$end),
    midpoint  = base::sapply(loci, function(x) (x$start + x$end) / 2),
    label     = base::sapply(loci, function(x)
      base::paste(base::sapply(base::strsplit(x$hits, " \\("), `[`, 1), collapse = "/")),
    substrate = base::sapply(loci, function(x) {
      idx <- base::match(x$hits[1], sig$hit_id)
      if (!base::is.na(idx)) sig$substrate[idx] else "PET"
    }),
    stringsAsFactors = FALSE
  )
  locus_positions$side <- base::ifelse(
    base::seq_len(base::nrow(locus_positions)) %% 2 == 1, 0.3, -0.3
  )

  p_overview <- ggplot2::ggplot(locus_positions) +
    ggplot2::labs(tag = "A") +
    ggplot2::annotate("rect", xmin = 0, xmax = genome_len,
                      ymin = -0.04, ymax = 0.04,
                      fill = "#D6D3CC", color = "#B8B4AB", linewidth = 0.3) +
    ggplot2::annotate("rect",
                      xmin = base::seq(500000, genome_len, by = 1e6),
                      xmax = base::pmin(base::seq(1e6, genome_len + 5e5, by = 1e6),
                                        genome_len),
                      ymin = -0.04, ymax = 0.04,
                      fill = "#C9C5BC", alpha = 0.5) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$start, xmax = .data$end,
                   ymin = -0.06, ymax = 0.06, fill = .data$substrate),
      color = "grey30", linewidth = 0.3, alpha = 0.8
    ) +
    ggplot2::scale_fill_manual(values = sub_colors, guide = "none") +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$midpoint, xend = .data$midpoint,
                   y = 0, yend = .data$side, color = .data$substrate),
      linewidth = 0.3, alpha = 0.6
    ) +
    ggplot2::scale_color_manual(values = sub_colors_dark, guide = "none") +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$midpoint,
                   y = .data$side + base::sign(.data$side) * 0.08,
                   label = .data$label, color = .data$substrate),
      size = 2.2, fontface = "bold"
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) base::paste0(x / 1e6, " Mb"),
      breaks = base::seq(0, genome_len, by = 5e5),
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::coord_cartesian(ylim = base::c(-0.55, 0.55)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_pub +
    ggplot2::theme(
      panel.grid         = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#EEEDEA", linewidth = 0.15),
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(size = 7, color = "grey50"),
      plot.tag           = ggplot2::element_text(face = "bold", size = 14,
                                                 color = "grey15"),
      plot.margin        = ggplot2::margin(2, 5, 0, 5)
    )

  # Assemble Panel A
  n_loci   <- base::length(locus_plots)
  ncol_val <- 4L

  p_legend_src <- ggplot2::ggplot(
    base::data.frame(
      x = base::seq_along(all_fill_colors), y = 1,
      s = base::factor(base::names(all_fill_colors),
                       levels = base::c(base::names(sub_colors), "Other"))
    ),
    ggplot2::aes(.data$x, .data$y, fill = .data$s)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = all_fill_colors, name = "Substrate") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1,
                                                 title.position = "left")) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.title   = ggplot2::element_text(face = "bold", size = 8),
      legend.text    = ggplot2::element_text(size = 7),
      legend.key.size = grid::unit(0.3, "cm"),
      plot.background = ggplot2::element_rect(fill = bg_col, color = NA)
    )
  p_legend <- patchwork::wrap_elements(cowplot::get_legend(p_legend_src))

  nrow_val   <- base::ceiling(n_loci / ncol_val)
  locus_rows <- base::list()
  for (r in base::seq_len(nrow_val)) {
    idx <- ((r - 1) * ncol_val + 1):base::min(r * ncol_val, n_loci)
    row_plots <- locus_plots[idx]
    while (base::length(row_plots) < ncol_val)
      row_plots <- base::c(row_plots, base::list(patchwork::plot_spacer()))
    locus_rows <- base::c(locus_rows, base::list(
      base::Reduce(`|`, row_plots) +
        patchwork::plot_layout(widths = base::rep(1, ncol_val))
    ))
  }
  p_genome <- base::Reduce(`/`, locus_rows)

  hit_order <- sig |> dplyr::arrange(.data$PAZy_bitscore) |> dplyr::pull(.data$hit_id)

  # ================================================================
  # Panel B: Bitscore + CLEAN EC class
  # ================================================================
  sig_b <- sig |>
    dplyr::mutate(
      hit_id     = base::factor(.data$hit_id, levels = hit_order),
      eval_label = base::ifelse(
        .data$PAZy_evalue < 1e-100,
        base::formatC(.data$PAZy_evalue, format = "e", digits = 0),
        base::formatC(.data$PAZy_evalue, format = "e", digits = 1)
      ),
      annot   = base::paste0(base::round(.data$PAZy_pident, 0),
                              "% id  |  E = ", .data$eval_label),
      ec_label = base::ifelse(!base::is.na(.data$clean_ec), .data$clean_ec, "")
    )

  p_score <- ggplot2::ggplot(sig_b, ggplot2::aes(y = .data$hit_id)) +
    ggplot2::geom_tile(
      ggplot2::aes(x = -18, fill = .data$ec_specific),
      width = 12, height = 0.75, color = "white", linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(
      values = ec_specific_colors,
      breaks = base::intersect(base::names(ec_specific_colors),
                               base::unique(sig_b$ec_specific)),
      name = "CLEAN EC"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      nrow = 1, title.position = "left",
      override.aes = base::list(alpha = 1), order = 1
    )) +
    ggplot2::geom_text(
      ggplot2::aes(x = -5, label = .data$ec_label),
      size = 2.0, color = "grey30", hjust = 1, family = "mono", fontface = "bold"
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_col(
      ggplot2::aes(x = .data$PAZy_bitscore, fill = .data$substrate,
                   alpha = .data$PAZy_pident),
      width = 0.6, color = NA
    ) +
    ggplot2::scale_fill_manual(values = sub_colors, name = "Substrate") +
    ggplot2::scale_alpha_continuous(range = base::c(0.25, 1), guide = "none") +
    ggplot2::guides(fill = ggplot2::guide_legend(
      nrow = 1, title.position = "left",
      override.aes = base::list(alpha = 1), order = 2
    )) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
      ggplot2::aes(x = -999, fill = .data$PAZy_pident),
      shape = 22, size = 0, show.legend = TRUE
    ) +
    ggplot2::scale_fill_gradient(
      low = "white", high = "grey20", name = "Identity (%)",
      limits = base::c(0, 100), breaks = base::c(0, 50, 100),
      guide = ggplot2::guide_colorbar(barwidth = 3, barheight = 0.4,
                                      title.position = "left", order = 3)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$PAZy_bitscore, label = .data$annot),
      hjust = -0.03, size = 2.1, color = "grey35", family = "mono"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = base::c(0.08, 0.45)),
      breaks = base::seq(0, 800, 200)
    ) +
    ggplot2::coord_cartesian(
      xlim = base::c(-30, base::max(sig_b$PAZy_bitscore, na.rm = TRUE) * 1.6)
    ) +
    ggplot2::labs(x = "Bit score", y = NULL) +
    theme_pub +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#EEEDEA", linewidth = 0.2),
      axis.text.y        = ggplot2::element_text(face = "bold", size = 7.5,
                                                 color = "grey30", hjust = 0),
      legend.position     = "bottom",
      legend.box          = "vertical",
      legend.title        = ggplot2::element_text(face = "bold", size = 6.5),
      legend.text         = ggplot2::element_text(size = 5.5),
      legend.key.size     = grid::unit(0.25, "cm"),
      legend.spacing.y    = grid::unit(0.05, "cm"),
      legend.margin       = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin   = ggplot2::margin(-3, 0, 0, 0),
      plot.margin         = ggplot2::margin(5, 8, 2, 5)
    )

  # ================================================================
  # Panel C: Catalytic site analysis
  # ================================================================

  # --- Read query FASTA or fall back to genbank_table translations ---
  query_fasta_path <- base::file.path(output_dir, "dnmb_module_pazy",
                                      "pazy_query_proteins.faa")
  if (!base::file.exists(query_fasta_path)) {
    query_fasta_path <- base::file.path(output_dir, "pazy_query_proteins.faa")
  }
  seqs_aa <- base::list()
  if (base::file.exists(query_fasta_path)) {
    fasta_lines  <- base::readLines(query_fasta_path)
    fasta_headers <- base::grep("^>", fasta_lines)
    for (i in base::seq_along(fasta_headers)) {
      h   <- fasta_headers[i]
      tag <- base::strsplit(base::sub("^>", "", fasta_lines[h]), " ")[[1]][1]
      end_idx <- if (i < base::length(fasta_headers)) fasta_headers[i + 1] - 1L
                 else base::length(fasta_lines)
      seqs_aa[[tag]] <- base::paste(fasta_lines[(h + 1):end_idx], collapse = "")
    }
  } else if ("translation" %in% base::names(genbank_table)) {
    for (j in base::seq_len(base::nrow(genbank_table))) {
      tag <- genbank_table$locus_tag[j]
      seq <- genbank_table$translation[j]
      if (!base::is.na(tag) && !base::is.na(seq) && base::nzchar(seq)) {
        seqs_aa[[tag]] <- base::as.character(seq)
      }
    }
  }

  flank_aa   <- 6L
  motif_data <- base::lapply(base::seq_len(base::nrow(sig)), function(i) {
    tag <- sig$locus_tag[i]
    s   <- seqs_aa[[tag]]
    if (base::is.null(s)) return(NULL)
    matches <- base::gregexpr("G.S.G", s)[[1]]
    if (matches[1] == -1L)
      matches <- base::gregexpr("[GA].S.[GA]", s)[[1]]
    if (matches[1] == -1L) return(NULL)
    m <- matches[1]
    ctx_start <- base::max(1L, m - flank_aa)
    ctx_end   <- base::min(base::nchar(s), m + 4L + flank_aa)
    base::data.frame(
      hit_id         = sig$hit_id[i],
      substrate      = sig$substrate[i],
      motif          = base::substr(s, m, m + 4L),
      motif_position = m,
      context        = base::substr(s, ctx_start, ctx_end),
      stringsAsFactors = FALSE
    )
  }) |> dplyr::bind_rows()

  sig <- dplyr::left_join(sig, motif_data, by = base::c("hit_id", "substrate"))

  # --- Read reference FASTA (needed for Panel C catalytic triad) ---
  pazy_module_dir <- .dnmb_db_module_dir("pazy", "current", cache_root = cache_root)
  ref_fasta_path  <- base::file.path(pazy_module_dir, "pazy_reference.faa")
  # Auto-install PAZy module if reference FASTA is missing
  if (!base::file.exists(ref_fasta_path)) {
    base::message("[DNMB] PAZy reference FASTA not found; attempting auto-install...")
    install_res <- base::tryCatch(
      dnmb_pazy_install_module(cache_root = cache_root),
      error = function(e) base::list(ok = FALSE)
    )
    if (base::isTRUE(install_res$ok)) {
      pazy_module_dir <- .dnmb_db_module_dir("pazy", "current", cache_root = cache_root)
      ref_fasta_path  <- base::file.path(pazy_module_dir, "pazy_reference.faa")
    }
  }
  has_ref_fasta <- base::file.exists(ref_fasta_path)
  ref_seqs <- base::list()
  if (has_ref_fasta) {
    ref_fasta   <- base::readLines(ref_fasta_path)
    ref_headers <- base::grep("^>", ref_fasta)
    for (i in base::seq_along(ref_headers)) {
      h       <- ref_headers[i]
      hdr     <- base::sub("^>", "", ref_fasta[h])
      pazy_id <- base::strsplit(hdr, " ")[[1]][1]
      end_idx <- if (i < base::length(ref_headers)) ref_headers[i + 1] - 1L
                 else base::length(ref_fasta)
      ref_seqs[[pazy_id]] <- base::paste(ref_fasta[(h + 1):end_idx], collapse = "")
    }
  } else {
    base::message("[DNMB] PAZy reference FASTA unavailable after install attempt; skipping Panel C")
  }

  # --- Alignment helper ---
  map_ref_to_query <- function(aln, ref_pos) {
    aln_q  <- base::as.character(Biostrings::alignedPattern(aln))
    aln_s  <- base::as.character(Biostrings::alignedSubject(aln))
    q_chars <- base::strsplit(aln_q, "")[[1]]
    s_chars <- base::strsplit(aln_s, "")[[1]]
    q_idx <- 0L; s_idx <- 0L
    for (k in base::seq_along(q_chars)) {
      if (s_chars[k] != "-") s_idx <- s_idx + 1L
      if (q_chars[k] != "-") q_idx <- q_idx + 1L
      if (s_idx == ref_pos) {
        if (q_chars[k] != "-") return(q_idx) else return(NA)
      }
    }
    NA
  }

  # --- Catalytic triad detection (requires reference FASTA) ---
  if (!has_ref_fasta || !base::length(ref_seqs)) {
    # Skip Panel C — assemble with Panels A + B only
    p_genome_full <- p_overview / p_genome / p_legend +
      patchwork::plot_layout(heights = base::c(0.35, 1, 0.08))
    p_genome_grob <- patchwork::wrap_elements(
      full = patchwork::patchworkGrob(p_genome_full)
    )

    combined <- p_genome_grob / p_score +
      patchwork::plot_layout(heights = base::c(1.05, 1)) +
      patchwork::plot_annotation(
        title    = "PAZy Substrate Map",
        subtitle = base::unique(d$contig[!base::is.na(d$contig)])[1],
        theme    = ggplot2::theme(
          plot.title    = ggplot2::element_text(face = "bold", size = 15,
                                               hjust = 0, color = "grey15"),
          plot.subtitle = ggplot2::element_text(size = 10, color = "#AAA79E",
                                               face = "italic", hjust = 0,
                                               margin = ggplot2::margin(b = 3)),
          plot.background = ggplot2::element_rect(fill = bg_col, color = NA)
        )
      )

    out_path <- base::file.path(plot_dir, "PAZy_overview.pdf")
    ggplot2::ggsave(out_path, combined, width = 17, height = 13,
                    device = grDevices::cairo_pdf)
    return(base::list(pdf = out_path))
  }

  triad_data <- base::lapply(base::seq_len(base::nrow(sig)), function(i) {
    tag    <- sig$locus_tag[i]
    q_seq  <- seqs_aa[[tag]]
    ref_id <- sig$PAZy_pazy_id[i]
    if (base::is.null(q_seq) || base::is.null(ref_seqs[[ref_id]])) return(NULL)
    r_seq   <- ref_seqs[[ref_id]]
    q_len   <- base::nchar(q_seq)
    r_len   <- base::nchar(r_seq)
    q_chars <- base::strsplit(q_seq, "")[[1]]
    r_chars <- base::strsplit(r_seq, "")[[1]]

    aln <- base::tryCatch(
      Biostrings::pairwiseAlignment(
        Biostrings::AAString(q_seq), Biostrings::AAString(r_seq),
        type = "local", substitutionMatrix = "BLOSUM62",
        gapOpening = 10, gapExtension = 0.5
      ),
      error = function(e) NULL
    )
    if (base::is.null(aln)) return(NULL)

    ref_gxsxg <- base::gregexpr("G.S.G", r_seq)[[1]]
    if (ref_gxsxg[1] == -1L)
      ref_gxsxg <- base::gregexpr("[GA].S.[GA]", r_seq)[[1]]

    q_gxsxg   <- base::gregexpr("G.S.G", q_seq)[[1]]
    if (q_gxsxg[1] == -1L)
      q_gxsxg <- base::gregexpr("[GA].S.[GA]", q_seq)[[1]]
    q_ser_pos <- if (q_gxsxg[1] != -1L) q_gxsxg[1] + 2L else NA

    ser_pos <- NA; ser_validated <- FALSE
    if (ref_gxsxg[1] != -1L) {
      for (m in base::as.integer(ref_gxsxg)) {
        ref_ser  <- m + 2L
        q_mapped <- map_ref_to_query(aln, ref_ser)
        if (!base::is.na(q_mapped) && q_mapped <= q_len &&
            q_chars[q_mapped] == "S") {
          ser_pos <- q_mapped; ser_validated <- TRUE; break
        }
      }
    }
    if (base::is.na(ser_pos) && !base::is.na(q_ser_pos)) ser_pos <- q_ser_pos

    asp_pos <- NA; asp_res <- NA_character_; asp_validated <- FALSE
    his_pos <- NA; his_validated <- FALSE
    rj <- NA  # track reference position for His search

    if (ref_gxsxg[1] != -1L) {
      ref_ser <- base::as.integer(ref_gxsxg[1]) + 2L
      for (offset in 20:base::min(250L, r_len - ref_ser)) {
        rj_cand <- ref_ser + offset
        if (rj_cand > r_len) break
        if (r_chars[rj_cand] %in% base::c("D", "E")) {
          q_mapped <- map_ref_to_query(aln, rj_cand)
          if (!base::is.na(q_mapped) && q_mapped <= q_len &&
              q_chars[q_mapped] == r_chars[rj_cand]) {
            asp_pos <- q_mapped; asp_res <- q_chars[q_mapped]
            asp_validated <- TRUE; rj <- rj_cand; break
          }
        }
      }
      if (base::is.na(asp_pos)) {
        for (offset in 20:base::min(250L, ref_ser - 1L)) {
          rj_cand <- ref_ser - offset
          if (rj_cand < 1L) break
          if (r_chars[rj_cand] %in% base::c("D", "E")) {
            q_mapped <- map_ref_to_query(aln, rj_cand)
            if (!base::is.na(q_mapped) && q_mapped <= q_len &&
                q_chars[q_mapped] == r_chars[rj_cand]) {
              asp_pos <- q_mapped; asp_res <- q_chars[q_mapped]
              asp_validated <- TRUE; rj <- rj_cand; break
            }
          }
        }
      }
      if (!base::is.na(asp_pos) && !base::is.na(rj)) {
        search_range <- (rj + 5L):base::min(r_len, rj + 200L)
        for (rj2 in search_range) {
          if (rj2 > r_len) break
          if (r_chars[rj2] == "H") {
            q_mapped <- map_ref_to_query(aln, rj2)
            if (!base::is.na(q_mapped) && q_mapped <= q_len &&
                q_chars[q_mapped] == "H") {
              his_pos <- q_mapped; his_validated <- TRUE; break
            }
          }
        }
      }
    }

    n_validated <- base::sum(base::c(ser_validated, asp_validated, his_validated))
    triad_type <- if (n_validated == 3L) {
      if (!base::is.na(asp_pos) && asp_pos < ser_pos) "D\u2013H\u2013S" else "S\u2013D\u2013H"
    } else if (n_validated == 2L) "partial"
    else if (n_validated == 1L) "Ser only"
    else "unverified"

    validation_str <- base::paste0(n_validated, "/3 aligned")

    get_ctx <- function(pos, seq_chars, fl = 2L) {
      if (base::is.na(pos)) return(NA_character_)
      base::paste(seq_chars[base::max(1L, pos - fl):base::min(base::length(seq_chars), pos + fl)],
                  collapse = "")
    }

    base::data.frame(
      hit_id        = sig$hit_id[i],
      locus_tag     = tag,
      seq_length    = q_len,
      ser_pos       = ser_pos,
      ser_ctx       = get_ctx(ser_pos, q_chars),
      ser_validated = ser_validated,
      asp_pos       = asp_pos,
      asp_res       = asp_res,
      asp_ctx       = get_ctx(asp_pos, q_chars),
      asp_validated = asp_validated,
      his_pos       = his_pos,
      his_ctx       = get_ctx(his_pos, q_chars),
      his_validated = his_validated,
      triad_type    = triad_type,
      validation    = validation_str,
      stringsAsFactors = FALSE
    )
  }) |> dplyr::bind_rows()

  if (base::nrow(triad_data) > 0L && "locus_tag" %in% base::names(triad_data)) {
    sig <- dplyr::left_join(sig, triad_data |> dplyr::select(-"locus_tag"),
                            by = "hit_id")
  }

  if (!"ser_pos" %in% base::names(sig)) {
    sig$ser_pos <- NA_integer_
    sig$asp_pos <- NA_integer_
    sig$his_pos <- NA_integer_
    sig$triad_type <- NA_character_
    sig$validation <- NA_character_
  }

  cat_df <- sig |>
    dplyr::filter(!base::is.na(.data$ser_pos))

  if (!base::nrow(cat_df)) {
    cat_df <- NULL
  } else {
    cat_df <- cat_df |>
      dplyr::mutate(hit_id = base::factor(.data$hit_id, levels = hit_order))
  }

  # --- Motif heatmap (skip if no catalytic triad data) ---
  if (base::is.null(cat_df)) {
    p_motif <- ggplot2::ggplot() + ggplot2::theme_void()
    p_triad <- ggplot2::ggplot() + ggplot2::theme_void()
  } else {
  motif_plot_data <- base::lapply(base::seq_len(base::nrow(cat_df)), function(i) {
    ctx <- cat_df$context[i]
    if (base::is.na(ctx)) return(NULL)
    chars <- base::strsplit(ctx, "")[[1]]
    motif_start_in_ctx <- flank_aa + 1L
    ser_pos_in_ctx     <- motif_start_in_ctx + 2L
    base::data.frame(
      hit_id    = cat_df$hit_id[i],
      substrate = cat_df$substrate[i],
      residue   = chars,
      pos       = base::seq_along(chars),
      is_motif  = base::seq_along(chars) >= motif_start_in_ctx &
                  base::seq_along(chars) <= (motif_start_in_ctx + 4L),
      is_ser    = base::seq_along(chars) == ser_pos_in_ctx,
      is_gly    = (base::seq_along(chars) == motif_start_in_ctx |
                   base::seq_along(chars) == (motif_start_in_ctx + 4L)),
      stringsAsFactors = FALSE
    )
  }) |> dplyr::bind_rows() |>
    dplyr::mutate(
      hit_id    = base::factor(.data$hit_id, levels = hit_order),
      res_color = dplyr::case_when(
        .data$is_ser ~ "#C0392B",
        .data$is_gly ~ "#2980B9",
        .data$is_motif ~ "#E07A5F",
        TRUE ~ "#B8B4AB"
      ),
      tile_fill = dplyr::case_when(
        .data$is_ser ~ "#FADBD8",
        .data$is_gly ~ "#D6EAF8",
        .data$is_motif ~ "#FEF0EC",
        TRUE ~ bg_col
      )
    )

  motif_labels <- cat_df |>
    dplyr::mutate(
      hit_id    = base::factor(.data$hit_id, levels = hit_order),
      pos_label = base::paste0(.data$motif, " @", .data$motif_position)
    )

  p_motif <- ggplot2::ggplot(motif_plot_data,
                             ggplot2::aes(x = .data$pos, y = .data$hit_id)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$tile_fill),
                       color = NA, width = 1, height = 0.6) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$residue, color = .data$res_color),
      size = 2.8, family = "mono",
      fontface = base::ifelse(motif_plot_data$is_motif, "bold", "plain")
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_text(
      data = motif_labels,
      ggplot2::aes(x = base::max(motif_plot_data$pos) + 1.2,
                   y = .data$hit_id, label = .data$pos_label),
      size = 1.8, color = "grey50", hjust = 0, family = "mono"
    ) +
    ggplot2::scale_x_continuous(breaks = NULL,
                                expand = ggplot2::expansion(mult = base::c(0.02, 0.28))) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = "G-x-S-x-G motif") +
    theme_pub +
    ggplot2::theme(
      panel.grid    = ggplot2::element_blank(),
      axis.text     = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 7.5,
                                            color = "grey35", hjust = 0.5,
                                            margin = ggplot2::margin(b = 2)),
      plot.margin   = ggplot2::margin(5, 5, 2, 2)
    )

  # --- Catalytic triad map ---
  triad_long <- cat_df |>
    dplyr::select("hit_id", "locus_tag", "seq_length",
                  "ser_pos", "ser_validated", "asp_pos", "asp_res",
                  "asp_validated", "his_pos", "his_validated",
                  "triad_type", "validation") |>
    tidyr::pivot_longer(
      cols      = base::c("ser_pos", "asp_pos", "his_pos"),
      names_to  = "residue_type",
      values_to = "position"
    ) |>
    dplyr::mutate(
      residue_label = dplyr::case_when(
        .data$residue_type == "ser_pos" ~ "S",
        .data$residue_type == "asp_pos" ~
          base::ifelse(!base::is.na(.data$asp_res) & .data$asp_res == "E", "E", "D"),
        .data$residue_type == "his_pos" ~ "H"
      ),
      validated = dplyr::case_when(
        .data$residue_type == "ser_pos" ~ .data$ser_validated,
        .data$residue_type == "asp_pos" ~ .data$asp_validated,
        .data$residue_type == "his_pos" ~ .data$his_validated
      ),
      pt_shape     = base::ifelse(.data$validated, 18L, 5L),
      pt_size      = base::ifelse(.data$validated, 3.8, 3.2),
      residue_type = base::factor(.data$residue_type,
                                  levels = base::c("ser_pos", "asp_pos", "his_pos"))
    ) |>
    dplyr::filter(!base::is.na(.data$position))

  triad_summary <- cat_df |>
    dplyr::mutate(
      summary_label = base::paste0(
        "S", .data$ser_pos,
        base::ifelse(!base::is.na(.data$asp_pos),
                     base::paste0(" ", base::ifelse(
                       !base::is.na(.data$asp_res) & .data$asp_res == "E", "E", "D"
                     ), .data$asp_pos), ""),
        base::ifelse(!base::is.na(.data$his_pos),
                     base::paste0(" H", .data$his_pos), "")
      ),
      info_label = base::paste0(.data$seq_length, " aa  ", .data$validation)
    )

  max_seq <- base::max(cat_df$seq_length, na.rm = TRUE)

  p_triad <- ggplot2::ggplot(cat_df, ggplot2::aes(y = .data$hit_id)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 1, xend = .data$seq_length, yend = .data$hit_id),
      color = "#E8E6E1", linewidth = 5.5, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = triad_long |> dplyr::filter(.data$validated),
      ggplot2::aes(x = .data$position, y = .data$hit_id,
                   color = .data$residue_type),
      size = 3.8, shape = 18
    ) +
    ggplot2::geom_text(
      data = triad_long |> dplyr::filter(.data$validated),
      ggplot2::aes(x = .data$position, y = .data$hit_id,
                   label = .data$residue_label),
      size = 2.0, color = "white", fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = triad_long |> dplyr::filter(!.data$validated),
      ggplot2::aes(x = .data$position, y = .data$hit_id,
                   color = .data$residue_type),
      size = 3.2, shape = 5, stroke = 0.7
    ) +
    ggplot2::geom_text(
      data = triad_long |> dplyr::filter(!.data$validated),
      ggplot2::aes(x = .data$position, y = .data$hit_id,
                   label = .data$residue_label),
      size = 1.8, color = "grey30", fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = triad_summary,
      ggplot2::aes(x = .data$seq_length + 5, y = .data$hit_id,
                   label = .data$summary_label),
      size = 2.0, color = "grey25", hjust = 0, fontface = "bold", family = "mono"
    ) +
    ggplot2::geom_text(
      data = triad_summary,
      ggplot2::aes(x = max_seq + 70, y = .data$hit_id, label = .data$info_label),
      size = 1.7, color = "grey50", hjust = 0, family = "mono"
    ) +
    ggplot2::scale_color_manual(
      values = base::c("ser_pos" = "#C0392B", "asp_pos" = "#2980B9",
                        "his_pos" = "#27AE60"),
      labels = base::c("ser_pos" = "Ser (nucleophile)",
                        "asp_pos" = "Asp/Glu (acid)",
                        "his_pos" = "His (base)"),
      name = "Catalytic residue"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = base::c(0.02, 0.22)),
      breaks = base::seq(0, 700, 100)
    ) +
    ggplot2::labs(
      x = "Residue position", y = NULL,
      subtitle = "Catalytic triad (alignment-validated)",
      caption = base::expression(
        base::phantom(0) * scriptstyle("filled = aligned with reference") *
        "      " * scriptstyle("open = unverified")
      )
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = base::list(size = 5, shape = 18),
      direction = "horizontal"
    )) +
    theme_pub +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#EEEDEA", linewidth = 0.2),
      axis.text.y        = ggplot2::element_text(face = "bold", size = 7,
                                                 color = "grey30", hjust = 0),
      axis.text.x        = ggplot2::element_text(size = 6, color = "grey50"),
      plot.subtitle      = ggplot2::element_text(face = "bold", size = 7.5,
                                                 color = "grey35", hjust = 0.5,
                                                 margin = ggplot2::margin(b = 2)),
      plot.caption       = ggplot2::element_text(size = 6, color = "grey45",
                                                 hjust = 0.5),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold", size = 7),
      legend.text      = ggplot2::element_text(size = 6.5),
      legend.key.size  = grid::unit(0.3, "cm"),
      legend.margin    = ggplot2::margin(0, 0, 0, 0),
      plot.margin      = ggplot2::margin(5, 5, 2, 2)
    )

  } # end if (!is.null(cat_df))

  # ================================================================
  # Assemble: A / (B | C)
  # ================================================================
  p_genome_full <- p_overview / p_genome / p_legend +
    patchwork::plot_layout(heights = base::c(0.35, 1, 0.08))
  p_genome_grob <- patchwork::wrap_elements(
    full = patchwork::patchworkGrob(p_genome_full)
  )

  p_score_tagged <- p_score +
    ggplot2::labs(tag = "B") +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14,
                                                    color = "grey15"))

  p_cat_combined <- (p_motif | p_triad) +
    patchwork::plot_layout(widths = base::c(0.45, 1))
  p_cat_tagged <- patchwork::wrap_elements(
    full = patchwork::patchworkGrob(p_cat_combined)
  )
  p_cat_tagged <- p_cat_tagged +
    ggplot2::labs(tag = "C") +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14,
                                                    color = "grey15"))

  bottom_row <- (p_score_tagged | p_cat_tagged) +
    patchwork::plot_layout(widths = base::c(1, 1.3))

  combined <- p_genome_grob / bottom_row +
    patchwork::plot_layout(heights = base::c(1.05, 1)) +
    patchwork::plot_annotation(
      title    = "PAZy Substrate Map",
      subtitle = base::unique(d$contig[!base::is.na(d$contig)])[1],
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 15,
                                             hjust = 0, color = "grey15"),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#AAA79E",
                                             face = "italic", hjust = 0,
                                             margin = ggplot2::margin(b = 3)),
        plot.background = ggplot2::element_rect(fill = bg_col, color = NA)
      )
    )

  out_path <- base::file.path(plot_dir, "PAZy_overview.pdf")
  ggplot2::ggsave(out_path, combined, width = 17, height = 13,
                  device = grDevices::cairo_pdf)

  base::list(pdf = out_path)
}


.dnmb_plot_pazy_module <- function(genbank_table, output_dir, cache_root = NULL) {
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  if (!base::nrow(tbl) || !"PAZy_family_id" %in% base::names(tbl)) {
    return(NULL)
  }
  tbl <- tbl[!is.na(tbl$PAZy_family_id) & base::nzchar(tbl$PAZy_family_id), , drop = FALSE]
  if (!base::nrow(tbl)) {
    return(NULL)
  }
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "PAZy_overview.pdf")

  .dnmb_plot_pazy_pub(genbank_table, output_dir = output_dir, cache_root = cache_root)

  list(pdf = pdf_path)
}

.dnmb_plot_merops_linear_module <- function(genbank_table, output_dir) {
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  if (!base::nrow(tbl) || !"MEROPS_family_id" %in% base::names(tbl)) {
    return(NULL)
  }
  tbl <- .dnmb_filter_merops_plot_hits(tbl)
  tbl <- tbl[!is.na(tbl$MEROPS_family_id) & base::nzchar(tbl$MEROPS_family_id), , drop = FALSE]
  if (!base::nrow(tbl)) {
    return(NULL)
  }
  tbl$merops_class <- .dnmb_merops_class_from_family(tbl$MEROPS_family_id)
  label_tbl <- .dnmb_merops_label_subset(tbl, max_labels_per_class = 4L)
  p <- ggplot2::ggplot(tbl) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = 0.82, ymax = 1.18, fill = .data$merops_class),
      color = NA,
      alpha = 0.85
    ) +
    ggplot2::geom_text(
      data = label_tbl,
      ggplot2::aes(x = .data$midpoint, y = 1.28, label = .data$locus_tag),
      angle = 90,
      size = 2.4,
      hjust = 0,
      show.legend = FALSE
    ) +
    ggplot2::facet_grid(merops_class ~ contig, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = .dnmb_merops_class_palette()) +
    ggplot2::labs(
      title = "MEROPS Class Tracks",
      subtitle = "Representative significant hits (evalue <= 1e-4)",
      x = "Genome coordinate (bp)",
      y = NULL,
      fill = "Class"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(12, 12, 18, 12)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "merops_class_tracks.pdf")
  .dnmb_module_plot_save(p, pdf_path, width = 14, height = 8)
  list(pdf = pdf_path)
}

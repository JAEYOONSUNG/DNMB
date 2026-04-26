.dnmb_integrated_defense_short_label <- function(x, max_len = 20L) {
  x <- as.character(x)
  ifelse(nchar(x) <= max_len, x, paste0(substr(x, 1, max_len - 1L), "..."))
}

.dnmb_integrated_defense_palette <- function(values, sources = NULL) {
  values <- as.character(values)
  sources <- as.character(sources %||% rep(NA_character_, length(values)))
  keep <- !is.na(values) & nzchar(values)
  values <- values[keep]
  sources <- sources[keep]
  if (!length(values)) {
    return(character())
  }
  key <- !duplicated(values)
  values_u <- values[key]
  sources_u <- sources[key]
  base_cols <- grDevices::hcl.colors(length(values_u), palette = "Dark 3")
  alpha_vals <- ifelse(
    sources_u == "DefenseFinder",
    1.0,
    ifelse(sources_u == "DefensePredictor", 0.25, 0.5)
  )
  cols <- vapply(seq_along(base_cols), function(i) {
    grDevices::adjustcolor(base_cols[[i]], alpha.f = alpha_vals[[i]])
  }, character(1))
  stats::setNames(cols, values_u)
}

.dnmb_integrated_defense_max_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(is.na(x))) {
    return(NA_real_)
  }
  max(x, na.rm = TRUE)
}

.dnmb_integrated_defense_clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("^DP_candidate::", "DP:", x)
  x <- gsub("^DefensePredictor_only:", "DP:", x)
  x <- gsub("\\[.*?\\]", "", x)
  x <- gsub("domain-containing", "dc", x, ignore.case = TRUE)
  x <- gsub("family protein", "", x, ignore.case = TRUE)
  x <- gsub(" family", "", x, ignore.case = TRUE)
  x <- gsub(" protein", "", x, ignore.case = TRUE)
  x <- gsub("CRISPR-associated", "Cas", x, ignore.case = TRUE)
  x <- gsub("restriction-modification", "RM", x, ignore.case = TRUE)
  x <- gsub("restriction endonuclease", "REase", x, ignore.case = TRUE)
  x <- gsub("methyltransferase", "MTase", x, ignore.case = TRUE)
  x <- gsub("ATP-dependent", "ATP", x, ignore.case = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

.dnmb_integrated_defense_display_name <- function(name, source = NULL, max_len = 20L) {
  x <- .dnmb_integrated_defense_clean_text(name)
  src <- as.character(source %||% rep("", length(x)))
  if (length(src) != length(x)) {
    src <- rep(src, length.out = length(x))
  }
  x[src == "DefensePredictor"] <- "DP candidate"
  .dnmb_integrated_defense_short_label(x, max_len = max_len)
}

.dnmb_integrated_defense_gene_name <- function(hit_rows) {
  n <- nrow(hit_rows)
  out <- rep(NA_character_, n)
  if ("DefenseFinder_gene_name" %in% names(hit_rows)) {
    val <- as.character(hit_rows$DefenseFinder_gene_name)
    keep <- !is.na(val) & nzchar(val)
    out[keep] <- val[keep]
  }
  if ("PADLOC_protein_name" %in% names(hit_rows)) {
    val <- as.character(hit_rows$PADLOC_protein_name)
    keep <- (is.na(out) | !nzchar(out)) & !is.na(val) & nzchar(val)
    out[keep] <- val[keep]
  }
  if ("REBASEfinder_enzyme_role" %in% names(hit_rows)) {
    val <- as.character(hit_rows$REBASEfinder_enzyme_role)
    keep <- (is.na(out) | !nzchar(out)) & !is.na(val) & nzchar(val)
    out[keep] <- val[keep]
  }
  if ("annotation_name" %in% names(hit_rows)) {
    val <- as.character(hit_rows$annotation_name)
    keep <- (is.na(out) | !nzchar(out)) & !is.na(val) & nzchar(val)
    out[keep] <- val[keep]
  }
  if ("locus_tag" %in% names(hit_rows)) {
    val <- as.character(hit_rows$locus_tag)
    keep <- is.na(out) | !nzchar(out)
    out[keep] <- val[keep]
  }
  .dnmb_integrated_defense_display_name(out, source = hit_rows$system_group_source, max_len = 14L)
}

.dnmb_integrated_defense_euler_label_positions <- function(ve,
                                                           set_sizes,
                                                           set_colors,
                                                           set_short,
                                                           venn_x,
                                                           venn_y,
                                                           venn_w,
                                                           venn_h) {
  label_tbl <- data.frame(
    tool = names(set_sizes),
    label = paste0(unname(set_short[names(set_sizes)]), " (", set_sizes, ")"),
    color = unname(set_colors[names(set_sizes)]),
    anchor_x = NA_real_,
    anchor_y = NA_real_,
    x = NA_real_,
    y = NA_real_,
    stringsAsFactors = FALSE
  )

  preferred_dirs <- list(
    DefenseFinder = c(0.00, 1.00),
    REBASE = c(1.00, 0.35),
    PADLOC = c(-1.00, -0.05),
    DefensePredictor = c(-0.75, -0.95),
    AntiDefenseFinder = c(0.00, 1.00),
    dbAPIS = c(1.00, -0.20),
    AcrFinder = c(-1.00, -0.20)
  )

  if (!is.null(ve$centers) && !is.null(ve$diameters)) {
    for (tool in label_tbl$tool) {
      if (!(tool %in% rownames(ve$centers)) || !(tool %in% names(ve$diameters))) {
        next
      }
      cx <- ve$centers[tool, "x"]
      cy <- ve$centers[tool, "y"]
      d <- ve$diameters[[tool]]
      vec <- preferred_dirs[[tool]] %||% c(0, 1)
      vec <- vec / sqrt(sum(vec^2))
      edge_x <- cx + vec[[1]] * (d / 2 + 0.02)
      edge_y <- cy + vec[[2]] * (d / 2 + 0.02)
      posx <- cx + vec[[1]] * (d / 2 + 0.18)
      posy <- cy + vec[[2]] * (d / 2 + 0.18)
      label_tbl$anchor_x[label_tbl$tool == tool] <- venn_x + venn_w * edge_x
      label_tbl$anchor_y[label_tbl$tool == tool] <- venn_y + venn_h * edge_y
      label_tbl$x[label_tbl$tool == tool] <- venn_x + venn_w * posx
      label_tbl$y[label_tbl$tool == tool] <- venn_y + venn_h * posy
    }
  } else if (!is.null(ve$ellipses)) {
    ell <- as.data.frame(ve$ellipses)
    if (all(c("h", "k", "a", "b") %in% names(ell)) && nrow(ell)) {
      x_min <- min(ell$h - ell$a, na.rm = TRUE)
      x_max <- max(ell$h + ell$a, na.rm = TRUE)
      y_min <- min(ell$k - ell$b, na.rm = TRUE)
      y_max <- max(ell$k + ell$b, na.rm = TRUE)
      x_span <- max(1e-6, x_max - x_min)
      y_span <- max(1e-6, y_max - y_min)
      pad <- 0.06
      for (tool in label_tbl$tool) {
        if (!(tool %in% rownames(ell))) {
          next
        }
        cx <- ell[tool, "h"]
        cy <- ell[tool, "k"]
        a <- max(1e-6, ell[tool, "a"])
        b <- max(1e-6, ell[tool, "b"])
        vec <- preferred_dirs[[tool]] %||% c(0, 1)
        vec <- vec / sqrt(sum(vec^2))
        radius <- 1 / sqrt((vec[[1]] / a)^2 + (vec[[2]] / b)^2)
        edge_x <- (cx + vec[[1]] * (radius + pad * x_span) - x_min) / x_span
        edge_y <- (cy + vec[[2]] * (radius + pad * y_span) - y_min) / y_span
        posx <- (cx + vec[[1]] * (radius + 0.20 * x_span) - x_min) / x_span
        posy <- (cy + vec[[2]] * (radius + 0.20 * y_span) - y_min) / y_span
        label_tbl$anchor_x[label_tbl$tool == tool] <- venn_x + venn_w * edge_x
        label_tbl$anchor_y[label_tbl$tool == tool] <- venn_y + venn_h * edge_y
        label_tbl$x[label_tbl$tool == tool] <- venn_x + venn_w * posx
        label_tbl$y[label_tbl$tool == tool] <- venn_y + venn_h * posy
      }
    }
  }

  miss <- is.na(label_tbl$x) | is.na(label_tbl$y)
  if (any(miss)) {
    fallback <- data.frame(
      tool = c("DefenseFinder", "REBASE", "PADLOC", "DefensePredictor", "AntiDefenseFinder", "dbAPIS", "AcrFinder"),
      x = c(0.12, 0.23, 0.09, 0.21, 0.12, 0.23, 0.08),
      y = c(0.79, 0.79, 0.50, 0.49, 0.79, 0.52, 0.52),
      stringsAsFactors = FALSE
    )
    idx <- match(label_tbl$tool[miss], fallback$tool)
    label_tbl$anchor_x[miss] <- fallback$x[idx]
    label_tbl$anchor_y[miss] <- fallback$y[idx]
    label_tbl$x[miss] <- fallback$x[idx]
    label_tbl$y[miss] <- fallback$y[idx]
  }

  # Mild pairwise repulsion in normalized page space.
  if (nrow(label_tbl) >= 2L) {
    w_est <- pmax(0.060, nchar(label_tbl$label) * 0.0095)
    h_est <- rep(0.028, nrow(label_tbl))
    for (pass in seq_len(8L)) {
      for (i in seq_len(nrow(label_tbl) - 1L)) {
        for (j in seq.int(i + 1L, nrow(label_tbl))) {
          dx <- label_tbl$x[[j]] - label_tbl$x[[i]]
          dy <- label_tbl$y[[j]] - label_tbl$y[[i]]
          overlap_x <- (w_est[[i]] + w_est[[j]]) / 2 - abs(dx)
          overlap_y <- (h_est[[i]] + h_est[[j]]) / 2 - abs(dy)
          if (overlap_x > 0 && overlap_y > 0) {
            if (abs(dx) >= abs(dy)) {
              shift <- overlap_x / 2 + 0.014
              signx <- ifelse(dx >= 0, 1, -1)
              label_tbl$x[[i]] <- label_tbl$x[[i]] - signx * shift / 2
              label_tbl$x[[j]] <- label_tbl$x[[j]] + signx * shift / 2
            } else {
              shift <- overlap_y / 2 + 0.012
              signy <- ifelse(dy >= 0, 1, -1)
              label_tbl$y[[i]] <- label_tbl$y[[i]] - signy * shift / 2
              label_tbl$y[[j]] <- label_tbl$y[[j]] + signy * shift / 2
            }
          }
        }
      }
    }
  }

  label_tbl$x <- pmin(0.34, pmax(0.02, label_tbl$x))
  label_tbl$y <- pmin(0.88, pmax(0.38, label_tbl$y))
  label_tbl
}

.dnmb_plot_integrated_defense_inventory <- function(system_summary, palette) {
  inventory_tbl <- system_summary |>
    dplyr::group_by(.data$display_name, .data$system_group_source) |>
    dplyr::summarise(
      n_systems = dplyr::n(),
      n_genes = sum(.data$genes_count, na.rm = TRUE),
      best_score = max(.data$weight, na.rm = TRUE),
      mean_wholeness = mean(.data$wholeness, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_systems), dplyr::desc(.data$n_genes), .data$display_name)

  inventory_tbl$display_name <- factor(inventory_tbl$display_name, levels = rev(inventory_tbl$display_name))
  label_w <- max(3.0, max(nchar(as.character(inventory_tbl$display_name)), na.rm = TRUE) * 0.22)
  inventory_tbl$label_x <- -label_w / 2

  ggplot2::ggplot(inventory_tbl, ggplot2::aes(y = .data$display_name)) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = -label_w,
        xmax = 0,
        ymin = as.numeric(.data$display_name) - 0.31,
        ymax = as.numeric(.data$display_name) + 0.31,
        fill = .data$display_name
      ),
      color = "grey35",
      linewidth = 0.22,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$label_x, label = .data$display_name),
      hjust = 0.5,
      size = 2.75,
      color = "white"
    ) +
    ggplot2::geom_col(
      ggplot2::aes(x = .data$n_systems, fill = .data$display_name),
      width = 0.62,
      color = "grey35",
      linewidth = 0.22,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$n_systems - 0.06, label = .data$n_systems),
      hjust = 1,
      size = 2.8,
      color = "white",
      fontface = "bold"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$n_systems + 0.10,
        label = paste0(.data$n_genes, " genes | score ", format(round(.data$best_score, 1), nsmall = 1), " | wholeness ", format(round(.data$mean_wholeness, 2), nsmall = 2))
      ),
      hjust = 0,
      size = 2.45,
      color = "grey30"
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.00, 0.34)),
      breaks = seq(0, max(inventory_tbl$n_systems), by = 1)
    ) +
    ggplot2::labs(title = NULL, x = "System loci detected", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 4, 4, 14)
    )
}

.dnmb_integrated_defense_assign_label_lanes <- function(df,
                                                        chars_per_full_width = 95,
                                                        min_bp = 28000,
                                                        pad_bp = 6000,
                                                        label_bp_scale = 1.15) {
  if (!nrow(df)) {
    df$label_lane <- integer()
    return(df)
  }
  df <- df[order(df$midpoint, df$start, df$end), , drop = FALSE]
  contig_len <- suppressWarnings(as.numeric(df$contig_length[[1]]))
  if (!is.finite(contig_len) || contig_len <= 0) {
    contig_len <- max(df$end, na.rm = TRUE)
  }
  label_chars <- pmax(nchar(as.character(df$context_label)), 4)
  char_bp <- contig_len / chars_per_full_width
  label_width_bp <- pmax(
    (label_chars + 2) * char_bp * label_bp_scale,
    (df$end - df$start) * 0.90,
    min_bp
  )
  interval_start <- pmax(0, df$midpoint - label_width_bp / 2 - pad_bp)
  interval_end <- df$midpoint + label_width_bp / 2 + pad_bp
  lane_end <- numeric()
  lane <- integer(nrow(df))
  for (i in seq_len(nrow(df))) {
    placed <- FALSE
    if (length(lane_end)) {
      for (k in seq_along(lane_end)) {
        if (interval_start[[i]] > lane_end[[k]]) {
          lane[[i]] <- k
          lane_end[[k]] <- interval_end[[i]]
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      lane[[i]] <- length(lane_end) + 1L
      lane_end <- c(lane_end, interval_end[[i]])
    }
  }
  df$label_lane <- lane
  df
}

.dnmb_plot_integrated_defense_context <- function(genbank_table,
                                                  output_dir,
                                                  defense_palette,
                                                  legend_position = "none") {
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  if (!nrow(tbl)) {
    return(NULL)
  }
  defense_tbl <- tbl[!is.na(tbl$DefenseFinder_system_id) & nzchar(tbl$DefenseFinder_system_id), , drop = FALSE]
  if (!nrow(defense_tbl)) {
    return(NULL)
  }

  defense_windows <- defense_tbl |>
    dplyr::group_by(.data$contig, .data$DefenseFinder_system_id, .data$DefenseFinder_system_subtype) |>
    dplyr::summarise(
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      genes_count = dplyr::n(),
      score = max(.data$DefenseFinder_system_score, na.rm = TRUE),
      wholeness = max(suppressWarnings(as.numeric(.data$DefenseFinder_system_wholeness)), na.rm = TRUE),
      profile_cov = if ("DefenseFinder_hit_profile_cov" %in% names(defense_tbl)) mean(suppressWarnings(as.numeric(.data$DefenseFinder_hit_profile_cov)), na.rm = TRUE) else NA_real_,
      seq_cov = if ("DefenseFinder_hit_seq_cov" %in% names(defense_tbl)) mean(suppressWarnings(as.numeric(.data$DefenseFinder_hit_seq_cov)), na.rm = TRUE) else NA_real_,
      context_label = dplyr::first(.data$context_label),
      system_group_source = dplyr::first(.data$system_group_source),
      .groups = "drop"
    )
  if (!nrow(defense_windows)) {
    return(NULL)
  }

  contig_lengths <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  contig_lengths$track <- 1
  defense_windows$track <- 1
  defense_windows$midpoint <- (defense_windows$start + defense_windows$end) / 2
  defense_windows$coverage_mean <- rowMeans(cbind(defense_windows$profile_cov, defense_windows$seq_cov), na.rm = TRUE)
  defense_windows$coverage_mean[!is.finite(defense_windows$coverage_mean)] <- NA_real_
  defense_windows$replicon_class <- ifelse(grepl("plasmid", defense_windows$contig, ignore.case = TRUE), "Plasmid", "Chromosome")
  defense_windows$color_value <- unname(defense_palette[match(defense_windows$DefenseFinder_system_subtype, names(defense_palette))])
  defense_windows$contig_length <- contig_lengths$length_bp[match(defense_windows$contig, contig_lengths$contig)]
  defense_windows <- dplyr::bind_rows(lapply(split(defense_windows, defense_windows$contig), .dnmb_integrated_defense_assign_label_lanes))
  defense_windows <- defense_windows[order(defense_windows$contig, defense_windows$midpoint), , drop = FALSE]

  dup_n <- ave(defense_windows$context_label, defense_windows$context_label, FUN = length)
  dup_i <- ave(seq_len(nrow(defense_windows)), defense_windows$context_label, FUN = seq_along)
  defense_windows$context_label <- ifelse(dup_n > 1, paste0(defense_windows$context_label, " (", dup_i, ")"), defense_windows$context_label)
  defense_windows$label_y <- 1.14 + 0.065 * (defense_windows$label_lane - 1L)
  max_lane <- max(defense_windows$label_lane, na.rm = TRUE)
  type_breaks <- names(defense_palette)
  n_type <- length(type_breaks)
  fill_rows <- max(1L, ceiling(n_type / 10L))

  contig_map <- contig_lengths$length_bp[match(defense_windows$contig, contig_lengths$contig)]
  offset_bp <- pmax(1200, contig_map * 0.0100)
  coverage_bg_tbl <- rbind(
    data.frame(contig = defense_windows$contig, x = defense_windows$midpoint - offset_bp / 2, xend = defense_windows$midpoint + offset_bp / 2, y = 0.915, yend = 0.915, stringsAsFactors = FALSE),
    data.frame(contig = defense_windows$contig, x = defense_windows$midpoint - offset_bp / 2, xend = defense_windows$midpoint + offset_bp / 2, y = 0.875, yend = 0.875, stringsAsFactors = FALSE)
  )
  coverage_tbl <- rbind(
    data.frame(contig = defense_windows$contig, x = defense_windows$midpoint - offset_bp / 2, xend = defense_windows$midpoint - offset_bp / 2 + offset_bp * pmax(0, pmin(1, defense_windows$profile_cov)), y = 0.915, yend = 0.915, metric = "P cov", stringsAsFactors = FALSE),
    data.frame(contig = defense_windows$contig, x = defense_windows$midpoint - offset_bp / 2, xend = defense_windows$midpoint - offset_bp / 2 + offset_bp * pmax(0, pmin(1, defense_windows$seq_cov)), y = 0.875, yend = 0.875, metric = "S cov", stringsAsFactors = FALSE)
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contig_lengths,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = .data$track, yend = .data$track),
      linewidth = 3.0,
      color = "grey88",
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = defense_windows,
      ggplot2::aes(x = .data$start, xend = .data$end, y = .data$track, yend = .data$track),
      linewidth = 2.4,
      lineend = "butt",
      color = defense_windows$color_value
    ) +
    ggplot2::geom_point(
      data = defense_windows,
      ggplot2::aes(x = .data$midpoint, y = .data$track, fill = .data$DefenseFinder_system_subtype, size = .data$score, shape = .data$replicon_class),
      color = "grey20",
      stroke = 0.35,
      show.legend = TRUE
    ) +
    ggplot2::geom_segment(
      data = defense_windows,
      ggplot2::aes(x = .data$midpoint, xend = .data$midpoint, y = 1.05, yend = .data$label_y - 0.02),
      linewidth = 0.35,
      color = grDevices::adjustcolor(defense_windows$color_value, alpha.f = 0.7),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = defense_windows,
      ggplot2::aes(x = .data$midpoint, y = .data$label_y, label = .data$context_label),
      size = 2.55,
      color = defense_windows$color_value,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = coverage_bg_tbl,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 1.5,
      lineend = "round",
      color = "grey82",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = coverage_tbl,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend, color = .data$metric),
      linewidth = 1.5,
      lineend = "round",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~contig, scales = "free_x", ncol = 1) +
    ggplot2::scale_fill_manual(values = defense_palette, name = "Defense subtype") +
    ggplot2::scale_color_manual(values = c("P cov" = "#2563EB", "S cov" = "#D97706"), name = "Coverage") +
    ggplot2::scale_shape_manual(values = c(Chromosome = 21, Plasmid = 24), name = "Replicon") +
    ggplot2::scale_size_continuous(range = c(4.0, 8.0), breaks = sort(unique(defense_windows$score)), name = "Score") +
    ggplot2::labs(
      title = "Integrated defense genome layout",
      x = "Genome coordinate (bp)",
      y = NULL
    ) +
    ggplot2::scale_x_continuous(labels = scales::label_comma(), expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::scale_y_continuous(
      limits = c(0.74, 1.18 + 0.065 * max_lane + 0.03),
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 0, b = 0)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.direction = "horizontal",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(0.08, "cm"),
      legend.spacing.x = grid::unit(0.14, "cm"),
      legend.box.spacing = grid::unit(0.10, "cm"),
      legend.key.width = grid::unit(0.42, "cm"),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.title = ggplot2::element_text(size = 9, face = "bold"),
      legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(4, 8, 8, 2)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 4,
        nrow = fill_rows,
        byrow = TRUE,
        override.aes = list(shape = 22, size = 5, color = NA)
      ),
      size = ggplot2::guide_legend(order = 1, nrow = 1, byrow = TRUE),
      shape = ggplot2::guide_legend(order = 2, nrow = 1, byrow = TRUE),
      color = ggplot2::guide_legend(order = 3, nrow = 1, byrow = TRUE)
    )
}

.dnmb_integrated_defense_overlap_data <- function(tbl,
                                                  dp_threshold = .dnmb_defensepredictor_default_threshold(),
                                                  defensefinder_activity = NULL) {
  gt <- as.data.frame(tbl, stringsAsFactors = FALSE)
  gt$DefenseFinder_system_activity <- if ("DefenseFinder_system_activity" %in% names(gt)) {
    as.character(gt$DefenseFinder_system_activity)
  } else {
    rep(NA_character_, nrow(gt))
  }
  gt$DefensePredictor_mean_log_odds <- if ("DefensePredictor_mean_log_odds" %in% names(gt)) {
    suppressWarnings(as.numeric(gt$DefensePredictor_mean_log_odds))
  } else {
    rep(NA_real_, nrow(gt))
  }
  gt$in_defensepredictor <- !is.na(gt$DefensePredictor_mean_log_odds) & gt$DefensePredictor_mean_log_odds >= as.numeric(dp_threshold)[1]
  gt$in_padloc <- if ("PADLOC_system" %in% names(gt)) {
    !is.na(gt$PADLOC_system) & nzchar(as.character(gt$PADLOC_system))
  } else {
    rep(FALSE, nrow(gt))
  }
  gt$in_defensefinder <- if ("DefenseFinder_system_id" %in% names(gt)) {
    keep <- !is.na(gt$DefenseFinder_system_id) & nzchar(as.character(gt$DefenseFinder_system_id))
    if (!is.null(defensefinder_activity)) {
      keep <- keep &
        !is.na(gt$DefenseFinder_system_activity) &
        gt$DefenseFinder_system_activity == as.character(defensefinder_activity)[1]
    }
    keep
  } else {
    rep(FALSE, nrow(gt))
  }
  gt$in_rebase <- if ("REBASEfinder_family_id" %in% names(gt)) {
    !is.na(gt$REBASEfinder_family_id) & nzchar(as.character(gt$REBASEfinder_family_id))
  } else {
    rep(FALSE, nrow(gt))
  }
  gt$in_dbapis <- if ("dbAPIS_family_id" %in% names(gt)) {
    !is.na(gt$dbAPIS_family_id) & nzchar(as.character(gt$dbAPIS_family_id))
  } else {
    rep(FALSE, nrow(gt))
  }
  gt$in_acrfinder <- if ("AcrFinder_family_id" %in% names(gt)) {
    !is.na(gt$AcrFinder_family_id) & nzchar(as.character(gt$AcrFinder_family_id))
  } else {
    rep(FALSE, nrow(gt))
  }
  if (!is.null(defensefinder_activity) && identical(as.character(defensefinder_activity)[1], "Anti-defense")) {
    gt$in_defensepredictor <- rep(FALSE, nrow(gt))
    gt$in_padloc <- rep(FALSE, nrow(gt))
    gt$in_rebase <- rep(FALSE, nrow(gt))
  }
  gt$accession <- as.character(gt$protein_id %||% NA_character_)
  if ("DefensePredictor_product_accession" %in% names(gt)) {
    missing_acc <- is.na(gt$accession) | !nzchar(gt$accession)
    gt$accession[missing_acc] <- as.character(gt$DefensePredictor_product_accession[missing_acc])
  }
  missing_acc <- is.na(gt$accession) | !nzchar(gt$accession)
  gt$accession[missing_acc] <- as.character(gt$locus_tag[missing_acc])
  gt <- gt[!is.na(gt$accession) & nzchar(gt$accession), , drop = FALSE]

  anti_mode <- !is.null(defensefinder_activity) && identical(as.character(defensefinder_activity)[1], "Anti-defense")
  tool_map <- if (anti_mode) {
    c(
      in_defensefinder = "AntiDefenseFinder",
      in_dbapis = "dbAPIS",
      in_acrfinder = "AcrFinder"
    )
  } else {
    c(
      in_defensepredictor = "DefensePredictor",
      in_padloc = "PADLOC",
      in_defensefinder = "DefenseFinder",
      in_rebase = "REBASE"
    )
  }
  long_parts <- lapply(names(tool_map), function(col) {
    keep <- gt[[col]] %in% TRUE
    if (!any(keep)) {
      return(NULL)
    }
    data.frame(
      accession = unique(as.character(gt$accession[keep])),
      tool = tool_map[[col]],
      stringsAsFactors = FALSE
    )
  })
  membership <- dplyr::bind_rows(long_parts)
  membership <- unique(membership)
  list(genbank = gt, membership = membership)
}

.dnmb_integrated_defense_build_members <- function(tbl,
                                                   dp_threshold = .dnmb_defensepredictor_default_threshold(),
                                                   region_gap = 25000L,
                                                   defensefinder_activity = NULL) {
  overlap <- .dnmb_integrated_defense_overlap_data(
    tbl,
    dp_threshold = dp_threshold,
    defensefinder_activity = defensefinder_activity
  )
  gt <- overlap$genbank
  anti_mode <- !is.null(defensefinder_activity) && identical(as.character(defensefinder_activity)[1], "Anti-defense")
  hit_rows <- gt[
    if (anti_mode) {
      gt$in_defensefinder | gt$in_dbapis | gt$in_acrfinder
    } else {
      gt$in_defensepredictor | gt$in_padloc | gt$in_defensefinder | gt$in_rebase
    },
    ,
    drop = FALSE
  ]
  if (!nrow(hit_rows)) {
    return(list(members = data.frame(), membership = overlap$membership))
  }

  hit_rows$seqid <- as.character(hit_rows$contig)
  hit_rows$strand <- as.character(hit_rows$direction)
  hit_rows$start <- suppressWarnings(as.numeric(hit_rows$start))
  hit_rows$end <- suppressWarnings(as.numeric(hit_rows$end))
  hit_rows$annotation_name <- as.character(hit_rows$product)
  support_cols <- if (anti_mode) {
    c("in_defensefinder", "in_dbapis", "in_acrfinder")
  } else {
    c("in_defensepredictor", "in_padloc", "in_defensefinder", "in_rebase")
  }
  support_tool_names <- if (anti_mode) {
    c("AntiDefenseFinder", "dbAPIS", "AcrFinder")
  } else {
    c("DefensePredictor", "PADLOC", "DefenseFinder", "REBASE")
  }
  hit_rows$support_tools <- apply(
    hit_rows[, support_cols, drop = FALSE],
    1,
    function(z) {
      paste(support_tool_names[as.logical(z)], collapse = ";")
    }
  )
  hit_rows$support_count <- rowSums(hit_rows[, support_cols, drop = FALSE], na.rm = TRUE)
  hit_rows$defensepredictor_hit_log_odds <- hit_rows$DefensePredictor_mean_log_odds
  hit_rows$defensepredictor_score_band <- as.character(hit_rows$DefensePredictor_score_band %||% .dnmb_defensepredictor_score_band(hit_rows$DefensePredictor_mean_log_odds))

  hit_rows$representative_source <- ifelse(
    hit_rows$in_defensefinder,
    if (anti_mode) "AntiDefenseFinder" else "DefenseFinder",
    ifelse(hit_rows$in_dbapis, "dbAPIS", ifelse(hit_rows$in_acrfinder, "AcrFinder", ifelse(hit_rows$in_padloc, "PADLOC", ifelse(hit_rows$in_rebase, "REBASE", "DefensePredictor"))))
  )
  df_subtype <- if ("DefenseFinder_system_subtype" %in% names(hit_rows)) as.character(hit_rows$DefenseFinder_system_subtype) else rep(NA_character_, nrow(hit_rows))
  df_type <- if ("DefenseFinder_system_type" %in% names(hit_rows)) as.character(hit_rows$DefenseFinder_system_type) else rep(NA_character_, nrow(hit_rows))
  df_label <- ifelse(!is.na(df_subtype) & nzchar(df_subtype), df_subtype, df_type)
  rb_family <- if ("REBASEfinder_family_id" %in% names(hit_rows)) as.character(hit_rows$REBASEfinder_family_id) else rep(NA_character_, nrow(hit_rows))
  rb_hit <- if ("REBASEfinder_hit_label" %in% names(hit_rows)) as.character(hit_rows$REBASEfinder_hit_label) else rep(NA_character_, nrow(hit_rows))
  rb_label <- ifelse(!is.na(rb_family) & nzchar(rb_family), rb_family, rb_hit)
  hit_rows$representative_system <- paste0("DefensePredictor_only:", hit_rows$defensepredictor_score_band)
  hit_rows$representative_system[hit_rows$in_rebase] <- rb_label[hit_rows$in_rebase]
  hit_rows$representative_system[hit_rows$in_padloc] <- as.character(hit_rows$PADLOC_system[hit_rows$in_padloc])
  if ("AcrFinder_hit_label" %in% names(hit_rows)) {
    acr_lab <- ifelse(
      !is.na(hit_rows$AcrFinder_hit_label) & nzchar(as.character(hit_rows$AcrFinder_hit_label)),
      as.character(hit_rows$AcrFinder_hit_label),
      as.character(hit_rows$AcrFinder_family_id)
    )
    hit_rows$representative_system[hit_rows$in_acrfinder] <- acr_lab[hit_rows$in_acrfinder]
  }
  if ("dbAPIS_hit_label" %in% names(hit_rows)) {
    dbapis_lab <- ifelse(
      !is.na(hit_rows$dbAPIS_hit_label) & nzchar(as.character(hit_rows$dbAPIS_hit_label)),
      as.character(hit_rows$dbAPIS_hit_label),
      as.character(hit_rows$dbAPIS_family_id)
    )
    hit_rows$representative_system[hit_rows$in_dbapis] <- dbapis_lab[hit_rows$in_dbapis]
  }
  hit_rows$representative_system[hit_rows$in_defensefinder] <- df_label[hit_rows$in_defensefinder]

  hit_rows$system_group_source <- NA_character_
  hit_rows$system_group_name <- NA_character_
  hit_rows$system_group_seed <- NA_character_
  hit_rows$system_locus_id <- NA_character_

  df_keep <- hit_rows$in_defensefinder & !is.na(hit_rows$DefenseFinder_system_id) & nzchar(as.character(hit_rows$DefenseFinder_system_id))
  hit_rows$system_group_source[df_keep] <- if (anti_mode) "AntiDefenseFinder" else "DefenseFinder"
  hit_rows$system_group_name[df_keep] <- df_label[df_keep]
  hit_rows$system_group_seed[df_keep] <- as.character(hit_rows$DefenseFinder_system_id[df_keep])
  hit_rows$system_locus_id[df_keep] <- paste0("DF::", hit_rows$DefenseFinder_system_id[df_keep])

  dbapis_keep <- is.na(hit_rows$system_locus_id) & hit_rows$in_dbapis
  hit_rows$system_group_source[dbapis_keep] <- "dbAPIS"
  hit_rows$system_group_name[dbapis_keep] <- ifelse(
    "dbAPIS_hit_label" %in% names(hit_rows) & !is.na(hit_rows$dbAPIS_hit_label[dbapis_keep]) & nzchar(as.character(hit_rows$dbAPIS_hit_label[dbapis_keep])),
    as.character(hit_rows$dbAPIS_hit_label[dbapis_keep]),
    as.character(hit_rows$dbAPIS_family_id[dbapis_keep])
  )
  hit_rows$system_group_seed[dbapis_keep] <- paste(hit_rows$seqid[dbapis_keep], hit_rows$locus_tag[dbapis_keep], hit_rows$dbAPIS_family_id[dbapis_keep], sep = "::")
  hit_rows$system_locus_id[dbapis_keep] <- paste0("DBA::", hit_rows$locus_tag[dbapis_keep])

  acrfinder_keep <- is.na(hit_rows$system_locus_id) & hit_rows$in_acrfinder
  hit_rows$system_group_source[acrfinder_keep] <- "AcrFinder"
  hit_rows$system_group_name[acrfinder_keep] <- ifelse(
    "AcrFinder_hit_label" %in% names(hit_rows) & !is.na(hit_rows$AcrFinder_hit_label[acrfinder_keep]) & nzchar(as.character(hit_rows$AcrFinder_hit_label[acrfinder_keep])),
    as.character(hit_rows$AcrFinder_hit_label[acrfinder_keep]),
    as.character(hit_rows$AcrFinder_family_id[acrfinder_keep])
  )
  hit_rows$system_group_seed[acrfinder_keep] <- paste(hit_rows$seqid[acrfinder_keep], hit_rows$locus_tag[acrfinder_keep], hit_rows$AcrFinder_family_id[acrfinder_keep], sep = "::")
  hit_rows$system_locus_id[acrfinder_keep] <- paste0("ACR::", hit_rows$locus_tag[acrfinder_keep])

  pd_keep <- !df_keep & hit_rows$in_padloc & !is.na(hit_rows$PADLOC_system_number)
  hit_rows$system_group_source[pd_keep] <- "PADLOC"
  hit_rows$system_group_name[pd_keep] <- as.character(hit_rows$PADLOC_system[pd_keep])
  hit_rows$system_group_seed[pd_keep] <- paste(hit_rows$seqid[pd_keep], hit_rows$PADLOC_system_number[pd_keep], sep = "::")
  hit_rows$system_locus_id[pd_keep] <- paste0("PD::", hit_rows$seqid[pd_keep], "::", hit_rows$PADLOC_system_number[pd_keep])

  rb_keep <- is.na(hit_rows$system_locus_id) & hit_rows$in_rebase
  hit_rows$system_group_source[rb_keep] <- "REBASE"
  hit_rows$system_group_name[rb_keep] <- rb_label[rb_keep]
  hit_rows$system_group_seed[rb_keep] <- hit_rows$system_group_name[rb_keep]

  dp_keep <- is.na(hit_rows$system_locus_id) & hit_rows$in_defensepredictor
  hit_rows$system_group_source[dp_keep] <- "DefensePredictor"
  hit_rows$system_group_name[dp_keep] <- ifelse(
    !is.na(hit_rows$annotation_name[dp_keep]) & nzchar(hit_rows$annotation_name[dp_keep]),
    paste0("DP_candidate::", hit_rows$annotation_name[dp_keep]),
    paste0("DefensePredictor_only:", hit_rows$defensepredictor_score_band[dp_keep])
  )
  hit_rows$system_group_seed[dp_keep] <- hit_rows$system_group_name[dp_keep]

  gap_idx <- which(is.na(hit_rows$system_locus_id))
  if (length(gap_idx)) {
    tmp <- hit_rows[gap_idx, , drop = FALSE]
    tmp$start <- suppressWarnings(as.numeric(tmp$start))
    tmp$end <- suppressWarnings(as.numeric(tmp$end))
    tmp$order_index <- seq_len(nrow(tmp))
    tmp <- tmp[order(tmp$system_group_source, tmp$seqid, tmp$start, tmp$end, tmp$locus_tag), , drop = FALSE]
    current_id <- 0L
    prev_source <- prev_seqid <- prev_seed <- NA_character_
    prev_end <- NA_real_
    locus_ids <- character(nrow(tmp))
    for (i in seq_len(nrow(tmp))) {
      same_block <- identical(prev_source, tmp$system_group_source[[i]]) &&
        identical(prev_seqid, tmp$seqid[[i]]) &&
        identical(prev_seed, tmp$system_group_seed[[i]]) &&
        !is.na(prev_end) &&
        !is.na(tmp$start[[i]]) &&
        (tmp$start[[i]] - prev_end) <= as.numeric(region_gap)
      if (!same_block) {
        current_id <- current_id + 1L
      }
      prefix <- ifelse(tmp$system_group_source[[i]] == "REBASE", "RB", "DP")
        locus_ids[[i]] <- sprintf("%s::%03d", prefix, current_id)
        prev_source <- tmp$system_group_source[[i]]
        prev_seqid <- tmp$seqid[[i]]
        prev_seed <- tmp$system_group_seed[[i]]
        prev_end <- tmp$end[[i]]
      }
      tmp$system_locus_id <- locus_ids
      hit_rows$system_locus_id[gap_idx] <- tmp$system_locus_id[order(tmp$order_index)]
  }

  hit_rows <- hit_rows[!is.na(hit_rows$system_locus_id) & nzchar(hit_rows$system_locus_id), , drop = FALSE]
  hit_rows$display_name <- .dnmb_integrated_defense_display_name(hit_rows$system_group_name, source = hit_rows$system_group_source, max_len = 18L)
  hit_rows$display_name_full <- .dnmb_integrated_defense_display_name(hit_rows$system_group_name, source = hit_rows$system_group_source, max_len = 999L)
  hit_rows$display_gene_name <- .dnmb_integrated_defense_gene_name(hit_rows)
  rownames(hit_rows) <- NULL
  list(members = hit_rows, membership = overlap$membership)
}

.dnmb_integrated_defense_overlap_panel <- function(tool_membership) {
  has_eulerr <- requireNamespace("eulerr", quietly = TRUE)
  has_venneuler <- requireNamespace("venneuler", quietly = TRUE)
  if ((!has_eulerr && !has_venneuler) || !requireNamespace("ggplotify", quietly = TRUE)) {
    return(NULL)
  }
  preferred_order <- c("DefensePredictor", "PADLOC", "DefenseFinder", "REBASE", "AntiDefenseFinder", "dbAPIS", "AcrFinder")
  set_colors <- c(
    "DefensePredictor" = "#D81B60",
    "PADLOC" = "#F39C12",
    "DefenseFinder" = "#1F77B4",
    "REBASE" = "#1B9E77",
    "AntiDefenseFinder" = "#6A5ACD",
    "dbAPIS" = "#2CA02C",
    "AcrFinder" = "#B22222"
  )
  if (is.null(tool_membership) || !is.data.frame(tool_membership) || !nrow(tool_membership)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  tool_membership <- as.data.frame(tool_membership, stringsAsFactors = FALSE)
  tool_order <- preferred_order[preferred_order %in% unique(as.character(tool_membership$tool))]
  if (!length(tool_order)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  membership_split <- split(tool_membership$accession, tool_membership$tool)
  all_tools_sets <- lapply(tool_order, function(tool) unique(as.character(membership_split[[tool]])))
  names(all_tools_sets) <- tool_order

  venn_keys <- vapply(
    split(tool_membership$tool, tool_membership$accession),
    function(x) paste(sort(unique(as.character(x))), collapse = "&"),
    character(1)
  )
  venn_tab <- table(venn_keys)
  venn_expr <- stats::setNames(as.numeric(venn_tab), names(venn_tab))
  set_sizes <- vapply(all_tools_sets, function(x) length(unique(x)), numeric(1))
  set_short <- c(
    "DefensePredictor" = "DP",
    "PADLOC" = "PD",
    "DefenseFinder" = "DF",
    "REBASE" = "RB",
    "AntiDefenseFinder" = "ADF",
    "dbAPIS" = "DBA",
    "AcrFinder" = "ACR"
  )
  ve <- NULL
  p_venn <- NULL
  if (has_eulerr) {
    ve <- tryCatch(
      eulerr::euler(venn_expr, shape = "ellipse"),
      error = function(e) NULL
    )
    if (!is.null(ve)) {
      venn_names <- rownames(ve$ellipses)
      venn_cols <- unname(set_colors[venn_names])
      p_venn <- ggplotify::as.ggplot(
        graphics::plot(
          ve,
          labels = FALSE,
          quantities = FALSE,
          fills = list(fill = venn_cols, alpha = 0.32),
          edges = list(col = grDevices::adjustcolor("grey45", alpha.f = 0.45), lwd = 0.7)
        )
      ) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
    }
  }
  if (is.null(ve) && has_venneuler) {
    ve <- tryCatch(
      venneuler::venneuler(venn_expr),
      error = function(e) NULL
    )
  }
  if (is.null(ve)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  if (is.null(p_venn)) {
    ve$labels <- rep("", length(set_sizes))
    ve_local <- ve
    venn_names <- rownames(ve$centers)
    venn_cols <- unname(set_colors[venn_names])
    p_venn <- ggplotify::as.ggplot(function() {
      par(mar = c(0, 0, 0, 0))
      plot(
        ve_local,
        main = NULL,
        col = venn_cols,
        border = grDevices::adjustcolor("grey45", alpha.f = 0.28),
        lwd = 0.5,
        cex = 0.82,
        col.txt = "black",
        alpha = 0.32
      )
    }) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  }

  accessions <- sort(unique(tool_membership$accession))
  incidence <- data.frame(accession = accessions, stringsAsFactors = FALSE)
  for (tool in tool_order) {
    incidence[[tool]] <- accessions %in% unique(all_tools_sets[[tool]])
  }
  subset_list <- list()
  for (k in seq_along(tool_order)) {
    cmb <- utils::combn(tool_order, k, simplify = FALSE)
    subset_list <- c(subset_list, cmb)
  }
  combo_counts <- do.call(rbind, lapply(subset_list, function(members) {
    keep <- rowSums(incidence[, members, drop = FALSE]) == length(members)
    data.frame(
      combo_key = paste(members, collapse = "&"),
      count = sum(keep),
      n_sets = length(members),
      DefensePredictor = "DefensePredictor" %in% members,
      PADLOC = "PADLOC" %in% members,
      DefenseFinder = "DefenseFinder" %in% members,
      REBASE = "REBASE" %in% members,
      stringsAsFactors = FALSE
    )
  }))
  combo_counts$members <- lapply(subset_list, identity)
  combo_counts <- combo_counts[combo_counts$count > 0, , drop = FALSE]
  combo_counts <- combo_counts[order(-combo_counts$count, combo_counts$n_sets, combo_counts$combo_key), , drop = FALSE]
  top_combo <- combo_counts
  top_combo$rank <- seq_len(nrow(top_combo))
  singleton_last_rank <- if (any(top_combo$n_sets == 1L)) max(top_combo$rank[top_combo$n_sets == 1L]) else NA_integer_
  rank_offset <- 0.00

  row_positions <- stats::setNames(1.00 + 0.55 * seq.int(0L, length(tool_order) - 1L), tool_order)
  matrix_df <- do.call(rbind, lapply(seq_len(nrow(top_combo)), function(i) {
    members <- unlist(top_combo$members[[i]])
    data.frame(
      rank = top_combo$rank[[i]] + rank_offset,
      tool = factor(tool_order, levels = rev(tool_order)),
      y = unname(row_positions[tool_order]),
      present = tool_order %in% members,
      stringsAsFactors = FALSE
    )
  }))
  connector_df <- do.call(rbind, lapply(seq_len(nrow(top_combo)), function(i) {
    members <- unlist(top_combo$members[[i]])
    yvals <- unname(row_positions[members])
    if (length(yvals) < 2) return(NULL)
    data.frame(rank = top_combo$rank[[i]] + rank_offset, y_min = min(yvals), y_max = max(yvals), stringsAsFactors = FALSE)
  }))
  if (is.null(connector_df)) {
    connector_df <- data.frame(rank = numeric(), y_min = numeric(), y_max = numeric())
  }

  set_df <- data.frame(tool = factor(tool_order, levels = rev(tool_order)), size = as.numeric(set_sizes[tool_order]), stringsAsFactors = FALSE)
  set_df$y <- unname(row_positions[tool_order])
  max_count <- max(top_combo$count)
  bar_base <- 3.45
  bar_height_available <- 2.85
  bar_scale <- bar_height_available / max_count
  left_scale <- 5.0 / max(set_df$size)

  bar_df <- top_combo
  bar_df$plot_rank <- bar_df$rank + rank_offset
  bar_df$xmin <- bar_df$plot_rank - 0.24
  bar_df$xmax <- bar_df$plot_rank + 0.24
  bar_df$ymin <- bar_base
  bar_df$ymax <- bar_base + bar_df$count * bar_scale
  y_upper <- max(bar_df$ymax, na.rm = TRUE) + 0.32

  set_bar_df <- set_df
  set_bar_df$xmin <- -set_bar_df$size * left_scale
  set_bar_df$xmax <- 0
  set_bar_df$ymin <- set_bar_df$y - 0.18
  set_bar_df$ymax <- set_bar_df$y + 0.18
  set_bar_df$label_x <- -0.05

  p_upset <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = set_bar_df,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax, fill = .data$tool),
      color = NA,
      alpha = 0.95
    ) +
    ggplot2::geom_text(data = set_df, ggplot2::aes(x = 0.10, y = .data$y, label = .data$size), hjust = 0, size = 2.5, color = "#333333") +
    ggplot2::geom_text(data = set_bar_df, ggplot2::aes(x = .data$label_x, y = .data$y, label = .data$tool), hjust = 1, size = 2.75, color = "#333333") +
    ggplot2::geom_rect(data = bar_df, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax), fill = "#6B7280", color = NA) +
    ggplot2::geom_text(data = bar_df, ggplot2::aes(x = .data$plot_rank, y = .data$ymax + 0.10, label = .data$count), size = 2.5, color = "#333333") +
    {
      if (!is.na(singleton_last_rank) && singleton_last_rank < nrow(top_combo)) {
        ggplot2::geom_segment(
          data = data.frame(x = singleton_last_rank + rank_offset + 0.5, y = 0.70, yend = 5.55),
          ggplot2::aes(x = .data$x, xend = .data$x, y = .data$y, yend = .data$yend),
          inherit.aes = FALSE,
          linewidth = 0.55,
          linetype = "22",
          color = "#9CA3AF"
        )
      }
    } +
    ggplot2::geom_segment(data = connector_df, ggplot2::aes(x = .data$rank, xend = .data$rank, y = .data$y_min, yend = .data$y_max), color = "#111827", linewidth = 0.8) +
    ggplot2::geom_point(data = matrix_df, ggplot2::aes(x = .data$rank, y = .data$y), color = "#D1D5DB", size = 2.0) +
    ggplot2::geom_point(data = matrix_df[matrix_df$present, , drop = FALSE], ggplot2::aes(x = .data$rank, y = .data$y, color = .data$tool), size = 2.6) +
    ggplot2::scale_fill_manual(values = set_colors) +
    ggplot2::scale_color_manual(values = set_colors) +
    ggplot2::scale_x_continuous(breaks = top_combo$rank + rank_offset, labels = top_combo$rank, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.55, y_upper), breaks = NULL, expand = c(0, 0)) +
    ggplot2::labs(x = "Inclusive overlap rank", y = NULL) +
    ggplot2::coord_cartesian(
      xlim = c(min(set_bar_df$xmin) - 0.25, max(top_combo$rank + rank_offset) + 0.45),
      ylim = c(0.55, y_upper),
      clip = "off",
      expand = FALSE
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.4),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(0, 0, 0, 14)
    )

  venn_x <- 0.07
  venn_y <- 0.54
  venn_w <- 0.20
  venn_h <- 0.30
  label_tbl <- .dnmb_integrated_defense_euler_label_positions(
    ve = ve,
    set_sizes = set_sizes,
    set_colors = set_colors,
    set_short = set_short,
    venn_x = venn_x,
    venn_y = venn_y,
    venn_w = venn_w,
    venn_h = venn_h
  )

  panel <- cowplot::ggdraw() +
    cowplot::draw_plot(p_upset, x = 0.00, y = 0.00, width = 1.00, height = 0.88) +
    cowplot::draw_plot(p_venn, x = venn_x, y = venn_y, width = venn_w, height = venn_h)
  for (i in seq_len(nrow(label_tbl))) {
    label_half_w <- max(0.028, nchar(label_tbl$label[[i]]) * 0.0055)
    label_half_h <- 0.012
    dx <- label_tbl$x[[i]] - label_tbl$anchor_x[[i]]
    dy <- label_tbl$y[[i]] - label_tbl$anchor_y[[i]]
    shrink <- 1
    if (abs(dx) >= abs(dy) && abs(dx) > 1e-6) {
      shrink <- max(0, 1 - label_half_w / abs(dx))
    } else if (abs(dy) > 1e-6) {
      shrink <- max(0, 1 - label_half_h / abs(dy))
    }
    line_x1 <- label_tbl$anchor_x[[i]] + dx * shrink
    line_y1 <- label_tbl$anchor_y[[i]] + dy * shrink
    panel <- panel +
      cowplot::draw_grob(
        grid::segmentsGrob(
          x0 = grid::unit(label_tbl$anchor_x[[i]], "npc"),
          y0 = grid::unit(label_tbl$anchor_y[[i]], "npc"),
          x1 = grid::unit(line_x1, "npc"),
          y1 = grid::unit(line_y1, "npc"),
          gp = grid::gpar(col = grDevices::adjustcolor(label_tbl$color[[i]], alpha.f = 0.55), lwd = 0.7)
        )
      ) +
      cowplot::draw_label(
        label_tbl$label[[i]],
        x = label_tbl$x[[i]],
        y = label_tbl$y[[i]],
        hjust = 0.5,
        vjust = 0.5,
        size = 4.1,
        color = label_tbl$color[[i]]
      )
  }
  panel
}

.dnmb_integrated_defense_plot_slug <- function(defensefinder_activity = NULL) {
  if (is.null(defensefinder_activity) || identical(as.character(defensefinder_activity)[1], "Defense")) {
    return("Defense")
  }
  if (identical(as.character(defensefinder_activity)[1], "Anti-defense")) {
    return("AntiDefense")
  }
  "Defense"
}

.dnmb_plot_integrated_defense_module <- function(genbank_table,
                                                 output_dir,
                                                 dp_threshold = .dnmb_defensepredictor_default_threshold(),
                                                 defensefinder_activity = "Defense") {
  built <- .dnmb_integrated_defense_build_members(
    genbank_table,
    dp_threshold = dp_threshold,
    defensefinder_activity = defensefinder_activity
  )
  system_members <- built$members
  tool_membership <- built$membership
  if (!is.data.frame(system_members) || !nrow(system_members)) {
    return(NULL)
  }
  plot_slug <- .dnmb_integrated_defense_plot_slug(defensefinder_activity)

  system_summary <- system_members |>
    dplyr::group_by(.data$system_locus_id, .data$system_group_source, .data$system_group_name, .data$display_name, .data$seqid) |>
    dplyr::summarise(
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      genes_count = dplyr::n(),
      support_count = max(.data$support_count, na.rm = TRUE),
      max_log_odds = .dnmb_integrated_defense_max_or_na(.data$defensepredictor_hit_log_odds),
      .groups = "drop"
    )
  system_summary$start <- suppressWarnings(as.numeric(system_summary$start))
  system_summary$end <- suppressWarnings(as.numeric(system_summary$end))
  system_summary$max_log_odds[!is.finite(system_summary$max_log_odds)] <- NA_real_
  system_summary$weight <- ifelse(is.finite(system_summary$max_log_odds), pmax(system_summary$support_count, system_summary$max_log_odds), system_summary$support_count)
  system_summary$wholeness <- 1

  plot_members <- system_members
  plot_members$contig <- plot_members$seqid
  plot_members$direction <- plot_members$strand
  plot_members$DefenseFinder_system_id <- plot_members$system_locus_id
  plot_members$DefenseFinder_system_subtype <- plot_members$display_name
  plot_members$DefenseFinder_system_score <- plot_members$support_count
  plot_members$DefenseFinder_system_wholeness <- 1
  plot_members$DefenseFinder_hit_profile_cov <- pmin(1, pmax(0.25, suppressWarnings(as.numeric(plot_members$support_count)) / 4))
  plot_members$DefenseFinder_hit_seq_cov <- ifelse(
    is.na(plot_members$defensepredictor_hit_log_odds),
    plot_members$DefenseFinder_hit_profile_cov,
    pmin(1, pmax(0.25, suppressWarnings(as.numeric(plot_members$defensepredictor_hit_log_odds)) / 12))
  )
  plot_members$DefenseFinder_gene_name <- plot_members$display_gene_name
  plot_members$context_label <- plot_members$display_name_full
  palette <- .dnmb_integrated_defense_palette(
    values = c(system_summary$display_name, plot_members$DefenseFinder_system_subtype),
    sources = c(system_summary$system_group_source, plot_members$system_group_source)
  )
  fill_rows <- max(1L, ceiling(length(palette) / 10L))

  overview_windows <- system_summary |>
    dplyr::transmute(
      contig = .data$seqid,
      DefenseFinder_system_id = .data$system_locus_id,
      DefenseFinder_system_subtype = .data$display_name,
      start = .data$start,
      end = .data$end,
      genes_count = .data$genes_count,
      weight = .data$weight,
      wholeness = .data$wholeness
    )

  p_inventory <- .dnmb_plot_integrated_defense_inventory(
    system_summary = system_summary,
    palette = palette
  )

  p_context <- .dnmb_plot_integrated_defense_context(
    plot_members,
    output_dir = output_dir,
    defense_palette = palette,
    legend_position = "none"
  ) + ggplot2::labs(title = NULL) + ggplot2::theme(plot.margin = ggplot2::margin(6, 8, 18, 2))

  p_context_legend <- .dnmb_plot_integrated_defense_context(
    plot_members,
    output_dir = output_dir,
    defense_palette = palette,
    legend_position = "bottom"
  ) + ggplot2::labs(title = NULL) + ggplot2::theme(plot.margin = ggplot2::margin(6, 8, 18, 2))

  p_context_meta_legend <- p_context_legend +
    ggplot2::guides(
      fill = "none",
      size = ggplot2::guide_legend(order = 1, nrow = 1, byrow = TRUE),
      shape = ggplot2::guide_legend(order = 2, nrow = 1, byrow = TRUE),
      color = ggplot2::guide_legend(order = 3, nrow = 1, byrow = TRUE)
    )
  p_context_fill_legend <- p_context_legend +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        nrow = fill_rows,
        byrow = TRUE,
        override.aes = list(shape = 22, size = 5, color = NA)
      ),
      size = "none",
      shape = "none",
      color = "none"
    )

  legend_context_meta <- cowplot::get_legend(p_context_meta_legend)
  legend_context_fill <- cowplot::get_legend(p_context_fill_legend)

  top_systems <- system_summary |>
    dplyr::transmute(
      DefenseFinder_system_id = .data$system_locus_id,
      DefenseFinder_system_subtype = .data$display_name,
      contig = .data$seqid,
      weight = .data$weight,
      region_start = .data$start,
      region_end = .data$end
    ) |>
    dplyr::arrange(.data$DefenseFinder_system_subtype, .data$contig, .data$region_start, dplyr::desc(.data$weight), .data$DefenseFinder_system_id)
  top_systems$region_start <- suppressWarnings(as.numeric(top_systems$region_start))
  top_systems$region_end <- suppressWarnings(as.numeric(top_systems$region_end))

  contig_lengths <- .dnmb_contig_lengths_for_plot(plot_members, output_dir = output_dir)
  sector_layout <- .dnmb_sector_layout(contig_lengths)
  p_detail <- .dnmb_defensefinder_radial_detail_plot(
    system_summary = top_systems,
    full_tbl = .dnmb_contig_ordered_table(plot_members),
    sector_layout = sector_layout,
    palette = palette
  )
  p_detail <- p_detail + ggplot2::theme(legend.position = "none")

  p_overlap <- .dnmb_integrated_defense_overlap_panel(tool_membership)
  if (is.null(p_overlap)) {
    p_overlap <- ggplot2::ggplot() + ggplot2::theme_void()
  }

  add_panel_header <- function(plot_obj, label, title, x = 0.0005, y = 0.968, size = 12, plot_y = 0.00, plot_h = 0.90) {
    cowplot::ggdraw() +
      cowplot::draw_plot(plot_obj, x = 0, y = plot_y, width = 1, height = plot_h) +
      cowplot::draw_label(paste0(label, "  ", title), x = x, y = y, hjust = 0, vjust = 1, fontface = "bold", size = size)
  }

  d_context_panel <- cowplot::plot_grid(p_context, ncol = 1, rel_heights = c(1.0), align = "v")
  top_row <- cowplot::plot_grid(
    add_panel_header(p_overlap, "A", "Overlap summary", plot_y = 0.06, plot_h = 0.92),
    add_panel_header(
      p_inventory,
      "B",
      if (identical(plot_slug, "AntiDefense")) "Integrated anti-defense inventory" else "Integrated defense inventory",
      plot_y = 0.06,
      plot_h = 0.84
    ),
    ncol = 2,
    rel_widths = c(1.04, 0.96)
  )

  p_detail_wrapped <- cowplot::ggdraw() +
    cowplot::draw_plot(p_detail, x = 0.005, y = 0.05, width = 0.99, height = 0.92)
  d_context_wrapped <- cowplot::ggdraw() +
    cowplot::draw_plot(d_context_panel, x = 0.01, y = 0.00, width = 0.98, height = 0.98)
  common_legend <- cowplot::ggdraw() +
    cowplot::draw_plot(
      cowplot::plot_grid(
        cowplot::ggdraw() +
          cowplot::draw_grob(legend_context_meta, x = 0.5, y = 0.80, width = 0.86, height = 0.22, hjust = 0.5, vjust = 0.5),
        cowplot::ggdraw() +
          cowplot::draw_grob(legend_context_fill, x = 0.5, y = 0.42, width = 0.995, height = 0.26, hjust = 0.5, vjust = 0.5),
        ncol = 1,
        rel_heights = c(0.40, 0.60)
      ),
      x = 0,
      y = 0,
      width = 1,
      height = 1
    )

  composite <- cowplot::plot_grid(
    top_row,
    add_panel_header(
      p_detail_wrapped,
      "C",
      if (identical(plot_slug, "AntiDefense")) "Integrated anti-defense radial detail" else "Integrated defense radial detail",
      plot_y = 0.01,
      plot_h = 0.90
    ),
    add_panel_header(
      d_context_wrapped,
      "D",
      if (identical(plot_slug, "AntiDefense")) "Integrated anti-defense genome layout" else "Integrated defense genome layout",
      plot_y = 0.11,
      plot_h = 0.82
    ),
    common_legend,
    ncol = 1,
    rel_heights = c(2.05, 4.00, 1.50, 0.22),
    align = "v"
  )

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(plot_slug, "_overview.pdf"))
  .dnmb_module_plot_save(composite, pdf_path, width = 10.5, height = 15.3)
  list(pdf = pdf_path)
}

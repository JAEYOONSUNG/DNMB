.dnmb_defensepredictor_plot_data <- function(genbank_table, output_dir) {
  raw_path <- file.path(
    output_dir,
    "dnmb_module_defensepredictor",
    "defense_predictor_output.csv"
  )

  raw <- NULL
  if (file.exists(raw_path)) {
    raw <- tryCatch(
      utils::read.csv(raw_path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
  }

  required_raw <- c(
    "mean_log_odds", "sd_log_odds", "genomic_accession",
    "start", "end", "locus_tag"
  )
  if (!is.null(raw) && nrow(raw) && all(required_raw %in% names(raw))) {
    tbl <- data.frame(
      contig_id = as.character(raw$genomic_accession),
      start = suppressWarnings(as.numeric(raw$start)),
      end = suppressWarnings(as.numeric(raw$end)),
      locus_tag = as.character(raw$locus_tag),
      gene = if ("symbol" %in% names(raw)) as.character(raw$symbol) else NA_character_,
      product = if ("name" %in% names(raw)) as.character(raw$name) else NA_character_,
      DefensePredictor_mean_log_odds = suppressWarnings(as.numeric(raw$mean_log_odds)),
      DefensePredictor_sd_log_odds = suppressWarnings(as.numeric(raw$sd_log_odds)),
      stringsAsFactors = FALSE
    )
  } else {
    # Older saved runs only contain stringent calls in the combined table.
    # Retain that fallback, while preferring the raw score file whenever it is
    # available so that the histogram and review tiers are not pre-filtered.
    base_tbl <- as.data.frame(genbank_table, stringsAsFactors = FALSE)
    required_base <- c(
      "start", "end", "locus_tag", "DefensePredictor_mean_log_odds"
    )
    if (!nrow(base_tbl) || !all(required_base %in% names(base_tbl))) {
      return(NULL)
    }
    keep <- is.finite(suppressWarnings(as.numeric(
      base_tbl$DefensePredictor_mean_log_odds
    )))
    base_tbl <- base_tbl[keep, , drop = FALSE]
    if (!nrow(base_tbl)) {
      return(NULL)
    }
    contig_number <- if ("contig_number" %in% names(base_tbl)) {
      suppressWarnings(as.integer(base_tbl$contig_number))
    } else {
      as.integer(factor(base_tbl$contig))
    }
    tbl <- data.frame(
      contig_id = sprintf("DNMB_CONTIG_%03d", contig_number),
      start = suppressWarnings(as.numeric(base_tbl$start)),
      end = suppressWarnings(as.numeric(base_tbl$end)),
      locus_tag = as.character(base_tbl$locus_tag),
      gene = if ("gene" %in% names(base_tbl)) as.character(base_tbl$gene) else NA_character_,
      product = if ("product" %in% names(base_tbl)) as.character(base_tbl$product) else NA_character_,
      DefensePredictor_mean_log_odds = suppressWarnings(as.numeric(
        base_tbl$DefensePredictor_mean_log_odds
      )),
      DefensePredictor_sd_log_odds = if ("DefensePredictor_sd_log_odds" %in% names(base_tbl)) {
        suppressWarnings(as.numeric(base_tbl$DefensePredictor_sd_log_odds))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }

  tbl <- tbl[
    is.finite(tbl$DefensePredictor_mean_log_odds) &
      is.finite(tbl$start) & is.finite(tbl$end) &
      !is.na(tbl$contig_id) & nzchar(tbl$contig_id),
    , drop = FALSE
  ]
  if (!nrow(tbl)) {
    return(NULL)
  }
  tbl$midpoint <- (tbl$start + tbl$end) / 2
  tbl$gene[is.na(tbl$gene)] <- ""
  tbl$product[is.na(tbl$product) | !nzchar(tbl$product)] <- "hypothetical protein"

  gbff_records <- .dnmb_parse_gbff_records(.dnmb_find_gbff_for_plot(output_dir))
  raw_order <- suppressWarnings(as.integer(sub(
    ".*?([0-9]+)$", "\\1", unique(tbl$contig_id)
  )))
  contigs <- data.frame(
    contig_id = unique(tbl$contig_id),
    order = raw_order,
    stringsAsFactors = FALSE
  )
  contigs$order[is.na(contigs$order)] <- seq_len(nrow(contigs))[is.na(contigs$order)]
  contigs <- contigs[order(contigs$order, contigs$contig_id), , drop = FALSE]

  if (nrow(gbff_records)) {
    idx <- match(contigs$order, gbff_records$order)
    contigs$length_bp <- suppressWarnings(as.numeric(gbff_records$length_bp[idx]))
    contigs$accession <- as.character(gbff_records$accession[idx])
    contigs$definition <- as.character(gbff_records$definition[idx])
  } else {
    contigs$length_bp <- NA_real_
    contigs$accession <- NA_character_
    contigs$definition <- NA_character_
  }
  max_end <- tapply(tbl$end, tbl$contig_id, max, na.rm = TRUE)
  missing_len <- !is.finite(contigs$length_bp) | contigs$length_bp <= 0
  contigs$length_bp[missing_len] <- unname(max_end[contigs$contig_id[missing_len]])

  plasmid_name <- sub(
    ".*\\bplasmid\\s+([^,.;]+).*$", "\\1", contigs$definition,
    ignore.case = TRUE
  )
  is_plasmid <- grepl("\\bplasmid\\b", contigs$definition, ignore.case = TRUE)
  accession_label <- ifelse(
    !is.na(contigs$accession) & nzchar(contigs$accession),
    paste0(" | ", contigs$accession),
    ""
  )
  contigs$replicon_short <- ifelse(
    is_plasmid,
    paste("Plasmid", plasmid_name),
    ifelse(contigs$order == 1L, "Chromosome", paste("Replicon", contigs$order))
  )
  contigs$sector_label <- paste0(contigs$replicon_short, accession_label)
  contigs$facet_label <- paste0(
    contigs$sector_label, " (", scales::label_comma()(contigs$length_bp), " bp)"
  )

  list(scores = tbl, contigs = contigs)
}

.dnmb_defensepredictor_score_band <- function(score, threshold = 4) {
  score <- suppressWarnings(as.numeric(score))
  factor(
    dplyr::case_when(
      score >= 10 ~ "DP_10plus",
      score >= 8 ~ "DP_8to10",
      score >= 6 ~ "DP_6to8",
      score >= threshold ~ "DP_4to6",
      score >= 2 ~ "DP_2to4",
      score >= 0 ~ "DP_0to2",
      TRUE ~ "DP_below0"
    ),
    levels = c(
      "DP_10plus", "DP_8to10", "DP_6to8", "DP_4to6",
      "DP_2to4", "DP_0to2", "DP_below0"
    )
  )
}

.dnmb_defensepredictor_pack_labels <- function(tbl, contigs,
                                                panel_width_in = 7.8,
                                                font_size_pt = 6.2) {
  if (!nrow(tbl)) {
    tbl$label_tier <- integer()
    tbl$label_x <- numeric()
    return(tbl)
  }

  locus_suffix <- sub("^.*_", "", tbl$locus_tag)
  has_gene <- !is.na(tbl$gene) & nzchar(trimws(tbl$gene))
  tbl$feature_label <- ifelse(
    has_gene,
    paste0(trimws(tbl$gene), " (", locus_suffix, ")"),
    tbl$locus_tag
  )
  tbl$label_tier <- 0L
  tbl$label_x <- tbl$midpoint
  tbl$label_half_width <- 0

  for (contig_id in unique(as.character(tbl$contig_id))) {
    ii <- which(as.character(tbl$contig_id) == contig_id)
    contig_len <- contigs$length_bp[match(contig_id, contigs$contig_id)]
    if (!length(contig_len) || !is.finite(contig_len) || contig_len <= 0) {
      contig_len <- max(tbl$end[ii], na.rm = TRUE)
    }
    label_chars <- vapply(
      strsplit(tbl$feature_label[ii], "\n", fixed = TRUE),
      function(x) max(nchar(x), 1L),
      numeric(1)
    )
    # Arial's average glyph is roughly 0.52 em. Converting that screen width
    # to genomic units lets the greedy interval packing remain stable across
    # replicon sizes and keeps edge labels inside the panel.
    width_fraction <- (
      label_chars * font_size_pt * 0.52 / 72 / panel_width_in
    ) + 0.012
    half_width <- pmin(0.24, width_fraction / 2) * contig_len
    label_x <- pmin(pmax(tbl$midpoint[ii], half_width), contig_len - half_width)
    left <- label_x - half_width
    right <- label_x + half_width

    ord <- order(label_x, -tbl$DefensePredictor_mean_log_odds[ii])
    tier_ends <- numeric()
    assigned <- integer(length(ii))
    for (jj in ord) {
      available <- which(left[jj] >= tier_ends)
      tier <- if (length(available)) available[[1]] else length(tier_ends) + 1L
      if (tier > length(tier_ends)) {
        tier_ends <- c(tier_ends, -Inf)
      }
      assigned[jj] <- tier - 1L
      tier_ends[tier] <- right[jj] + 0.006 * contig_len
    }
    tbl$label_tier[ii] <- assigned
    tbl$label_x[ii] <- label_x
    tbl$label_half_width[ii] <- half_width
  }
  tbl
}

.dnmb_plot_defensepredictor_module <- function(
    genbank_table,
    output_dir,
    threshold = .dnmb_defensepredictor_default_threshold(),
    reference_threshold = 0) {
  add_panel_header <- function(plot_obj, label, title, x = 0.0005, y = 0.975,
                               size = 12, plot_y = 0.035, plot_h = 0.90) {
    cowplot::ggdraw() +
      cowplot::draw_plot(plot_obj, x = 0, y = plot_y, width = 1, height = plot_h) +
      cowplot::draw_label(
        paste0(label, "  ", title), x = x, y = y,
        hjust = 0, vjust = 1, fontface = "bold", size = size,
        fontfamily = .dnmb_plot_font_family()
      )
  }

  plot_data <- .dnmb_defensepredictor_plot_data(genbank_table, output_dir)
  if (is.null(plot_data)) {
    return(NULL)
  }
  all_scores <- plot_data$scores
  contigs <- plot_data$contigs
  threshold <- as.numeric(threshold)[1]
  reference_threshold <- as.numeric(reference_threshold)[1]

  all_scores$score_band_plot <- .dnmb_defensepredictor_score_band(
    all_scores$DefensePredictor_mean_log_odds,
    threshold = threshold
  )
  candidates <- all_scores[
    all_scores$DefensePredictor_mean_log_odds >= reference_threshold,
    , drop = FALSE
  ]
  if (!nrow(candidates)) {
    return(NULL)
  }

  band_levels <- c(
    "DP_10plus", "DP_8to10", "DP_6to8", "DP_4to6",
    "DP_2to4", "DP_0to2"
  )
  band_palette <- c(
    DP_10plus = "#C2185B",
    DP_8to10 = "#D84A8A",
    DP_6to8 = "#E67BB4",
    DP_4to6 = "#F2A7CF",
    DP_2to4 = "#B59ACB",
    DP_0to2 = "#D9D3E3"
  )
  band_labels <- c(
    DP_10plus = "\u226510",
    DP_8to10 = "8\u2013<10",
    DP_6to8 = "6\u2013<8",
    DP_4to6 = "4\u2013<6 (stringent)",
    DP_2to4 = "2\u2013<4 (review)",
    DP_0to2 = "0\u2013<2 (reference only)"
  )

  base_theme <- ggplot2::theme_bw(base_size = 10.5, base_family = .dnmb_plot_font_family()) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(family = .dnmb_plot_font_family())
    )

  max_count <- max(
    graphics::hist(
      all_scores$DefensePredictor_mean_log_odds,
      breaks = 45,
      plot = FALSE
    )$counts,
    na.rm = TRUE
  )
  threshold_labels <- data.frame(
    x = c(reference_threshold, 2, threshold),
    y = c(max_count * 0.98, max_count * 0.78, max_count * 0.98),
    label = c("model-positive 0", "review 2", "stringent 4"),
    hjust = c(1.03, 0.5, -0.03),
    stringsAsFactors = FALSE
  )
  p_hist <- ggplot2::ggplot(
    all_scores,
    ggplot2::aes(x = .data$DefensePredictor_mean_log_odds)
  ) +
    ggplot2::geom_histogram(
      bins = 45, fill = "#D84A8A", color = "white",
      linewidth = 0.2, alpha = 0.9
    ) +
    ggplot2::geom_vline(
      xintercept = c(reference_threshold, 2, threshold),
      linetype = c("dotted", "dotdash", "dashed"),
      linewidth = c(0.55, 0.55, 0.72),
      color = c("#777777", "#7A6294", "#333333")
    ) +
    ggplot2::geom_text(
      data = threshold_labels,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, hjust = .data$hjust),
      inherit.aes = FALSE, vjust = 1, size = 2.55, color = "#333333",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::labs(x = "Mean log-odds", y = "Scored proteins") +
    base_theme +
    ggplot2::theme(legend.position = "none")

  candidate_contigs <- contigs[contigs$contig_id %in% candidates$contig_id, , drop = FALSE]
  candidates$contig_facet <- candidate_contigs$facet_label[
    match(candidates$contig_id, candidate_contigs$contig_id)
  ]
  candidates$contig_facet <- factor(
    candidates$contig_facet,
    levels = candidate_contigs$facet_label
  )
  candidate_contigs$contig_facet <- factor(
    candidate_contigs$facet_label,
    levels = candidate_contigs$facet_label
  )

  plasmid_ids <- candidate_contigs$contig_id[
    grepl("^Plasmid", candidate_contigs$replicon_short)
  ]
  layout_labels <- candidates[
    candidates$DefensePredictor_mean_log_odds >= 2 |
      candidates$contig_id %in% plasmid_ids,
    , drop = FALSE
  ]
  if (!nrow(layout_labels)) {
    top_by_contig <- unlist(lapply(
      split(seq_len(nrow(candidates)), candidates$contig_id),
      function(ii) ii[which.max(candidates$DefensePredictor_mean_log_odds[ii])]
    ), use.names = FALSE)
    layout_labels <- candidates[top_by_contig, , drop = FALSE]
  }
  layout_labels <- .dnmb_defensepredictor_pack_labels(
    layout_labels, candidate_contigs
  )
  label_base_y <- 0.52
  label_step <- 0.16
  layout_labels$label_y <- label_base_y + layout_labels$label_tier * label_step
  max_label_tier <- max(layout_labels$label_tier, na.rm = TRUE)
  layout_ymax <- label_base_y + max_label_tier * label_step + 0.12
  p_layout <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = candidate_contigs,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = 0.18, yend = 0.18),
      linewidth = 1.15, color = "grey80"
    ) +
    ggplot2::geom_segment(
      data = layout_labels,
      ggplot2::aes(
        x = .data$midpoint, xend = .data$label_x,
        y = 0.25, yend = .data$label_y - 0.045,
        group = interaction(.data$contig_id, .data$locus_tag)
      ),
      linewidth = 0.25, color = "grey58", alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = candidates,
      ggplot2::aes(
        x = .data$midpoint, y = 0.18,
        fill = .data$score_band_plot,
        size = pmax(.data$DefensePredictor_mean_log_odds, 0)
      ),
      shape = 21, color = "grey20", stroke = 0.22, alpha = 0.92
    ) +
    ggplot2::geom_text(
      data = layout_labels,
      ggplot2::aes(
        x = .data$label_x, y = .data$label_y,
        label = .data$feature_label
      ),
      size = 2.18, lineheight = 0.88, color = "grey20",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::facet_wrap(
      ~contig_facet, ncol = 1, scales = "free_x", strip.position = "top"
    ) +
    ggplot2::scale_fill_manual(
      values = band_palette, breaks = band_levels, labels = band_labels,
      drop = FALSE, name = "Mean log-odds tier"
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 6.5), guide = "none") +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0.012, 0.012))
    ) +
    ggplot2::coord_cartesian(ylim = c(0.08, layout_ymax), clip = "on") +
    ggplot2::labs(x = "Replicon coordinate (bp)", y = NULL) +
    base_theme +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold", size = 8.5, margin = ggplot2::margin(t = 1, b = 1)
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey97", colour = "grey85", linewidth = 0.3
      ),
      panel.spacing.y = grid::unit(0.10, "lines"),
      legend.position = "none"
    )

  candidates$replicon_short <- candidate_contigs$replicon_short[
    match(candidates$contig_id, candidate_contigs$contig_id)
  ]
  candidates$feature_id <- ifelse(
    !is.na(candidates$gene) & nzchar(trimws(candidates$gene)),
    paste0(trimws(candidates$gene), " | ", candidates$locus_tag),
    candidates$locus_tag
  )
  candidates$annotation <- paste0(
    ifelse(grepl("^Plasmid", candidates$replicon_short),
           paste0("[", candidates$replicon_short, "]  "), ""),
    candidates$product
  )
  candidates <- candidates[
    order(
      -candidates$DefensePredictor_mean_log_odds,
      candidates$contig_id,
      candidates$locus_tag
    ),
    , drop = FALSE
  ]
  candidates$candidate_rank <- seq_len(nrow(candidates))
  candidates$feature_id <- factor(
    candidates$feature_id,
    levels = rev(candidates$feature_id)
  )

  annotation_x <- max(22, ceiling(max(candidates$DefensePredictor_mean_log_odds, na.rm = TRUE)) + 3)
  table_xmax <- annotation_x + 34
  p_top <- ggplot2::ggplot(
    candidates,
    ggplot2::aes(x = .data$DefensePredictor_mean_log_odds, y = .data$feature_id)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$score_band_plot),
      width = 0.64, color = "grey35", linewidth = 0.22
    ) +
    ggplot2::geom_vline(
      xintercept = threshold, linetype = "dashed",
      linewidth = 0.55, color = "#333333"
    ) +
    ggplot2::geom_vline(
      xintercept = 2, linetype = "dotdash",
      linewidth = 0.45, color = "#7A6294"
    ) +
    ggplot2::geom_vline(
      xintercept = annotation_x - 0.8,
      linewidth = 0.35, color = "grey82"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = pmin(.data$DefensePredictor_mean_log_odds + 0.22, annotation_x - 1.25),
        label = ifelse(
          is.finite(.data$DefensePredictor_sd_log_odds),
          paste0("\u00b1", sprintf("%.1f", .data$DefensePredictor_sd_log_odds)),
          ""
        )
      ),
      hjust = 0, size = 1.82, color = "grey30",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = annotation_x, label = .data$annotation),
      hjust = 0, size = 1.84, color = "grey27",
      family = .dnmb_plot_font_family()
    ) +
    ggplot2::scale_fill_manual(
      values = band_palette, breaks = band_levels, labels = band_labels,
      drop = FALSE, name = "Mean log-odds tier"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 2, 4, 8, 12, 16),
      limits = c(0, table_xmax),
      expand = ggplot2::expansion(mult = c(0, 0.008))
    ) +
    ggplot2::labs(
      x = "Mean log-odds  |  call cutoff = 4; scores 0\u2013<4 are visual review only",
      y = NULL
    ) +
    base_theme +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6.2, color = "grey20"),
      axis.text.x = ggplot2::element_text(size = 7),
      axis.title.x = ggplot2::element_text(size = 8.3),
      plot.margin = ggplot2::margin(5, 8, 4, 4)
    )

  legend_plot <- p_layout +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(face = "bold", size = 9.2),
      legend.text = ggplot2::element_text(size = 8.1),
      legend.key.width = grid::unit(1.05, "lines"),
      legend.key.height = grid::unit(0.82, "lines"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))
  common_legend <- .dnmb_with_plot_pdf_device(
    cowplot::get_legend(legend_plot), width = 10, height = 1
  )
  legend_panel <- cowplot::ggdraw() + cowplot::draw_grob(
    common_legend, x = 0.5, y = 0.5, width = 0.98, height = 0.95,
    hjust = 0.5, vjust = 0.5
  )

  hist_h <- 1.55
  layout_h <- max(1.75, 0.58 * nrow(candidate_contigs) + 0.16 * (max_label_tier + 1) + 0.55)
  candidate_h <- max(3.35, 0.135 * nrow(candidates) + 0.90)
  legend_h <- 0.48
  composite <- cowplot::plot_grid(
    add_panel_header(
      p_hist, "A", "DefensePredictor score distribution",
      plot_y = 0.04, plot_h = 0.88
    ),
    add_panel_header(
      p_layout, "B", "Candidate loci across chromosome and plasmids",
      plot_y = 0.03, plot_h = 0.90
    ),
    add_panel_header(
      p_top, "C", "Ranked candidates through the model-positive boundary",
      plot_y = 0.025, plot_h = 0.925
    ),
    legend_panel,
    ncol = 1,
    rel_heights = c(hist_h, layout_h, candidate_h, legend_h)
  )
  total_h <- hist_h + layout_h + candidate_h + legend_h
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "DefensePredictor_overview.pdf")
  .dnmb_module_plot_save(composite, pdf_path, width = 10, height = total_h)
  list(
    pdf = pdf_path,
    stringent_calls = sum(candidates$DefensePredictor_mean_log_odds >= threshold),
    review_candidates = sum(
      candidates$DefensePredictor_mean_log_odds >= 2 &
        candidates$DefensePredictor_mean_log_odds < threshold
    ),
    reference_candidates = sum(
      candidates$DefensePredictor_mean_log_odds >= reference_threshold &
        candidates$DefensePredictor_mean_log_odds < 2
    )
  )
}

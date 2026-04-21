.dnmb_defensefinder_activity_slug <- function(activity) {
  activity <- trimws(as.character(activity)[1])
  if (!nzchar(activity) || is.na(activity)) {
    return("DefenseFinder")
  }
  if (identical(activity, "Anti-defense")) {
    return("AntiDefenseFinder")
  }
  "DefenseFinder"
}

.dnmb_defensefinder_plot_table <- function(genbank_table, activity = NULL) {
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  req <- c("DefenseFinder_system_id", "DefenseFinder_system_subtype")
  if (!base::nrow(tbl) || !base::all(req %in% base::names(tbl))) {
    return(data.frame())
  }
  tbl <- tbl[!is.na(tbl$DefenseFinder_system_id) & base::nzchar(tbl$DefenseFinder_system_id), , drop = FALSE]
  if (!base::nrow(tbl)) {
    return(tbl)
  }
  if (is.null(activity) || !"DefenseFinder_system_activity" %in% base::names(tbl)) {
    return(tbl)
  }
  keep <- !is.na(tbl$DefenseFinder_system_activity) &
    base::as.character(tbl$DefenseFinder_system_activity) == base::as.character(activity)[1]
  tbl[keep, , drop = FALSE]
}

.dnmb_plot_defensefinder_activity_module <- function(genbank_table, output_dir, activity = NULL) {
  tbl <- .dnmb_defensefinder_plot_table(genbank_table, activity = activity)
  if (!base::nrow(tbl)) {
    return(NULL)
  }
  activity_label <- .dnmb_defensefinder_activity_slug(activity)
  contig_lengths <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  overview_windows <- tbl |>
    dplyr::group_by(.data$contig, .data$DefenseFinder_system_id, .data$DefenseFinder_system_subtype) |>
    dplyr::summarise(
      start = min(.data$start),
      end = max(.data$end),
      genes_count = dplyr::n(),
      weight = max(.data$DefenseFinder_system_score, na.rm = TRUE),
      wholeness = max(suppressWarnings(as.numeric(.data$DefenseFinder_system_wholeness)), na.rm = TRUE),
      .groups = "drop"
    )
  defense_palette <- .dnmb_defensefinder_palette(tbl$DefenseFinder_system_subtype)
  p_inventory <- .dnmb_plot_defensefinder_inventory(
    overview_windows,
    defense_palette,
    legend_position = "none"
  ) +
    ggplot2::labs(title = paste0("A   ", activity_label, " inventory"))
  p_context <- .dnmb_plot_defensefinder_context(
    genbank_table,
    output_dir = output_dir,
    defense_palette = defense_palette,
    system_ids = unique(tbl$DefenseFinder_system_id),
    legend_position = "none"
  ) +
    ggplot2::labs(title = paste0("B   ", activity_label, " genome layout"))
  p_context_legend <- .dnmb_plot_defensefinder_context(
    genbank_table,
    output_dir = output_dir,
    defense_palette = defense_palette,
    system_ids = unique(tbl$DefenseFinder_system_id),
    legend_position = "bottom"
  )
  legend_context <- cowplot::get_legend(p_context_legend)
  legend_row <- cowplot::ggdraw() +
    cowplot::draw_grob(legend_context, x = 0.5, y = 0.5, width = 0.92, height = 0.92, hjust = 0.5, vjust = 0.5)

  top_systems <- tbl |>
    dplyr::group_by(.data$DefenseFinder_system_id, .data$DefenseFinder_system_subtype, .data$contig) |>
    dplyr::summarise(
      weight = max(.data$DefenseFinder_system_score, na.rm = TRUE),
      region_start = min(.data$start, na.rm = TRUE),
      region_end = max(.data$end, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$DefenseFinder_system_subtype, .data$contig, .data$region_start, dplyr::desc(.data$weight), .data$DefenseFinder_system_id)
  full_tbl <- .dnmb_contig_ordered_table(genbank_table)
  sector_layout <- .dnmb_sector_layout(contig_lengths)
  p_detail <- .dnmb_defensefinder_radial_detail_plot(
    system_summary = top_systems,
    full_tbl = full_tbl,
    sector_layout = sector_layout,
    palette = defense_palette
  )
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(activity_label, "_overview.pdf"))
  composite <- cowplot::plot_grid(
    p_inventory,
    p_context,
    legend_row,
    p_detail,
    labels = c("", "", "", "C"),
    label_size = 14,
    label_fontface = "bold",
    label_x = 0,
    label_y = c(1, 1, 1, 1),
    hjust = 0,
    ncol = 1,
    align = "v",
    rel_heights = c(0.8, 1.35, 0.25, 4.0)
  )
  .dnmb_module_plot_save(composite, pdf_path, width = 9.5, height = 13)
  list(pdf = pdf_path)
}

.dnmb_plot_defensefinder_module <- function(genbank_table, output_dir) {
  plots <- list()
  defense_plot <- .dnmb_plot_defensefinder_activity_module(
    genbank_table = genbank_table,
    output_dir = output_dir,
    activity = "Defense"
  )
  if (is.list(defense_plot)) {
    plots$defense_pdf <- defense_plot$pdf
  }
  antidefense_plot <- .dnmb_plot_defensefinder_activity_module(
    genbank_table = genbank_table,
    output_dir = output_dir,
    activity = "Anti-defense"
  )
  if (is.list(antidefense_plot)) {
    plots$antidefense_pdf <- antidefense_plot$pdf
  }
  if (!base::length(plots)) {
    return(NULL)
  }
  plots
}

.dnmb_defensefinder_palette <- function(values) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) {
    return(character())
  }
  pal <- grDevices::hcl.colors(length(values), palette = "Dark 3")
  stats::setNames(pal, values)
}

.dnmb_defensefinder_category <- function(subtype) {
  subtype <- as.character(subtype)
  categories <- dplyr::case_when(
    grepl("^RM_", subtype)        ~ "Restriction-modification",
    grepl("^Abi", subtype)        ~ "Abortive infection",
    grepl("^Maz|^Vap|^Hic|^Doc|^Rel|^Par|^Pas|^Hip|^Phd|^Ccd|^YoeB|^MqsR|^GraT|^Toxin|^TA_", subtype, ignore.case = TRUE) ~ "Toxin-antitoxin",
    grepl("^CRISPR|^Cas", subtype) ~ "CRISPR-Cas",
    grepl("^BREX", subtype)       ~ "BREX",
    grepl("^DISARM", subtype)     ~ "DISARM",
    grepl("^Wadjet", subtype)     ~ "Wadjet",
    grepl("^Zorya", subtype)      ~ "Zorya",
    grepl("^Thoeris", subtype)    ~ "Thoeris",
    grepl("^Gabija", subtype)     ~ "Gabija",
    grepl("^Lamassu", subtype)    ~ "Lamassu",
    grepl("^Druantia", subtype)   ~ "Druantia",
    grepl("^Kiwa", subtype)       ~ "Kiwa",
    grepl("^Septu", subtype)      ~ "Septu",
    grepl("^Hachiman", subtype)   ~ "Hachiman",
    grepl("^Shedu", subtype)      ~ "Shedu",
    grepl("^Retron", subtype)     ~ "Retron",
    grepl("^CBASS|^cyclic", subtype, ignore.case = TRUE) ~ "CBASS",
    grepl("^Pycsar", subtype)     ~ "Pycsar",
    grepl("^DND|^Dnd", subtype)   ~ "DND",
    TRUE                          ~ "Other"
  )
  categories
}

.dnmb_plot_defensefinder_inventory <- function(system_summary,
                                              palette,
                                              legend_position = "none") {
  inventory_tbl <- system_summary |>
    dplyr::group_by(.data$DefenseFinder_system_subtype) |>
    dplyr::summarise(
      n_systems = dplyr::n(),
      n_genes = sum(.data$genes_count, na.rm = TRUE),
      best_score = max(.data$weight, na.rm = TRUE),
      mean_wholeness = mean(.data$wholeness, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_systems), dplyr::desc(.data$n_genes), .data$DefenseFinder_system_subtype)
  inventory_tbl$DefenseFinder_system_subtype <- factor(
    inventory_tbl$DefenseFinder_system_subtype,
    levels = rev(inventory_tbl$DefenseFinder_system_subtype)
  )

  ggplot2::ggplot(inventory_tbl, ggplot2::aes(
    x = .data$n_systems,
    y = .data$DefenseFinder_system_subtype
  )) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$DefenseFinder_system_subtype),
      width = 0.6, alpha = 0.9, color = "grey40", linewidth = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 0.05, label = .data$DefenseFinder_system_subtype),
      hjust = 0, size = 2.8, color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$n_systems - 0.05, label = .data$n_systems),
      hjust = 1, size = 3.0, color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$n_systems + 0.08,
        label = paste0(.data$n_genes, " genes | score ",
                       format(round(.data$best_score, 1), nsmall = 1),
                       " | wholeness ",
                       format(round(.data$mean_wholeness, 2), nsmall = 2))
      ),
      hjust = 0, size = 2.5, color = "grey30"
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.52)),
      breaks = seq(0, max(inventory_tbl$n_systems), by = 1)
    ) +
    ggplot2::labs(
      title = "A   Defensefinder inventory",
      x = "Systems detected",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.title.position = "plot",
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 4, 4, 18)
    )
}

.dnmb_plot_defensefinder_context <- function(genbank_table,
                                            output_dir,
                                            defense_palette,
                                            system_ids = NULL,
                                            legend_position = "none") {
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  if (!nrow(tbl)) {
    return(NULL)
  }

  defense_tbl <- tbl[!is.na(tbl$DefenseFinder_system_id) & nzchar(tbl$DefenseFinder_system_id), , drop = FALSE]
  if (!is.null(system_ids)) {
    system_ids <- unique(as.character(system_ids))
    defense_tbl <- defense_tbl[base::as.character(defense_tbl$DefenseFinder_system_id) %in% system_ids, , drop = FALSE]
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
      .groups = "drop"
    )
  if (!nrow(defense_windows)) {
    return(NULL)
  }
  contig_lengths <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  # Only show contigs with defense system hits
  hit_contigs <- unique(defense_windows$contig)
  contig_lengths <- contig_lengths[contig_lengths$contig %in% hit_contigs, , drop = FALSE]
  contig_lengths$track <- 1
  defense_windows$track <- 1
  defense_windows$midpoint <- (defense_windows$start + defense_windows$end) / 2
  defense_windows <- defense_windows[order(defense_windows$contig, defense_windows$midpoint), , drop = FALSE]
  defense_windows$label_tier <- ave(seq_len(nrow(defense_windows)), defense_windows$contig, FUN = function(x) ((seq_along(x) - 1) %% 4))
  defense_windows$label_y <- 1.04 + defense_windows$label_tier * 0.045
  defense_windows$label <- ifelse(
    duplicated(defense_windows$DefenseFinder_system_subtype) | duplicated(defense_windows$DefenseFinder_system_subtype, fromLast = TRUE),
    paste0(defense_windows$DefenseFinder_system_subtype, " (", ave(seq_len(nrow(defense_windows)), defense_windows$DefenseFinder_system_subtype, FUN = seq_along), ")"),
    as.character(defense_windows$DefenseFinder_system_subtype)
  )
  defense_windows$coverage_mean <- rowMeans(
    cbind(defense_windows$profile_cov, defense_windows$seq_cov),
    na.rm = TRUE
  )
  defense_windows$coverage_mean[!is.finite(defense_windows$coverage_mean)] <- NA_real_
  defense_windows$replicon_class <- ifelse(
    grepl("plasmid", defense_windows$contig, ignore.case = TRUE),
    "Plasmid",
    "Chromosome"
  )
  defense_windows$color_value <- unname(defense_palette[match(defense_windows$DefenseFinder_system_subtype, names(defense_palette))])
  contig_map <- contig_lengths$length_bp[match(defense_windows$contig, contig_lengths$contig)]
  offset_bp <- pmax(1200, contig_map * 0.0100)
  coverage_bg_tbl <- rbind(
    data.frame(
      contig = defense_windows$contig,
      x = defense_windows$midpoint - offset_bp / 2,
      xend = defense_windows$midpoint + offset_bp / 2,
      y = 0.915,
      yend = 0.915,
      stringsAsFactors = FALSE
    ),
    data.frame(
      contig = defense_windows$contig,
      x = defense_windows$midpoint - offset_bp / 2,
      xend = defense_windows$midpoint + offset_bp / 2,
      y = 0.875,
      yend = 0.875,
      stringsAsFactors = FALSE
    )
  )
  coverage_tbl <- rbind(
    data.frame(
      contig = defense_windows$contig,
      x = defense_windows$midpoint - offset_bp / 2,
      xend = defense_windows$midpoint - offset_bp / 2 + offset_bp * pmax(0, pmin(1, defense_windows$profile_cov)),
      y = 0.915,
      yend = 0.915,
      metric = "P cov",
      stringsAsFactors = FALSE
    ),
    data.frame(
      contig = defense_windows$contig,
      x = defense_windows$midpoint - offset_bp / 2,
      xend = defense_windows$midpoint - offset_bp / 2 + offset_bp * pmax(0, pmin(1, defense_windows$seq_cov)),
      y = 0.875,
      yend = 0.875,
      metric = "S cov",
      stringsAsFactors = FALSE
    )
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
      ggplot2::aes(
        x = .data$midpoint,
        y = .data$track,
        fill = .data$DefenseFinder_system_subtype,
        size = .data$score,
        shape = .data$replicon_class
      ),
      color = "grey20",
      stroke = 0.35,
      show.legend = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = defense_windows,
      ggplot2::aes(x = .data$midpoint, y = 1.02, label = .data$label),
      size = 2.3,
      color = defense_windows$color_value,
      show.legend = FALSE,
      inherit.aes = FALSE,
      direction = "y", nudge_y = 0.04,
      segment.size = 0.2, segment.color = "grey55",
      max.overlaps = Inf, seed = 42,
      min.segment.length = 0,
      force = 1.5, force_pull = 0.3,
      box.padding = 0.15,
      ylim = c(1.04, 1.21)
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
    ggplot2::scale_fill_manual(values = defense_palette, guide = "none") +
    ggplot2::scale_color_manual(values = c("P cov" = "#2563EB", "S cov" = "#D97706"), name = "Coverage") +
    ggplot2::scale_shape_manual(values = c(Chromosome = 21, Plasmid = 24), name = "Replicon") +
    ggplot2::scale_size_continuous(range = c(4.0, 8.0), breaks = sort(unique(defense_windows$score)), name = "Score") +
    ggplot2::labs(
      title = "B   Defensefinder genome layout",
      subtitle = "Segment = system span | point size = score | lower mini-bars = profile/sequence coverage",
      x = "Genome coordinate (bp)",
      y = NULL
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0.02, 0.08)),
      breaks = function(x) {
        # Match the closest contig length to this facet's axis range
        data_max <- contig_lengths$length_bp[which.min(abs(contig_lengths$length_bp - x[2]))]
        b <- scales::breaks_extended()(c(0, data_max))
        thresh <- diff(range(x)) * 0.07
        b <- b[abs(b - data_max) > thresh]
        sort(unique(c(b, data_max)))
      }
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(4, 4, 4, 18),
      legend.position = legend_position,
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.84, 1.22),
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    ggplot2::guides(
      fill = "none",
      size = ggplot2::guide_legend(order = 1, nrow = 1, byrow = TRUE),
      shape = ggplot2::guide_legend(order = 2, nrow = 1, byrow = TRUE),
      color = ggplot2::guide_legend(order = 3, nrow = 1, byrow = TRUE)
    )
}


.dnmb_defensefinder_radial_detail_plot <- function(system_summary,
                                                   full_tbl,
                                                   sector_layout,
                                                   palette) {
  system_summary <- as.data.frame(system_summary, stringsAsFactors = FALSE)
  full_tbl <- as.data.frame(full_tbl, stringsAsFactors = FALSE)
  if (!nrow(system_summary) || !nrow(sector_layout)) {
    return(NULL)
  }

  center_x <- 0
  center_y <- 0
  backbone_mid <- 1.42
  backbone_half <- 0.09
  panel_ring_start <- 3.25
  panel_ring_step <- 2.80
  panel_span <- 0.84
  panel_body <- 0.42
  panel_head <- 0.75

  angle_dist <- function(a, b) {
    abs(atan2(sin(a - b), cos(a - b)))
  }
  ring_conflict_count <- function(layout_df, gap_rad) {
    if (!nrow(layout_df)) {
      return(0L)
    }
    total <- 0L
    for (ring in unique(layout_df$ring_id)) {
      sub <- layout_df[layout_df$ring_id == ring, , drop = FALSE]
      if (nrow(sub) < 2L) {
        next
      }
      for (i in seq_len(nrow(sub) - 1L)) {
        for (j in seq.int(i + 1L, nrow(sub))) {
          overlap <- angle_dist(sub$display_angle[[i]], sub$display_angle[[j]]) -
            ((sub$span_rad[[i]] + sub$span_rad[[j]]) / 2 + gap_rad)
          if (is.na(overlap) || overlap < 0) {
            total <- total + 1L
          }
        }
      }
    }
    total
  }

  backbone_bg <- dplyr::bind_rows(lapply(seq_len(nrow(sector_layout)), function(i) {
    .dnmb_arc_band_polygon(
      xmin = 0,
      xmax = sector_layout$length_bp[[i]],
      start_angle = sector_layout$theta_start[[i]],
      end_angle = sector_layout$theta_end[[i]],
      r_inner = backbone_mid - backbone_half,
      r_outer = backbone_mid + backbone_half,
      center_x = center_x,
      center_y = center_y,
      n = 220L
    )
  }))

  highlight_arcs <- dplyr::bind_rows(lapply(seq_len(nrow(system_summary)), function(i) {
    subtype <- as.character(system_summary$DefenseFinder_system_subtype[[i]])
    fill_value <- palette[[subtype]]
    if (is.null(fill_value) || is.na(fill_value)) {
      fill_value <- "#F87171"
    }
    contig_idx <- match(system_summary$contig[[i]], sector_layout$contig)
    path <- .dnmb_arc_path_segments(
      start = system_summary$region_start[[i]],
      end = system_summary$region_end[[i]],
      xmin = 0,
      xmax = sector_layout$length_bp[[contig_idx]],
      start_angle = sector_layout$theta_start[[contig_idx]],
      end_angle = sector_layout$theta_end[[contig_idx]],
      radius = backbone_mid,
      center_x = center_x,
      center_y = center_y,
      n = 90L
    )
    path$fill_value <- fill_value
    path$system_id <- as.character(system_summary$DefenseFinder_system_id[[i]])
    path
  }))

  system_summary$angle_mid <- .dnmb_sector_position_angle(
    contig = system_summary$contig,
    position = (system_summary$region_start + system_summary$region_end) / 2,
    sector_layout = sector_layout
  )
  system_summary <- system_summary[order(system_summary$angle_mid, decreasing = TRUE), , drop = FALSE]
  system_summary$subtype_n <- ave(system_summary$DefenseFinder_system_subtype, system_summary$DefenseFinder_system_subtype, FUN = length)
  system_summary$subtype_idx <- ave(seq_len(nrow(system_summary)), system_summary$DefenseFinder_system_subtype, FUN = seq_along)
  system_summary$display_label <- ifelse(
    system_summary$subtype_n > 1,
    paste0(system_summary$DefenseFinder_system_subtype, " (", system_summary$subtype_idx, ")"),
    as.character(system_summary$DefenseFinder_system_subtype)
  )
  system_summary$contig_len <- sector_layout$length_bp[match(system_summary$contig, sector_layout$contig)]
  system_summary$zoom_pad <- pmax((system_summary$region_end - system_summary$region_start) * 0.85, 4500)
  system_summary$zoom_start <- pmax(0, system_summary$region_start - system_summary$zoom_pad)
  system_summary$zoom_end <- pmin(system_summary$contig_len, system_summary$region_end + system_summary$zoom_pad)
  system_summary$zoom_width <- system_summary$zoom_end - system_summary$zoom_start
  bp_per_radian <- stats::median(system_summary$zoom_width, na.rm = TRUE) / 0.84
  system_summary$panel_span <- system_summary$zoom_width / bp_per_radian
  gap_rad <- 0.09
  panel_layout <- NULL
  panel_rings <- NULL
  for (n_rings in seq.int(4L, 7L)) {
    candidate_rings <- seq(panel_ring_start, by = panel_ring_step, length.out = n_rings)
    candidate_layout <- .dnmb_assign_panel_layout(
      system_summary$angle_mid,
      span_rad = system_summary$panel_span,
      ring_radii = candidate_rings,
      gap_rad = gap_rad
    )
    if (ring_conflict_count(candidate_layout, gap_rad = gap_rad) == 0L) {
      panel_layout <- candidate_layout
      panel_rings <- candidate_rings
      break
    }
    panel_layout <- candidate_layout
    panel_rings <- candidate_rings
  }
  system_summary$display_angle <- panel_layout$display_angle
  system_summary$ring_id <- panel_layout$ring_id
  system_summary$r_mid <- panel_layout$r_mid

  panel_gene_list <- list()
  panel_bg_list <- list()
  stem_list <- list()
  cap_list <- list()
  title_list <- list()
  endpoint_list <- list()
  gene_label_list <- list()

  for (i in seq_len(nrow(system_summary))) {
    one <- system_summary[i, , drop = FALSE]
    contig_id <- as.character(one$contig[[1]])
    contig_len <- sector_layout$length_bp[match(contig_id, sector_layout$contig)]
    subtype <- as.character(one$DefenseFinder_system_subtype[[1]])
    subtype_color <- palette[[subtype]]
    if (is.null(subtype_color) || is.na(subtype_color)) {
      subtype_color <- "#F87171"
    }
    region_start <- suppressWarnings(as.numeric(one$region_start[[1]]))
    region_end <- suppressWarnings(as.numeric(one$region_end[[1]]))
    zoom_start <- suppressWarnings(as.numeric(one$zoom_start[[1]]))
    zoom_end <- suppressWarnings(as.numeric(one$zoom_end[[1]]))
    theta_mid <- suppressWarnings(as.numeric(one$display_angle[[1]]))
    panel_span_i <- suppressWarnings(as.numeric(one$panel_span[[1]]))
    theta_start <- theta_mid + panel_span_i / 2
    theta_end <- theta_mid - panel_span_i / 2
    r_mid <- suppressWarnings(as.numeric(one$r_mid[[1]]))

    panel_bg_list[[length(panel_bg_list) + 1L]] <- transform(
      .dnmb_arc_band_polygon(
        xmin = zoom_start,
        xmax = zoom_end,
        start_angle = theta_start,
        end_angle = theta_end,
        r_inner = r_mid - panel_head / 2,
        r_outer = r_mid + panel_head / 2,
        center_x = center_x,
        center_y = center_y,
        n = 120L
      ),
      group_id = one$DefenseFinder_system_id[[1]],
      border_value = subtype_color
    )

    zoom_tbl <- full_tbl[
      full_tbl$contig == contig_id &
        !is.na(full_tbl$start) &
        !is.na(full_tbl$end) &
        full_tbl$end >= zoom_start &
        full_tbl$start <= zoom_end,
      ,
      drop = FALSE
    ]
    zoom_tbl$defense_role <- "neighbor"
    sys_loci <- full_tbl$locus_tag[full_tbl$DefenseFinder_system_id == one$DefenseFinder_system_id[[1]]]
    match_idx <- match(zoom_tbl$locus_tag, sys_loci)
    zoom_tbl$defense_role[!is.na(match_idx)] <- subtype

    if ("DefenseFinder_gene_name" %in% names(zoom_tbl)) {
      zoom_tbl$label <- NA_character_
      sys_tbl <- full_tbl[full_tbl$DefenseFinder_system_id == one$DefenseFinder_system_id[[1]], , drop = FALSE]
      gene_idx <- match(zoom_tbl$locus_tag, sys_tbl$locus_tag)
      ok <- !is.na(gene_idx)
      short_locus <- gsub("^.*?(RS[0-9]+)$", "\\1", zoom_tbl$locus_tag[ok])
      zoom_tbl$label[ok] <- paste0(
        sub("^.*__", "", as.character(sys_tbl$DefenseFinder_gene_name[gene_idx[ok]])),
        "\n",
        short_locus
      )
    } else {
      zoom_tbl$label <- NA_character_
    }

    for (j in seq_len(nrow(zoom_tbl))) {
      poly <- .dnmb_arc_gene_polygon(
        start = zoom_tbl$start[[j]],
        end = zoom_tbl$end[[j]],
        direction = zoom_tbl$direction[[j]],
        xmin = zoom_start,
        xmax = zoom_end,
        start_angle = theta_start,
        end_angle = theta_end,
        r_mid = r_mid,
        body_thickness = panel_body,
        head_thickness = panel_head,
        center_x = center_x,
        center_y = center_y
      )
      if (!nrow(poly)) {
        next
      }
      poly$feature_id <- paste(one$DefenseFinder_system_id[[1]], j, sep = "::")
      poly$fill_value <- if (zoom_tbl$defense_role[[j]] == "neighbor") "neighbor" else subtype
      panel_gene_list[[length(panel_gene_list) + 1L]] <- poly

      if (!is.na(zoom_tbl$label[[j]]) && nzchar(zoom_tbl$label[[j]])) {
        mid_ang <- .dnmb_arc_map_angle((zoom_tbl$start[[j]] + zoom_tbl$end[[j]]) / 2, zoom_start, zoom_end, theta_start, theta_end)
        gene_label_list[[length(gene_label_list) + 1L]] <- data.frame(
          system_id = one$DefenseFinder_system_id[[1]],
          base_angle = mid_ang,
          lower_bound = min(theta_start, theta_end) + panel_span_i * 0.06,
          upper_bound = max(theta_start, theta_end) - panel_span_i * 0.06,
          r_mid = r_mid,
          label = zoom_tbl$label[[j]],
          stringsAsFactors = FALSE
        )
      }
    }

    actual_mid <- .dnmb_sector_position_angle(contig_id, (region_start + region_end) / 2, sector_layout)
    neck_radius <- backbone_mid + backbone_half + 0.015
    cap_radius <- r_mid
    stem_top <- .dnmb_arc_xy(theta_mid, cap_radius, center_x = center_x, center_y = center_y)
    stem_bottom <- .dnmb_arc_xy(actual_mid, neck_radius, center_x = center_x, center_y = center_y)
    stem_list[[length(stem_list) + 1L]] <- data.frame(
      x = stem_bottom$x,
      y = stem_bottom$y,
      xend = stem_top$x,
      yend = stem_top$y,
      group_id = one$DefenseFinder_system_id[[1]],
      color_value = subtype_color,
      stringsAsFactors = FALSE
    )
    cap_list[[length(cap_list) + 1L]] <- transform(
      .dnmb_arc_angle_path(
        angle_start = theta_mid + panel_span_i * 0.48,
        angle_end = theta_mid - panel_span_i * 0.48,
        radius = cap_radius,
        center_x = center_x,
        center_y = center_y,
        n = 90L
      ),
      group_id = one$DefenseFinder_system_id[[1]],
      color_value = subtype_color
    )

    title_list[[length(title_list) + 1L]] <- data.frame(
      base_angle = theta_mid,
      lower_bound = min(theta_start, theta_end) + panel_span_i * 0.12,
      upper_bound = max(theta_start, theta_end) - panel_span_i * 0.12,
      title_radius = r_mid + panel_head / 2 + 0.42 + 0.07 * ((one$ring_id[[1]] - 1) %% 2),
      label = one$display_label[[1]],
      stringsAsFactors = FALSE
    )

    endpoint_offset <- min(0.08, panel_span_i * 0.09)
    endpoint_theta <- c(theta_start + endpoint_offset, theta_end - endpoint_offset)
    endpoint_xy <- .dnmb_arc_xy(endpoint_theta, r_mid - panel_head / 2 - 0.18, center_x = center_x, center_y = center_y)
    endpoint_facing <- .dnmb_label_radial_facing(endpoint_theta * 180 / pi)
    endpoint_list[[length(endpoint_list) + 1L]] <- data.frame(
      x = endpoint_xy$x,
      y = endpoint_xy$y,
      label = c(.dnmb_fmt_bp_exact(zoom_start), .dnmb_fmt_bp_exact(zoom_end)),
      angle = endpoint_facing$angle,
      hjust = endpoint_facing$hjust,
      stringsAsFactors = FALSE
    )
  }

  panel_gene_tbl <- dplyr::bind_rows(panel_gene_list)
  panel_bg_tbl <- dplyr::bind_rows(panel_bg_list)
  stem_tbl <- dplyr::bind_rows(stem_list)
  cap_tbl <- dplyr::bind_rows(cap_list)
  title_tbl <- dplyr::bind_rows(title_list)
  endpoint_tbl <- dplyr::bind_rows(endpoint_list)
  gene_label_tbl <- dplyr::bind_rows(gene_label_list)
  if (nrow(title_tbl)) {
    title_tbl$adjusted_angle <- .dnmb_repel_label_angles_bounded(
      angle_rad = title_tbl$base_angle,
      label_text = title_tbl$label,
      lower_bound = title_tbl$lower_bound,
      upper_bound = title_tbl$upper_bound,
      base_sep = 0.078
    )
    title_xy <- .dnmb_arc_xy(title_tbl$adjusted_angle, title_tbl$title_radius, center_x = center_x, center_y = center_y)
    title_facing <- .dnmb_label_angle_facing(title_tbl$adjusted_angle * 180 / pi - 90)
    title_tbl$x <- title_xy$x
    title_tbl$y <- title_xy$y
    title_tbl$angle <- title_facing$angle
    title_tbl$hjust <- title_facing$hjust
  }
  if (nrow(gene_label_tbl)) {
    gene_label_tbl <- gene_label_tbl |>
      dplyr::group_by(.data$system_id) |>
      dplyr::arrange(.data$base_angle, .by_group = TRUE) |>
      dplyr::mutate(
        adjusted_angle = .dnmb_repel_label_angles_bounded(
          angle_rad = .data$base_angle,
          label_text = .data$label,
          lower_bound = .data$lower_bound,
          upper_bound = .data$upper_bound,
          base_sep = 0.085
        ),
        label_radius = .data$r_mid - panel_head / 2 - 0.22,
        x = center_x + .data$label_radius * cos(.data$adjusted_angle),
        y = center_y + .data$label_radius * sin(.data$adjusted_angle)
      ) |>
      dplyr::ungroup()
    gene_facing <- .dnmb_label_radial_facing(gene_label_tbl$adjusted_angle * 180 / pi)
    gene_label_tbl$angle <- gene_facing$angle
    gene_label_tbl$hjust <- ifelse(cos(gene_label_tbl$adjusted_angle) >= 0, 1, 0)
  }
  scale_candidates <- c(500, 1000, 2000, 5000, 10000)
  target_scale <- stats::median(system_summary$zoom_width, na.rm = TRUE) * 0.22
  scale_bp <- max(scale_candidates[scale_candidates <= target_scale], na.rm = TRUE)
  if (!is.finite(scale_bp)) {
    scale_bp <- min(scale_candidates)
  }
  scale_span <- scale_bp / bp_per_radian
  max_ring <- max(panel_rings, na.rm = TRUE)
  # Draw scale bar as horizontal line, positioned above Type legend
  sb_arc_len <- scale_span * stats::median(panel_rings)
  sb_x <- max_ring - 1.3
  sb_y <- -(max_ring - 1.8)
  scale_bar <- data.frame(x = c(sb_x - sb_arc_len / 2, sb_x + sb_arc_len / 2), y = c(sb_y, sb_y))
  scale_label_xy <- list(x = sb_x, y = sb_y + 0.72)
  scale_label <- paste0(format(scale_bp, big.mark = ","), " bp")

  legend_breaks <- names(palette)
  fill_values <- c("neighbor" = "grey96", palette)

  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = backbone_bg, ggplot2::aes(x = .data$x, y = .data$y), fill = "grey97", color = "grey78", linewidth = 0.25) +
    ggplot2::geom_path(data = highlight_arcs, ggplot2::aes(x = .data$x, y = .data$y, group = .data$system_id), color = highlight_arcs$fill_value, linewidth = 1.55, lineend = "round") +
    ggplot2::geom_segment(
      data = stem_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = grDevices::adjustcolor(stem_tbl$color_value, alpha.f = 0.5),
      linewidth = 0.58,
      lineend = "round"
    ) +
    ggplot2::geom_polygon(
      data = panel_bg_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group_id),
      fill = grDevices::adjustcolor("#F8FAFC", alpha.f = 0.7),
      color = grDevices::adjustcolor(panel_bg_tbl$border_value, alpha.f = 0.55),
      linewidth = 0.24,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_path(
      data = cap_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group_id),
      color = cap_tbl$color_value,
      linewidth = 0.72,
      lineend = "round"
    ) +
    ggplot2::geom_polygon(data = panel_gene_tbl, ggplot2::aes(x = .data$x, y = .data$y, group = .data$feature_id, fill = .data$fill_value), color = "grey25", linewidth = 0.16, linejoin = "mitre") +
    ggplot2::geom_text(data = title_tbl, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$angle, hjust = .data$hjust), size = 2.8, fontface = "bold") +
    ggplot2::geom_text(data = endpoint_tbl, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$angle, hjust = .data$hjust), size = 1.75, color = "grey25", inherit.aes = FALSE) +
    ggplot2::geom_text(data = gene_label_tbl, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, angle = .data$angle, hjust = .data$hjust), size = 1.55, lineheight = 0.80, vjust = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = scale_bar$x[1], xend = scale_bar$x[2], y = scale_bar$y[1], yend = scale_bar$y[2]), linewidth = 1.15, color = "grey25", lineend = "round", inherit.aes = FALSE) +
    ggplot2::geom_text(ggplot2::aes(x = scale_label_xy$x, y = scale_label_xy$y, label = scale_label), size = 2.8, color = "grey20", inherit.aes = FALSE) +
    ggplot2::coord_equal(
      xlim = c(-(max_ring + 2.4), max_ring + 2.6),
      ylim = c(-(max_ring + 2.0), max_ring + 0.2),
      clip = "off",
      expand = FALSE
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 18),
      legend.position = c(0.88, 0.03),
      legend.justification = c(1, 0),
      legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.70), color = NA),
      legend.key.height = grid::unit(0.32, "cm"),
      legend.key.width = grid::unit(0.42, "cm"),
      legend.spacing.y = grid::unit(0.02, "cm"),
      legend.title = ggplot2::element_text(size = 9, face = "bold"),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::scale_fill_manual(
      values = fill_values,
      breaks = legend_breaks,
      name = "Defense subtype"
    )
}

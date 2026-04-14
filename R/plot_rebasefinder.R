.dnmb_plot_rebasefinder_overview <- function(genbank_table, output_dir) {
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  req <- "REBASEfinder_family_id"
  if (!nrow(tbl) || !all(req %in% names(tbl))) return(NULL)

  tbl <- tbl[!is.na(tbl$REBASEfinder_family_id) & nzchar(tbl$REBASEfinder_family_id), , drop = FALSE]
  if (!nrow(tbl)) return(NULL)

  # Stage cache may store numeric columns as character — cast upfront
  for (col in c("REBASEfinder_blast_identity", "REBASEfinder_blast_bitscore",
                "REBASEfinder_blast_length", "start", "end")) {
    if (col %in% names(tbl)) tbl[[col]] <- suppressWarnings(as.numeric(tbl[[col]]))
  }
  if ("REBASEfinder_typing_eligible" %in% names(tbl)) {
    tbl$REBASEfinder_typing_eligible <- as.logical(tbl$REBASEfinder_typing_eligible)
  }
  # Infer enzyme_role from hit name when missing or inconsistent
  if ("REBASEfinder_hit_label" %in% names(tbl) && "REBASEfinder_enzyme_role" %in% names(tbl)) {
    hit <- as.character(tbl$REBASEfinder_hit_label)
    inferred <- ifelse(grepl("^M[0-9]?\\.", hit), "M",
                ifelse(grepl("^R[0-9]?\\.", hit), "R",
                ifelse(grepl("^S[0-9]?\\.", hit), "S", NA_character_)))
    fix <- !is.na(inferred) & (is.na(tbl$REBASEfinder_enzyme_role) |
           tbl$REBASEfinder_enzyme_role != inferred)
    tbl$REBASEfinder_enzyme_role[fix] <- inferred[fix]
  }

  rm_palette <- .dnmb_rebasefinder_palette(tbl$REBASEfinder_family_id)
  role_palette <- .dnmb_rebasefinder_role_palette(tbl)

  # Motif verification
  motif_verified <- .dnmb_rebasefinder_motif_verified(tbl)

  # Confidence filter
  has_eligible <- "REBASEfinder_typing_eligible" %in% names(tbl)
  has_ev <- "REBASEfinder_evidence_mode" %in% names(tbl)
  is_hc <- if (has_eligible) {
    !is.na(tbl$REBASEfinder_typing_eligible) & tbl$REBASEfinder_typing_eligible == TRUE
  } else if (has_ev) {
    !is.na(tbl$REBASEfinder_evidence_mode) & tbl$REBASEfinder_evidence_mode == "high_confidence"
  } else {
    rep(TRUE, nrow(tbl))
  }

  n_hc <- sum(is_hc)
  tbl_verified <- tbl[is_hc & motif_verified, , drop = FALSE]

  # If no genes pass the motif+confidence filter, fall back to showing
  # ALL high-confidence hits so the overview is never empty when there
  # ARE REBASE matches. The user expects to see whatever data exists.
  tbl_show <- if (nrow(tbl_verified) > 0) tbl_verified else tbl[is_hc, , drop = FALSE]
  if (!nrow(tbl_show)) tbl_show <- tbl

  p_inventory <- .dnmb_plot_rebasefinder_inventory(
    tbl_show, rm_palette,
    n_total = nrow(tbl), n_hc = n_hc, n_verified = nrow(tbl_verified)
  )

  p_context <- .dnmb_plot_rebasefinder_context(
    genbank_table, output_dir, rm_palette, legend_position = "none"
  )
  p_context_leg <- .dnmb_plot_rebasefinder_context(
    genbank_table, output_dir, rm_palette, legend_position = "bottom"
  )
  legend_context <- cowplot::get_legend(p_context_leg)
  legend_row <- cowplot::ggdraw() +
    cowplot::draw_grob(legend_context, x = 0.5, y = 0.5,
                       width = 0.92, height = 0.92, hjust = 0.5, vjust = 0.5)

  # Shared display labels with operon grouping
  display_info <- .dnmb_rebasefinder_display_labels(tbl)

  # Each subplot wrapped in tryCatch — if any fails, use a placeholder
  # so the composite never crashes entirely.
  empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  p_blast <- tryCatch(
    .dnmb_plot_rebasefinder_blast_quality(tbl, role_palette, display_info, rm_palette),
    error = function(e) empty_plot
  )
  uniprot_doms <- tryCatch(
    .dnmb_rebasefinder_uniprot_domains(tbl, output_dir),
    error = function(e) NULL
  )
  p_domain <- tryCatch(
    .dnmb_plot_rebasefinder_domain_map(tbl, display_info, uniprot_doms),
    error = function(e) empty_plot
  )
  p_motif <- tryCatch(
    .dnmb_plot_rebasefinder_motif_verification(tbl, display_info),
    error = function(e) empty_plot
  )

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "REBASE_overview.pdf")

  n_hits <- nrow(display_info)
  cd_inch <- max(2.0, 0.30 * n_hits + 0.6)

  n_rm_types <- length(unique(tbl_show$REBASEfinder_family_id))
  a_inch <- max(0.7, 0.45 * n_rm_types + 0.5)

  # Build bottom row (C+D+E) with tryCatch
  bottom_row <- tryCatch(
    cowplot::plot_grid(
      p_blast, p_domain, p_motif,
      ncol = 3, rel_widths = c(1.8, 0.6, 0.8),
      align = "h", axis = "tb"
    ),
    error = function(e) NULL
  )

  # Full composite: A + B + legend + C/D/E (or A+B if bottom fails)
  composite <- tryCatch({
    if (!is.null(bottom_row)) {
      cowplot::plot_grid(
        p_inventory, p_context, legend_row, bottom_row,
        labels = c("A", "B", "", "C"),
        label_size = 14, label_fontface = "bold",
        label_x = 0, label_y = c(1.02, 1.02, 1, 1.02),
        hjust = 0, ncol = 1,
        rel_heights = c(a_inch, 2.0, 0.35, cd_inch)
      )
    } else {
      cowplot::plot_grid(
        p_inventory, p_context, legend_row,
        labels = c("A", "B", ""),
        label_size = 14, label_fontface = "bold",
        label_x = 0, hjust = 0, ncol = 1,
        rel_heights = c(a_inch, 2.0, 0.35)
      )
    }
  }, error = function(e) {
    cowplot::plot_grid(
      p_inventory, p_context, legend_row,
      labels = c("A", "B", ""),
      label_size = 14, label_fontface = "bold",
      label_x = 0, hjust = 0, ncol = 1,
      rel_heights = c(a_inch, 2.0, 0.35)
    )
  })
  total_height <- a_inch + 2.0 + 0.35 + cd_inch
  .dnmb_module_plot_save(composite, pdf_path, width = 12, height = min(18, max(9, total_height)))
  list(pdf = pdf_path)
}


# ====================================================================
# Palettes
# ====================================================================
.dnmb_rebasefinder_palette <- function(values) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return(character())
  pal <- grDevices::hcl.colors(max(length(values), 3), palette = "Dark 3")
  stats::setNames(pal[seq_along(values)], values)
}

.dnmb_rebasefinder_role_palette <- function(tbl) {
  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  if (!has_role) return(c(Other = "grey50"))
  roles <- unique(as.character(stats::na.omit(tbl$REBASEfinder_enzyme_role)))
  roles <- roles[nzchar(roles)]
  if (!length(roles)) return(c(Other = "grey50"))
  pal <- grDevices::hcl.colors(max(length(roles), 3), palette = "Dark 3")
  out <- stats::setNames(pal[seq_along(roles)], roles)
  c(out, Other = "grey60")
}

# Use full contig name (no truncation)
.dnmb_rebasefinder_short_contig <- function(x) {
  as.character(x)
}


# ====================================================================
# Shared display label builder — operon grouping + C/D y-axis sync
# ====================================================================
.dnmb_rebasefinder_display_labels <- function(tbl) {
  has_hit      <- "REBASEfinder_hit_label" %in% names(tbl)
  has_rec      <- "REBASEfinder_rec_seq" %in% names(tbl)
  has_eligible <- "REBASEfinder_typing_eligible" %in% names(tbl)
  has_identity <- "REBASEfinder_blast_identity" %in% names(tbl)
  has_operon   <- "REBASEfinder_operon_id" %in% names(tbl)
  has_blastlen <- "REBASEfinder_blast_length" %in% names(tbl)
  has_translation <- "translation" %in% names(tbl)

  hit_label <- if (has_hit) as.character(tbl$REBASEfinder_hit_label) else tbl$locus_tag
  rec_seq   <- if (has_rec) as.character(tbl$REBASEfinder_rec_seq) else NA_character_
  identity  <- if (has_identity) tbl$REBASEfinder_blast_identity else rep(NA_real_, nrow(tbl))
  operon    <- if (has_operon) as.character(tbl$REBASEfinder_operon_id) else seq_len(nrow(tbl))
  confidence <- if (has_eligible) {
    ifelse(tbl$REBASEfinder_typing_eligible == TRUE, "High", "Low")
  } else if (has_identity) {
    ifelse(!is.na(identity) & identity >= 0.5, "High", "Low")
  } else {
    rep("High", nrow(tbl))
  }

  # Coverage: blast_length / protein_length
  coverage <- rep(NA_real_, nrow(tbl))
  if (has_blastlen && has_translation) {
    prot_len <- nchar(as.character(tbl$translation))
    blast_len <- as.numeric(tbl$REBASEfinder_blast_length)
    coverage <- ifelse(prot_len > 0 & !is.na(blast_len), blast_len / prot_len, NA_real_)
  }

  # hit_label (rec_seq) [id%/cov%] | locus_tag
  base_label <- ifelse(
    !is.na(rec_seq) & nzchar(rec_seq) & rec_seq != "?",
    paste0(hit_label, " (", rec_seq, ")"),
    hit_label
  )
  id_str <- ifelse(!is.na(identity), paste0(round(identity * 100, 1), "%"), "")
  cov_str <- ifelse(!is.na(coverage), paste0(round(coverage * 100, 0), "%"), "")
  qual_tag <- ifelse(
    nzchar(id_str) & nzchar(cov_str),
    paste0(" [", id_str, "/", cov_str, "]"),
    ifelse(nzchar(id_str), paste0(" [", id_str, "]"), "")
  )
  display <- paste0(base_label, qual_tag, " | ", tbl$locus_tag)

  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  has_rmtype <- "REBASEfinder_family_id" %in% names(tbl)
  role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role) else rep(NA_character_, nrow(tbl))
  rm_type <- if (has_rmtype) as.character(tbl$REBASEfinder_family_id) else rep(NA_character_, nrow(tbl))

  out <- data.frame(
    locus_tag = tbl$locus_tag,
    display = display,
    identity = identity,
    confidence = confidence,
    operon = operon,
    role = role,
    rm_type = rm_type,
    stringsAsFactors = FALSE
  )
  # Sort: by operon group, then identity descending within group
  operon_max_id <- tapply(out$identity, out$operon, max, na.rm = TRUE)
  out$operon_rank <- operon_max_id[out$operon]
  out <- out[order(-out$operon_rank, out$operon, -out$identity), , drop = FALSE]
  out$display <- factor(out$display, levels = rev(unique(out$display)))
  out$operon_rank <- NULL
  out
}


# ====================================================================
# Panel A
# ====================================================================
.dnmb_plot_rebasefinder_inventory <- function(tbl, palette,
                                              n_total = nrow(tbl),
                                              n_hc = nrow(tbl),
                                              n_verified = nrow(tbl)) {
  if (!nrow(tbl)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = paste0("No motif-verified hits (",
                                              n_hc, " high-confidence / ", n_total, " total)"),
                               size = 4, color = "grey50"))
  }

  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  has_rec  <- "REBASEfinder_rec_seq" %in% names(tbl)
  has_id   <- "REBASEfinder_blast_identity" %in% names(tbl)

  inv <- tbl |>
    dplyr::group_by(.data$REBASEfinder_family_id) |>
    dplyr::summarise(
      n_genes = dplyr::n(),
      roles = if (has_role) paste(sort(unique(stats::na.omit(.data$REBASEfinder_enzyme_role))), collapse = "/") else "",
      rec_seqs = if (has_rec) {
        seqs <- unique(stats::na.omit(.data$REBASEfinder_rec_seq))
        seqs <- seqs[nzchar(seqs) & seqs != "?"]
        if (length(seqs)) paste(seqs[seq_len(min(3, length(seqs)))], collapse = ", ") else ""
      } else "",
      mean_identity = if (has_id) mean(.data$REBASEfinder_blast_identity, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_genes), .data$REBASEfinder_family_id)

  inv$REBASEfinder_family_id <- factor(inv$REBASEfinder_family_id, levels = rev(inv$REBASEfinder_family_id))

  inv$annot <- vapply(seq_len(nrow(inv)), function(i) {
    ng <- inv$n_genes[i]
    line1 <- paste(c(
      paste0(ng, ifelse(ng == 1, " gene", " genes")),
      if (nzchar(inv$roles[i])) paste0("roles: ", inv$roles[i]) else NULL,
      if (!is.na(inv$mean_identity[i])) paste0("identity: ", round(inv$mean_identity[i] * 100, 1), "%") else NULL
    ), collapse = " | ")
    if (nzchar(inv$rec_seqs[i])) paste0(line1, "\nrec: ", inv$rec_seqs[i]) else line1
  }, character(1))

  ggplot2::ggplot(inv, ggplot2::aes(x = .data$n_genes, y = .data$REBASEfinder_family_id)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$REBASEfinder_family_id),
                      width = 0.3, alpha = 0.9, color = "grey40", linewidth = 0.2, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(x = 0.05, label = .data$REBASEfinder_family_id),
                       hjust = 0, size = 3.2, fontface = "bold", color = "white") +
    ggplot2::geom_text(ggplot2::aes(x = .data$n_genes + 0.15, label = .data$annot),
                       hjust = 0, size = 2.5, color = "grey30") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.65)),
                                breaks = seq(0, max(inv$n_genes), by = max(1, ceiling(max(inv$n_genes) / 10)))) +
    ggplot2::labs(title = paste0("R-M system inventory (", n_verified, " verified / ",
                                  n_hc, " high-conf / ", n_total, " total)"),
                  x = "Genes detected", y = NULL) +
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


# ====================================================================
# Panel B
# ====================================================================
.dnmb_plot_rebasefinder_context <- function(genbank_table, output_dir,
                                            rm_palette, legend_position = "none") {
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  if (!nrow(tbl)) return(ggplot2::ggplot() + ggplot2::theme_void())
  for (col in c("REBASEfinder_blast_identity", "REBASEfinder_blast_bitscore",
                "REBASEfinder_blast_length", "start", "end")) {
    if (col %in% names(tbl)) tbl[[col]] <- suppressWarnings(as.numeric(tbl[[col]]))
  }
  if ("REBASEfinder_typing_eligible" %in% names(tbl)) {
    tbl$REBASEfinder_typing_eligible <- as.logical(tbl$REBASEfinder_typing_eligible)
  }
  rm_tbl <- tbl[!is.na(tbl$REBASEfinder_family_id) & nzchar(tbl$REBASEfinder_family_id), , drop = FALSE]
  if (!nrow(rm_tbl)) return(ggplot2::ggplot() + ggplot2::theme_void())

  contig_lengths <- .dnmb_contig_lengths_for_plot(tbl, output_dir = output_dir)
  # Only show contigs with RM gene hits
  hit_contigs <- unique(rm_tbl$contig)
  contig_lengths <- contig_lengths[contig_lengths$contig %in% hit_contigs, , drop = FALSE]
  contig_lengths$track <- 1
  contig_map <- stats::setNames(.dnmb_rebasefinder_short_contig(contig_lengths$contig), contig_lengths$contig)
  contig_lengths$contig_short <- contig_map[contig_lengths$contig]

  has_eligible <- "REBASEfinder_typing_eligible" %in% names(rm_tbl)
  is_hc <- if (has_eligible) !is.na(rm_tbl$REBASEfinder_typing_eligible) & rm_tbl$REBASEfinder_typing_eligible == TRUE else rep(TRUE, nrow(rm_tbl))

  has_identity <- "REBASEfinder_blast_identity" %in% names(rm_tbl)
  blast_id <- if (has_identity) rm_tbl$REBASEfinder_blast_identity else rep(NA_real_, nrow(rm_tbl))

  rm_genes <- data.frame(contig = rm_tbl$contig, contig_short = contig_map[rm_tbl$contig],
                          start = rm_tbl$start, end = rm_tbl$end,
                          rm_type = as.character(rm_tbl$REBASEfinder_family_id),
                          identity = blast_id,
                          confidence = ifelse(is_hc, "High", "Low"), stringsAsFactors = FALSE)

  for (i in seq_len(nrow(rm_genes))) {
    clen <- contig_lengths$length_bp[match(rm_genes$contig[i], contig_lengths$contig)]
    min_w <- clen * 0.008
    if ((rm_genes$end[i] - rm_genes$start[i]) < min_w) {
      mid <- (rm_genes$start[i] + rm_genes$end[i]) / 2
      rm_genes$start[i] <- mid - min_w / 2; rm_genes$end[i] <- mid + min_w / 2
    }
  }
  rm_genes$midpoint <- (rm_genes$start + rm_genes$end) / 2
  rm_genes$track <- 1

  # Compute per-gene fill color: rm_type base color blended toward white by identity
  id_frac <- ifelse(is.na(rm_genes$identity), 0.5, rm_genes$identity)
  base_cols <- rm_palette[rm_genes$rm_type]
  rm_genes$fill_color <- vapply(seq_len(nrow(rm_genes)), function(i) {
    bc <- grDevices::col2rgb(base_cols[i]) / 255
    w <- id_frac[i]  # 1 = full color, 0 = white
    blended <- bc * w + 1 * (1 - w)
    grDevices::rgb(blended[1], blended[2], blended[3])
  }, character(1))

  has_hit <- "REBASEfinder_hit_label" %in% names(rm_tbl)
  has_rec <- "REBASEfinder_rec_seq" %in% names(rm_tbl)
  short_label <- if (has_hit) { sl <- sub("ORF[0-9]+P$", "", as.character(rm_tbl$REBASEfinder_hit_label)); ifelse(nzchar(sl), sl, as.character(rm_tbl$REBASEfinder_hit_label)) } else as.character(rm_tbl$REBASEfinder_family_id)
  # Add locus_tag, then recognition sequence at the bottom
  short_label <- paste0(short_label, "\n", rm_tbl$locus_tag)
  if (has_rec) { rec <- as.character(rm_tbl$REBASEfinder_rec_seq); short_label <- ifelse(!is.na(rec) & nzchar(rec) & rec != "?", paste0(short_label, "\n(", rec, ")"), short_label) }
  rm_genes$label <- short_label
  rm_hc <- rm_genes[rm_genes$confidence == "High", , drop = FALSE]
  rm_lc <- rm_genes[rm_genes$confidence == "Low", , drop = FALSE]

  # Split labels into isolated (geom_text) vs crowded (geom_text_repel)
  # based on proximity of midpoints within each contig
  .split_crowded <- function(df, genome_frac = 0.06) {
    if (!nrow(df)) return(list(isolated = df[0, , drop = FALSE], crowded = df[0, , drop = FALSE]))
    is_crowded <- rep(FALSE, nrow(df))
    for (ctg in unique(df$contig)) {
      idx <- which(df$contig == ctg)
      if (length(idx) < 2) next
      mids <- df$midpoint[idx]
      clen <- max(df$end[idx]) - min(df$start[idx])
      if (clen <= 0) clen <- max(mids)
      thresh <- clen * genome_frac
      for (k in seq_along(idx)) {
        dists <- abs(mids[k] - mids[-k])
        if (any(dists < thresh)) is_crowded[idx[k]] <- TRUE
      }
    }
    list(isolated = df[!is_crowded, , drop = FALSE], crowded = df[is_crowded, , drop = FALSE])
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = contig_lengths, ggplot2::aes(x = 0, xend = .data$length_bp, y = .data$track, yend = .data$track), linewidth = 1.2, color = "grey85", lineend = "round")
  if (nrow(rm_lc)) {
    p <- p + ggplot2::geom_rect(data = rm_lc, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = .data$track - 0.12, ymax = .data$track + 0.12), fill = rm_lc$fill_color, color = "grey50", linewidth = 0.3, linetype = "dashed")
    p <- p + ggrepel::geom_text_repel(
      data = rm_lc,
      ggplot2::aes(x = .data$midpoint, y = .data$track + 0.18, label = .data$label),
      size = 1.7, color = "grey55", fontface = "italic", lineheight = 0.9,
      direction = "both", nudge_y = 0.10,
      segment.size = 0.12, segment.color = "grey70",
      max.overlaps = Inf, seed = 42, min.segment.length = 0.1,
      force = 2, force_pull = 0.5, box.padding = 0.25
    )
  }
  if (nrow(rm_hc)) {
    p <- p + ggplot2::geom_rect(data = rm_hc, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = .data$track - 0.12, ymax = .data$track + 0.12), fill = rm_hc$fill_color, color = "grey20", linewidth = 0.4) +
      ggplot2::geom_segment(data = rm_hc, ggplot2::aes(x = .data$midpoint, xend = .data$midpoint, y = .data$track + 0.12, yend = .data$track + 0.20, color = .data$rm_type), linewidth = 0.5, show.legend = FALSE)
    p <- p + ggrepel::geom_text_repel(
      data = rm_hc,
      ggplot2::aes(x = .data$midpoint, y = .data$track + 0.22, label = .data$label),
      size = 1.9, color = "grey20", fontface = "bold", lineheight = 0.9,
      direction = "both", nudge_y = 0.12,
      segment.size = 0.15, segment.color = "grey55",
      max.overlaps = Inf, seed = 42,
      min.segment.length = 0.1,
      force = 2, force_pull = 0.5,
      box.padding = 0.3, point.padding = 0.1
    )
  }
  # R-M Type legend (full opacity) + Identity gradient legend
  legend_df <- data.frame(rm_type = names(rm_palette), x = NA_real_, y = NA_real_, stringsAsFactors = FALSE)
  p <- p + ggplot2::geom_point(data = legend_df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$rm_type), shape = 22, size = 3, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = rm_palette, name = "R-M Type") +
    ggplot2::scale_color_manual(values = rm_palette, guide = "none") +
    # Identity gradient colorbar via ggnewscale
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(data = rm_genes, ggplot2::aes(x = .data$midpoint, y = .data$track, color = .data$identity), size = 0, na.rm = TRUE) +
    ggplot2::scale_color_gradient(low = "grey90", high = "grey20", name = "Identity",
                                   labels = scales::percent,
                                   guide = ggplot2::guide_colorbar(
                                     barwidth = 8, barheight = 0.8,
                                     title.position = "left", order = 1,
                                     nbin = 300, frame.colour = "grey40",
                                     frame.linewidth = 0.3,
                                     ticks.colour = "grey40"
                                   )) +
    ggplot2::facet_wrap(~ contig_short, ncol = 1, scales = "free_x") +
    ggplot2::scale_y_continuous(limits = c(0.4, 1.9), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(labels = function(x) ifelse(x >= 1e6, paste0(round(x / 1e6, 1), " Mb"), ifelse(x >= 1e3, paste0(round(x / 1e3), " kb"), x)), expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(title = "R-M system genome context", subtitle = "solid = high-confidence, dashed = low-confidence | color intensity = BLAST identity", x = "Position", y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", size = 9), strip.background = ggplot2::element_rect(fill = "grey95"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.title.position = "plot",
                   plot.subtitle = ggplot2::element_text(size = 8, color = "grey40"),
                   legend.position = legend_position,
                   legend.box = "horizontal",
                   legend.box.just = "left",
                   legend.spacing.x = ggplot2::unit(0.8, "cm"),
                   plot.margin = ggplot2::margin(4, 8, 4, 18))
}


# ====================================================================
# Panel C: BLAST Match Quality — legend horizontal at bottom
# ====================================================================
.dnmb_plot_rebasefinder_blast_quality <- function(tbl, role_palette, display_info, rm_palette = NULL) {
  has_identity <- "REBASEfinder_blast_identity" %in% names(tbl)
  has_bitscore <- "REBASEfinder_blast_bitscore" %in% names(tbl)
  has_role     <- "REBASEfinder_enzyme_role" %in% names(tbl)

  if (!has_identity) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No BLAST identity data", size = 4, color = "grey50"))

  plot_tbl <- merge(display_info, data.frame(
    locus_tag = tbl$locus_tag,
    bitscore = if (has_bitscore) tbl$REBASEfinder_blast_bitscore else NA_real_,
    stringsAsFactors = FALSE
  ), by = "locus_tag")
  # Use role from display_info (already included)
  if (!"role" %in% names(plot_tbl) || all(is.na(plot_tbl$role))) {
    plot_tbl$role <- "Other"
  }
  plot_tbl$role[is.na(plot_tbl$role) | !nzchar(plot_tbl$role)] <- "Other"
  plot_tbl <- plot_tbl[!is.na(plot_tbl$identity), , drop = FALSE]
  if (!nrow(plot_tbl)) return(ggplot2::ggplot() + ggplot2::theme_void())

  # Infer methylation type for M subunits + recognition sequence
  meth_type_all <- .dnmb_rebasefinder_infer_methylation_type(tbl)
  plot_tbl$meth_type <- meth_type_all[match(plot_tbl$locus_tag, tbl$locus_tag)]
  has_rec <- "REBASEfinder_rec_seq" %in% names(tbl)
  if (has_rec) {
    rec_all <- as.character(tbl$REBASEfinder_rec_seq)
    plot_tbl$rec_seq <- rec_all[match(plot_tbl$locus_tag, tbl$locus_tag)]
  } else {
    plot_tbl$rec_seq <- NA_character_
  }

  # Ensure rm_palette is available
  if (is.null(rm_palette)) rm_palette <- .dnmb_rebasefinder_palette(plot_tbl$rm_type)
  all_rm <- unique(plot_tbl$rm_type[!is.na(plot_tbl$rm_type)])
  miss_rm <- setdiff(all_rm, names(rm_palette))
  if (length(miss_rm)) rm_palette <- c(rm_palette, stats::setNames(rep("grey60", length(miss_rm)), miss_rm))

  # Build y-axis labels with rm_type prefix
  type_label <- ifelse(!is.na(plot_tbl$rm_type), paste0("[", plot_tbl$rm_type, "] "), "")
  plot_tbl$display_typed <- paste0(type_label, as.character(plot_tbl$display))
  # Preserve factor order from display
  lvls <- levels(plot_tbl$display)
  type_map <- stats::setNames(plot_tbl$display_typed, as.character(plot_tbl$display))
  new_lvls <- type_map[lvls]
  plot_tbl$display_typed <- factor(plot_tbl$display_typed, levels = new_lvls)

  # Color methylated bases in y-axis labels (rec_seq portion)
  # Uses REBASE bairoch metadata for exact position when available,
  # falls back to coloring ALL matching bases when position is unknown.
  meth_pal <- c("N6A" = "#D32F2F", "N5C" = "#1565C0", "N4C" = "#FF8F00")
  meth_target <- c("N6A" = "A", "N5C" = "C", "N4C" = "C")
  has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
  bairoch <- tryCatch(.dnmb_rebasefinder_download_bairoch(), error = function(e) NULL)
  if (has_ggtext) {
    old_lvls <- levels(plot_tbl$display_typed)
    label_vec <- as.character(plot_tbl$display_typed)
    lvl_map <- stats::setNames(old_lvls, old_lvls)
    has_hit <- "REBASEfinder_hit_label" %in% names(plot_tbl)
    for (i in seq_len(nrow(plot_tbl))) {
      mt <- plot_tbl$meth_type[i]
      rs <- plot_tbl$rec_seq[i]
      if (is.na(mt) || is.na(rs) || !nzchar(rs) || rs == "?") next
      col <- meth_pal[mt]
      if (is.na(col)) next

      # Try position-specific coloring from bairoch
      colored_rs <- NULL
      if (!is.null(bairoch) && has_hit) {
        hit_name <- as.character(plot_tbl$REBASEfinder_hit_label[i])
        if (!is.na(hit_name) && nzchar(hit_name)) {
          b_idx <- match(hit_name, bairoch$enzyme_name)
          if (!is.na(b_idx)) {
            mpos <- bairoch$meth_pos[b_idx]
            if (!is.na(mpos) && mpos != "?" && grepl("^[0-9]+$", mpos)) {
              pos <- as.integer(mpos)
              if (pos >= 1 && pos <= nchar(rs)) {
                chars <- strsplit(rs, "")[[1]]
                chars[pos] <- paste0(
                  "<span style='color:", col, ";font-weight:bold'>",
                  chars[pos], "</span>")
                colored_rs <- paste(chars, collapse = "")
              }
            }
          }
        }
      }

      # Fallback: color ALL matching target bases
      if (is.null(colored_rs)) {
        target <- meth_target[mt]
        if (!is.na(target)) {
          colored_rs <- gsub(
            target,
            paste0("<span style='color:", col, ";font-weight:bold'>", target, "</span>"),
            rs
          )
        }
      }

      if (!is.null(colored_rs)) {
        old_label <- label_vec[i]
        new_label <- sub(rs, colored_rs, old_label, fixed = TRUE)
        label_vec[i] <- new_label
        lvl_map[lvl_map == old_label] <- new_label
      }
    }
    new_lvls <- unname(lvl_map[old_lvls])
    plot_tbl$display_typed <- factor(label_vec, levels = new_lvls)
  }

  # Operon group separators
  operon_bounds <- .dnmb_rebasefinder_operon_separators(plot_tbl)

  all_roles <- unique(plot_tbl$role)
  missing_roles <- setdiff(all_roles, names(role_palette))
  if (length(missing_roles)) role_palette <- c(role_palette, stats::setNames(rep("grey60", length(missing_roles)), missing_roles))

  shade_df <- data.frame(xmin = -Inf, xmax = 50, ymin = -Inf, ymax = Inf)

  p <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$identity * 100, y = .data$display_typed)) +
    ggplot2::geom_rect(data = shade_df, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax), fill = "#FFF3E0", alpha = 0.35, inherit.aes = FALSE) +
    ggplot2::geom_vline(xintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.4)

  # Operon group separators (horizontal lines)
  if (nrow(operon_bounds)) {
    p <- p + ggplot2::geom_hline(yintercept = operon_bounds$y, color = "grey75", linewidth = 0.3, linetype = "dotted")
  }

  if (has_bitscore && any(!is.na(plot_tbl$bitscore))) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(size = .data$bitscore, fill = .data$role,
                   color = .data$rm_type, shape = .data$confidence),
      alpha = 0.8, stroke = 1.0
    ) +
      ggplot2::scale_size_continuous(name = "Bitscore", range = c(2.5, 7),
                                     breaks = scales::pretty_breaks(n = 3))
  } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(fill = .data$role, color = .data$rm_type,
                   shape = .data$confidence),
      size = 4, alpha = 0.8, stroke = 1.0
    )
  }

  x_left_limit <- max(0, min(plot_tbl$identity * 100, na.rm = TRUE) - 5)

  meth_annot_label <- if (has_ggtext && any(!is.na(plot_tbl$meth_type))) {
    paste0("<span style='color:#D32F2F'>N6A</span>",
           " / <span style='color:#1565C0'>N5C</span>",
           " / <span style='color:#FF8F00'>N4C</span>")
  } else NULL

  p + ggplot2::scale_shape_manual(values = c(High = 21, Low = 1), name = "Confidence") +
    ggplot2::annotate("text", x = 49, y = Inf, label = "low", hjust = 1, vjust = 1.5,
                      size = 2.3, color = "#BF360C", fontface = "italic") +
    ggplot2::annotate("text", x = 51, y = Inf, label = "high confidence", hjust = 0, vjust = 1.5,
                      size = 2.3, color = "#2E7D32", fontface = "italic") +
    { if (!is.null(meth_annot_label) && has_ggtext)
      ggtext::geom_richtext(
        data = data.frame(x = 80, y = Inf, label = paste0("Methylated base: ", meth_annot_label)),
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        inherit.aes = FALSE, size = 2.5, hjust = 0, vjust = 1.5,
        fill = NA, label.colour = NA, label.padding = grid::unit(0, "pt")
      )
    else NULL } +
    ggplot2::scale_fill_manual(values = role_palette, name = "Enzyme Role") +
    ggplot2::scale_color_manual(values = rm_palette, name = "R-M Type") +
    ggplot2::scale_x_continuous(limits = c(x_left_limit, 115),
                                breaks = seq(0, 100, by = 10)) +
    ggplot2::labs(
      title = "REBASE BLAST match quality",
      subtitle = NULL,
      x = "Sequence identity (%)", y = NULL
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, nrow = 1,
                                    override.aes = list(shape = 21, size = 3)),
      color = ggplot2::guide_legend(order = 2, nrow = 1,
                                     override.aes = list(shape = 21, fill = "grey80", size = 3)),
      shape = ggplot2::guide_legend(order = 3, nrow = 1,
                                     override.aes = list(size = 3, fill = "grey50", color = "grey30")),
      size = ggplot2::guide_legend(order = 4, nrow = 1, title.position = "left")
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey92"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = if (has_ggtext) ggtext::element_markdown(size = 8) else ggplot2::element_text(size = 8),
      axis.text.y = if (has_ggtext) ggtext::element_markdown(size = 7.5) else ggplot2::element_text(size = 7.5),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.key.size = ggplot2::unit(0.3, "cm"),
      legend.text = ggplot2::element_text(size = 6.5),
      legend.title = ggplot2::element_text(size = 7),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.spacing.y = ggplot2::unit(0.05, "cm"),
      legend.spacing.x = ggplot2::unit(0.1, "cm"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

# Compute operon group separator y-positions (between adjacent different operons)
.dnmb_rebasefinder_operon_separators <- function(plot_tbl) {
  if (!"operon" %in% names(plot_tbl) || nrow(plot_tbl) < 2) return(data.frame(y = numeric(0)))
  lvls <- levels(plot_tbl$display)
  ordered_operons <- plot_tbl$operon[match(lvls, plot_tbl$display)]
  seps <- which(ordered_operons[-1] != ordered_operons[-length(ordered_operons)])
  if (!length(seps)) return(data.frame(y = numeric(0)))
  data.frame(y = seps + 0.5)
}


# ====================================================================
# Motif definitions & scanning
# Literature-based patterns with positional constraints & co-occurrence
# Refs: Malone+ 1995 (amino-MTase motifs), Posfai+ 1989 / Kumar+ 1994
#   (C5-MTase I-X), Pingoud & Jeltsch 2001 (PD-ExK), Dunin-Horkawicz+ 2006
#   (GIY-YIG), Schluckebier+ 1995 (SAM fold), Murray 2000 (Type I)
# ====================================================================
.dnmb_rebasefinder_motif_definitions <- function() {
  list(
    # -- Methyltransferase motifs (M / S subunits) --
    "SAM"     = list(pattern = "[FLIVM].G[TSAG]G",
                     expected_role = "M",
                     full = "SAM binding (Motif I) FxGxG",
                     pos_range = c(0.01, 0.85),
                     weight = 3L),
    "MTase"   = list(pattern = "[DN]PP[YFW]",
                     expected_role = "M",
                     full = "N4C/N6A amino-MTase catalytic (Motif IV) [DN]PP[YFW]",
                     pos_range = c(0.20, 0.70),
                     weight = 3L),
    "C5-IV"   = list(pattern = "[PE]C[QG]",
                     expected_role = "M",
                     full = "C5-MTase catalytic Cys (Motif IV) [PE]C[QG]",
                     pos_range = c(0.30, 0.55),
                     weight = 3L),
    "C5-VI"   = list(pattern = "E[NHQS][VIL]",
                     expected_role = "M",
                     full = "C5-MTase Glu proton donor (Motif VI) E[NHQS][VIL]",
                     pos_range = c(0.40, 0.65),
                     weight = 2L),
    "C5-VII"  = list(pattern = "Q.?R.R",
                     expected_role = "M",
                     full = "C5-MTase DNA binding (Motif VII) QxRxR",
                     pos_range = c(0.50, 0.75),
                     weight = 1L),
    # -- REase / nuclease motifs (R subunits) --
    "PD-ExK"  = list(pattern = "PD.{8,20}[DE][LIVMFY]K",
                     expected_role = "R",
                     full = "REase PD-(D/E)xK catalytic triad",
                     pos_range = c(0.05, 0.65),
                     weight = 3L),
    "HNH"     = list(pattern = "H.{1,3}N.{5,40}H",
                     expected_role = "R",
                     full = "HNH nuclease (His-Asn-His)",
                     pos_range = NULL,
                     weight = 3L),
    "GIY-YIG" = list(pattern = "G[LIVMA]Y.{2,4}Y[IVLA]G",
                     expected_role = "R",
                     full = "GIY-YIG endonuclease",
                     pos_range = c(0.0, 0.40),
                     weight = 3L),
    "P-loop"  = list(pattern = "[AG].{4}GK[ST]",
                     expected_role = "R",
                     full = "Walker A / P-loop ATPase GxxxxGK[ST]",
                     pos_range = c(0.20, 0.55),
                     weight = 2L),
    "DEAD"    = list(pattern = "DE[AHCF][DHQ]",
                     expected_role = "R",
                     full = "Walker B / DEAD-box helicase (Type I HsdR)",
                     pos_range = c(0.25, 0.60),
                     weight = 2L),
    "PLD"     = list(pattern = "H[LIVMF]K.{4}D",
                     expected_role = "R",
                     full = "PLD/HKD phosphodiesterase nuclease",
                     pos_range = NULL,
                     weight = 2L),
    "Mrr"     = list(pattern = "D.{8,15}[EQ].[KR].{20,60}[DE].{0,5}[KR]",
                     expected_role = "R",
                     full = "Mrr-like Type IV REase (modified-DNA restriction)",
                     pos_range = c(0.10, 0.70),
                     weight = 2L)
  )
}

# Backward-compatible alias for legacy panel code that references old names
.dnmb_rebasefinder_motif_names_display <- function() {
  defs <- .dnmb_rebasefinder_motif_definitions()
  # Group C5 sub-motifs for display
  nms <- names(defs)
  nms[nms == "C5-IV"]  <- "C5-PC"
  nms
}

# ====================================================================
# Positional constraint check — does a hit fall within expected region?
# ====================================================================
.dnmb_rebasefinder_pos_in_range <- function(pos, prot_len, pos_range) {
  if (is.null(pos_range) || is.na(pos) || prot_len <= 0) return(TRUE)
  frac <- pos / prot_len
  frac >= pos_range[1] && frac <= pos_range[2]
}

# ====================================================================
# Co-occurrence scoring — confidence level per gene
# ====================================================================
.dnmb_rebasefinder_motif_score <- function(gene_result, role, prot_len,
                                            motif_defs = .dnmb_rebasefinder_motif_definitions()) {
  motif_names <- names(motif_defs)
  score <- 0L
  has <- character(0)
  for (j in seq_along(motif_names)) {
    mn <- motif_names[j]
    info <- gene_result[[j]]
    if (is.null(info) || info$n_hits == 0L) next
    # Check positional constraint for best hit
    in_range <- .dnmb_rebasefinder_pos_in_range(info$pos, prot_len, motif_defs[[mn]]$pos_range)
    w <- motif_defs[[mn]]$weight
    score <- score + w + if (in_range) 2L else 0L
    has <- c(has, mn)
  }
  # Co-occurrence bonuses
  if (all(c("SAM", "MTase") %in% has)) score <- score + 3L      # amino-MTase pair
  if (all(c("SAM", "C5-IV", "C5-VI") %in% has)) score <- score + 3L  # C5-MTase triad
  if (all(c("P-loop", "DEAD") %in% has)) score <- score + 3L    # Type I HsdR helicase
  # Protein length penalty for very short ORFs
  min_len <- if (!is.na(role) && role %in% c("M", "S")) 200L else if (!is.na(role) && role == "R") 150L else 100L
  if (prot_len < min_len) score <- max(0L, score - 3L)
  list(score = score, motifs_found = has)
}

# Returns a list per gene: for each motif, ALL hits (gregexpr)
# Each hit: list(status, hits = list of list(pos, end, match))
.dnmb_rebasefinder_scan_motifs_detailed <- function(tbl) {
  has_translation <- "translation" %in% names(tbl)
  if (!has_translation) return(NULL)

  motif_defs <- .dnmb_rebasefinder_motif_definitions()
  motif_names <- names(motif_defs)
  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)

  results <- lapply(seq_len(nrow(tbl)), function(i) {
    seq <- toupper(as.character(tbl$translation[i]))
    prot_len <- nchar(seq)
    if (is.na(seq) || !nzchar(seq)) {
      return(lapply(motif_names, function(mn) list(
        status = NA, hits = list(), n_hits = 0L, prot_len = 0L,
        pos = NA_integer_, match = NA_character_, match_len = NA_integer_
      )))
    }
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_

    lapply(motif_names, function(mn) {
      m <- gregexpr(motif_defs[[mn]]$pattern, seq)[[1]]
      if (m[1] > 0) {
        mlens <- attr(m, "match.length")
        all_hits <- lapply(seq_along(m), function(k) {
          p <- as.integer(m[k])
          e <- as.integer(m[k] + mlens[k] - 1L)
          in_range <- .dnmb_rebasefinder_pos_in_range(p, prot_len, motif_defs[[mn]]$pos_range)
          list(pos = p, end = e,
               match = substr(seq, m[k], m[k] + mlens[k] - 1L),
               in_range = in_range)
        })
        # Prefer hits within positional range
        range_hits <- Filter(function(h) h$in_range, all_hits)
        best_hits <- if (length(range_hits)) range_hits else all_hits
        expected <- motif_defs[[mn]]$expected_role
        in_pos <- length(range_hits) > 0
        status <- if (!is.na(role) && role == expected && in_pos) "present"
                  else if (!is.na(role) && role == expected) "present~"
                  else if (in_pos) "present*"
                  else "present*~"
        list(status = status, hits = all_hits, n_hits = length(all_hits),
             prot_len = prot_len, n_in_range = length(range_hits),
             pos = best_hits[[1]]$pos, match = best_hits[[1]]$match,
             match_len = nchar(best_hits[[1]]$match))
      } else {
        list(status = "absent", hits = list(), n_hits = 0L,
             prot_len = prot_len, n_in_range = 0L,
             pos = NA_integer_, match = NA_character_, match_len = NA_integer_)
      }
    })
  })
  stats::setNames(results, tbl$locus_tag)
}

# Simplified scan for motif_verified flag
.dnmb_rebasefinder_scan_motifs <- function(tbl) {
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  if (is.null(detailed)) return(NULL)
  motif_names <- names(.dnmb_rebasefinder_motif_definitions())
  do.call(rbind, lapply(names(detailed), function(lt) {
    row <- vapply(motif_names, function(mn) detailed[[lt]][[which(motif_names == mn)]]$status %||% NA_character_, character(1))
    as.data.frame(as.list(row), stringsAsFactors = FALSE)
  })) -> df
  df$locus_tag <- names(detailed)
  df
}

# Methylation type inference from catalytic motifs
# C5-PC motif → N5C; amino-MTase motif → N6A (most common); fallback NA
.dnmb_rebasefinder_infer_methylation_type <- function(tbl) {
  # Strategy (in priority order):
  # 1. REBASE bairoch metadata — authoritative, per-enzyme lookup
  # 2. Protein motif detection — heuristic fallback
  # 3. Recognition sequence heuristic — last resort
  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  has_hit <- "REBASEfinder_hit_label" %in% names(tbl)
  has_rec <- "REBASEfinder_rec_seq" %in% names(tbl)

  # Pre-load bairoch lookup (cached after first download)
  bairoch <- tryCatch(.dnmb_rebasefinder_download_bairoch(), error = function(e) NULL)

  # Motif detection (fallback)
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  motif_names <- if (!is.null(detailed)) names(.dnmb_rebasefinder_motif_definitions()) else character()
  c5_idx <- which(motif_names == "C5-IV")
  mtase_idx <- which(motif_names == "MTase")

  vapply(seq_len(nrow(tbl)), function(i) {
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_
    if (is.na(role) || !role %in% c("M", "S")) return(NA_character_)

    # --- 1. REBASE bairoch lookup ---
    if (!is.null(bairoch) && has_hit) {
      hit_name <- as.character(tbl$REBASEfinder_hit_label[i])
      if (!is.na(hit_name) && nzchar(hit_name)) {
        idx <- match(hit_name, bairoch$enzyme_name)
        if (!is.na(idx)) {
          mt <- bairoch$meth_type[idx]
          if (!is.na(mt) && nzchar(mt)) {
            if (grepl("m6A", mt)) return("N6A")
            if (grepl("m5C", mt)) return("N5C")
            if (grepl("m4C|Nm4C", mt)) return("N4C")
          }
        }
      }
    }

    # --- 2. Motif detection ---
    if (!is.null(detailed) && length(c5_idx) && length(mtase_idx)) {
      c5_hit <- !is.na(detailed[[i]][[c5_idx]]$status) &&
        detailed[[i]][[c5_idx]]$status %in% c("present", "present*")
      mtase_hit <- !is.na(detailed[[i]][[mtase_idx]]$status) &&
        detailed[[i]][[mtase_idx]]$status %in% c("present", "present*")
      if (c5_hit) return("N5C")
      if (mtase_hit) return("N6A")
    }

    # --- 3. Recognition sequence heuristic ---
    if (has_rec) {
      rs <- as.character(tbl$REBASEfinder_rec_seq[i])
      if (!is.na(rs) && nzchar(rs) && rs != "?") {
        bases <- toupper(gsub("[^ACGTWSMKRYBDHVN]", "", rs))
        if (grepl("A", bases)) return("N6A")
        if (grepl("C", bases) && !grepl("A", bases)) return("N5C")
      }
    }

    NA_character_
  }, character(1))
}

# Motif-verified flag — uses co-occurrence scoring
# "present" = correct role + in positional range
# "present~" = correct role but outside expected position
# "present*" / "present*~" = wrong role
.dnmb_rebasefinder_motif_verified <- function(tbl) {
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  if (is.null(detailed)) return(rep(TRUE, nrow(tbl)))

  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  has_translation <- "translation" %in% names(tbl)
  motif_defs <- .dnmb_rebasefinder_motif_definitions()
  motif_names <- names(motif_defs)

  vapply(seq_len(nrow(tbl)), function(i) {
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_
    if (is.na(role) || !nzchar(role)) return(TRUE)
    prot_len <- if (has_translation) nchar(as.character(tbl$translation[i])) else 0L
    # Co-occurrence score
    sc <- .dnmb_rebasefinder_motif_score(detailed[[i]], role, prot_len, motif_defs)
    # Require at least one role-matched motif in expected position
    expected_idx <- which(vapply(motif_names, function(mn) motif_defs[[mn]]$expected_role == role, logical(1)))
    if (!length(expected_idx)) return(TRUE)
    has_match <- any(vapply(expected_idx, function(j) {
      s <- detailed[[i]][[j]]$status
      !is.na(s) && grepl("^present", s)
    }, logical(1)))
    has_match && sc$score >= 2L
  }, logical(1))
}


# ====================================================================
# Panel D: Motif Verification Tile — with residue + position
# ====================================================================
.dnmb_plot_rebasefinder_motif_verification <- function(tbl, display_info) {
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  if (is.null(detailed)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No protein sequences", size = 3.5, color = "grey50"))
  }

  motif_defs <- .dnmb_rebasefinder_motif_definitions()
  motif_names <- names(motif_defs)
  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)

  # Build long table: only show motifs relevant to each gene's role
  # M/S subunit → SAM binding + MTase catalytic (indices 1,2)
  # R subunit → REase + HNH + GIY-YIG + P-loop (indices 3,4,5,6)
  has_hit <- "REBASEfinder_hit_label" %in% names(tbl)
  long <- do.call(rbind, lapply(seq_along(detailed), function(i) {
    lt <- names(detailed)[i]
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_
    # Fallback: infer role from hit name prefix
    if ((is.na(role) || !nzchar(role)) && has_hit) {
      hn <- as.character(tbl$REBASEfinder_hit_label[i])
      if (!is.na(hn) && grepl("^M[0-9]?\\.", hn)) role <- "M"
      else if (!is.na(hn) && grepl("^R[0-9]?\\.", hn)) role <- "R"
      else if (!is.na(hn) && grepl("^S[0-9]?\\.", hn)) role <- "S"
    }

    do.call(rbind, lapply(seq_along(motif_names), function(j) {
      info <- detailed[[i]][[j]]
      expected <- motif_defs[[motif_names[j]]]$expected_role

      is_relevant <- if (is.na(role) || !nzchar(role)) TRUE else role == expected
      if (!is.na(role) && role == "S") is_relevant <- expected == "M"

      if (is_relevant) {
        status <- if (!is.na(info$status) && grepl("^present", info$status)) "found" else "missing"
        # Distinguish positional: "found" (in range) vs "found_oop" (out of position)
        if (status == "found" && !is.null(info$n_in_range) && info$n_in_range == 0L) {
          status <- "found_oop"
        }
      } else {
        status <- "na"
      }

      # Build label from ALL hits
      if (status %in% c("found", "found_oop") && info$n_hits > 0) {
        hit_labels <- vapply(info$hits, function(h) {
          trunc <- if (nchar(h$match) > 5) paste0(substr(h$match, 1, 4), "..") else h$match
          paste0(trunc, "\n", h$pos, "-", h$end)
        }, character(1))
        # Show up to 3 hits
        if (length(hit_labels) > 3) {
          tile_label <- paste0(paste(hit_labels[1:3], collapse = "\n"), "\n+", length(hit_labels) - 3, " more")
        } else {
          tile_label <- paste(hit_labels, collapse = "\n")
        }
        n_found <- info$n_hits
      } else {
        tile_label <- ""
        n_found <- 0L
      }

      data.frame(
        locus_tag = lt,
        motif = motif_names[j],
        status = status,
        expected_role = expected,
        tile_label = tile_label,
        n_found = n_found,
        stringsAsFactors = FALSE
      )
    }))
  }))

  long <- merge(long, display_info[, c("locus_tag", "display", "operon")], by = "locus_tag")
  long$motif <- factor(long$motif, levels = motif_names)

  # Status palette: role-specific colors, with out-of-position variants
  long$fill_status <- ifelse(
    long$status == "found" & long$expected_role == "R", "found_R",
    ifelse(long$status == "found", "found_M",
    ifelse(long$status == "found_oop" & long$expected_role == "R", "oop_R",
    ifelse(long$status == "found_oop", "oop_M", long$status)))
  )

  status_pal <- c(
    "found_M" = "#2E7D32",   # green — MTase motif in expected position
    "found_R" = "#C2185B",   # rose — REase motif in expected position
    "oop_M"   = "#A5D6A7",   # light green — MTase motif out of expected position
    "oop_R"   = "#F48FB1",   # light pink — REase motif out of expected position
    "missing" = "#BDBDBD",   # grey
    "na"      = "#F5F5F5"    # very light grey
  )

  # Operon separators
  operon_bounds <- .dnmb_rebasefinder_operon_separators(
    merge(display_info, data.frame(locus_tag = tbl$locus_tag), by = "locus_tag")
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$motif, y = .data$display)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill_status), color = "white", linewidth = 0.6)

  if (nrow(operon_bounds)) {
    p <- p + ggplot2::geom_hline(yintercept = operon_bounds$y, color = "grey75", linewidth = 0.3, linetype = "dotted")
  }

  p +
    ggplot2::geom_text(ggplot2::aes(label = .data$tile_label),
                       size = 1.5, color = "grey10", lineheight = 0.75) +
    ggplot2::scale_fill_manual(
      values = status_pal, name = "Motif",
      labels = c("found_M" = "MTase", "found_R" = "REase",
                 "oop_M" = "MTase (oop)", "oop_R" = "REase (oop)",
                 "missing" = "Missing", "na" = "N/A"),
      breaks = c("found_M", "found_R", "oop_M", "oop_R", "missing", "na")
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::labs(title = "  Functional motif verification", x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(size = 7, angle = 30, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      legend.key.size = ggplot2::unit(0.3, "cm"),
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 6.5),
      legend.spacing.x = ggplot2::unit(0.1, "cm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(4, 4, 4, 2)
    )
}


# ====================================================================
# UniProt domain query (optional, requires internet + jsonlite)
# Returns data.frame: locus_tag, type, description, start, end
# Cached in output_dir to avoid repeated API calls
# ====================================================================
.dnmb_rebasefinder_uniprot_domains <- function(tbl, output_dir = NULL) {
  if (!"protein_id" %in% names(tbl)) return(NULL)
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)

  # Check cache first
  cache_path <- NULL
  if (!is.null(output_dir)) {
    cache_dir <- file.path(output_dir, "dnmb_module_rebasefinder")
    if (dir.exists(cache_dir)) {
      cache_path <- file.path(cache_dir, "uniprot_domains_cache.rds")
      if (file.exists(cache_path)) {
        cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
        if (!is.null(cached)) return(cached)
      }
    }
  }

  results <- tryCatch({
    do.call(rbind, lapply(seq_len(nrow(tbl)), function(i) {
      pid <- as.character(tbl$protein_id[i])
      lt <- tbl$locus_tag[i]
      if (is.na(pid) || !nzchar(pid)) return(NULL)

      url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=xref:refseq-", pid,
                    "&format=json&fields=accession,ft_act_site,ft_binding,ft_domain,ft_motif,ft_site")
      resp <- tryCatch({
        Sys.sleep(0.2)
        jsonlite::fromJSON(url, simplifyDataFrame = TRUE)
      }, error = function(e) NULL)

      if (is.null(resp) || !length(resp$results) || !nrow(resp$results)) return(NULL)
      feats <- resp$results$features[[1]]
      if (is.null(feats) || !nrow(feats)) return(NULL)

      data.frame(locus_tag = lt, type = feats$type, description = feats$description,
                 start = feats$location$start$value, end = feats$location$end$value,
                 stringsAsFactors = FALSE)
    }))
  }, error = function(e) NULL)

  # Save cache
  if (!is.null(cache_path) && !is.null(results)) {
    tryCatch(saveRDS(results, cache_path), error = function(e) NULL)
  }
  results
}


# ====================================================================
# Panel D (new): Protein domain map — motif positions along protein
# ====================================================================
.dnmb_plot_rebasefinder_domain_map <- function(tbl, display_info, uniprot_doms = NULL) {
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  has_translation <- "translation" %in% names(tbl)
  if (is.null(detailed) || !has_translation) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  motif_defs <- .dnmb_rebasefinder_motif_definitions()
  motif_names <- names(motif_defs)

  # Motif color palette — M motifs distinguished by methylation type, R motifs red tones
  motif_pal <- c(
    "SAM"     = "#2E7D32",   # dark green (M - SAM binding)
    "MTase"   = "#9CCC65",   # yellow-green / 연두 (M - N4C/N6A amino-MTase)
    "C5-IV"   = "#1565C0",   # blue (M - N5C catalytic Cys)
    "C5-VI"   = "#42A5F5",   # light blue (M - N5C Glu donor)
    "C5-VII"  = "#90CAF9",   # pale blue (M - N5C DNA binding)
    "PD-ExK"  = "#C62828",   # dark red (R)
    "HNH"     = "#EF5350",   # light red (R)
    "GIY-YIG" = "#F48FB1",   # pink (R)
    "P-loop"  = "#E91E63",   # rose (R)
    "DEAD"    = "#AD1457",   # magenta (R - Walker B)
    "PLD"     = "#880E4F"    # dark magenta (R - PLD nuclease)
  )

  # Build protein backbones
  prot_lengths <- nchar(as.character(tbl$translation))

  backbone <- merge(
    data.frame(locus_tag = tbl$locus_tag, prot_len = prot_lengths, stringsAsFactors = FALSE),
    display_info[, c("locus_tag", "display")],
    by = "locus_tag"
  )

  # Motif markers — ALL hits from gregexpr
  markers <- do.call(rbind, lapply(seq_along(detailed), function(i) {
    lt <- names(detailed)[i]
    do.call(rbind, lapply(seq_along(motif_names), function(j) {
      info <- detailed[[i]][[j]]
      if (info$n_hits > 0) {
        do.call(rbind, lapply(info$hits, function(h) {
          data.frame(locus_tag = lt, pos = h$pos,
                     motif_short = motif_names[j],
                     stringsAsFactors = FALSE)
        }))
      } else NULL
    }))
  }))

  if (is.null(markers) || !nrow(markers)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  markers <- merge(markers, display_info[, c("locus_tag", "display")], by = "locus_tag")
  markers <- merge(markers, backbone[, c("locus_tag", "prot_len")], by = "locus_tag")

  # Normalize position to fraction of protein length
  markers$frac <- markers$pos / markers$prot_len

  p <- ggplot2::ggplot() +
    # Protein backbone — geom_tile with exact height matching domain rects
    ggplot2::geom_tile(
      data = backbone,
      ggplot2::aes(x = 0.5, y = .data$display, width = 1, height = 0.30),
      fill = "grey80", color = NA
    )

  # UniProt annotation overlay
  if (!is.null(uniprot_doms) && nrow(uniprot_doms)) {
    dom_plot <- merge(uniprot_doms, backbone[, c("locus_tag", "display", "prot_len")], by = "locus_tag")
    dom_plot$frac_start <- dom_plot$start / dom_plot$prot_len
    dom_plot$frac_end <- dom_plot$end / dom_plot$prot_len
    dom_plot$label <- ifelse(
      !is.na(dom_plot$description) & nzchar(dom_plot$description),
      dom_plot$description, dom_plot$type
    )

    # Separate domains (regions) vs residue-level features (sites)
    is_region <- (dom_plot$end - dom_plot$start) > 5
    dom_regions <- dom_plot[is_region, , drop = FALSE]
    dom_sites   <- dom_plot[!is_region, , drop = FALSE]

    # Domain regions — semi-transparent rectangles
    if (nrow(dom_regions)) {
      p <- p +
        ggplot2::geom_rect(
          data = dom_regions,
          ggplot2::aes(xmin = .data$frac_start, xmax = .data$frac_end,
                       ymin = as.numeric(.data$display) - 0.15,
                       ymax = as.numeric(.data$display) + 0.15),
          fill = "#1565C0", alpha = 0.12, color = "#1565C0",
          linewidth = 0.3, linetype = "solid"
        ) +
        ggplot2::geom_text(
          data = dom_regions,
          ggplot2::aes(x = (.data$frac_start + .data$frac_end) / 2,
                       y = as.numeric(.data$display) - 0.18,
                       label = .data$label),
          size = 1.3, color = "#0D47A1", vjust = 1, fontface = "italic",
          check_overlap = TRUE
        )
    }

    # Residue-level features (active site, binding site) — vertical lines
    if (nrow(dom_sites)) {
      dom_sites$frac_mid <- (dom_sites$frac_start + dom_sites$frac_end) / 2
      dom_sites$site_label <- paste0(substr(dom_sites$type, 1, 1), "@", dom_sites$start)
      p <- p +
        ggplot2::geom_segment(
          data = dom_sites,
          ggplot2::aes(x = .data$frac_mid, xend = .data$frac_mid,
                       y = as.numeric(.data$display) - 0.15,
                       yend = as.numeric(.data$display) + 0.15),
          color = "#FF6F00", linewidth = 0.5, alpha = 0.8
        ) +
        ggplot2::geom_text(
          data = dom_sites,
          ggplot2::aes(x = .data$frac_mid,
                       y = as.numeric(.data$display) + 0.18,
                       label = .data$site_label),
          size = 1.3, color = "#E65100", fontface = "bold",
          vjust = 0, check_overlap = TRUE
        )
    }
  }

  p <- p +
    # Motif position markers
    ggplot2::geom_point(
      data = markers,
      ggplot2::aes(x = .data$frac, y = .data$display, color = .data$motif_short),
      size = 2.5, alpha = 0.9, shape = 18
    )

  p <- p +
    # Protein length annotation
    ggplot2::geom_text(
      data = backbone,
      ggplot2::aes(x = 1.08, y = .data$display, label = paste0(.data$prot_len, " aa")),
      hjust = 0, size = 2.0, color = "grey50"
    ) +
    ggplot2::scale_color_manual(values = motif_pal, name = "Motif") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::scale_x_continuous(
      limits = c(-0.02, 1.25),
      breaks = c(0, 0.5, 1),
      labels = c("N", "50%", "C")
    ) +
    ggplot2::labs(title = "  Protein domain map", x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 7),
      legend.position = "bottom",
      legend.justification = "left",
      legend.key.size = ggplot2::unit(0.25, "cm"),
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 6),
      legend.spacing.x = ggplot2::unit(0.08, "cm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(4, 2, 4, 2)
    )
}

.dnmb_iselement_plot_module_dir <- function(output_dir) {
  candidates <- c(
    base::file.path(output_dir, "dnmb_module_iselement"),
    base::file.path(output_dir, "dnmb_module_ISelement")
  )
  existing <- candidates[base::dir.exists(candidates)]
  if (base::length(existing)) {
    return(existing[[1]])
  }
  candidates[[1]]
}

.dnmb_plot_iselement_module <- function(genbank_table, output_dir) {
  genbank_table <- .dnmb_backfill_iselement_columns(genbank_table, output_dir = output_dir)
  tbl <- .dnmb_contig_ordered_table(genbank_table)
  if (!base::nrow(tbl) || !"ISelement_element_id" %in% base::names(tbl)) {
    return(NULL)
  }

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  module_dir <- .dnmb_iselement_plot_module_dir(output_dir)

  # --- Check for full analysis data ---
  census_path <- base::file.path(module_dir, "iselement_census.tsv")
  lp_path <- base::file.path(module_dir, "iselement_landing_pads.tsv")
  elements_path <- base::file.path(module_dir, "iselement_elements.tsv")
  tm_path <- base::file.path(module_dir, "iselement_target_models.tsv")
  has_full <- base::file.exists(census_path) && base::file.exists(lp_path)

  if (has_full) {
    result <- tryCatch(
      .dnmb_plot_iselement_comprehensive(tbl, output_dir, plot_dir,
                                          census_path, lp_path, elements_path, tm_path),
      error = function(e) {
        message("[DNMB plot] Comprehensive IS plot failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(result)) {
      # Also generate separate landing pad linear map
      tryCatch(
        .dnmb_plot_iselement_landing_pads(tbl, output_dir),
        error = function(e) {
          message("[DNMB plot] Landing pad map failed: ", conditionMessage(e))
        }
      )
      return(result)
    }
  }

  NULL
}

# ---------------------------------------------------------------------------
# Comprehensive IS Element Overview: circlize genome map + ggplot2 panels
# Two-page PDF: Page 1 = circlize circular map with IS links,
#               Page 2 = ggplot2 census + recognition + landing pad heatmap
# ---------------------------------------------------------------------------
.dnmb_plot_iselement_comprehensive <- function(genbank_table, output_dir, plot_dir,
                                                census_path, lp_path, elements_path, tm_path) {
  census <- utils::read.delim(census_path, check.names = FALSE)
  landing_pads <- utils::read.delim(lp_path, check.names = FALSE)
  elements <- if (base::file.exists(elements_path)) utils::read.delim(elements_path, check.names = FALSE) else data.frame()
  target_models <- if (base::file.exists(tm_path)) utils::read.delim(tm_path, check.names = FALSE) else data.frame()
  if (!nrow(census) || !nrow(landing_pads)) return(NULL)

  genbank_candidates <- list.files(output_dir, pattern = "\\.(gbff|gbk|gb)$", full.names = TRUE)
  if (!length(genbank_candidates)) genbank_candidates <- list.files(dirname(output_dir), pattern = "\\.(gbff|gbk|gb)$", full.names = TRUE)
  genes <- NULL
  if (length(genbank_candidates)) {
    genes <- tryCatch({
      parsed <- .dnmb_parse_genbank_features(genbank_candidates[1])
      .dnmb_predict_gene_essentiality(.dnmb_build_gene_table(parsed$features))
    }, error = function(e) NULL)
  }
  if (is.null(genes) || !nrow(genes)) return(NULL)

  genome_len <- max(genes$end, na.rm = TRUE)
  main_contig <- genes$contig[which.max(genes$end)]

  ess_track <- genes[genes$contig == main_contig, , drop = FALSE]
  is_track <- if (nrow(elements)) {
    el <- elements[elements$contig == main_contig, , drop = FALSE]
    el$fam <- ifelse(is.na(el$element_family) | !nzchar(el$element_family), "unknown", el$element_family)
    el
  } else data.frame()
  lp_track <- landing_pads[landing_pads$contig == main_contig, , drop = FALSE]

  fam_cols <- c(IS110="#7E57C2",IS982="#26A69A",IS701="#EF5350",IS630="#42A5F5",
                IS66="#FFA726",IS4="#AB47BC",IS21="#26C6DA",IS3="#EC407A",
                IS481="#5C6BC0",ISL3="#66BB6A","IS200/IS605"="#8D6E63",
                IS256="#78909C",IS5="#D4E157",IS1595="#FF7043",unknown="#BDBDBD")

  pdf_path <- base::file.path(plot_dir, "ISelement_overview.pdf")

  # Use optimized circlize script if available
  ctg <- main_contig
  out_file <- pdf_path
  opt_script <- base::system.file("scripts", "iselement_circlize_optimized.R", package = "DNMB")
  if (!nzchar(opt_script) || !base::file.exists(opt_script)) {
    # Dev mode: look in R/ directory
    candidate_dirs <- c(
      base::file.path(base::dirname(output_dir), "R"),
      base::file.path(base::dirname(base::dirname(output_dir)), "R"),
      base::file.path(base::dirname(base::dirname(plot_dir)), "R")
    )
    for (d in candidate_dirs) {
      p <- base::file.path(d, "_iselement_circlize_optimized.R")
      if (base::file.exists(p)) { opt_script <- p; break }
    }
  }
  if (nzchar(opt_script) && base::file.exists(opt_script)) {
    # Set context for the optimized script
    gbff_files <- list.files(output_dir, pattern = "\\.(gbff|gbk|gb)$", full.names = TRUE)
    if (!length(gbff_files)) gbff_files <- list.files(dirname(output_dir), pattern = "\\.(gbff|gbk|gb)$", full.names = TRUE)
    dnmb_iselement_context <- list(
      pkg_dir = base::system.file(package = "DNMB"),
      output_dir = output_dir,
      plot_dir = plot_dir,
      genbank_path = if (length(gbff_files)) gbff_files[1] else "",
      parsed = list(features = genes),
      genes = genes,
      elements = elements,
      census = census,
      landing_pads = landing_pads,
      target_models = target_models,
      full_result = list(elements = elements, census = census, target_models = target_models, landing_pads = landing_pads),
      main_contig = main_contig,
      out_file = pdf_path
    )
    base::source(opt_script, local = TRUE)
    return(list(pdf = pdf_path, comprehensive = TRUE))
  }

  # Fallback: inline circlize (simplified)
  grDevices::pdf(pdf_path, width = 14, height = 14)

  # ---- Page 1: Circlize genome map ----
  graphics::par(mar=c(2, 2, 3, 2))
  circlize::circos.clear()
  circlize::circos.par(start.degree=90, gap.degree=3, cell.padding=c(0,0,0,0),
                       track.margin=c(0.005,0.005), points.overflow.warning=FALSE)
  genome_bed <- data.frame(chr=main_contig, start=0, end=genome_len)
  circlize::circos.genomicInitialize(genome_bed, plotType=NULL, major.by=200000, axis.labels.cex=0.55)
  graphics::title(paste0("IS Element Genome Distribution & Landing Pad Map\n", main_contig,
                          "  (", round(genome_len/1e6,2), " Mb)  |  ", nrow(is_track),
                          " IS elements, ", length(unique(is_track$fam)), " families"),
                  cex.main=1.15, font.main=2, line=0.8)

  # Coordinate axis
  circlize::circos.track(ylim=c(0,1), track.height=0.03, bg.border=NA, panel.fun=function(x,y) {
    circlize::circos.genomicAxis(h="top", major.by=500000, labels.cex=0.5, labels.facing="clockwise")
  })

  # Track 1: Gene Essentiality (colorRamp2 gradient)
  ess_bed <- data.frame(chr=main_contig, start=ess_track$start, end=ess_track$end, value=ess_track$essentiality_score)
  ess_bed <- ess_bed[!is.na(ess_bed$start) & !is.na(ess_bed$end) & !is.na(ess_bed$value), , drop=FALSE]
  ess_col_fun <- circlize::colorRamp2(c(0,0.15,0.35,0.55,0.75,1.0), c("#E3F2FD","#90CAF9","#42A5F5","#1E88E5","#1565C0","#0D47A1"))
  if (nrow(ess_bed)) circlize::circos.genomicTrack(ess_bed, track.height=0.055, bg.border="grey60", bg.lwd=0.4,
                                panel.fun=function(region, value, ...) {
    circlize::circos.genomicRect(region, value, ytop.column=1, ybottom=0, col=ess_col_fun(value[[1]]), border=NA)
  })

  # Track 2: IS Elements (genomicTrack, family-colored)
  is_bed <- data.frame(chr=main_contig, start=pmin(is_track$start, is_track$end), end=pmax(is_track$start, is_track$end), family=is_track$fam)
  is_bed <- is_bed[!is.na(is_bed$start) & !is.na(is_bed$end), , drop=FALSE]
  if (!nrow(is_bed)) { grDevices::dev.off(); return(NULL) }
  circlize::circos.genomicTrack(is_bed, ylim=c(0,1), track.height=0.10, bg.border="grey60", bg.lwd=0.4,
                                panel.fun=function(region, value, ...) {
    for (i in seq_len(nrow(region))) {
      col <- fam_cols[value$family[i]]; if (is.na(col)) col <- "#BDBDBD"
      circlize::circos.rect(region$start[i], 0.05, region$end[i], 0.95, col=col, border=grDevices::adjustcolor("grey30",0.6), lwd=0.3)
    }
  })

  # Track 3: Landing Pad Score (bar chart with top-15 highlighting)
  lp_bed <- data.frame(chr=main_contig, start=lp_track$region_start, end=lp_track$region_end, score=lp_track$landing_pad_score)
  lp_bed <- lp_bed[!is.na(lp_bed$start) & !is.na(lp_bed$end) & !is.na(lp_bed$score), , drop=FALSE]
  lp_col_fun <- circlize::colorRamp2(c(0.55,0.65,0.73,0.80,0.87), c("#E53935","#FF8F00","#FDD835","#66BB6A","#1B5E20"))
  if (nrow(lp_bed)) circlize::circos.genomicTrack(lp_bed, track.height=0.085, bg.border="grey60", bg.lwd=0.4,
                                panel.fun=function(region, value, ...) {
    scores <- value[[1]]
    norm_h <- pmin(0.95, pmax(0.08, (scores - 0.55) / (0.87 - 0.55) * 0.87 + 0.08))
    for (i in seq_len(nrow(region))) {
      circlize::circos.rect(region$start[i], 0, region$end[i], norm_h[i], col=lp_col_fun(scores[i]), border=NA)
    }
    top_idx <- order(scores, decreasing=TRUE)[seq_len(min(15, length(scores)))]
    for (i in top_idx) circlize::circos.rect(region$start[i], 0, region$end[i], norm_h[i], col=NA, border="grey20", lwd=0.4)
    ref_y <- (0.80 - 0.55) / (0.87 - 0.55) * 0.87 + 0.08
    circlize::circos.lines(circlize::CELL_META$cell.xlim, c(ref_y, ref_y), lty=3, col="grey30", lwd=0.5)
  })

  # Track 4: IS density (circos.genomicDensity)
  is_density_bed <- data.frame(chr=main_contig, start=pmin(is_track$start, is_track$end), end=pmax(is_track$start, is_track$end))
  is_density_bed <- is_density_bed[!is.na(is_density_bed$start) & !is.na(is_density_bed$end), , drop=FALSE]
  if (nrow(is_density_bed) > 1) {
    circlize::circos.genomicDensity(is_density_bed, track.height=0.065, bg.border="grey60", bg.lwd=0.4,
                                    col="#7B1FA2", window.size=50000, overlap=TRUE)
  }

  # Track 5: Essential gene density
  ess_high_idx <- !is.na(ess_track$essentiality_score) & ess_track$essentiality_score >= 0.45
  ess_high_bed <- data.frame(chr=main_contig,
    start=pmin(ess_track$start[ess_high_idx], ess_track$end[ess_high_idx]),
    end=pmax(ess_track$start[ess_high_idx], ess_track$end[ess_high_idx]))
  ess_high_bed <- ess_high_bed[!is.na(ess_high_bed$start) & !is.na(ess_high_bed$end), , drop=FALSE]
  if (nrow(ess_high_bed) > 2) {
    circlize::circos.genomicDensity(ess_high_bed, track.height=0.055, bg.border="grey60", bg.lwd=0.4,
                                    col="#C62828", window.size=60000, overlap=TRUE)
  }

  # IS family links — connect same-family IS elements
  link_families <- is_track %>% dplyr::count(.data$fam, sort=TRUE) %>%
    dplyr::filter(.data$n >= 2, .data$fam != "unknown") %>% dplyr::pull(.data$fam)
  for (fam_name in link_families) {
    fam_el <- is_track %>% dplyr::filter(.data$fam == fam_name) %>% dplyr::arrange(.data$start)
    if (nrow(fam_el) < 2) next
    col <- fam_cols[fam_name]; if (is.na(col)) col <- "#9E9E9E"
    link_col <- grDevices::adjustcolor(col, alpha.f=0.18)
    border_col <- grDevices::adjustcolor(col, alpha.f=0.40)
    for (i in seq_len(nrow(fam_el)-1)) {
      dist <- fam_el$start[i+1] - fam_el$end[i]
      if (dist > 30000) {
        h <- if (dist > 500000) 0.5 else if (dist > 200000) 0.4 else 0.3
        circlize::circos.link(main_contig, c(fam_el$start[i], fam_el$end[i]),
                              main_contig, c(fam_el$start[i+1], fam_el$end[i+1]),
                              col=link_col, border=border_col, lwd=0.4, h.ratio=h)
      }
    }
  }

  circlize::circos.track(ylim=c(0,1), track.height=0.001, bg.border=NA)

  # Center legends
  graphics::text(0, 0.35, "Track legend (outer to inner):", cex=0.85, font=2)
  graphics::text(0, 0.28, "1. Gene Essentiality  2. IS Elements  3. Landing Pad Score  4. IS Density", cex=0.7)

  fam_counts <- sort(table(is_track$fam), decreasing=TRUE)
  legend_fams <- names(fam_counts)[seq_len(min(12, length(fam_counts)))]
  for (i in seq_along(legend_fams)) {
    y_pos <- 0.2 - (i-1)*0.045
    graphics::rect(-0.35, y_pos-0.015, -0.28, y_pos+0.015, col=fam_cols[legend_fams[i]], border=NA)
    graphics::text(-0.26, y_pos, paste0(legend_fams[i]," (",fam_counts[legend_fams[i]],")"), cex=0.55, adj=c(0,0.5))
  }
  graphics::text(0.15, 0.2, "Essentiality:", cex=0.65, font=2)
  for (i in seq_along(c("#81C784","#FFF9C4","#FB8C00","#E53935"))) {
    y_pos <- 0.15 - (i-1)*0.04
    graphics::rect(0.15, y_pos-0.013, 0.20, y_pos+0.013, col=c("#81C784","#FFF9C4","#FB8C00","#E53935")[i], border=NA)
    graphics::text(0.22, y_pos, c("Low","Medium","High","Essential")[i], cex=0.55, adj=c(0,0.5))
  }
  graphics::text(0.15, -0.05, "LP Score:", cex=0.65, font=2)
  for (i in seq_along(c("#EF5350","#FF8F00","#FDD835","#66BB6A","#2E7D32"))) {
    y_pos <- -0.10 - (i-1)*0.04
    graphics::rect(0.15, y_pos-0.013, 0.20, y_pos+0.013, col=c("#EF5350","#FF8F00","#FDD835","#66BB6A","#2E7D32")[i], border=NA)
    graphics::text(0.22, y_pos, c("<0.65","0.65-0.72","0.72-0.78","0.78-0.82",">0.82")[i], cex=0.55, adj=c(0,0.5))
  }
  circlize::circos.clear()

  # ---- Panel 2: Census ----
  graphics::par(mar=c(3,5,3,1))
  census_ordered <- census[order(census$element_count), , drop=FALSE]
  bar_labels <- ifelse(!is.na(census_ordered$recognition_motif) & nzchar(census_ordered$recognition_motif),
                       paste0(census_ordered$family, " [", census_ordered$recognition_motif, "]"),
                       census_ordered$family)
  graphics::barplot(
    rbind(census_ordered$low_confidence_count, census_ordered$medium_confidence_count, census_ordered$high_confidence_count),
    beside=FALSE, horiz=TRUE, names.arg=bar_labels,
    col=c("#FFC107","#42A5F5","#66BB6A"), border=NA, las=1, cex.names=0.7, cex.axis=0.7,
    xlab="Count", main="B) IS Census")
  graphics::legend("bottomright", legend=c("Low","Medium","High"),
                   fill=c("#FFC107","#42A5F5","#66BB6A"), cex=0.65, bty="n", title="Confidence")

  # ---- Panel 3: Recognition sequences ----
  graphics::par(mar=c(3,7,3,1))
  tm <- target_models[target_models$n_elements >= 2, , drop=FALSE]
  tm <- tm[order(tm$n_elements), , drop=FALSE]
  if (nrow(tm)) {
    plot(0, 0, type="n", xlim=c(0,3), ylim=c(0.5,nrow(tm)+0.5), axes=FALSE, xlab="", ylab="", main="C) Recognition Sequences")
    for (i in seq_len(nrow(tm))) {
      conf_col <- switch(tm$model_confidence[i], high="#66BB6A", medium="#42A5F5", low="#FFC107", "#E0E0E0")
      graphics::rect(0, i-0.35, 0.3, i+0.35, col=conf_col, border=NA)
      motif <- if (!is.na(tm$model_motif[i])) tm$model_motif[i] else "-"
      tsd_info <- if (!is.na(tm$dominant_tsd_len[i])) paste0("TSD: ",tm$dominant_tsd_len[i],"bp") else "TSD: -"
      graphics::text(0.4, i, motif, adj=c(0,0.5), cex=0.8, font=2)
      graphics::text(1.6, i, tsd_info, adj=c(0,0.5), cex=0.65, col="grey30")
    }
    graphics::axis(2, at=seq_len(nrow(tm)), labels=paste0(tm$family," (n=",tm$n_elements,")"),
                   las=1, cex.axis=0.7, tick=FALSE, line=-0.5)
    graphics::legend("bottomright", legend=c("High","Medium","Low"),
                     fill=c("#66BB6A","#42A5F5","#FFC107"), cex=0.6, bty="n", title="Model Conf.")
  }

  # ---- Panel 4: Top landing pads table ----
  graphics::par(mar=c(2,1,3,1))
  top_n <- min(10L, nrow(lp_track))
  top_lp <- utils::head(lp_track, top_n)
  plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,top_n+2), axes=FALSE, xlab="", ylab="", main="D) Top Landing Pads")
  graphics::text(0, top_n+1.5, "Rank", adj=c(0,0.5), cex=0.6, font=2)
  graphics::text(1, top_n+1.5, "Position", adj=c(0,0.5), cex=0.6, font=2)
  graphics::text(3, top_n+1.5, "Size", adj=c(0,0.5), cex=0.6, font=2)
  graphics::text(4.2, top_n+1.5, "Score", adj=c(0,0.5), cex=0.6, font=2)
  graphics::text(5.5, top_n+1.5, "Flanking genes", adj=c(0,0.5), cex=0.6, font=2)
  graphics::segments(0, top_n+1.1, 10, top_n+1.1, lwd=0.5)
  for (i in seq_len(top_n)) {
    y <- top_n + 1 - i
    sc <- top_lp$landing_pad_score[i]
    bg_col <- if (sc >= 0.82) "#E8F5E9" else if (sc >= 0.78) "#F1F8E9" else "#FFFDE7"
    graphics::rect(-0.2, y-0.4, 10.2, y+0.4, col=bg_col, border=NA)
    graphics::text(0, y, top_lp$landing_pad_id[i], adj=c(0,0.5), cex=0.55)
    graphics::text(1, y, paste0(round(top_lp$region_start[i]/1000)," kb"), adj=c(0,0.5), cex=0.55)
    graphics::text(3, y, paste0(top_lp$region_size_bp[i]," bp"), adj=c(0,0.5), cex=0.55)
    graphics::text(4.2, y, sprintf("%.3f", sc), adj=c(0,0.5), cex=0.55, font=2)
    left_p <- substr(top_lp$left_gene_product[i], 1, 18)
    right_p <- substr(top_lp$right_gene_product[i], 1, 18)
    graphics::text(5.5, y, paste0(left_p, " | ", right_p), adj=c(0,0.5), cex=0.45)
  }

  grDevices::dev.off()
  list(pdf = pdf_path, comprehensive = TRUE)
}

# ---------------------------------------------------------------------------
# IS Census Overview: Family-wise element counts with recognition info
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Landing Pad Genome Map: Linear genome with landing pads, IS elements, essentiality
# ---------------------------------------------------------------------------
.dnmb_plot_iselement_landing_pads <- function(genbank_table, output_dir) {
  lp_path <- base::file.path(.dnmb_iselement_plot_module_dir(output_dir), "iselement_landing_pads.tsv")
  if (!base::file.exists(lp_path)) {
    return(NULL)
  }
  landing_pads <- utils::read.delim(lp_path, check.names = FALSE)
  if (!base::nrow(landing_pads)) {
    return(NULL)
  }

  tbl <- .dnmb_contig_ordered_table(genbank_table)
  if (!base::nrow(tbl)) {
    return(NULL)
  }

  # Get contig lengths from data
  contig_info <- tbl %>%
    dplyr::group_by(.data$contig) %>%
    dplyr::summarise(contig_len = max(.data$end, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$contig_len))

  # Filter to main contigs
  total_len <- sum(contig_info$contig_len)
  contig_info <- contig_info %>%
    dplyr::mutate(frac = .data$contig_len / total_len) %>%
    dplyr::filter(.data$frac >= 0.05 | dplyr::row_number() <= 5)
  main_contigs <- contig_info$contig

  landing_pads <- landing_pads %>% dplyr::filter(.data$contig %in% main_contigs)
  if (!nrow(landing_pads)) return(NULL)

  # Essentiality track from genbank_table
  ess_data <- tbl %>%
    dplyr::filter(.data$contig %in% main_contigs) %>%
    dplyr::select(dplyr::any_of(c("contig", "start", "end", "midpoint", "ISelement_essentiality_class"))) %>%
    dplyr::filter(!is.na(.data$midpoint))
  if ("ISelement_essentiality_class" %in% names(ess_data)) {
    ess_data$ess_class <- ess_data$ISelement_essentiality_class
  } else {
    ess_data$ess_class <- "unknown"
  }

  # IS element positions from genbank_table
  is_data <- tbl %>%
    dplyr::filter(.data$contig %in% main_contigs) %>%
    dplyr::filter(!is.na(.data$ISelement_element_id) & nzchar(.data$ISelement_element_id))

  plots <- list()

  for (ctg in main_contigs) {
    ctg_len <- contig_info$contig_len[contig_info$contig == ctg]
    lp_ctg <- landing_pads %>% dplyr::filter(.data$contig == ctg)
    is_ctg <- is_data %>% dplyr::filter(.data$contig == ctg)
    ess_ctg <- ess_data %>% dplyr::filter(.data$contig == ctg)

    # Build multi-track plot
    track_data <- list()

    # Track 1: Landing pads
    if (nrow(lp_ctg)) {
      lp_track <- lp_ctg %>%
        dplyr::mutate(
          track = "Landing Pads",
          ymin = 2, ymax = 3,
          fill_score = .data$landing_pad_score
        )
      track_data$lp <- lp_track
    }

    # Track 2: IS elements
    if (nrow(is_ctg)) {
      is_track <- is_ctg %>%
        dplyr::mutate(
          track = "IS Elements",
          ymin = 0.5, ymax = 1.5
        )
      track_data$is <- is_track
    }

    p <- ggplot2::ggplot() +
      ggplot2::xlim(0, ctg_len)

    # Draw landing pads
    if (!is.null(track_data$lp)) {
      p <- p + ggplot2::geom_rect(
        data = track_data$lp,
        ggplot2::aes(
          xmin = .data$region_start, xmax = .data$region_end,
          ymin = .data$ymin, ymax = .data$ymax,
          fill = .data$fill_score
        ),
        alpha = 0.9
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#EF5350", mid = "#FFC107", high = "#66BB6A",
        midpoint = 0.55, limits = c(0, 1),
        name = "Landing Pad\nScore"
      )
    }

    # Draw IS elements
    if (!is.null(track_data$is)) {
      family_col_name <- if ("ISelement_element_type" %in% names(is_ctg)) "ISelement_element_type" else "ISelement_feature_type"
      is_ctg$is_family <- if (family_col_name %in% names(is_ctg)) as.character(is_ctg[[family_col_name]]) else "IS_element"
      p <- p + ggplot2::geom_rect(
        data = is_ctg,
        ggplot2::aes(
          xmin = .data$start, xmax = .data$end,
          ymin = 0.5, ymax = 1.5
        ),
        fill = "#7E57C2", alpha = 0.8
      )
    }

    # Track labels
    p <- p +
      ggplot2::annotate("text", x = 0, y = 2.5, label = "Landing Pads", hjust = 1.1, size = 3, fontface = "bold") +
      ggplot2::annotate("text", x = 0, y = 1.0, label = "IS Elements", hjust = 1.1, size = 3, fontface = "bold") +
      ggplot2::labs(
        title = paste0("Landing Pad Map: ", ctg),
        x = "Genome coordinate (bp)"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", size = 11)
      ) +
      ggplot2::coord_cartesian(ylim = c(-0.5, 4))

    plots[[ctg]] <- p
  }

  if (!length(plots)) return(NULL)

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- base::file.path(plot_dir, "iselement_landing_pad_map.pdf")
  n_contigs <- length(plots)
  plot_height <- max(4, n_contigs * 3.5)
  .dnmb_module_plot_save_multi(plots, pdf_path, labels = NULL, ncol = 1,
                                width = 16, height = plot_height)
  list(pdf = pdf_path)
}

# --- Placeholder plot when no prophage regions are detected ---
#' Internal ISelement plotting helpers
#'
#' Plot-construction and fallback helpers for ISelement overview figures.
#'
#' @name dnmb_internal_plot_iselement
#' @keywords internal
#' @noRd
NULL

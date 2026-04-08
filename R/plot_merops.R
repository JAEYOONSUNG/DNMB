.dnmb_merops_class_from_family <- function(values) {
  values <- as.character(values)
  class_id <- substr(values, 1, 1)
  class_id[!class_id %in% c("M", "C", "S", "A")] <- "other"
  class_id
}

.dnmb_merops_class_palette <- function() {
  c(
    M = "#22C55E",
    C = "#EF4444",
    S = "#0EA5E9",
    A = "#FB923C",
    other = "#A1A1AA"
  )
}

.dnmb_filter_merops_plot_hits <- function(tbl,
                                          evalue_threshold = 1e-10,
                                          pident_threshold = 30,
                                          qcov_threshold = 0.3) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!"MEROPS_evalue" %in% names(tbl)) tbl$MEROPS_evalue <- NA_real_
  if (!"MEROPS_pident" %in% names(tbl)) tbl$MEROPS_pident <- NA_real_
  if (!"MEROPS_qcov" %in% names(tbl)) tbl$MEROPS_qcov <- NA_real_
  tbl$MEROPS_evalue <- suppressWarnings(as.numeric(tbl$MEROPS_evalue))
  tbl$MEROPS_pident <- suppressWarnings(as.numeric(tbl$MEROPS_pident))
  tbl$MEROPS_qcov <- suppressWarnings(as.numeric(tbl$MEROPS_qcov))
  keep <- !is.na(tbl$MEROPS_family_id) & nzchar(tbl$MEROPS_family_id) &
    !is.na(tbl$MEROPS_evalue) & tbl$MEROPS_evalue <= evalue_threshold
  # Apply pident filter if data available
  has_pident <- !all(is.na(tbl$MEROPS_pident))
  if (has_pident) keep <- keep & (is.na(tbl$MEROPS_pident) | tbl$MEROPS_pident >= pident_threshold)
  # Apply qcov filter if data available
  has_qcov <- !all(is.na(tbl$MEROPS_qcov))
  if (has_qcov) keep <- keep & (is.na(tbl$MEROPS_qcov) | tbl$MEROPS_qcov >= qcov_threshold)
  tbl <- tbl[keep, , drop = FALSE]
  if (!nrow(tbl)) return(tbl)
  # Mark confirmed vs unassigned subfamilies
  # 1) MEROPS hit_label based: assigned subfamily = confirmed
  if ("MEROPS_hit_label" %in% names(tbl)) {
    tbl$merops_confirmed <- !grepl("unassigned|UPW|UPB|UNB|UPA|UPC|UPD|UNA|UNC|non-peptidase", tbl$MEROPS_hit_label, ignore.case = TRUE)
  } else {
    tbl$merops_confirmed <- rep(TRUE, nrow(tbl))
  }
  # 2) Product keyword rescue: if product contains protease/peptidase keywords, confirm it
  if ("product" %in% names(tbl)) {
    product_match <- grepl("peptidase|protease|proteinase|endopeptidase|exopeptidase|carboxypeptidase|aminopeptidase|dipeptidase", tbl$product, ignore.case = TRUE)
    tbl$merops_confirmed <- tbl$merops_confirmed | product_match
  }
  # 3) Known false positive demotion: demote confirmed hits whose product
  #    indicates non-protease function (shared domain homologues)
  if ("product" %in% names(tbl)) {
    fp_product <- grepl(
      "synthase|synthetase|transferase|transaminase|dehydrogenase|kinase|reductase|polymerase|DNA repair|hydratase|isomerase|amidotransferase",
      tbl$product, ignore.case = TRUE
    )
    # Only demote if product does NOT also contain peptidase/protease keywords
    is_protease_product <- grepl("peptidase|protease|proteinase", tbl$product, ignore.case = TRUE)
    tbl$merops_confirmed[fp_product & !is_protease_product] <- FALSE
  }
  # 4) Known non-protease MEROPS families: I87 (inhibitor), S09X (alpha/beta hydrolase fold)
  if ("MEROPS_family_id" %in% names(tbl)) {
    non_protease_fam <- tbl$MEROPS_family_id %in% c("I87")
    s09x_hydrolase <- tbl$MEROPS_family_id == "S09X" &
      grepl("hydrolase|esterase|lipase", tbl$product, ignore.case = TRUE) &
      !grepl("peptidase|protease", tbl$product, ignore.case = TRUE)
    tbl$merops_confirmed[non_protease_fam | s09x_hydrolase] <- FALSE
  }
  rownames(tbl) <- NULL
  tbl
}

.dnmb_merops_track_levels <- function(values) {
  classes <- unique(.dnmb_merops_class_from_family(values))
  levels <- c("M", "C", "S", "A", "other")
  levels[levels %in% classes]
}

.dnmb_genome_axis_breaks <- function(genome_length) {
  genome_length <- as.numeric(genome_length)[1]
  if (is.na(genome_length) || genome_length <= 0) {
    return(list(brk = numeric(), labels = character()))
  }

  step_size <- genome_length / 24
  if (genome_length >= 1e6) {
    label_unit <- "Mb"
    scale_factor <- 1e6
  } else {
    label_unit <- "kb"
    scale_factor <- 1e3
  }
  brk <- pretty(c(0, genome_length), n = max(4, floor(genome_length / max(step_size, 1))))
  brk <- brk[brk >= 0 & brk <= genome_length]
  labels <- paste0(round(brk / scale_factor, 2), label_unit)
  list(brk = brk, labels = labels)
}

.dnmb_merops_label_subset <- function(tbl, max_labels_per_class = 12L) {
  tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if ("MEROPS_bitscore" %in% names(tbl)) {
    tbl$MEROPS_bitscore <- suppressWarnings(as.numeric(tbl$MEROPS_bitscore))
  } else {
    tbl$MEROPS_bitscore <- NA_real_
  }
  tbl$MEROPS_evalue <- suppressWarnings(as.numeric(tbl$MEROPS_evalue))
  tbl$label_rank <- seq_len(nrow(tbl))
  out <- lapply(split(tbl, tbl$merops_class), function(df) {
    df <- df[order(df$MEROPS_evalue, -df$MEROPS_bitscore, df$label_rank), , drop = FALSE]
    utils::head(df, max_labels_per_class)
  })
  out <- dplyr::bind_rows(out)
  rownames(out) <- NULL
  out
}

.dnmb_truncate_label_text <- function(text, width = 56L) {
  text <- as.character(text)
  width <- as.integer(width)[1]
  if (is.na(width) || width < 8L) {
    return(text)
  }
  long <- !is.na(text) & nchar(text, type = "width") > width
  text[long] <- paste0(substr(text[long], 1, width - 3L), "...")
  text
}

.dnmb_merops_class_label <- function(family_id) {
  family_id <- as.character(family_id)
  class_id <- substr(family_id, 1, 1)
  label <- c(
    A = "aspartic peptidase",
    C = "cysteine peptidase",
    G = "glutamic peptidase",
    I = "peptidase inhibitor",
    M = "metallopeptidase",
    N = "asparagine peptide lyase",
    P = "mixed peptidase",
    S = "serine peptidase",
    T = "threonine peptidase",
    U = "unknown peptidase"
  )
  out <- unname(label[match(class_id, names(label))])
  out[is.na(out) | !nzchar(out)] <- "other peptidase"
  out
}

.dnmb_merops_family_tier_label <- function(family_id, hit_label) {
  family_id <- as.character(family_id)
  hit_label <- as.character(hit_label)
  extract_first <- function(pattern, x) {
    vapply(x, function(one) {
      if (is.na(one) || !nzchar(one)) {
        return(NA_character_)
      }
      m <- regexpr(pattern, one, perl = TRUE, ignore.case = TRUE)
      if (is.na(m) || m < 1L) {
        return(NA_character_)
      }
      regmatches(one, m)
    }, character(1))
  }

  tier <- extract_first("subfamily [A-Z0-9]+", hit_label)
  family_tier <- extract_first("family [A-Z0-9]+", hit_label)
  missing_tier <- is.na(tier) | !nzchar(tier)
  tier[missing_tier] <- family_tier[missing_tier]
  missing_tier <- is.na(tier) | !nzchar(tier)
  inferred <- ifelse(
    grepl("^[A-Z][0-9]+[A-Z]+$", family_id),
    family_id,
    ifelse(grepl("^[A-Z][0-9]+$", family_id), family_id, family_id)
  )
  tier[missing_tier] <- inferred[missing_tier]
  tier <- sub("^subfamily\\s+", "", tier, ignore.case = TRUE)
  tier <- sub("^family\\s+", "", tier, ignore.case = TRUE)
  tier
}

.dnmb_merops_label_text <- function(locus_tag, product, family_id, hit_label) {
  n <- max(length(locus_tag), length(product), length(family_id), length(hit_label))
  rep_len_chr <- function(x, n) rep_len(as.character(x), n)
  locus_tag <- rep_len_chr(locus_tag, n)
  product <- rep_len_chr(product, n)
  family_id <- rep_len_chr(family_id, n)
  hit_label <- rep_len_chr(hit_label, n)
  product <- .dnmb_truncate_label_text(gsub("\\s+", " ", trimws(product)), width = 56L)
  missing_product <- is.na(product) | !nzchar(product)
  tier_label <- .dnmb_merops_family_tier_label(family_id, hit_label)
  class_label <- .dnmb_merops_class_label(family_id)
  header <- paste0(tier_label, " | ", class_label, " | ", locus_tag)
  annotation <- product
  annotation[missing_product] <- "annotation unavailable"
  out <- paste(header, annotation, sep = "\n")
  out
}

.dnmb_read_merops_best_hits_for_plot <- function(output_dir) {
  blast_path <- file.path(output_dir, "dnmb_module_merops", "merops_blastp.tsv")
  if (!file.exists(blast_path)) {
    return(data.frame())
  }

  raw <- tryCatch(
    read.delim(blast_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(raw) || ncol(raw) < 12L || !nrow(raw)) {
    return(data.frame())
  }

  raw <- raw[, seq_len(min(13L, ncol(raw))), drop = FALSE]
  colnames(raw) <- c(
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    if (ncol(raw) >= 13L) "stitle" else NULL
  )
  raw$evalue <- suppressWarnings(as.numeric(raw$evalue))
  raw$bitscore <- suppressWarnings(as.numeric(raw$bitscore))
  raw$ord <- seq_len(nrow(raw))
  raw <- raw[order(raw$qseqid, raw$evalue, -raw$bitscore, raw$ord), , drop = FALSE]
  best <- raw[!duplicated(raw$qseqid), c("qseqid", "evalue", "bitscore"), drop = FALSE]
  colnames(best) <- c("locus_tag", "MEROPS_evalue", "MEROPS_bitscore")
  rownames(best) <- NULL
  best
}

.dnmb_plot_merops_module <- function(genbank_table, output_dir) {
  required <- c("locus_tag", "contig", "start", "end", "MEROPS_family_id")
  if (!all(required %in% names(genbank_table))) {
    return(NULL)
  }
  all_tbl <- as.data.frame(genbank_table, stringsAsFactors = FALSE)
  all_tbl$start <- suppressWarnings(as.numeric(all_tbl$start))
  all_tbl$end <- suppressWarnings(as.numeric(all_tbl$end))
  if (!"MEROPS_evalue" %in% names(all_tbl)) {
    all_tbl$MEROPS_evalue <- NA_real_
  }
  if (!"MEROPS_bitscore" %in% names(all_tbl)) {
    all_tbl$MEROPS_bitscore <- NA_real_
  }
  best_metrics <- .dnmb_read_merops_best_hits_for_plot(output_dir)
  if (nrow(best_metrics)) {
    idx <- match(all_tbl$locus_tag, best_metrics$locus_tag)
    fill_evalue <- is.na(all_tbl$MEROPS_evalue) & !is.na(idx)
    fill_bitscore <- is.na(all_tbl$MEROPS_bitscore) & !is.na(idx)
    all_tbl$MEROPS_evalue[fill_evalue] <- best_metrics$MEROPS_evalue[idx[fill_evalue]]
    all_tbl$MEROPS_bitscore[fill_bitscore] <- best_metrics$MEROPS_bitscore[idx[fill_bitscore]]
  }
  tbl <- all_tbl
  tbl$MEROPS_evalue <- suppressWarnings(as.numeric(tbl$MEROPS_evalue))
  tbl <- .dnmb_filter_merops_plot_hits(tbl)
  tbl <- tbl[!is.na(tbl$start) & !is.na(tbl$end), , drop = FALSE]
  if (!nrow(tbl)) {
    return(NULL)
  }

  contigs <- .dnmb_contig_lengths_for_plot(all_tbl, output_dir = output_dir)
  contigs <- contigs[!is.na(contigs$length_bp) & contigs$length_bp > 0, , drop = FALSE]
  contigs <- contigs[order(contigs$contig_number, contigs$contig), , drop = FALSE]
  class_palette <- .dnmb_merops_class_palette()
  tbl$merops_class <- .dnmb_merops_class_from_family(tbl$MEROPS_family_id)
  tbl$plot_color <- unname(class_palette[match(tbl$merops_class, names(class_palette))])
  # Thin labels so they have >= min_deg angular spacing per contig
  gap_sizes <- ifelse(contigs$length_bp > 100000, 5, 2)
  total_bp <- sum(contigs$length_bp)
  total_gap <- sum(gap_sizes)
  avail_deg <- 360 - total_gap
  min_deg <- 8
  label_tbl <- do.call(rbind, lapply(split(tbl, tbl$contig), function(ct) {
    clen <- contigs$length_bp[contigs$contig == ct$contig[1]]
    if (!length(clen) || clen == 0) return(ct)
    sector_deg <- (clen / total_bp) * avail_deg
    ct <- ct[order(ct$start), , drop = FALSE]
    # Only thin labels on small contigs where labels would overlap
    if (nrow(ct) <= 2 || sector_deg >= 30) return(ct)
    deg_per_bp <- sector_deg / clen
    # Greedily select labels with >= min_deg spacing; always keep first & last
    keep <- rep(FALSE, nrow(ct))
    keep[1] <- TRUE; keep[nrow(ct)] <- TRUE
    last_pos <- ct$start[1]
    for (i in seq(2, nrow(ct) - 1)) {
      if ((ct$start[i] - last_pos) * deg_per_bp >= min_deg) {
        keep[i] <- TRUE
        last_pos <- ct$start[i]
      }
    }
    ct[keep, , drop = FALSE]
  }))
  label_key <- paste(label_tbl$contig, label_tbl$start, label_tbl$end, label_tbl$locus_tag, sep = "::")
  tbl_key <- paste(tbl$contig, tbl$start, tbl$end, tbl$locus_tag, sep = "::")
  tbl$is_labeled <- tbl_key %in% label_key
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "Protease_overview.pdf")

  # Precompute bar colors
  tbl$neg_log_e <- pmin(50, -log10(pmax(1e-300, suppressWarnings(as.numeric(tbl$MEROPS_evalue)))))
  tbl$pident_val <- suppressWarnings(as.numeric(tbl$MEROPS_pident))
  tbl$pident_val[is.na(tbl$pident_val)] <- 50
  tbl$bar_color <- vapply(seq_len(nrow(tbl)), function(i) {
    base_rgb <- grDevices::col2rgb(tbl$plot_color[i]) / 255
    blend <- scales::rescale(tbl$pident_val[i], to = c(0.6, 1.0), from = c(30, 100))
    grDevices::rgb(base_rgb[1]*blend + 0.92*(1-blend), base_rgb[2]*blend + 0.92*(1-blend), base_rgb[3]*blend + 0.92*(1-blend))
  }, character(1))

  # Prepare label data (shared between circlize draws)
  label_bed <- label_tbl[, c("contig", "start", "end"), drop = FALSE]
  label_product <- if ("product" %in% names(label_tbl)) label_tbl$product else NA_character_
  label_hit <- if ("MEROPS_hit_label" %in% names(label_tbl)) label_tbl$MEROPS_hit_label else NA_character_
  label_text <- .dnmb_merops_label_text(label_tbl$locus_tag, label_product, label_tbl$MEROPS_family_id, label_hit)
  label_cols <- unname(class_palette[match(label_tbl$merops_class, names(class_palette))])
  label_cols[is.na(label_cols)] <- class_palette[["other"]]
  is_confirmed <- if ("merops_confirmed" %in% names(label_tbl)) label_tbl$merops_confirmed else rep(TRUE, nrow(label_tbl))
  label_text_cols <- vapply(seq_along(label_cols), function(i) {
    if (isTRUE(is_confirmed[i])) label_cols[i] else grDevices::adjustcolor(label_cols[i], alpha.f = 0.35)
  }, character(1))
  label_line_cols <- vapply(seq_along(label_cols), function(i) {
    if (isTRUE(is_confirmed[i])) label_cols[i] else grDevices::adjustcolor(label_cols[i], alpha.f = 0.25)
  }, character(1))

  # Family clustering links
  link_candidates <- tbl[order(tbl$contig, tbl$MEROPS_family_id, tbl$start), , drop = FALSE]
  link_candidates <- link_candidates |>
    dplyr::group_by(.data$contig, .data$MEROPS_family_id) |>
    dplyr::filter(dplyr::n() >= 3) |>
    dplyr::mutate(next_start = dplyr::lead(.data$start), next_end = dplyr::lead(.data$end)) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$next_start), abs(.data$next_start - .data$end) <= 50000)

  # --- Panel B: Family distribution by class (faceted vertical bars) ---
  family_detail <- tbl |>
    dplyr::mutate(status = ifelse(.data$merops_confirmed, "Confirmed", "Unassigned")) |>
    dplyr::count(.data$MEROPS_family_id, .data$merops_class, .data$status, name = "n") |>
    dplyr::arrange(.data$merops_class, dplyr::desc(.data$n))
  family_totals <- family_detail |>
    dplyr::group_by(.data$MEROPS_family_id, .data$merops_class) |>
    dplyr::summarise(total = sum(.data$n), confirmed = sum(.data$n[.data$status == "Confirmed"]), .groups = "drop") |>
    dplyr::arrange(.data$merops_class, dplyr::desc(.data$confirmed), dplyr::desc(.data$total))
  family_detail$merops_class[is.na(family_detail$merops_class)] <- "other"
  # Class summary rows (first bar in each group = class total)
  class_n <- family_totals |> dplyr::count(.data$merops_class, wt = .data$total, name = "class_total")
  class_conf <- family_totals |>
    dplyr::group_by(.data$merops_class) |>
    dplyr::summarise(class_confirmed = sum(.data$confirmed), .groups = "drop")
  class_summary <- dplyr::left_join(class_n, class_conf, by = "merops_class")
  # Create class summary rows as pseudo-families
  class_rows_conf <- data.frame(
    MEROPS_family_id = paste0("ALL_", class_summary$merops_class),
    merops_class = class_summary$merops_class,
    status = "Confirmed", n = class_summary$class_confirmed, stringsAsFactors = FALSE)
  class_rows_unassigned <- data.frame(
    MEROPS_family_id = paste0("ALL_", class_summary$merops_class),
    merops_class = class_summary$merops_class,
    status = "Unassigned", n = class_summary$class_total - class_summary$class_confirmed, stringsAsFactors = FALSE)
  class_rows <- rbind(class_rows_conf, class_rows_unassigned[class_rows_unassigned$n > 0, ])
  # Add class totals for labels
  class_total_labels <- data.frame(
    MEROPS_family_id = paste0("ALL_", class_summary$merops_class),
    merops_class = class_summary$merops_class,
    total = class_summary$class_total,
    confirmed = class_summary$class_confirmed, stringsAsFactors = FALSE)
  # Merge family + class summary
  all_detail <- rbind(family_detail, class_rows)
  all_totals <- rbind(family_totals, class_total_labels)
  # Factor levels: within each class, ALL_ first, then families by confirmed desc
  ordered_ids <- c()
  cls_order <- class_n$merops_class[order(-class_n$class_total)]
  # Move "other" to the end so it appears as the rightmost facet
  cls_order <- c(cls_order[cls_order != "other"], cls_order[cls_order == "other"])
  for (cls in cls_order) {
    cls_families <- family_totals$MEROPS_family_id[family_totals$merops_class == cls]
    ordered_ids <- c(ordered_ids, paste0("ALL_", cls), cls_families)
  }
  all_detail$MEROPS_family_id <- factor(all_detail$MEROPS_family_id, levels = ordered_ids)
  all_totals$MEROPS_family_id <- factor(all_totals$MEROPS_family_id, levels = ordered_ids)
  # Facet labels
  facet_labels <- stats::setNames(
    dplyr::case_when(class_n$merops_class == "M" ~ "Metallo (M)",
      class_n$merops_class == "C" ~ "Cysteine (C)",
      class_n$merops_class == "S" ~ "Serine (S)",
      class_n$merops_class == "A" ~ "Aspartic (A)",
      TRUE ~ "Other"),
    class_n$merops_class)
  all_detail$class_label <- facet_labels[all_detail$merops_class]
  all_totals$class_label <- facet_labels[all_totals$merops_class]
  all_detail$class_label <- factor(all_detail$class_label, levels = facet_labels[cls_order])
  all_totals$class_label <- factor(all_totals$class_label, levels = facet_labels[cls_order])
  # Fill palette
  bc_fill_vals <- c()
  for (cls in names(class_palette)) {
    bc_fill_vals[paste0(cls, "::Confirmed")] <- class_palette[[cls]]
    bc_fill_vals[paste0(cls, "::Unassigned")] <- grDevices::adjustcolor(class_palette[[cls]], alpha.f = 0.35)
  }
  all_detail$status <- factor(all_detail$status, levels = c("Confirmed", "Unassigned"))
  all_detail$fill_key <- paste0(all_detail$merops_class, "::", all_detail$status)
  # Display labels: strip "ALL_" prefix for class summary bars
  x_labels <- levels(all_detail$MEROPS_family_id)
  x_display <- ifelse(grepl("^ALL_", x_labels), "ALL", x_labels)

  p_bc <- ggplot2::ggplot(all_detail, ggplot2::aes(x = .data$MEROPS_family_id, y = .data$n, fill = .data$fill_key)) +
    ggplot2::geom_col(width = 0.7, position = ggplot2::position_stack(reverse = TRUE), show.legend = FALSE) +
    ggplot2::geom_text(data = all_totals,
      ggplot2::aes(x = .data$MEROPS_family_id, y = .data$total,
        label = ifelse(.data$confirmed == .data$total, as.character(.data$total),
          paste0(.data$confirmed, "/", .data$total))),
      inherit.aes = FALSE, vjust = -0.3, size = 2.2, color = "#4B5563") +
    ggplot2::scale_fill_manual(values = bc_fill_vals) +
    ggplot2::scale_x_discrete(labels = stats::setNames(x_display, x_labels)) +
    ggplot2::facet_grid(. ~ class_label, scales = "free_x", space = "free_x") +
    ggplot2::labs(title = "B. Family distribution by protease class", x = NULL, y = "Hits") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
      plot.title.position = "plot",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8, color = "#374151", angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9, color = "#374151"),
      axis.title.y = ggplot2::element_text(size = 11),
      strip.text.x = ggplot2::element_text(face = "bold", size = 10)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)))

  # Compose single-page PDF: circlize (Panel A, top 72%) + cowplot (B+C, bottom 28%)
  grDevices::pdf(pdf_path, width = 14, height = 18, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)

  # --- Panel A: circlize via par(fig) — top portion ---
  graphics::par(fig = c(0, 1, 0.30, 1), mar = c(0, 0, 2, 0), bg = "white")
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.after = gap_sizes,
    track.margin = c(0.01, 0.01), cell.padding = c(0, 0, 0, 0),
    points.overflow.warning = FALSE, canvas.xlim = c(-1.0, 1.0), canvas.ylim = c(-1.0, 1.0))
  circlize::circos.initialize(factors = contigs$contig, xlim = cbind(rep(0, nrow(contigs)), contigs$length_bp))
  # Axis
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.04,
    bg.border = "#94A3B8", bg.lwd = 0.3, bg.col = "#F8FAFC",
    panel.fun = function(x, y) {
      sector <- circlize::CELL_META$sector.index; max_x <- circlize::CELL_META$xlim[2]
      # Fixed 0.2Mb interval for large contigs; start+end only for small ones
      if (max_x >= 200000) {
        axis_brk <- seq(0, max_x, by = 200000)
        axis_brk <- axis_brk[axis_brk >= 0 & axis_brk <= max_x]
        axis_labels <- paste0(round(axis_brk / 1e6, 1), "Mb")
      } else {
        axis_brk <- c(0, max_x)
        scale_f <- if (max_x >= 1e6) 1e6 else 1e3
        unit <- if (max_x >= 1e6) "Mb" else "kb"
        axis_labels <- paste0(round(axis_brk / scale_f, 1), unit)
      }
      # Draw ticks only, then add labels vertically (clockwise facing)
      circlize::circos.axis(h = "top", major.at = axis_brk, labels = FALSE,
        col = "#94A3B8", lwd = 0.3)
      for (ai in seq_along(axis_brk)) {
        circlize::circos.text(x = axis_brk[ai], y = 1.3, labels = axis_labels[ai],
          facing = "clockwise", niceFacing = TRUE, cex = 0.7, col = "#64748B", adj = c(0, 0.5))
      }
      circlize::circos.text(x = max_x / 2, y = 0.4, labels = contigs$sector_label[match(sector, contigs$contig)],
        cex = 0.45, facing = "bending.inside", col = "#475569", font = 2)
    })
  # Spacer
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.08, bg.border = NA, bg.col = NA)
  # Labels
  circlize::circos.genomicLabels(bed = label_bed, labels = label_text, side = "outside", niceFacing = TRUE,
    col = label_text_cols, cex = 0.42, font = 2, padding = 0.15,
    connection_height = circlize::mm_h(14), line_col = label_line_cols, line_lwd = 0.4,
    labels_height = circlize::cm_h(3.2),
    track.margin = c(0.005, 0.005))
  # Bars
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.09,
    bg.border = "#CBD5E1", bg.lwd = 0.3, bg.col = "white",
    panel.fun = function(x, y) {
      sector <- circlize::CELL_META$sector.index
      hits <- tbl[tbl$contig == sector, , drop = FALSE]
      if (!nrow(hits)) return()
      for (i in seq_len(nrow(hits))) {
        bar_alpha <- if ("merops_confirmed" %in% names(hits) && !isTRUE(hits$merops_confirmed[i])) 0.3 else 0.85
        circlize::circos.rect(xleft = hits$start[i], ybottom = 0.05,
          xright = max(hits$end[i], hits$start[i] + circlize::CELL_META$xlim[2] * 0.001),
          ytop = 0.95, col = grDevices::adjustcolor(hits$plot_color[i], alpha.f = bar_alpha), border = NA)
      }
    })
  # Links (disabled — too few to be informative at this scale)
  graphics::title(main = "A. MEROPS protease distribution", cex.main = 1.5, col.main = "#1E293B", font.main = 2, adj = 0)
  graphics::legend("topleft", legend = c("Metallo (M)", "Cysteine (C)", "Serine (S)", "Aspartic (A)", "Other"),
    fill = unname(class_palette[c("M", "C", "S", "A", "other")]), border = NA, cex = 0.9, bty = "n",
    title = "Protease class", title.col = "#374151", title.font = 2)
  circlize::circos.clear()

  # --- Panel B: faceted family chart in bottom viewport ---
  print(p_bc, vp = grid::viewport(x = 0.5, y = 0.14, width = 0.97, height = 0.26))

  list(pdf = pdf_path)
}


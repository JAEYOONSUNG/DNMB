#' Extract organism-preferred TIR design rules from mRNAcal output
#'
#' Stratifies CDS by within-genome CAI (codon adaptation index) — a strong
#' abundance proxy validated against PaxDB on three species (Spearman 0.34-0.58)
#' — and contrasts the top decile (presumed highly-translated) against the
#' bottom decile (presumed lowly-translated). The differential profile reveals
#' the design choices the organism's translation machinery favours: SD motif
#' variant, spacer length, start codon, A/G content of the 70S IC footprint,
#' codon usage per amino acid, internal SD avoidance, etc.
#'
#' Useful for non-model bacteria where no proteomics data is available — the
#' genome's own highly-CAI gene set serves as the endogenous reference.
#'
#' @param results mRNAcal results table (data.frame from
#'   `dnmb_run_mrnacal_module(...)$results`) or path to
#'   `mrnacal_translation_efficiency.tsv`. Genbank-table merged variants
#'   (with `mRNAcal_*` prefixes) are also accepted.
#' @param top_quantile Top fraction of CAI to treat as preferred set
#'   (default 0.10 = top decile).
#' @param bottom_quantile Bottom fraction of CAI for anti-preferred set
#'   (default 0.10).
#' @param output_dir Where the design-rule outputs are written. Defaults to
#'   `mrnacal_design_rules/` in the current working directory.
#' @param organism Optional label for the title of the summary PDF.
#'
#' @return A list with components: `summary` (per-feature comparison),
#'   `codon_table` (preferred-codon table per AA), `rbs_pwm` (top vs bottom
#'   PWM over the RBS region) and file paths in `files`.
#'
#' @export
dnmb_extract_design_rules <- function(results,
                                      top_quantile = 0.10,
                                      bottom_quantile = 0.10,
                                      output_dir = NULL,
                                      organism = NULL) {
  results <- .dnmb_validate_load_results(results)
  if (!base::nrow(results)) {
    stop("mRNAcal results are empty.", call. = FALSE)
  }
  if (!"cai" %in% base::names(results)) {
    stop("`results` must contain a `cai` column.", call. = FALSE)
  }
  if (base::is.null(output_dir)) {
    output_dir <- base::file.path(base::getwd(), "mrnacal_design_rules")
  }
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cai <- suppressWarnings(base::as.numeric(results$cai))
  ok <- !base::is.na(cai)
  if (base::sum(ok) < 50L) {
    stop("Not enough CAI-scored genes (n<50) to extract design rules.", call. = FALSE)
  }
  q_top <- stats::quantile(cai[ok], probs = 1 - top_quantile, na.rm = TRUE, names = FALSE)
  q_bot <- stats::quantile(cai[ok], probs = bottom_quantile, na.rm = TRUE, names = FALSE)
  group <- base::rep(NA_character_, base::length(cai))
  group[ok & cai >= q_top] <- "top"
  group[ok & cai <= q_bot] <- "bottom"
  results$design_group <- group

  top_idx <- base::which(group == "top")
  bot_idx <- base::which(group == "bottom")
  if (!base::length(top_idx) || !base::length(bot_idx)) {
    stop("Design rule extraction failed: empty top or bottom group.", call. = FALSE)
  }

  summary_tbl <- .dnmb_design_feature_summary(results, top_idx, bot_idx)
  summary_path <- base::file.path(output_dir, "mrnacal_design_summary.tsv")
  readr::write_tsv(summary_tbl, summary_path)

  codon_table <- .dnmb_design_preferred_codons(results, top_idx)
  codon_path <- base::file.path(output_dir, "mrnacal_design_codon_table.tsv")
  readr::write_tsv(codon_table, codon_path)

  spacer_tbl <- .dnmb_design_spacer_distribution(results, top_idx, bot_idx)
  spacer_path <- base::file.path(output_dir, "mrnacal_design_spacer.tsv")
  readr::write_tsv(spacer_tbl, spacer_path)

  start_tbl <- .dnmb_design_start_distribution(results, top_idx, bot_idx)
  start_path <- base::file.path(output_dir, "mrnacal_design_start.tsv")
  readr::write_tsv(start_tbl, start_path)

  rbs_pwm <- .dnmb_design_rbs_pwm(results, top_idx, bot_idx)

  pdf_path <- base::file.path(output_dir, "mrnacal_design_summary.pdf")
  ok_plot <- tryCatch(
    .dnmb_design_summary_plot(
      pdf_path, summary_tbl, codon_table, spacer_tbl, start_tbl, rbs_pwm,
      organism = organism, n_top = base::length(top_idx), n_bot = base::length(bot_idx),
      cai_top = q_top, cai_bot = q_bot
    ),
    error = function(e) {
      base::warning("Design summary plot failed: ", conditionMessage(e), call. = FALSE)
      FALSE
    }
  )

  files <- list(
    summary = summary_path,
    codon_table = codon_path,
    spacer = spacer_path,
    start = start_path
  )
  if (base::isTRUE(ok_plot)) files$plot <- pdf_path

  list(
    summary = summary_tbl,
    codon_table = codon_table,
    spacer = spacer_tbl,
    start = start_tbl,
    rbs_pwm = rbs_pwm,
    n_top = base::length(top_idx),
    n_bottom = base::length(bot_idx),
    cai_top_threshold = q_top,
    cai_bottom_threshold = q_bot,
    files = files
  )
}

.dnmb_design_feature_summary <- function(results, top_idx, bot_idx) {
  candidates <- c(
    "rbs_score", "rbs_spacer", "rbs_mismatches",
    "duplex_score", "duplex_energy",
    "accessibility_score", "plfold_accessibility_score",
    "tir_plfold_unpaired_probability", "start_plfold_unpaired_probability",
    "downstream_plfold_unpaired_probability",
    "upstream_au_score", "upstream20_at_fraction",
    "tir_core_a_fraction", "tir_core_g_fraction", "tir_core_score",
    "ncs_at_fraction", "ncs_lysine_codon_count",
    "ncs_mfe_dg", "ncs_fold_score",
    "internal_sd_count", "internal_sd_penalty",
    "fold_mfe", "fold_mfe_per_nt", "fold_score",
    "early_k_score", "start_codon_score",
    "cai", "tai", "codon_efficiency_score"
  )
  feats <- base::intersect(candidates, base::names(results))
  rows <- list()
  for (f in feats) {
    x_top <- suppressWarnings(base::as.numeric(results[[f]][top_idx]))
    x_bot <- suppressWarnings(base::as.numeric(results[[f]][bot_idx]))
    x_top <- x_top[!base::is.na(x_top)]
    x_bot <- x_bot[!base::is.na(x_bot)]
    if (!base::length(x_top) || !base::length(x_bot)) next
    p <- tryCatch(
      suppressWarnings(stats::wilcox.test(x_top, x_bot)$p.value),
      error = function(e) NA_real_
    )
    fc <- if (base::abs(median(x_bot, na.rm = TRUE)) > 1e-9) {
      median(x_top, na.rm = TRUE) / median(x_bot, na.rm = TRUE)
    } else {
      NA_real_
    }
    rows[[base::length(rows) + 1L]] <- data.frame(
      feature = f,
      n_top = base::length(x_top),
      n_bottom = base::length(x_bot),
      top_median = base::round(median(x_top, na.rm = TRUE), 4),
      bottom_median = base::round(median(x_bot, na.rm = TRUE), 4),
      diff_top_minus_bot = base::round(
        median(x_top, na.rm = TRUE) - median(x_bot, na.rm = TRUE), 4
      ),
      ratio_top_over_bot = base::round(fc, 3),
      wilcox_p = base::signif(p, 3),
      stringsAsFactors = FALSE
    )
  }
  out <- if (base::length(rows)) dplyr::bind_rows(rows) else data.frame()
  if (base::nrow(out)) {
    out <- out[base::order(-base::abs(out$diff_top_minus_bot)), , drop = FALSE]
  }
  out
}

.dnmb_design_preferred_codons <- function(results, top_idx) {
  if (!"sequence_dna" %in% base::names(results)) {
    return(data.frame())
  }
  cds_seqs <- character()
  if ("rearranged_nt_seq" %in% base::names(results)) {
    cds_seqs <- base::as.character(results$rearranged_nt_seq[top_idx])
  } else if ("sequence_dna" %in% base::names(results)) {
    win_up <- suppressWarnings(base::as.integer(results$window_upstream[top_idx]))
    seq_full <- base::as.character(results$sequence_dna[top_idx])
    cds_seqs <- base::vapply(base::seq_along(seq_full), function(k) {
      up <- if (base::is.na(win_up[k])) 0L else win_up[k]
      base::substr(seq_full[k], up + 1L, base::nchar(seq_full[k]))
    }, character(1))
  }
  cds_seqs <- cds_seqs[!base::is.na(cds_seqs) & base::nzchar(cds_seqs)]
  if (!base::length(cds_seqs)) return(data.frame())
  total <- .dnmb_mrnacal_codon_counts_total(cds_seqs)
  code <- .dnmb_mrnacal_genetic_code()
  rows <- list()
  for (aa in base::unique(code)) {
    if (aa %in% c("M", "W", "*")) next
    family <- base::names(code)[code == aa]
    counts <- total[family]
    s <- base::sum(counts, na.rm = TRUE)
    if (s == 0L) next
    freq <- counts / s
    preferred <- base::names(counts)[base::which.max(counts)]
    rows[[base::length(rows) + 1L]] <- data.frame(
      aa = aa,
      family = base::paste(family, collapse = ","),
      preferred_codon = preferred,
      preferred_freq = base::round(freq[preferred], 3),
      counts = base::paste(base::sprintf("%s:%d", family, counts), collapse = ","),
      total = s,
      stringsAsFactors = FALSE
    )
  }
  if (base::length(rows)) dplyr::bind_rows(rows) else data.frame()
}

.dnmb_design_spacer_distribution <- function(results, top_idx, bot_idx) {
  if (!"rbs_spacer" %in% base::names(results)) return(data.frame())
  build <- function(idx, label) {
    sp <- suppressWarnings(base::as.integer(results$rbs_spacer[idx]))
    sp <- sp[!base::is.na(sp)]
    if (!base::length(sp)) return(data.frame())
    tab <- base::as.data.frame(base::table(sp), stringsAsFactors = FALSE)
    base::names(tab) <- c("spacer_nt", "n")
    tab$spacer_nt <- base::as.integer(tab$spacer_nt)
    tab$group <- label
    tab$fraction <- base::round(tab$n / base::sum(tab$n), 4)
    tab
  }
  dplyr::bind_rows(build(top_idx, "top"), build(bot_idx, "bottom"))
}

.dnmb_design_start_distribution <- function(results, top_idx, bot_idx) {
  if (!"start_codon" %in% base::names(results)) return(data.frame())
  build <- function(idx, label) {
    sc <- base::toupper(base::as.character(results$start_codon[idx]))
    sc <- sc[!base::is.na(sc) & base::nzchar(sc)]
    if (!base::length(sc)) return(data.frame())
    tab <- base::as.data.frame(base::table(sc), stringsAsFactors = FALSE)
    base::names(tab) <- c("start_codon", "n")
    tab$group <- label
    tab$fraction <- base::round(tab$n / base::sum(tab$n), 4)
    tab
  }
  dplyr::bind_rows(build(top_idx, "top"), build(bot_idx, "bottom"))
}

.dnmb_design_rbs_pwm <- function(results, top_idx, bot_idx, win_left = -20L, win_right = 5L) {
  if (!all(c("sequence_dna", "window_upstream") %in% base::names(results))) {
    return(NULL)
  }
  extract <- function(idx) {
    seqs <- character()
    for (k in idx) {
      seq_full <- base::as.character(results$sequence_dna[k])
      up <- suppressWarnings(base::as.integer(results$window_upstream[k]))
      if (base::is.na(up) || !base::nzchar(seq_full)) next
      start_coord <- up + 1L
      left <- start_coord + win_left
      right <- start_coord + win_right
      if (left < 1L || right > base::nchar(seq_full)) next
      seqs <- c(seqs, base::substr(seq_full, left, right))
    }
    seqs
  }
  build_pwm <- function(seqs) {
    if (!base::length(seqs)) return(NULL)
    L <- base::nchar(seqs[[1]])
    seqs <- seqs[base::nchar(seqs) == L]
    if (!base::length(seqs)) return(NULL)
    mat <- base::matrix(0, nrow = 4, ncol = L,
                        dimnames = list(c("A", "C", "G", "T"), base::seq.int(win_left, win_right)))
    for (s in seqs) {
      chars <- base::strsplit(s, "", fixed = TRUE)[[1]]
      for (j in base::seq_len(L)) {
        nt <- chars[j]
        if (nt %in% rownames(mat)) mat[nt, j] <- mat[nt, j] + 1L
      }
    }
    base::sweep(mat, 2, base::pmax(1L, base::colSums(mat)), "/")
  }
  list(top = build_pwm(extract(top_idx)), bottom = build_pwm(extract(bot_idx)),
       positions = base::seq.int(win_left, win_right))
}

.dnmb_design_summary_plot <- function(path, summary_tbl, codon_table, spacer_tbl,
                                      start_tbl, rbs_pwm, organism, n_top, n_bot,
                                      cai_top, cai_bot) {
  grDevices::pdf(path, width = 12, height = 13)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(base::matrix(c(1, 1, 2, 3, 4, 4, 5, 5, 6, 6), ncol = 2, byrow = TRUE),
                   heights = c(1.7, 1.3, 1.3, 1.7, 0.9))

  # Panel 1: PWM heatmap (top set)
  if (!base::is.null(rbs_pwm) && !base::is.null(rbs_pwm$top)) {
    graphics::par(mar = c(4, 4, 3, 1))
    mat <- rbs_pwm$top
    image_pal <- grDevices::colorRampPalette(c("#FFFFFF", "#0EA5E9", "#1E3A8A"))(40)
    graphics::image(
      x = base::as.integer(colnames(mat)),
      y = 1:4,
      z = base::t(mat),
      col = image_pal, zlim = c(0, 1),
      xlab = "Position relative to AUG",
      ylab = "",
      yaxt = "n",
      main = base::sprintf(
        "Preferred-set RBS region (top %d genes, CAI >= %.3f)%s",
        n_top, cai_top, if (!base::is.null(organism)) base::paste0(" — ", organism) else ""
      )
    )
    graphics::axis(2, at = 1:4, labels = c("A", "C", "G", "T"), las = 1)
    graphics::abline(v = 0, col = "#DC2626", lty = 2)
    for (i in 1:4) {
      for (j in base::seq_len(ncol(mat))) {
        v <- mat[i, j]
        if (v > 0.05) {
          graphics::text(
            x = base::as.integer(colnames(mat))[j], y = i,
            labels = base::sprintf("%.2f", v),
            cex = 0.55, col = if (v > 0.5) "white" else "#1F2937"
          )
        }
      }
    }
  } else {
    graphics::plot.new(); graphics::title("RBS PWM unavailable")
  }

  # Panel 2: spacer distribution
  graphics::par(mar = c(4, 4, 3, 1))
  if (base::nrow(spacer_tbl)) {
    sp_top <- spacer_tbl[spacer_tbl$group == "top", , drop = FALSE]
    sp_bot <- spacer_tbl[spacer_tbl$group == "bottom", , drop = FALSE]
    rng <- base::range(spacer_tbl$spacer_nt, na.rm = TRUE)
    pos <- base::seq.int(rng[1], rng[2])
    h_top <- stats::setNames(base::rep(0, base::length(pos)), pos)
    h_bot <- h_top
    h_top[base::as.character(sp_top$spacer_nt)] <- sp_top$fraction
    h_bot[base::as.character(sp_bot$spacer_nt)] <- sp_bot$fraction
    bar_mat <- base::rbind(top = h_top, bottom = h_bot)
    graphics::barplot(
      bar_mat, beside = TRUE, names.arg = pos,
      col = c("#15803D", "#B91C1C"),
      ylim = c(0, base::max(bar_mat) * 1.1),
      ylab = "Fraction", xlab = "SD-AUG spacer (nt)",
      main = "Spacer length: top vs bottom CAI"
    )
    graphics::legend("topright", legend = c("top", "bottom"),
                     fill = c("#15803D", "#B91C1C"), bty = "n")
  } else {
    graphics::plot.new(); graphics::title("spacer data missing")
  }

  # Panel 3: start codon distribution
  graphics::par(mar = c(4, 4, 3, 1))
  if (base::nrow(start_tbl)) {
    starts <- base::sort(base::unique(start_tbl$start_codon))
    sc_top <- stats::setNames(base::rep(0, base::length(starts)), starts)
    sc_bot <- sc_top
    sc_top[start_tbl$start_codon[start_tbl$group == "top"]] <- start_tbl$fraction[start_tbl$group == "top"]
    sc_bot[start_tbl$start_codon[start_tbl$group == "bottom"]] <- start_tbl$fraction[start_tbl$group == "bottom"]
    bar_mat <- base::rbind(top = sc_top, bottom = sc_bot)
    graphics::barplot(
      bar_mat, beside = TRUE, names.arg = starts,
      col = c("#15803D", "#B91C1C"),
      ylim = c(0, 1), ylab = "Fraction", xlab = "Start codon",
      main = "Start codon: top vs bottom CAI"
    )
    graphics::legend("topright", legend = c("top", "bottom"),
                     fill = c("#15803D", "#B91C1C"), bty = "n")
  } else {
    graphics::plot.new(); graphics::title("start codon data missing")
  }

  # Panel 4: top features by |diff|
  graphics::par(mar = c(4, 9, 3, 1))
  if (base::nrow(summary_tbl)) {
    top12 <- utils::head(summary_tbl, 12)
    cols <- ifelse(top12$diff_top_minus_bot >= 0, "#15803D", "#B91C1C")
    graphics::barplot(
      base::rev(top12$diff_top_minus_bot),
      names.arg = base::rev(top12$feature),
      horiz = TRUE, las = 1,
      col = base::rev(cols), border = NA,
      xlab = "median(top) − median(bottom)",
      main = "Differential features (top vs bottom CAI)"
    )
  } else {
    graphics::plot.new(); graphics::title("summary missing")
  }

  # Panel 5: codon table preferred summary (full width)
  graphics::par(mar = c(4, 4, 3, 1))
  if (base::nrow(codon_table)) {
    aa_order <- base::order(codon_table$aa)
    barvec <- codon_table$preferred_freq[aa_order]
    aa_lab <- base::paste0(codon_table$aa[aa_order], "\n", codon_table$preferred_codon[aa_order])
    bp <- graphics::barplot(
      barvec, names.arg = aa_lab, las = 1, col = "#0EA5E9",
      ylim = c(0, 1.05), ylab = "Top-codon freq within AA family",
      main = "Preferred codon per AA (top-CAI set)",
      cex.names = 0.85
    )
    graphics::abline(h = 0.5, col = "#9CA3AF", lty = 3)
    graphics::text(bp, barvec + 0.04, labels = base::sprintf("%.2f", barvec), cex = 0.7)
  } else {
    graphics::plot.new(); graphics::title("codon table missing")
  }

  # Panel 6: text legend / methodology
  graphics::par(mar = c(2, 2, 2, 2))
  graphics::plot.new()
  graphics::title(main = "Method", adj = 0, line = 0)
  msg <- base::paste0(
    "CAI is a validated abundance proxy (E. coli rho=0.57, B. subtilis 0.38, C. jejuni 0.35).\n",
    "Top set: CAI >= ", base::round(cai_top, 3), " (n=", n_top, ").\n",
    "Bottom set: CAI <= ", base::round(cai_bot, 3), " (n=", n_bot, ").\n",
    "Differential features that are positive (green) are enriched in the top set;\n",
    "negative (red) features are avoided. Use these as design hints for synthetic\n",
    "constructs targeted at this organism. No proteomics data needed."
  )
  graphics::text(0.02, 0.95, msg, adj = c(0, 1), cex = 0.85, family = "mono")
  TRUE
}

#' Validate mRNAcal scores against measured protein abundance
#'
#' Cross-checks the mRNAcal output against an external abundance reference,
#' typically PaxDB integrated whole-organism estimates. Reports per-feature
#' Spearman/Pearson correlations between abundance and mRNAcal columns and
#' renders a multi-panel scatter plot.
#'
#' @param results mRNAcal results table (data.frame) or path to
#'   `mrnacal_translation_efficiency.tsv`. May also be the
#'   genbank_table-style merged output containing `mRNAcal_*` columns.
#' @param abundance Optional data.frame with at least an id column and an
#'   `abundance` column, or a path to a TSV/CSV. When NULL the function
#'   tries to fetch PaxDB data using `taxid`.
#' @param taxid Optional NCBI taxonomy ID. When NULL the function tries to
#'   read it from `genbank` (looks for `/db_xref="taxon:NNNN"`).
#' @param genbank Optional GenBank file path used to detect taxid and
#'   build a locus_tag <-> protein_id map.
#' @param paxdb_url Optional direct URL to a PaxDB integrated TXT file.
#'   When provided, overrides the auto-constructed URL.
#' @param paxdb_dataset PaxDB dataset slug. Defaults to
#'   `WHOLE_ORGANISM-integrated`.
#' @param gene_metadata Optional data.frame with at least
#'   `locus_tag` and any of `gene`/`protein_id`/`old_locus_tag`,
#'   used to broaden ID matching when the mRNAcal results table
#'   lacks gene names. Falls back to the global `genbank_table`.
#' @param cache_dir Where to cache downloaded PaxDB files.
#' @param output_dir Directory where the validation TSVs and PDF go.
#' @param generate_plot Whether to render the scatter PDF.
#'
#' @return A list with `summary` (per-feature correlation table),
#'   `merged` (mRNAcal + abundance per gene), `taxid`, `paxdb_source`,
#'   and file paths in `files`.
#'
#' @export
dnmb_validate_mrnacal <- function(results,
                                  abundance = NULL,
                                  taxid = NULL,
                                  genbank = NULL,
                                  paxdb_url = NULL,
                                  paxdb_dataset = "WHOLE_ORGANISM-integrated",
                                  gene_metadata = NULL,
                                  cache_dir = NULL,
                                  output_dir = NULL,
                                  generate_plot = TRUE) {
  results <- .dnmb_validate_load_results(results)
  if (!base::nrow(results)) {
    stop("mRNAcal results are empty.", call. = FALSE)
  }
  if (base::is.null(gene_metadata)) {
    gene_metadata <- base::get0("genbank_table", envir = base::.GlobalEnv, inherits = FALSE)
  }
  if (base::is.data.frame(gene_metadata) && base::nrow(gene_metadata)) {
    add_cols <- base::intersect(
      c("gene", "protein_id", "old_locus_tag", "gene_synonym"),
      base::names(gene_metadata)
    )
    if (base::length(add_cols) && "locus_tag" %in% base::names(gene_metadata)) {
      meta <- base::as.data.frame(
        gene_metadata[, c("locus_tag", add_cols), drop = FALSE],
        stringsAsFactors = FALSE
      )
      meta$locus_tag <- base::as.character(meta$locus_tag)
      meta <- meta[!base::duplicated(meta$locus_tag), , drop = FALSE]
      idx <- base::match(base::as.character(results$locus_tag), meta$locus_tag)
      for (col in add_cols) {
        if (!col %in% base::names(results)) {
          results[[col]] <- meta[[col]][idx]
        }
      }
    }
  }
  if (base::is.null(output_dir)) {
    output_dir <- base::file.path(base::getwd(), "mrnacal_validation")
  }
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  paxdb_source <- "user_supplied"
  if (base::is.null(abundance)) {
    if (base::is.null(taxid) && !base::is.null(genbank)) {
      taxid <- .dnmb_validate_taxid_from_genbank(genbank)
    }
    if (base::is.null(taxid)) {
      stop("Provide either `abundance` or a `taxid` (or a `genbank` file containing /db_xref=\"taxon:NNN\").", call. = FALSE)
    }
    fetched <- .dnmb_paxdb_fetch(
      taxid = taxid,
      paxdb_url = paxdb_url,
      dataset = paxdb_dataset,
      cache_dir = cache_dir
    )
    if (base::is.null(fetched$file)) {
      stop(base::sprintf(
        "Could not retrieve PaxDB data for taxid %s (dataset %s). %s",
        taxid, paxdb_dataset, fetched$error %||% ""
      ), call. = FALSE)
    }
    abundance <- .dnmb_paxdb_parse(fetched$file)
    paxdb_source <- fetched$url %||% fetched$source
  } else {
    abundance <- .dnmb_validate_load_abundance(abundance)
  }
  if (!base::nrow(abundance)) {
    stop("Abundance table is empty.", call. = FALSE)
  }

  merged <- .dnmb_validate_match(results, abundance, genbank = genbank)
  if (!base::nrow(merged)) {
    stop("No mRNAcal genes matched the abundance table.", call. = FALSE)
  }
  merged_path <- base::file.path(output_dir, "mrnacal_validation_merged.tsv")
  readr::write_tsv(merged, merged_path)

  summary_tbl <- .dnmb_validate_correlations(merged)
  summary_path <- base::file.path(output_dir, "mrnacal_validation_summary.tsv")
  readr::write_tsv(summary_tbl, summary_path)

  files <- list(merged = merged_path, summary = summary_path)
  if (base::isTRUE(generate_plot)) {
    plot_path <- base::file.path(output_dir, "mrnacal_validation_scatter.pdf")
    ok <- tryCatch(
      .dnmb_validate_plot(merged, summary_tbl, plot_path, taxid = taxid, source = paxdb_source),
      error = function(e) {
        base::warning("Validation plot failed: ", conditionMessage(e), call. = FALSE)
        FALSE
      }
    )
    if (base::isTRUE(ok)) {
      files$plot <- plot_path
    }
  }

  list(
    summary = summary_tbl,
    merged = merged,
    taxid = taxid,
    paxdb_source = paxdb_source,
    files = files,
    n_matched = base::nrow(merged)
  )
}

.dnmb_validate_load_results <- function(x) {
  if (base::is.character(x) && base::length(x) == 1L && base::file.exists(x)) {
    x <- readr::read_tsv(x, show_col_types = FALSE)
  }
  if (!base::is.data.frame(x)) {
    stop("`results` must be a data.frame or a TSV path.", call. = FALSE)
  }
  x <- base::as.data.frame(x, stringsAsFactors = FALSE)
  if (!"locus_tag" %in% base::names(x)) {
    stop("`results` must contain a locus_tag column.", call. = FALSE)
  }
  if ("mRNAcal_tir_score" %in% base::names(x) && !"tir_score" %in% base::names(x)) {
    drop_prefix <- function(name) base::sub("^mRNAcal_", "", name)
    candidates <- base::grep("^mRNAcal_", base::names(x), value = TRUE)
    for (col in candidates) {
      bare <- drop_prefix(col)
      if (!bare %in% base::names(x)) {
        x[[bare]] <- x[[col]]
      }
    }
  }
  x
}

.dnmb_validate_load_abundance <- function(x) {
  if (base::is.character(x) && base::length(x) == 1L && base::file.exists(x)) {
    sep <- if (base::grepl("\\.csv$", x, ignore.case = TRUE)) "," else "\t"
    x <- utils::read.table(x, header = TRUE, sep = sep, stringsAsFactors = FALSE, comment.char = "#", quote = "\"")
  }
  if (!base::is.data.frame(x)) {
    stop("`abundance` must be a data.frame or a TSV/CSV path.", call. = FALSE)
  }
  x <- base::as.data.frame(x, stringsAsFactors = FALSE)
  ab_col <- base::intersect(c("abundance", "ppm", "value"), base::names(x))
  if (!base::length(ab_col)) {
    stop("`abundance` must contain an `abundance`, `ppm`, or `value` column.", call. = FALSE)
  }
  x$abundance <- suppressWarnings(base::as.numeric(x[[ab_col[[1]]]]))
  id_col <- base::intersect(
    c("locus_tag", "protein_id", "external_id", "string_id", "id", "name"),
    base::names(x)
  )
  if (!base::length(id_col)) {
    stop("`abundance` must contain at least one id column (locus_tag/protein_id/external_id/string_id/id/name).", call. = FALSE)
  }
  x$external_id <- base::as.character(x[[id_col[[1]]]])
  x[!base::is.na(x$abundance) & x$abundance > 0, , drop = FALSE]
}

.dnmb_validate_taxid_from_genbank <- function(genbank) {
  if (base::is.null(genbank) || !base::file.exists(genbank)) {
    return(NULL)
  }
  lines <- base::readLines(genbank, n = 5000L, warn = FALSE)
  hit <- base::grep("/db_xref=\"taxon:[0-9]+\"", lines, value = TRUE)
  if (!base::length(hit)) {
    return(NULL)
  }
  m <- stringr::str_match(hit[[1]], "taxon:([0-9]+)")[, 2]
  if (base::is.na(m) || !base::nzchar(m)) {
    return(NULL)
  }
  base::as.integer(m)
}

.dnmb_paxdb_fetch <- function(taxid,
                              paxdb_url = NULL,
                              dataset = "WHOLE_ORGANISM-integrated",
                              cache_dir = NULL) {
  if (base::is.null(cache_dir) || !base::nzchar(cache_dir)) {
    cache_dir <- tools::R_user_dir("DNMB", which = "cache")
  }
  base::dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- base::file.path(
    cache_dir,
    base::sprintf("paxdb_%s_%s.txt", taxid, dataset)
  )
  if (base::file.exists(cache_file) && base::file.size(cache_file) > 100L) {
    return(list(file = cache_file, url = paste0("cache:", cache_file), source = "cache"))
  }
  candidates <- character()
  if (!base::is.null(paxdb_url) && base::nzchar(paxdb_url)) {
    candidates <- c(candidates, paxdb_url)
  }
  candidates <- c(
    candidates,
    base::sprintf("https://pax-db.org/downloads/latest/datasets/%s/%s-%s.txt", taxid, taxid, dataset),
    base::sprintf("https://pax-db.org/downloads/6.1/datasets/%s/%s-%s.txt", taxid, taxid, dataset),
    base::sprintf("https://pax-db.org/downloads/6.0/datasets/%s/%s-%s.txt", taxid, taxid, dataset)
  )
  candidates <- base::unique(candidates)
  errors <- character()
  for (url in candidates) {
    ok <- tryCatch({
      utils::download.file(url, cache_file, quiet = TRUE, mode = "wb")
      base::file.exists(cache_file) && base::file.size(cache_file) > 100L
    }, error = function(e) {
      errors <<- c(errors, base::sprintf("%s: %s", url, conditionMessage(e)))
      FALSE
    }, warning = function(w) {
      errors <<- c(errors, base::sprintf("%s: %s", url, conditionMessage(w)))
      FALSE
    })
    if (base::isTRUE(ok)) {
      first_lines <- tryCatch(base::readLines(cache_file, n = 3L, warn = FALSE), error = function(e) character())
      looks_html <- base::any(base::grepl("<html|<!doctype", first_lines, ignore.case = TRUE))
      if (looks_html) {
        base::file.remove(cache_file)
        next
      }
      return(list(file = cache_file, url = url, source = "download"))
    }
  }
  list(file = NULL, error = base::paste(errors, collapse = " | "))
}

.dnmb_paxdb_parse <- function(file) {
  lines <- base::readLines(file, warn = FALSE)
  data_lines <- lines[!base::grepl("^#", lines) & base::nzchar(base::trimws(lines))]
  if (!base::length(data_lines)) {
    return(data.frame())
  }
  tbl <- tryCatch(
    utils::read.table(text = data_lines, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = ""),
    error = function(e) data.frame()
  )
  if (!base::nrow(tbl) || base::ncol(tbl) < 2L) {
    return(data.frame())
  }
  if (base::ncol(tbl) >= 3L) {
    base::names(tbl)[1:3] <- c("gene_name", "string_id", "abundance")
  } else {
    base::names(tbl)[1:2] <- c("string_id", "abundance")
    tbl$gene_name <- NA_character_
  }
  tbl$abundance <- suppressWarnings(base::as.numeric(tbl$abundance))
  tbl$string_id <- base::as.character(tbl$string_id)
  tbl$gene_name <- base::as.character(tbl$gene_name)
  tbl$external_id <- base::sub("^[0-9]+\\.", "", tbl$string_id)
  tbl <- tbl[!base::is.na(tbl$abundance) & tbl$abundance > 0, , drop = FALSE]
  tbl
}

.dnmb_validate_match <- function(results, abundance, genbank = NULL) {
  abundance <- base::as.data.frame(abundance, stringsAsFactors = FALSE)
  abundance$external_id <- base::as.character(abundance$external_id)
  if (!"gene_name" %in% base::names(abundance)) {
    abundance$gene_name <- NA_character_
  }
  abundance$gene_name <- base::as.character(abundance$gene_name)
  results <- base::as.data.frame(results, stringsAsFactors = FALSE)
  results$locus_tag <- base::as.character(results$locus_tag)
  if ("protein_id" %in% base::names(results)) {
    results$protein_id <- base::as.character(results$protein_id)
  }
  if ("gene" %in% base::names(results)) {
    results$gene <- base::as.character(results$gene)
  }

  norm_id <- function(x) {
    base::tolower(base::gsub("[^A-Za-z0-9]", "", base::as.character(x)))
  }

  ext_raw <- abundance$external_id
  ext_strip_ver <- base::sub("\\.[0-9]+$", "", ext_raw)
  ext_norm <- norm_id(ext_raw)
  gene_norm <- norm_id(abundance$gene_name)

  matches <- base::rep(NA_integer_, base::nrow(results))

  try_keys <- function(values) {
    if (!base::length(values)) return(invisible(NULL))
    raw <- base::as.character(values)
    raw_strip <- base::sub("\\.[0-9]+$", "", raw)
    norm <- norm_id(raw)
    candidates <- list(
      base::match(raw, ext_raw),
      base::match(raw, ext_strip_ver),
      base::match(raw_strip, ext_raw),
      base::match(raw_strip, ext_strip_ver),
      base::match(norm, ext_norm),
      if (base::any(base::nzchar(gene_norm), na.rm = TRUE)) base::match(norm, gene_norm) else integer(0)
    )
    for (m in candidates) {
      if (!base::length(m)) next
      take <- base::is.na(matches) & !base::is.na(m)
      matches[take] <<- m[take]
    }
  }

  if ("locus_tag" %in% base::names(results)) try_keys(results$locus_tag)
  if ("protein_id" %in% base::names(results)) try_keys(results$protein_id)
  if ("gene" %in% base::names(results)) try_keys(results$gene)
  if ("old_locus_tag" %in% base::names(results)) try_keys(results$old_locus_tag)
  if ("gene_synonym" %in% base::names(results)) try_keys(results$gene_synonym)

  keep <- !base::is.na(matches)
  if (!base::any(keep)) {
    return(results[0, , drop = FALSE])
  }
  out <- results[keep, , drop = FALSE]
  out$abundance <- abundance$abundance[matches[keep]]
  out$abundance_external_id <- abundance$external_id[matches[keep]]
  if ("gene_name" %in% base::names(abundance)) {
    out$abundance_gene_name <- abundance$gene_name[matches[keep]]
  }
  out$abundance_log10 <- base::log10(out$abundance)
  out
}

.dnmb_validate_correlations <- function(merged) {
  feature_candidates <- c(
    "tir_score", "tir_score_percentile",
    "rbs_score", "duplex_score", "duplex_energy",
    "accessibility_score", "plfold_accessibility_score",
    "rbs_plfold_unpaired_probability", "start_plfold_unpaired_probability",
    "tir_plfold_unpaired_probability", "downstream_plfold_unpaired_probability",
    "standby_plfold_unpaired_probability",
    "upstream_au_score", "upstream20_at_fraction",
    "start_codon_score", "early_k_score",
    "ncs_at_fraction", "ncs_lysine_codon_count",
    "tir_core_a_fraction", "tir_core_g_fraction",
    "internal_sd_count", "internal_sd_penalty",
    "fold_mfe", "fold_mfe_per_nt", "fold_score",
    "cai", "cai_score", "tai", "tai_score", "codon_efficiency_score"
  )
  feats <- base::intersect(feature_candidates, base::names(merged))
  if (!base::length(feats)) {
    return(data.frame())
  }
  abundance_log <- merged$abundance_log10
  rows <- vector("list", base::length(feats))
  for (i in base::seq_along(feats)) {
    f <- feats[[i]]
    x <- suppressWarnings(base::as.numeric(merged[[f]]))
    keep <- !base::is.na(x) & !base::is.na(abundance_log)
    n <- base::sum(keep)
    if (n < 5L) {
      rows[[i]] <- data.frame(feature = f, n = n, spearman_rho = NA_real_, spearman_p = NA_real_, pearson_r = NA_real_, pearson_p = NA_real_, stringsAsFactors = FALSE)
      next
    }
    sp <- suppressWarnings(stats::cor.test(x[keep], abundance_log[keep], method = "spearman", exact = FALSE))
    pe <- suppressWarnings(stats::cor.test(x[keep], abundance_log[keep], method = "pearson"))
    rows[[i]] <- data.frame(
      feature = f,
      n = n,
      spearman_rho = base::round(base::as.numeric(sp$estimate), 4),
      spearman_p = base::signif(base::as.numeric(sp$p.value), 3),
      pearson_r = base::round(base::as.numeric(pe$estimate), 4),
      pearson_p = base::signif(base::as.numeric(pe$p.value), 3),
      stringsAsFactors = FALSE
    )
  }
  out <- dplyr::bind_rows(rows)
  out <- out[base::order(-base::abs(out$spearman_rho)), , drop = FALSE]
  out
}

.dnmb_validate_plot <- function(merged, summary_tbl, path, taxid = NULL, source = NULL) {
  if (!base::nrow(merged)) {
    return(FALSE)
  }
  panel_features <- c(
    "tir_score", "codon_efficiency_score",
    "cai", "tai",
    "plfold_accessibility_score", "downstream_plfold_unpaired_probability"
  )
  panel_features <- base::intersect(panel_features, base::names(merged))
  if (!base::length(panel_features)) {
    return(FALSE)
  }
  grDevices::pdf(path, width = 10, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
  for (f in panel_features) {
    x <- suppressWarnings(base::as.numeric(merged[[f]]))
    y <- merged$abundance_log10
    keep <- !base::is.na(x) & !base::is.na(y)
    if (base::sum(keep) < 5L) {
      graphics::plot.new()
      graphics::title(main = base::paste0(f, " (insufficient data)"))
      next
    }
    rho <- suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
    graphics::plot(
      x[keep], y[keep],
      xlab = f, ylab = "log10(abundance)",
      pch = 16, cex = 0.4, col = grDevices::adjustcolor("steelblue", alpha.f = 0.5),
      main = base::sprintf("%s  rho=%.3f  n=%d", f, rho, base::sum(keep))
    )
    fit <- stats::lm(y[keep] ~ x[keep])
    graphics::abline(fit, col = "firebrick", lwd = 1.5)
  }
  graphics::par(mfrow = c(1, 1), mar = c(8, 4, 3, 1))
  cor_tbl <- summary_tbl[base::order(summary_tbl$spearman_rho), , drop = FALSE]
  cor_tbl <- cor_tbl[!base::is.na(cor_tbl$spearman_rho), , drop = FALSE]
  if (base::nrow(cor_tbl)) {
    cols <- ifelse(cor_tbl$spearman_rho >= 0, "steelblue", "firebrick")
    graphics::barplot(
      height = cor_tbl$spearman_rho,
      names.arg = cor_tbl$feature,
      las = 2,
      col = cols,
      ylim = base::range(c(-1, 1, cor_tbl$spearman_rho), na.rm = TRUE),
      ylab = "Spearman rho vs log10(abundance)",
      main = base::sprintf(
        "mRNAcal validation%s%s",
        if (!base::is.null(taxid)) base::paste0(" (taxid=", taxid, ")") else "",
        if (!base::is.null(source)) base::paste0("\n", source) else ""
      )
    )
    graphics::abline(h = 0, lty = 2, col = "grey50")
  }
  TRUE
}

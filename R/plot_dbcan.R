.dnmb_dbcan_plot_palette <- function() {
  c(
    CAZyme = "#087F8C",
    TC = "#D97706",
    TF = "#7C3AED",
    STP = "#DC2626",
    Sulfatase = "#2563EB",
    Peptidase = "#C026D3",
    other = "#9CA3AF",
    null = "#E5E7EB"
  )
}

.dnmb_dbcan_class_palette <- function() {
  c(
    GH = "#0072B2",
    GT = "#009E73",
    CE = "#D55E00",
    PL = "#CC79A7",
    AA = "#E69F00",
    CBM = "#56B4E9",
    other = "#7A7A7A"
  )
}

.dnmb_dbcan_evidence_palette <- function() {
  c(
    very_high = "#0B4F6C",
    high = "#2A7F9E",
    medium = "#74A9CF",
    low = "#F2B134",
    audit = "#D95F59",
    legacy = "#A0A0A0"
  )
}

.dnmb_dbcan_plot_text <- function(x, max_chars = 78L) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  long <- nchar(x) > max_chars
  x[long] <- paste0(substr(x[long], 1L, max_chars - 3L), "...")
  x
}

.dnmb_dbcan_expand_plot_families <- function(tbl, family_col) {
  if (is.null(family_col) || !family_col %in% names(tbl) || !nrow(tbl)) {
    return(tibble::tibble())
  }
  rows <- lapply(seq_len(nrow(tbl)), function(i) {
    families <- .dnmb_dbcan_family_tokens(tbl[[family_col]][[i]])
    if (!length(families)) return(NULL)
    data.frame(
      dbcan_gene_row = if ("dbcan_gene_row" %in% names(tbl)) tbl$dbcan_gene_row[[i]] else i,
      locus_tag = if ("locus_tag" %in% names(tbl)) as.character(tbl$locus_tag[[i]]) else NA_character_,
      family = families,
      class = sub("^([A-Za-z]+).*$", "\\1", families),
      evidence_tier = if ("dbcan_evidence_tier" %in% names(tbl)) as.character(tbl$dbcan_evidence_tier[[i]]) else "legacy",
      substrate = if ("dbcan_resolved_substrate" %in% names(tbl)) as.character(tbl$dbcan_resolved_substrate[[i]]) else NA_character_,
      substrate_source = if ("dbcan_substrate_source" %in% names(tbl)) as.character(tbl$dbcan_substrate_source[[i]]) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- tibble::as_tibble(dplyr::bind_rows(rows))
  if (!nrow(out)) return(out)
  out$class[!out$class %in% c("GH", "GT", "CE", "PL", "AA", "CBM")] <- "other"
  out$evidence_tier[is.na(out$evidence_tier) | !nzchar(out$evidence_tier)] <- "legacy"
  out$substrate[is.na(out$substrate) | !nzchar(out$substrate)] <- "unknown"
  out$substrate_source[is.na(out$substrate_source) | !nzchar(out$substrate_source)] <- "legacy"
  out$evidence_tier <- factor(out$evidence_tier, levels = c("very_high", "high", "medium", "low", "audit", "legacy"))
  out$class <- factor(out$class, levels = c("GH", "GT", "CE", "PL", "AA", "CBM", "other"))
  unique(out)
}

.dnmb_dbcan_read_plot_inputs <- function(genbank_table, output_dir) {
  tbl <- as.data.frame(genbank_table, stringsAsFactors = FALSE)
  if (!nrow(tbl)) return(list(genes = tbl, cgc = data.frame(), domains = data.frame()))
  row_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_gene_row", "dbcan_gene_row"))
  tbl$dbcan_gene_row <- if (!is.null(row_col)) {
    suppressWarnings(as.integer(tbl[[row_col]]))
  } else {
    seq_len(nrow(tbl))
  }

  family_col <- .dnmb_pick_column(
    tbl,
    c(
      "dbCAN_dbcan_all_families",
      "dbCAN_dbcan_hit",
      "dbCAN_family_id",
      "dbcan_all_families",
      "dbcan_hit",
      "family_id"
    )
  )
  evidence_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_evidence_tier", "dbcan_evidence_tier"))
  substrate_source_col <- .dnmb_pick_column(
    tbl,
    c("dbCAN_dbcan_substrate_source", "dbcan_substrate_source")
  )
  substrate_candidates <- c(
    "dbCAN_substrate_label",
    "dbCAN_dbcan_pul_substrate",
    "dbCAN_dbcan_sub_substrate",
    "dbCAN_dbcan_overview_substrate",
    "substrate_label",
    "dbcan_pul_substrate",
    "dbcan_sub_substrate",
    "dbcan_overview_substrate"
  )
  substrate_candidates <- substrate_candidates[substrate_candidates %in% names(tbl)]

  tbl$dbcan_all_families <- if (!is.null(family_col)) as.character(tbl[[family_col]]) else rep(NA_character_, nrow(tbl))
  tbl$dbcan_evidence_tier <- if (!is.null(evidence_col)) as.character(tbl[[evidence_col]]) else rep("legacy", nrow(tbl))
  tbl$dbcan_substrate_source <- if (!is.null(substrate_source_col)) {
    as.character(tbl[[substrate_source_col]])
  } else {
    rep(NA_character_, nrow(tbl))
  }
  tbl$dbcan_resolved_substrate <- rep(NA_character_, nrow(tbl))
  source_by_column <- c(
    dbCAN_dbcan_pul_substrate = "dbCAN-PUL",
    dbCAN_dbcan_sub_substrate = "dbCAN-sub",
    dbCAN_dbcan_overview_substrate = "overview",
    dbcan_pul_substrate = "dbCAN-PUL",
    dbcan_sub_substrate = "dbCAN-sub",
    dbcan_overview_substrate = "overview"
  )
  for (column_name in substrate_candidates) {
    value <- as.character(tbl[[column_name]])
    take <- (is.na(tbl$dbcan_resolved_substrate) | !nzchar(tbl$dbcan_resolved_substrate)) &
      !is.na(value) & nzchar(value) & value != "-"
    tbl$dbcan_resolved_substrate[take] <- value[take]
    inferred_source <- unname(source_by_column[column_name])
    if (!length(inferred_source) || is.na(inferred_source) || !nzchar(inferred_source)) {
      inferred_source <- "legacy"
    }
    source_missing <- is.na(tbl$dbcan_substrate_source) | !nzchar(tbl$dbcan_substrate_source)
    tbl$dbcan_substrate_source[take & source_missing] <- inferred_source
  }

  unique_match <- function(target, source) {
    target <- .dnmb_dbcan_clean_query_id(target)
    source <- .dnmb_dbcan_clean_query_id(source)
    source_unique <- !duplicated(source) & !duplicated(source, fromLast = TRUE) & !is.na(source)
    target_unique <- !duplicated(target) & !duplicated(target, fromLast = TRUE) & !is.na(target)
    out <- match(target, ifelse(source_unique, source, NA_character_))
    out[!target_unique] <- NA_integer_
    out
  }
  id_map_path <- file.path(output_dir, "dnmb_module_dbcan", "dbcan_gene_id_map.tsv")
  id_map <- if (file.exists(id_map_path)) {
    tryCatch(utils::read.delim(id_map_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  } else NULL

  overview_path <- file.path(output_dir, "dnmb_module_dbcan", "run_dbcan", "overview.tsv")
  if (file.exists(overview_path)) {
    overview <- dnmb_dbcan_parse_overview(overview_path, id_map = id_map)
    idx <- rep(NA_integer_, nrow(tbl))
    if ("dbcan_gene_row" %in% names(overview) && any(!is.na(overview$dbcan_gene_row)) &&
        any(!is.na(tbl$dbcan_gene_row))) {
      idx <- match(tbl$dbcan_gene_row, suppressWarnings(as.integer(overview$dbcan_gene_row)))
    } else {
      idx <- unique_match(tbl$locus_tag, overview$query)
      if ("protein_id" %in% names(tbl)) {
        unresolved <- is.na(idx)
        idx[unresolved] <- unique_match(tbl$protein_id, overview$query)[unresolved]
      }
    }
    mapped <- !is.na(idx)
    tbl$dbcan_all_families[mapped] <- overview$dbcan_all_families[idx[mapped]]
    tbl$dbcan_evidence_tier[mapped] <- overview$dbcan_evidence_tier[idx[mapped]]
    overview_substrate <- overview$dbcan_overview_substrate[idx]
    fill_substrate <- mapped & (
      is.na(tbl$dbcan_resolved_substrate) | !nzchar(tbl$dbcan_resolved_substrate) |
        tbl$dbcan_substrate_source == "family_prior"
    ) &
      !is.na(overview_substrate) & nzchar(overview_substrate)
    tbl$dbcan_resolved_substrate[fill_substrate] <- overview_substrate[fill_substrate]
    tbl$dbcan_substrate_source[fill_substrate] <- "overview"
  }

  tbl$dbcan_evidence_tier[is.na(tbl$dbcan_evidence_tier) | !nzchar(tbl$dbcan_evidence_tier)] <- "legacy"
  tbl$dbcan_resolved_substrate <- vapply(
    tbl$dbcan_resolved_substrate,
    .dnmb_dbcan_normalize_substrate,
    character(1)
  )
  tbl$dbcan_substrate_source[
    is.na(tbl$dbcan_substrate_source) | !nzchar(tbl$dbcan_substrate_source)
  ] <- "legacy"

  contig_key_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_contig_key", "dbcan_contig_key"))
  if (!is.null(contig_key_col)) {
    tbl$dbcan_contig_key <- as.character(tbl[[contig_key_col]])
  } else if ("contig_number" %in% names(tbl)) {
    number <- suppressWarnings(as.integer(tbl$contig_number))
    tbl$dbcan_contig_key <- ifelse(is.na(number), as.character(tbl$contig), sprintf("replicon_%04d", number))
  } else {
    tbl$dbcan_contig_key <- as.character(tbl$contig)
  }
  tbl$dbcan_contig_key[is.na(tbl$dbcan_contig_key) | !nzchar(tbl$dbcan_contig_key)] <- "unknown_replicon"

  cgc_id_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_id", "dbCAN_cgc_id", "dbcan_cgc_id"))
  cgc_type_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_gene_type", "dbCAN_cgc_gene_type", "dbcan_cgc_gene_type"))
  cgc_family_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_protein_family", "dbCAN_cgc_protein_family", "dbcan_cgc_protein_family"))

  if (!is.null(cgc_id_col) && !is.null(cgc_type_col)) {
    cgc <- tbl
    cgc$dbcan_cgc_id <- as.character(cgc[[cgc_id_col]])
    cgc$dbcan_cgc_gene_type <- as.character(cgc[[cgc_type_col]])
    cgc$dbcan_cgc_protein_family <- if (!is.null(cgc_family_col)) as.character(cgc[[cgc_family_col]]) else NA_character_
    cgc <- cgc[!is.na(cgc$dbcan_cgc_id) & nzchar(cgc$dbcan_cgc_id), , drop = FALSE]
  } else {
    cgc_path <- file.path(output_dir, "dnmb_module_dbcan", "run_dbcan", "cgc_standard_out.tsv")
    cgc <- data.frame()
    if (file.exists(cgc_path)) {
      cgc_raw <- tryCatch(dnmb_dbcan_parse_cgc_standard(cgc_path), error = function(e) NULL)
      if (is.data.frame(cgc_raw) && nrow(cgc_raw)) {
        if (is.data.frame(id_map) && nrow(id_map)) {
          contig_keys <- unique(as.character(id_map$dbcan_contig_key))
          contig_map <- tibble::tibble(
            contig = contig_keys,
            safe_contig = paste0("ctg", seq_along(contig_keys))
          )
          cgc_raw <- .dnmb_dbcan_restore_cgc_ids(cgc_raw, contig_map)
          cgc_raw <- .dnmb_dbcan_restore_query_map(cgc_raw, id_map)
        }
        if ("dbcan_gene_row" %in% names(cgc_raw) && any(!is.na(cgc_raw$dbcan_gene_row))) {
          idx <- match(suppressWarnings(as.integer(cgc_raw$dbcan_gene_row)), tbl$dbcan_gene_row)
        } else {
          idx <- match(unique_match(cgc_raw$query, tbl$locus_tag), seq_len(nrow(tbl)))
          if ("protein_id" %in% names(tbl)) {
            unresolved <- is.na(idx)
            idx[unresolved] <- match(
              unique_match(cgc_raw$query, tbl$protein_id)[unresolved],
              seq_len(nrow(tbl))
            )
          }
        }
        keep <- !is.na(idx)
        cgc <- tbl[idx[keep], , drop = FALSE]
        cgc$dbcan_cgc_id <- as.character(cgc_raw$dbcan_cgc_id[keep])
        cgc$dbcan_cgc_gene_type <- as.character(cgc_raw$dbcan_cgc_gene_type[keep])
        cgc$dbcan_cgc_protein_family <- as.character(cgc_raw$dbcan_cgc_protein_family[keep])
      }
    }
  }

  if (nrow(cgc)) {
    type <- trimws(as.character(cgc$dbcan_cgc_gene_type))
    type[tolower(type) %in% c("prodoric", "tf")] <- "TF"
    type[tolower(type) == "cazyme"] <- "CAZyme"
    type[tolower(type) == "tc"] <- "TC"
    type[tolower(type) == "stp"] <- "STP"
    type[tolower(type) == "sulfatase"] <- "Sulfatase"
    type[tolower(type) == "peptidase"] <- "Peptidase"
    type[tolower(type) %in% c("", "na", "null")] <- "null"
    type[!type %in% names(.dnmb_dbcan_plot_palette())] <- "other"
    cgc$gene_type <- type
    cgc$start <- suppressWarnings(as.numeric(cgc$start))
    cgc$end <- suppressWarnings(as.numeric(cgc$end))
    cgc <- cgc[is.finite(cgc$start) & is.finite(cgc$end), , drop = FALSE]
  }

  cazy <- tbl[!is.na(tbl$dbcan_all_families) & nzchar(tbl$dbcan_all_families), , drop = FALSE]
  cazy$start <- suppressWarnings(as.numeric(cazy$start))
  cazy$end <- suppressWarnings(as.numeric(cazy$end))
  cazy <- cazy[is.finite(cazy$start) & is.finite(cazy$end), , drop = FALSE]
  cazy$primary_family <- vapply(
    cazy$dbcan_all_families,
    function(value) .dnmb_dbcan_primary_family(.dnmb_dbcan_family_tokens(value)),
    character(1)
  )
  cazy$class <- sub("^([A-Za-z]+).*$", "\\1", cazy$primary_family)
  cazy$class[!cazy$class %in% c("GH", "GT", "CE", "PL", "AA", "CBM")] <- "other"
  cazy$class <- factor(cazy$class, levels = c("GH", "GT", "CE", "PL", "AA", "CBM", "other"))
  cazy$dbcan_evidence_tier <- factor(cazy$dbcan_evidence_tier, levels = c("very_high", "high", "medium", "low", "audit", "legacy"))
  domains <- .dnmb_dbcan_expand_plot_families(tbl, "dbcan_all_families")
  list(genes = tbl, cazy = cazy, cgc = cgc, domains = domains)
}

.dnmb_dbcan_contig_summary <- function(tbl) {
  if (!nrow(tbl)) return(tibble::tibble())
  display <- if ("contig" %in% names(tbl)) as.character(tbl$contig) else tbl$dbcan_contig_key
  data.frame(
    dbcan_contig_key = tbl$dbcan_contig_key,
    display = display,
    end = suppressWarnings(as.numeric(tbl$end)),
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(.data$dbcan_contig_key) |>
    dplyr::summarise(
      display = dplyr::first(.data$display),
      length_bp = max(.data$end, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      display = ifelse(is.na(.data$display) | !nzchar(.data$display), .data$dbcan_contig_key, .data$display),
      label = paste0(.dnmb_dbcan_plot_text(.data$display, 62L), " [", .data$dbcan_contig_key, "]")
    )
}

.dnmb_dbcan_hit_contig_groups <- function(inputs, contigs_per_page = 8L) {
  gene_order <- unique(as.character(inputs$genes$dbcan_contig_key))
  hit_keys <- unique(c(
    as.character(inputs$cazy$dbcan_contig_key),
    as.character(inputs$cgc$dbcan_contig_key)
  ))
  hit_keys <- hit_keys[!is.na(hit_keys) & nzchar(hit_keys)]
  keys <- c(gene_order[gene_order %in% hit_keys], setdiff(hit_keys, gene_order))
  if (!length(keys)) return(list())
  n_pages <- ceiling(length(keys) / max(1L, as.integer(contigs_per_page)))
  if (n_pages <= 1L) return(list(keys))
  unname(split(keys, cut(seq_along(keys), breaks = n_pages, labels = FALSE)))
}

.dnmb_dbcan_summary_page <- function(inputs, contig_keys = NULL,
                                     map_page = 1L, map_pages = 1L,
                                     include_details = TRUE,
                                     contig_ncol = 1L) {
  cazy <- inputs$cazy
  domains <- inputs$domains
  contigs <- .dnmb_dbcan_contig_summary(inputs$genes)
  if (!is.null(contig_keys)) {
    contigs <- contigs[contigs$dbcan_contig_key %in% contig_keys, , drop = FALSE]
    cazy <- cazy[cazy$dbcan_contig_key %in% contig_keys, , drop = FALSE]
  }
  class_palette <- .dnmb_dbcan_class_palette()
  evidence_palette <- .dnmb_dbcan_evidence_palette()

  contigs$contig_label <- factor(contigs$label, levels = contigs$label)
  cazy$midpoint <- (cazy$start + cazy$end) / 2
  cazy$contig_label <- factor(
    contigs$label[match(cazy$dbcan_contig_key, contigs$dbcan_contig_key)],
    levels = contigs$label
  )
  if (nrow(cazy)) {
    cazy <- cazy[!is.na(cazy$contig_label), , drop = FALSE]
  }
  p_map <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = contigs,
      ggplot2::aes(x = 0, xend = .data$length_bp, y = 0, yend = 0),
      linewidth = 0.8,
      color = "#C7CDD4"
    ) +
    ggplot2::facet_wrap(
      ~contig_label,
      ncol = max(1L, as.integer(contig_ncol)[1L]),
      scales = "free_x"
    ) +
    ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ggplot2::scale_y_continuous(limits = c(-0.25, 0.25), breaks = NULL) +
    ggplot2::labs(
      title = paste0(
        "A  Genome-wide CAZyme positions",
        if (map_pages > 1L) paste0("  ", map_page, "/", map_pages) else ""
      ),
      x = "Genome coordinate",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.title.position = "plot",
      strip.text = ggplot2::element_text(face = "bold", size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.key.height = grid::unit(0.28, "cm"),
      legend.spacing.x = grid::unit(0.10, "cm"),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 2)
    )
  if (nrow(cazy)) {
    p_map <- p_map +
      ggplot2::geom_rect(
        data = cazy,
        ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = -0.16, ymax = 0.16, fill = .data$class, alpha = .data$dbcan_evidence_tier),
        color = "#2F3B45",
        linewidth = 0.15
      ) +
      ggplot2::geom_point(
        data = cazy,
        ggplot2::aes(x = .data$midpoint, y = 0, fill = .data$class, alpha = .data$dbcan_evidence_tier),
        shape = 21,
        size = 2.2,
        color = "#263238",
        stroke = 0.2
      ) +
      ggplot2::scale_fill_manual(values = class_palette, drop = FALSE, name = "CAZy class") +
      ggplot2::scale_alpha_manual(values = c(very_high = 1, high = 0.9, medium = 0.72, low = 0.5, audit = 0.28, legacy = 0.65), drop = FALSE, name = "Evidence") +
      ggplot2::guides(
        alpha = ggplot2::guide_legend(nrow = 1L, byrow = TRUE, order = 1L),
        fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE, order = 2L)
      )
  }

  class_counts <- if (nrow(domains)) {
    domains |>
      dplyr::count(.data$class, .data$evidence_tier, name = "genes")
  } else tibble::tibble(class = character(), evidence_tier = character(), genes = integer())
  p_class <- ggplot2::ggplot(class_counts, ggplot2::aes(x = .data$class, y = .data$genes, fill = .data$evidence_tier)) +
    ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    ggplot2::labs(title = "B  CAZy evidence by class", x = NULL, y = "Gene-family calls") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.title.position = "plot",
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 2)
    )
  if (nrow(class_counts)) {
    p_class <- p_class + ggplot2::scale_fill_manual(values = evidence_palette, drop = FALSE, name = "Evidence")
  }

  substrate_counts <- if (nrow(domains)) {
    domains |>
      dplyr::filter(
        .data$substrate != "unknown",
        .data$substrate_source != "family_prior"
      ) |>
      tidyr::separate_rows("substrate", sep = ";\\s*") |>
      dplyr::count(.data$substrate, name = "genes", sort = TRUE) |>
      utils::head(15L)
  } else tibble::tibble(substrate = character(), genes = integer())
  if (nrow(substrate_counts)) {
    substrate_counts$substrate <- factor(substrate_counts$substrate, levels = rev(substrate_counts$substrate))
  }
  p_substrate <- ggplot2::ggplot(substrate_counts, ggplot2::aes(x = .data$genes, y = .data$substrate)) +
    ggplot2::geom_col(width = 0.68, fill = "#2A9D8F") +
    ggplot2::geom_text(ggplot2::aes(label = .data$genes), hjust = -0.2, size = 3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = "C  Supported substrates", x = "Gene-family calls", y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.title.position = "plot",
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 2)
    )

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(
        "dbCAN CAZyme/CGC overview  |  ",
        nrow(inputs$cazy), " mapped genes  |  ",
        length(unique(inputs$cgc$dbcan_cgc_id)), " CGCs"
      ),
      x = 0,
      hjust = 0,
      fontface = "bold",
      size = 15
    )
  if (!isTRUE(include_details)) {
    return(cowplot::plot_grid(title, p_map, ncol = 1, rel_heights = c(0.06, 0.94)))
  }
  bottom <- cowplot::plot_grid(p_class, p_substrate, ncol = 2, rel_widths = c(1, 1.15))
  page <- cowplot::plot_grid(p_map, bottom, ncol = 1, rel_heights = c(0.58, 0.42))
  cowplot::plot_grid(title, page, ncol = 1, rel_heights = c(0.06, 0.94))
}

.dnmb_dbcan_cgc_pages <- function(cgc, clusters_per_page = 12L) {
  if (!nrow(cgc)) return(list())
  palette <- .dnmb_dbcan_plot_palette()
  summary <- cgc |>
    dplyr::group_by(.data$dbcan_cgc_id) |>
    dplyr::summarise(
      dbcan_contig_key = dplyr::first(.data$dbcan_contig_key),
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      families = paste(sort(unique(unlist(lapply(.data$dbcan_all_families[.data$gene_type == "CAZyme"], .dnmb_dbcan_family_tokens)))), collapse = ", "),
      substrate = dplyr::first(stats::na.omit(.data$dbcan_resolved_substrate), default = NA_character_),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$dbcan_contig_key, .data$start)
  summary$local_id <- sub("^.*\\|", "", summary$dbcan_cgc_id)
  summary$facet_label <- paste0(summary$dbcan_contig_key, " - ", summary$local_id)
  summary$facet_label <- ifelse(
    nzchar(summary$families),
    paste0(summary$facet_label, " | ", summary$families),
    summary$facet_label
  )
  summary$facet_label <- ifelse(
    !is.na(summary$substrate) & nzchar(summary$substrate),
    paste0(summary$facet_label, " | ", summary$substrate),
    summary$facet_label
  )
  summary$facet_label <- .dnmb_dbcan_plot_text(summary$facet_label, 105L)

  n_pages <- ceiling(nrow(summary) / max(1L, as.integer(clusters_per_page)))
  page_groups <- if (n_pages <= 1L) {
    list(seq_len(nrow(summary)))
  } else {
    split(seq_len(nrow(summary)), cut(seq_len(nrow(summary)), breaks = n_pages, labels = FALSE))
  }
  lapply(seq_along(page_groups), function(page_index) {
    selected <- summary[page_groups[[page_index]], , drop = FALSE]
    track <- cgc[cgc$dbcan_cgc_id %in% selected$dbcan_cgc_id, , drop = FALSE]
    track$facet_label <- selected$facet_label[match(track$dbcan_cgc_id, selected$dbcan_cgc_id)]
    track$facet_label <- factor(track$facet_label, levels = rev(selected$facet_label))
    offsets <- stats::setNames(selected$start, selected$dbcan_cgc_id)
    track$xmin_rel <- track$start - offsets[track$dbcan_cgc_id]
    track$xmax_rel <- track$end - offsets[track$dbcan_cgc_id]
    track$forward <- !grepl("complement|minus|^-|^-1$", as.character(track$direction), ignore.case = TRUE)
    family_label <- vapply(track$dbcan_all_families, function(value) {
      family <- .dnmb_dbcan_family_tokens(value)
      if (length(family)) paste(family, collapse = "/") else ""
    }, character(1))
    width <- abs(track$xmax_rel - track$xmin_rel)
    track$gene_label <- ifelse(track$gene_type == "CAZyme" & width >= 450, .dnmb_dbcan_plot_text(family_label, 18L), "")

    ggplot2::ggplot(
      track,
      ggplot2::aes(
        xmin = .data$xmin_rel,
        xmax = .data$xmax_rel,
        y = .data$facet_label,
        fill = .data$gene_type,
        forward = .data$forward
      )
    ) +
      gggenes::geom_gene_arrow(
        arrowhead_height = grid::unit(4, "mm"),
        arrowhead_width = grid::unit(2, "mm"),
        arrow_body_height = grid::unit(4, "mm"),
        color = "#334155",
        size = 0.2
      ) +
      ggplot2::geom_text(
        data = track[nzchar(track$gene_label), , drop = FALSE],
        ggplot2::aes(x = (.data$xmin_rel + .data$xmax_rel) / 2, y = .data$facet_label, label = .data$gene_label),
        inherit.aes = FALSE,
        size = 2.3,
        color = "white",
        fontface = "bold",
        check_overlap = TRUE
      ) +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE, name = "Gene type") +
      ggplot2::guides(
        fill = ggplot2::guide_legend(nrow = 1L, byrow = TRUE)
      ) +
      ggplot2::scale_x_continuous(labels = function(x) paste0(round(x / 1000, 1), " kb")) +
      gggenes::theme_genes() +
      ggplot2::labs(
        title = paste0(
          "D  CGC gene context",
          if (length(page_groups) > 1L) {
            paste0("  ", page_index, "/", length(page_groups))
          } else {
            ""
          }
        ),
        subtitle = "Full replicon key is retained; CAZyme labels show all mapped families",
        x = "Coordinate relative to CGC start",
        y = NULL
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.title.position = "plot",
        plot.subtitle = ggplot2::element_text(size = 9, color = "#475569"),
        axis.text.y = ggplot2::element_text(size = 7.2),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.key.height = grid::unit(0.28, "cm"),
        legend.key.width = grid::unit(0.50, "cm"),
        legend.spacing.x = grid::unit(0.10, "cm"),
        plot.margin = ggplot2::margin(8, 12, 8, 2)
      )
  })
}

.dnmb_dbcan_matrix_pages <- function(domains, families_per_page = 30L,
                                     substrates_per_page = 12L,
                                     substrate_mode = c("supported", "family_prior")) {
  substrate_mode <- match.arg(substrate_mode)
  if (!nrow(domains)) return(list())
  data <- as.data.frame(domains, stringsAsFactors = FALSE)
  if (!"substrate_source" %in% names(data)) data$substrate_source <- "legacy"
  if (identical(substrate_mode, "family_prior")) {
    data <- data[
      data$substrate_source == "family_prior" & data$substrate != "unknown",
      ,
      drop = FALSE
    ]
  } else {
    data$substrate[data$substrate_source == "family_prior"] <- "unknown"
  }
  if (!nrow(data)) return(list())

  links <- data |>
    tidyr::separate_rows("substrate", sep = ";\\s*") |>
    dplyr::mutate(substrate = ifelse(is.na(.data$substrate) | !nzchar(.data$substrate), "unknown", .data$substrate)) |>
    dplyr::count(.data$family, .data$substrate, name = "genes")
  if (identical(substrate_mode, "family_prior")) {
    links <- links[links$substrate != "unknown", , drop = FALSE]
  }
  if (!nrow(links)) return(list())

  balanced_groups <- function(values, max_per_page) {
    n_pages <- ceiling(length(values) / max(1L, as.integer(max_per_page)))
    if (n_pages <= 1L) return(list(values))
    unname(split(values, cut(seq_along(values), breaks = n_pages, labels = FALSE)))
  }
  substrate_order <- links |>
    dplyr::group_by(.data$substrate) |>
    dplyr::summarise(total = sum(.data$genes), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$total), .data$substrate)
  substrate_groups <- balanced_groups(substrate_order$substrate, substrates_per_page)
  specs <- list()
  for (substrates in substrate_groups) {
    subset <- links[links$substrate %in% substrates, , drop = FALSE]
    family_order <- subset |>
      dplyr::group_by(.data$family) |>
      dplyr::summarise(total = sum(.data$genes), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(.data$total), .data$family)
    family_groups <- balanced_groups(family_order$family, families_per_page)
    for (families in family_groups) {
      specs[[length(specs) + 1L]] <- list(
        data = subset[subset$family %in% families, , drop = FALSE],
        families = families,
        substrates = substrates
      )
    }
  }

  title_text <- if (identical(substrate_mode, "family_prior")) {
    "F  Family-level substrate prior (hypothesis only)"
  } else {
    "E  Substrate-family support"
  }
  subtitle_text <- if (identical(substrate_mode, "family_prior")) {
    "Broad family associations only; these are not gene-specific substrate evidence"
  } else {
    "PUL/dbCAN-sub/overview support; unknown has no gene- or CGC-level substrate assignment"
  }
  colors <- if (identical(substrate_mode, "family_prior")) {
    c(low = "#F8D9A0", high = "#D97706")
  } else {
    c(low = "#BFE3DA", high = "#087F8C")
  }

  lapply(seq_along(specs), function(page_index) {
    spec <- specs[[page_index]]
    plot_data <- spec$data
    plot_data$family <- factor(plot_data$family, levels = rev(spec$families))
    plot_data$substrate <- factor(plot_data$substrate, levels = spec$substrates)
    ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$substrate, y = .data$family)) +
      ggplot2::geom_point(ggplot2::aes(size = .data$genes, fill = .data$genes), shape = 21, color = "#334155", stroke = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = .data$genes), size = 2.5, color = "#0F172A") +
      ggplot2::scale_size_continuous(range = c(3, 11), guide = "none") +
      ggplot2::scale_fill_gradient(low = colors[["low"]], high = colors[["high"]], guide = "none") +
      ggplot2::labs(
        title = paste0(
          title_text,
          if (length(specs) > 1L) {
            paste0("  ", page_index, "/", length(specs))
          } else {
            ""
          }
        ),
        subtitle = subtitle_text,
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.title.position = "plot",
        plot.subtitle = ggplot2::element_text(size = 9, color = "#475569"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = ggplot2::element_text(size = 8),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(8, 12, 12, 2)
      )
  })
}

.dnmb_dbcan_one_page_layout <- function(inputs) {
  hit_groups <- .dnmb_dbcan_hit_contig_groups(
    inputs,
    contigs_per_page = .Machine$integer.max
  )
  contig_keys <- unique(unlist(hit_groups, use.names = FALSE))
  n_contigs <- length(contig_keys)
  contig_ncol <- if (n_contigs <= 4L) {
    1L
  } else if (n_contigs <= 12L) {
    2L
  } else {
    3L
  }
  summary_plot <- .dnmb_dbcan_summary_page(
    inputs,
    contig_keys = contig_keys,
    include_details = TRUE,
    contig_ncol = contig_ncol
  )

  n_cgc <- length(unique(as.character(inputs$cgc$dbcan_cgc_id)))
  cgc_plot <- if (n_cgc > 0L) {
    .dnmb_dbcan_cgc_pages(
      inputs$cgc,
      clusters_per_page = max(1L, n_cgc)
    )[[1L]]
  } else {
    NULL
  }

  matrix_plots <- c(
    .dnmb_dbcan_matrix_pages(
      inputs$domains,
      families_per_page = .Machine$integer.max,
      substrates_per_page = .Machine$integer.max,
      substrate_mode = "supported"
    ),
    .dnmb_dbcan_matrix_pages(
      inputs$domains,
      families_per_page = .Machine$integer.max,
      substrates_per_page = .Machine$integer.max,
      substrate_mode = "family_prior"
    )
  )
  matrix_plots <- Filter(Negate(is.null), matrix_plots)
  matrix_rows <- if (length(matrix_plots)) {
    vapply(matrix_plots, function(plot) {
      if (!is.data.frame(plot$data) || !"family" %in% names(plot$data)) return(1L)
      max(1L, length(unique(as.character(plot$data$family))))
    }, integer(1))
  } else {
    integer(0)
  }
  matrix_plot <- if (length(matrix_plots) == 1L) {
    matrix_plots[[1L]]
  } else if (length(matrix_plots) > 1L) {
    cowplot::plot_grid(
      plotlist = matrix_plots,
      ncol = 1L,
      rel_heights = pmax(1, matrix_rows)
    )
  } else {
    NULL
  }

  if (!is.null(cgc_plot) && !is.null(matrix_plot)) {
    detail_plot <- cowplot::plot_grid(
      cgc_plot, matrix_plot,
      ncol = 2L,
      rel_widths = c(1.55, 1)
    )
    page_width <- 28
  } else if (!is.null(cgc_plot)) {
    detail_plot <- cgc_plot
    page_width <- 20
  } else if (!is.null(matrix_plot)) {
    detail_plot <- matrix_plot
    page_width <- 18
  } else {
    detail_plot <- NULL
    page_width <- 14
  }

  contig_rows <- max(1L, ceiling(max(1L, n_contigs) / contig_ncol))
  summary_height <- max(7, 5.8 + 0.55 * max(0L, contig_rows - 2L))
  detail_rows <- max(c(n_cgc, sum(matrix_rows), 0L))
  detail_height <- if (is.null(detail_plot)) 0 else max(10, 2.4 + 0.21 * detail_rows)

  if (is.null(detail_plot)) {
    page_plot <- summary_plot
    page_height <- max(9, summary_height)
  } else {
    page_plot <- cowplot::plot_grid(
      summary_plot, detail_plot,
      ncol = 1L,
      rel_heights = c(summary_height, detail_height)
    )
    page_height <- summary_height + detail_height
  }

  list(
    plot = page_plot,
    width = page_width,
    height = page_height,
    n_contigs = n_contigs,
    n_cgc = n_cgc,
    n_matrix_rows = sum(matrix_rows)
  )
}

.dnmb_plot_dbcan_module <- function(genbank_table, output_dir) {
  inputs <- .dnmb_dbcan_read_plot_inputs(genbank_table, output_dir)
  if (!nrow(inputs$cazy) && !nrow(inputs$cgc)) {
    stop("dbCAN plot: no mapped CAZyme or CGC data found.", call. = FALSE)
  }

  layout <- .dnmb_dbcan_one_page_layout(inputs)
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "dbcan_cgc_overview.pdf")
  .dnmb_plot_pdf_device(
    pdf_path,
    width = layout$width,
    height = layout$height,
    onefile = FALSE
  )
  device_open <- TRUE
  on.exit({
    if (device_open) grDevices::dev.off()
  }, add = TRUE)
  print(layout$plot)
  grDevices::dev.off()
  device_open <- FALSE

  list(
    pdf = pdf_path,
    pages = 1L,
    width = layout$width,
    height = layout$height
  )
}

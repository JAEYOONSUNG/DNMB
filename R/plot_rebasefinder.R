.dnmb_rebasefinder_cleanup_native_pdfs <- function(output_dir) {
  output_dir <- base::normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )
  run_root <- if (base::basename(output_dir) == "dnmb_module_rebasefinder") {
    base::dirname(output_dir)
  } else {
    output_dir
  }
  legacy_pdfs <- base::unique(c(
    base::file.path(run_root, "RM_system_dotplot.pdf"),
    base::file.path(run_root, "dnmb_module_rebasefinder", "RM_system_dotplot.pdf")
  ))
  existing <- legacy_pdfs[base::file.exists(legacy_pdfs)]
  if (base::length(existing)) {
    base::unlink(existing, force = TRUE)
  }
  base::invisible(existing)
}

.dnmb_plot_rebasefinder_overview <- function(genbank_table, output_dir, cache_root = NULL) {
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
  tbl <- .dnmb_rebasefinder_ensure_partial_cols(tbl)
  if ("REBASEfinder_typing_eligible" %in% names(tbl)) {
    tbl$REBASEfinder_typing_eligible <- as.logical(tbl$REBASEfinder_typing_eligible)
  }
  # Curated roles are authoritative. Hit-name parsing is a legacy fallback only.
  if ("REBASEfinder_hit_label" %in% names(tbl) && "REBASEfinder_enzyme_role" %in% names(tbl)) {
    hit <- as.character(tbl$REBASEfinder_hit_label)
    inferred <- .dnmb_rebasefinder_role_from_hit(hit)
    fix <- !is.na(inferred) & (is.na(tbl$REBASEfinder_enzyme_role) |
           !nzchar(as.character(tbl$REBASEfinder_enzyme_role)))
    tbl$REBASEfinder_enzyme_role[fix] <- inferred[fix]
  }

  rm_palette <- .dnmb_rebasefinder_palette(tbl$REBASEfinder_family_id)
  role_palette <- .dnmb_rebasefinder_role_palette(tbl)

  # Motif verification
  motif_verified <- .dnmb_rebasefinder_motif_verified(tbl)
  motif_files <- .dnmb_rebasefinder_write_motif_hits(output_dir, tbl)
  motif_geometry <- if (!base::is.null(motif_files) &&
                        base::is.data.frame(motif_files$qualified)) {
    tryCatch(
      .dnmb_rebasefinder_write_motif_geometry(
        output_dir = output_dir,
        tbl = tbl,
        motif_hits = motif_files$qualified
      ),
      error = function(e) list(pairs = data.frame(), summary = data.frame())
    )
  } else {
    list(pairs = data.frame(), summary = data.frame())
  }
  motif_contacts <- tryCatch(
    .dnmb_rebasefinder_read_motif_contacts(output_dir),
    error = function(e) list(pairs = data.frame(), summary = data.frame())
  )

  # Confidence filter
  has_eligible <- "REBASEfinder_typing_eligible" %in% names(tbl)
  has_tier <- "REBASEfinder_curation_tier" %in% names(tbl)
  has_ev <- "REBASEfinder_evidence_mode" %in% names(tbl)
  is_hc <- if (has_tier) {
    !is.na(tbl$REBASEfinder_curation_tier) & tbl$REBASEfinder_curation_tier == "high"
  } else if (has_eligible) {
    !is.na(tbl$REBASEfinder_typing_eligible) & tbl$REBASEfinder_typing_eligible == TRUE
  } else if (has_ev) {
    !is.na(tbl$REBASEfinder_evidence_mode) & tbl$REBASEfinder_evidence_mode == "high_confidence"
  } else {
    rep(TRUE, nrow(tbl))
  }

  inventory_scope <- .dnmb_rebasefinder_inventory_scope(tbl)
  inventory_tbl <- tbl[inventory_scope$system_mask, , drop = FALSE]
  inventory_hc <- is_hc[inventory_scope$system_mask]
  inventory_verified <- motif_verified[inventory_scope$system_mask]
  n_hc <- sum(inventory_hc, na.rm = TRUE)
  n_verified <- sum(inventory_hc & inventory_verified, na.rm = TRUE)

  p_inventory <- .dnmb_plot_rebasefinder_inventory(
    inventory_tbl, rm_palette,
    n_total = nrow(inventory_tbl), n_hc = n_hc, n_verified = n_verified,
    n_orphan_mtase = inventory_scope$n_orphan_mtase,
    n_unassociated = inventory_scope$n_unassociated
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
    .dnmb_plot_rebasefinder_blast_quality(tbl, role_palette, display_info, rm_palette, cache_root = cache_root),
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
	    .dnmb_plot_rebasefinder_motif_verification(
	      tbl,
	      display_info,
	      geometry_pairs = motif_geometry$pairs,
	      geometry_summary = motif_geometry$summary,
	      contact_pairs = motif_contacts$pairs,
	      contact_summary = motif_contacts$summary
	    ),
	    error = function(e) empty_plot
	  )

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, "REBASE_overview.pdf")

  n_hits <- nrow(display_info)
  cd_inch <- max(2.8, 0.42 * n_hits + 0.8)

  n_rm_types <- length(unique(inventory_tbl$REBASEfinder_family_id))
  a_inch <- max(0.7, 0.45 * n_rm_types + 0.5)

  # Build bottom row (C): BLAST, protein domain map, square motif heatmap.
  bottom_row <- tryCatch(
    cowplot::plot_grid(
      p_blast, p_domain, p_motif,
      ncol = 3, rel_widths = c(1.55, 0.50, 1.70),
      align = "h", axis = "tb"
    ),
    error = function(e) NULL
  )
  if (is.null(bottom_row)) {
    bottom_row <- tryCatch(
      cowplot::plot_grid(
        p_blast, p_motif,
        ncol = 2, rel_widths = c(1.62, 1.0),
        align = "h", axis = "tb"
      ),
      error = function(e) NULL
    )
  }
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
  .dnmb_module_plot_save(composite, pdf_path, width = 18, height = min(22, max(10, total_height)))
  if (base::file.exists(pdf_path) && base::file.info(pdf_path)$size > 500) {
    .dnmb_rebasefinder_cleanup_native_pdfs(output_dir)
  }
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

.dnmb_rebasefinder_pretty_hit_label <- function(x) {
  x <- as.character(x)
  out <- x
  out <- ifelse(grepl("^typeI_context:R:", out), "Type I R (operon)", out)
  out <- ifelse(grepl("^typeI_context:M:", out), "Type I M (operon)", out)
  out <- ifelse(grepl("^typeI_context:S:", out), "Type I S (operon)", out)
  out <- ifelse(grepl("^typeIII_context:R:", out), "Type III R (operon)", out)
  out <- ifelse(grepl("^typeIII_context:M:", out), "Type III M (operon)", out)
  out <- ifelse(grepl("^typeII_context:R:", out), "Type II R (operon)", out)
  out <- ifelse(grepl("^typeII_context:M:", out), "Type II M (operon)", out)
  out <- ifelse(grepl("^typeIV_context:", out), "Type IV MDRE", out)
  out
}

.dnmb_rebasefinder_ensure_partial_cols <- function(tbl) {
  tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!base::nrow(tbl)) return(tbl)
  if (all(c("REBASEfinder_aa_len", "REBASEfinder_expected_min_aa",
            "REBASEfinder_partial_status", "REBASEfinder_partial_reason") %in% base::names(tbl))) {
    return(tbl)
  }
  partial <- .dnmb_rebasefinder_sequence_partial_table(
    tbl,
    family_col = "REBASEfinder_family_id",
    role_col = "REBASEfinder_enzyme_role"
  )
  if (!"REBASEfinder_aa_len" %in% base::names(tbl)) tbl$REBASEfinder_aa_len <- partial$aa_len
  if (!"REBASEfinder_expected_min_aa" %in% base::names(tbl)) tbl$REBASEfinder_expected_min_aa <- partial$expected_min_aa
  if (!"REBASEfinder_partial_status" %in% base::names(tbl)) tbl$REBASEfinder_partial_status <- partial$partial_status
  if (!"REBASEfinder_partial_reason" %in% base::names(tbl)) tbl$REBASEfinder_partial_reason <- partial$partial_reason
  tbl
}

.dnmb_rebasefinder_structure_motif_status <- function(row) {
  structure_status <- if ("REBASEfinder_structure_status" %in% names(row)) {
    as.character(row$REBASEfinder_structure_status[[1]])
  } else {
    NA_character_
  }
  if (is.na(structure_status) || !nzchar(structure_status)) return("structure_not_available")
  if (!grepl("^structure_supported", structure_status)) return(structure_status)
  role <- if ("REBASEfinder_enzyme_role" %in% names(row)) as.character(row$REBASEfinder_enzyme_role[[1]]) else NA_character_
  structure_role <- if ("REBASEfinder_structure_role" %in% names(row)) as.character(row$REBASEfinder_structure_role[[1]]) else NA_character_
  chain_role <- if ("REBASEfinder_structure_chain_role" %in% names(row)) as.character(row$REBASEfinder_structure_chain_role[[1]]) else NA_character_
  ref_role <- stats::na.omit(c(structure_role, chain_role))
  if (!length(ref_role) || is.na(role) || !nzchar(role)) return("structure_supported")
  if (any(grepl(role, ref_role, fixed = TRUE))) return("structure_supported_role_consistent")
  "structure_supported_role_check"
}

.dnmb_rebasefinder_motif_hits_table <- function(tbl, role_relevant_only = TRUE) {
  tbl <- .dnmb_rebasefinder_ensure_partial_cols(tbl)
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  if (is.null(detailed)) return(data.frame())
  motif_defs <- .dnmb_rebasefinder_motif_definitions()
  motif_names <- names(motif_defs)
  annotations <- .dnmb_rebasefinder_gene_annotation_text(tbl)

  rows <- list()
  n <- 0L
  for (i in seq_len(nrow(tbl))) {
    locus <- as.character(tbl$locus_tag[[i]])
    row <- tbl[i, , drop = FALSE]
    seq_len <- if ("REBASEfinder_aa_len" %in% names(tbl)) tbl$REBASEfinder_aa_len[[i]] else NA_integer_
    role <- if ("REBASEfinder_enzyme_role" %in% names(tbl)) as.character(tbl$REBASEfinder_enzyme_role[[i]]) else NA_character_
    family <- if ("REBASEfinder_family_id" %in% names(tbl)) as.character(tbl$REBASEfinder_family_id[[i]]) else NA_character_
    hit_label <- if ("REBASEfinder_hit_label" %in% names(tbl)) as.character(tbl$REBASEfinder_hit_label[[i]]) else NA_character_
    structure_check <- .dnmb_rebasefinder_structure_motif_status(row)
    assessment <- .dnmb_rebasefinder_motif_assessment(
      detailed[[i]], role = role, family = family, annotation = annotations[[i]]
    )

    for (mn in motif_names) {
      info <- detailed[[i]][[which(motif_names == mn)]]
      def <- motif_defs[[mn]]
      expected_role <- .dnmb_rebasefinder_motif_role_label(def)
      role_relevant <- if (is.na(role) || !nzchar(role)) {
        NA
      } else {
        .dnmb_rebasefinder_motif_role_match(def, role, family)
      }
      motif_evidence_level <- if (isTRUE(assessment$supported[[mn]])) {
        if (identical(def$evidence_kind, "domain_architecture")) "architecture_supported" else "supported"
      } else if (isTRUE(assessment$hint[[mn]])) {
        if (identical(def$evidence_kind, "domain_architecture")) "architecture_candidate" else "sequence_hint"
      } else if (isTRUE(role_relevant)) {
        "raw_unqualified"
      } else {
        "raw_non_role"
      }
      hit_rows <- if (!is.null(info) && !is.null(info$hits) && length(info$hits)) {
        info$hits
      } else if (motif_evidence_level %in% c("architecture_supported", "architecture_candidate")) {
        list(list(
          match = assessment$derived_label[[mn]], pos = NA_integer_, end = NA_integer_,
          in_range = TRUE
        ))
      } else {
        list()
      }
      if (!length(hit_rows)) next
      for (h in hit_rows) {
        evidence_level <- motif_evidence_level
        if (evidence_level == "supported" && !is.null(def$pos_range) && !isTRUE(h$in_range)) {
          evidence_level <- "raw_unqualified"
        }
        n <- n + 1L
        rows[[n]] <- data.frame(
          locus_tag = locus,
          family_id = family,
          enzyme_role = role,
          hit_label = hit_label,
          motif = mn,
          motif_description = def$full,
          regex = def$pattern,
          evidence_kind = def$evidence_kind,
          evidence_level = evidence_level,
          coherent_architecture = evidence_level %in% c("supported", "architecture_supported"),
          gene_functional_verified = assessment$verified,
          mtase_architecture = assessment$mtase_architecture,
          typeiv_architecture = assessment$typeiv_architecture,
          expected_role = expected_role,
          role_relevant = role_relevant,
          match = h$match,
          start_aa = h$pos,
          end_aa = h$end,
          aa_len = seq_len,
          relative_pos = if (!is.na(seq_len) && seq_len > 0) h$pos / seq_len else NA_real_,
          in_expected_position = h$in_range,
          partial_status = if ("REBASEfinder_partial_status" %in% names(tbl)) tbl$REBASEfinder_partial_status[[i]] else NA_character_,
          partial_reason = if ("REBASEfinder_partial_reason" %in% names(tbl)) tbl$REBASEfinder_partial_reason[[i]] else NA_character_,
          structure_status = if ("REBASEfinder_structure_status" %in% names(tbl)) as.character(tbl$REBASEfinder_structure_status[[i]]) else NA_character_,
          structure_reference_id = if ("REBASEfinder_structure_reference_id" %in% names(tbl)) as.character(tbl$REBASEfinder_structure_reference_id[[i]]) else NA_character_,
          structure_role = if ("REBASEfinder_structure_role" %in% names(tbl)) as.character(tbl$REBASEfinder_structure_role[[i]]) else NA_character_,
          structure_chain_role = if ("REBASEfinder_structure_chain_role" %in% names(tbl)) as.character(tbl$REBASEfinder_structure_chain_role[[i]]) else NA_character_,
          structural_verification = structure_check,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  if (isTRUE(role_relevant_only) && "role_relevant" %in% names(out)) {
    qualified <- out$evidence_level %in% c(
      "supported", "architecture_supported", "sequence_hint", "architecture_candidate"
    )
    out <- out[(is.na(out$role_relevant) | out$role_relevant) & qualified, , drop = FALSE]
  }
  if (!nrow(out)) return(out)
  out[order(out$locus_tag, out$start_aa, out$motif), , drop = FALSE]
}

.dnmb_rebasefinder_write_motif_hits <- function(output_dir, tbl) {
  raw <- .dnmb_rebasefinder_motif_hits_table(tbl, role_relevant_only = FALSE)
  if (!nrow(raw)) return(invisible(NULL))
  out <- .dnmb_rebasefinder_motif_hits_table(tbl, role_relevant_only = TRUE)
  module_dir <- file.path(output_dir, "dnmb_module_rebasefinder")
  dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
  # Remove stray copies left in the run root by older versions of this writer.
  for (stray in c("DNMB_REBASEfinder_motif_hits.tsv",
                  "DNMB_REBASEfinder_motif_hits.xlsx",
                  "DNMB_REBASEfinder_motif_hits_raw.tsv",
                  "DNMB_REBASEfinder_motif_hits_raw.xlsx")) {
    stray_path <- file.path(output_dir, stray)
    if (file.exists(stray_path)) unlink(stray_path, force = TRUE)
  }
  tsv <- file.path(module_dir, "DNMB_REBASEfinder_motif_hits.tsv")
  utils::write.table(out, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  xlsx <- file.path(module_dir, "DNMB_REBASEfinder_motif_hits.xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    tryCatch(openxlsx::write.xlsx(out, xlsx, overwrite = TRUE), error = function(e) NULL)
  }
  raw_tsv <- file.path(module_dir, "DNMB_REBASEfinder_motif_hits_raw.tsv")
  utils::write.table(raw, raw_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  raw_xlsx <- file.path(module_dir, "DNMB_REBASEfinder_motif_hits_raw.xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    tryCatch(openxlsx::write.xlsx(raw, raw_xlsx, overwrite = TRUE), error = function(e) NULL)
  }
  invisible(list(
    tsv = tsv,
    xlsx = if (file.exists(xlsx)) xlsx else NA_character_,
    raw_tsv = raw_tsv,
    raw_xlsx = if (file.exists(raw_xlsx)) raw_xlsx else NA_character_,
    qualified = out,
    raw = raw
  ))
}


# ====================================================================
# Shared display label builder — operon grouping + C/D y-axis sync
# ====================================================================
.dnmb_rebasefinder_partner_loci <- function(x, known_loci) {
  x <- base::as.character(x)
  x <- x[!base::is.na(x) & base::nzchar(x)]
  if (!base::length(x)) return(character())
  parts <- base::unlist(base::strsplit(x, "\\s*\\|\\s*", perl = TRUE), use.names = FALSE)
  parts <- base::trimws(parts)
  parts <- parts[base::nzchar(parts)]
  loci <- base::sub(":.*$", "", parts)
  loci <- loci[loci %in% known_loci]
  base::unique(loci)
}

.dnmb_rebasefinder_operon_groups_for_plot <- function(tbl) {
  tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  n <- base::nrow(tbl)
  if (!n || !"locus_tag" %in% base::names(tbl)) {
    return(base::data.frame(
      locus_tag = character(), operon_group = character(),
      group_contig_rank = integer(), group_start = numeric(),
      member_start = numeric(), member_end = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  loci <- base::as.character(tbl$locus_tag)
  parent <- base::seq_len(n)
  find <- function(i) {
    while (!base::identical(parent[[i]], i)) {
      parent[[i]] <<- parent[[parent[[i]]]]
      i <- parent[[i]]
    }
    i
  }
  union <- function(a, b) {
    ia <- base::match(a, loci)
    ib <- base::match(b, loci)
    if (base::is.na(ia) || base::is.na(ib)) return(invisible(NULL))
    ra <- find(ia); rb <- find(ib)
    if (!base::identical(ra, rb)) parent[[rb]] <<- ra
    invisible(NULL)
  }

  if ("REBASEfinder_operon_id" %in% base::names(tbl)) {
    op <- base::as.character(tbl$REBASEfinder_operon_id)
    op[base::is.na(op) | !base::nzchar(op)] <- NA_character_
    for (oid in base::unique(stats::na.omit(op))) {
      members <- loci[op == oid]
      if (base::length(members) > 1L) {
        base::invisible(base::lapply(members[-1L], function(m) union(members[[1]], m)))
      }
    }
  }

  partner_cols <- base::intersect(
    c("REBASEfinder_typei_context_partners", "REBASEfinder_typeiii_context_partners"),
    base::names(tbl)
  )
  if (base::length(partner_cols)) {
    for (col in partner_cols) {
      for (i in base::seq_len(n)) {
        partners <- .dnmb_rebasefinder_partner_loci(tbl[[col]][[i]], loci)
        if (base::length(partners)) {
          base::invisible(base::lapply(partners, function(p) union(loci[[i]], p)))
        }
      }
    }
  }

  root <- base::vapply(base::seq_len(n), find, integer(1))
  group_members <- base::split(loci, root)
  group_label <- base::vapply(root, function(r) {
    base::paste(base::sort(group_members[[base::as.character(r)]]), collapse = "+")
  }, character(1))

  contig <- if ("contig" %in% base::names(tbl)) base::as.character(tbl$contig) else base::rep("contig", n)
  contig[base::is.na(contig) | !base::nzchar(contig)] <- "contig"
  contig_rank <- base::match(contig, base::unique(contig))
  start <- if ("start" %in% base::names(tbl)) suppressWarnings(base::as.numeric(tbl$start)) else base::seq_len(n)
  end <- if ("end" %in% base::names(tbl)) suppressWarnings(base::as.numeric(tbl$end)) else start
  start[base::is.na(start)] <- base::seq_len(n)[base::is.na(start)]
  end[base::is.na(end)] <- start[base::is.na(end)]

  group_start <- stats::ave(start, group_label, FUN = function(z) min(z, na.rm = TRUE))
  group_contig_rank <- stats::ave(contig_rank, group_label, FUN = function(z) min(z, na.rm = TRUE))
  base::data.frame(
    locus_tag = loci,
    operon_group = group_label,
    group_contig_rank = group_contig_rank,
    group_start = group_start,
    member_start = start,
    member_end = end,
    stringsAsFactors = FALSE
  )
}

.dnmb_rebasefinder_display_labels <- function(tbl) {
  has_hit      <- "REBASEfinder_hit_label" %in% names(tbl)
  has_rec      <- "REBASEfinder_rec_seq" %in% names(tbl)
  has_ref_rec  <- "REBASEfinder_reference_rec_seq" %in% names(tbl)
  has_eligible <- "REBASEfinder_typing_eligible" %in% names(tbl)
  has_tier <- "REBASEfinder_curation_tier" %in% names(tbl)
  has_identity <- "REBASEfinder_blast_identity" %in% names(tbl)
  has_blastlen <- "REBASEfinder_blast_length" %in% names(tbl)
  has_translation <- "translation" %in% names(tbl)
  has_partial <- "REBASEfinder_partial_status" %in% names(tbl)
  has_aa_len <- "REBASEfinder_aa_len" %in% names(tbl)
  group_info <- .dnmb_rebasefinder_operon_groups_for_plot(tbl)

  hit_label <- if (has_hit) .dnmb_rebasefinder_pretty_hit_label(tbl$REBASEfinder_hit_label) else tbl$locus_tag
  rec_seq   <- if (has_rec) as.character(tbl$REBASEfinder_rec_seq) else NA_character_
  ref_rec_seq <- if (has_ref_rec) as.character(tbl$REBASEfinder_reference_rec_seq) else NA_character_
  valid_rec <- !is.na(rec_seq) & nzchar(rec_seq) & rec_seq != "?"
  valid_ref_rec <- !is.na(ref_rec_seq) & nzchar(ref_rec_seq) & ref_rec_seq != "?"
  identity  <- if (has_identity) tbl$REBASEfinder_blast_identity else rep(NA_real_, nrow(tbl))
  operon    <- group_info$operon_group[match(tbl$locus_tag, group_info$locus_tag)]
  operon[is.na(operon) | !nzchar(operon)] <- as.character(tbl$locus_tag[is.na(operon) | !nzchar(operon)])
  confidence <- if (has_tier) {
    ifelse(tbl$REBASEfinder_curation_tier == "high", "High",
      ifelse(tbl$REBASEfinder_curation_tier == "medium", "Medium", "Low"))
  } else if (has_eligible) {
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

  # Keep recognition sequence and locus on a second line so long labels do
  # not consume the C panel width.
  base_label <- hit_label
  id_str <- ifelse(!is.na(identity), paste0(round(identity * 100, 1), "%"), "")
  cov_str <- ifelse(!is.na(coverage), paste0(round(coverage * 100, 0), "%"), "")
  qual_tag <- ifelse(
    nzchar(id_str) & nzchar(cov_str),
    paste0(" [", id_str, "/", cov_str, "]"),
    ifelse(nzchar(id_str), paste0(" [", id_str, "]"), "")
  )
  partial_tag <- if (has_partial) {
    is_partial <- !is.na(tbl$REBASEfinder_partial_status) &
      tbl$REBASEfinder_partial_status == "partial_or_short"
    len <- if (has_aa_len) suppressWarnings(as.integer(tbl$REBASEfinder_aa_len)) else NA_integer_
    ifelse(is_partial,
           ifelse(!is.na(len), paste0(" [partial/short ", len, " aa]"), " [partial/short]"),
           "")
  } else {
    rep("", nrow(tbl))
  }
  second_line <- ifelse(
    valid_rec,
    paste0("rec: ", rec_seq, " | ", tbl$locus_tag),
    ifelse(
      valid_ref_rec,
      paste0("ref rec: ", ref_rec_seq, " | ", tbl$locus_tag),
      as.character(tbl$locus_tag)
    )
  )
  display <- paste0(base_label, qual_tag, partial_tag, "\n", second_line)

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
  group_info <- group_info[match(out$locus_tag, group_info$locus_tag), , drop = FALSE]
  role_order <- base::match(out$role, c("R", "M", "S", "Other"))
  role_order[base::is.na(role_order)] <- 99L
  # Sort by genomic operon block first, then by the member's genomic order.
  # This keeps M/R/S partners adjacent across all C subpanels.
  out <- out[base::order(group_info$group_contig_rank, group_info$group_start,
                         group_info$member_start, group_info$member_end,
                         role_order, out$locus_tag), , drop = FALSE]
	  out$display <- factor(out$display, levels = rev(unique(out$display)))
	  out
	}

.dnmb_rebasefinder_display_levels <- function(display_info) {
  lvls <- levels(display_info$display)
  if (is.null(lvls)) lvls <- unique(as.character(display_info$display))
  lvls[!is.na(lvls) & nzchar(lvls)]
}


# ====================================================================
# Panel A
# ====================================================================
.dnmb_rebasefinder_inventory_scope <- function(tbl) {
  tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  n <- base::nrow(tbl)
  if (!n || !"REBASEfinder_rm_association_class" %in% base::names(tbl)) {
    return(list(
      system_mask = base::rep(TRUE, n),
      n_orphan_mtase = 0L,
      n_unassociated = 0L
    ))
  }

  association <- base::as.character(tbl$REBASEfinder_rm_association_class)
  supported <- association %in% c(
    "classic_RM_supported",
    "fused_RM_supported",
    "modification_dependent_restriction_supported"
  )
  orphan_mtase <- association == "DNA_MTase_RM_association_unproven"
  orphan_mtase[base::is.na(orphan_mtase)] <- FALSE
  supported[base::is.na(supported)] <- FALSE
  list(
    system_mask = supported,
    n_orphan_mtase = base::sum(orphan_mtase),
    n_unassociated = base::sum(!supported)
  )
}

.dnmb_plot_rebasefinder_inventory <- function(tbl, palette,
                                              n_total = nrow(tbl),
                                              n_hc = nrow(tbl),
                                              n_verified = nrow(tbl),
                                              n_orphan_mtase = 0L,
                                              n_unassociated = 0L) {
  if (!nrow(tbl)) {
    detail <- if (n_orphan_mtase > 0L) {
      paste0(n_orphan_mtase, " orphan DNA MTase locus/loci retained in detailed panels")
    } else if (n_unassociated > 0L) {
      paste0(n_unassociated, " unassociated candidate locus/loci retained in detailed panels")
    } else {
      "No curated R-M loci"
    }
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = paste0("No association-supported R-M loci\n", detail),
                               size = 4, color = "grey50"))
  }

	  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
	  has_rec  <- "REBASEfinder_rec_seq" %in% names(tbl)
	  has_ref_rec <- "REBASEfinder_reference_rec_seq" %in% names(tbl)
	  has_id   <- "REBASEfinder_blast_identity" %in% names(tbl)
	  has_eligible <- "REBASEfinder_typing_eligible" %in% names(tbl)
	  has_tier <- "REBASEfinder_curation_tier" %in% names(tbl)
	  has_struct <- "REBASEfinder_structure_status" %in% names(tbl)
	  has_typeiii_ctx <- "REBASEfinder_typeiii_context_status" %in% names(tbl)
	  identity_vec <- if (has_id) suppressWarnings(as.numeric(tbl$REBASEfinder_blast_identity)) else rep(NA_real_, nrow(tbl))
	  eligible_vec <- if (has_tier) {
	    tbl$REBASEfinder_curation_tier %in% c("high", "medium")
	  } else if (has_eligible) {
	    as.logical(tbl$REBASEfinder_typing_eligible)
	  } else {
	    rep(TRUE, nrow(tbl))
	  }
	  tbl$.inventory_meaningful <- eligible_vec %in% TRUE & (is.na(identity_vec) | identity_vec > 0.40)
	  rec <- if (has_rec) as.character(tbl$REBASEfinder_rec_seq) else rep(NA_character_, nrow(tbl))
	  ref_rec <- if (has_ref_rec) as.character(tbl$REBASEfinder_reference_rec_seq) else rep(NA_character_, nrow(tbl))
	  valid_rec <- !is.na(rec) & nzchar(rec) & rec != "?"
	  valid_ref_rec <- !is.na(ref_rec) & nzchar(ref_rec) & ref_rec != "?"
	  tbl$.inventory_rec_label <- ifelse(
	    valid_rec,
	    paste0("rec:", rec),
	    ifelse(valid_ref_rec, paste0("ref:", ref_rec), NA_character_)
	  )

	  inv <- tbl |>
	    dplyr::group_by(.data$REBASEfinder_family_id) |>
	    dplyr::summarise(
	      n_genes = dplyr::n(),
	      n_meaningful = sum(.data$.inventory_meaningful, na.rm = TRUE),
	      n_weak = .data$n_genes - .data$n_meaningful,
	      roles = if (has_role) paste(sort(unique(stats::na.omit(.data$REBASEfinder_enzyme_role))), collapse = "/") else "",
      n_structure = if (has_struct) sum(grepl("^structure_supported", .data$REBASEfinder_structure_status), na.rm = TRUE) else 0L,
      n_typeiii_complete = if (has_typeiii_ctx) sum(.data$REBASEfinder_typeiii_context_status == "complete_mod_res", na.rm = TRUE) else 0L,
	      rec_seqs = {
	        seqs <- unique(stats::na.omit(.data$.inventory_rec_label))
	        if (length(seqs)) paste(seqs[seq_len(min(3, length(seqs)))], collapse = ", ") else ""
	      },
      mean_identity = if (has_id) mean(.data$REBASEfinder_blast_identity, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  rm_order <- c("Type I", "Type II", "Type III", "Type IV")
  present_order <- c(
    rm_order[rm_order %in% inv$REBASEfinder_family_id],
    setdiff(as.character(inv$REBASEfinder_family_id), rm_order)
  )
	  inv <- inv[match(present_order, inv$REBASEfinder_family_id), , drop = FALSE]
	  inv$REBASEfinder_family_id <- factor(inv$REBASEfinder_family_id, levels = rev(present_order))
	  inv_bar <- rbind(
	    data.frame(REBASEfinder_family_id = inv$REBASEfinder_family_id,
	               segment = "meaningful", n = inv$n_meaningful, stringsAsFactors = FALSE),
	    data.frame(REBASEfinder_family_id = inv$REBASEfinder_family_id,
	               segment = "weak", n = inv$n_weak, stringsAsFactors = FALSE)
	  )
	  inv_bar <- inv_bar[inv_bar$n > 0, , drop = FALSE]
	  inv_bar$REBASEfinder_family_id <- factor(inv_bar$REBASEfinder_family_id, levels = levels(inv$REBASEfinder_family_id))
	  inv_bar$segment <- factor(inv_bar$segment, levels = c("weak", "meaningful"))

	  inv$annot <- vapply(seq_len(nrow(inv)), function(i) {
	    ng <- inv$n_genes[i]
	    line1 <- paste(c(
	      paste0(ng, ifelse(ng == 1, " gene", " genes")),
	      paste0("meaningful: ", inv$n_meaningful[i]),
	      if (nzchar(inv$roles[i])) paste0("roles: ", inv$roles[i]) else NULL,
      if (!is.na(inv$mean_identity[i])) paste0("identity: ", round(inv$mean_identity[i] * 100, 1), "%") else NULL,
      if (inv$n_structure[i] > 0) paste0("structure: ", inv$n_structure[i]) else NULL,
      if (inv$n_typeiii_complete[i] > 0) paste0("Type III operon: ", inv$n_typeiii_complete[i]) else NULL
    ), collapse = " | ")
	    if (nzchar(inv$rec_seqs[i])) paste0(line1, "\n", inv$rec_seqs[i]) else line1
  }, character(1))

	  unassociated_note <- if (n_unassociated > 0L) {
	    paste0(
	      " | ", n_unassociated, " unassociated candidate",
	      ifelse(n_unassociated == 1L, "", "s"), " shown in detailed panels",
	      if (n_orphan_mtase > 0L) paste0(" (", n_orphan_mtase, " orphan DNA MTase)") else ""
	    )
	  } else {
	    ""
	  }

	  ggplot2::ggplot(inv, ggplot2::aes(y = .data$REBASEfinder_family_id)) +
	    ggplot2::geom_col(
	      data = inv_bar,
	      ggplot2::aes(x = .data$n, fill = .data$REBASEfinder_family_id, alpha = .data$segment),
	      width = 0.3, color = "grey40", linewidth = 0.2, show.legend = FALSE
	    ) +
	    ggplot2::geom_text(ggplot2::aes(x = 0.05, label = .data$REBASEfinder_family_id),
	                       hjust = 0, size = 3.2, fontface = "bold", color = "white") +
    ggplot2::geom_text(ggplot2::aes(x = .data$n_genes + 0.15, label = .data$annot),
                       hjust = 0, size = 2.5, color = "grey30") +
	    ggplot2::scale_fill_manual(values = palette) +
	    ggplot2::scale_alpha_manual(values = c(meaningful = 0.9, weak = 0.25)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.65)),
                                breaks = seq(0, max(inv$n_genes), by = max(1, ceiling(max(inv$n_genes) / 10)))) +
	    ggplot2::labs(title = paste0("R-M system inventory (", n_verified, " associated loci with complete functional suite / ",
	                                  n_hc, " high-conf / ", n_total, " associated loci)",
	                                  unassociated_note),
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
  has_tier <- "REBASEfinder_curation_tier" %in% names(rm_tbl)
  confidence <- if (has_tier) {
    ifelse(rm_tbl$REBASEfinder_curation_tier == "high", "High",
      ifelse(rm_tbl$REBASEfinder_curation_tier == "medium", "Medium", "Low"))
  } else if (has_eligible) {
    ifelse(!is.na(rm_tbl$REBASEfinder_typing_eligible) & rm_tbl$REBASEfinder_typing_eligible == TRUE, "High", "Low")
  } else {
    rep("High", nrow(rm_tbl))
  }

  has_identity <- "REBASEfinder_blast_identity" %in% names(rm_tbl)
  blast_id <- if (has_identity) rm_tbl$REBASEfinder_blast_identity else rep(NA_real_, nrow(rm_tbl))

  rm_genes <- data.frame(contig = rm_tbl$contig, contig_short = contig_map[rm_tbl$contig],
                          start = rm_tbl$start, end = rm_tbl$end,
                          rm_type = as.character(rm_tbl$REBASEfinder_family_id),
                          identity = blast_id,
                          confidence = confidence, stringsAsFactors = FALSE)

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
  has_ref_rec <- "REBASEfinder_reference_rec_seq" %in% names(rm_tbl)
  short_label <- if (has_hit) {
    hit_lab <- .dnmb_rebasefinder_pretty_hit_label(rm_tbl$REBASEfinder_hit_label)
    sl <- sub("ORF[0-9]+P$", "", hit_lab)
    ifelse(nzchar(sl), sl, hit_lab)
  } else as.character(rm_tbl$REBASEfinder_family_id)
  # Add locus_tag, then recognition sequence at the bottom
  short_label <- paste0(short_label, "\n", rm_tbl$locus_tag)
  rec <- if (has_rec) as.character(rm_tbl$REBASEfinder_rec_seq) else rep(NA_character_, nrow(rm_tbl))
  ref_rec <- if (has_ref_rec) as.character(rm_tbl$REBASEfinder_reference_rec_seq) else rep(NA_character_, nrow(rm_tbl))
  valid_rec <- !is.na(rec) & nzchar(rec) & rec != "?"
  valid_ref_rec <- !is.na(ref_rec) & nzchar(ref_rec) & ref_rec != "?"
  short_label <- ifelse(valid_rec, paste0(short_label, "\n(", rec, ")"), short_label)
  short_label <- ifelse(!valid_rec & valid_ref_rec, paste0(short_label, "\n(ref: ", ref_rec, ")"), short_label)
  if ("REBASEfinder_partial_status" %in% names(rm_tbl)) {
    is_partial <- !is.na(rm_tbl$REBASEfinder_partial_status) &
      rm_tbl$REBASEfinder_partial_status == "partial_or_short"
    len <- if ("REBASEfinder_aa_len" %in% names(rm_tbl)) suppressWarnings(as.integer(rm_tbl$REBASEfinder_aa_len)) else NA_integer_
    partial_label <- ifelse(!is.na(len), paste0("partial/short ", len, " aa"), "partial/short")
    short_label <- ifelse(is_partial, paste0(short_label, "\n", partial_label), short_label)
  }
  rm_genes$label <- short_label
  rm_hc <- rm_genes[rm_genes$confidence == "High", , drop = FALSE]
  rm_lc <- rm_genes[rm_genes$confidence != "High", , drop = FALSE]

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
	    p <- p + ggplot2::geom_rect(data = rm_lc, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = .data$track - 0.12, ymax = .data$track + 0.12), fill = rm_lc$fill_color, color = "grey65", linewidth = 0.15, linetype = "dashed")
	    p <- p + ggrepel::geom_label_repel(
	      data = rm_lc,
	      ggplot2::aes(x = .data$midpoint, y = .data$track + 0.18, label = .data$label),
	      size = 1.7, color = "grey55", fontface = "italic", lineheight = 0.9,
	      fill = scales::alpha("white", 0.88), label.size = 0.08,
	      label.padding = grid::unit(0.06, "lines"),
	      direction = "both", nudge_y = 0.10,
	      segment.size = 0.10, segment.color = "grey65",
	      max.overlaps = Inf, seed = 42, min.segment.length = 0.1,
	      force = 2, force_pull = 0.5, box.padding = 0.25
	    )
	  }
	  if (nrow(rm_hc)) {
	    p <- p + ggplot2::geom_rect(data = rm_hc, ggplot2::aes(xmin = .data$start, xmax = .data$end, ymin = .data$track - 0.12, ymax = .data$track + 0.12), fill = rm_hc$fill_color, color = NA, linewidth = 0) +
	      ggplot2::geom_segment(data = rm_hc, ggplot2::aes(x = .data$midpoint, xend = .data$midpoint, y = .data$track + 0.12, yend = .data$track + 0.20, color = .data$rm_type), linewidth = 0.5, show.legend = FALSE)
	    p <- p + ggrepel::geom_label_repel(
	      data = rm_hc,
	      ggplot2::aes(x = .data$midpoint, y = .data$track + 0.22, label = .data$label),
	      size = 1.9, color = "grey20", fontface = "bold", lineheight = 0.9,
	      fill = scales::alpha("white", 0.90), label.size = 0.08,
	      label.padding = grid::unit(0.06, "lines"),
	      direction = "both", nudge_y = 0.12,
	      segment.size = 0.12, segment.color = "grey55",
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
    ggplot2::labs(title = "R-M system genome context", subtitle = "solid = high-confidence, dashed = medium/review | color intensity = BLAST identity", x = "Position", y = NULL) +
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
.dnmb_plot_rebasefinder_blast_quality <- function(tbl, role_palette, display_info, rm_palette = NULL, cache_root = NULL) {
  has_identity <- "REBASEfinder_blast_identity" %in% names(tbl)
  has_bitscore <- "REBASEfinder_blast_bitscore" %in% names(tbl)

  if (!has_identity) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No BLAST identity data", size = 4, color = "grey50"))

  plot_tbl <- merge(display_info, data.frame(
    locus_tag = tbl$locus_tag,
	    REBASEfinder_hit_label = if ("REBASEfinder_hit_label" %in% names(tbl)) as.character(tbl$REBASEfinder_hit_label) else NA_character_,
	    bitscore = if (has_bitscore) tbl$REBASEfinder_blast_bitscore else NA_real_,
	    structure_status = if ("REBASEfinder_structure_status" %in% names(tbl)) as.character(tbl$REBASEfinder_structure_status) else NA_character_,
	    structure_pass = if ("REBASEfinder_structure_pass" %in% names(tbl)) as.logical(tbl$REBASEfinder_structure_pass) else NA,
	    structure_file_exists = if ("REBASEfinder_structure_file_exists" %in% names(tbl)) as.logical(tbl$REBASEfinder_structure_file_exists) else NA,
	    foldseek_hit_present = if ("REBASEfinder_foldseek_hit_present" %in% names(tbl)) as.logical(tbl$REBASEfinder_foldseek_hit_present) else NA,
	    structure_coverage_status = if ("REBASEfinder_structure_coverage_status" %in% names(tbl)) as.character(tbl$REBASEfinder_structure_coverage_status) else NA_character_,
	    typeiii_context_status = if ("REBASEfinder_typeiii_context_status" %in% names(tbl)) as.character(tbl$REBASEfinder_typeiii_context_status) else NA_character_,
	    stringsAsFactors = FALSE
	  ), by = "locus_tag", sort = FALSE)
  display_lvls <- .dnmb_rebasefinder_display_levels(display_info)
  plot_tbl$display <- factor(as.character(plot_tbl$display), levels = display_lvls)
  # Use role from display_info (already included)
  if (!"role" %in% names(plot_tbl) || all(is.na(plot_tbl$role))) {
    plot_tbl$role <- "Other"
  }
  plot_tbl$role[is.na(plot_tbl$role) | !nzchar(plot_tbl$role)] <- "Other"
  if (!nrow(plot_tbl)) return(ggplot2::ggplot() + ggplot2::theme_void())

  # Infer methylation type for M subunits + recognition sequence
  meth_annot <- .dnmb_rebasefinder_methylation_annotations(tbl, cache_root = cache_root)
  ma_idx <- match(plot_tbl$locus_tag, meth_annot$locus_tag)
  plot_tbl$meth_type <- meth_annot$meth_type[ma_idx]
  plot_tbl$meth_pos <- meth_annot$meth_pos[ma_idx]
  plot_tbl$meth_all <- meth_annot$meth_all[ma_idx]
  plot_tbl$meth_source <- meth_annot$meth_source[ma_idx]
  plot_tbl$rec_seq <- meth_annot$rec_seq[ma_idx]

  # Ensure rm_palette is available
  if (is.null(rm_palette)) rm_palette <- .dnmb_rebasefinder_palette(plot_tbl$rm_type)
  all_rm <- unique(plot_tbl$rm_type[!is.na(plot_tbl$rm_type)])
  miss_rm <- setdiff(all_rm, names(rm_palette))
  if (length(miss_rm)) rm_palette <- c(rm_palette, stats::setNames(rep("grey60", length(miss_rm)), miss_rm))

  # Build y-axis labels with rm_type prefix
  type_label <- ifelse(!is.na(plot_tbl$rm_type), paste0("[", plot_tbl$rm_type, "] "), "")
  plot_tbl$display_typed <- paste0(type_label, as.character(plot_tbl$display))
  # Preserve factor order from display
  lvls <- display_lvls
  type_map <- stats::setNames(plot_tbl$display_typed, as.character(plot_tbl$display))
  new_lvls <- type_map[lvls]
  plot_tbl$display_typed <- factor(plot_tbl$display_typed, levels = new_lvls)

  # Color methylated bases in y-axis labels (rec_seq portion)
  # Uses REBASE bairoch metadata for exact position when available,
  # falls back to coloring ALL matching bases when position is unknown.
  meth_pal <- c("N6A" = "#D32F2F", "N5C" = "#1565C0", "N4C" = "#FF8F00")
  meth_target <- c("N6A" = "A", "N5C" = "C", "N4C" = "C")
  has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
  if (has_ggtext) {
    old_lvls <- levels(plot_tbl$display_typed)
    label_vec <- as.character(plot_tbl$display_typed)
    lvl_map <- stats::setNames(old_lvls, old_lvls)
    for (i in seq_len(nrow(plot_tbl))) {
      mt <- plot_tbl$meth_type[i]
      rs <- plot_tbl$rec_seq[i]
      col <- meth_pal[mt]
      if (is.na(col)) next
      # Try position-specific coloring from bairoch
      colored_rs <- NULL
      if (!is.na(rs) && nzchar(rs) && rs != "?") {
        mpos <- plot_tbl$meth_pos[i]
        if (!is.na(mpos) && mpos != "?" && grepl("^-?[0-9]+$", mpos)) {
          pos <- as.integer(mpos)
          if (pos < 0) pos <- nchar(rs) + pos + 1L
          if (pos >= 1 && pos <= nchar(rs)) {
            chars <- strsplit(rs, "")[[1]]
            chars[pos] <- paste0(
              "<span style='color:", col, ";font-weight:bold'>",
              chars[pos], "</span>")
            colored_rs <- paste(chars, collapse = "")
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
      }

      old_label <- label_vec[i]
      new_label <- old_label
      if (!is.null(colored_rs) && !is.na(rs) && nzchar(rs)) {
        new_label <- sub(rs, colored_rs, new_label, fixed = TRUE)
      }
      label_vec[i] <- new_label
      lvl_map[lvl_map == old_label] <- new_label
	    }
	    new_lvls <- unname(lvl_map[old_lvls])
	    label_vec <- gsub("\n", "<br>", label_vec, fixed = TRUE)
	    new_lvls <- gsub("\n", "<br>", new_lvls, fixed = TRUE)
	    plot_tbl$display_typed <- factor(label_vec, levels = new_lvls)
	  }

  # Operon group separators
  operon_bounds <- .dnmb_rebasefinder_operon_separators(plot_tbl, y_col = "display_typed")
  operon_blocks <- .dnmb_rebasefinder_operon_blocks(plot_tbl, y_col = "display_typed")

  all_roles <- unique(plot_tbl$role)
  missing_roles <- setdiff(all_roles, names(role_palette))
  if (length(missing_roles)) role_palette <- c(role_palette, stats::setNames(rep("grey60", length(missing_roles)), missing_roles))

  shade_df <- data.frame(xmin = -Inf, xmax = 50, ymin = -Inf, ymax = Inf)

  blast_tbl <- plot_tbl[!is.na(plot_tbl$identity), , drop = FALSE]
  context_only_tbl <- plot_tbl[is.na(plot_tbl$identity), , drop = FALSE]

  p <- ggplot2::ggplot(plot_tbl, ggplot2::aes(y = .data$display_typed)) +
    ggplot2::geom_blank(ggplot2::aes(x = 0))
  if (nrow(operon_blocks)) {
    p <- p + ggplot2::geom_rect(
      data = operon_blocks,
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = .data$ymin, ymax = .data$ymax),
      fill = operon_blocks$fill, alpha = 0.55, color = NA,
      inherit.aes = FALSE
    )
  }
  p <- p +
    ggplot2::geom_rect(data = shade_df, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax), fill = "#FFF3E0", alpha = 0.35, inherit.aes = FALSE) +
    ggplot2::geom_vline(xintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.4)

  # Operon group separators (horizontal lines)
  if (nrow(operon_bounds)) {
    p <- p + ggplot2::geom_hline(yintercept = operon_bounds$y, color = "grey75", linewidth = 0.3, linetype = "dotted")
  }

  if (has_bitscore && nrow(blast_tbl) && any(!is.na(blast_tbl$bitscore))) {
    p <- p + ggplot2::geom_point(
      data = blast_tbl,
      ggplot2::aes(size = .data$bitscore, fill = .data$role,
                   color = .data$rm_type, shape = .data$confidence,
                   x = .data$identity * 100),
      alpha = 0.8, stroke = 1.0
    ) +
      ggplot2::scale_size_continuous(name = "Bitscore", range = c(2.5, 7),
                                     breaks = scales::pretty_breaks(n = 3))
  } else if (nrow(blast_tbl)) {
    p <- p + ggplot2::geom_point(
      data = blast_tbl,
      ggplot2::aes(fill = .data$role, color = .data$rm_type,
                   shape = .data$confidence, x = .data$identity * 100),
      size = 4, alpha = 0.8, stroke = 1.0
    )
  }
  if (nrow(context_only_tbl)) {
    p <- p + ggplot2::geom_point(
      data = context_only_tbl,
      ggplot2::aes(x = 0, y = .data$display_typed),
      inherit.aes = FALSE,
      shape = 124, size = 2.5, stroke = 0.7, color = "grey60",
      show.legend = FALSE
    )
  }
	  plot_tbl$foldseek_supported <- grepl("^foldseek_supported", plot_tbl$structure_coverage_status) |
	    (grepl("^structure_supported", plot_tbl$structure_status) & plot_tbl$structure_pass %in% TRUE)
	  plot_tbl$foldseek_supported[is.na(plot_tbl$foldseek_supported)] <- FALSE
	  struct_tbl <- plot_tbl[plot_tbl$foldseek_supported & !is.na(plot_tbl$identity), , drop = FALSE]
  if (nrow(struct_tbl)) {
    p <- p + ggplot2::geom_point(
      data = struct_tbl,
      ggplot2::aes(x = .data$identity * 100, y = .data$display_typed),
      inherit.aes = FALSE,
      shape = 8, size = 2.2, stroke = 0.55, color = "#111827",
      show.legend = FALSE
    )
  }

  x_left_limit <- if (nrow(context_only_tbl)) {
    0
  } else if (nrow(blast_tbl)) {
    max(0, min(blast_tbl$identity * 100, na.rm = TRUE) - 5)
  } else {
    0
  }

	  meth_annot_label <- if (any(!is.na(plot_tbl$meth_type))) {
	    if (has_ggtext) {
	      paste0("Methylated base: <span style='color:#D32F2F'>N6A</span>",
	             " / <span style='color:#1565C0'>N5C</span>",
	             " / <span style='color:#FF8F00'>N4C</span>")
	    } else {
	      "Methylated base: N6A / N5C / N4C"
	    }
	  } else NULL
	  has_structure_tracking <- "structure_file_exists" %in% names(plot_tbl) &&
	    any(!is.na(plot_tbl$structure_file_exists) | !is.na(plot_tbl$foldseek_hit_present) |
	          (!is.na(plot_tbl$structure_coverage_status) & nzchar(plot_tbl$structure_coverage_status)))
	  foldseek_summary <- if (has_structure_tracking) {
	    n_files <- sum(plot_tbl$structure_file_exists %in% TRUE, na.rm = TRUE)
	    n_foldseek <- sum(plot_tbl$foldseek_hit_present %in% TRUE, na.rm = TRUE)
	    n_supported <- sum(plot_tbl$foldseek_supported %in% TRUE, na.rm = TRUE)
	    paste0("Foldseek: ", n_foldseek, "/", nrow(plot_tbl),
	           " checked, ", n_supported, " supported")
	  } else if (nrow(struct_tbl)) {
	    "star = Foldseek-supported"
	  } else {
	    NULL
	  }
		  subtitle_parts <- c(
		    meth_annot_label,
		    foldseek_summary,
		    if (nrow(context_only_tbl)) "tick = operon/annotation-only" else NULL
		  )
	  struct_subtitle <- if (length(subtitle_parts)) {
	    paste(subtitle_parts, collapse = if (has_ggtext) "<br>" else "\n")
	  } else NULL

  panel_core <- p + ggplot2::scale_shape_manual(values = c(High = 21, Medium = 22, Low = 1), name = "Confidence") +
    ggplot2::annotate("text", x = 49, y = Inf, label = "low", hjust = 1, vjust = 1.5,
                      size = 2.3, color = "#BF360C", fontface = "italic") +
    ggplot2::annotate("text", x = 51, y = Inf, label = "high confidence", hjust = 0, vjust = 1.5,
                      size = 2.3, color = "#2E7D32", fontface = "italic") +
    ggplot2::scale_fill_manual(values = role_palette, name = "Enzyme Role") +
    ggplot2::scale_color_manual(values = rm_palette, name = "R-M Type") +
    ggplot2::scale_x_continuous(limits = c(x_left_limit, 100),
                                breaks = seq(0, 100, by = 10)) +
    ggplot2::scale_y_discrete(limits = new_lvls, drop = FALSE) +
    ggplot2::labs(
      title = "REBASE BLAST",
      subtitle = struct_subtitle,
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
      axis.text.y = if (has_ggtext) ggtext::element_markdown(size = 8.3) else ggplot2::element_text(size = 8.3),
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

  panel_core
}

# Compute operon group separator y-positions (between adjacent different operons)
.dnmb_rebasefinder_operon_separators <- function(plot_tbl, y_col = "display") {
  if (!"operon" %in% names(plot_tbl) || !y_col %in% names(plot_tbl) || nrow(plot_tbl) < 2) {
    return(data.frame(y = numeric(0)))
  }
  y_vals <- plot_tbl[[y_col]]
  lvls <- levels(y_vals)
  if (is.null(lvls)) lvls <- unique(as.character(y_vals))
  lvls <- lvls[!is.na(lvls)]
  ordered_operons <- plot_tbl$operon[match(lvls, as.character(y_vals))]
  seps <- which(ordered_operons[-1] != ordered_operons[-length(ordered_operons)])
  if (!length(seps)) return(data.frame(y = numeric(0)))
  data.frame(y = seps + 0.5)
}

.dnmb_rebasefinder_operon_blocks <- function(plot_tbl, y_col = "display") {
  if (!"operon" %in% names(plot_tbl) || !y_col %in% names(plot_tbl) || nrow(plot_tbl) < 2) {
    return(data.frame())
  }
  y_vals <- plot_tbl[[y_col]]
  lvls <- levels(y_vals)
  if (is.null(lvls)) lvls <- unique(as.character(y_vals))
  lvls <- lvls[!is.na(lvls)]
  ordered_operons <- plot_tbl$operon[match(lvls, as.character(y_vals))]
  if (!length(ordered_operons)) return(data.frame())
  run <- rle(ordered_operons)
  ends <- cumsum(run$lengths)
  starts <- ends - run$lengths + 1L
  out <- data.frame(
    operon = run$values,
    ymin = starts - 0.5,
    ymax = ends + 0.5,
    block = seq_along(run$values),
    stringsAsFactors = FALSE
  )
  out <- out[out$ymax - out$ymin > 1, , drop = FALSE]
  if (!nrow(out)) return(out)
  out$fill <- ifelse(out$block %% 2L == 0L, "#F7FAFC", "#EEF6FF")
  out
}

.dnmb_rebasefinder_meth_type_label <- function(meth_type) {
  meth_type <- as.character(meth_type)
  ifelse(grepl("m6A", meth_type), "N6A",
    ifelse(grepl("m5C", meth_type), "N5C",
      ifelse(grepl("m4C|Nm4C", meth_type), "N4C", NA_character_)
    )
  )
}

.dnmb_rebasefinder_curated_mtase_motif_verified <- function(tbl) {
  n <- base::nrow(tbl)
  out <- base::rep(NA, n)
  if ("REBASEfinder_mtase_motif_verified" %in% base::names(tbl)) {
    out <- suppressWarnings(base::as.logical(tbl$REBASEfinder_mtase_motif_verified))
  }
  missing <- base::which(base::is.na(out))
  if (!base::length(missing)) return(out)

  translations <- if ("translation" %in% base::names(tbl)) {
    base::as.character(tbl$translation)
  } else {
    base::rep(NA_character_, n)
  }
  annotations <- .dnmb_rebasefinder_gene_annotation_text(tbl)
  out[missing] <- base::vapply(missing, function(i) {
    .dnmb_rebasefinder_mtase_sequence_signals(translations[[i]], annotations[[i]])$verified
  }, logical(1))
  out
}

.dnmb_rebasefinder_methylation_annotations <- function(tbl, cache_root = NULL) {
  n <- nrow(tbl)
  rec_seq <- if ("REBASEfinder_rec_seq" %in% names(tbl)) as.character(tbl$REBASEfinder_rec_seq) else rep(NA_character_, n)
  reference_rec_seq <- if ("REBASEfinder_reference_rec_seq" %in% names(tbl)) {
    as.character(tbl$REBASEfinder_reference_rec_seq)
  } else {
    rep(NA_character_, n)
  }
  rec_seq[rec_seq %in% c("", "?", "NA")] <- NA_character_
  reference_rec_seq[reference_rec_seq %in% c("", "?", "NA")] <- NA_character_
  use_reference <- is.na(rec_seq) & !is.na(reference_rec_seq)
  rec_seq[use_reference] <- reference_rec_seq[use_reference]
  out <- data.frame(
    locus_tag = if ("locus_tag" %in% names(tbl)) as.character(tbl$locus_tag) else rep(NA_character_, n),
    meth_type = rep(NA_character_, n),
    meth_pos = rep(NA_character_, n),
    meth_all = rep(NA_character_, n),
    rec_seq = rec_seq,
    rec_source = ifelse(use_reference, "rebase_reference", ifelse(!is.na(rec_seq), "query_prediction", NA_character_)),
    meth_source = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
  if (!n) return(out)

  bairoch <- tryCatch(.dnmb_rebasefinder_download_bairoch(cache_root = cache_root), error = function(e) NULL)
  hit <- if ("REBASEfinder_hit_label" %in% names(tbl)) as.character(tbl$REBASEfinder_hit_label) else rep(NA_character_, n)
  if (!is.null(bairoch) && nrow(bairoch)) {
    donor <- if ("REBASEfinder_recognition_donor" %in% names(tbl)) {
      as.character(tbl$REBASEfinder_recognition_donor)
    } else {
      rep(NA_character_, n)
    }
    recognition_match <- if ("REBASEfinder_recognition_match" %in% names(tbl)) {
      as.character(tbl$REBASEfinder_recognition_match)
    } else {
      rep(NA_character_, n)
    }
    bairoch_name <- trimws(as.character(bairoch$enzyme_name))
    bairoch_rec <- as.character(bairoch$rec_seq)
    bairoch_type <- as.character(bairoch$meth_type)
    bairoch_pos <- as.character(bairoch$meth_pos)
    bairoch_all <- as.character(bairoch$meth_all)

    for (i in seq_len(n)) {
      candidate_names <- c(donor[[i]], recognition_match[[i]], hit[[i]])
      candidate_sources <- c("recognition_donor", "recognition_match", "clean_hit")
      candidate_names <- .dnmb_rebasefinder_clean_rebase_subject(candidate_names)
      valid_candidate <- !is.na(candidate_names) & nzchar(candidate_names) &
        !grepl("_context:", candidate_names, fixed = TRUE)
      candidate_names <- candidate_names[valid_candidate]
      candidate_sources <- candidate_sources[valid_candidate]
      if (!length(candidate_names)) next

      for (j in seq_along(candidate_names)) {
        exact <- which(bairoch_name == candidate_names[[j]])
        if (!length(exact)) next
        exact_rec <- exact[!is.na(bairoch_rec[exact]) & nzchar(bairoch_rec[exact])]
        if (length(exact_rec) && (is.na(out$rec_seq[[i]]) || !nzchar(out$rec_seq[[i]]))) {
          out$rec_seq[[i]] <- bairoch_rec[[exact_rec[[1]]]]
          out$rec_source[[i]] <- paste0("rebase_", candidate_sources[[j]])
        }
        exact_meth <- exact[!is.na(bairoch_type[exact]) & nzchar(bairoch_type[exact])]
        if (!length(exact_meth)) next
        exact_meth <- exact_meth[[1]]
        out$meth_type[[i]] <- .dnmb_rebasefinder_meth_type_label(bairoch_type[[exact_meth]])
        out$meth_pos[[i]] <- bairoch_pos[[exact_meth]]
        out$meth_all[[i]] <- bairoch_all[[exact_meth]]
        out$meth_source[[i]] <- paste0("bairoch_", candidate_sources[[j]])
        break
      }
    }

    rec_known <- !is.na(out$rec_seq) & nzchar(out$rec_seq)
    rec_missing_type <- rec_known & is.na(out$meth_type)
    if (any(rec_missing_type)) {
      for (i in which(rec_missing_type)) {
        idx <- which(!is.na(bairoch_rec) & bairoch_rec == out$rec_seq[[i]] &
                       !is.na(bairoch_type) & nzchar(bairoch_type))
        if (!length(idx)) next
        inferred_type <- .dnmb_rebasefinder_meth_type_label(bairoch_type[idx])
        inferred_pos <- trimws(bairoch_pos[idx])
        inferred_pos[is.na(inferred_pos) | !nzchar(inferred_pos)] <- NA_character_
        valid <- !is.na(inferred_type) & nzchar(inferred_type)
        if (!any(valid)) next
        idx <- idx[valid]
        inferred_type <- inferred_type[valid]
        inferred_pos <- inferred_pos[valid]
        signatures <- paste(inferred_type, ifelse(is.na(inferred_pos), "<NA>", inferred_pos), sep = "\r")
        if (length(unique(signatures)) != 1L) next
        chosen <- idx[[1]]
        out$meth_type[[i]] <- inferred_type[[1]]
        out$meth_pos[[i]] <- inferred_pos[[1]]
        out$meth_all[[i]] <- bairoch_all[[chosen]]
        out$meth_source[[i]] <- "bairoch_rec_seq_unique"
      }
    }
  }

  group_info <- .dnmb_rebasefinder_operon_groups_for_plot(tbl)
  group <- group_info$operon_group[match(out$locus_tag, group_info$locus_tag)]
  if (length(group) == n) {
    for (g in unique(group[!is.na(group)])) {
      idx <- which(group == g)
      donors <- idx[!is.na(out$meth_type[idx])]
      if (!length(donors)) next
      donor <- donors[[1]]
      missing_type <- idx[is.na(out$meth_type[idx])]
      if (length(missing_type)) {
        out$meth_type[missing_type] <- out$meth_type[[donor]]
        out$meth_pos[missing_type] <- out$meth_pos[[donor]]
        out$meth_all[missing_type] <- out$meth_all[[donor]]
        out$meth_source[missing_type] <- paste0("operon_", out$meth_source[[donor]])
      }
      missing_rec <- idx[is.na(out$rec_seq[idx]) | !nzchar(out$rec_seq[idx])]
      if (length(missing_rec) && !is.na(out$rec_seq[[donor]]) && nzchar(out$rec_seq[[donor]])) {
        out$rec_seq[missing_rec] <- out$rec_seq[[donor]]
        out$rec_source[missing_rec] <- paste0("operon_", out$rec_source[[donor]])
      }
    }
  }

  curated_mtase_pair <- .dnmb_rebasefinder_curated_mtase_motif_verified(tbl)
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  if (!is.null(detailed)) {
    motif_names <- names(.dnmb_rebasefinder_motif_definitions())
    c5_idx <- which(motif_names == "N5C-PC")
    amino_idx <- which(motif_names == "Amino-IV")
    for (i in seq_len(n)) {
      if (!is.na(out$meth_type[[i]])) next
      c5_hit <- length(c5_idx) && detailed[[i]][[c5_idx]]$status %in% c("present", "present*")
      amino_hit <- length(amino_idx) && detailed[[i]][[amino_idx]]$status %in% c("present", "present*")
      mtase_pair <- curated_mtase_pair[[i]]
      if (isTRUE(c5_hit) && isTRUE(mtase_pair)) {
        out$meth_type[[i]] <- "N5C"
        out$meth_pos[[i]] <- "?"
        out$meth_source[[i]] <- "motif_fallback"
      } else if (isTRUE(amino_hit) && !is.na(out$rec_seq[[i]]) && nzchar(out$rec_seq[[i]])) {
        bases <- toupper(gsub("[^ACGTWSMKRYBDHVN]", "", out$rec_seq[[i]]))
        has_a <- grepl("A", bases)
        has_c <- grepl("C", bases)
        if (has_a && !has_c) {
          out$meth_type[[i]] <- "N6A"; out$meth_pos[[i]] <- "?"
          out$meth_source[[i]] <- "motif_recseq_fallback"
        } else if (has_c && !has_a) {
          out$meth_type[[i]] <- "N4C"; out$meth_pos[[i]] <- "?"
          out$meth_source[[i]] <- "motif_recseq_fallback"
        }
      }
    }
  }
  out
}

.dnmb_rebasefinder_methylation_tag <- function(meth_type, meth_pos) {
  meth_type <- as.character(meth_type)
  meth_pos <- as.character(meth_pos)
  has <- !is.na(meth_type) & nzchar(meth_type)
  out <- rep("", length(meth_type))
  out[has] <- paste0(meth_type[has], "@", ifelse(!is.na(meth_pos[has]) & nzchar(meth_pos[has]), meth_pos[has], "?"))
  out
}


# ====================================================================
# Motif definitions & scanning
# Literature-based patterns with positional constraints & co-occurrence
# Refs: Malone+ 1995 (amino-MTase motifs), Posfai+ 1989 / Kumar+ 1994
#   (C5-MTase I-X), Pingoud & Jeltsch 2001 (PD-ExK), Dunin-Horkawicz+ 2006
#   (GIY-YIG), Schluckebier+ 1995 (SAM fold), Murray 2000 (Type I)
# ====================================================================
.dnmb_rebasefinder_motif_definitions <- function() {
  .dnmb_rebasefinder_motif_catalog()
}

.dnmb_rebasefinder_motif_roles <- function(def) {
  roles <- def$expected_role
  roles <- base::as.character(roles)
  roles[!base::is.na(roles) & base::nzchar(roles)]
}

.dnmb_rebasefinder_motif_role_label <- function(def) {
  roles <- .dnmb_rebasefinder_motif_roles(def)
  if (!base::length(roles)) NA_character_ else base::paste(roles, collapse = "/")
}

.dnmb_rebasefinder_motif_families <- function(def) {
  fam <- def$expected_family
  fam <- base::as.character(fam)
  fam[!base::is.na(fam) & base::nzchar(fam)]
}

.dnmb_rebasefinder_motif_family_match <- function(def, family) {
  families <- .dnmb_rebasefinder_motif_families(def)
  if (!base::length(families)) return(TRUE)
  family <- base::as.character(family)[1]
  if (base::is.na(family) || !base::nzchar(family)) return(FALSE)
  family %in% families
}

.dnmb_rebasefinder_motif_role_match <- function(def, role, family = NA_character_) {
  role <- base::as.character(role)[1]
  if (base::is.na(role) || !base::nzchar(role)) return(FALSE)
  expected <- .dnmb_rebasefinder_motif_roles(def)
  role_match <- role %in% expected ||
    (role == "RM" && base::any(expected %in% c("M", "R")))
  role_match && .dnmb_rebasefinder_motif_family_match(def, family)
}

.dnmb_rebasefinder_motif_primary_role <- function(def) {
  if (!base::is.null(def$axis) && !base::is.na(def$axis) && base::nzchar(def$axis)) {
    return(base::as.character(def$axis)[1])
  }
  roles <- .dnmb_rebasefinder_motif_roles(def)
  if (!base::length(roles)) NA_character_ else roles[[1]]
}

# Backward-compatible alias for legacy panel code that references old names
.dnmb_rebasefinder_motif_names_display <- function() {
  names(.dnmb_rebasefinder_motif_definitions())
}

# ====================================================================
# Positional constraint check — does a hit fall within expected region?
# ====================================================================
.dnmb_rebasefinder_pos_in_range <- function(pos, prot_len, pos_range) {
  if (is.null(pos_range) || is.na(pos) || prot_len <= 0) return(TRUE)
  frac <- pos / prot_len
  frac >= pos_range[1] && frac <= pos_range[2]
}

.dnmb_rebasefinder_motif_info_positions <- function(gene_result, name, in_range_only = FALSE) {
  info <- gene_result[[name]]
  if (base::is.null(info) || !base::length(info$hits)) return(integer())
  hits <- info$hits
  if (isTRUE(in_range_only)) hits <- base::Filter(function(hit) isTRUE(hit$in_range), hits)
  if (!base::length(hits)) integer() else base::vapply(hits, `[[`, integer(1), "pos")
}

.dnmb_rebasefinder_motif_text_count <- function(text, pattern) {
  text <- base::tolower(base::as.character(text)[1])
  if (base::is.na(text) || !base::nzchar(text)) return(0L)
  hit <- base::gregexpr(pattern, text, perl = TRUE)[[1]]
  if (!base::length(hit) || hit[[1]] < 1L) 0L else base::length(hit)
}

.dnmb_rebasefinder_motif_assessment <- function(gene_result,
                                                 role,
                                                 family,
                                                 annotation = NA_character_) {
  definitions <- .dnmb_rebasefinder_motif_definitions()
  motif_names <- base::names(definitions)
  supported <- stats::setNames(base::rep(FALSE, base::length(motif_names)), motif_names)
  hint <- supported
  derived_label <- stats::setNames(base::rep(NA_character_, base::length(motif_names)), motif_names)
  applicable <- stats::setNames(base::vapply(definitions, function(definition) {
    .dnmb_rebasefinder_motif_role_match(definition, role, family)
  }, logical(1)), motif_names)
  present <- stats::setNames(base::vapply(motif_names, function(name) {
    info <- gene_result[[name]]
    !base::is.null(info) && !base::is.null(info$n_hits) && info$n_hits > 0L
  }, logical(1)), motif_names)

  annotation <- base::as.character(annotation)[1]
  if (base::is.na(annotation)) annotation <- ""
  annotation_lower <- base::tolower(annotation)
  flags <- .dnmb_rebasefinder_annotation_flags(annotation)
  dna_mtase <- isTRUE(flags$dna_mtase[[1]]) && !isTRUE(flags$explicit_non_dna_mtase[[1]])
  c5_domain <- base::grepl(
    "pf00145|c5[-_ ]?dna[-_ ]?methyl|cytosine-5.*dna.*methyl",
    annotation_lower,
    perl = TRUE
  )
  rm_rease <- isTRUE(flags$rm_rease[[1]]) &&
    !base::any(c(flags$repair_nuclease[[1]], flags$unrelated_helicase[[1]],
                 flags$mobile_nuclease[[1]], flags$other_defense[[1]]))

  ordered <- function(names, max_span = Inf, in_range_only = FALSE) {
    .dnmb_rebasefinder_motif_ordered(
      base::lapply(names, function(name) {
        .dnmb_rebasefinder_motif_info_positions(gene_result, name, in_range_only = in_range_only)
      }),
      max_span = max_span
    )
  }

  amino_pair <- present[["SAM"]] && present[["Amino-IV"]]
  c5_pair <- ordered(c("N5C-PC", "N5C-ENV"), max_span = 250L)
  c5_extended <- ordered(c("N5C-PC", "N5C-ENV", "N5C-QxRxR"), max_span = 450L)
  mmei_sequence <- ordered(c("MmeI-X", "MmeI-I", "Amino-IV"), max_span = 650L)
  mmei_domain_count <- base::sum(base::vapply(
    c("pf20465", "pf20466", "pf20467", "pf20473"),
    function(pattern) base::grepl(pattern, annotation_lower, fixed = TRUE),
    logical(1)
  ))
  mmei_architecture <- mmei_domain_count >= 3L || (
    base::grepl("mmei", annotation_lower, fixed = TRUE) &&
      base::grepl("target recognition", annotation_lower, fixed = TRUE) &&
      base::grepl("dna-methyltransferase domain", annotation_lower, fixed = TRUE)
  )

  if (dna_mtase && amino_pair) supported[c("SAM", "Amino-IV")] <- TRUE
  if (dna_mtase && c5_pair) {
    supported[c("N5C-PC", "N5C-ENV")] <- TRUE
    if (c5_extended) supported[["N5C-QxRxR"]] <- TRUE
    if (present[["SAM"]]) supported[["SAM"]] <- TRUE
  }
  if (mmei_architecture) {
    supported[["MmeI-architecture"]] <- TRUE
    derived_label[["MmeI-architecture"]] <- paste0(mmei_domain_count, "-dom")
    for (name in c("MmeI-X", "MmeI-I", "Amino-IV", "MmeI-PDExK")) {
      if (present[[name]]) supported[[name]] <- TRUE
    }
  }
  c5_names <- c("N5C-PC", "N5C-ENV", "N5C-QxRxR")
  if (role %in% c("M", "RM")) {
    if (c5_domain || c5_pair) {
      applicable[["Amino-IV"]] <- FALSE
    } else {
      applicable[c5_names] <- FALSE
    }
    if (mmei_architecture) {
      applicable[c5_names] <- FALSE
      applicable[c("PD-ExK", "HNH", "GIY-YIG", "PLD-HKD")] <- FALSE
    }
  }
  if (dna_mtase) {
    hint[["Amino-IV"]] <- present[["Amino-IV"]] && !supported[["Amino-IV"]]
    hint[["SAM"]] <- present[["SAM"]] && present[["Amino-IV"]] && !supported[["SAM"]]
    if (c5_domain && (present[["N5C-PC"]] || present[["N5C-ENV"]])) {
      hint[["N5C-PC"]] <- present[["N5C-PC"]] && !supported[["N5C-PC"]]
      hint[["N5C-ENV"]] <- present[["N5C-ENV"]] && !supported[["N5C-ENV"]]
    }
  }

  hsd_s_count <- base::max(c(
    .dnmb_rebasefinder_motif_text_count(annotation, "pf01420"),
    .dnmb_rebasefinder_motif_text_count(annotation, "type i restriction modification dna specificity domain"),
    .dnmb_rebasefinder_motif_text_count(annotation, "rmtype1_s_[a-z0-9_-]+")
  ))
  if (identical(role, "S") && identical(family, "Type I")) {
    if (hsd_s_count >= 2L) {
      supported[["HsdS-2TRD"]] <- TRUE
      derived_label[["HsdS-2TRD"]] <- "2x TRD"
    } else if (isTRUE(flags$rm_specificity[[1]])) {
      hint[["HsdS-2TRD"]] <- TRUE
      derived_label[["HsdS-2TRD"]] <- if (hsd_s_count == 1L) "1x TRD" else "HsdS profile"
    }
  }

  hsd_motor <- ordered(c("P-loop", "HsdR-WB", "HsdR-MIII"), max_span = 750L, in_range_only = TRUE)
  typeiii_motor <- ordered(c("ResIII-WA", "ResIII-WB", "ResIII-MIII"), max_span = 750L, in_range_only = TRUE)
  if (identical(role, "R") && identical(family, "Type I") && rm_rease) {
    if (hsd_motor) supported[c("P-loop", "HsdR-WB", "HsdR-MIII")] <- TRUE
    motor_present <- present[c("P-loop", "HsdR-WB", "HsdR-MIII")]
    if (base::sum(motor_present) >= 2L && !hsd_motor) {
      hint[base::names(motor_present)[motor_present]] <- TRUE
    }
    hsd_pd_in_range <- base::length(.dnmb_rebasefinder_motif_info_positions(
      gene_result, "HsdR-PD", in_range_only = TRUE
    )) > 0L
    if (hsd_pd_in_range) supported[["HsdR-PD"]] <- TRUE
  }
  if (identical(role, "R") && identical(family, "Type III") && rm_rease) {
    if (typeiii_motor) supported[c("ResIII-WA", "ResIII-WB", "ResIII-MIII")] <- TRUE
    motor_present <- present[c("ResIII-WA", "ResIII-WB", "ResIII-MIII")]
    if (base::sum(motor_present) >= 2L && !typeiii_motor) {
      hint[base::names(motor_present)[motor_present]] <- TRUE
    }
    typeiii_pd_in_range <- base::length(.dnmb_rebasefinder_motif_info_positions(
      gene_result, "ResIII-PD", in_range_only = TRUE
    )) > 0L
    if (typeiii_pd_in_range) supported[["ResIII-PD"]] <- TRUE
    if (present[["ResIII-PD"]] && !typeiii_pd_in_range) hint[["ResIII-PD"]] <- TRUE
  }

  hnh_domain <- base::grepl("\\bhnh\\b|h-n-h|pf01844|pf13391", annotation_lower, perl = TRUE)
  giy_domain <- base::grepl("giy[-_ ]?yig|pf01541|pf02820", annotation_lower, perl = TRUE)
  pld_domain <- base::grepl("pf13091|pldc_2|pld[- ]like|phospholipase d/nuclease|hkd", annotation_lower, perl = TRUE)
  pld_count <- base::length(.dnmb_rebasefinder_motif_info_positions(gene_result, "PLD-HKD"))
  typeii_role <- role %in% c("R", "RM") && identical(family, "Type II")
  if (typeii_role && rm_rease) {
    if (present[["PD-ExK"]]) supported[["PD-ExK"]] <- TRUE
    if (present[["HNH"]] && hnh_domain) supported[["HNH"]] <- TRUE
    if (present[["GIY-YIG"]] && (giy_domain || !base::any(c(flags$repair_nuclease[[1]], flags$mobile_nuclease[[1]])))) {
      supported[["GIY-YIG"]] <- TRUE
    }
    if (pld_count >= 2L && pld_domain) supported[["PLD-HKD"]] <- TRUE
    if (pld_count == 1L && pld_domain) hint[["PLD-HKD"]] <- TRUE
  }

  typeiv_label <- NA_character_
  typeiv_display <- NA_character_
  if (identical(role, "R") && identical(family, "Type IV")) {
    resiii_pld <- base::grepl("pf04851|resiii|dexhc_re_i_iii_res", annotation_lower, perl = TRUE) && pld_domain
    mspji <- base::grepl("mspji|pvu?rts1i|sauusi", annotation_lower, perl = TRUE)
    mcr <- base::grepl("\\bmcr[abc]\\b", annotation_lower, perl = TRUE)
    if (isTRUE(flags$typeiv_gmrsd[[1]])) {
      typeiv_label <- "GmrSD domains"; typeiv_display <- "GmrSD"
    } else if (resiii_pld) {
      typeiv_label <- "ResIII+PLD"; typeiv_display <- "Res+PLD"
    } else if (isTRUE(flags$typeiv_mrr[[1]])) {
      typeiv_label <- "Mrr PF04471"; typeiv_display <- "Mrr"
    } else if (mspji) {
      typeiv_label <- "MspJI/PvuRts1I"; typeiv_display <- "MspJI"
    } else if (mcr) {
      typeiv_label <- "McrBC profile"; typeiv_display <- "McrBC"
    }
    if (!base::is.na(typeiv_label)) {
      supported[["TypeIV-profile"]] <- TRUE
      derived_label[["TypeIV-profile"]] <- typeiv_display
    } else if (isTRUE(flags$typeiv_specific[[1]])) {
      hint[["TypeIV-profile"]] <- TRUE
      derived_label[["TypeIV-profile"]] <- "Type IV profile"
    }
  }

  m_complete <- (dna_mtase && (amino_pair || c5_pair)) || (mmei_architecture && mmei_sequence)
  typeii_nuclease <- base::any(supported[c("PD-ExK", "HNH", "GIY-YIG", "PLD-HKD")])
  mmei_nuclease <- mmei_architecture && present[["MmeI-PDExK"]]
  sequence_available <- base::any(base::vapply(gene_result, function(info) {
    !base::is.null(info$prot_len) && !base::is.na(info$prot_len) && info$prot_len > 0L
  }, logical(1)))
  architecture_supported <- base::any(supported[c(
    "MmeI-architecture", "HsdS-2TRD", "TypeIV-profile"
  )])
  verified <- if (base::is.na(role) || !base::nzchar(role) ||
                  (!sequence_available && !architecture_supported)) {
    NA
  } else if (role == "M") {
    m_complete
  } else if (role == "S") {
    if (identical(family, "Type I")) supported[["HsdS-2TRD"]] else NA
  } else if (role == "RM") {
    m_complete && (mmei_nuclease || typeii_nuclease) && mmei_architecture
  } else if (role == "R" && identical(family, "Type I")) {
    hsd_motor && supported[["HsdR-PD"]]
  } else if (role == "R" && identical(family, "Type III")) {
    typeiii_motor && supported[["ResIII-PD"]]
  } else if (role == "R" && identical(family, "Type II")) {
    typeii_nuclease
  } else if (role == "R" && identical(family, "Type IV")) {
    supported[["TypeIV-profile"]]
  } else {
    NA
  }

  supported <- supported & applicable
  hint <- hint & applicable & !supported
  list(
    applicable = applicable,
    present = present,
    supported = supported,
    hint = hint,
    derived_label = derived_label,
    verified = verified,
    mtase_architecture = if (mmei_sequence) "MmeI_X-I-IV" else if (c5_extended) "C5_PCQ-ENV-QxRxR" else if (c5_pair) "C5_PCQ-ENV" else if (amino_pair) "amino_I-IV" else NA_character_,
    typeiv_architecture = typeiv_label
  )
}

# ====================================================================
# Co-occurrence scoring — confidence level per gene
# ====================================================================
.dnmb_rebasefinder_motif_score <- function(gene_result, role, prot_len,
                                            family = NA_character_,
                                            motif_defs = .dnmb_rebasefinder_motif_definitions()) {
  motif_names <- names(motif_defs)
  score <- 0L
  has <- character(0)
  for (j in seq_along(motif_names)) {
    mn <- motif_names[j]
    info <- gene_result[[j]]
    if (is.null(info) || info$n_hits == 0L) next
    if (!.dnmb_rebasefinder_motif_role_match(motif_defs[[mn]], role, family)) next
    # Check positional constraint for best hit
    in_range <- .dnmb_rebasefinder_pos_in_range(info$pos, prot_len, motif_defs[[mn]]$pos_range)
    w <- motif_defs[[mn]]$weight
    score <- score + w + if (in_range) 2L else 0L
    has <- c(has, mn)
  }
  # Co-occurrence bonuses
  sam_hit <- "SAM" %in% has
  if (sam_hit && "Amino-IV" %in% has) score <- score + 3L      # amino-MTase pair
  if (all(c("N5C-PC", "N5C-ENV") %in% has)) score <- score + 3L
  if (all(c("MmeI-X", "MmeI-I", "Amino-IV") %in% has)) score <- score + 4L
  hsd_r_motifs <- c("P-loop", "HsdR-WB", "HsdR-MIII")
  if (all(hsd_r_motifs %in% has)) score <- score + 4L            # Type I HsdR helicase
  else if (sum(hsd_r_motifs %in% has) >= 2L) score <- score + 2L
  if (all(c("ResIII-WA", "ResIII-WB", "ResIII-MIII", "ResIII-PD") %in% has)) score <- score + 4L
  else if (sum(c("ResIII-WA", "ResIII-WB", "ResIII-MIII", "ResIII-PD") %in% has) >= 2L) score <- score + 2L
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
  has_family <- "REBASEfinder_family_id" %in% names(tbl)

  results <- lapply(seq_len(nrow(tbl)), function(i) {
    seq <- .dnmb_rebasefinder_normalize_protein(tbl$translation[i])
    prot_len <- if (is.na(seq)) 0L else nchar(seq)
    if (is.na(seq)) {
      return(stats::setNames(lapply(motif_names, function(mn) list(
        status = NA, hits = list(), n_hits = 0L, prot_len = 0L,
        n_in_range = 0L, pos = NA_integer_, match = NA_character_, match_len = NA_integer_
      )), motif_names))
    }
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_
    family <- if (has_family) as.character(tbl$REBASEfinder_family_id[i]) else NA_character_

    stats::setNames(lapply(motif_names, function(mn) {
      pattern <- motif_defs[[mn]]$pattern
      if (is.null(pattern) || is.na(pattern) || !nzchar(pattern)) {
        return(list(status = "absent", hits = list(), n_hits = 0L,
                    prot_len = prot_len, n_in_range = 0L,
                    pos = NA_integer_, match = NA_character_, match_len = NA_integer_))
      }
      m <- gregexpr(pattern, seq, perl = TRUE)[[1]]
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
        role_match <- .dnmb_rebasefinder_motif_role_match(motif_defs[[mn]], role, family)
        in_pos <- length(range_hits) > 0
        status <- if (role_match && in_pos) "present"
                  else if (role_match) "present~"
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
    }), motif_names)
  })
  stats::setNames(results, tbl$locus_tag)
}

# Simplified scan for motif_verified flag
.dnmb_rebasefinder_scan_motifs <- function(tbl) {
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  if (is.null(detailed)) return(NULL)
  motif_names <- names(.dnmb_rebasefinder_motif_definitions())
  do.call(rbind, lapply(names(detailed), function(lt) {
    row <- vapply(motif_names, function(mn) {
      status <- detailed[[lt]][[which(motif_names == mn)]]$status
      if (is.null(status) || is.na(status)) NA_character_ else as.character(status)
    }, character(1))
    as.data.frame(as.list(row), stringsAsFactors = FALSE, check.names = FALSE)
  })) -> df
  df$locus_tag <- names(detailed)
  df
}

# Methylation type inference from catalytic motifs
# C5-PC motif -> N5C; amino-MTase motif needs REBASE metadata or a
# recognition-sequence heuristic and is not forced to N6A.
.dnmb_rebasefinder_infer_methylation_type <- function(tbl, cache_root = NULL) {
  # Strategy (in priority order):
  # 1. REBASE bairoch metadata — authoritative, per-enzyme lookup
  # 2. Protein motif detection — heuristic fallback
  # 3. Recognition sequence heuristic — last resort
  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  has_hit <- "REBASEfinder_hit_label" %in% names(tbl)
  has_rec <- "REBASEfinder_rec_seq" %in% names(tbl)

  # Pre-load bairoch lookup (cached after first download)
  bairoch <- tryCatch(.dnmb_rebasefinder_download_bairoch(cache_root = cache_root), error = function(e) NULL)

  # Motif detection (fallback)
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  annot <- .dnmb_rebasefinder_methylation_annotations(tbl, cache_root = cache_root)

  vapply(seq_len(nrow(tbl)), function(i) {
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_
    if (is.na(role) || !role %in% c("M", "S", "RM")) return(NA_character_)
    annot$meth_type[[i]]
  }, character(1))
}

# Motif-verified flag — uses co-occurrence scoring
# "present" = correct role + in positional range
# "present~" = correct role but outside expected position
# "present*" / "present*~" = wrong role
.dnmb_rebasefinder_motif_verified <- function(tbl) {
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  if (is.null(detailed)) return(rep(NA, nrow(tbl)))

  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  has_family <- "REBASEfinder_family_id" %in% names(tbl)
  annotations <- .dnmb_rebasefinder_gene_annotation_text(tbl)

  vapply(seq_len(nrow(tbl)), function(i) {
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_
    family <- if (has_family) as.character(tbl$REBASEfinder_family_id[i]) else NA_character_
    assessment <- .dnmb_rebasefinder_motif_assessment(
      detailed[[i]], role = role, family = family, annotation = annotations[[i]]
    )
    assessment$verified
  }, logical(1))
}


# ====================================================================
# Panel D: Motif Verification Tile — with residue + position
# ====================================================================
.dnmb_rebasefinder_core_catalytic_pairs <- function() {
  pairs <- base::data.frame(
    pair_group = c(
      "amino_mtase", "c5_mtase", "mmei_mtase", "pld_nuclease",
      "hsdr_nuclease_motor", "resiii_nuclease_motor"
    ),
    motif_a = c("SAM", "N5C-PC", "MmeI-I", "PLD-HKD", "HsdR-PD", "ResIII-PD"),
    motif_b = c("Amino-IV", "N5C-ENV", "Amino-IV", "PLD-HKD", "P-loop", "ResIII-WA"),
    stringsAsFactors = FALSE
  )
  pairs$pair_key <- base::vapply(base::seq_len(base::nrow(pairs)), function(i) {
    base::paste(base::sort(c(pairs$motif_a[[i]], pairs$motif_b[[i]])), collapse = "|")
  }, character(1))
  pairs
}

.dnmb_rebasefinder_pair_key <- function(motif_a, motif_b) {
  motif_a <- base::as.character(motif_a)
  motif_b <- base::as.character(motif_b)
  base::mapply(function(a, b) {
    if (base::is.na(a) || base::is.na(b) || !base::nzchar(a) || !base::nzchar(b)) {
      return(NA_character_)
    }
    base::paste(base::sort(c(a, b)), collapse = "|")
  }, motif_a, motif_b, USE.NAMES = FALSE)
}

.dnmb_rebasefinder_core_pair_for_locus <- function(locus,
                                                    tbl,
                                                    long,
                                                    specs) {
  locus <- .dnmb_module_clean_annotation_key(locus)
  tbl_locus <- .dnmb_module_clean_annotation_key(tbl$locus_tag)
  idx <- base::match(locus, tbl_locus)
  if (base::is.na(idx)) return(specs[0, , drop = FALSE])

  value_at <- function(candidates) {
    col <- base::intersect(candidates, base::names(tbl))
    if (!base::length(col)) return(NA_character_)
    base::as.character(tbl[[col[[1]]]][[idx]])
  }
  family <- value_at(c("REBASEfinder_family_id", "family_id"))
  role <- value_at(c("REBASEfinder_enzyme_role", "enzyme_role"))
  canonical_family <- .dnmb_rebasefinder_canonical_structure_family(family)
  locus_long <- long[.dnmb_module_clean_annotation_key(long$locus_tag) == locus, , drop = FALSE]
  motif_status <- stats::setNames(base::as.character(locus_long$status), base::as.character(locus_long$motif))
  status_for <- function(motif) {
    if (!motif %in% base::names(motif_status)) return(NA_character_)
    value <- motif_status[motif]
    if (!base::length(value)) NA_character_ else base::as.character(value[[1]])
  }
  relevant <- function(motifs) {
    status <- base::vapply(motifs, status_for, character(1))
    base::all(!base::is.na(status) & status != "na")
  }
  detected <- function(motifs) {
    status <- base::vapply(motifs, status_for, character(1))
    base::any(status %in% c("supported", "hint"))
  }

  groups <- character()
  if (!base::is.na(role) && role %in% c("M", "RM")) {
    mmei_signal <- relevant(c("MmeI-I", "Amino-IV")) &&
      detected(c("MmeI-I", "MmeI-architecture"))
    c5_signal <- relevant(c("N5C-PC", "N5C-ENV")) &&
      detected(c("N5C-PC", "N5C-ENV"))
    if (mmei_signal) {
      groups <- c(groups, "mmei_mtase")
    } else if (c5_signal) {
      groups <- c(groups, "c5_mtase")
    } else if (relevant(c("SAM", "Amino-IV"))) {
      groups <- c(groups, "amino_mtase")
    }
  }
  if (!base::is.na(role) && role %in% c("R", "RM")) {
    if (base::identical(canonical_family, "Type I")) {
      groups <- c(groups, "hsdr_nuclease_motor")
    } else if (base::identical(canonical_family, "Type III")) {
      groups <- c(groups, "resiii_nuclease_motor")
    } else if (base::identical(canonical_family, "Type II") &&
               relevant(c("PLD-HKD")) && detected(c("PLD-HKD"))) {
      groups <- c(groups, "pld_nuclease")
    }
  }
  specs[specs$pair_group %in% base::unique(groups), , drop = FALSE]
}

.dnmb_rebasefinder_catalytic_pair_boxes <- function(tbl,
                                                     long,
                                                     display_info,
                                                     display_lvls,
                                                     motif_levels,
                                                     geometry_pairs = NULL,
                                                     geometry_summary = NULL,
                                                     contact_pairs = NULL) {
  empty <- base::data.frame(
    locus_tag = character(), pair_group = character(), motif_a = character(),
    motif_b = character(), pair_state = character(), xmin = numeric(),
    xmax = numeric(), ymin = numeric(), ymax = numeric(), stringsAsFactors = FALSE
  )
  specs <- .dnmb_rebasefinder_core_catalytic_pairs()
  geometry_pairs <- base::as.data.frame(geometry_pairs, stringsAsFactors = FALSE)
  geometry_summary <- base::as.data.frame(geometry_summary, stringsAsFactors = FALSE)
  contact_pairs <- base::as.data.frame(contact_pairs, stringsAsFactors = FALSE)

  observed <- specs[0, c("pair_group", "motif_a", "motif_b", "pair_key"), drop = FALSE]
  observed$locus_tag <- character()
  observed <- observed[, c("locus_tag", "pair_group", "motif_a", "motif_b", "pair_key"), drop = FALSE]
  core_rows <- function(x) {
    if (!base::nrow(x) || !base::all(c("locus_tag", "motif_a", "motif_b") %in% base::names(x))) {
      return(observed[0, , drop = FALSE])
    }
    x$pair_key <- .dnmb_rebasefinder_pair_key(x$motif_a, x$motif_b)
    idx <- base::match(x$pair_key, specs$pair_key)
    keep <- !base::is.na(idx)
    if (!base::any(keep)) return(observed[0, , drop = FALSE])
    base::data.frame(
      locus_tag = .dnmb_module_clean_annotation_key(x$locus_tag[keep]),
      pair_group = specs$pair_group[idx[keep]],
      motif_a = specs$motif_a[idx[keep]],
      motif_b = specs$motif_b[idx[keep]],
      pair_key = specs$pair_key[idx[keep]],
      stringsAsFactors = FALSE
    )
  }
  observed <- base::rbind(core_rows(geometry_pairs), core_rows(contact_pairs))

  summary_loci <- if (base::nrow(geometry_summary) && "locus_tag" %in% base::names(geometry_summary)) {
    .dnmb_module_clean_annotation_key(geometry_summary$locus_tag)
  } else {
    character()
  }
  homology_col <- base::intersect(
    c("REBASEfinder_homology_geometry_status", "homology_geometry_status"),
    base::names(tbl)
  )
  homology_loci <- if (base::length(homology_col)) {
    status <- base::as.character(tbl[[homology_col[[1]]]])
    .dnmb_module_clean_annotation_key(tbl$locus_tag[
      !base::is.na(status) & base::nzchar(status) &
        !status %in% c("not_run", "not_applicable", "no_candidate")
    ])
  } else {
    character()
  }
  fallback_loci <- base::unique(c(summary_loci, homology_loci))
  fallback <- if (base::length(fallback_loci)) {
    base::do.call(base::rbind, base::lapply(fallback_loci, function(locus) {
      inferred <- .dnmb_rebasefinder_core_pair_for_locus(locus, tbl, long, specs)
      if (!base::nrow(inferred)) return(NULL)
      base::data.frame(
        locus_tag = locus,
        pair_group = inferred$pair_group,
        motif_a = inferred$motif_a,
        motif_b = inferred$motif_b,
        pair_key = inferred$pair_key,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    NULL
  }
  candidates <- base::rbind(observed, fallback)
  if (base::is.null(candidates) || !base::nrow(candidates)) return(empty)
  candidates <- candidates[
    !base::is.na(candidates$locus_tag) & base::nzchar(candidates$locus_tag),
    , drop = FALSE
  ]
  candidates <- candidates[!base::duplicated(candidates[c("locus_tag", "pair_group")]), , drop = FALSE]

  geometry_status <- function(locus, pair_key) {
    if (!base::nrow(geometry_pairs) ||
        !base::all(c("locus_tag", "motif_a", "motif_b") %in% base::names(geometry_pairs))) {
      return(character())
    }
    key <- .dnmb_rebasefinder_pair_key(geometry_pairs$motif_a, geometry_pairs$motif_b)
    rows <- .dnmb_module_clean_annotation_key(geometry_pairs$locus_tag) == locus & key == pair_key
    col <- base::intersect(c("combined_status", "geometry_status"), base::names(geometry_pairs))
    if (!base::any(rows) || !base::length(col)) character() else base::as.character(geometry_pairs[[col[[1]]]][rows])
  }
  summary_status <- function(locus) {
    if (!base::nrow(geometry_summary) ||
        !base::all(c("locus_tag", "structural_adjacency_status") %in% base::names(geometry_summary))) {
      return(character())
    }
    idx <- base::match(locus, .dnmb_module_clean_annotation_key(geometry_summary$locus_tag))
    if (base::is.na(idx)) character() else base::as.character(geometry_summary$structural_adjacency_status[[idx]])
  }
  contact_status <- function(locus, pair_key) {
    if (!base::nrow(contact_pairs) ||
        !base::all(c("locus_tag", "motif_a", "motif_b", "contact_status") %in% base::names(contact_pairs))) {
      return(character())
    }
    key <- .dnmb_rebasefinder_pair_key(contact_pairs$motif_a, contact_pairs$motif_b)
    rows <- .dnmb_module_clean_annotation_key(contact_pairs$locus_tag) == locus & key == pair_key
    base::as.character(contact_pairs$contact_status[rows])
  }
  homology_status <- function(locus) {
    if (!base::length(homology_col)) return(character())
    idx <- base::match(locus, .dnmb_module_clean_annotation_key(tbl$locus_tag))
    if (base::is.na(idx)) character() else base::as.character(tbl[[homology_col[[1]]]][[idx]])
  }
  supported <- c("homology_model_supported", "3d_fold_supported", "3d_geometry_supported")
  neutral <- c("", "no_structure", "not_applicable", "unavailable", "structure_missing")
  candidates$pair_state <- base::vapply(base::seq_len(base::nrow(candidates)), function(i) {
    locus <- candidates$locus_tag[[i]]
    pair_key <- candidates$pair_key[[i]]
    pair_geometry <- geometry_status(locus, pair_key)
    locus_geometry <- summary_status(locus)
    homology <- homology_status(locus)
    contacts <- contact_status(locus, pair_key)
    explicit_negative <- base::any(!base::is.na(c(locus_geometry, pair_geometry, homology)) &
      !c(locus_geometry, pair_geometry, homology) %in% c(neutral, supported))
    if (!explicit_negative &&
        base::any(c(locus_geometry, pair_geometry) %in% supported)) {
      return("verified")
    }
    if (explicit_negative) return("not_verified")
    if (base::any(contacts %in% c("contact_strong", "contact_supportive"))) {
      return("contact_supported")
    }
    "not_verified"
  }, character(1))

  x_a <- base::match(candidates$motif_a, motif_levels)
  x_b <- base::match(candidates$motif_b, motif_levels)
  motif_keep <- !base::is.na(x_a) & !base::is.na(x_b)
  if (!base::any(motif_keep)) return(empty)
  candidates <- candidates[motif_keep, , drop = FALSE]
  candidates <- base::do.call(base::rbind, base::lapply(
    base::split(base::seq_len(base::nrow(candidates)), candidates$locus_tag),
    function(idx) {
      rows <- candidates[idx, , drop = FALSE]
      states <- base::as.character(rows$pair_state)
      state <- if (base::any(states == "not_verified")) {
        "not_verified"
      } else if (base::any(states == "verified")) {
        "verified"
      } else {
        "contact_supported"
      }
      base::data.frame(
        locus_tag = rows$locus_tag[[1]],
        pair_group = base::paste(base::unique(rows$pair_group), collapse = ","),
        motif_a = base::paste(base::unique(rows$motif_a), collapse = ","),
        motif_b = base::paste(base::unique(rows$motif_b), collapse = ","),
        pair_key = base::paste(base::unique(rows$pair_key), collapse = ","),
        pair_state = state,
        stringsAsFactors = FALSE
      )
    }
  ))
  display_idx <- base::match(candidates$locus_tag, .dnmb_module_clean_annotation_key(display_info$locus_tag))
  display <- base::as.character(display_info$display[display_idx])
  y <- base::match(display, display_lvls)
  keep <- !base::is.na(y)
  if (!base::any(keep)) return(empty)
  candidates <- candidates[keep, , drop = FALSE]
  y <- y[keep]
  candidates$xmin <- 0.53
  candidates$xmax <- base::length(motif_levels) + 0.47
  candidates$ymin <- y - 0.45
  candidates$ymax <- y + 0.45
  candidates$pair_state <- base::factor(
    candidates$pair_state,
    levels = c("verified", "contact_supported", "not_verified")
  )
  candidates[, base::names(empty), drop = FALSE]
}

.dnmb_plot_rebasefinder_motif_verification <- function(tbl,
                                                       display_info,
                                                       geometry_pairs = NULL,
                                                       geometry_summary = NULL,
                                                       contact_pairs = NULL,
                                                       contact_summary = NULL) {
  detailed <- .dnmb_rebasefinder_scan_motifs_detailed(tbl)
  if (is.null(detailed)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No protein sequences", size = 3.5, color = "grey50"))
  }

  motif_defs <- .dnmb_rebasefinder_motif_definitions()
  motif_names <- names(motif_defs)
  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  has_family <- "REBASEfinder_family_id" %in% names(tbl)
  annotations <- .dnmb_rebasefinder_gene_annotation_text(tbl)

  # Build long table: only show motifs relevant to each gene's type and role.
  has_hit <- "REBASEfinder_hit_label" %in% names(tbl)
  long <- do.call(rbind, lapply(seq_along(detailed), function(i) {
    lt <- names(detailed)[i]
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_
    family <- if (has_family) as.character(tbl$REBASEfinder_family_id[i]) else NA_character_
    # Fallback: infer role from hit name prefix
    if ((is.na(role) || !nzchar(role)) && has_hit) {
      hn <- as.character(tbl$REBASEfinder_hit_label[i])
      inferred <- .dnmb_rebasefinder_role_from_hit(hn)
      if (!is.na(inferred)) role <- inferred
    }
    assessment <- .dnmb_rebasefinder_motif_assessment(
      detailed[[i]], role = role, family = family, annotation = annotations[[i]]
    )

    do.call(rbind, lapply(seq_along(motif_names), function(j) {
      info <- detailed[[i]][[j]]
      def <- motif_defs[[motif_names[j]]]
      expected <- .dnmb_rebasefinder_motif_role_label(def)
      expected_primary <- .dnmb_rebasefinder_motif_primary_role(def)

      is_relevant <- isTRUE(assessment$applicable[[motif_names[j]]])
      status <- if (!is_relevant) "na" else if (isTRUE(assessment$supported[[motif_names[j]]])) {
        "supported"
      } else if (isTRUE(assessment$hint[[motif_names[j]]])) {
        "hint"
      } else {
        "missing"
      }

		      motif_tile_label <- function(h, max_residues) {
		        residue <- h$match
		        if (nchar(residue) > max_residues) {
		          residue <- substr(residue, 1, max_residues)
		        }
		        start <- suppressWarnings(as.integer(h$pos))
		        end <- suppressWarnings(as.integer(h$end))
		        if (is.na(end) && !is.na(start)) end <- start + nchar(h$match) - 1L
		        pos_label <- if (!is.na(start) && !is.na(end)) {
		          paste0(start, "-", end)
		        } else {
		          "NA-NA"
		        }
		        paste0(residue, "\n", pos_label)
		      }
		      if (status %in% c("supported", "hint") && info$n_hits > 0) {
		        if (identical(expected_primary, "R")) {
		          plot_hits <- Filter(function(h) isTRUE(h$in_range), info$hits)
		          if (!length(plot_hits)) plot_hits <- info$hits
		          hit_labels <- vapply(plot_hits, function(h) {
		            motif_tile_label(h, max_residues = 4L)
			          }, character(1))
	        if (length(hit_labels) > 1L) {
	            tile_label <- paste0(hit_labels[[1]], " +", length(hit_labels) - 1L)
	          } else {
	            tile_label <- hit_labels
	          }
		        } else {
		          hit_labels <- vapply(info$hits, function(h) {
		            motif_tile_label(h, max_residues = 6L)
			          }, character(1))
	          if (length(hit_labels) > 1) {
	            tile_label <- paste0(hit_labels[[1]], " +", length(hit_labels) - 1)
	          } else {
	            tile_label <- hit_labels
	          }
	        }
	        n_found <- info$n_hits
	      } else if (status %in% c("supported", "hint") &&
	                 !is.na(assessment$derived_label[[motif_names[j]]])) {
	        tile_label <- gsub(" ", "\n", assessment$derived_label[[motif_names[j]]], fixed = TRUE)
	        n_found <- 1L
      } else {
        tile_label <- ""
        n_found <- 0L
      }

      data.frame(
        locus_tag = lt,
        motif = motif_names[j],
        status = status,
        expected_role = expected,
        expected_primary_role = expected_primary,
        tile_label = tile_label,
        n_found = n_found,
        stringsAsFactors = FALSE
      )
    }))
  }))

	  display_lvls <- .dnmb_rebasefinder_display_levels(display_info)
	  long <- merge(long, display_info[, c("locus_tag", "display", "operon")], by = "locus_tag", sort = FALSE)
	  long$display <- factor(as.character(long$display), levels = display_lvls)
		  motif_short <- c(
		    "SAM" = "SAM",
			    "Amino-IV" = "MTase\nIV",
		    "N5C-PC" = "m5C\nPCQ",
		    "N5C-ENV" = "m5C\nENV",
		    "N5C-QxRxR" = "m5C\nQxRxR",
		    "MmeI-X" = "MmeI\nX",
		    "MmeI-I" = "MmeI\nSAM-I",
		    "MmeI-PDExK" = "MmeI\nREase",
		    "MmeI-architecture" = "MmeI\narch",
		    "HsdS-2TRD" = "HsdS\n2xTRD",
		    "HsdR-PD" = "HsdR\nPD",
		    "PD-ExK" = "PD\nExK",
	    "HNH" = "HNH",
		    "GIY-YIG" = "GIY",
		    "PLD-HKD" = "PLD\n2xHKD",
			    "ResIII-WA" = "Res\nWA",
			    "ResIII-WB" = "Res\nWB",
			    "ResIII-MIII" = "Res\nMIII",
			    "ResIII-PD" = "Res\nPD",
			    "P-loop" = "P\nloop",
		    "HsdR-WB" = "HsdR\nWB",
		    "HsdR-MIII" = "HsdR\nMIII",
			    "TypeIV-profile" = "Type IV\nprofile"
			  )
			  motif_levels <- motif_names[motif_names %in% unique(long$motif[long$status != "na"])]
			  long <- long[long$motif %in% motif_levels, , drop = FALSE]
			  long$motif_label <- ifelse(long$motif %in% names(motif_short), motif_short[long$motif], long$motif)
		  motif_label_levels <- ifelse(motif_levels %in% names(motif_short), motif_short[motif_levels], motif_levels)
		  long$motif_label <- factor(long$motif_label, levels = motif_label_levels)
		  pair_boxes <- .dnmb_rebasefinder_catalytic_pair_boxes(
		    tbl = tbl,
		    long = long,
		    display_info = display_info,
		    display_lvls = display_lvls,
		    motif_levels = motif_levels,
		    geometry_pairs = geometry_pairs,
		    geometry_summary = geometry_summary,
		    contact_pairs = contact_pairs
		  )

  # Status palette: green/rose/teal identify coherent role-specific evidence;
  # amber is an isolated sequence/profile hint and never verifies a protein.
	  long$fill_status <- ifelse(
	    long$status == "supported" & long$expected_primary_role == "R", "supported_R",
	    ifelse(long$status == "supported" & long$expected_primary_role == "S", "supported_S",
	    ifelse(long$status == "supported" & long$expected_primary_role == "RM", "supported_RM",
	    ifelse(long$status == "supported", "supported_M",
	    ifelse(long$status == "hint", "hint", long$status))))
	  )

  status_pal <- c(
    "supported_M" = "#2E7D32",
    "supported_R" = "#C2185B",
    "supported_S" = "#00897B",
	    "supported_RM" = "#6A4C93",
	    "hint" = "#F9A825",
	    "missing" = "#BDBDBD",   # grey
	    "na"      = "#F5F5F5"    # very light grey
  )

  # Operon separators
  operon_bounds <- .dnmb_rebasefinder_operon_separators(
    merge(display_info, data.frame(locus_tag = tbl$locus_tag), by = "locus_tag", sort = FALSE)
  )

			  long$label_color <- ifelse(grepl("^supported_", long$fill_status), "found", "other")

	  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$motif_label, y = .data$display)) +
	    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill_status), color = "white", linewidth = 0.6)
	  if (base::nrow(pair_boxes)) {
	    p <- p + ggplot2::geom_rect(
	      data = pair_boxes,
	      ggplot2::aes(
	        xmin = .data$xmin, xmax = .data$xmax,
	        ymin = .data$ymin, ymax = .data$ymax,
	        color = .data$pair_state, linetype = .data$pair_state
	      ),
	      fill = NA,
	      linewidth = 0.525,
	      inherit.aes = FALSE,
	      show.legend = TRUE
	    )
	  }

  if (nrow(operon_bounds)) {
    p <- p + ggplot2::geom_hline(yintercept = operon_bounds$y, color = "grey75", linewidth = 0.3, linetype = "dotted")
  }

	  found_labels <- long[long$label_color == "found", , drop = FALSE]
	  other_labels <- long[long$label_color != "found", , drop = FALSE]
	  p <- p +
	    ggplot2::geom_text(
	      data = found_labels,
	      ggplot2::aes(label = .data$tile_label),
	      color = "white", size = 2.05, fontface = "bold", lineheight = 0.82
	    ) +
	    ggplot2::geom_text(
	      data = other_labels,
	      ggplot2::aes(label = .data$tile_label),
	      color = "grey15", size = 2.05, fontface = "bold", lineheight = 0.82
	    ) +
		    ggplot2::scale_fill_manual(
		      values = status_pal, name = "Evidence",
	      labels = c("supported_M" = "MTase set", "supported_R" = "Restriction evidence",
	                 "supported_S" = "Specificity architecture", "supported_RM" = "RM fusion architecture",
	                 "hint" = "Sequence/profile hint", "missing" = "Not detected", "na" = "N/A"),
			      breaks = c("supported_M", "supported_R", "supported_S", "supported_RM",
			                 "hint", "missing", "na")
	    ) +
	    ggplot2::scale_color_manual(
	      values = c(
	        verified = "#2166AC",
	        contact_supported = "#67A9CF",
	        not_verified = "#B2182B"
	      ),
	      name = "Catalytic pair",
	      breaks = c("verified", "contact_supported", "not_verified"),
	      labels = c(
	        verified = "3D verified",
	        contact_supported = "Contact-supported",
	        not_verified = "Not 3D-verified"
	      )
	    ) +
	    ggplot2::scale_linetype_manual(
	      values = c(verified = "solid", contact_supported = "dashed", not_verified = "solid"),
	      name = "Catalytic pair",
	      breaks = c("verified", "contact_supported", "not_verified"),
	      labels = c(
	        verified = "3D verified",
	        contact_supported = "Contact-supported",
	        not_verified = "Not 3D-verified"
	      )
	    ) +
	    ggplot2::scale_y_discrete(limits = display_lvls, drop = FALSE) +
	    ggplot2::coord_fixed(ratio = 1, clip = "off") +
	    ggplot2::guides(
	      fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 2),
	      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 1),
	      linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 1)
	    ) +
			    ggplot2::labs(title = "  Functional evidence: motifs and catalytic-pair proximity", x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
	      axis.text.x = ggplot2::element_text(size = 7.5, angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.82),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
	      legend.position = "bottom",
	      legend.justification = "left",
	      legend.box = "vertical",
	      legend.box.just = "left",
	      legend.key.size = ggplot2::unit(0.3, "cm"),
      legend.title = ggplot2::element_text(size = 7),
	      legend.text = ggplot2::element_text(size = 6.5),
	      legend.spacing.x = ggplot2::unit(0.1, "cm"),
	      legend.spacing.y = ggplot2::unit(0.02, "cm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
	      plot.margin = ggplot2::margin(4, 4, 4, 2)
	    )
	  base::attr(p, "dnmb_catalytic_pair_boxes") <- pair_boxes
	  p
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
  allow_remote <- base::tolower(base::Sys.getenv("DNMB_ALLOW_REMOTE_ANNOTATION", unset = "")) %in%
    c("1", "true", "yes", "on")
  if (!allow_remote) return(NULL)

  results <- tryCatch({
    do.call(rbind, lapply(seq_len(nrow(tbl)), function(i) {
      pid <- as.character(tbl$protein_id[i])
      lt <- tbl$locus_tag[i]
      valid_refseq <- !is.na(pid) && base::grepl(
        "^(WP|NP|XP|YP|ZP|AP)_?[0-9]+([.][0-9]+)?$", pid, ignore.case = TRUE
      )
      if (!valid_refseq) return(NULL)

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
	    "SAM"      = "#2E7D32",  # dark green (M - SAM binding)
	    "Amino-IV" = "#9CCC65",  # amino-MTase catalytic motif IV
    "N5C-PC"   = "#1565C0",  # N5C C5-MTase catalytic PC motif
    "N5C-ENV"  = "#00838F",
    "N5C-QxRxR" = "#26A69A",
    "MmeI-X" = "#558B2F",
    "MmeI-I" = "#33691E",
    "MmeI-PDExK" = "#B71C1C",
    "HsdR-PD" = "#B71C1C",
    "PD-ExK"  = "#C62828",   # dark red (R)
    "HNH"     = "#EF5350",   # light red (R)
    "GIY-YIG" = "#F48FB1",   # pink (R)
	    "P-loop"  = "#E91E63",   # rose (R)
	    "HsdR-WB" = "#AD1457",
	    "HsdR-MIII" = "#C2185B",
	    "PLD-HKD" = "#880E4F",
	    "ResIII-WA" = "#6A1B9A",
	    "ResIII-WB" = "#8E24AA",
	    "ResIII-MIII" = "#9C27B0",
	    "ResIII-PD" = "#AB47BC"
  )

  # Build protein backbones
  prot_lengths <- vapply(tbl$translation, function(sequence) {
    normalized <- .dnmb_rebasefinder_normalize_protein(sequence)
    if (is.na(normalized)) NA_integer_ else nchar(normalized)
  }, integer(1))
  display_lvls <- .dnmb_rebasefinder_display_levels(display_info)

  backbone <- merge(
    data.frame(locus_tag = tbl$locus_tag, prot_len = prot_lengths, stringsAsFactors = FALSE),
    display_info[, c("locus_tag", "display")],
    by = "locus_tag",
    sort = FALSE
  )
  backbone$display <- factor(as.character(backbone$display), levels = display_lvls)

  has_role <- "REBASEfinder_enzyme_role" %in% names(tbl)
  has_family <- "REBASEfinder_family_id" %in% names(tbl)
  annotations <- .dnmb_rebasefinder_gene_annotation_text(tbl)

  # Motif markers: only motifs matching each gene's R-M type and subunit role.
  markers <- do.call(rbind, lapply(seq_along(detailed), function(i) {
    lt <- names(detailed)[i]
    role <- if (has_role) as.character(tbl$REBASEfinder_enzyme_role[i]) else NA_character_
    family <- if (has_family) as.character(tbl$REBASEfinder_family_id[i]) else NA_character_
    assessment <- .dnmb_rebasefinder_motif_assessment(
      detailed[[i]], role = role, family = family, annotation = annotations[[i]]
    )
    do.call(rbind, lapply(seq_along(motif_names), function(j) {
      info <- detailed[[i]][[j]]
      def <- motif_defs[[motif_names[j]]]
      if (info$n_hits <= 0) return(NULL)
      if (!isTRUE(assessment$supported[[motif_names[j]]])) return(NULL)
      plot_hits <- Filter(function(h) isTRUE(h$in_range), info$hits)
      if (!length(plot_hits)) plot_hits <- info$hits
      do.call(rbind, lapply(plot_hits, function(h) {
        data.frame(locus_tag = lt, pos = h$pos,
                   motif_short = motif_names[j],
                   stringsAsFactors = FALSE)
      }))
    }))
  }))

  if (is.null(markers) || !nrow(markers)) {
    markers <- data.frame(
      locus_tag = character(), pos = integer(), motif_short = character(),
      display = factor(levels = display_lvls), prot_len = integer(), frac = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    markers <- merge(markers, display_info[, c("locus_tag", "display")], by = "locus_tag", sort = FALSE)
    markers <- merge(markers, backbone[, c("locus_tag", "prot_len")], by = "locus_tag", sort = FALSE)
    markers$display <- factor(as.character(markers$display), levels = display_lvls)
    markers$frac <- markers$pos / markers$prot_len
  }
  operon_blocks <- .dnmb_rebasefinder_operon_blocks(display_info)

  p <- ggplot2::ggplot()
  if (nrow(operon_blocks)) {
    p <- p + ggplot2::geom_rect(
      data = operon_blocks,
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = .data$ymin, ymax = .data$ymax),
      fill = operon_blocks$fill, alpha = 0.55, color = NA,
      inherit.aes = FALSE
    )
  }
  p <- p +
    # Protein backbone: geom_tile with exact height matching domain rects.
    ggplot2::geom_tile(
      data = backbone,
      ggplot2::aes(x = 0.5, y = .data$display, width = 1, height = 0.30),
      fill = "grey80", color = NA
    )

  # UniProt annotation overlay
  if (!is.null(uniprot_doms) && nrow(uniprot_doms)) {
    dom_plot <- merge(uniprot_doms, backbone[, c("locus_tag", "display", "prot_len")], by = "locus_tag", sort = FALSE)
    dom_plot$display <- factor(as.character(dom_plot$display), levels = display_lvls)
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

  if (nrow(markers)) {
    p <- p + ggplot2::geom_point(
      data = markers,
      ggplot2::aes(x = .data$frac, y = .data$display, color = .data$motif_short),
      size = 2.5, alpha = 0.9, shape = 18
    )
  }

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
    ggplot2::scale_y_discrete(limits = display_lvls, drop = FALSE) +
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

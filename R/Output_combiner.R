#' Result DNMB table
#'
#' @param gb_table A data frame containing the data.
#' @param InterPro_search Logical, whether to include p-values for means.
#' @param InterPro_site Logical, whether to include p-values for medians.
#' @param codon Logical, whether to use chi-square test for categorical variables.
#' @return A data frame containing the baseline table.
#' @export
#'

DNMB_table <- function(
    genbank_table = TRUE,
    EggNOG_table = FALSE,
    InterProScan_table = FALSE,
    InterProScan_site = FALSE,
    codon_usage = TRUE,
    tRNA_anticodon = TRUE,
    tRNA_distribution = TRUE,
    RBS_table = TRUE,
    CRISPR_by_spacer = TRUE,
    save_dir = NULL
) {
  # Set the save directory to the current working directory if not provided
  if (is.null(save_dir)) {
    save_dir <- getwd()
  }

  # Check for genbank_table
  if (isTRUE(genbank_table)) {
    if (exists("genbank_table", envir = .GlobalEnv)) {
      genbank_table <- get("genbank_table", envir = .GlobalEnv)
    } else {
      stop("genbank_table is set to TRUE but not found in the environment.")
    }
  }

  # Check for EggNOG_table
  if (isTRUE(EggNOG_table)) {
    if (exists("EggNOG_table", envir = .GlobalEnv)) {
      InterPro_table <- get("EggNOG_table", envir = .GlobalEnv)
    } else {
      stop("EggNOG_table is set to TRUE but not found in the environment.")
    }
  }

  # Check for codon_usage
  if (isTRUE(codon_usage)) {
    if (exists("codon_usage", envir = .GlobalEnv)) {
      codon_usage <- get("codon_usage", envir = .GlobalEnv)
    } else {
      stop("codon_usage is set to TRUE but not found in the environment.")
    }
  }

  # Check for tRNA_anticodon
  if (isTRUE(tRNA_anticodon)) {
    if (exists("tRNA_anticodon", envir = .GlobalEnv)) {
      tRNA_anticodon <- get("tRNA_anticodon", envir = .GlobalEnv)
    } else {
      stop("tRNA_anticodon is set to TRUE but not found in the environment.")
    }
  }

  # Check for tRNA_distribution
  if (isTRUE(tRNA_distribution)) {
    if (exists("tRNA_distribution", envir = .GlobalEnv)) {
      tRNA_distribution <- get("tRNA_distribution", envir = .GlobalEnv)
    } else {
      stop("tRNA_distribution is set to TRUE but not found in the environment.")
    }
  }

  # Check for RBS_table
  if (isTRUE(RBS_table)) {
    if (exists("RBS_table", envir = .GlobalEnv)) {
      RBS_table <- get("RBS_table", envir = .GlobalEnv)
    } else {
      stop("RBS_table is set to TRUE but not found in the environment.")
    }
  }

  # Check for CRISPR_by_spacer
  if (isTRUE(CRISPR_by_spacer)) {
    if (exists("CRISPR_by_spacer", envir = .GlobalEnv)) {
      CRISPR_by_spacer <- get("CRISPR_by_spacer", envir = .GlobalEnv)
    } else {
      warning("CRISPR_by_spacer is set to TRUE but not found in the environment. Creating an empty sheet.")
      CRISPR_by_spacer <- NULL
    }
  }

  # InterProScan_site is actually driven by .GlobalEnv, not by the
  # caller flag. The parameter default is FALSE for historical reasons,
  # but that should NOT suppress the sheet when a fresh run has parsed
  # `.tsv.sites` and populated `InterProScan_site` in the global env.
  # Unless the caller passed an actual data.frame, always look up the
  # global env and use it; otherwise leave as NULL so the sheet is
  # skipped. The downstream guard requires is.data.frame() anyway.
  if (!is.data.frame(InterProScan_site)) {
    if (exists("InterProScan_site", envir = .GlobalEnv, inherits = FALSE)) {
      InterProScan_site <- get("InterProScan_site", envir = .GlobalEnv, inherits = FALSE)
      message("InterProScan_site found in the environment.")
    } else {
      InterProScan_site <- NULL
    }
  }

  genbank_table <- dnmb_prepare_genbank_table_for_output(genbank_table)
  module_detail_table <- dnmb_build_module_details_table(genbank_table)
  genbank_table <- dnmb_order_genbank_table_for_output(genbank_table)

  # Create a new workbook
  wb <- openxlsx::createWorkbook()

  # Add worksheets
  openxlsx::addWorksheet(wb, "1.GenBank_table")
  openxlsx::addWorksheet(wb, "2.Codon_usage")
  openxlsx::addWorksheet(wb, "3.tRNA_anticodon")
  openxlsx::addWorksheet(wb, "4.tRNA_distribution")
  openxlsx::addWorksheet(wb, "5.RBS")
  openxlsx::addWorksheet(wb, "6.CRISPR_table")

  no_wrap_style <- openxlsx::createStyle(wrapText = FALSE)

  # Write data to worksheets
  openxlsx::writeData(wb, "1.GenBank_table", genbank_table, startRow = 1, startCol = 1)
  dnmb_apply_group_header_styles(wb, "1.GenBank_table", genbank_table, startRow = 1, startCol = 1)
  dnmb_write_hyperlink_columns(wb, "1.GenBank_table", genbank_table, startRow = 1, startCol = 1)
  dnmb_apply_numeric_display_styles(wb, "1.GenBank_table", genbank_table, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "2.Codon_usage", codon_usage, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "3.tRNA_anticodon", tRNA_anticodon, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "4.tRNA_distribution", tRNA_distribution, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "5.RBS", RBS_table, startRow = 1, startCol = 1)

  # Handle CRISPR table (empty or not)
  if (is.null(CRISPR_by_spacer)) {
    openxlsx::writeData(wb, "6.CRISPR_table", data.frame(), startRow = 1, startCol = 1)
  } else {
    openxlsx::writeData(wb, "6.CRISPR_table", CRISPR_by_spacer, startRow = 1, startCol = 1)
  }

  openxlsx::addStyle(wb, sheet = "1.GenBank_table", style = no_wrap_style, rows = 1:nrow(genbank_table) + 1, cols = 1:ncol(genbank_table), gridExpand = TRUE)

  # Adding the worksheet only if InterProScan_site holds an actual table
  if (is.data.frame(InterProScan_site) && nrow(InterProScan_site) > 0) {
    openxlsx::addWorksheet(wb, "7.InterPro_site")
    openxlsx::writeData(wb, "7.InterPro_site", InterProScan_site, startRow = 1, startCol = 1)
    message("InterProScan_site sheet added (", nrow(InterProScan_site), " rows).")
  }

  if (is.data.frame(module_detail_table) && nrow(module_detail_table)) {
    openxlsx::addWorksheet(wb, "8.Module_details")
    openxlsx::writeData(wb, "8.Module_details", module_detail_table, startRow = 1, startCol = 1)
    dnmb_apply_group_header_styles(wb, "8.Module_details", module_detail_table, startRow = 1, startCol = 1)
    dnmb_apply_module_category_styles(wb, "8.Module_details", module_detail_table, startRow = 1, startCol = 1)
    dnmb_write_hyperlink_columns(wb, "8.Module_details", module_detail_table, startRow = 1, startCol = 1)
    dnmb_apply_numeric_display_styles(wb, "8.Module_details", module_detail_table, startRow = 1, startCol = 1)
    openxlsx::addStyle(
      wb,
      sheet = "8.Module_details",
      style = no_wrap_style,
      rows = seq_len(nrow(module_detail_table)) + 1L,
      cols = seq_len(ncol(module_detail_table)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  # Prophage summary sheet
  prophage_summary <- tryCatch({
    gbff_candidates <- list.files(getwd(), pattern = "\\.(gbff|gbk|gb)$", full.names = TRUE)
    gbff_for_gc <- if (length(gbff_candidates)) gbff_candidates[[1]] else NULL
    .dnmb_prophage_build_summary_table(genbank_table, gbff_path = gbff_for_gc)
  }, error = function(e) data.frame())
  if (is.data.frame(prophage_summary) && nrow(prophage_summary)) {
    openxlsx::addWorksheet(wb, "9.Prophage_summary")
    openxlsx::writeData(wb, "9.Prophage_summary", prophage_summary, startRow = 1, startCol = 1)
    prophage_header_style <- openxlsx::createStyle(
      fgFill = "#FDE8E8", fontColour = "#1F2937",
      textDecoration = "bold", halign = "center",
      border = "Bottom", borderColour = "#D1D5DB"
    )
    openxlsx::addStyle(wb, "9.Prophage_summary", prophage_header_style,
                       rows = 1, cols = seq_len(ncol(prophage_summary)), gridExpand = TRUE)
    openxlsx::setColWidths(wb, "9.Prophage_summary", cols = seq_len(ncol(prophage_summary)), widths = "auto")
    message("Prophage_summary sheet added (", nrow(prophage_summary), " regions).")
  }

  # IS element census sheet
  is_census <- tryCatch({
    census_path <- file.path(getwd(), "dnmb_module_iselement", "iselement_census.tsv")
    if (file.exists(census_path)) utils::read.delim(census_path, check.names = FALSE) else data.frame()
  }, error = function(e) data.frame())
  if (is.data.frame(is_census) && nrow(is_census)) {
    openxlsx::addWorksheet(wb, "10.IS_census")
    openxlsx::writeData(wb, "10.IS_census", is_census, startRow = 1, startCol = 1)
    is_census_header <- openxlsx::createStyle(
      fgFill = "#E8F5E9", fontColour = "#1F2937",
      textDecoration = "bold", halign = "center",
      border = "Bottom", borderColour = "#D1D5DB"
    )
    openxlsx::addStyle(wb, "10.IS_census", is_census_header,
                       rows = 1, cols = seq_len(ncol(is_census)), gridExpand = TRUE)
    openxlsx::setColWidths(wb, "10.IS_census", cols = seq_len(ncol(is_census)), widths = "auto")
    message("IS_census sheet added (", nrow(is_census), " families).")
  }

  # Landing pad sheet
  is_landing_pads <- tryCatch({
    lp_path <- file.path(getwd(), "dnmb_module_iselement", "iselement_landing_pads.tsv")
    if (file.exists(lp_path)) utils::read.delim(lp_path, check.names = FALSE) else data.frame()
  }, error = function(e) data.frame())
  if (is.data.frame(is_landing_pads) && nrow(is_landing_pads)) {
    openxlsx::addWorksheet(wb, "11.Landing_pads")
    openxlsx::writeData(wb, "11.Landing_pads", is_landing_pads, startRow = 1, startCol = 1)
    lp_header <- openxlsx::createStyle(
      fgFill = "#E3F2FD", fontColour = "#1F2937",
      textDecoration = "bold", halign = "center",
      border = "Bottom", borderColour = "#D1D5DB"
    )
    openxlsx::addStyle(wb, "11.Landing_pads", lp_header,
                       rows = 1, cols = seq_len(ncol(is_landing_pads)), gridExpand = TRUE)
    openxlsx::setColWidths(wb, "11.Landing_pads", cols = seq_len(ncol(is_landing_pads)), widths = "auto")
    # Conditional formatting for landing pad score
    score_col <- match("landing_pad_score", names(is_landing_pads))
    if (!is.na(score_col)) {
      openxlsx::conditionalFormatting(
        wb, "11.Landing_pads",
        cols = score_col,
        rows = 2:(nrow(is_landing_pads) + 1),
        type = "colourScale",
        style = c("#EF5350", "#FFC107", "#66BB6A"),
        rule = c(0, 0.55, 1)
      )
    }
    message("Landing_pads sheet added (", nrow(is_landing_pads), " pads).")
  }

  # Construct the file name
  gb_dir <- getwd()
  gb_files <- list.files(gb_dir, pattern = "\\.gbk$|\\.gb$|\\.gbff$", full.names = FALSE)
  gb_files <- gb_files[!grepl("^~\\$", gb_files)]

  file_name <- paste0(.dnmb_beg2char(gb_files, "."), "_total.xlsx")
  save_path <- file.path(save_dir, file_name)

  # Save the workbook
  dnmb_apply_workbook_sheet_order(wb)
  openxlsx::saveWorkbook(wb, file = save_path, overwrite = TRUE)

  message(paste("Workbook saved as", save_path))
}

.dnmb_backfill_iselement_columns <- function(genbank_table, output_dir = getwd()) {
  if (!is.data.frame(genbank_table) || !nrow(genbank_table) || "ISelement_element_id" %in% names(genbank_table)) {
    return(genbank_table)
  }

  module_dir <- if (exists(".dnmb_iselement_plot_module_dir", mode = "function", inherits = TRUE)) {
    .dnmb_iselement_plot_module_dir(output_dir)
  } else {
    candidates <- c(
      file.path(output_dir, "dnmb_module_iselement"),
      file.path(output_dir, "dnmb_module_ISelement")
    )
    existing <- candidates[dir.exists(candidates)]
    if (length(existing)) existing[[1]] else candidates[[1]]
  }

  elements_path <- file.path(module_dir, "iselement_elements.tsv")
  if (!file.exists(elements_path)) {
    return(genbank_table)
  }

  elements <- tryCatch(utils::read.delim(elements_path, check.names = FALSE), error = function(e) data.frame())
  if (!is.data.frame(elements) || !nrow(elements)) {
    return(genbank_table)
  }

  hits <- tryCatch(.dnmb_iselement_normalize_hits(elements), error = function(e) data.frame())
  if (!is.data.frame(hits) || !nrow(hits)) {
    return(genbank_table)
  }

  module_table <- tryCatch(.dnmb_iselement_output_table(genbank_table, hits), error = function(e) data.frame())
  if (!is.data.frame(module_table) || !nrow(module_table)) {
    return(genbank_table)
  }

  if ("essentiality_class" %in% names(genbank_table) && "locus_tag" %in% names(module_table)) {
    ess_map <- stats::setNames(genbank_table$essentiality_class, genbank_table$locus_tag)
    module_table$iselement_essentiality_class <- ess_map[module_table$locus_tag]
  }

  landing_pads_path <- file.path(module_dir, "iselement_landing_pads.tsv")
  if (file.exists(landing_pads_path)) {
    landing_pads <- tryCatch(utils::read.delim(landing_pads_path, check.names = FALSE), error = function(e) data.frame())
    if (is.data.frame(landing_pads) && nrow(landing_pads)) {
      module_table$iselement_nearest_landing_pad_rank <- .dnmb_nearest_landing_pad_rank(module_table, landing_pads, genbank_table)
    }
  }

  run <- structure(
    list(database = "ISelement", output_table = module_table, hits = hits),
    class = "dnmb_module_run"
  )
  append_module_results(genbank_table, list(ISelement = run))
}

dnmb_prepare_genbank_table_for_output <- function(genbank_table) {
  if (!is.data.frame(genbank_table) || !nrow(genbank_table)) {
    return(genbank_table)
  }

  out <- as.data.frame(genbank_table, stringsAsFactors = FALSE, check.names = FALSE)
  drop_cols <- intersect("PAZy_verified", names(out))
  if (length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  out <- .dnmb_backfill_iselement_columns(out, output_dir = getwd())

  character_cols <- names(out)[vapply(out, is.character, logical(1))]
  if (length(character_cols)) {
    for (column_name in character_cols) {
      out[[column_name]] <- dnmb_sanitize_output_text(out[[column_name]])
    }
  }

  if ("note" %in% names(out)) {
    out$note <- dnmb_strip_automated_note_phrase(out$note)
  }

  if ("CLEAN_best_hit_label" %in% names(out)) {
    if (!"CLEAN_expasy_link" %in% names(out)) {
      out$CLEAN_expasy_link <- NA_character_
    }
    if (!"CLEAN_brenda_link" %in% names(out)) {
      out$CLEAN_brenda_link <- NA_character_
    }
    missing_expasy <- is.na(out$CLEAN_expasy_link) | !nzchar(out$CLEAN_expasy_link)
    missing_brenda <- is.na(out$CLEAN_brenda_link) | !nzchar(out$CLEAN_brenda_link)
    out$CLEAN_expasy_link[missing_expasy] <- .dnmb_clean_expasy_url(out$CLEAN_best_hit_label[missing_expasy])
    out$CLEAN_brenda_link[missing_brenda] <- .dnmb_clean_brenda_url(out$CLEAN_best_hit_label[missing_brenda])
  }

  out
}

dnmb_order_genbank_table_for_output <- function(genbank_table) {
  if (!is.data.frame(genbank_table) || !nrow(genbank_table)) {
    return(genbank_table)
  }

  out <- as.data.frame(genbank_table, stringsAsFactors = FALSE, check.names = FALSE)
  separator_cols <- paste0(dnmb_supported_module_prefixes(), "_SECTION")
  present_separators <- intersect(separator_cols, names(out))
  if (length(present_separators)) {
    out[present_separators] <- NULL
  }

  module_prefixes <- dnmb_module_prefix_order(out)
  module_cols <- unlist(lapply(module_prefixes, function(prefix) grep(paste0("^", prefix, "_"), names(out), value = TRUE)), use.names = FALSE)
  base_cols <- setdiff(names(out), module_cols)
  if ("locus_tag" %in% base_cols) {
    base_cols <- c("locus_tag", setdiff(base_cols, "locus_tag"))
  }
  if ("old_locus_tag" %in% base_cols) {
    base_cols <- c("locus_tag", "old_locus_tag", setdiff(base_cols, c("locus_tag", "old_locus_tag")))
  }

  ordered_cols <- base_cols
  for (prefix in module_prefixes) {
    block_cols <- dnmb_module_block_columns(out, prefix = prefix)
    if (!length(block_cols)) {
      next
    }
    separator_name <- paste0(prefix, "_SECTION")
    out[[separator_name]] <- ""
    ordered_cols <- c(ordered_cols, separator_name, block_cols)
  }

  ordered_cols <- unique(c(ordered_cols, names(out)))
  out[, ordered_cols, drop = FALSE]
}

dnmb_supported_module_prefixes <- function() {
  c("CLEAN", "GapMindAA", "GapMindCarbon", "GapMind", "DefenseFinder", "dbAPIS", "AcrFinder", "Promotech", "mRNAcal", "PADLOC", "DefensePredictor", "ISelement", "Prophage", "dbCAN", "MEROPS", "PAZy")
}

dnmb_module_prefix_order <- function(genbank_table) {
  if (!is.data.frame(genbank_table) || !ncol(genbank_table)) {
    return(character())
  }
  prefixes <- dnmb_supported_module_prefixes()
  # Prefer more specific names first so GapMindAA is not collapsed into GapMind.
  prefixes <- prefixes[order(nchar(prefixes), decreasing = TRUE)]
  found <- character()
  for (column_name in names(genbank_table)) {
    matches <- prefixes[startsWith(column_name, paste0(prefixes, "_"))]
    if (length(matches)) {
      found <- c(found, matches[[1]])
    }
  }
  unique(found)
}

dnmb_module_prefix_kind <- function(prefix) {
  prefix <- as.character(prefix)[1]
  if (prefix %in% c("GapMindAA", "GapMindCarbon", "GapMind")) {
    return("GapMind")
  }
  if (prefix %in% c("DefenseFinder")) {
    return("DefenseFinder")
  }
  prefix
}

dnmb_module_block_columns <- function(genbank_table, prefix) {
  prefix <- as.character(prefix)[1]
  prefix_pattern <- paste0("^", gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", prefix), "_")
  block_cols <- names(genbank_table)[grepl(prefix_pattern, names(genbank_table))]
  if (!length(block_cols)) {
    return(character())
  }

  desired_suffixes <- switch(
    dnmb_module_prefix_kind(prefix),
    CLEAN = c("best_hit_label", "best_distance", "pvalue_best_pvalue", "predictions", "expasy_link", "brenda_link"),
    dbCAN = c(
      "dbcan_hit",
      "family_id",
      "support",
      "profile_id",
      "evalue",
      "coverage",
      "profile_length",
      "gene_length",
      "profile_start",
      "profile_end",
      "gene_start",
      "gene_end",
      "dbcan_cgc_id",
      "dbcan_cgc_gene_type",
      "dbcan_cgc_protein_family",
      "dbcan_pul_id",
      "dbcan_pul_substrate",
      "dbcan_pul_bitscore",
      "dbcan_pul_signature_pairs",
      "dbcan_sub_substrate",
      "dbcan_sub_score"
    ),
    MEROPS = c(
      "family_id",
      "hit_label",
      "subject_accession",
      "evalue",
      "bitscore",
      "pident",
      "qcov",
      "scov",
      "alignment_length",
      "query_length",
      "subject_length",
      "query_start",
      "query_end",
      "subject_start",
      "subject_end"
    ),
    PAZy = c(
      "family_id",
      "hit_label",
      "substrate_label",
      "pazy_id",
      "accession_list",
      "organism_name",
      "evalue",
      "bitscore",
      "pident",
      "qcov",
      "scov",
      "alignment_length",
      "query_length",
      "subject_length"
    ),
    GapMind = c(
      "pathway_id",
      "pathway_desc",
      "step_id",
      "step_score",
      "confidence",
      "on_best_path",
      "reference_id",
      "candidate_type",
      "support",
      "locus_tag2",
      "curated_ids",
      "curated_desc",
      "blast_bits",
      "identity",
      "blast_coverage",
      "blast_score",
      "hmm_bits",
      "hmm_id",
      "hmm_coverage",
      "hmm_score",
      "hmm_desc",
      "other_ids",
      "other_bits",
      "other_identity",
      "other_coverage"
    ),
    DefenseFinder = c(
      "system_type",
      "system_subtype",
      "system_activity",
      "system_id",
      "gene_name",
      "hit_status",
      "system_wholeness",
      "system_score",
      "genes_count",
      "system_begin",
      "system_end",
      "profiles_in_system",
      "hit_i_eval",
      "hit_score",
      "hit_profile_cov",
      "hit_seq_cov",
      "support"
    ),
    dbAPIS = c(
      "family_id",
      "hit_label",
      "defense_type",
      "clan_id",
      "clan_defense_type",
      "representative_protein",
      "description",
      "i_evalue",
      "hit_score",
      "hmm_coverage",
      "query_coverage",
      "support",
      "evidence_mode",
      "substrate_label",
      "profile_length",
      "gene_length",
      "hmm_from",
      "hmm_to",
      "ali_from",
      "ali_to"
    ),
    AcrFinder = c(
      "family_id",
      "hit_label",
      "enzyme_role",
      "evidence_mode",
      "substrate_label",
      "classification",
      "pident",
      "support"
    ),
    Promotech = c(
      "family_id",
      "hit_label",
      "promoter_id",
      "promoter_score",
      "promoter_start",
      "promoter_end",
      "promoter_strand",
      "distance_to_gene",
      "promoter_sequence",
      "support"
    ),
    mRNAcal = c(
      "family_id",
      "hit_label",
      "tir_score",
      "tir_score_band",
      "tir_score_percentile",
      "rbs_motif",
      "rbs_spacer",
      "rbs_score",
      "anti_sd_sequence",
      "duplex_energy",
      "duplex_score",
      "duplex_motif",
      "start_codon",
      "start_codon_score",
      "early_k_score",
      "lysine_codon_count_2_8",
      "aaa_count_2_8",
      "aag_count_2_8",
      "early_coding_at_fraction",
      "ncs45_sequence",
      "ncs45_at_fraction",
      "ncs45_lysine_codon_count",
      "ncs45_aaa_count",
      "ncs45_aag_count",
      "early_aaa_run",
      "early_poly_a_penalty",
      "skik_like",
      "internal_sd_count",
      "internal_sd_motifs",
      "internal_sd_min_position",
      "internal_sd_penalty",
      "upstream20_sequence",
      "upstream20_at_fraction",
      "upstream_au_score",
      "fold_mfe",
      "fold_mfe_per_nt",
      "fold_score",
      "rbs_unpaired_fraction",
      "start_unpaired_fraction",
      "rbs_plfold_unpaired_probability",
      "start_plfold_unpaired_probability",
      "tir_plfold_unpaired_probability",
      "standby_plfold_unpaired_probability",
      "downstream_plfold_unpaired_probability",
      "plfold_accessibility_score",
      "accessibility_score",
      "accessibility_method",
      "support"
    ),
    PADLOC = c(
      "system",
      "system_number",
      "protein_name",
      "target_description",
      "hmm_accession",
      "hmm_name",
      "full_seq_evalue",
      "domain_ievalue",
      "target_coverage",
      "hmm_coverage",
      "relative_position",
      "seqid",
      "support"
    ),
    DefensePredictor = c(
      "score_band",
      "mean_log_odds",
      "sd_log_odds",
      "min_log_odds",
      "max_log_odds",
      "protein_context_id",
      "product_accession",
      "genomic_accession",
      "name",
      "symbol",
      "feature_class",
      "feature_interval_length",
      "product_length",
      "assembly",
      "assembly_unit",
      "seq_type",
      "chromosome",
      "attributes",
      "support"
    ),
    ISelement = c(
      "element_id",
      "element_type",
      "feature_type",
      "confidence",
      "confidence_score",
      "family_system",
      "evidence_label",
      "supporting_gene_count",
      "supporting_genes",
      "support"
    ),
    Prophage = c(
      "prophage_id",
      "prophage_rank",
      "prophage_my_status",
      "prophage_pp",
      "prophage_start",
      "prophage_end",
      "support"
    ),
    character()
  )

  desired <- paste0(prefix, "_", desired_suffixes)
  c(intersect(desired, block_cols), setdiff(block_cols, desired))
}

dnmb_build_module_details_table <- function(genbank_table) {
  if (!is.data.frame(genbank_table) || !nrow(genbank_table)) {
    return(data.frame())
  }

  out <- as.data.frame(genbank_table, stringsAsFactors = FALSE, check.names = FALSE)
  base_cols <- intersect(
    c("locus_tag", "old_locus_tag", "gene", "product", "protein_id", "contig", "start", "end", "direction"),
    names(out)
  )
  detail_parts <- list(
    dnmb_module_details_clean(out, base_cols),
    dnmb_module_details_gapmind_full(out, base_cols, version = "aa", label = "GapMindAA"),
    dnmb_module_details_gapmind_full(out, base_cols, version = "carbon", label = "GapMindCarbon"),
    dnmb_module_details_gapmind(out, base_cols, prefix = "GapMind", label = "GapMind"),
    dnmb_module_details_defensefinder(out, base_cols),
    dnmb_module_details_dbapis(out, base_cols),
    dnmb_module_details_acrfinder(out, base_cols),
    dnmb_module_details_promotech(out, base_cols),
    dnmb_module_details_mrnacal(out, base_cols),
    dnmb_module_details_padloc(out, base_cols),
    dnmb_module_details_defensepredictor(out, base_cols),
    dnmb_module_details_iselement(out, base_cols),
    dnmb_module_details_prophage(out, base_cols),
    dnmb_module_details_dbcan(out, base_cols),
    dnmb_module_details_merops(out, base_cols),
    dnmb_module_details_pazy(out, base_cols),
    dnmb_module_details_rebasefinder(out, base_cols)
  )
  detail_parts <- Filter(function(x) is.data.frame(x) && nrow(x), detail_parts)
  if (!length(detail_parts)) {
    return(data.frame())
  }

  detail <- dplyr::bind_rows(detail_parts)
  dynamic_levels <- dnmb_module_detail_order(out)
  detail_levels <- unique(as.character(detail$module_category))
  detail_levels <- detail_levels[!is.na(detail_levels) & nzchar(detail_levels)]
  dynamic_levels <- c(dynamic_levels, setdiff(detail_levels, dynamic_levels))
  detail$module_category <- factor(detail$module_category, levels = dynamic_levels)
  locus_rank <- suppressWarnings(as.numeric(sub("^.*?(\\d+)$", "\\1", detail$locus_tag)))
  detail <- detail[order(detail$module_category, locus_rank, detail$locus_tag), , drop = FALSE]
  detail$module_category <- as.character(detail$module_category)
  ordered_cols <- c(
    "module_category",
    "locus_tag",
    "old_locus_tag",
    "gene",
    "product",
    "protein_id",
    "location",
    "primary_call",
    "reference_id",
    "significance",
    "score_summary",
    "alignment_summary",
    "context_summary",
    "resource_links"
  )
  detail <- detail[, intersect(ordered_cols, names(detail)), drop = FALSE]
  rownames(detail) <- NULL
  detail
}

dnmb_module_detail_order <- function(genbank_table) {
  prefix_order <- dnmb_module_prefix_order(genbank_table)
  label_map <- c(
    CLEAN = "CLEAN",
    GapMindAA = "GapMindAA",
    GapMindCarbon = "GapMindCarbon",
    GapMind = "GapMind",
    DefenseFinder = "DefenseFinder",
    dbAPIS = "dbAPIS",
    AcrFinder = "AcrFinder",
    Promotech = "Promotech",
    mRNAcal = "mRNAcal",
    PADLOC = "PADLOC",
    DefensePredictor = "DefensePredictor",
    REBASEfinder = "REBASE",
    ISelement = "ISelement",
    Prophage = "Prophage",
    dbCAN = "dbCAN",
    MEROPS = "MEROPS",
    PAZy = "PAZy"
  )
  out <- unname(label_map[prefix_order])
  out[!is.na(out) & nzchar(out)]
}

dnmb_module_detail_base <- function(genbank_table, base_cols, keep) {
  base <- genbank_table[keep, base_cols, drop = FALSE]
  for (column_name in setdiff(c("locus_tag", "old_locus_tag", "gene", "product", "protein_id", "contig", "start", "end", "direction"), names(base))) {
    base[[column_name]] <- NA
  }
  base <- base[, c("locus_tag", "old_locus_tag", "gene", "product", "protein_id", "contig", "start", "end", "direction"), drop = FALSE]
  base$location <- dnmb_module_location_summary(base$contig, base$start, base$end, base$direction)
  base[, c("locus_tag", "old_locus_tag", "gene", "product", "protein_id", "location"), drop = FALSE]
}

dnmb_module_detail_template <- function(base) {
  cbind(
    data.frame(
      module_category = character(nrow(base)),
      primary_call = character(nrow(base)),
      reference_id = character(nrow(base)),
      significance = character(nrow(base)),
      score_summary = character(nrow(base)),
      alignment_summary = character(nrow(base)),
      context_summary = character(nrow(base)),
      resource_links = character(nrow(base)),
      stringsAsFactors = FALSE
    ),
    base,
    stringsAsFactors = FALSE
  )
}

dnmb_module_details_clean <- function(genbank_table, base_cols) {
  if (!all(c("CLEAN_best_hit_label", "CLEAN_best_distance", "CLEAN_pvalue_best_pvalue") %in% names(genbank_table))) {
    return(data.frame())
  }
  keep <- !is.na(genbank_table$CLEAN_best_hit_label) & nzchar(genbank_table$CLEAN_best_hit_label)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "CLEAN"
  out$primary_call <- as.character(genbank_table$CLEAN_best_hit_label[keep])
  out$significance <- dnmb_detail_key_value("pvalue", genbank_table$CLEAN_pvalue_best_pvalue[keep], digits = 3)
  out$score_summary <- dnmb_detail_key_value("distance", genbank_table$CLEAN_best_distance[keep], digits = 3)
  out$context_summary <- if ("CLEAN_predictions" %in% names(genbank_table)) as.character(genbank_table$CLEAN_predictions[keep]) else NA_character_
  out$resource_links <- dnmb_detail_links(
    if ("CLEAN_expasy_link" %in% names(genbank_table)) as.character(genbank_table$CLEAN_expasy_link[keep]) else NA_character_,
    if ("CLEAN_brenda_link" %in% names(genbank_table)) as.character(genbank_table$CLEAN_brenda_link[keep]) else NA_character_
  )
  out
}

dnmb_module_details_dbapis <- function(genbank_table, base_cols) {
  label_col <- if ("dbAPIS_family_id" %in% names(genbank_table)) genbank_table$dbAPIS_family_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "dbAPIS"
  out$primary_call <- if ("dbAPIS_hit_label" %in% names(genbank_table)) {
    label <- as.character(genbank_table$dbAPIS_hit_label[keep])
    fallback <- as.character(genbank_table$dbAPIS_family_id[keep])
    ifelse(!is.na(label) & nzchar(label), label, fallback)
  } else {
    as.character(genbank_table$dbAPIS_family_id[keep])
  }
  out$reference_id <- dnmb_detail_join(
    if ("dbAPIS_family_id" %in% names(genbank_table)) paste0("family=", as.character(genbank_table$dbAPIS_family_id[keep])) else NA_character_,
    if ("dbAPIS_clan_id" %in% names(genbank_table)) paste0("clan=", as.character(genbank_table$dbAPIS_clan_id[keep])) else NA_character_,
    if ("dbAPIS_representative_protein" %in% names(genbank_table)) paste0("rep=", as.character(genbank_table$dbAPIS_representative_protein[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_key_value("evalue", if ("dbAPIS_i_evalue" %in% names(genbank_table)) genbank_table$dbAPIS_i_evalue[keep] else NA, digits = 3)
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("score", if ("dbAPIS_hit_score" %in% names(genbank_table)) genbank_table$dbAPIS_hit_score[keep] else NA, digits = 2),
    dnmb_detail_key_value("hmm_cov", if ("dbAPIS_hmm_coverage" %in% names(genbank_table)) genbank_table$dbAPIS_hmm_coverage[keep] else NA, digits = 3),
    dnmb_detail_key_value("query_cov", if ("dbAPIS_query_coverage" %in% names(genbank_table)) genbank_table$dbAPIS_query_coverage[keep] else NA, digits = 3)
  )
  out$alignment_summary <- dnmb_detail_join(
    dnmb_detail_span("hmm", if ("dbAPIS_hmm_from" %in% names(genbank_table)) genbank_table$dbAPIS_hmm_from[keep] else NA, if ("dbAPIS_hmm_to" %in% names(genbank_table)) genbank_table$dbAPIS_hmm_to[keep] else NA, if ("dbAPIS_profile_length" %in% names(genbank_table)) genbank_table$dbAPIS_profile_length[keep] else NA),
    dnmb_detail_span("gene", if ("dbAPIS_ali_from" %in% names(genbank_table)) genbank_table$dbAPIS_ali_from[keep] else NA, if ("dbAPIS_ali_to" %in% names(genbank_table)) genbank_table$dbAPIS_ali_to[keep] else NA, if ("dbAPIS_gene_length" %in% names(genbank_table)) genbank_table$dbAPIS_gene_length[keep] else NA)
  )
  out$context_summary <- dnmb_detail_join(
    if ("dbAPIS_clan_defense_type" %in% names(genbank_table)) paste0("target=", as.character(genbank_table$dbAPIS_clan_defense_type[keep])) else NA_character_,
    if ("dbAPIS_defense_type" %in% names(genbank_table)) paste0("type=", as.character(genbank_table$dbAPIS_defense_type[keep])) else NA_character_,
    if ("dbAPIS_description" %in% names(genbank_table)) as.character(genbank_table$dbAPIS_description[keep]) else NA_character_,
    if ("dbAPIS_support" %in% names(genbank_table)) as.character(genbank_table$dbAPIS_support[keep]) else NA_character_
  )
  out
}

dnmb_module_details_acrfinder <- function(genbank_table, base_cols) {
  label_col <- if ("AcrFinder_family_id" %in% names(genbank_table)) genbank_table$AcrFinder_family_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "AcrFinder"
  out$primary_call <- if ("AcrFinder_hit_label" %in% names(genbank_table)) {
    label <- as.character(genbank_table$AcrFinder_hit_label[keep])
    fallback <- as.character(genbank_table$AcrFinder_family_id[keep])
    ifelse(!is.na(label) & nzchar(label), label, fallback)
  } else {
    as.character(genbank_table$AcrFinder_family_id[keep])
  }
  out$reference_id <- if ("AcrFinder_family_id" %in% names(genbank_table)) paste0("acr=", as.character(genbank_table$AcrFinder_family_id[keep])) else NA_character_
  out$significance <- dnmb_detail_key_value("pident", if ("AcrFinder_pident" %in% names(genbank_table)) genbank_table$AcrFinder_pident[keep] else NA, digits = 2)
  out$score_summary <- dnmb_detail_join(
    if ("AcrFinder_evidence_mode" %in% names(genbank_table)) paste0("evidence=", as.character(genbank_table$AcrFinder_evidence_mode[keep])) else NA_character_,
    if ("AcrFinder_classification" %in% names(genbank_table)) paste0("classification=", as.character(genbank_table$AcrFinder_classification[keep])) else NA_character_
  )
  out$context_summary <- dnmb_detail_join(
    if ("AcrFinder_enzyme_role" %in% names(genbank_table)) paste0("target=", as.character(genbank_table$AcrFinder_enzyme_role[keep])) else NA_character_,
    if ("AcrFinder_substrate_label" %in% names(genbank_table)) paste0("label=", as.character(genbank_table$AcrFinder_substrate_label[keep])) else NA_character_,
    if ("AcrFinder_support" %in% names(genbank_table)) as.character(genbank_table$AcrFinder_support[keep]) else NA_character_
  )
  out
}

dnmb_module_details_promotech <- function(genbank_table, base_cols) {
  label_col <- if ("Promotech_promoter_id" %in% names(genbank_table)) genbank_table$Promotech_promoter_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "Promotech"
  out$primary_call <- if ("Promotech_hit_label" %in% names(genbank_table)) {
    label <- as.character(genbank_table$Promotech_hit_label[keep])
    fallback <- if ("Promotech_family_id" %in% names(genbank_table)) as.character(genbank_table$Promotech_family_id[keep]) else "promoter"
    ifelse(!is.na(label) & nzchar(label), label, fallback)
  } else {
    "Promotech promoter"
  }
  out$reference_id <- paste0("promoter=", as.character(genbank_table$Promotech_promoter_id[keep]))
  out$significance <- dnmb_detail_key_value("score", if ("Promotech_promoter_score" %in% names(genbank_table)) genbank_table$Promotech_promoter_score[keep] else NA, digits = 3)
  out$alignment_summary <- dnmb_detail_join(
    if (all(c("Promotech_promoter_start", "Promotech_promoter_end") %in% names(genbank_table))) {
      paste0("interval=", genbank_table$Promotech_promoter_start[keep], "-", genbank_table$Promotech_promoter_end[keep])
    } else {
      NA_character_
    },
    if ("Promotech_promoter_strand" %in% names(genbank_table)) paste0("strand=", as.character(genbank_table$Promotech_promoter_strand[keep])) else NA_character_,
    dnmb_detail_key_value("distance_to_gene", if ("Promotech_distance_to_gene" %in% names(genbank_table)) genbank_table$Promotech_distance_to_gene[keep] else NA, digits = 0)
  )
  out$context_summary <- dnmb_detail_join(
    if ("Promotech_promoter_sequence" %in% names(genbank_table)) paste0("sequence=", as.character(genbank_table$Promotech_promoter_sequence[keep])) else NA_character_,
    if ("Promotech_support" %in% names(genbank_table)) as.character(genbank_table$Promotech_support[keep]) else NA_character_
  )
  out
}

dnmb_module_details_mrnacal <- function(genbank_table, base_cols) {
  label_col <- if ("mRNAcal_tir_score" %in% names(genbank_table)) genbank_table$mRNAcal_tir_score else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(suppressWarnings(as.numeric(label_col)))
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "mRNAcal"
  out$primary_call <- if ("mRNAcal_hit_label" %in% names(genbank_table)) {
    as.character(genbank_table$mRNAcal_hit_label[keep])
  } else {
    paste0("Translation efficiency: ", as.character(genbank_table$mRNAcal_tir_score_band[keep]))
  }
  out$reference_id <- if ("mRNAcal_rbs_motif" %in% names(genbank_table)) {
    paste0("RBS=", as.character(genbank_table$mRNAcal_rbs_motif[keep]))
  } else {
    NA_character_
  }
  out$significance <- dnmb_detail_join(
    dnmb_detail_key_value("TIR_score", genbank_table$mRNAcal_tir_score[keep], digits = 2),
    dnmb_detail_key_value("pct", if ("mRNAcal_tir_score_percentile" %in% names(genbank_table)) genbank_table$mRNAcal_tir_score_percentile[keep] else NA, digits = 1)
  )
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("RBS", if ("mRNAcal_rbs_score" %in% names(genbank_table)) genbank_table$mRNAcal_rbs_score[keep] else NA, digits = 1),
    dnmb_detail_key_value("antiSD", if ("mRNAcal_duplex_score" %in% names(genbank_table)) genbank_table$mRNAcal_duplex_score[keep] else NA, digits = 1),
    dnmb_detail_key_value("access", if ("mRNAcal_accessibility_score" %in% names(genbank_table)) genbank_table$mRNAcal_accessibility_score[keep] else NA, digits = 1),
    dnmb_detail_key_value("downA", if ("mRNAcal_downstream_plfold_unpaired_probability" %in% names(genbank_table)) 100 * suppressWarnings(as.numeric(genbank_table$mRNAcal_downstream_plfold_unpaired_probability[keep])) else NA, digits = 1),
    dnmb_detail_key_value("upAU", if ("mRNAcal_upstream_au_score" %in% names(genbank_table)) genbank_table$mRNAcal_upstream_au_score[keep] else NA, digits = 1),
    dnmb_detail_key_value("earlyK", if ("mRNAcal_lysine_codon_count_2_8" %in% names(genbank_table)) genbank_table$mRNAcal_lysine_codon_count_2_8[keep] else NA, digits = 0),
    dnmb_detail_key_value("NCS45_K", if ("mRNAcal_ncs45_lysine_codon_count" %in% names(genbank_table)) genbank_table$mRNAcal_ncs45_lysine_codon_count[keep] else NA, digits = 0),
    dnmb_detail_key_value("intSD", if ("mRNAcal_internal_sd_count" %in% names(genbank_table)) genbank_table$mRNAcal_internal_sd_count[keep] else NA, digits = 0),
    dnmb_detail_key_value("MFE", if ("mRNAcal_fold_mfe" %in% names(genbank_table)) genbank_table$mRNAcal_fold_mfe[keep] else NA, digits = 2)
  )
  out$alignment_summary <- dnmb_detail_join(
    dnmb_detail_key_value("spacer", if ("mRNAcal_rbs_spacer" %in% names(genbank_table)) genbank_table$mRNAcal_rbs_spacer[keep] else NA, digits = 0),
    dnmb_detail_key_value("duplex_dG", if ("mRNAcal_duplex_energy" %in% names(genbank_table)) genbank_table$mRNAcal_duplex_energy[keep] else NA, digits = 2),
    if ("mRNAcal_start_codon" %in% names(genbank_table)) paste0("start=", as.character(genbank_table$mRNAcal_start_codon[keep])) else NA_character_,
    if ("mRNAcal_upstream20_sequence" %in% names(genbank_table)) paste0("up20=", as.character(genbank_table$mRNAcal_upstream20_sequence[keep])) else NA_character_,
    if ("mRNAcal_skik_like" %in% names(genbank_table)) paste0("SKIK_like=", as.character(genbank_table$mRNAcal_skik_like[keep])) else NA_character_,
    if ("mRNAcal_second_codon" %in% names(genbank_table)) paste0("second=", as.character(genbank_table$mRNAcal_second_codon[keep])) else NA_character_
  )
  out$context_summary <- dnmb_detail_join(
    if ("mRNAcal_accessibility_method" %in% names(genbank_table)) paste0("access_method=", as.character(genbank_table$mRNAcal_accessibility_method[keep])) else NA_character_,
    if ("mRNAcal_early_codons" %in% names(genbank_table)) paste0("early_codons=", as.character(genbank_table$mRNAcal_early_codons[keep])) else NA_character_,
    if ("mRNAcal_support" %in% names(genbank_table)) as.character(genbank_table$mRNAcal_support[keep]) else NA_character_
  )
  out
}

dnmb_module_details_dbcan <- function(genbank_table, base_cols) {
  label_col <- if ("dbCAN_dbcan_hit" %in% names(genbank_table)) genbank_table$dbCAN_dbcan_hit else if ("dbCAN_family_id" %in% names(genbank_table)) genbank_table$dbCAN_family_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "dbCAN"
  out$primary_call <- if ("dbCAN_dbcan_hit" %in% names(genbank_table)) as.character(genbank_table$dbCAN_dbcan_hit[keep]) else as.character(genbank_table$dbCAN_family_id[keep])
  out$reference_id <- dnmb_detail_join(
    if ("dbCAN_profile_id" %in% names(genbank_table)) paste0("profile=", as.character(genbank_table$dbCAN_profile_id[keep])) else NA_character_,
    if ("dbCAN_family_id" %in% names(genbank_table)) paste0("family=", as.character(genbank_table$dbCAN_family_id[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_key_value("evalue", if ("dbCAN_evalue" %in% names(genbank_table)) genbank_table$dbCAN_evalue[keep] else NA, digits = 3)
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("coverage", if ("dbCAN_coverage" %in% names(genbank_table)) genbank_table$dbCAN_coverage[keep] else NA, digits = 3),
    if ("dbCAN_support" %in% names(genbank_table)) paste0("support=", as.character(genbank_table$dbCAN_support[keep])) else NA_character_
  )
  out$alignment_summary <- dnmb_detail_join(
    dnmb_detail_span("profile", if ("dbCAN_profile_start" %in% names(genbank_table)) genbank_table$dbCAN_profile_start[keep] else NA, if ("dbCAN_profile_end" %in% names(genbank_table)) genbank_table$dbCAN_profile_end[keep] else NA, if ("dbCAN_profile_length" %in% names(genbank_table)) genbank_table$dbCAN_profile_length[keep] else NA),
    dnmb_detail_span("gene", if ("dbCAN_gene_start" %in% names(genbank_table)) genbank_table$dbCAN_gene_start[keep] else NA, if ("dbCAN_gene_end" %in% names(genbank_table)) genbank_table$dbCAN_gene_end[keep] else NA, if ("dbCAN_gene_length" %in% names(genbank_table)) genbank_table$dbCAN_gene_length[keep] else NA)
  )
  out$context_summary <- dnmb_detail_join(
    if ("dbCAN_dbcan_cgc_id" %in% names(genbank_table)) paste0("cgc=", as.character(genbank_table$dbCAN_dbcan_cgc_id[keep])) else NA_character_,
    if ("dbCAN_dbcan_cgc_gene_type" %in% names(genbank_table)) paste0("cgc_type=", as.character(genbank_table$dbCAN_dbcan_cgc_gene_type[keep])) else NA_character_,
    if ("dbCAN_dbcan_cgc_protein_family" %in% names(genbank_table)) paste0("cgc_family=", as.character(genbank_table$dbCAN_dbcan_cgc_protein_family[keep])) else NA_character_,
    if ("dbCAN_dbcan_pul_id" %in% names(genbank_table)) paste0("pul=", as.character(genbank_table$dbCAN_dbcan_pul_id[keep])) else NA_character_,
    if ("dbCAN_dbcan_pul_substrate" %in% names(genbank_table)) paste0("pul_substrate=", as.character(genbank_table$dbCAN_dbcan_pul_substrate[keep])) else NA_character_,
    dnmb_detail_key_value("pul_bitscore", if ("dbCAN_dbcan_pul_bitscore" %in% names(genbank_table)) genbank_table$dbCAN_dbcan_pul_bitscore[keep] else NA, digits = 3),
    if ("dbCAN_dbcan_pul_signature_pairs" %in% names(genbank_table)) paste0("signature=", as.character(genbank_table$dbCAN_dbcan_pul_signature_pairs[keep])) else NA_character_,
    if ("dbCAN_dbcan_sub_substrate" %in% names(genbank_table)) paste0("dbcan_sub=", as.character(genbank_table$dbCAN_dbcan_sub_substrate[keep])) else NA_character_,
    if ("dbCAN_dbcan_sub_score" %in% names(genbank_table)) paste0("dbcan_sub_score=", as.character(genbank_table$dbCAN_dbcan_sub_score[keep])) else NA_character_
  )
  out
}

dnmb_module_details_merops <- function(genbank_table, base_cols) {
  label_col <- if ("MEROPS_hit_label" %in% names(genbank_table)) genbank_table$MEROPS_hit_label else if ("MEROPS_family_id" %in% names(genbank_table)) genbank_table$MEROPS_family_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "MEROPS"
  out$primary_call <- if ("MEROPS_family_id" %in% names(genbank_table)) as.character(genbank_table$MEROPS_family_id[keep]) else as.character(genbank_table$MEROPS_hit_label[keep])
  out$reference_id <- dnmb_detail_join(
    if ("MEROPS_subject_accession" %in% names(genbank_table)) paste0("subject=", as.character(genbank_table$MEROPS_subject_accession[keep])) else NA_character_,
    if ("MEROPS_hit_label" %in% names(genbank_table)) paste0("label=", as.character(genbank_table$MEROPS_hit_label[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_key_value("evalue", if ("MEROPS_evalue" %in% names(genbank_table)) genbank_table$MEROPS_evalue[keep] else NA, digits = 3)
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("bitscore", if ("MEROPS_bitscore" %in% names(genbank_table)) genbank_table$MEROPS_bitscore[keep] else NA, digits = 3),
    dnmb_detail_key_value("pident", if ("MEROPS_pident" %in% names(genbank_table)) genbank_table$MEROPS_pident[keep] else NA, digits = 3),
    dnmb_detail_key_value("qcov", if ("MEROPS_qcov" %in% names(genbank_table)) genbank_table$MEROPS_qcov[keep] else NA, digits = 3),
    dnmb_detail_key_value("scov", if ("MEROPS_scov" %in% names(genbank_table)) genbank_table$MEROPS_scov[keep] else NA, digits = 3)
  )
  out$alignment_summary <- dnmb_detail_join(
    dnmb_detail_key_value("aln", if ("MEROPS_alignment_length" %in% names(genbank_table)) genbank_table$MEROPS_alignment_length[keep] else NA, digits = 0),
    dnmb_detail_span("query", if ("MEROPS_query_start" %in% names(genbank_table)) genbank_table$MEROPS_query_start[keep] else NA, if ("MEROPS_query_end" %in% names(genbank_table)) genbank_table$MEROPS_query_end[keep] else NA, if ("MEROPS_query_length" %in% names(genbank_table)) genbank_table$MEROPS_query_length[keep] else NA),
    dnmb_detail_span("subject", if ("MEROPS_subject_start" %in% names(genbank_table)) genbank_table$MEROPS_subject_start[keep] else NA, if ("MEROPS_subject_end" %in% names(genbank_table)) genbank_table$MEROPS_subject_end[keep] else NA, if ("MEROPS_subject_length" %in% names(genbank_table)) genbank_table$MEROPS_subject_length[keep] else NA)
  )
  out
}

dnmb_module_details_pazy <- function(genbank_table, base_cols) {
  label_col <- if ("PAZy_hit_label" %in% names(genbank_table)) genbank_table$PAZy_hit_label else if ("PAZy_family_id" %in% names(genbank_table)) genbank_table$PAZy_family_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "PAZy"
  out$primary_call <- if ("PAZy_family_id" %in% names(genbank_table)) as.character(genbank_table$PAZy_family_id[keep]) else as.character(genbank_table$PAZy_hit_label[keep])
  out$reference_id <- dnmb_detail_join(
    if ("PAZy_pazy_id" %in% names(genbank_table)) paste0("pazy=", as.character(genbank_table$PAZy_pazy_id[keep])) else NA_character_,
    if ("PAZy_accession_list" %in% names(genbank_table)) paste0("accession=", as.character(genbank_table$PAZy_accession_list[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_key_value("evalue", if ("PAZy_evalue" %in% names(genbank_table)) genbank_table$PAZy_evalue[keep] else NA, digits = 3)
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("bitscore", if ("PAZy_bitscore" %in% names(genbank_table)) genbank_table$PAZy_bitscore[keep] else NA, digits = 3),
    dnmb_detail_key_value("pident", if ("PAZy_pident" %in% names(genbank_table)) genbank_table$PAZy_pident[keep] else NA, digits = 3),
    dnmb_detail_key_value("qcov", if ("PAZy_qcov" %in% names(genbank_table)) genbank_table$PAZy_qcov[keep] else NA, digits = 3),
    dnmb_detail_key_value("scov", if ("PAZy_scov" %in% names(genbank_table)) genbank_table$PAZy_scov[keep] else NA, digits = 3)
  )
  out$alignment_summary <- dnmb_detail_join(
    dnmb_detail_key_value("aln", if ("PAZy_alignment_length" %in% names(genbank_table)) genbank_table$PAZy_alignment_length[keep] else NA, digits = 0),
    dnmb_detail_key_value("query_len", if ("PAZy_query_length" %in% names(genbank_table)) genbank_table$PAZy_query_length[keep] else NA, digits = 0),
    dnmb_detail_key_value("subject_len", if ("PAZy_subject_length" %in% names(genbank_table)) genbank_table$PAZy_subject_length[keep] else NA, digits = 0)
  )
  out$context_summary <- dnmb_detail_join(
    if ("PAZy_substrate_label" %in% names(genbank_table)) paste0("substrate=", as.character(genbank_table$PAZy_substrate_label[keep])) else NA_character_,
    if ("PAZy_organism_name" %in% names(genbank_table)) paste0("organism=", as.character(genbank_table$PAZy_organism_name[keep])) else NA_character_
  )
  out
}

dnmb_module_details_rebasefinder <- function(genbank_table, base_cols) {
  label_col <- if ("REBASEfinder_family_id" %in% names(genbank_table)) genbank_table$REBASEfinder_family_id else NULL
  if (is.null(label_col)) return(data.frame())
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) return(data.frame())
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "REBASE"
  out$primary_call <- as.character(genbank_table$REBASEfinder_family_id[keep])
  out$reference_id <- dnmb_detail_join(
    if ("REBASEfinder_hit_label" %in% names(genbank_table)) as.character(genbank_table$REBASEfinder_hit_label[keep]) else NA_character_,
    if ("REBASEfinder_rec_seq" %in% names(genbank_table)) paste0("rec=", as.character(genbank_table$REBASEfinder_rec_seq[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_join(
    dnmb_detail_key_value("evalue", if ("REBASEfinder_blast_evalue" %in% names(genbank_table)) genbank_table$REBASEfinder_blast_evalue[keep] else NA, digits = 3),
    dnmb_detail_key_value("identity", if ("REBASEfinder_blast_identity" %in% names(genbank_table)) genbank_table$REBASEfinder_blast_identity[keep] else NA, digits = 1)
  )
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("bitscore", if ("REBASEfinder_blast_bitscore" %in% names(genbank_table)) genbank_table$REBASEfinder_blast_bitscore[keep] else NA, digits = 1),
    dnmb_detail_key_value("aln_len", if ("REBASEfinder_blast_length" %in% names(genbank_table)) genbank_table$REBASEfinder_blast_length[keep] else NA, digits = 0)
  )
  out$context_summary <- dnmb_detail_join(
    if ("REBASEfinder_enzyme_role" %in% names(genbank_table)) paste0("role=", as.character(genbank_table$REBASEfinder_enzyme_role[keep])) else NA_character_,
    if ("REBASEfinder_operon_id" %in% names(genbank_table)) paste0("operon=", as.character(genbank_table$REBASEfinder_operon_id[keep])) else NA_character_
  )
  out
}

dnmb_module_details_padloc <- function(genbank_table, base_cols) {
  label_col <- if ("PADLOC_system" %in% names(genbank_table)) genbank_table$PADLOC_system else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "PADLOC"
  out$primary_call <- as.character(genbank_table$PADLOC_system[keep])
  out$reference_id <- dnmb_detail_join(
    if ("PADLOC_system_number" %in% names(genbank_table)) paste0("system=", as.character(genbank_table$PADLOC_system_number[keep])) else NA_character_,
    if ("PADLOC_hmm_accession" %in% names(genbank_table)) paste0("hmm=", as.character(genbank_table$PADLOC_hmm_accession[keep])) else NA_character_,
    if ("PADLOC_protein_name" %in% names(genbank_table)) paste0("protein=", as.character(genbank_table$PADLOC_protein_name[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_join(
    dnmb_detail_key_value("full_e", if ("PADLOC_full_seq_evalue" %in% names(genbank_table)) genbank_table$PADLOC_full_seq_evalue[keep] else NA, digits = 3),
    dnmb_detail_key_value("domain_iE", if ("PADLOC_domain_ievalue" %in% names(genbank_table)) genbank_table$PADLOC_domain_ievalue[keep] else NA, digits = 3)
  )
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("target_cov", if ("PADLOC_target_coverage" %in% names(genbank_table)) genbank_table$PADLOC_target_coverage[keep] else NA, digits = 3),
    dnmb_detail_key_value("hmm_cov", if ("PADLOC_hmm_coverage" %in% names(genbank_table)) genbank_table$PADLOC_hmm_coverage[keep] else NA, digits = 3),
    dnmb_detail_key_value("rel_pos", if ("PADLOC_relative_position" %in% names(genbank_table)) genbank_table$PADLOC_relative_position[keep] else NA, digits = 0)
  )
  out$context_summary <- dnmb_detail_join(
    if ("PADLOC_target_description" %in% names(genbank_table)) as.character(genbank_table$PADLOC_target_description[keep]) else NA_character_,
    if ("PADLOC_support" %in% names(genbank_table)) as.character(genbank_table$PADLOC_support[keep]) else NA_character_
  )
  out
}

dnmb_module_details_defensepredictor <- function(genbank_table, base_cols) {
  label_col <- if ("DefensePredictor_mean_log_odds" %in% names(genbank_table)) genbank_table$DefensePredictor_mean_log_odds else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "DefensePredictor"
  out$primary_call <- if ("DefensePredictor_score_band" %in% names(genbank_table)) as.character(genbank_table$DefensePredictor_score_band[keep]) else "DefensePredictor"
  out$reference_id <- dnmb_detail_join(
    if ("DefensePredictor_product_accession" %in% names(genbank_table)) paste0("product_accession=", as.character(genbank_table$DefensePredictor_product_accession[keep])) else NA_character_,
    if ("DefensePredictor_protein_context_id" %in% names(genbank_table)) paste0("context=", as.character(genbank_table$DefensePredictor_protein_context_id[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_join(
    dnmb_detail_key_value("mean", if ("DefensePredictor_mean_log_odds" %in% names(genbank_table)) genbank_table$DefensePredictor_mean_log_odds[keep] else NA, digits = 3),
    dnmb_detail_key_value("sd", if ("DefensePredictor_sd_log_odds" %in% names(genbank_table)) genbank_table$DefensePredictor_sd_log_odds[keep] else NA, digits = 3)
  )
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("min", if ("DefensePredictor_min_log_odds" %in% names(genbank_table)) genbank_table$DefensePredictor_min_log_odds[keep] else NA, digits = 3),
    dnmb_detail_key_value("max", if ("DefensePredictor_max_log_odds" %in% names(genbank_table)) genbank_table$DefensePredictor_max_log_odds[keep] else NA, digits = 3)
  )
  out$alignment_summary <- dnmb_detail_join(
    dnmb_detail_key_value("feature_nt", if ("DefensePredictor_feature_interval_length" %in% names(genbank_table)) genbank_table$DefensePredictor_feature_interval_length[keep] else NA, digits = 0),
    dnmb_detail_key_value("aa_len", if ("DefensePredictor_product_length" %in% names(genbank_table)) genbank_table$DefensePredictor_product_length[keep] else NA, digits = 0)
  )
  out$context_summary <- dnmb_detail_join(
    if ("DefensePredictor_name" %in% names(genbank_table)) paste0("name=", as.character(genbank_table$DefensePredictor_name[keep])) else NA_character_,
    if ("DefensePredictor_symbol" %in% names(genbank_table)) paste0("symbol=", as.character(genbank_table$DefensePredictor_symbol[keep])) else NA_character_,
    if ("DefensePredictor_support" %in% names(genbank_table)) as.character(genbank_table$DefensePredictor_support[keep]) else NA_character_
  )
  out
}

dnmb_module_details_gapmind <- function(genbank_table, base_cols, prefix = "GapMind", label = prefix) {
  label_col <- if (paste0(prefix, "_pathway_id") %in% names(genbank_table)) genbank_table[[paste0(prefix, "_pathway_id")]] else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- label
  out$primary_call <- as.character(genbank_table[[paste0(prefix, "_pathway_id")]][keep])
  out$reference_id <- dnmb_detail_join(
    if (paste0(prefix, "_reference_id") %in% names(genbank_table)) paste0("ref=", as.character(genbank_table[[paste0(prefix, "_reference_id")]][keep])) else NA_character_,
    if (paste0(prefix, "_step_id") %in% names(genbank_table)) paste0("step=", as.character(genbank_table[[paste0(prefix, "_step_id")]][keep])) else NA_character_
  )
  out$significance <- dnmb_detail_join(
    dnmb_detail_key_value("score", if (paste0(prefix, "_step_score") %in% names(genbank_table)) genbank_table[[paste0(prefix, "_step_score")]][keep] else NA, digits = 2),
    if (paste0(prefix, "_confidence") %in% names(genbank_table)) paste0("confidence=", as.character(genbank_table[[paste0(prefix, "_confidence")]][keep])) else NA_character_,
    if (paste0(prefix, "_on_best_path") %in% names(genbank_table)) paste0("best_path=", as.character(genbank_table[[paste0(prefix, "_on_best_path")]][keep])) else NA_character_
  )
  out$score_summary <- dnmb_detail_join(
    if (paste0(prefix, "_candidate_type") %in% names(genbank_table)) paste0("type=", as.character(genbank_table[[paste0(prefix, "_candidate_type")]][keep])) else NA_character_,
    if (paste0(prefix, "_pathway_desc") %in% names(genbank_table)) paste0("pathway_desc=", as.character(genbank_table[[paste0(prefix, "_pathway_desc")]][keep])) else NA_character_
  )
  out$context_summary <- dnmb_detail_join(
    if (paste0(prefix, "_support") %in% names(genbank_table)) as.character(genbank_table[[paste0(prefix, "_support")]][keep]) else NA_character_,
    if (paste0(prefix, "_locus_tag2") %in% names(genbank_table)) paste0("split_locus=", as.character(genbank_table[[paste0(prefix, "_locus_tag2")]][keep])) else NA_character_
  )
  out
}

dnmb_module_details_defensefinder <- function(genbank_table, base_cols) {
  label_col <- if ("DefenseFinder_system_id" %in% names(genbank_table)) genbank_table$DefenseFinder_system_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "DefenseFinder"
  out$primary_call <- if ("DefenseFinder_system_subtype" %in% names(genbank_table)) as.character(genbank_table$DefenseFinder_system_subtype[keep]) else as.character(genbank_table$DefenseFinder_system_type[keep])
  out$reference_id <- dnmb_detail_join(
    if ("DefenseFinder_system_id" %in% names(genbank_table)) paste0("sys=", as.character(genbank_table$DefenseFinder_system_id[keep])) else NA_character_,
    if ("DefenseFinder_gene_name" %in% names(genbank_table)) paste0("gene=", as.character(genbank_table$DefenseFinder_gene_name[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_join(
    dnmb_detail_key_value("system_score", if ("DefenseFinder_system_score" %in% names(genbank_table)) genbank_table$DefenseFinder_system_score[keep] else NA, digits = 3),
    dnmb_detail_key_value("wholeness", if ("DefenseFinder_system_wholeness" %in% names(genbank_table)) genbank_table$DefenseFinder_system_wholeness[keep] else NA, digits = 3),
    if ("DefenseFinder_hit_status" %in% names(genbank_table)) paste0("status=", as.character(genbank_table$DefenseFinder_hit_status[keep])) else NA_character_
  )
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("hit_score", if ("DefenseFinder_hit_score" %in% names(genbank_table)) genbank_table$DefenseFinder_hit_score[keep] else NA, digits = 3),
    dnmb_detail_key_value("i_eval", if ("DefenseFinder_hit_i_eval" %in% names(genbank_table)) genbank_table$DefenseFinder_hit_i_eval[keep] else NA, digits = 3),
    dnmb_detail_key_value("profile_cov", if ("DefenseFinder_hit_profile_cov" %in% names(genbank_table)) genbank_table$DefenseFinder_hit_profile_cov[keep] else NA, digits = 3),
    dnmb_detail_key_value("seq_cov", if ("DefenseFinder_hit_seq_cov" %in% names(genbank_table)) genbank_table$DefenseFinder_hit_seq_cov[keep] else NA, digits = 3)
  )
  out$context_summary <- dnmb_detail_join(
    if ("DefenseFinder_system_type" %in% names(genbank_table)) paste0("type=", as.character(genbank_table$DefenseFinder_system_type[keep])) else NA_character_,
    if ("DefenseFinder_system_activity" %in% names(genbank_table)) paste0("activity=", as.character(genbank_table$DefenseFinder_system_activity[keep])) else NA_character_,
    if ("DefenseFinder_profiles_in_system" %in% names(genbank_table)) paste0("profiles=", as.character(genbank_table$DefenseFinder_profiles_in_system[keep])) else NA_character_,
    if ("DefenseFinder_support" %in% names(genbank_table)) as.character(genbank_table$DefenseFinder_support[keep]) else NA_character_
  )
  out
}

dnmb_module_details_iselement <- function(genbank_table, base_cols) {
  label_col <- if ("ISelement_element_id" %in% names(genbank_table)) genbank_table$ISelement_element_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "ISelement"
  out$primary_call <- if ("ISelement_element_type" %in% names(genbank_table)) as.character(genbank_table$ISelement_element_type[keep]) else as.character(genbank_table$ISelement_element_id[keep])
  out$reference_id <- dnmb_detail_join(
    if ("ISelement_element_id" %in% names(genbank_table)) paste0("element=", as.character(genbank_table$ISelement_element_id[keep])) else NA_character_,
    if ("ISelement_family_system" %in% names(genbank_table)) paste0("family_system=", as.character(genbank_table$ISelement_family_system[keep])) else NA_character_
  )
  out$significance <- dnmb_detail_join(
    if ("ISelement_confidence" %in% names(genbank_table)) paste0("confidence=", as.character(genbank_table$ISelement_confidence[keep])) else NA_character_,
    dnmb_detail_key_value("score", if ("ISelement_confidence_score" %in% names(genbank_table)) genbank_table$ISelement_confidence_score[keep] else NA, digits = 2)
  )
  out$context_summary <- dnmb_detail_join(
    if ("ISelement_feature_type" %in% names(genbank_table)) paste0("feature_type=", as.character(genbank_table$ISelement_feature_type[keep])) else NA_character_,
    if ("ISelement_evidence_label" %in% names(genbank_table)) paste0("evidence=", as.character(genbank_table$ISelement_evidence_label[keep])) else NA_character_,
    dnmb_detail_key_value("supporting_genes", if ("ISelement_supporting_gene_count" %in% names(genbank_table)) genbank_table$ISelement_supporting_gene_count[keep] else NA, digits = 0),
    if ("ISelement_supporting_genes" %in% names(genbank_table)) as.character(genbank_table$ISelement_supporting_genes[keep]) else NA_character_,
    if ("ISelement_support" %in% names(genbank_table)) as.character(genbank_table$ISelement_support[keep]) else NA_character_
  )
  out
}

dnmb_module_details_prophage <- function(genbank_table, base_cols) {
  label_col <- if ("Prophage_prophage_id" %in% names(genbank_table)) genbank_table$Prophage_prophage_id else NULL
  if (is.null(label_col)) {
    return(data.frame())
  }
  keep <- !is.na(label_col) & nzchar(label_col)
  if (!any(keep)) {
    return(data.frame())
  }
  base <- dnmb_module_detail_base(genbank_table, base_cols, keep)
  out <- dnmb_module_detail_template(base)
  out$module_category <- "Prophage"
  out$primary_call <- as.character(genbank_table$Prophage_prophage_id[keep])
  out$reference_id <- if ("Prophage_prophage_id" %in% names(genbank_table)) paste0("region=", as.character(genbank_table$Prophage_prophage_id[keep])) else NA_character_
  out$significance <- dnmb_detail_join(
    dnmb_detail_key_value("rank", if ("Prophage_prophage_rank" %in% names(genbank_table)) genbank_table$Prophage_prophage_rank[keep] else NA, digits = 2),
    dnmb_detail_key_value("my_status", if ("Prophage_prophage_my_status" %in% names(genbank_table)) genbank_table$Prophage_prophage_my_status[keep] else NA, digits = 2),
    dnmb_detail_key_value("pp", if ("Prophage_prophage_pp" %in% names(genbank_table)) genbank_table$Prophage_prophage_pp[keep] else NA, digits = 2)
  )
  out$alignment_summary <- dnmb_detail_span(
    "region",
    if ("Prophage_prophage_start" %in% names(genbank_table)) genbank_table$Prophage_prophage_start[keep] else NA,
    if ("Prophage_prophage_end" %in% names(genbank_table)) genbank_table$Prophage_prophage_end[keep] else NA,
    NA
  )
  out$context_summary <- if ("Prophage_support" %in% names(genbank_table)) as.character(genbank_table$Prophage_support[keep]) else NA_character_
  out
}

dnmb_gapmind_candidate_file <- function(version = "aa", cwd = getwd()) {
  version <- tolower(as.character(version)[1])
  candidates <- if (identical(version, "aa")) {
    c(
      file.path(cwd, "dnmb_module_gapmindaa", "aa.sum.cand"),
      file.path(cwd, "dnmb_module_gapmind_aa", "aa.sum.cand"),
      file.path(cwd, "dnmb_module_gapmind", "aa.sum.cand")
    )
  } else {
    c(
      file.path(cwd, "dnmb_module_gapmindcarbon", "aa.sum.cand"),
      file.path(cwd, "dnmb_module_gapmind_carbon", "aa.sum.cand")
    )
  }
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) {
    return(NULL)
  }
  hits[[1]]
}

dnmb_gapmind_step_file <- function(version = "aa", cwd = getwd()) {
  version <- tolower(as.character(version)[1])
  candidates <- if (identical(version, "aa")) {
    c(
      file.path(cwd, "dnmb_module_gapmindaa", "aa.sum.steps"),
      file.path(cwd, "dnmb_module_gapmind_aa", "aa.sum.steps"),
      file.path(cwd, "dnmb_module_gapmind", "aa.sum.steps")
    )
  } else {
    c(
      file.path(cwd, "dnmb_module_gapmindcarbon", "aa.sum.steps"),
      file.path(cwd, "dnmb_module_gapmind_carbon", "aa.sum.steps")
    )
  }
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) {
    return(NULL)
  }
  hits[[1]]
}

dnmb_gapmind_pathway_table_file <- function(version = "aa", cwd = getwd()) {
  version <- tolower(as.character(version)[1])
  candidates <- c(
    file.path(dirname(cwd), "cache", "shared-module-cache", "db_modules", "gapmind", version, "PaperBLAST", "gaps", version, paste0(version, ".table")),
    file.path(cwd, "cache", "shared-module-cache", "db_modules", "gapmind", version, "PaperBLAST", "gaps", version, paste0(version, ".table")),
    file.path(dnmb_db_home(cache_root = NULL, create = FALSE), "gapmind", version, "PaperBLAST", "gaps", version, paste0(version, ".table"))
  )
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) {
    return(NULL)
  }
  hits[[1]]
}

dnmb_module_details_gapmind_full <- function(genbank_table, base_cols, version = "aa", label = "GapMind") {
  cand_path <- dnmb_gapmind_candidate_file(version = version)
  if (is.null(cand_path)) {
    return(data.frame())
  }
  step_path <- dnmb_gapmind_step_file(version = version)
  pathway_table_path <- dnmb_gapmind_pathway_table_file(version = version)
  pathway_tbl <- if (!is.null(pathway_table_path) && file.exists(pathway_table_path)) dnmb_gapmind_parse_pathway_table(pathway_table_path) else data.frame()
  candidates <- dnmb_gapmind_parse_candidates(cand_path, steps_path = step_path, pathway_table = pathway_tbl)
  if (!is.data.frame(candidates) || !nrow(candidates)) {
    return(data.frame())
  }

  candidates$locus_tag <- .dnmb_module_clean_annotation_key(candidates$locus_tag)
  keep <- !is.na(candidates$locus_tag) & candidates$locus_tag %in% genbank_table$locus_tag
  candidates <- candidates[keep, , drop = FALSE]
  if (!nrow(candidates)) {
    return(data.frame())
  }

  index <- match(candidates$locus_tag, genbank_table$locus_tag)
  base <- dnmb_module_detail_base(genbank_table, base_cols, index)
  out <- dnmb_module_detail_template(base)
  out$module_category <- label
  out$primary_call <- as.character(candidates$pathway)
  out$reference_id <- dnmb_detail_join(
    if ("curatedIds" %in% names(candidates)) paste0("ref=", as.character(candidates$curatedIds)) else NA_character_,
    if ("hmmId" %in% names(candidates)) paste0("hmm=", as.character(candidates$hmmId)) else NA_character_,
    paste0("step=", as.character(candidates$step))
  )
  out$significance <- dnmb_detail_join(
    dnmb_detail_key_value("score", candidates$score, digits = 2),
    if ("on_best_path" %in% names(candidates)) paste0("best_path=", as.character(candidates$on_best_path)) else NA_character_
  )
  out$score_summary <- dnmb_detail_join(
    dnmb_detail_key_value("blastBits", if ("blastBits" %in% names(candidates)) candidates$blastBits else NA, digits = 3),
    dnmb_detail_key_value("identity", if ("identity" %in% names(candidates)) candidates$identity else NA, digits = 3),
    dnmb_detail_key_value("blastCoverage", if ("blastCoverage" %in% names(candidates)) candidates$blastCoverage else NA, digits = 3),
    dnmb_detail_key_value("hmmBits", if ("hmmBits" %in% names(candidates)) candidates$hmmBits else NA, digits = 3),
    dnmb_detail_key_value("hmmCoverage", if ("hmmCoverage" %in% names(candidates)) candidates$hmmCoverage else NA, digits = 3)
  )
  out$alignment_summary <- dnmb_detail_join(
    dnmb_detail_key_value("blastScore", if ("blastScore" %in% names(candidates)) candidates$blastScore else NA, digits = 2),
    dnmb_detail_key_value("hmmScore", if ("hmmScore" %in% names(candidates)) candidates$hmmScore else NA, digits = 2),
    if ("locus_tag2" %in% names(candidates)) paste0("split_locus=", as.character(candidates$locus_tag2)) else NA_character_
  )
  out$context_summary <- dnmb_detail_join(
    if ("pathway_desc" %in% names(candidates)) paste0("pathway_desc=", as.character(candidates$pathway_desc)) else NA_character_,
    if ("step" %in% names(candidates)) paste0("step=", as.character(candidates$step)) else NA_character_,
    if ("curatedDesc" %in% names(candidates)) paste0("curated_desc=", as.character(candidates$curatedDesc)) else NA_character_,
    if ("otherIds" %in% names(candidates)) paste0("other=", as.character(candidates$otherIds)) else NA_character_
  )
  out
}

dnmb_module_location_summary <- function(contig, start, end, direction) {
  contig <- as.character(contig)
  start <- suppressWarnings(as.numeric(start))
  end <- suppressWarnings(as.numeric(end))
  direction <- as.character(direction)
  vapply(seq_along(contig), function(i) {
    parts <- character()
    if (!is.na(contig[[i]]) && nzchar(contig[[i]])) parts <- c(parts, contig[[i]])
    coord <- NULL
    if (!is.na(start[[i]]) || !is.na(end[[i]])) {
      coord <- paste0(
        ifelse(is.na(start[[i]]), "?", format(start[[i]], scientific = FALSE, trim = TRUE)),
        "-",
        ifelse(is.na(end[[i]]), "?", format(end[[i]], scientific = FALSE, trim = TRUE))
      )
    }
    if (!is.null(coord)) parts <- c(parts, coord)
    if (!is.na(direction[[i]]) && nzchar(direction[[i]])) parts <- c(parts, paste0("dir=", direction[[i]]))
    if (!length(parts)) return(NA_character_)
    paste(parts, collapse = " | ")
  }, character(1))
}

dnmb_detail_key_value <- function(label, value, digits = 3) {
  value <- suppressWarnings(as.numeric(value))
  vapply(seq_along(value), function(i) {
    if (is.na(value[[i]])) {
      return(NA_character_)
    }
    paste0(label, "=", format(signif(value[[i]], digits), trim = TRUE))
  }, character(1))
}

dnmb_detail_span <- function(label, start, end, total = NA) {
  start <- suppressWarnings(as.numeric(start))
  end <- suppressWarnings(as.numeric(end))
  total <- suppressWarnings(as.numeric(total))
  n <- max(length(start), length(end), length(total))
  start <- rep_len(start, n)
  end <- rep_len(end, n)
  total <- rep_len(total, n)
  vapply(seq_along(start), function(i) {
    if (is.na(start[[i]]) && is.na(end[[i]]) && is.na(total[[i]])) {
      return(NA_character_)
    }
    span <- paste0(
      ifelse(is.na(start[[i]]), "?", format(start[[i]], scientific = FALSE, trim = TRUE)),
      "-",
      ifelse(is.na(end[[i]]), "?", format(end[[i]], scientific = FALSE, trim = TRUE))
    )
    if (!is.na(total[[i]])) {
      span <- paste0(span, "/", format(total[[i]], scientific = FALSE, trim = TRUE))
    }
    paste0(label, "=", span)
  }, character(1))
}

dnmb_detail_join <- function(...) {
  pieces <- list(...)
  if (!length(pieces)) {
    return(character())
  }
  n <- max(vapply(pieces, length, integer(1)))
  pieces <- lapply(pieces, function(x) rep_len(as.character(x), n))
  vapply(seq_len(n), function(i) {
    vals <- vapply(pieces, function(x) x[[i]], character(1))
    vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
    if (!length(vals)) {
      return(NA_character_)
    }
    paste(vals, collapse = "; ")
  }, character(1))
}

dnmb_detail_links <- function(expasy_link, brenda_link) {
  expasy_link <- as.character(expasy_link)
  brenda_link <- as.character(brenda_link)
  n <- max(length(expasy_link), length(brenda_link))
  expasy_link <- rep_len(expasy_link, n)
  brenda_link <- rep_len(brenda_link, n)
  vapply(seq_len(n), function(i) {
    vals <- c(
      if (!is.na(expasy_link[[i]]) && nzchar(trimws(expasy_link[[i]]))) paste0("ExPASy=", expasy_link[[i]]) else NA_character_,
      if (!is.na(brenda_link[[i]]) && nzchar(trimws(brenda_link[[i]]))) paste0("BRENDA=", brenda_link[[i]]) else NA_character_
    )
    vals <- vals[!is.na(vals)]
    if (!length(vals)) {
      return(NA_character_)
    }
    paste(vals, collapse = "; ")
  }, character(1))
}

dnmb_sanitize_output_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x <- gsub("[\r\n\t]+", " ", x, perl = TRUE)
  x <- gsub("[[:space:]]{2,}", " ", x, perl = TRUE)
  trimws(x)
}

dnmb_strip_automated_note_phrase <- function(x) {
  x <- dnmb_sanitize_output_text(x)
  x <- gsub("Derived by automated computational analysis using\\s*", "", x, perl = TRUE)
  x <- gsub("\\s*;\\s*;\\s*", "; ", x, perl = TRUE)
  x <- gsub("\\s*,\\s*,\\s*", ", ", x, perl = TRUE)
  x <- gsub("[[:space:]]{2,}", " ", x, perl = TRUE)
  trimws(x)
}

dnmb_apply_workbook_sheet_order <- function(wb) {
  desired <- c(
    "1.GenBank_table",
    "2.Codon_usage",
    "3.tRNA_anticodon",
    "4.tRNA_distribution",
    "5.RBS",
    "6.CRISPR_table",
    "7.InterPro_site",
    "8.Module_details",
    "9.Prophage_summary"
  )
  existing <- names(wb)
  ordered <- c(intersect(desired, existing), setdiff(existing, desired))
  if (length(ordered) == length(existing)) {
    openxlsx::worksheetOrder(wb) <- match(ordered, existing)
  }
  invisible(wb)
}

dnmb_excel_hyperlink_label <- function(column_name) {
  column_name <- tolower(as.character(column_name)[1])
  if (grepl("expasy", column_name, fixed = TRUE)) {
    return("ExPASy")
  }
  if (grepl("brenda", column_name, fixed = TRUE)) {
    return("BRENDA")
  }
  "Link"
}

dnmb_header_group_for_column <- function(column_name) {
  column_name <- as.character(column_name)[1]
  if (grepl("^CLEAN_", column_name)) {
    return("CLEAN")
  }
  if (grepl("^dbCAN_", column_name)) {
    return("dbCAN")
  }
  if (grepl("^MEROPS_", column_name)) {
    return("MEROPS")
  }
  if (grepl("^PAZy_", column_name)) {
    return("PAZy")
  }
  if (grepl("^GapMind", column_name)) {
    return("GapMind")
  }
  if (grepl("^DefenseFinder_", column_name)) {
    return("DefenseFinder")
  }
  if (grepl("^dbAPIS_", column_name)) {
    return("dbAPIS")
  }
  if (grepl("^AcrFinder_", column_name)) {
    return("AcrFinder")
  }
  if (grepl("^Promotech_", column_name)) {
    return("Promotech")
  }
  if (grepl("^mRNAcal_", column_name)) {
    return("mRNAcal")
  }
  if (grepl("^PADLOC_", column_name)) {
    return("PADLOC")
  }
  if (grepl("^DefensePredictor_", column_name)) {
    return("DefensePredictor")
  }
  if (grepl("^ISelement_", column_name)) {
    return("ISelement")
  }
  if (grepl("^Prophage_", column_name)) {
    return("Prophage")
  }
  "DNMB"
}

dnmb_header_group_styles <- function() {
  list(
    DNMB = openxlsx::createStyle(
      fgFill = "#DCEAF7",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#8BA7C1"
    ),
    CLEAN = openxlsx::createStyle(
      fgFill = "#E5F3E1",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#8FB27E"
    ),
    dbCAN = openxlsx::createStyle(
      fgFill = "#FBE5D6",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#C98D61"
    ),
    MEROPS = openxlsx::createStyle(
      fgFill = "#F5DDE3",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#B56D7E"
    ),
    PAZy = openxlsx::createStyle(
      fgFill = "#F7F0C8",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#B49B3C"
    ),
    GapMind = openxlsx::createStyle(
      fgFill = "#D8F1F0",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#4D8E92"
    ),
    DefenseFinder = openxlsx::createStyle(
      fgFill = "#E9E2F7",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#7A67AF"
    ),
    dbAPIS = openxlsx::createStyle(
      fgFill = "#DDF6D2",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#52B846"
    ),
    AcrFinder = openxlsx::createStyle(
      fgFill = "#D7F0E6",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#3B8F72"
    ),
    Promotech = openxlsx::createStyle(
      fgFill = "#EAF3D8",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#7B9B38"
    ),
    mRNAcal = openxlsx::createStyle(
      fgFill = "#E1EEF7",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#3A79A8"
    ),
    PADLOC = openxlsx::createStyle(
      fgFill = "#FDE9D9",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#CC7A29"
    ),
    DefensePredictor = openxlsx::createStyle(
      fgFill = "#F8DDEA",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#B34B7F"
    ),
    REBASE = openxlsx::createStyle(
      fgFill = "#E8F2EA",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#4F8A61"
    ),
    ISelement = openxlsx::createStyle(
      fgFill = "#E2F1EC",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#5F8F80"
    ),
    Prophage = openxlsx::createStyle(
      fgFill = "#F2E5D7",
      textDecoration = "bold",
      halign = "center",
      valign = "center",
      border = "Bottom",
      borderColour = "#B7814F"
    )
  )
}

dnmb_apply_group_header_styles <- function(wb, sheet, table, startRow = 1L, startCol = 1L) {
  if (!is.data.frame(table) || !ncol(table)) {
    return(invisible(NULL))
  }

  styles <- dnmb_header_group_styles()
  column_names <- names(table)
  for (index in seq_along(column_names)) {
    group_name <- dnmb_header_group_for_column(column_names[[index]])
    style <- styles[[group_name]]
    if (is.null(style)) {
      next
    }
    openxlsx::addStyle(
      wb,
      sheet = sheet,
      style = style,
      rows = startRow,
      cols = startCol + index - 1L,
      gridExpand = FALSE,
      stack = TRUE
    )
  }

  invisible(NULL)
}

dnmb_apply_module_category_styles <- function(wb, sheet, table, startRow = 1L, startCol = 1L) {
  if (!is.data.frame(table) || !nrow(table) || !"module_category" %in% names(table)) {
    return(invisible(NULL))
  }
  styles <- dnmb_header_group_styles()
  module_col <- startCol + match("module_category", names(table)) - 1L
  categories <- as.character(table$module_category)
  for (group_name in unique(categories[!is.na(categories) & nzchar(categories)])) {
    style_group <- gsub("(AA|Carbon)$", "", group_name)
    style <- styles[[style_group]]
    if (is.null(style)) {
      next
    }
    rows <- startRow + which(categories == group_name)
    openxlsx::addStyle(
      wb,
      sheet = sheet,
      style = style,
      rows = rows,
      cols = module_col,
      gridExpand = FALSE,
      stack = TRUE
    )
  }
  invisible(NULL)
}

dnmb_write_hyperlink_columns <- function(wb, sheet, table, startRow = 1L, startCol = 1L) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(invisible(NULL))
  }

  link_cols <- grep("_link$", names(table), value = TRUE)
  if (!length(link_cols)) {
    return(invisible(NULL))
  }

  for (column_name in link_cols) {
    values <- as.character(table[[column_name]])
    values <- trimws(values)
    formulas <- rep(NA_character_, length(values))
    keep <- !is.na(values) & nzchar(values)
    if (any(keep)) {
      label <- dnmb_excel_hyperlink_label(column_name)
      escaped <- gsub("\"", "\"\"", values[keep], fixed = TRUE)
      formulas[keep] <- paste0('=HYPERLINK("', escaped, '","', label, '")')
    }
    if (any(!is.na(formulas))) {
      col_index <- startCol + match(column_name, names(table)) - 1L
      for (row_index in which(!is.na(formulas))) {
        openxlsx::writeFormula(
          wb,
          sheet = sheet,
          x = formulas[[row_index]],
          startCol = col_index,
          startRow = startRow + row_index
        )
      }
    }
  }

  invisible(NULL)
}

dnmb_apply_numeric_display_styles <- function(wb, sheet, table, startRow = 1L, startCol = 1L) {
  if (!is.data.frame(table) || !nrow(table)) {
    return(invisible(NULL))
  }

  numeric_cols <- intersect(
    c("CLEAN_best_distance", "CLEAN_pvalue_best_pvalue", "distance", "pvalue"),
    names(table)
  )
  if (!length(numeric_cols)) {
    return(invisible(NULL))
  }

  three_decimal_style <- openxlsx::createStyle(numFmt = "0.000")
  for (column_name in numeric_cols) {
    col_index <- startCol + match(column_name, names(table)) - 1L
    openxlsx::addStyle(
      wb,
      sheet = sheet,
      style = three_decimal_style,
      rows = seq.int(startRow + 1L, startRow + nrow(table)),
      cols = col_index,
      gridExpand = TRUE,
      stack = TRUE
    )
  }

  invisible(NULL)
}

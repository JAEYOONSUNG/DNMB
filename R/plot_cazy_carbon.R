.dnmb_cct_cazy_substrate_map <- function() {
  # Mapping: CAZy family -> substrate acted upon -> degradation product
  # Codex-verified 2026-03-26:
  #   GH13_31: starch/alpha-glucan only (removed fructan/sucrose dual mapping)
  #   GH18: endo-chitinase only (chitin->chitooligosaccharides), NOT exo
  #   GH20: exo-hexosaminidase (GlcNAc release) — separate from GH18

  #   GH130_4: phosphorylase (Man-1P release), not simple hydrolase
  #   GH5/GH9: cellulase (endo), only shown if genome has actual hits
  rbind(
    # Starch/maltose pathway (alpha-glucans)
    data.frame(family = c("GH13","GH13_5","GH13_20","GH13_31","GH13_39","GH13_45"),
               substrate = "maltose", product = "glucose", stringsAsFactors = FALSE, row.names = NULL),
    # GH4: 6-phospho-beta-glucosidase (maltose-6P pathway)
    data.frame(family = "GH4",
               substrate = "maltose", product = "glucose-6P", stringsAsFactors = FALSE, row.names = NULL),
    # Cellobiose (beta-glucans)
    data.frame(family = c("GH1","GH5","GH9","GH44","GH48"),
               substrate = "cellobiose", product = "glucose", stringsAsFactors = FALSE, row.names = NULL),
    # Sucrose (GH32=invertase, GH68=levansucrase)
    data.frame(family = c("GH32","GH68"),
               substrate = "sucrose", product = "fructose+glucose", stringsAsFactors = FALSE, row.names = NULL),
    # Chitin: GH18/GH19 are endo-chitinases (chitin->chitooligosaccharides)
    data.frame(family = c("GH18","GH19"),
               substrate = "NAG", product = "chitooligosaccharides", stringsAsFactors = FALSE, row.names = NULL),
    # GH20: exo-hexosaminidase (chitobiose->GlcNAc)
    data.frame(family = "GH20",
               substrate = "NAG", product = "GlcNAc", stringsAsFactors = FALSE, row.names = NULL),
    # Peptidoglycan
    data.frame(family = c("GH23","GH24","GH25","GH73"),
               substrate = "peptidoglycan", product = "MurNAc+GlcNAc", stringsAsFactors = FALSE, row.names = NULL),
    # Galactose
    data.frame(family = c("GH36","GH42","GH2"),
               substrate = "galactose", product = "galactose", stringsAsFactors = FALSE, row.names = NULL),
    # Mannose: GH92=exo-mannosidase, GH130/GH130_4=phosphorylase (Man-1P)
    data.frame(family = "GH92",
               substrate = "mannose", product = "mannose", stringsAsFactors = FALSE, row.names = NULL),
    data.frame(family = c("GH130","GH130_4"),
               substrate = "mannose", product = "Man-1P", stringsAsFactors = FALSE, row.names = NULL),
    # Xylose (xylan degradation)
    data.frame(family = c("GH10","GH11","GH43"),
               substrate = "xylose", product = "xylose", stringsAsFactors = FALSE, row.names = NULL),
    # Trehalose
    data.frame(family = c("GH37","GH65"),
               substrate = "trehalose", product = "glucose", stringsAsFactors = FALSE, row.names = NULL),
    # Arabinose
    data.frame(family = c("GH51","GH54"),
               substrate = "arabinose", product = "arabinose", stringsAsFactors = FALSE, row.names = NULL),
    # Glucuronate (pectin-related)
    data.frame(family = c("GH28","GH105"),
               substrate = "glucuronate", product = "glucuronate", stringsAsFactors = FALSE, row.names = NULL),
    # Glucosamine
    data.frame(family = "GH3",
               substrate = "glucosamine", product = "glucosamine-6P", stringsAsFactors = FALSE, row.names = NULL)
  )
}

# --- Extract ALL transport data (rules + steps + PTS) from GapMind Carbon ---
.dnmb_cct_transport_extract <- function(output_dir) {
  base_dir <- output_dir

  # --- 1. Read aa.sum.rules to identify transport rules per pathway ---
  rules_file <- NULL
  for (f in c(file.path(base_dir, "dnmb_module_gapmindcarbon", "aa.sum.rules"),
              file.path(base_dir, "dnmb_module_gapmind_carbon", "aa.sum.rules"),
              file.path(base_dir, "aa.sum.rules"))) {
    if (file.exists(f)) { rules_file <- f; break }
  }
  steps_file <- NULL
  for (f in c(file.path(base_dir, "dnmb_module_gapmindcarbon", "aa.sum.steps"),
              file.path(base_dir, "dnmb_module_gapmind_carbon", "aa.sum.steps"),
              file.path(base_dir, "aa.sum.steps"))) {
    if (file.exists(f)) { steps_file <- f; break }
  }
  if (is.null(rules_file) && is.null(steps_file)) return(NULL)

  # --- 2. From rules: collect step names ONLY from *-transport rules ---
  # Do NOT include *-utilization rules (they mix metabolic enzymes with PTS)
  rule_transport_steps <- character(0)
  if (!is.null(rules_file)) {
    rules <- utils::read.delim(rules_file, stringsAsFactors = FALSE, row.names = NULL)
    is_tr <- grepl("-transport$", rules$rule)
    tr_rules <- rules[is_tr, , drop = FALSE]
    if (nrow(tr_rules) > 0) {
      ep <- unlist(strsplit(tr_rules$expandedPath, "\\s+"))
      p2 <- unlist(strsplit(tr_rules$path2, "\\s+"))
      rule_transport_steps <- unique(c(ep, p2))
      rule_transport_steps <- rule_transport_steps[nzchar(rule_transport_steps) &
                                                    !is.na(rule_transport_steps)]
    }
  }

  # --- 3. From steps: get locus_tag assignments ---
  if (is.null(steps_file)) return(NULL)
  st <- utils::read.delim(steps_file, stringsAsFactors = FALSE, row.names = NULL)
  req <- c("pathway", "step", "onBestPath", "score")
  if (!all(req %in% names(st))) return(NULL)

  # Transport/permease gene regex — actual membrane transport proteins only
  # PTS components (phosphotransferase system)
  # ABC/MFS/SSS permeases and binding proteins
  # Do NOT match metabolic enzymes (galK, nagA, pgmA, glk, etc.)
  pts_transport_re <- paste0(
    "^ptsG-crr|^ptsG$|^ptsH$|^ptsI$|^ptsS$|^crr$|",
    "^nagEIIA|^nagEcb|^nagEcba|^nagPcb|^nagP$|^nagF$|^nag3$|",
    "^malEII|^malE[_1FG]|^malF|^malG|^malK|^musE|^musF|^musG|^musK|",
    "^manP$|^manX$|^manY$|^manZ$|^manMFS|^mapP$|",
    "^bglF$|^bglG$|^bglT$|^ascB$|",
    "^celEII|",
    "^sacP$|^scrT$|^sut$|^sut1$|^cscB$|",
    "^fruII|^fruA$|^fruB$|^fruD$|^fruE$|^fruF$|^fruG$|^fruI$|^fruK$|^fruP$|",
    "^levD$|^levDE$|^levE$|^levF$|^levG$|",
    "^treB$|^treEII|^TRET1$|",
    "^galP$|^gal2$|^ptcA$|^ptcB$|^ptcEIIC$|^lacP$|^lacA\\'$|^lacB\\'$|^lacC\\'$|",
    "^xylT$|^araE$|^araS$|^araT$|^araU$|^araV$|^exuT$|",
    "^cdt$|^cbt[A-F]$|^ceb[EFG]$|^cbpB$|^cbpC$|",
    "^MFS-|^SSS-|^SWEET|^SemiSWEET|^PAST-|^SGLT|",
    "^glcU|^glcS$|^glcT$|^glcV$|^glcP$|^gluP$|",
    "^thu[EFGK]$|^agl[EFGK]|^gts[ABCD]$|",
    "^mgl[ABC]$|^msiK$|^msdB|^msdC|",
    "^STP6$|^MAL11$|^Slc2a5$|^SLC45|^uhpT$|",
    "^kguT$|^scrT$|^larD$|^dctA$|^alsT$|^cycA$|^kgtP$|^pcaK$")
  is_transport <- (st$step %in% rule_transport_steps) |
                  grepl(pts_transport_re, st$step)
  tr <- st[is_transport, , drop = FALSE]
  tr <- tr[!is.na(tr$locusId) & nzchar(tr$locusId), , drop = FALSE]
  # Keep score == 0 assignments as low-confidence transport evidence.
  # GapMind often places weak but still meaningful transporter assignments on best path with score 0.
  tr <- tr[!is.na(tr$score) & tr$score >= 0, , drop = FALSE]
  tr <- tr[!duplicated(tr[, c("pathway", "step", "locusId")]), , drop = FALSE]
  if (nrow(tr) == 0) return(NULL)

  # PTS classification (order matters: ptsG-crr before crr)
  tr$is_pts <- grepl(paste0(
    "^ptsG-crr|^ptsG$|^ptsH$|^ptsI$|",
    "^nagEIIA|^nagEcb|^nagEcba|^nagF$|^nagPcb|",
    "^malEII|^manP$|",
    "^bgl[FG]$|^ascB$|^celEII|^sac[P]$|",
    "^fru[ABDI]|^fruII|^levD|^levE|^levF|^levG|",
    "^tre[AB]|^treEII|^crr$"), tr$step, ignore.case = TRUE)

  tr$confidence <- dplyr::case_when(
    tr$score >= 2 ~ "high",  tr$score >= 1 ~ "medium",
    tr$score >  0 ~ "low",   TRUE ~ "none")
  tr$locus_tag <- tr$locusId
  tr$gene_name <- if ("sysName" %in% names(tr)) tr$sysName else tr$step

  # --- 4. Include locusId2 alternatives ---
  if ("locusId2" %in% names(tr)) {
    alt <- tr[!is.na(tr$locusId2) & nzchar(as.character(tr$locusId2)), , drop = FALSE]
    if (nrow(alt) > 0) {
      alt2 <- alt
      alt2$locus_tag <- alt2$locusId2
      alt2$gene_name <- if ("sysName2" %in% names(alt2)) alt2$sysName2 else alt2$step
      alt2$score <- if ("score2" %in% names(alt2)) as.numeric(alt2$score2) else 0
      alt2 <- alt2[!is.na(alt2$score) & alt2$score >= 0, , drop = FALSE]
      if (nrow(alt2) > 0) {
        alt2$confidence <- dplyr::case_when(
          alt2$score >= 2 ~ "high",  alt2$score >= 1 ~ "medium",
          alt2$score >  0 ~ "low",   TRUE ~ "none")
        tr <- rbind(tr, alt2)
      }
    }
  }
  tr <- tr[!duplicated(tr[, c("pathway", "step", "locus_tag")]), , drop = FALSE]
  tr
}

# --- Extract PTS components from aa.sum.cand for richer detail ---
.dnmb_cct_pts_detail <- function(output_dir) {
  cand_file <- NULL
  for (f in c(file.path(output_dir, "dnmb_module_gapmindcarbon", "aa.sum.cand"),
              file.path(output_dir, "dnmb_module_gapmind_carbon", "aa.sum.cand"),
              file.path(output_dir, "aa.sum.cand"))) {
    if (file.exists(f)) { cand_file <- f; break }
  }
  if (is.null(cand_file)) return(NULL)
  ca <- utils::read.delim(cand_file, stringsAsFactors = FALSE, row.names = NULL)
  pts_re <- paste0(
    "pts[GHIS]|ptsG-crr|crr$|",
    "nag[EFPK]|nagEIIA|nagEcb|nagPcb|nagEcba|nagF$|",
    "malEII|manP$|",
    "bgl[FG]$|ascB|celEII|sac[P]$|",
    "fruII|fru[ABDI]|levD|levE|levF|levG|",
    "tre[AB]|treEII|",
    "PTS|phosphotransferase")
  pts_rows <- ca[grepl(pts_re, ca$step, ignore.case = TRUE) |
                 grepl(pts_re, ca$desc, ignore.case = TRUE), , drop = FALSE]
  if (nrow(pts_rows) == 0) return(NULL)
  pts_rows$locus_tag <- if ("locusId" %in% names(pts_rows)) pts_rows$locusId else NA_character_
  # Component classification (ORDER MATTERS: multi-domain fusions first)
  pts_rows$pts_component <- dplyr::case_when(
    grepl("ptsG-crr|nagEcba|malEIICBA|fruII-ABC", pts_rows$step) ~ "EIICBA",
    grepl("malEIICB", pts_rows$step)                              ~ "EIICB",
    grepl("nagF", pts_rows$step)                                  ~ "EI-HPr-EIIA",
    grepl("^ptsI$|\\bEI$", pts_rows$step)                        ~ "EI",
    grepl("^ptsH$|\\bHPr\\b", pts_rows$step)                     ~ "HPr",
    grepl("EIIA|^crr$|nagEIIA|malEIIA|fruII-A$", pts_rows$step)  ~ "EIIA",
    grepl("EIIB|nagPcb|nagEcb|fruII-B$|celEIIB|sacP$|ascB$|bglF$|malEIIC", pts_rows$step) ~ "EIIB/C",
    grepl("EIIC|^ptsG$|fruII-C$|celEIIC|manP$|levG$|treB$", pts_rows$step) ~ "EIIC",
    TRUE ~ "PTS"
  )
  # Include locusId2
  if ("locusId2" %in% names(pts_rows)) {
    alt <- pts_rows[!is.na(pts_rows$locusId2) & nzchar(as.character(pts_rows$locusId2)), , drop = FALSE]
    if (nrow(alt) > 0) {
      alt2 <- alt; alt2$locus_tag <- alt2$locusId2
      pts_rows <- rbind(pts_rows, alt2)
    }
  }
  pts_rows <- pts_rows[!is.na(pts_rows$locus_tag) & nzchar(pts_rows$locus_tag), , drop = FALSE]
  pts_rows <- pts_rows[!duplicated(pts_rows[, c("pathway","step","locus_tag")]), , drop = FALSE]
  pts_rows
}

# --- Extract CAZy hits with full annotation ---
.dnmb_cct_cazy_extract <- function(genbank_table, output_dir = NULL) {
  if (is.list(genbank_table) && !is.data.frame(genbank_table) &&
      "features" %in% names(genbank_table))
    genbank_table <- genbank_table$features
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE, row.names = NULL)
  family_col <- .dnmb_pick_column(tbl, c("dbCAN_family_id", "family_id"))
  # Fallback: read dbCAN merged results from module output directory
  if (is.null(family_col) && !is.null(output_dir)) {
    dbcan_tsv <- file.path(output_dir, "dnmb_module_dbcan", "dbcan_merged_fixed.tsv")
    if (file.exists(dbcan_tsv)) {
      tbl <- utils::read.delim(dbcan_tsv, stringsAsFactors = FALSE, row.names = NULL)
      family_col <- .dnmb_pick_column(tbl, c("dbCAN_family_id", "family_id"))
    }
  }
  if (is.null(family_col)) return(NULL)
  tbl <- tbl[!is.na(tbl[[family_col]]) & base::nzchar(tbl[[family_col]]), , drop = FALSE]
  if (!base::nrow(tbl)) return(NULL)
  cgc_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_id"))
  cgc_type_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_gene_type"))
  data.frame(
    locus_tag = tbl$locus_tag,
    gene = if ("gene" %in% names(tbl)) tbl$gene else NA_character_,
    gene_product = if ("product" %in% names(tbl)) tbl$product else NA_character_,
    family = as.character(tbl[[family_col]]),
    contig = if ("contig" %in% names(tbl)) tbl$contig else NA_character_,
    start = if ("start" %in% names(tbl)) as.numeric(tbl$start) else NA_real_,
    end = if ("end" %in% names(tbl)) as.numeric(tbl$end) else NA_real_,
    direction = if ("direction" %in% names(tbl)) tbl$direction else NA_character_,
    cgc_id = if (!is.null(cgc_col)) tbl[[cgc_col]] else NA_character_,
    cgc_type = if (!is.null(cgc_type_col)) tbl[[cgc_type_col]] else NA_character_,
    stringsAsFactors = FALSE
  )
}

# --- Infer operons: CGC + proximity-based ---
.dnmb_cct_infer_operons <- function(genbank_table, cazy_hits, transport_data) {
  if (is.list(genbank_table) && !is.data.frame(genbank_table) &&
      "features" %in% names(genbank_table))
    genbank_table <- genbank_table$features
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE, row.names = NULL)
  # Alias: genbank parser returns "strand"; operon logic expects "direction"
  if (!"direction" %in% names(tbl) && "strand" %in% names(tbl))
    tbl$direction <- tbl$strand
  if (!all(c("locus_tag","contig","start","end","direction") %in% names(tbl))) return(NULL)
  # Collect all relevant locus_tags
  relevant_lt <- unique(c(
    cazy_hits$locus_tag,
    if (!is.null(transport_data)) transport_data$locus_tag else character(0)
  ))
  relevant_lt <- relevant_lt[!is.na(relevant_lt) & nzchar(relevant_lt)]
  if (length(relevant_lt) == 0) return(NULL)

  # Method 1: CGC clusters from dbCAN
  cgc_col <- .dnmb_pick_column(tbl, c("dbCAN_dbcan_cgc_id"))
  cgc_groups <- list()
  if (!is.null(cgc_col)) {
    cgc_tbl <- tbl[!is.na(tbl[[cgc_col]]) & nzchar(tbl[[cgc_col]]), , drop = FALSE]
    if (nrow(cgc_tbl) > 0) {
      for (cid in unique(cgc_tbl[[cgc_col]])) {
        members <- cgc_tbl$locus_tag[cgc_tbl[[cgc_col]] == cid]
        # Only keep if at least one relevant gene
        if (any(members %in% relevant_lt)) {
          short_id <- sub(".*\\|", "", cid)  # "CGC1" from "contig|CGC1"
          cgc_groups[[short_id]] <- members
        }
      }
    }
  }

  # Method 2: Proximity-based (same strand, <200bp gap)
  sub_tbl <- tbl[tbl$locus_tag %in% relevant_lt, , drop = FALSE]
  sub_tbl <- sub_tbl[order(sub_tbl$contig, as.numeric(sub_tbl$start)), , drop = FALSE]
  prox_groups <- list()
  if (nrow(sub_tbl) >= 2) {
    grp_id <- 1; current_grp <- c(sub_tbl$locus_tag[1])
    for (i in 2:nrow(sub_tbl)) {
      same_contig <- sub_tbl$contig[i] == sub_tbl$contig[i-1]
      same_strand <- sub_tbl$direction[i] == sub_tbl$direction[i-1]
      gap <- as.numeric(sub_tbl$start[i]) - as.numeric(sub_tbl$end[i-1])
      if (same_contig && same_strand && !is.na(gap) && gap >= 0 && gap < 200) {
        current_grp <- c(current_grp, sub_tbl$locus_tag[i])
      } else {
        if (length(current_grp) >= 2) {
          prox_groups[[paste0("prox_", grp_id)]] <- current_grp
          grp_id <- grp_id + 1
        }
        current_grp <- c(sub_tbl$locus_tag[i])
      }
    }
    if (length(current_grp) >= 2) {
      prox_groups[[paste0("prox_", grp_id)]] <- current_grp
    }
  }

  # Merge: CGC takes priority
  all_groups <- c(cgc_groups, prox_groups)
  # Build locus_tag -> operon_id map
  lt_to_operon <- character(0)
  for (gname in names(all_groups)) {
    for (lt in all_groups[[gname]]) {
      if (is.na(lt_to_operon[lt]) || !nzchar(lt_to_operon[lt])) {
        lt_to_operon[lt] <- gname
      }
    }
  }
  list(groups = all_groups, lt_map = lt_to_operon)
}

# --- Multi-evidence confidence scoring ---
.dnmb_cct_evidence_score <- function(cazy_match, gapmind_match, operon_match, product_match) {
  score <- sum(c(cazy_match, gapmind_match, operon_match, product_match), na.rm = TRUE)
  dplyr::case_when(score >= 3 ~ "very_high", score == 2 ~ "high",
                    score == 1 ~ "medium", TRUE ~ "low")
}

# --- Build comprehensive pathway rows for grid layout ---
.dnmb_cct_build_rows <- function(cazy_hits, transport_data, pts_detail,
                                  operon_info, cazy_map, genbank_table) {
  # Join CAZy hits to substrate map
  cazy_anno <- merge(cazy_hits, cazy_map, by = "family", all.x = FALSE)
  if (nrow(cazy_anno) > 0) {
    cazy_anno$enzyme_label <- ifelse(
      !is.na(cazy_anno$gene) & nzchar(cazy_anno$gene),
      paste0(cazy_anno$gene, " (", cazy_anno$family, ")"),
      cazy_anno$family)
  }

  # --- Bug 3: Operon-based PTS<->CAZy cross-matching ---
  # Before substrate assignment: if PTS and CAZy share an operon,
  # propagate substrate info bidirectionally
  if (!is.null(operon_info) && !is.null(pts_detail) && nrow(cazy_anno) > 0) {
    cazy_lts <- unique(cazy_anno$locus_tag)
    pts_lts <- unique(pts_detail$locus_tag[!is.na(pts_detail$locus_tag)])
    for (gname in names(operon_info$groups)) {
      members <- operon_info$groups[[gname]]
      cazy_in_op <- cazy_lts[cazy_lts %in% members]
      pts_in_op <- pts_lts[pts_lts %in% members]
      if (length(cazy_in_op) > 0 && length(pts_in_op) > 0) {
        # Get CAZy substrates for this operon
        cazy_subs_op <- unique(cazy_anno$substrate[cazy_anno$locus_tag %in% cazy_in_op])
        # Get PTS pathways for this operon
        pts_pathways_op <- unique(pts_detail$pathway[pts_detail$locus_tag %in% pts_in_op])
        # Cross-assign: add PTS pathways to CAZy substrates (if PTS pathway != CAZy substrate)
        for (ps in pts_pathways_op) {
          if (!ps %in% cazy_subs_op && ps %in% cazy_map$substrate) {
            # PTS has a pathway that CAZy doesn't cover in this operon
            # -> add PTS locus_tags as transport evidence for the CAZy substrate
          }
        }
        # Cross-assign: give PTS the CAZy substrate if PTS pathway is generic
        for (cs in cazy_subs_op) {
          pts_missing <- pts_detail$locus_tag %in% pts_in_op & pts_detail$pathway != cs
          if (any(pts_missing)) {
            # Duplicate PTS rows with the CAZy substrate as pathway
            extra_pts <- pts_detail[pts_detail$locus_tag %in% pts_in_op, , drop = FALSE]
            extra_pts$pathway <- cs
            extra_pts <- extra_pts[!duplicated(extra_pts[, c("pathway","locus_tag","pts_component")]), , drop = FALSE]
            pts_detail <- dplyr::bind_rows(pts_detail, extra_pts)
          }
        }
        # Also propagate: give transport_data the CAZy substrate
        if (!is.null(transport_data)) {
          tr_in_op <- transport_data$locus_tag[transport_data$locus_tag %in% members]
          for (cs in cazy_subs_op) {
            tr_missing <- tr_in_op[!tr_in_op %in% transport_data$locus_tag[transport_data$pathway == cs]]
            if (length(tr_missing) > 0) {
              extra_tr <- transport_data[transport_data$locus_tag %in% tr_missing, , drop = FALSE]
              extra_tr$pathway <- cs
              transport_data <- dplyr::bind_rows(transport_data, extra_tr[!duplicated(extra_tr$locus_tag), , drop = FALSE])
            }
          }
        }
      }
    }
  }

  # --- Bug 4: Neighbor gene context analysis ---
  # Build a lookup of product annotations for ±5 genes around each PTS locus_tag
  if (is.list(genbank_table) && !is.data.frame(genbank_table) &&
      "features" %in% names(genbank_table))
    genbank_table <- genbank_table$features
  gt <- as.data.frame(genbank_table, stringsAsFactors = FALSE, row.names = NULL)
  neighbor_substrate_hints <- list()  # locus_tag -> vector of substrate hints
  if (!is.null(pts_detail) && nrow(pts_detail) > 0 &&
      all(c("locus_tag","contig","start","product") %in% names(gt))) {
    gt_sorted <- gt[order(gt$contig, as.numeric(gt$start)), , drop = FALSE]
    gt_sorted$row_idx <- seq_len(nrow(gt_sorted))
    pts_lts_all <- unique(pts_detail$locus_tag[!is.na(pts_detail$locus_tag)])
    # Substrate hint patterns from neighbor genes
    neighbor_kw <- list(
      maltose = "maltose|maltos|amylas|alpha.glucos|starch|pullulan",
      cellobiose = "cellobiose|cellulas|beta.glucos|cellulose",
      sucrose = "sucrose|fructos|levan|invert",
      NAG = "N-acetylglucosamin|chitin|GlcNAc|nag[A-Z]",
      galactose = "galactos|galactokinas",
      mannose = "mannos|man[A-Z].*kinas",
      xylose = "xylos|xylan",
      trehalose = "trehalos",
      glucose = "glucose.*kinas|glucokinas|glk",
      fructose = "fructose.*kinas|fructokinas|fruK",
      mannitol = "mannitol|mtl",
      sorbitol = "sorbitol|glucitol|srl",
      lactose = "lactose|beta.galactosid|lac[ZYA]",
      arabinose = "arabinos",
      glucuronate = "glucuronate|pectin",
      glucosamine = "glucosamine",
      ribose = "ribose|ribokinas",
      glycerol = "glycerol|glp[KFDA]"
    )
    for (lt in pts_lts_all) {
      idx <- which(gt_sorted$locus_tag == lt)
      if (length(idx) == 0) next
      idx <- idx[1]
      # ±5 genes
      lo <- max(1, idx - 5)
      hi <- min(nrow(gt_sorted), idx + 5)
      neighbors <- gt_sorted[lo:hi, , drop = FALSE]
      neighbor_products <- paste(neighbors$product[!is.na(neighbors$product)], collapse = " | ")
      hints <- character(0)
      for (sub_name in names(neighbor_kw)) {
        if (grepl(neighbor_kw[[sub_name]], neighbor_products, ignore.case = TRUE)) {
          hints <- c(hints, sub_name)
        }
      }
      if (length(hints) > 0) {
        neighbor_substrate_hints[[lt]] <- unique(hints)
        # Also: add PTS rows for the hinted substrates if not already present
        for (h_sub in hints) {
          existing <- pts_detail$pathway == h_sub & pts_detail$locus_tag == lt
          if (!any(existing)) {
            extra <- pts_detail[pts_detail$locus_tag == lt, , drop = FALSE]
            if (nrow(extra) > 0) {
              extra <- extra[!duplicated(extra$pts_component), , drop = FALSE]
              extra$pathway <- h_sub
              pts_detail <- dplyr::bind_rows(pts_detail, extra)
            }
          }
        }
      }
    }
  }

  # Unique substrates from CAZy + ALL transport (not just those matching CAZy)
  cazy_subs <- if (nrow(cazy_anno) > 0) unique(cazy_anno$substrate) else character(0)
  tr_subs <- if (!is.null(transport_data)) unique(transport_data$pathway) else character(0)
  pts_subs <- if (!is.null(pts_detail)) unique(pts_detail$pathway) else character(0)
  substrates <- sort(unique(c(cazy_subs, tr_subs, pts_subs)))
  if (length(substrates) == 0) return(NULL)

  # Product keyword patterns for evidence scoring
  product_kw <- list(
    maltose = "maltose|starch|amylase|alpha-glucos",
    cellobiose = "cellobiose|cellulase|beta-glucos",
    sucrose = "sucrose|fructos|levansucrase",
    NAG = "N-acetylglucosamine|chitinase|chitin|GlcNAc|nagE|nagP",
    galactose = "galactos",
    mannose = "mannos",
    xylose = "xylan|xylos",
    trehalose = "trehalos",
    glucose = "glucose|glucos",
    peptidoglycan = "peptidoglycan|lysozyme|muramidase",
    glucosamine = "glucosamine",
    arabinose = "arabino",
    glucuronate = "glucuronate|pectin",
    fructose = "fructose|fructokinas",
    lactose = "lactose|beta-galactosid",
    ribose = "ribose|ribokinas",
    mannitol = "mannitol",
    sorbitol = "sorbitol|glucitol",
    rhamnose = "rhamnos",
    fucose = "fucose|fucosid",
    galacturonate = "galacturonate|pectin",
    gluconate = "gluconate",
    glycerol = "glycerol|glycerokinas",
    acetate = "acetate|acetyl-CoA synthet",
    citrate = "citrate",
    succinate = "succinate",
    `L-lactate` = "lactate|L-lactate",
    `D-lactate` = "D-lactate",
    ethanol = "ethanol|alcohol dehydrogen"
  )

  rows <- list()
  for (sub in substrates) {
    # Enzymes for this substrate (may be empty if only transport evidence)
    enz <- if (nrow(cazy_anno) > 0) {
      e <- cazy_anno[cazy_anno$substrate == sub, , drop = FALSE]
      e[!duplicated(e$locus_tag), , drop = FALSE]
    } else {
      data.frame(locus_tag = character(0), family = character(0),
                 enzyme_label = character(0), gene = character(0),
                 gene_product = character(0), substrate = character(0),
                 product = character(0), stringsAsFactors = FALSE, row.names = NULL)
    }

    # Transporters for this substrate
    tr_rows <- NULL
    if (!is.null(transport_data)) {
      tr_rows <- transport_data[transport_data$pathway == sub, , drop = FALSE]
      tr_rows <- tr_rows[!is.na(tr_rows$locus_tag) & nzchar(tr_rows$locus_tag), , drop = FALSE]
      tr_rows <- tr_rows[!duplicated(tr_rows$locus_tag), , drop = FALSE]
    }

    # PTS detail for this substrate
    pts_rows <- NULL
    if (!is.null(pts_detail)) {
      pts_rows <- pts_detail[pts_detail$pathway == sub, , drop = FALSE]
      pts_rows <- pts_rows[!is.na(pts_rows$locus_tag), , drop = FALSE]
      pts_rows <- pts_rows[!duplicated(pts_rows[, c("locus_tag","pts_component")]), , drop = FALSE]
    }

    # Bug 2: Track fused PTS components per locus_tag
    # e.g., nagEcba -> one locus_tag covers EIIC+EIIB+EIIA
    pts_fused <- NULL
    if (!is.null(pts_rows) && nrow(pts_rows) > 0) {
      pts_fused <- pts_rows |>
        dplyr::group_by(.data$locus_tag) |>
        dplyr::summarise(
          components = paste(unique(.data$pts_component), collapse = "|"),
          n_components = dplyr::n_distinct(.data$pts_component),
          is_fused = dplyr::n_distinct(.data$pts_component) > 1,
          .groups = "drop"
        )
    }

    # Bug 1: Include PTS locus_tags in operon matching
    all_lt <- c(
      enz$locus_tag,
      if (!is.null(tr_rows)) tr_rows$locus_tag else NULL,
      if (!is.null(pts_rows)) pts_rows$locus_tag else NULL
    )
    all_lt <- unique(all_lt[!is.na(all_lt) & nzchar(all_lt)])
    operon_ids <- character(0)
    if (!is.null(operon_info) && length(all_lt) > 0) {
      operon_ids <- unique(operon_info$lt_map[all_lt[all_lt %in% names(operon_info$lt_map)]])
      operon_ids <- operon_ids[!is.na(operon_ids)]
    }

    # Bug 4: Check if any PTS locus_tags have neighbor-based substrate hints
    has_neighbor_evidence <- FALSE
    if (length(neighbor_substrate_hints) > 0 && !is.null(pts_rows) && nrow(pts_rows) > 0) {
      for (lt in unique(pts_rows$locus_tag)) {
        if (lt %in% names(neighbor_substrate_hints) && sub %in% neighbor_substrate_hints[[lt]]) {
          has_neighbor_evidence <- TRUE
          break
        }
      }
    }

    # Evidence scoring
    if (nrow(enz) > 0) {
      enz$product_evidence <- FALSE
      kw_pattern <- product_kw[[sub]]
      if (!is.null(kw_pattern)) {
        enz$product_evidence <- grepl(kw_pattern, enz$gene_product, ignore.case = TRUE)
      }

      # Evidence score per enzyme
      enz$evidence_n <- 1L  # Always have CAZy family match
      enz$evidence_n <- enz$evidence_n +
        as.integer(!is.null(tr_rows) && nrow(tr_rows) > 0) +      # GapMind transport match
        as.integer(length(operon_ids) > 0) +                        # Operon co-localization
        as.integer(enz$product_evidence) +                          # Product keyword
        as.integer(has_neighbor_evidence)                           # Bug 4: neighbor gene context
      enz$evidence_conf <- dplyr::case_when(
        enz$evidence_n >= 4 ~ "very_high", enz$evidence_n == 3 ~ "very_high",
        enz$evidence_n == 2 ~ "high",
        enz$evidence_n == 1 ~ "medium", TRUE ~ "low")
    }

    products <- if (nrow(enz) > 0) unique(enz$product) else sub

    rows[[sub]] <- list(
      substrate = sub,
      enzymes = enz,
      transporters = tr_rows,
      pts_components = pts_rows,
      pts_fused = pts_fused,
      products = products,
      operon_ids = operon_ids,
      has_pts = !is.null(pts_rows) && nrow(pts_rows) > 0,
      has_neighbor_evidence = has_neighbor_evidence
    )
  }
  rows
}


# === Complex carbohydrate degradation cascade definitions ===
.dnmb_cct_complex_carb_cascades <- function() {
  s <- function(substrate, product, gh_families, cleavage_type, bond_desc) {
    list(substrate = substrate, product = product, gh_families = gh_families,
         cleavage_type = cleavage_type, bond_desc = bond_desc)
  }
  list(
    starch = list(
      id = "starch", polymer = "Starch / Pullulan (\u03b1-glucan)",
      monomer = "glucose", entry = "G6P",
      chain_color = "#E6550D", unit_label = "Glc",
      unit_snfg = c("glucose"),
      bond_label = "alpha-1,4",
      steps = list(
        s("starch", "maltodextrin", c("GH13","GH13_5","GH13_14","GH49","GH57"), "endo",
          "alpha-1,4/1,6-glucan cleavage (incl. pullulanase)"),
        s("maltodextrin", "maltose", c("GH13_20","GH13_39","GH13_45","GH15"), "exo",
          "alpha-1,4-glucan non-reducing end"),
        s("maltose", "glucose", c("GH13","GH13_31","GH4","GH65"), "exo",
          "alpha-1,4-glucoside hydrolysis"))),
    cellulose = list(
      id = "cellulose", polymer = "Cellulose",
      monomer = "glucose", entry = "G6P",
      chain_color = "#D94801", unit_label = "Glc",
      unit_snfg = c("glucose"),
      bond_label = "beta-1,4",
      steps = list(
        s("cellulose", "cellodextrin", c("GH5","GH9","GH44","GH48"), "endo",
          "beta-1,4-glucan internal cleavage"),
        s("cellodextrin", "cellobiose", c("GH5","GH48"), "exo",
          "beta-1,4-glucan reducing/non-reducing end"),
        s("cellobiose", "glucose", c("GH1","GH3"), "exo",
          "beta-1,4-glucoside hydrolysis"))),
    chitin = list(
      id = "chitin", polymer = "Chitin",
      monomer = "GlcNAc", entry = "F6P",
      chain_color = "#3182BD", unit_label = "GlcNAc",
      unit_snfg = c("GlcNAc"),
      bond_label = "beta-1,4",
      steps = list(
        s("chitin", "chitooligosaccharide", c("GH18","GH19"), "endo",
          "beta-1,4-GlcNAc internal cleavage"),
        s("chitooligosaccharide", "GlcNAc", c("GH20","GH3"), "exo",
          "beta-1,4-GlcNAc terminal hydrolysis"))),
    xylan = list(
      id = "xylan", polymer = "Xylan",
      monomer = "xylose", entry = "GAP",
      chain_color = "#FF8C00", unit_label = "Xyl",
      unit_snfg = c("xylose"),
      bond_label = "beta-1,4",
      steps = list(
        s("xylan", "xylooligosaccharide", c("GH10","GH11"), "endo",
          "beta-1,4-xylan internal cleavage"),
        s("xylooligosaccharide", "xylose", c("GH43","GH3"), "exo",
          "beta-1,4-xyloside hydrolysis"))),
    fructan = list(
      id = "fructan", polymer = "Inulin/Levan (fructan)",
      monomer = "fructose", entry = "F6P",
      chain_color = "#2CA02C", unit_label = "Fru",
      unit_snfg = c("fructose"),
      bond_label = "beta-2,1",
      steps = list(
        s("fructan", "fructooligosaccharide", c("GH32","GH68"), "endo",
          "beta-2,1/2,6-fructan internal cleavage"),
        s("fructooligosaccharide", "fructose", c("GH32"), "exo",
          "beta-fructoside hydrolysis"))),
    mannan = list(
      id = "mannan", polymer = "beta-Mannan",
      monomer = "mannose", entry = "F6P",
      chain_color = "#D6604D", unit_label = "Man",
      unit_snfg = c("mannose"),
      bond_label = "beta-1,4",
      steps = list(
        s("mannan", "manno-oligosaccharide", c("GH26","GH113"), "endo",
          "beta-1,4-mannan internal cleavage"),
        s("manno-oligosaccharide", "mannose", c("GH130","GH130_4","GH92"), "exo",
          "beta-mannoside hydrolysis"))),
    galactan = list(
      id = "galactan", polymer = "Galactan/Raffinose",
      monomer = "galactose", entry = "G6P",
      chain_color = "#9467BD", unit_label = "Gal",
      unit_snfg = c("galactose"),
      bond_label = "alpha-1,6",
      steps = list(
        s("galactan", "galactose", c("GH36","GH42","GH2"), "exo",
          "alpha/beta-galactoside hydrolysis"))),
    peptidoglycan = list(
      id = "peptidoglycan", polymer = "Peptidoglycan",
      monomer = "MurNAc+GlcNAc", entry = "F6P",
      chain_color = "#8C564B", unit_label = "MurNAc",
      unit_snfg = c("MurNAc", "GlcNAc"),
      bond_label = "beta-1,4",
      steps = list(
        s("peptidoglycan", "muropeptide", c("GH23","GH24","GH25","GH73"), "endo",
          "beta-1,4-MurNAc-GlcNAc cleavage"),
        s("muropeptide", "MurNAc+GlcNAc", c("GH18","GH20"), "exo",
          "N-acetylhexosaminidase hydrolysis")))
  )
}

# === Sugar chain drawing with SNFG symbols and GH scissor marks ===
.dnmb_cct_sugar_chain_layers <- function(x_start, y, n_units, unit_size = 0.35,
                                          cut_positions, gh_labels,
                                          cleavage_types, chain_color,
                                          unit_label = NULL,
                                          unit_snfg = NULL,
                                          bond_label = NULL) {
  layers <- list()
  spacing <- unit_size * 1.4

  # Resolve SNFG sugar types for each chain unit (cycling for heteropolymers)
  if (is.null(unit_snfg)) unit_snfg <- "Glc"
  snfg_types <- rep_len(unit_snfg, n_units)
  snfg_info  <- .dnmb_snfg_lookup(snfg_types)
  # Resolve abbreviations: canonical names are already short
  abbrs <- vapply(snfg_types, function(s) {
    cn <- .dnmb_snfg_resolve_synonym(s)
    if (nchar(cn) <= 6) cn else substr(cn, 1, 3)
  }, character(1))

  # Unit center positions
  ux <- x_start + (seq_len(n_units) - 1) * spacing
  uy <- rep(y, n_units)

  # Bond lines FIRST (center-to-center, no gap — symbols paint over ends)
  # Determine linetype from bond_label: beta = dashed, alpha = solid
  if (n_units > 1) {
    is_beta <- !is.null(bond_label) && grepl("\u03b2", bond_label)
    bond_lty <- if (is_beta) "dashed" else "solid"
    for (bi in seq_len(n_units - 1)) {
      seg_df <- data.frame(x = ux[bi], xend = ux[bi + 1],
                           y = y, yend = y, stringsAsFactors = FALSE, row.names = NULL)
      layers[[length(layers) + 1]] <- ggplot2::geom_segment(
        data = seg_df,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        color = "#888888", linewidth = 0.8, linetype = bond_lty,
        lineend = "round", inherit.aes = FALSE)
      # Bond label above each bond segment
      if (!is.null(bond_label) && nchar(bond_label) > 0) {
        mid_x <- (ux[bi] + ux[bi + 1]) / 2
        layers[[length(layers) + 1]] <- ggplot2::annotate(
          "text", x = mid_x, y = y + unit_size * 0.65,
          label = bond_label, size = 1.3, color = "grey40", fontface = "italic")
      }
    }
  }

  # Draw SNFG symbols ON TOP of bond lines
  for (i in seq_len(n_units)) {
    col_i <- snfg_info$snfg_color[i]
    sym_layers <- .dnmb_snfg_render_symbol(ux[i], uy[i], snfg_types[i],
      r = 0.27, point_size = 7, active = TRUE)
    for (sl in sym_layers) layers[[length(layers) + 1]] <- sl
    layers[[length(layers) + 1]] <- ggplot2::annotate(
      "text", x = ux[i], y = uy[i] - unit_size * 0.7,
      label = abbrs[i], size = 1.1, color = colorspace::darken(col_i, 0.2),
      fontface = "bold")
  }

  # GH cleavage scissor marks at cut positions
  if (length(cut_positions) > 0) {
    for (ci in seq_along(cut_positions)) {
      pos <- cut_positions[ci]
      if (pos < 1 || pos >= n_units) next
      cx <- (ux[pos] + ux[pos + 1]) / 2
      is_endo <- grepl("endo", cleavage_types[ci], ignore.case = TRUE)
      sc_col <- if (is_endo) "#D32F2F" else "#1565C0"
      sc_half <- unit_size * 0.4
      sc_df <- data.frame(
        x = c(cx - sc_half * 0.4, cx + sc_half * 0.4),
        xend = c(cx + sc_half * 0.4, cx - sc_half * 0.4),
        y = c(y - sc_half, y - sc_half),
        yend = c(y + sc_half, y + sc_half), stringsAsFactors = FALSE, row.names = NULL)
      layers[[length(layers) + 1]] <- ggplot2::geom_segment(
        data = sc_df,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        color = sc_col, linewidth = 0.6, lineend = "round")
      ctype_lab <- if (is_endo) "(endo)" else "(exo)"
      sc_lab_df <- data.frame(
        x = cx, y = y + sc_half + unit_size * 0.4,
        label = paste0(gh_labels[ci], "\n", ctype_lab), stringsAsFactors = FALSE, row.names = NULL)
      layers[[length(layers) + 1]] <- ggplot2::geom_text(
        data = sc_lab_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = 1.2, color = sc_col, fontface = "bold", lineheight = 0.8)
    }
  }

  # Return chain width for auto-placement calculations
  chain_width <- (n_units - 1) * spacing + unit_size
  attr(layers, "chain_width") <- chain_width
  attr(layers, "unit_xs") <- ux
  layers
}

# === Cell shape layers: cytoplasm/extracellular zones + membrane bilayer ===
.dnmb_cct_cell_shape_layers <- function(xmin, xmax, ymin, ymax,
                                          membrane_x, membrane_width = 0.4) {
  layers <- list()
  # No background fill rects — clean white background
  mem_inner <- membrane_x - membrane_width / 4
  mem_outer <- membrane_x + membrane_width / 4
  mem_seg <- data.frame(
    x = c(mem_inner, mem_outer), xend = c(mem_inner, mem_outer),
    y = c(ymin, ymin), yend = c(ymax, ymax), stringsAsFactors = FALSE, row.names = NULL)
  layers[[length(layers) + 1]] <- ggplot2::geom_segment(
    data = mem_seg,
    ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    color = "#D4A574", linewidth = 1.0, alpha = 0.8)
  n_dots <- max(3, round((ymax - ymin) / 0.8))
  dot_ys <- seq(ymin + 0.3, ymax - 0.3, length.out = n_dots)
  head_df <- data.frame(
    x = c(rep(mem_inner, n_dots), rep(mem_outer, n_dots)),
    y = c(dot_ys, dot_ys), stringsAsFactors = FALSE, row.names = NULL)
  layers[[length(layers) + 1]] <- ggplot2::geom_point(
    data = head_df, ggplot2::aes(x = .data$x, y = .data$y),
    shape = 21, size = 1.5, fill = "#FFCC80", color = "#FFCC80", stroke = 0)
  layers[[length(layers) + 1]] <- ggplot2::annotate("text",
    x = (xmin + membrane_x) / 2, y = ymax + 0.4,
    label = "Cytoplasm", size = 2.8, fontface = "bold.italic", color = "#546E7A")
  layers[[length(layers) + 1]] <- ggplot2::annotate("text",
    x = (membrane_x + xmax) / 2, y = ymax + 0.4,
    label = "Extracellular", size = 2.8, fontface = "bold.italic", color = "#795548")
  layers[[length(layers) + 1]] <- ggplot2::annotate("text",
    x = membrane_x, y = ymax + 0.4,
    label = "Membrane", size = 2.5, fontface = "bold", color = "#8D6E63")
  layers
}

# === Match active cascades against genome's CAZy hits ===
.dnmb_cct_match_cascades <- function(cazy_hits, cascades) {
  if (is.null(cazy_hits) || nrow(cazy_hits) == 0) return(list())
  hit_families <- unique(cazy_hits$family)
  active <- list()
  for (cname in names(cascades)) {
    casc <- cascades[[cname]]
    step_hits <- list()
    any_hit <- FALSE
    for (si in seq_along(casc$steps)) {
      step <- casc$steps[[si]]
      matched <- intersect(step$gh_families, hit_families)
      if (length(matched) > 0) {
        any_hit <- TRUE
        lt_rows <- cazy_hits[cazy_hits$family %in% matched, , drop = FALSE]
        step_hits[[si]] <- list(
          matched_families = matched, locus_tags = unique(lt_rows$locus_tag),
          cleavage_type = step$cleavage_type, bond_desc = step$bond_desc,
          substrate = step$substrate, product = step$product)
      } else {
        step_hits[[si]] <- list(
          matched_families = character(0), locus_tags = character(0),
          cleavage_type = step$cleavage_type, bond_desc = step$bond_desc,
          substrate = step$substrate, product = step$product)
      }
    }
    if (any_hit) active[[cname]] <- list(cascade = casc, step_hits = step_hits)
  }
  active
}

# === Main plotting function (v4: cell shape + complex carb cascades) ===
.dnmb_plot_cazy_carbon_transport_map <- function(genbank_table, output_dir,
                                                  file_stub = "CAZy_overview") {
  # Direct call to 3zone version (no fallback)
  return(.dnmb_plot_cazy_carbon_3zone_v2(genbank_table, output_dir = output_dir,
                                          file_stub = file_stub))
}

# Legacy pane 5 function kept for reference but not called
.dnmb_plot_cazy_carbon_transport_map_legacy <- function(genbank_table, output_dir,
                                                  file_stub = "CAZy_overview") {
  # Guard: accept raw genbank parse list (extract features data.frame)
  if (is.list(genbank_table) && !is.data.frame(genbank_table) &&
      "features" %in% names(genbank_table))
    genbank_table <- genbank_table$features

  # ---- Data extraction ----
  cazy_map   <- .dnmb_cct_cazy_substrate_map()
  cazy_hits  <- .dnmb_cct_cazy_extract(genbank_table, output_dir)
  if (is.null(cazy_hits)) cazy_hits <- data.frame(
    locus_tag = character(0), family = character(0), stringsAsFactors = FALSE, row.names = NULL)
  transport  <- .dnmb_cct_transport_extract(output_dir)
  pts_detail <- .dnmb_cct_pts_detail(output_dir)
  operon_info <- .dnmb_cct_infer_operons(genbank_table, cazy_hits, transport)
  rows <- .dnmb_cct_build_rows(cazy_hits, transport, pts_detail, operon_info,
                                cazy_map, genbank_table)
  if (is.null(rows)) rows <- list()  # empty is OK

  # ---- v4: Complex carbohydrate cascade detection ----
  cascades <- .dnmb_cct_complex_carb_cascades()
  active_cascades <- .dnmb_cct_match_cascades(cazy_hits, cascades)
  n_cascades <- length(active_cascades)

  # ---- Layout v6: vertical zones, better proportions ----
  # Y-axis = depth (top->bottom): substrate -> membrane -> enzymes -> glycolysis/TCA
  # X-axis = substrate spread (each substrate = one column)
  #
  # Zone Y coordinates:
  y_substrate   <- 20.0   # top: substrate nodes (extracellular)
  y_cascade_top <- 21.5   # polysaccharide cascades above substrates
  y_transport   <- 18.0   # transporter segments
  y_membrane    <- 17.0   # membrane line (horizontal)
  y_enz_start   <- 16.0   # enzyme segments (cytoplasmic side)
  y_enz_step    <- 0.7    # vertical spacing per enzyme
  y_product     <- 13.0   # degradation product nodes
  y_backbone    <- 10.0   # glycolysis backbone (horizontal, bottom zone)
  y_tca         <- 6.0    # TCA cycle center (below backbone)
  y_ppp         <- 11.5   # PPP branch (parallel to backbone)

  # Glycolysis backbone runs LEFT->RIGHT (horizontal)
  bx <- 1.0   # legacy compatibility (used for entry_map lookups)

  # TCA cycle: perfect circle, 8 nodes at 45 deg
  tca_cx <- 12.0; tca_cy <- y_tca; tca_r <- 1.8
  tca_n <- 8L
  # Start from top (OAA receives AcCoA), go clockwise
  tca_ang <- pi/2 - (0:(tca_n - 1)) * (2 * pi / tca_n)
  tca_x <- tca_cx + tca_r * cos(tca_ang)
  tca_y_pts <- tca_cy + tca_r * sin(tca_ang)
  tca_ids    <- c("OAA","Citrate","Isocit","AKG","SucCoA","Succinate","Fumarate","Malate")
  tca_labels <- c("OAA","Citrate","Isocitrate","alpha-KG","Succinyl-CoA","Succinate","Fumarate","Malate")
  tca_arrow_color <- c("#66BB6A","#43A047","#388E3C","#2E7D32",
                        "#558B2F","#689F38","#7CB342","#8BC34A")

  # Glycolysis backbone: HORIZONTAL (left->right) at y_backbone
  # Scale to span full substrate width (calculated after sub_x assignment)
  bb_x_start <- 0.0
  bb_x_step  <- 2.8  # provisional; recalculated after sub_x
  bb_ids <- c("Glucose","G6P","F6P","FBP","DHAP","GAP","3PG","PEP","Pyruvate","AcCoA")
  bb_labels <- c("Glucose","G6P","F6P","FBP","DHAP","GAP","3-PG","PEP","Pyruvate","AcCoA")
  bb_nodes <- data.frame(
    id    = bb_ids,
    x     = bb_x_start + (seq_along(bb_ids) - 1) * bb_x_step,
    y     = rep(y_backbone, length(bb_ids)),
    label = bb_labels,
    stringsAsFactors = FALSE, row.names = NULL)
  # DHAP branches down from FBP
  bb_nodes$y[bb_nodes$id == "DHAP"] <- y_backbone - 0.8
  bb_nodes$x[bb_nodes$id == "DHAP"] <- bb_nodes$x[bb_nodes$id == "FBP"] + 0.0

  # Color gradient along glycolysis: warm (left) -> cool (right)
  bb_gradient <- c("#E57373","#EF5350","#F44336","#E53935","#D32F2F",
                    "#C62828","#B71C1C","#880E4F","#6A1B9A")
  # Cofactor annotations at key glycolysis steps
  bb_cofactors <- list(
    list(from = "Glucose", to = "G6P",  label = "ATP -> ADP"),
    list(from = "F6P",     to = "FBP",  label = "ATP -> ADP"),
    list(from = "GAP",     to = "3PG",  label = "NAD+ -> NADH\n+ ATP"),
    list(from = "3PG",     to = "PEP",  label = "H2O"),
    list(from = "PEP",     to = "Pyruvate", label = "ADP -> ATP")
  )

  # PPP branch: runs parallel above backbone
  g6p_x <- bb_nodes$x[bb_nodes$id == "G6P"]
  ppp_nodes <- data.frame(
    id    = c("6PGL","6PG","Ru5P","R5P","Xu5P"),
    x     = g6p_x + c(1.4, 2.8, 4.2, 5.2, 5.2),
    y     = rep(y_ppp, 5),
    label = c("6PGL","6PG","Ru5P","R5P","Xu5P"),
    stringsAsFactors = FALSE, row.names = NULL)
  # R5P and Xu5P fan out
  ppp_nodes$y[ppp_nodes$id == "R5P"]  <- y_ppp + 0.5
  ppp_nodes$y[ppp_nodes$id == "Xu5P"] <- y_ppp - 0.5

  # TCA cycle nodes
  tca_nodes <- data.frame(
    id = tca_ids, x = tca_x, y = tca_y_pts, label = tca_labels,
    stringsAsFactors = FALSE, row.names = NULL)

  # Substrate -> backbone entry point mapping
  # Entry points reference the node IDs in bb_nodes, ppp_nodes, or tca_nodes
  entry_map <- c(
    # -> G6P (hexoses via glucose-6-phosphate)
    glucose = "G6P", maltose = "G6P", cellobiose = "G6P",
    trehalose = "G6P", galactose = "G6P", lactose = "G6P",
    # -> F6P (fructose-6-phosphate entry)
    fructose = "F6P", sucrose = "F6P", mannose = "F6P",
    NAG = "F6P", glucosamine = "F6P", peptidoglycan = "F6P",
    mannitol = "F6P", sorbitol = "F6P",
    # -> Xu5P/Ru5P via PPP
    xylose = "Xu5P", arabinose = "Xu5P", ribose = "R5P",
    # -> GAP (Entner-Doudoroff)
    glucuronate = "GAP", galacturonate = "GAP", gluconate = "GAP",
    rhamnose = "GAP", glycerol = "DHAP",
    # -> PEP
    fucose = "PEP", myoinositol = "PEP",
    # -> Pyruvate
    `L-lactate` = "Pyruvate", `D-lactate` = "Pyruvate", ethanol = "Pyruvate",
    # -> AcCoA
    acetate = "AcCoA",
    # -> TCA cycle nodes directly
    citrate = "Citrate", succinate = "Succinate", fumarate = "Fumarate",
    `L-malate` = "Malate"
  )

  # ---- Build node coordinate lookups ----
  all_nodes <- rbind(bb_nodes, ppp_nodes, tca_nodes)
  node_y <- setNames(all_nodes$y, all_nodes$id)
  node_x <- setNames(all_nodes$x, all_nodes$id)

  # ---- Show ALL substrates from entry_map ----
  all_substrates <- names(entry_map)
  has_evidence <- all_substrates %in% names(rows)
  substrates <- all_substrates
  sub_entry <- entry_map[substrates]

  # Group substrates by entry point, evidence-positive first
  unique_entries <- unique(sub_entry)
  entry_groups <- list()
  for (eid in unique_entries) {
    subs <- substrates[sub_entry == eid]
    subs_pos <- subs[subs %in% names(rows)]
    subs_neg <- subs[!subs %in% names(rows)]
    entry_groups[[eid]] <- c(subs_pos, subs_neg)
  }

  # ---- Assign substrate X (horizontal spread, grouped by entry) ----
  # Differentiated spacing: enzyme evidence > transport-only > no-evidence
  sub_x <- numeric(length(substrates)); names(sub_x) <- substrates
  cursor_x <- 0.0
  for (eid in names(entry_groups)) {
    subs <- entry_groups[[eid]]
    for (s in subs) {
      sub_x[s] <- cursor_x
      has_d <- s %in% names(rows)
      has_enz <- has_d && nrow(rows[[s]]$enzymes) > 0
      cursor_x <- cursor_x + if (has_enz) 1.6 else if (has_d) 1.1 else 0.5
    }
    cursor_x <- cursor_x + 0.3  # gap between entry groups
  }
  # sub_y is FIXED for all substrates (extracellular zone)
  sub_y <- setNames(rep(y_substrate, length(substrates)), substrates)

  # ---- Rescale backbone/TCA/PPP — compact, centered under substrates ----
  sub_x_range <- max(sub_x) - min(sub_x)
  sub_cx <- (min(sub_x) + max(sub_x)) / 2  # substrate center
  # Fixed backbone step (compact, readable) centered under substrates
  bb_x_step <- 2.5
  bb_total_w <- bb_x_step * (length(bb_ids) - 1)
  bb_x_start <- sub_cx - bb_total_w / 2
  bb_nodes$x <- bb_x_start + (seq_along(bb_ids) - 1) * bb_x_step
  bb_nodes$x[bb_nodes$id == "DHAP"] <- bb_nodes$x[bb_nodes$id == "FBP"]
  # TCA cycle: right of backbone (after AcCoA), shifted right by radius
  acoa_x <- bb_nodes$x[bb_nodes$id == "AcCoA"]
  tca_cx <- acoa_x + tca_r + 1.0
  tca_x <- tca_cx + tca_r * cos(tca_ang)
  tca_y_pts <- tca_cy + tca_r * sin(tca_ang)
  tca_nodes$x <- tca_x; tca_nodes$y <- tca_y_pts
  # PPP branch relative to G6P
  g6p_x <- bb_nodes$x[bb_nodes$id == "G6P"]
  ppp_step <- 1.4
  ppp_nodes$x <- g6p_x + c(1, 2, 3, 3.8, 3.8) * ppp_step
  ppp_nodes$y <- rep(y_ppp, 5)
  ppp_nodes$y[ppp_nodes$id == "R5P"]  <- y_ppp + 0.5
  ppp_nodes$y[ppp_nodes$id == "Xu5P"] <- y_ppp - 0.5
  # Refresh coordinate lookups
  all_nodes <- rbind(bb_nodes, ppp_nodes, tca_nodes)
  node_y <- setNames(all_nodes$y, all_nodes$id)
  node_x <- setNames(all_nodes$x, all_nodes$id)

  # ---- Y positions for vertical pathway chain (per substrate column) ----
  # Each substrate column goes: substrate(top) -> transport -> membrane -> enzymes -> product(bottom)
  # Actual positions set by zone constants above

  # ---- Confidence colors (same as gapmind_aa) ----
  line_col <- c(very_high = "#1A9641", high = "#2CA25F", medium = "#FEC44F",
                low = "#F03B20", none = "#CCCCCC")
  conf_levels <- c("very_high", "high", "medium", "low", "none")

  # ---- Pathway ribbon colors (unique per substrate, grouped by entry) ----
  pw_colors <- c(
    # G6P group
    glucose = "#E6550D", maltose = "#B8860B", cellobiose = "#D94801",
    trehalose = "#8C6BB1", galactose = "#9467BD", lactose = "#C49C94",
    # F6P group
    fructose = "#2CA02C", sucrose = "#31A354", mannose = "#D6604D",
    NAG = "#3182BD", glucosamine = "#E7298A", peptidoglycan = "#8C564B",
    mannitol = "#17BECF", sorbitol = "#9EDAE5",
    # GAP group
    xylose = "#FF8C00", arabinose = "#FD8D3C", ribose = "#FFBB78",
    glucuronate = "#636363", galacturonate = "#969696", gluconate = "#BCBD22",
    rhamnose = "#AEC7E8", glycerol = "#7F7F7F",
    # PEP group
    fucose = "#9C9EDE", myoinositol = "#CE93D8",
    # Pyruvate group
    `L-lactate` = "#6BAED6", `D-lactate` = "#4292C6", ethanol = "#2171B5",
    # AcCoA group
    acetate = "#756BB1",
    # TCA group
    citrate = "#006D2C", succinate = "#41AB5D", fumarate = "#74C476",
    `L-malate` = "#A1D99B"
  )

  # ---- Build segment data for all pathways (v5: vertical columns) ----
  seg_list <- list()    # enzyme/transport segments with confidence coloring
  ribbon_list <- list() # background ribbon gradients
  label_list <- list()  # gene name + locus_tag labels
  sub_nodes <- list()   # substrate nodes (top)
  prod_nodes <- list()  # degradation product nodes (bottom)
  pts_cascade_list <- list() # PTS cascade mini-segments
  pts_label_list <- list()   # summarized PTS transporter labels
  operon_boxes <- list()

  for (si in seq_along(substrates)) {
    sub <- substrates[si]
    sx <- sub_x[sub]  # column X for this substrate
    pw_col <- unname(pw_colors[sub])
    if (is.na(pw_col) || is.null(pw_col)) pw_col <- "#999999"
    entry_id <- entry_map[sub]; if (is.na(entry_id)) entry_id <- "G6P"
    has_data <- sub %in% names(rows)

    # --- Substrate node (top, extracellular) --- always shown
    sub_nodes[[length(sub_nodes) + 1]] <- data.frame(
      x = sx, y = y_substrate, label = sub, has_data = has_data,
      stringsAsFactors = FALSE, row.names = NULL)

    if (!has_data) {
      # ---- NO EVIDENCE: short grey vertical connector ----
      ribbon_list[[length(ribbon_list) + 1]] <- data.frame(
        x = sx, y = y_substrate, xend = sx, yend = y_membrane + 0.5,
        color = "#E0E0E0", stringsAsFactors = FALSE, row.names = NULL)
      next
    }

    # ---- HAS EVIDENCE: full vertical pathway chain ----
    r <- rows[[sub]]

    # --- Product node (near backbone) ---
    prods <- r$products
    n_prd <- length(prods)
    prd_xs <- if (n_prd <= 1) sx else seq(sx - 0.2 * (n_prd - 1), sx + 0.2 * (n_prd - 1), length.out = n_prd)
    entry_x <- if (entry_id %in% names(node_x)) node_x[entry_id] else bb_nodes$x[1]
    entry_yy <- if (entry_id %in% names(node_y)) node_y[entry_id] else y_backbone
    for (j in seq_len(n_prd)) {
      prod_nodes[[length(prod_nodes) + 1]] <- data.frame(
        x = prd_xs[j], y = y_product, label = prods[j], stringsAsFactors = FALSE, row.names = NULL)
      # Simple diagonal from product to entry point (thin, transparent)
      ribbon_list[[length(ribbon_list) + 1]] <- data.frame(
        x = prd_xs[j], y = y_product, xend = entry_x, yend = entry_yy,
        color = pw_col, stringsAsFactors = FALSE, row.names = NULL)
    }

    # --- CAZy enzyme segments (vertical, between membrane and product) ---
    enz <- r$enzymes
    n_enz <- nrow(enz)
    if (n_enz > 4) {
      enz <- enz[order(-enz$evidence_n), , drop = FALSE][1:4, , drop = FALSE]
      n_enz <- 4
    }
    enz_ys <- if (n_enz == 0) numeric(0) else {
      seq(y_enz_start, by = -y_enz_step, length.out = n_enz)
    }
    for (j in seq_len(n_enz)) {
      ev_conf <- enz$evidence_conf[j]
      ev_col <- unname(line_col[ev_conf])
      if (is.na(ev_col)) ev_col <- "#CCCCCC"
      y1 <- if (j == 1) y_membrane else enz_ys[j - 1]
      y2 <- enz_ys[j]
      seg_list[[length(seg_list) + 1]] <- data.frame(
        x = sx, y = y1, xend = sx, yend = y2,
        confidence = ev_conf, color = ev_col,
        stringsAsFactors = FALSE, row.names = NULL)
      # Gene label — only show family name (compact), max 2 labels
      if (j <= 2) {
        gene_lab <- enz$enzyme_label[j]
        label_list[[length(label_list) + 1]] <- data.frame(
          x = sx + 0.4, y = (y1 + y2) / 2,
          label = gene_lab,
          color = ifelse(ev_conf == "none", "#AAAAAA", "#333333"),
          stringsAsFactors = FALSE, row.names = NULL)
      }
    }

    # --- Transport segments (vertical, above membrane) ---
    tr_rows <- r$transporters
    pts_comp <- r$pts_components
    has_pts <- r$has_pts && !is.null(pts_comp) && nrow(pts_comp) > 0

    tr_items <- NULL
    if (!is.null(tr_rows) && nrow(tr_rows) > 0) {
      tr_items <- tr_rows[!duplicated(tr_rows$locus_tag), , drop = FALSE]
    }

    # Enzyme bottom -> product connection
    enz_end_y <- if (n_enz > 0) min(enz_ys) else y_membrane
    seg_list[[length(seg_list) + 1]] <- data.frame(
      x = sx, y = enz_end_y, xend = sx, yend = y_product,
      confidence = "none", color = "#D5D5D5", stringsAsFactors = FALSE, row.names = NULL)

    # PTS cascade — simplified: one colored bar at membrane + summary label
    if (has_pts) {
      # Single PTS indicator bar at membrane
      pts_cascade_list[[length(pts_cascade_list) + 1]] <- data.frame(
        x = sx - 0.5, y = y_membrane,
        xend = sx + 0.5, yend = y_membrane,
        comp = "PTS", lt = "PTS",
        is_fused = FALSE, stringsAsFactors = FALSE, row.names = NULL)
      # Summary label
      n_pts_comp <- length(unique(pts_comp$pts_component))
      pts_label_list[[length(pts_label_list) + 1]] <- data.frame(
        x = sx,
        y = y_membrane + 0.55,
        label = paste0("PTS(", n_pts_comp, ")"),
        color = "#4E342E",
        stringsAsFactors = FALSE
        )
    }

    # Non-PTS transporters — single vertical bar + count label
    nonpts <- if (!is.null(tr_items)) tr_items[!tr_items$is_pts, , drop = FALSE] else NULL
    n_nonpts <- if (!is.null(nonpts)) nrow(nonpts) else 0L
    if (n_nonpts > 0) {
      best_conf <- nonpts$confidence[which.max(match(nonpts$confidence, c("high","medium","low","none")))]
      best_col <- unname(line_col[best_conf])
      if (is.na(best_col)) best_col <- "#CCCCCC"
      # One vertical bar crossing membrane
      seg_list[[length(seg_list) + 1]] <- data.frame(
        x = sx, y = y_membrane - 0.4,
        xend = sx, yend = y_membrane + 0.4,
        confidence = best_conf, color = best_col, stringsAsFactors = FALSE, row.names = NULL)
      # Best gene name + count
      best_gene <- nonpts$gene_name[1]
      if (is.na(best_gene) || !nzchar(best_gene)) best_gene <- nonpts$step[1]
      tr_label <- if (n_nonpts == 1) best_gene else paste0(best_gene, " +", n_nonpts - 1)
      label_list[[length(label_list) + 1]] <- data.frame(
        x = sx + 0.45, y = y_membrane - 0.3,
        label = tr_label,
        color = "#555555",
        stringsAsFactors = FALSE, row.names = NULL)
    }

    # Substrate -> membrane ribbon (vertical, colored)
    ribbon_list[[length(ribbon_list) + 1]] <- data.frame(
      x = sx, y = y_substrate, xend = sx, yend = y_membrane + 0.5,
      color = pw_col, stringsAsFactors = FALSE, row.names = NULL)

    # Operon brackets (enzyme section)
    if (length(r$operon_ids) > 0 && !is.null(operon_info)) {
      all_op_ys <- enz_ys
      all_op_lts <- enz$locus_tag
      for (oid in r$operon_ids) {
        members <- operon_info$groups[[oid]]
        in_op <- which(all_op_lts %in% members)
        if (length(in_op) >= 2) {
          op_ys <- all_op_ys[in_op]
          operon_boxes[[length(operon_boxes) + 1]] <- data.frame(
            xmin = sx - 0.35, xmax = sx + 0.35,
            ymin = min(op_ys) - 0.2, ymax = max(op_ys) + 0.2,
            label = oid, stringsAsFactors = FALSE, row.names = NULL)
        }
      }
    }
  }

  # Combine data frames
  seg_df <- if (length(seg_list) > 0) do.call(rbind, seg_list) else NULL
  ribbon_df <- if (length(ribbon_list) > 0) do.call(rbind, ribbon_list) else NULL
  label_df <- if (length(label_list) > 0) do.call(rbind, label_list) else NULL
  sub_df <- if (length(sub_nodes) > 0) do.call(rbind, sub_nodes) else NULL
  prod_df <- if (length(prod_nodes) > 0) do.call(rbind, prod_nodes) else NULL
  pts_df <- if (length(pts_cascade_list) > 0) do.call(rbind, pts_cascade_list) else NULL
  pts_label_df <- if (length(pts_label_list) > 0) do.call(rbind, pts_label_list) else NULL
  op_df <- if (length(operon_boxes) > 0) do.call(rbind, operon_boxes) else NULL

  # ---- Summary stats ----
  n_sub_total <- length(substrates)
  n_sub_hit <- sum(has_evidence)
  n_enz <- sum(vapply(rows, function(r) as.integer(nrow(r$enzymes)), integer(1)))
  n_tr  <- sum(vapply(rows, function(r) {
    if (!is.null(r$transporters)) as.integer(nrow(r$transporters)) else 0L
  }, integer(1)))
  n_pts <- sum(vapply(rows, function(r) r$has_pts, logical(1)))
  n_op  <- if (!is.null(operon_info)) length(operon_info$groups) else 0L
  subtitle_text <- paste0(
    n_sub_hit, "/", n_sub_total, " substrates with evidence \u2022 ",
    n_enz, " CAZy enzymes \u2022 ",
    n_tr, " transporters (", n_pts, " w/ PTS) \u2022 ",
    n_op, " operon clusters",
    if (n_cascades > 0) paste0(" \u2022 ", n_cascades, " degradation cascades") else "")

  build_ribbon_gradient_chunks <- function(row, n_chunk = 14L) {
    t0 <- seq(0, 1 - 1 / n_chunk, length.out = n_chunk)
    t1 <- seq(1 / n_chunk, 1, length.out = n_chunk)
    cols <- grDevices::colorRampPalette(c(
      colorspace::lighten(row$color, 0.72),
      colorspace::lighten(row$color, 0.30),
      row$color
    ))(n_chunk)
    data.frame(
      x = row$x + (row$xend - row$x) * t0,
      y = row$y + (row$yend - row$y) * t0,
      xend = row$x + (row$xend - row$x) * t1,
      yend = row$y + (row$yend - row$y) * t1,
      seg_color = cols,
      seg_alpha = seq(0.35, 1.00, length.out = n_chunk),
      stringsAsFactors = FALSE
    )
  }
  # end build_ribbon_gradient_chunks

  # ---- Build ggplot (gapmind_aa style) ----
  # ---- Layout: gapmind_aa backbone + TCA circle + PPP branch ----
  # Backbone on LEFT (x=bx), pathways branch RIGHT to membrane & substrate
  bx <- 1.0   # backbone x (same as gapmind_aa)

  # TCA cycle: perfect circle, 8 nodes at 45 deg (same as gapmind_aa)
  tca_cx <- 0.0; tca_cy <- 6.5; tca_r <- 1.0
  tca_ang <- pi/2 - (0:7) * (2 * pi / 8)
  tca_x <- tca_cx + tca_r * cos(tca_ang)
  tca_y <- tca_cy + tca_r * sin(tca_ang)
  tca_ids    <- c("OAA","Citrate","Isocit","AKG","SucCoA","Succinate","Fumarate","Malate")
  tca_labels <- c("OAA","Citrate","Isocit","\u03b1-KG","Suc-CoA","Succinate","Fumarate","Malate")

  # Glycolysis backbone nodes (gapmind_aa coordinates)
  bb_nodes <- data.frame(
    id    = c("Glucose","G6P","F6P","FBP","DHAP","GAP","3PG","PEP","Pyruvate","AcCoA"),
    x     = c(bx, bx, bx, bx, bx-1.5, bx, bx, bx, bx, bx),
    y     = c(14.0, 13.0, 12.5, 12.0, 11.5, 11.5, 11.0, 10.5, 9.5, 9.0),
    label = c("Glucose","G6P","F6P","FBP","DHAP","GAP","3-PG","PEP","Pyruvate","AcCoA"),
    stringsAsFactors = FALSE, row.names = NULL)

  # PPP branch nodes
  ppp_nodes <- data.frame(
    id    = c("6PGL","6PG","Ru5P","R5P","Xu5P"),
    x     = c(bx+1.5, bx+1.5, bx+1.5, bx+2.0, bx+1.5),
    y     = c(13.0, 12.7, 12.4, 12.1, 12.1),
    label = c("6PGL","6PG","Ru5P","R5P","Xu5P"),
    stringsAsFactors = FALSE, row.names = NULL)

  # TCA cycle nodes
  tca_nodes <- data.frame(
    id = tca_ids, x = tca_x, y = tca_y, label = tca_labels,
    stringsAsFactors = FALSE, row.names = NULL)

  # Substrate -> backbone entry point mapping
  # Entry points reference the node IDs in bb_nodes, ppp_nodes, or tca_nodes
  entry_map <- c(
    # -> G6P (hexoses via glucose-6-phosphate)
    glucose = "G6P", maltose = "G6P", cellobiose = "G6P",
    trehalose = "G6P", galactose = "G6P", lactose = "G6P",
    # -> F6P (fructose-6-phosphate entry)
    fructose = "F6P", sucrose = "F6P", mannose = "F6P",
    NAG = "F6P", glucosamine = "F6P", peptidoglycan = "F6P",
    mannitol = "F6P", sorbitol = "F6P",
    # -> Xu5P/Ru5P via PPP
    xylose = "Xu5P", arabinose = "Xu5P", ribose = "R5P",
    # -> GAP (Entner-Doudoroff)
    glucuronate = "GAP", galacturonate = "GAP", gluconate = "GAP",
    rhamnose = "GAP", glycerol = "DHAP",
    # -> PEP
    fucose = "PEP", myoinositol = "PEP",
    # -> Pyruvate
    `L-lactate` = "Pyruvate", `D-lactate` = "Pyruvate", ethanol = "Pyruvate",
    # -> AcCoA
    acetate = "AcCoA",
    # -> TCA cycle nodes directly
    citrate = "Citrate", succinate = "Succinate", fumarate = "Fumarate",
    `L-malate` = "Malate"
  )

  # ---- Build node y-lookup from all backbone/ppp/tca nodes ----
  all_nodes <- rbind(bb_nodes, ppp_nodes, tca_nodes)
  node_y <- setNames(all_nodes$y, all_nodes$id)
  node_x <- setNames(all_nodes$x, all_nodes$id)

  # ---- Show ALL substrates from entry_map ----
  all_substrates <- names(entry_map)
  has_evidence <- all_substrates %in% names(rows)
  substrates <- all_substrates
  sub_entry <- entry_map[substrates]

  # Group substrates by entry point, evidence-positive first
  unique_entries <- unique(sub_entry)
  entry_groups <- list()
  for (eid in unique_entries) {
    subs <- substrates[sub_entry == eid]
    subs_pos <- subs[subs %in% names(rows)]
    subs_neg <- subs[!subs %in% names(rows)]
    entry_groups[[eid]] <- c(subs_pos, subs_neg)
  }

  # Assign substrate y: fan out from entry node y, spaced by 0.55 (evidence) / 0.25 (none)
  sub_y <- numeric(length(substrates)); names(sub_y) <- substrates
  for (eid in names(entry_groups)) {
    subs <- entry_groups[[eid]]
    entry_y <- if (eid %in% names(node_y)) node_y[eid] else 10.0
    n <- length(subs)
    row_h <- ifelse(subs %in% names(rows), 0.55, 0.25)
    cum_h <- cumsum(c(0, row_h[-n]))
    total_h <- sum(row_h) - row_h[n]
    offset <- entry_y + total_h / 2 - cum_h
    for (i in seq_along(subs)) sub_y[subs[i]] <- offset[i]
  }

  # ---- X positions for pathway chain ----
  x_product  <- 2.5   # degradation product node
  x_enz_start <- 4.0  # first enzyme segment start
  x_enz_step  <- 1.2  # spacing per enzyme
  x_membrane  <- 9.0  # membrane line
  x_tr_start  <- 10.0 # first transport segment
  x_tr_step   <- 1.0  # spacing per transporter
  x_substrate <- 14.0 # substrate node (right end)

  # ---- Confidence colors (same as gapmind_aa) ----
  line_col <- c(very_high = "#1A9641", high = "#2CA25F", medium = "#FEC44F",
                low = "#F03B20", none = "#CCCCCC")
  conf_levels <- c("very_high", "high", "medium", "low", "none")

  # ---- Pathway ribbon colors (unique per substrate, grouped by entry) ----
  pw_colors <- c(
    # G6P group
    glucose = "#E6550D", maltose = "#B8860B", cellobiose = "#D94801",
    trehalose = "#8C6BB1", galactose = "#9467BD", lactose = "#C49C94",
    # F6P group
    fructose = "#2CA02C", sucrose = "#31A354", mannose = "#D6604D",
    NAG = "#3182BD", glucosamine = "#E7298A", peptidoglycan = "#8C564B",
    mannitol = "#17BECF", sorbitol = "#9EDAE5",
    # GAP group
    xylose = "#FF8C00", arabinose = "#FD8D3C", ribose = "#FFBB78",
    glucuronate = "#636363", galacturonate = "#969696", gluconate = "#BCBD22",
    rhamnose = "#AEC7E8", glycerol = "#7F7F7F",
    # PEP group
    fucose = "#9C9EDE", myoinositol = "#CE93D8",
    # Pyruvate group
    `L-lactate` = "#6BAED6", `D-lactate` = "#4292C6", ethanol = "#2171B5",
    # AcCoA group
    acetate = "#756BB1",
    # TCA group
    citrate = "#006D2C", succinate = "#41AB5D", fumarate = "#74C476",
    `L-malate` = "#A1D99B"
  )

  # ---- Build segment data for all pathways (gapmind_aa style) ----
  seg_list <- list()    # enzyme/transport segments with confidence coloring
  ribbon_list <- list() # background ribbon gradients
  dot_list <- list()    # intermediate node circles
  label_list <- list()  # gene name + locus_tag labels
  sub_nodes <- list()   # substrate "product" nodes (bold border, right end)
  prod_nodes <- list()  # degradation product nodes
  pts_cascade_list <- list() # PTS cascade mini-segments
  operon_boxes <- list()

  for (si in seq_along(substrates)) {
    sub <- substrates[si]
    y <- sub_y[sub]
    pw_col <- unname(pw_colors[sub])
    if (is.na(pw_col) || is.null(pw_col)) pw_col <- "#999999"
    entry_id <- entry_map[sub]; if (is.na(entry_id)) entry_id <- "G6P"
    has_data <- sub %in% names(rows)

    # --- Substrate node (right end) --- always shown
    sub_nodes[[length(sub_nodes) + 1]] <- data.frame(
      x = x_substrate, y = y, label = sub, has_data = has_data,
      stringsAsFactors = FALSE, row.names = NULL)

    if (!has_data) {
      # ---- NO EVIDENCE: short grey connector from membrane area to substrate ----
      ribbon_list[[length(ribbon_list) + 1]] <- data.frame(
        x = x_membrane - 1, y = y, xend = x_substrate, yend = y,
        color = "#E0E0E0", stringsAsFactors = FALSE, row.names = NULL)
      next
    }

    # ---- HAS EVIDENCE: full pathway chain ----
    r <- rows[[sub]]

    # --- Entry connection (backbone -> product) ---
    prods <- r$products
    n_prd <- length(prods)
    prd_ys <- if (n_prd <= 1) y else seq(y + 0.15 * (n_prd - 1), y - 0.15 * (n_prd - 1), length.out = n_prd)
    for (j in seq_len(n_prd)) {
      prod_nodes[[length(prod_nodes) + 1]] <- data.frame(
        x = x_product, y = prd_ys[j], label = prods[j], stringsAsFactors = FALSE, row.names = NULL)
      # backbone_entry -> product (colored ribbon)
      entry_x <- if (entry_id %in% names(node_x)) node_x[entry_id] else bx
      entry_yy <- if (entry_id %in% names(node_y)) node_y[entry_id] else y
      ribbon_list[[length(ribbon_list) + 1]] <- data.frame(
        x = entry_x, y = entry_yy, xend = x_product, yend = prd_ys[j],
        color = pw_col, stringsAsFactors = FALSE, row.names = NULL)
    }

    # --- CAZy enzyme segments ---
    enz <- r$enzymes
    n_enz <- nrow(enz)
    if (n_enz > 6) {
      enz <- enz[order(-enz$evidence_n), , drop = FALSE][1:6, , drop = FALSE]
      n_enz <- 6
    }
    enz_xs <- if (n_enz == 0) numeric(0) else {
      seq(x_enz_start, by = x_enz_step, length.out = n_enz)
    }
    for (j in seq_len(n_enz)) {
      ev_conf <- enz$evidence_conf[j]
      ev_col <- unname(line_col[ev_conf])
      if (is.na(ev_col)) ev_col <- "#CCCCCC"
      x1 <- if (j == 1) x_product else enz_xs[j - 1]
      x2 <- enz_xs[j]
      seg_list[[length(seg_list) + 1]] <- data.frame(
        x = x1, y = y, xend = x2, yend = y,
        confidence = ev_conf, color = ev_col,
        stringsAsFactors = FALSE, row.names = NULL)
      # Gene label — on segment midpoint, offset above
      gene_lab <- enz$enzyme_label[j]
      lt <- enz$locus_tag[j]
      label_list[[length(label_list) + 1]] <- data.frame(
        x = (x1 + x2) / 2, y = y + 0.30,
        label = paste0(gene_lab, "\n", lt),
        color = ifelse(ev_conf == "none", "#AAAAAA", "#222222"),
        stringsAsFactors = FALSE, row.names = NULL)
      # Intermediate dot
      dot_list[[length(dot_list) + 1]] <- data.frame(x = x2, y = y, stringsAsFactors = FALSE, row.names = NULL)
    }

    # --- Membrane -> Transport -> Substrate chain ---
    tr_rows <- r$transporters
    pts_comp <- r$pts_components
    has_pts <- r$has_pts && !is.null(pts_comp) && nrow(pts_comp) > 0

    # Get non-PTS transporters
    tr_items <- NULL
    if (!is.null(tr_rows) && nrow(tr_rows) > 0) {
      tr_items <- tr_rows[!duplicated(tr_rows$locus_tag), , drop = FALSE]
    }

    # Enzyme end -> membrane connection
    enz_end_x <- if (n_enz > 0) max(enz_xs) else x_product
    seg_list[[length(seg_list) + 1]] <- data.frame(
      x = enz_end_x, y = y, xend = x_membrane, yend = y,
      confidence = "none", color = "#D5D5D5", stringsAsFactors = FALSE, row.names = NULL)

    # PTS cascade — biological vertical stacking around membrane
    # EIIC: on membrane, EIIB: cytoplasmic face, EIIA/HPr/EI: cytoplasm
    if (has_pts) {
      unique_lts <- unique(pts_comp$locus_tag)
      if (length(unique_lts) > 5) unique_lts <- unique_lts[1:5]
      # Vertical y-offsets by component (relative to pathway row y)
      pts_y_offset <- c(EIIC = 0, `EIIB/C` = 0, EIICB = 0, EIICBA = 0,
                        EIIB = -0.25, EIIA = -0.45, HPr = -0.60, EI = -0.75,
                        `EI-HPr-EIIA` = -0.50, PTS = -0.30)
      for (lt in unique_lts) {
        lt_rows <- pts_comp[pts_comp$locus_tag == lt, , drop = FALSE]
        comps <- unique(lt_rows$pts_component)
        for (comp_i in comps) {
          dy <- unname(pts_y_offset[comp_i])
          if (is.na(dy)) dy <- -0.30
          pts_cascade_list[[length(pts_cascade_list) + 1]] <- data.frame(
            x = x_membrane - 0.4, y = y + dy,
            xend = x_membrane + 0.4, yend = y + dy,
            comp = comp_i, lt = lt,
            is_fused = length(comps) > 1, stringsAsFactors = FALSE, row.names = NULL)
        }
      }
    }

    # Non-PTS transporters — vertical stacking on membrane
    # ABC: SBP above membrane, TMD on membrane, NBD below
    # MFS/other: on membrane line
    nonpts <- if (!is.null(tr_items)) tr_items[!tr_items$is_pts, , drop = FALSE] else NULL
    if (!is.null(nonpts) && nrow(nonpts) > 3) {
      conf_rank <- c(high = 1, medium = 2, low = 3, none = 4)
      nonpts$rank <- conf_rank[nonpts$confidence]
      nonpts <- nonpts[order(nonpts$rank), , drop = FALSE][1:3, ]
    }
    n_nonpts <- if (!is.null(nonpts)) nrow(nonpts) else 0L
    for (j in seq_len(n_nonpts)) {
      tr_conf <- nonpts$confidence[j]
      tr_col <- unname(line_col[tr_conf])
      if (is.na(tr_col)) tr_col <- "#CCCCCC"
      tr_gene <- nonpts$gene_name[j]
      if (is.na(tr_gene) || !nzchar(tr_gene)) tr_gene <- nonpts$locus_tag[j]
      # Detect ABC subunits by gene name pattern
      is_sbp <- grepl("SBP|[Bb]inding|[Pp]eri", tr_gene, ignore.case = TRUE)
      is_nbd <- grepl("NBD|ATP|[Aa]ase", tr_gene, ignore.case = TRUE)
      tr_dy <- if (is_sbp) 0.25 else if (is_nbd) -0.25 else 0
      # Place as short horizontal segment on/near membrane
      seg_list[[length(seg_list) + 1]] <- data.frame(
        x = x_membrane - 0.3, y = y + tr_dy,
        xend = x_membrane + 0.3, yend = y + tr_dy,
        confidence = tr_conf, color = tr_col, stringsAsFactors = FALSE, row.names = NULL)
      # Label on cytoplasmic side (below) — only gene name, no locus_tag
      label_list[[length(label_list) + 1]] <- data.frame(
        x = x_membrane + 0.6, y = y + tr_dy,
        label = tr_gene,
        color = ifelse(tr_conf == "none", "#AAAAAA", "#222222"),
        stringsAsFactors = FALSE, row.names = NULL)
    }

    # Connection: enzyme end → membrane, then membrane → substrate
    ribbon_list[[length(ribbon_list) + 1]] <- data.frame(
      x = x_membrane + 0.5, y = y, xend = x_substrate, yend = y,
      color = pw_col, stringsAsFactors = FALSE, row.names = NULL)

    # Operon brackets (enzyme side only — transporters now stacked on membrane)
    if (length(r$operon_ids) > 0 && !is.null(operon_info)) {
      all_op_xs <- enz_xs
      all_op_lts <- enz$locus_tag
      for (oid in r$operon_ids) {
        members <- operon_info$groups[[oid]]
        in_op <- which(all_op_lts %in% members)
        if (length(in_op) >= 2) {
          op_xs <- all_op_xs[in_op]
          operon_boxes[[length(operon_boxes) + 1]] <- data.frame(
            xmin = min(op_xs) - 0.3, xmax = max(op_xs) + 0.3,
            ymin = y - 0.35, ymax = y + 0.35,
            label = oid, stringsAsFactors = FALSE, row.names = NULL)
        }
      }
    }
  }

  # Combine data frames
  seg_df <- if (length(seg_list) > 0) do.call(rbind, seg_list) else NULL
  ribbon_df <- if (length(ribbon_list) > 0) do.call(rbind, ribbon_list) else NULL
  dot_df <- if (length(dot_list) > 0) do.call(rbind, dot_list) else NULL
  label_df <- if (length(label_list) > 0) do.call(rbind, label_list) else NULL
  sub_df <- if (length(sub_nodes) > 0) do.call(rbind, sub_nodes) else NULL
  prod_df <- if (length(prod_nodes) > 0) do.call(rbind, prod_nodes) else NULL
  pts_df <- if (length(pts_cascade_list) > 0) do.call(rbind, pts_cascade_list) else NULL
  op_df <- if (length(operon_boxes) > 0) do.call(rbind, operon_boxes) else NULL

  # ---- Summary stats ----
  n_sub_total <- length(substrates)
  n_sub_hit <- sum(has_evidence)
  n_enz <- sum(vapply(rows, function(r) as.integer(nrow(r$enzymes)), integer(1)))
  n_tr  <- sum(vapply(rows, function(r) {
    if (!is.null(r$transporters)) as.integer(nrow(r$transporters)) else 0L
  }, integer(1)))
  n_pts <- sum(vapply(rows, function(r) r$has_pts, logical(1)))
  n_op  <- if (!is.null(op_df)) nrow(op_df) else 0L
  subtitle_text <- paste0(
    n_sub_hit, "/", n_sub_total, " substrates with evidence \u2022 ",
    n_enz, " CAZy enzymes \u2022 ",
    n_tr, " transporters (", n_pts, " w/ PTS) \u2022 ",
    n_op, " operon clusters",
    if (n_cascades > 0) paste0(" \u2022 ", n_cascades, " degradation cascades") else "")


  p <- ggplot2::ggplot()

  # ---- L-1: Cell shape background (cytoplasm + membrane + extracellular) ----
  cell_xmin <- -1.0
  cell_xmax <- x_substrate + 2.0
  cell_ymin <- min(bb_nodes$y, sub_y, tca_y) - 1.0
  cell_ymax <- max(bb_nodes$y, sub_y, ppp_nodes$y) + 1.0
  cell_layers <- .dnmb_cct_cell_shape_layers(
    xmin = cell_xmin, xmax = cell_xmax,
    ymin = cell_ymin, ymax = cell_ymax,
    membrane_x = x_membrane, membrane_width = 0.4)
  for (cl in cell_layers) p <- p + cl

  # ---- L0: Ribbon gradients (simplified: 2-layer glow for evidence, thin for no-evidence) ----
  if (!is.null(ribbon_df) && nrow(ribbon_df) > 0) {
    # No-evidence ribbons (grey) = thin single line
    noev <- ribbon_df[ribbon_df$color == "#E0E0E0", , drop = FALSE]
    if (nrow(noev) > 0) {
      p <- p + ggplot2::geom_segment(data = noev,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        color = "#E0E0E0", linewidth = 0.4, alpha = 0.5, lineend = "butt")
    }
    # Evidence ribbons = 2-layer soft glow
    ev <- ribbon_df[ribbon_df$color != "#E0E0E0", , drop = FALSE]
    if (nrow(ev) > 0) {
      for (li in seq_along(c(4, 1.5))) {
        lw <- c(4, 1.5)[li]; a <- c(0.12, 0.25)[li]
        for (i in seq_len(nrow(ev))) {
          row <- ev[i, , drop = FALSE]
          p <- p + ggplot2::geom_segment(data = row,
            ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
            color = row$color, linewidth = lw, alpha = a, lineend = "butt")
        }
      }
    }
  }

  # ---- L1: Glycolysis backbone + TCA circle + PPP branch ----
  # Glycolysis: sequential segments
  bb_seg <- data.frame(
    x = bb_nodes$x[-nrow(bb_nodes)], y = bb_nodes$y[-nrow(bb_nodes)],
    xend = bb_nodes$x[-1], yend = bb_nodes$y[-1],
    stringsAsFactors = FALSE, row.names = NULL)
  p <- p +
    ggplot2::geom_segment(data = bb_seg,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "#D5D5D5", linewidth = 1.8)

  # AcCoA -> OAA connection (glycolysis -> TCA)
  p <- p + ggplot2::annotate("segment",
    x = bx, y = 9.0, xend = tca_x[1], yend = tca_y[1],
    color = "#D5D5D5", linewidth = 1.5)

  # TCA cycle: segments around the circle
  for (ti in seq_along(tca_ids)) {
    ti_next <- if (ti == length(tca_ids)) 1L else ti + 1L
    p <- p + ggplot2::annotate("segment",
      x = tca_x[ti], y = tca_y[ti],
      xend = tca_x[ti_next], yend = tca_y[ti_next],
      color = "#C8E6C9", linewidth = 1.2)
  }
  # TCA nodes
  p <- p +
    ggplot2::geom_point(data = tca_nodes,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 5, fill = "#E8F5E9", color = "#81C784", stroke = 0.3) +
    ggplot2::geom_text(data = tca_nodes,
      ggplot2::aes(x = .data$x, y = .data$y - 0.25, label = .data$label),
      size = 1.4, fontface = "bold", color = "#2E7D32")

  # PPP branch: G6P -> 6PGL -> 6PG -> Ru5P -> R5P / Xu5P
  ppp_route <- c("G6P","6PGL","6PG","Ru5P")
  for (pi in seq_along(ppp_route)[-1]) {
    p <- p + ggplot2::annotate("segment",
      x = node_x[ppp_route[pi-1]], y = node_y[ppp_route[pi-1]],
      xend = node_x[ppp_route[pi]], yend = node_y[ppp_route[pi]],
      color = "#B3E5FC", linewidth = 1.2)
  }
  # Ru5P -> R5P and Ru5P -> Xu5P
  p <- p +
    ggplot2::annotate("segment",
      x = node_x["Ru5P"], y = node_y["Ru5P"],
      xend = node_x["R5P"], yend = node_y["R5P"],
      color = "#B3E5FC", linewidth = 1.0) +
    ggplot2::annotate("segment",
      x = node_x["Ru5P"], y = node_y["Ru5P"],
      xend = node_x["Xu5P"], yend = node_y["Xu5P"],
      color = "#B3E5FC", linewidth = 1.0)
  # PPP nodes
  p <- p +
    ggplot2::geom_point(data = ppp_nodes,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 5, fill = "#E1F5FE", color = "#4FC3F7", stroke = 0.3) +
    ggplot2::geom_text(data = ppp_nodes,
      ggplot2::aes(x = .data$x + 0.3, y = .data$y, label = .data$label),
      size = 1.4, fontface = "bold", color = "#0277BD", hjust = 0)

  # Glycolysis backbone nodes + labels
  p <- p +
    ggplot2::geom_point(data = bb_nodes,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 7, fill = "#F0F0F0", color = "#AAAAAA", stroke = 0.3) +
    ggplot2::geom_text(data = bb_nodes,
      ggplot2::aes(x = .data$x - 0.4, y = .data$y, label = .data$label),
      size = 1.8, fontface = "bold", color = "#333333", hjust = 1)

  # ---- L2: Membrane line (vertical dashed) ----
  mem_y_range <- range(sub_y) + c(-1.0, 1.0)
  p <- p +
    ggplot2::annotate("segment",
      x = x_membrane, xend = x_membrane,
      y = mem_y_range[1], yend = mem_y_range[2],
      color = "#D4A574", linewidth = 1.2, linetype = "solid") +
    ggplot2::annotate("text", x = x_membrane, y = mem_y_range[2] + 0.3,
      label = "Membrane", size = 2.5, fontface = "bold", color = "#8D6E63")

  # ---- L3: Enzyme/transport segments (colored by confidence) ----
  if (!is.null(seg_df) && nrow(seg_df) > 0) {
    seg_df$confidence <- factor(seg_df$confidence, levels = conf_levels)
    p <- p +
      ggplot2::geom_segment(data = seg_df,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                     color = .data$confidence),
        linewidth = 0.8, lineend = "round", show.legend = TRUE) +
      ggplot2::scale_color_manual(
        name = "Confidence",
        values = c(very_high = "#1A9641", high = "#2CA25F", medium = "#FEC44F",
                   low = "#F03B20", none = "#CCCCCC"),
        labels = c(very_high = "Very High", high = "High", medium = "Medium",
                   low = "Low", none = "Not found"),
        drop = FALSE)
  }

  # ---- L4: Intermediate dots (small, unobtrusive — no white fill background) ----
  if (!is.null(dot_df) && nrow(dot_df) > 0) {
    dot_df <- dot_df[!duplicated(dot_df), , drop = FALSE]
    p <- p +
      ggplot2::geom_point(data = dot_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        shape = 21, size = 3, fill = "#F5F5F5", color = "#CCCCCC", stroke = 0.3)
  }

  # ---- L5: Product nodes (degradation products) ----
  if (!is.null(prod_df) && nrow(prod_df) > 0) {
    prod_df <- prod_df[!duplicated(prod_df), , drop = FALSE]
    p <- p +
      ggplot2::geom_point(data = prod_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        shape = 21, size = 9, fill = "#C8E6C9", color = "#CCCCCC", stroke = 0.3) +
      ggplot2::geom_text(data = prod_df,
        ggplot2::aes(x = .data$x, y = .data$y - 0.22, label = .data$label),
        size = 1.5, color = "#1B5E20", fontface = "italic")
  }

  # ---- L6: PTS/transporter — vertically stacked around membrane line ----
  if (!is.null(pts_df) && nrow(pts_df) > 0) {
    pts_comp_colors <- c(EI = "#FF8A65", HPr = "#FFB74D", EIIA = "#FFD54F",
                         `EIIB/C` = "#AED581", EIIC = "#4DB6AC",
                         EIICB = "#4FC3F7", EIICBA = "#7986CB",
                         `EI-HPr-EIIA` = "#CE93D8", PTS = "#BDBDBD")
    fused_brackets <- list()
    for (i in seq_len(nrow(pts_df))) {
      row <- pts_df[i, , drop = FALSE]
      cc <- unname(pts_comp_colors[row$comp])
      if (is.na(cc)) cc <- "#BDBDBD"
      is_fused_seg <- if ("is_fused" %in% names(row)) isTRUE(row$is_fused) else FALSE
      seg_lw <- if (is_fused_seg) 2.0 else 2.8
      p <- p +
        ggplot2::annotate("segment",
          x = row$x, y = row$y, xend = row$xend, yend = row$yend,
          color = cc, linewidth = seg_lw, alpha = 0.7, lineend = "round") +
        ggplot2::annotate("text",
          x = row$xend + 0.15, y = row$y,
          label = row$comp, size = 0.9, fontface = "bold", color = "#E65100",
          hjust = 0)
      if (is_fused_seg) {
        lt <- row$lt
        if (!lt %in% names(fused_brackets)) {
          fused_brackets[[lt]] <- list(xmin = row$x, xmax = row$xend,
            ymin = row$y, ymax = row$y, lt = lt)
        } else {
          fused_brackets[[lt]]$xmin <- min(fused_brackets[[lt]]$xmin, row$x)
          fused_brackets[[lt]]$xmax <- max(fused_brackets[[lt]]$xmax, row$xend)
          fused_brackets[[lt]]$ymin <- min(fused_brackets[[lt]]$ymin, row$y)
          fused_brackets[[lt]]$ymax <- max(fused_brackets[[lt]]$ymax, row$y)
        }
      }
    }
    for (fb in fused_brackets) {
      p <- p +
        ggplot2::annotate("rect",
          xmin = fb$xmin - 0.05, xmax = fb$xmax + 0.05,
          ymin = fb$ymin - 0.08, ymax = fb$ymax + 0.08,
          fill = NA, color = "#795548", linewidth = 0.3, linetype = "dashed")
    }
  }

  # ---- L7: Gene labels — cytoplasmic side only (no locus_tags in extracellular) ----
  # Remove any labels placed beyond the membrane (extracellular zone)
  if (!is.null(label_df) && nrow(label_df) > 0) {
    label_df <- label_df[label_df$x <= x_membrane + 0.1, , drop = FALSE]
  }
  if (!is.null(label_df) && nrow(label_df) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p +
        ggrepel::geom_text_repel(data = label_df,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
          size = 1.5, color = label_df$color, fontface = "bold", lineheight = 0.8,
          max.overlaps = Inf, segment.size = 0.2, segment.color = "#CCCCCC",
          box.padding = 0.15, point.padding = 0.1, force = 2,
          min.segment.length = 0, seed = 42)
    } else {
      p <- p +
        ggplot2::geom_text(data = label_df,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
          size = 1.5, color = label_df$color, fontface = "bold", lineheight = 0.8,
          check_overlap = TRUE)
    }
  }

  # ---- L8: Substrate nodes (composite SNFG: mono/di/tri + bonds) ----
  if (!is.null(sub_df) && nrow(sub_df) > 0) {
    for (i in seq_len(nrow(sub_df))) {
      sub_name <- sub_df$label[i]
      scol <- unname(pw_colors[sub_name])
      if (is.na(scol)) scol <- "#999999"
      has_d <- sub_df$has_data[i]

      lab_fill <- if (has_d) scol else "#DDDDDD"
      lab_col  <- if (has_d) "white" else "#999999"
      lab_face <- if (has_d) "bold" else "plain"

      # Render composite SNFG node (handles mono/di/tri/non-sugar)
      node_layers <- .dnmb_snfg_render_metabolite_node(
        cx = sub_df$x[i], cy = sub_df$y[i],
        metabolite = sub_name, active = has_d,
        symbol_r = 0.33, point_size = 9, spacing = 0.70)
      for (nl in node_layers) p <- p + nl

      # Substrate label (offset right of the rightmost symbol)
      comp <- .dnmb_snfg_metabolite_composition(sub_name)
      n_mono <- length(comp$monomers)
      label_offset <- if (n_mono > 1) (n_mono - 1) * 0.70 / 2 + 0.65 else 0.55
      p <- p + ggplot2::geom_label(
        data = sub_df[i, , drop = FALSE],
        ggplot2::aes(x = .data$x + label_offset, y = .data$y, label = .data$label),
        size = 2.5, fontface = lab_face, color = lab_col, fill = lab_fill,
        label.padding = grid::unit(0.12, "lines"), linewidth = 0)
    }
  }

  # ---- L9: Operon brackets (dashed boxes) ----
  if (!is.null(op_df) && nrow(op_df) > 0) {
    p <- p +
      ggplot2::geom_rect(data = op_df,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$ymin, ymax = .data$ymax),
        fill = "#FFF3E0", color = "#FF8F00", linewidth = 0.5,
        linetype = "dashed", alpha = 0.25) +
      ggplot2::geom_text(data = op_df,
        ggplot2::aes(x = (.data$xmin + .data$xmax) / 2,
                     y = .data$ymax + 0.12, label = .data$label),
        size = 1.5, color = "#E65100", fontface = "bold")
  }

  # Collect locus_tags already shown in enzyme pathway labels
  displayed_lts <- if (!is.null(label_df)) unique(unlist(strsplit(label_df$label, "\n"))) else character(0)

  # ---- L10: Extracellular complex carbohydrate SNFG cascade panel ----
  if (n_cascades > 0) {
    # Auto-placement in extracellular zone (right of substrate nodes)
    casc_x_start <- x_substrate + 2.0
    casc_unit_sz <- 0.32
    casc_row_h   <- 2.2   # vertical spacing between cascades
    casc_y_start <- max(sub_y)  # start from top of pathway area

    # Section title
    p <- p + ggplot2::annotate("text",
      x = casc_x_start + 2.0, y = casc_y_start + 1.0,
      label = "Extracellular Polysaccharide Degradation (SNFG)",
      size = 2.8, fontface = "bold", color = "#4E342E")

    for (ci in seq_along(active_cascades)) {
      ac <- active_cascades[[ci]]
      casc <- ac$cascade
      hits <- ac$step_hits
      cy <- casc_y_start - (ci - 1) * casc_row_h

      # Polymer name label (left of chain)
      p <- p + ggplot2::annotate("text",
        x = casc_x_start - 0.3, y = cy,
        label = casc$polymer, size = 1.6, fontface = "bold",
        color = casc$chain_color, hjust = 1)

      n_steps <- length(hits)
      n_units <- n_steps + 2
      cut_pos <- integer(0)
      gh_labs <- character(0)
      ct_labs <- character(0)
      lt_annotations <- list()

      for (si in seq_along(hits)) {
        hit <- hits[[si]]
        if (length(hit$matched_families) > 0) {
          cut_pos <- c(cut_pos, si)
          gh_labs <- c(gh_labs, paste(hit$matched_families, collapse = "/"))
          ct_labs <- c(ct_labs, hit$cleavage_type)
          if (length(hit$locus_tags) > 0) {
            lt_annotations[[length(lt_annotations) + 1]] <- list(
              pos = si, lts = hit$locus_tags)
          }
        }
      }

      # Draw SNFG sugar chain with proper monosaccharide symbols
      chain_layers <- .dnmb_cct_sugar_chain_layers(
        x_start = casc_x_start, y = cy,
        n_units = n_units, unit_size = casc_unit_sz,
        cut_positions = cut_pos, gh_labels = gh_labs,
        cleavage_types = ct_labs, chain_color = casc$chain_color,
        unit_label = casc$unit_label,
        unit_snfg = casc$unit_snfg,
        bond_label = casc$bond_label)
      for (cl in chain_layers) p <- p + cl

      # Locus tag annotations below cleavage sites
      spacing <- casc_unit_sz * 1.4
      unit_xs <- casc_x_start + (seq_len(n_units) - 1) * spacing
      # Locus_tags omitted from extracellular zone — GH family labels above
      # scissors are sufficient enzyme identification

      # Dashed arrow from chain end -> matching substrate node in pathway
      chain_end_x <- max(unit_xs) + casc_unit_sz * 0.5
      # Find the substrate node y for this cascade's monomer
      sub_target <- casc$monomer
      if (sub_target %in% names(sub_y)) {
        target_y <- sub_y[sub_target]
      } else {
        entry_id <- casc$entry
        target_y <- if (entry_id %in% names(node_y)) node_y[entry_id] else cy
      }
      p <- p + ggplot2::annotate("segment",
        x = casc_x_start - 0.1, y = cy,
        xend = x_substrate + 0.8, yend = target_y,
        color = casc$chain_color, linewidth = 0.4,
        linetype = "dashed", alpha = 0.5,
        arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "closed"))
      p <- p + ggplot2::annotate("text",
        x = chain_end_x + 0.3, y = cy,
        label = paste0("\u2192 ", casc$monomer),
        size = 1.3, color = casc$chain_color, fontface = "bold", hjust = 0)
    }
  }

  # ---- Zone labels (top) ----
  top_y <- max(sub_y) + 1.2
  if (n_cascades > 0) top_y <- max(sub_y) + 1.5
  p <- p +
    ggplot2::annotate("text", x = 0, y = top_y,
      label = "Glycolysis", size = 3.0, fontface = "bold", color = "#333333") +
    ggplot2::annotate("text", x = x_product, y = top_y,
      label = "Product", size = 3.0, fontface = "bold", color = "#2E7D32") +
    ggplot2::annotate("text", x = (x_enz_start + x_membrane) / 2, y = top_y,
      label = "CAZy Enzymes", size = 3.0, fontface = "bold", color = "#333333") +
    ggplot2::annotate("text", x = x_substrate, y = top_y,
      label = "Substrate", size = 3.0, fontface = "bold", color = "#333333")

  # ---- SNFG Shape Legend (bottom-right) ----
  snfg_legend <- data.frame(
    label = c("Hexose (circle)", "HexNAc (square)", "Hexosamine (crossed sq)",
              "Hexuronate (split diamond)", "Deoxyhexose (triangle)",
              "Pentose (star)", "Ketose (pentagon)"),
    representative = c("Glc", "GlcNAc", "GlcN", "GlcA", "Fuc", "Xyl", "Fru"),
    stringsAsFactors = FALSE, row.names = NULL)
  leg_x0 <- x_substrate + 1.5
  leg_y0 <- min(sub_y) - 0.5
  p <- p + ggplot2::annotate("text",
    x = leg_x0, y = leg_y0 + 0.5,
    label = "SNFG Shapes", size = 2.0, fontface = "bold", color = "#333333", hjust = 0)
  for (li in seq_len(nrow(snfg_legend))) {
    ly <- leg_y0 - (li - 1) * 0.5
    leg_layers <- .dnmb_snfg_render_symbol(leg_x0, ly,
      snfg_legend$representative[li], r = 0.27, point_size = 7, active = TRUE)
    for (ll in leg_layers) p <- p + ll
    p <- p + ggplot2::annotate("text", x = leg_x0 + 0.5, y = ly,
      label = snfg_legend$label[li], size = 1.4, color = "#333333", hjust = 0)
  }

  # ---- Coord, theme, save (gapmind_aa style) ----
  # Expand right for extracellular SNFG cascade chains
  max_x <- x_substrate + 3
  if (n_cascades > 0) {
    max_chain_units <- max(vapply(active_cascades, function(ac)
      length(ac$step_hits) + 2L, integer(1)))
    max_x <- max(max_x, x_substrate + 2.0 + max_chain_units * 0.32 * 1.4 + 3.0)
  }
  y_lo <- min(bb_nodes$y, sub_y, tca_y) - 1.5 - nrow(snfg_legend) * 0.4
  y_hi <- max(bb_nodes$y, sub_y, ppp_nodes$y) + 2.0
  if (n_cascades > 0) y_hi <- max(y_hi, max(sub_y) + 1.8)
  p <- p +
    ggplot2::coord_fixed(ratio = 1,
      xlim = c(-1.5, max_x), ylim = c(y_lo, y_hi), clip = "off") +
    ggplot2::labs(
      title = "CAZy + Carbon + Membrane Transport Integrated Pathway Map",
      subtitle = subtitle_text) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0,
                                         margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 7, hjust = 0, color = "#666666",
                                            margin = ggplot2::margin(b = 3)),
      legend.position = c(0.0, 0.95),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 7),
      legend.key.width = grid::unit(16, "pt"),
      legend.key.height = grid::unit(6, "pt"),
      plot.margin = ggplot2::margin(5, 5, 5, 5))

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(file_stub, ".pdf"))
  x_range <- diff(c(-1.5, max_x))
  y_range <- diff(c(y_lo, y_hi))
  plot_scale <- 1.0
  .dnmb_module_plot_save(p, pdf_path,
    width = x_range * plot_scale + 1,
    height = y_range * plot_scale + 1)
  list(pdf = pdf_path)
}

#' Match GH enzymes to cascade gh_prefixes (prefix matching for subfamilies)
#' @keywords internal
.dnmb_cct_3zone_match_gh <- function(gh_enzymes, gh_prefixes) {
  if (is.null(gh_enzymes) || !is.data.frame(gh_enzymes) || nrow(gh_enzymes) == 0) {
    return(data.frame(gh_family = character(0), locus_tag = character(0),
                      gene_name = character(0), stringsAsFactors = FALSE, row.names = NULL))
  }
  # Build regex: GH13 should match GH13, GH13_20, GH13_31 etc.
  patterns <- paste0("^(", paste(gh_prefixes, collapse = "|"), ")(\\b|_|$)")
  matched <- grepl(patterns, gh_enzymes$gh_family)
  gh_enzymes[matched, , drop = FALSE]
}

#' Classify GH family as endo-acting (cuts mid-chain) or exo-acting (cuts end)
#' @keywords internal
.dnmb_gh_is_endo <- function(gh_family) {
  base <- sub("(_.*|$)", "", gh_family)
  endo <- c("GH5","GH6","GH7","GH9","GH10","GH11","GH13","GH18","GH19",
            "GH23","GH28","GH44","GH45","GH48","GH26","GH8")
  base %in% endo
}

# Entry intermediate node builder — x positions derived from carbon source positions
.dnmb_cct_entry_intermediate_nodes <- function(cs_xs, bx, lx, by, s, n) {
  g <- function(v) round(v * 2) / 2  # snap to 0.5 grid

  # Map each entry intermediate to the mean x of its associated carbon sources
  # If the associated carbon source is absent, use lx (backbone left)
  cs_mean_x <- function(ids, fallback = lx) {
    present <- ids[ids %in% names(cs_xs)]
    if (length(present) > 0) g(mean(cs_xs[present])) else fallback
  }

  # Glucose family intermediates: Glc-1-P links to Maltose/Cellobiose/Trehalose
  glc_x <- cs_mean_x(c("Maltose", "Cellobiose", "Trehalose"))
  # Galactose family: Gal-1-P links to Galactose/Lactose
  gal_x <- cs_mean_x(c("Galactose", "Lactose"))
  # GlcNAc family: GlcNAc-6-P links to NAG/Glucosamine
  nac_x <- cs_mean_x(c("NAG", "Glucosamine"))
  # Mannose family: Man-6-P links to Mannose/Mannitol
  man_x <- cs_mean_x(c("Mannose", "Mannitol"))
  # Fructose family: Fru-1-P links to Fructose/Sucrose
  fru_x <- cs_mean_x(c("Fructose", "Sucrose"))
  # Glycerol: Glycerol-3-P links to Glycerol
  gly_x <- cs_mean_x(c("Glycerol"))

  # Spread x positions: offset from backbone to avoid overlap
  # Each entry intermediate gets x between its carbon source and backbone
  spread_x <- function(cs_x) g((cs_x + bx) / 2)

  rbind(
    n("Glc-1-P",      spread_x(glc_x), by,        "Glc-1-P",      "entry_intermediate", "glucose"),
    n("Gal-1-P",      spread_x(gal_x), by - 1*s,  "Gal-1-P",      "entry_intermediate", "galactose"),
    n("GlcNAc-6-P",   spread_x(nac_x), by - 2*s,  "GlcNAc-6-P",   "entry_intermediate", "glcnac"),
    n("Man-6-P",      spread_x(man_x), by - 2.5*s,"Man-6-P",       "entry_intermediate", "mannose"),
    n("Fru-1-P",      spread_x(fru_x), by - 3.5*s,"Fru-1-P",       "entry_intermediate", "fructose"),
    n("Glycerol-3-P", spread_x(gly_x), by - 4*s,  "Glycerol-3-P",  "entry_intermediate", "glycerol")
  )
}

# Carbon source node builder — dynamic, uses named cs_xs vector
.dnmb_cct_carbon_source_nodes <- function(cs_xs, cs_y, n) {
  cs_sugar_type <- c(
    Maltose = "glucose", Cellobiose = "glucose", Galactose = "galactose",
    Trehalose = "glucose", Lactose = "galactose", Mannose = "mannose",
    NAG = "glcnac", Glucosamine = "glcnac", Fructose = "fructose",
    Sucrose = "fructose", Mannitol = "mannose", Glycerol = "glycerol",
    Xylose = "xylose", Arabinose = "arabinose", Ribose = "ribose",
    Fucose = "fucose", Rhamnose = "rhamnose", Gluconate = "glca"
  )
  cs_ids <- names(cs_xs)
  rows <- lapply(cs_ids, function(id) {
    st <- if (id %in% names(cs_sugar_type)) cs_sugar_type[id] else "generic"
    n(id, unname(cs_xs[id]), cs_y, id, "carbon_source", unname(st))
  })
  do.call(rbind, rows)
}


# NOTE: Legacy definitions of .dnmb_cct_sugar_nodes and .dnmb_cct_3zone_cyto_edges
# were removed. Active definitions are below (~line 5100+).

# ---- 5-pointed star polygon for SNFG star shapes ----
# Returns data.frame(x, y) of 10 vertices (5 outer + 5 inner tips)

# ============================================================
# 3-ZONE VERSION: Vertical glycolysis + hierarchy + monomer hub
# ============================================================
.dnmb_snfg_normalize_sugar_type_v2 <- function(sugar_type) {
  if (length(sugar_type) == 0L) return(character(0))
  out <- as.character(sugar_type)
  na_idx <- is.na(out)
  key <- tolower(trimws(out))
  key <- sub("-p$", "", key)

  alias_map <- c(
    glcnac = "GlcNAc",
    nag = "GlcNAc",
    nacetylglucosamine = "GlcNAc",
    glca = "GlcA",
    gala = "GalA"
  )

  mapped <- unname(alias_map[key])
  replace_idx <- !is.na(mapped)
  out[replace_idx] <- mapped[replace_idx]
  out[na_idx] <- NA_character_
  out
}

 # Lookup SNFG style for a vector of sugar_type strings
 # Returns data.frame with columns: snfg_shape, snfg_color, snfg_fill
 .dnmb_snfg_lookup <- function(sugar_types) {
   ref <- .dnmb_snfg_style()
  sugar_types <- .dnmb_snfg_normalize_sugar_type_v2(sugar_types)
   idx <- match(sugar_types, ref$sugar_type)
   idx[is.na(idx)] <- match("generic", ref$sugar_type)
   ref[idx, c("snfg_shape", "snfg_color", "snfg_fill"), drop = FALSE]
 }
 
# ---------------------------------------------------------------------------
# SNFG polygon-based symbol rendering (ggplot2-only, no ggstar dependency)
# ---------------------------------------------------------------------------

#' Return polygon vertex coordinates for an SNFG shape
#' @param shape_type One of "circle", "square", "diamond", "triangle", "star"
#' @param x,y Centre coordinates
#' @param size Radius / half-width of the symbol (plot units)
#' @return data.frame with columns \code{x}, \code{y}
#' @keywords internal
.dnmb_snfg_polygon_coords_v2 <- function(shape_type, x = 0, y = 0, size = 0.075) {

  switch(shape_type,
    circle = {
      # 20-gon approximation
      theta <- seq(0, 2 * pi, length.out = 21L)[-21L]
      data.frame(x = x + size * cos(theta),
                 y = y + size * sin(theta))
    },
    square = {
      # axis-aligned square
      hs <- size  # half-side
      data.frame(x = x + c(-hs, hs, hs, -hs),
                 y = y + c(-hs, -hs, hs, hs))
    },
    diamond = {
      # 45-degree rotated square
      data.frame(x = x + size * c(0, 1, 0, -1),
                 y = y + size * c(1, 0, -1, 0))
    },
    triangle = {
      # equilateral, pointing up
      angles <- c(pi / 2, pi / 2 + 2 * pi / 3, pi / 2 + 4 * pi / 3)
      data.frame(x = x + size * cos(angles),
                 y = y + size * sin(angles))
    },
    star = {
      # 5-pointed star: 10 vertices, alternating outer / inner radius
      outer_r <- size
      inner_r <- size * 0.38  # classic pentagram ratio
      angles  <- seq(pi / 2, pi / 2 + 2 * pi, length.out = 11L)[-11L]
      radii   <- rep(c(outer_r, inner_r), 5L)
      data.frame(x = x + radii * cos(angles),
                 y = y + radii * sin(angles))
    },
    pentagon = {
      # regular pentagon, point up
      angles <- seq(pi / 2, pi / 2 + 2 * pi, length.out = 6L)[-6L]
      data.frame(x = x + size * cos(angles),
                 y = y + size * sin(angles))
    },
    # fallback: small circle
    {
      theta <- seq(0, 2 * pi, length.out = 21L)[-21L]
      data.frame(x = x + size * cos(theta),
                 y = y + size * sin(theta))
    }
  )
}

#' Map a sugar type to its SNFG shape name and colours
#' @param sugar_type character scalar (lowercase or common abbreviation)
#' @return list with elements \code{shape}, \code{fill}, \code{color}
#' @keywords internal
.dnmb_snfg_sugar_spec_v2 <- function(sugar_type) {
  # Direct polygon-shape + color lookup for SNFG rendering.
  # This is the SINGLE source of truth used by .dnmb_snfg_symbol_layers_v2()
  # and must match the SNFG symbol sheet exactly.

  # Border color is a darker shade of the fill for visual depth
  specs <- list(
    glucose      = list(shape = "circle",   fill = "#0072BC", color = "#004A7C"),
    galactose    = list(shape = "circle",   fill = "#FFD400", color = "#CCA800"),
    mannose      = list(shape = "circle",   fill = "#00A651", color = "#006B34"),
    fructose     = list(shape = "pentagon", fill = "#00A651", color = "#006B34"),
    glcnac       = list(shape = "square",   fill = "#0072BC", color = "#004A7C"),
    galnac       = list(shape = "square",   fill = "#FFD400", color = "#CCA800"),
    fucose       = list(shape = "triangle", fill = "#ED1C24", color = "#B51219"),
    rhamnose     = list(shape = "triangle", fill = "#00A651", color = "#006B34"),
    xylose       = list(shape = "star",     fill = "#F47920", color = "#C55A10"),
    arabinose    = list(shape = "star",     fill = "#00A651", color = "#006B34"),
    glca         = list(shape = "diamond",  fill = "#0072BC", color = "#004A7C", half_fill = TRUE),
    gala         = list(shape = "diamond",  fill = "#FFD400", color = "#CCA800", half_fill = TRUE),
    ribose       = list(shape = "star",     fill = "#F69EA1", color = "#C47078"),
    glycerol     = list(shape = "circle",   fill = "#3182BD", color = "#1B5A8A"),
    glucosamine  = list(shape = "square",   fill = "#0072BC", color = "#004A7C", crossed = TRUE),
    trehalose    = list(shape = "circle",   fill = "#0072BC", color = "#004A7C"),
    sucrose      = list(shape = "pentagon", fill = "#00A651", color = "#006B34"),
    maltose      = list(shape = "circle",   fill = "#0072BC", color = "#004A7C"),
    cellobiose   = list(shape = "circle",   fill = "#0072BC", color = "#004A7C"),
    lactose      = list(shape = "circle",   fill = "#FFD400", color = "#CCA800"),
    mannitol     = list(shape = "circle",   fill = "#3182BD", color = "#1B5A8A"),
    gluconate    = list(shape = "diamond",  fill = "#0072BC", color = "#004A7C", half_fill = TRUE),
    intermediate = list(shape = "circle",   fill = "#FFFFFF", color = "#000000"),
    phospho_sugar = list(shape = "circle",  fill = "#CCCCCC", color = "#888888"),
    organic_acid = list(shape = "diamond",  fill = "#E7298A", color = "#B51A6E"),
    generic      = list(shape = "circle",   fill = "#AAAAAA", color = "#666666")
  )

  base_type <- .dnmb_snfg_normalize_sugar_type_v2(sub("-[Pp]$", "", as.character(sugar_type)))
  key <- tolower(base_type)
  sp  <- specs[[key]]
  if (is.null(sp)) sp <- specs[["generic"]]
  if (is.null(sp$half_fill)) sp$half_fill <- FALSE
  sp
}

#' Build ggplot2 layers for a single SNFG symbol (polygon + optional badge)
#'
#' For half-filled diamonds (GlcA / GalA) a smaller white diamond is drawn
#' on top to mimic the SNFG convention.
#'
#' @param x,y Position in plot coordinates
#' @param sugar_type Character scalar (e.g. "glucose", "GlcNAc")
#' @param size Radius of the symbol
#' @return A list of ggplot2 layers
#' @keywords internal
.dnmb_snfg_symbol_layers_v2 <- function(x, y, sugar_type, size = 0.075, label = NULL) {

  sp    <- .dnmb_snfg_sugar_spec_v2(sugar_type)
  verts <- .dnmb_snfg_polygon_coords_v2(sp$shape, x, y, size)
  halo_verts <- .dnmb_snfg_polygon_coords_v2(sp$shape, x, y, size * 1.18)

  layers <- list(
    ggplot2::geom_polygon(
      data    = halo_verts,
      mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      fill    = "#FFFFFF",
      color   = NA,
      alpha   = 0.92,
      inherit.aes = FALSE
    )
  )

  shape_layers <- if (isTRUE(sp$half_fill) && identical(sp$shape, "diamond")) {
    .dnmb_snfg_split_diamond_layers(
      cx = x, cy = y, r = size,
      fill_color = sp$fill, border_color = sp$color, border_lw = 0.3
    )
  } else {
    list(
      ggplot2::geom_polygon(
        data    = verts,
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
        fill    = sp$fill,
        color   = sp$color,
        linewidth = 0.3,
        inherit.aes = FALSE
      )
    )
  }
  layers <- c(layers, shape_layers)

  # Abbreviation label inside the symbol
  if (!is.null(label) && nzchar(label)) {
    # White text for dark fills, black for light fills
    txt_col <- if (mean(grDevices::col2rgb(sp$fill)) < 180) "white" else "#333333"
    lbl_df <- data.frame(x = x, y = y, label = label)
    layers[[length(layers) + 1L]] <- ggplot2::geom_text(
      data = lbl_df,
      mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
      size = size * 10, fontface = "bold", color = txt_col, inherit.aes = FALSE)
  }

  # Crossed square (hexosamine): diagonal line from top-left to bottom-right
  # Upper triangle = colored, lower triangle = white
  if (isTRUE(sp$crossed)) {
    h <- size
    # White lower triangle overlay
    lower_tri <- data.frame(x = x + h * c(-1, -1, 1), y = y + h * c(1, -1, -1))
    layers[[length(layers) + 1L]] <- ggplot2::geom_polygon(
      data = lower_tri,
      mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      fill = "white", color = NA, inherit.aes = FALSE)
    # Diagonal line
    diag_df <- data.frame(x = x - h, xend = x + h, y = y + h, yend = y - h)
    layers[[length(layers) + 1L]] <- ggplot2::geom_segment(
      data = diag_df,
      mapping = ggplot2::aes(x = .data[["x"]], xend = .data[["xend"]],
                              y = .data[["y"]], yend = .data[["yend"]]),
      linewidth = 0.3, color = sp$color, inherit.aes = FALSE)
  }

  # Phosphorylation badge (optional, triggered by a trailing "P")
  if (grepl("-P$", sugar_type)) {
    badge_df <- data.frame(x = x + size * 0.75, y = y + size * 0.75,
                           label = "P")
    layers[[length(layers) + 1L]] <- ggplot2::geom_text(
      data    = badge_df,
      mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                             label = .data[["label"]]),
      size    = size * 12,
      fontface = "bold",
      color   = "#333333",
      inherit.aes = FALSE
    )
  }

  layers
}

#' Add an SNFG symbol to a ggplot — unified rendering entry point
#'
#' This is the single function that both the transport map and the symbol sheet
#' should call.  It resolves sugar_name → shape + color via
#' \code{.dnmb_snfg_sugar_spec_v2()} and appends polygon layers to \code{p}.
#'
#' @param p A ggplot2 object
#' @param cx,cy Centre coordinates
#' @param sugar_name Character scalar (e.g. "glucose", "GlcNAc", "fructose")
#' @param r Radius of the symbol (plot units)
#' @return The modified ggplot2 object with symbol layers added
#' @keywords internal
.dnmb_snfg_render_symbol_v2 <- function(p, cx, cy, sugar_name, r = 0.075, label = NULL) {
  layers <- .dnmb_snfg_symbol_layers_v2(cx, cy, sugar_name, size = r, label = label)
  for (ly in layers) p <- p + ly
  p
}

# Sugar abbreviation lookup
.dnmb_snfg_abbreviation_v2 <- function(sugar_type) {
  abbr <- c(
    glucose = "Glc", galactose = "Gal", mannose = "Man", fructose = "Fru",
    glcnac = "NAc", galnac = "NAc", fucose = "Fuc", rhamnose = "Rha",
    xylose = "Xyl", arabinose = "Ara", glca = "GA", gala = "GA",
    ribose = "Rib", glycerol = "Gly", glucosamine = "GlN",
    trehalose = "Tre", sucrose = "Suc", maltose = "Mal",
    cellobiose = "Cel", lactose = "Lac", mannitol = "Mtl",
    gluconate = "Gnt", intermediate = "", organic_acid = "",
    phospho_sugar = "", generic = "")
  key <- tolower(sugar_type)
  if (key %in% names(abbr)) unname(abbr[key]) else ""
}

.dnmb_cct_canonical_sugar_id <- function(sugar_type) {
  if (length(sugar_type) == 0L) return(character(0))
  raw <- tolower(trimws(as.character(sugar_type)))
  norm <- .dnmb_snfg_normalize_sugar_type_v2(sugar_type)
  out <- tolower(trimws(as.character(norm)))
  bad <- is.na(out) | !nzchar(out)
  out[bad] <- raw[bad]
  out
}

.dnmb_cct_sugar_route_color <- function(sugar_type, fallback = "#666666") {
  if (is.null(sugar_type) || length(sugar_type) == 0 || is.na(sugar_type) || !nzchar(sugar_type)) {
    return(fallback)
  }
  sp <- .dnmb_snfg_sugar_spec_v2(tolower(as.character(sugar_type)))
  if (!is.null(sp$color) && nzchar(sp$color)) sp$color else fallback
}

.dnmb_cct_reference_gray <- function(strength = c("light", "mid", "dark")) {
  strength <- match.arg(strength)
  switch(strength,
    light = "#D9DDE1",
    mid = "#C8CDD3",
    dark = "#B4BBC3"
  )
}

.dnmb_cct_short_step_label <- function(step_id, pathway_id = NA_character_) {
  x <- as.character(step_id)[1]
  if (is.na(x) || !nzchar(x)) return("")
  pw <- tolower(as.character(pathway_id)[1])
  alias_map <- c(
    nagA = "NagA", nagB = "NagB", galK = "GalK", galE = "GalE", galT = "GalT",
    pgmA = "PgmA", glk = "Glk", gnd = "Gnd", manA = "ManA", gntK = "GntK",
    kdgK = "KdgK", garK = "GarK", garL = "GarL", garR = "GarR", gci = "Gci",
    udh = "Udh", glpK = "GlpK", glpO = "GlpO", tpi = "Tpi", fba = "Fba",
    lacE = "LacE", lacF = "LacF", lacG = "LacG", lacK = "LacK", lacZ = "LacZ",
    malE1 = "MalE", malF1 = "MalF", malG1 = "MalG", malK1 = "MalK", susB = "SusB",
    thuE = "ThuE", thuF = "ThuF", thuG = "ThuG", thuK = "ThuK", PsTP = "Tre PTS",
    rbsA = "RbsA", rbsB = "RbsB", rbsC = "RbsC", rbsK = "RbsK",
    xylB = "XylB", xylF = "XylF", xylG = "XylG", xylH = "XylH",
    araA = "AraA", araB = "AraB", araD = "AraD", chvE = "Ara SBP",
    mglA = "MglA", mglB = "MglB", mglC = "MglC", srlD = "SrlD",
    ackA = "AckA", pta = "Pta", adh = "Adh", eda = "Eda", `L-LDH` = "L-LDH",
    `D-LDH` = "D-LDH", `1pfk` = "1-PFK", gguA = "GguA", gguB = "GguB",
    `deoxyribonate-dehyd` = "DR dehyd", `deoxyribonate-transport` = "DR transport",
    `ketodeoxyribonate-cleavage` = "DR cleavage", TM0027 = "PTS-A", TM0028 = "PTS-B",
    TM0029 = "PTS-C", TM0030 = "PTS-D", TM0031 = "Bgl"
  )
  if (x %in% names(alias_map)) return(unname(alias_map[x]))
  if (identical(pw, "trehalose") && x == "glk") return("Tre Glk")
  if (identical(pw, "sucrose") && x == "glk") return("ScrK/Glk")
  x <- sub("[-_](transport|utilization).*$", "", x, ignore.case = TRUE)
  x <- sub("\\s+(transport|utilization).*$", "", x, ignore.case = TRUE)
  x <- sub("\\s*\\(.*\\)$", "", x)
  x <- sub("/.*$", "", x)
  x <- sub("\\|.*$", "", x)
  x <- gsub("_", "-", x, fixed = TRUE)
  if (nchar(x) > 18) {
    parts <- unlist(strsplit(x, "[- ]+"))
    if (length(parts) > 0 && nzchar(parts[1])) x <- parts[1]
  }
  x
}

.dnmb_cct_is_generic_feature_id <- function(x) {
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) return(TRUE)
  grepl("^TM\\d+$", x, ignore.case = TRUE) ||
    grepl("^[A-Z0-9]+_RS\\d+$", x) ||
    grepl("^gene[_-]?\\d+$", x, ignore.case = TRUE) ||
    grepl("^locus[_-]?\\d+$", x, ignore.case = TRUE)
}

.dnmb_cct_normalize_locus_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("[\u2212\u2010\u2011\u2012\u2013\u2014\u2015]", "-", x, perl = TRUE)
  x <- gsub("\\s+", "", x, perl = TRUE)
  x
}

.dnmb_cct_filter_valid_locus_tags <- function(x, valid_tags) {
  x <- .dnmb_cct_normalize_locus_tag(x)
  valid_tags <- unique(.dnmb_cct_normalize_locus_tag(valid_tags))
  valid_tags <- valid_tags[!is.na(valid_tags) & nzchar(valid_tags)]
  out <- vapply(x, function(one) {
    if (is.na(one) || !nzchar(one)) return("")
    parts <- trimws(unlist(strsplit(one, ",", fixed = TRUE)))
    parts <- parts[parts %in% valid_tags]
    if (!length(parts)) return("")
    paste(unique(parts), collapse = ",")
  }, character(1))
  out
}

.dnmb_cct_prefer_human_label <- function(step_id, gene_name = NA_character_, pathway_id = NA_character_) {
  gn <- as.character(gene_name)[1]
  if (!is.na(gn) && nzchar(gn) && !.dnmb_cct_is_generic_feature_id(gn)) return(gn)
  .dnmb_cct_short_step_label(step_id, pathway_id = pathway_id)
}

.dnmb_cct_is_transport_like_step <- function(step_id, gene_name = NA_character_, is_pts = FALSE) {
  if (isTRUE(is_pts)) return(TRUE)
  # Generic keywords (match anywhere in step_id or gene_name)
  generic_re <- "transport|permease|porter|uptake|symporter|antiporter|binding.protein|\\bsbp\\b|\\babc\\b|\\bmfs\\b|\\bpts\\b|\\beiic|\\beiib|\\beiia"
  # Specific gene names (use word boundaries, check each field separately)
  gene_re <- "\\b(mgl[ABC]|msiK|malE|malF|malG|malK|thu[EFGK]|rbs[ABC]|xyl[FGH]|ggu[AB]|glc[TUVP]|galP|lac[PEFG]|chvE|gts[ABCD]|kguT|ara[ESTUV]|mtl[AEK]|fucP|nagE|nagF|nagPcb|fruA|fruB|fruD|scrT|sacP|cscB|treB|manP|bglF|celB|ascB|ptsG|ptsH|ptsI|crr)\\b"
  sid <- as.character(step_id)[1]
  gnm <- as.character(gene_name)[1]
  if (grepl(generic_re, sid, ignore.case = TRUE)) return(TRUE)
  if (!is.na(gnm) && grepl(generic_re, gnm, ignore.case = TRUE)) return(TRUE)
  if (grepl(gene_re, sid, ignore.case = TRUE)) return(TRUE)
  if (!is.na(gnm) && grepl(gene_re, gnm, ignore.case = TRUE)) return(TRUE)
  FALSE
}

.dnmb_cct_transporter_half_span <- function(step_id, gene_name = NA_character_, is_pts = FALSE) {
  txt <- paste(step_id, gene_name)
  if (isTRUE(is_pts)) return(0.34)
  if (grepl("abc|sbp|binding|malE|mgl[ABC]|rbs[ABC]|xyl[FGH]|thu[EFG]", txt, ignore.case = TRUE)) return(0.30)
  if (grepl("mfs|permease|symporter|transporter|porter|gnt|lac|gal|man|srl|iolT|kgtP", txt, ignore.case = TRUE)) return(0.26)
  0.22
}

.dnmb_cct_transporter_pref_regex <- function(cs_id) {
  switch(as.character(cs_id),
    Maltose = "mal|malt|sus|agl|msiK",
    Cellobiose = "cbt|ceb|cbp|cel|bgl|asc|mgl|agl|msiK|TM00",
    Galactose = "gal|ggu|chv|mgl|glc",
    Trehalose = "thu|tre|mal|agl|mgl|pts",
    Lactose = "lac|gal|ggu|mgl",
    Mannose = "man|glc|mfs|pts",
    NAG = "nag|pts|glc|ngc",
    Glucosamine = "nag|pts|glc|ngc",
    Fructose = "fru|pts|glc|ara",
    Sucrose = "scr|sac|csc|mgl|agl|pts",
    Mannitol = "mtl|man|pts",
    Glycerol = "glp|gyl",
    Xylose = "xyl|ara|ggu|mgl",
    Arabinose = "ara|ggu|xyl|chv",
    Ribose = "rbs|nup|plt|mgl",
    Fucose = "fuc|plt|pts",
    Gluconate = "gnt|kgu|mfs|glu",
    ""
  )
}

.dnmb_cct_pathway_context_regex <- function(pathway_id) {
  switch(tolower(as.character(pathway_id)[1]),
    maltose = "maltose|starch|alpha[ -]?gluc|amylo|pullulan",
    cellobiose = "cellobiose|cellulose|beta[ -]?gluc|cellodextr",
    galactose = "galactose|lactose|raffinose|galactan|beta-gal",
    trehalose = "trehalose|trehalose-6|alpha,alpha-trehalose",
    lactose = "lactose|galactose|beta-gal|galactan",
    mannose = "mannose|mannan|galactomannan|beta-mannan",
    mannitol = "mannitol|mannose|mannan|galactomannan|sorbitol",
    nag = "chitin|glcnac|n-acetylglucosamine|nag|peptidoglycan",
    glucosamine = "glucosamine|glcnac|chitin|nag",
    fructose = "fructose|fructan|inulin|levan|sucrose",
    sucrose = "sucrose|fructose|fructan|inulin|levan",
    glycerol = "glycerol|glycerol-3-phosphate|glp",
    xylose = "xylose|xylan|xylan",
    arabinose = "arabinose|arabinan|xylan|xylan",
    ribose = "ribose|ribonucleoside",
    fucose = "fucose|fuculose",
    rhamnose = "rhamnose|rhamnulose",
    gluconate = "gluconate|glucuronate|galacturonate|pectin|pectate|gala|glca|kdg",
    ""
  )
}

.dnmb_cct_pathway_dbcan_families <- function(pathway_id) {
  pid <- tolower(as.character(pathway_id)[1])
  fams <- switch(pid,
    maltose = c("GH13", "GH13_11", "GH13_20", "GH13_31", "GH13_36", "GH13_48", "GH4"),
    cellobiose = c("GH1", "GH3", "GH5", "GH6", "GH7", "GH9", "GH44", "GH45", "GH48"),
    galactose = c("GH35", "GH36", "GH42", "GH2"),
    trehalose = c("GH37", "GH65", "GH4"),
    lactose = c("GH35", "GH42", "GH2", "GH1"),
    mannose = c("GH26", "GH113", "GH92", "GH130", "GH130_2", "GH130_6", "GH130_4", "GH27", "GH36"),
    mannitol = c("GH26", "GH113", "GH92", "GH130", "GH130_2", "GH130_6", "GH130_4", "GH27", "GH36"),
    nag = c("GH18", "GH19", "GH20", "GH3"),
    glucosamine = c("GH18", "GH19", "GH20", "GH3"),
    fructose = c("GH32", "GH68", "GH91"),
    sucrose = c("GH32", "GH68", "GH91"),
    glycerol = character(0),
    xylose = c("GH10", "GH11", "GH30", "GH43", "GH8"),
    arabinose = c("GH43", "GH51", "GH54", "GH62"),
    ribose = character(0),
    fucose = c("GH29", "GH95", "GH151"),
    gluconate = c("GH28", "GH78", "GH106", "GH35", "GH42", "GH88", "GH105", "GH53"),
    character(0)
  )
  unique(fams)
}

.dnmb_cct_transporter_kind <- function(step_id, gene_name = NA_character_, is_pts = FALSE) {
  txt <- paste(step_id, gene_name)
  if (isTRUE(is_pts)) return("PTS")
  if (grepl("abc|sbp|binding|malE|mgl[ABC]|rbs[ABC]|xyl[FGH]|thu[EFG]|pot[ABCD]", txt, ignore.case = TRUE)) {
    return("ABC")
  }
  if (grepl("mfs|permease|symporter|transporter|porter|gnt|lac|gal|man|srl|iolT|kgtP|thuK|glpO", txt, ignore.case = TRUE)) {
    return("MFS")
  }
  "GEN"
}

.dnmb_cct_transporter_display_score <- function(cs_id, step_id, gene_name,
                                                confidence = "low", score = NA_real_,
                                                context_score = 0,
                                                is_pts = FALSE) {
  txt <- paste(step_id, gene_name)
  pref_re <- .dnmb_cct_transporter_pref_regex(cs_id)
  conf_score <- c(none = 0, low = 1, medium = 2, high = 3)[as.character(confidence)]
  conf_score[is.na(conf_score)] <- 0
  kind <- .dnmb_cct_transporter_kind(step_id, gene_name, is_pts = is_pts)
  kind_bonus <- c(PTS = 1.4, ABC = 1.1, MFS = 0.9, GEN = 0.4)[kind]
  path_bonus <- if (nzchar(pref_re) && grepl(pref_re, txt, ignore.case = TRUE)) 2.4 else 0
  generic_bonus <- if (grepl("transport|permease|porter|symporter|abc|binding|pts|mfs|tm00", txt, ignore.case = TRUE)) 0.8 else 0
  name_bonus <- if (!is.na(gene_name) && nzchar(gene_name) && gene_name != step_id) 0.15 else 0
  numeric_bonus <- ifelse(is.na(score), 0, pmax(0, score) * 0.35)
  conf_score * 2.0 + kind_bonus + path_bonus + generic_bonus + name_bonus + numeric_bonus + context_score
}

.dnmb_cct_select_transporters <- function(transporters, max_per_lane = 4L) {
  if (is.null(transporters) || nrow(transporters) == 0) return(transporters)
  transporters$display_score <- mapply(
    .dnmb_cct_transporter_display_score,
    cs_id = transporters$cs_id,
    step_id = transporters$step,
    gene_name = transporters$gene_name,
    confidence = transporters$confidence,
    score = transporters$step_score,
    context_score = transporters$context_score,
    is_pts = transporters$is_pts
  )
  out_idx <- integer(0)
  split_idx <- split(seq_len(nrow(transporters)), transporters$cs_id)
  for (grp in split_idx) {
    sub <- transporters[grp, , drop = FALSE]
    sub <- sub[order(-sub$display_score, -sub$conf_rank, sub$gene_name, sub$locus_tag), , drop = FALSE]
    keep <- integer(0)
    seen_kind <- character(0)
    for (i in seq_len(nrow(sub))) {
      if (length(keep) >= max_per_lane) break
      k <- sub$kind[i]
      if (!k %in% seen_kind || sub$display_score[i] >= median(sub$display_score, na.rm = TRUE)) {
        keep <- c(keep, i)
        seen_kind <- c(seen_kind, k)
      }
    }
    if (length(keep) < min(max_per_lane, nrow(sub))) {
      extra <- setdiff(seq_len(nrow(sub)), keep)
      keep <- c(keep, head(extra, min(max_per_lane, nrow(sub)) - length(keep)))
    }
    out_idx <- c(out_idx, grp[keep])
  }
  out <- transporters[sort(unique(out_idx)), , drop = FALSE]
  out <- out[order(out$cs_id, -out$display_score, -out$conf_rank, out$gene_name, out$locus_tag), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_cct_annotate_transport_context <- function(transporters, genbank_table, window = 10L) {
  if (is.null(transporters) || nrow(transporters) == 0) return(transporters)
  gt <- as.data.frame(genbank_table, stringsAsFactors = FALSE, row.names = NULL)
  req <- c("locus_tag", "contig", "start")
  if (!all(req %in% names(gt))) {
    transporters$context_score <- 0
    transporters$context_hits <- NA_character_
    return(transporters)
  }

  gt <- gt[order(gt$contig, as.numeric(gt$start)), , drop = FALSE]
  gt$row_idx <- seq_len(nrow(gt))
  lt_idx <- match(transporters$locus_tag, gt$locus_tag)
  cazy_map <- .dnmb_cct_cazy_substrate_map()
  fam_to_sub <- split(cazy_map$substrate, cazy_map$family)
  if (!"direction" %in% names(gt) && "strand" %in% names(gt)) gt$direction <- gt$strand
  if (!"end" %in% names(gt)) gt$end <- gt$start

  transporters$context_score <- 0
  transporters$context_hits <- NA_character_

  for (i in seq_len(nrow(transporters))) {
    gi <- lt_idx[i]
    if (is.na(gi)) next
    contig_i <- gt$contig[gi]
    lo <- max(1, gi - window)
    hi <- min(nrow(gt), gi + window)
    nb <- gt[lo:hi, , drop = FALSE]
    nb <- nb[nb$contig == contig_i & nb$locus_tag != transporters$locus_tag[i], , drop = FALSE]
    if (nrow(nb) == 0) next

    rx <- .dnmb_cct_pathway_context_regex(transporters$pathway[i])
    fam_pref <- .dnmb_cct_pathway_dbcan_families(transporters$cs_id[i])
    score <- 0
    hits <- character(0)
    tr_dir <- if ("direction" %in% names(gt)) as.character(gt$direction[gi]) else NA_character_
    tr_end <- if ("end" %in% names(gt)) suppressWarnings(as.numeric(gt$end[gi])) else NA_real_
    tr_start <- suppressWarnings(as.numeric(gt$start[gi]))

    same_strand_idx <- if ("direction" %in% names(nb) && !is.na(tr_dir) && nzchar(tr_dir)) {
      !is.na(nb$direction) & as.character(nb$direction) == tr_dir
    } else {
      rep(FALSE, nrow(nb))
    }
    near_bp_idx <- if (!is.na(tr_start)) {
      nb_start <- suppressWarnings(as.numeric(nb$start))
      nb_end <- suppressWarnings(as.numeric(nb$end))
      d1 <- abs(nb_start - tr_end)
      d2 <- abs(nb_end - tr_start)
      pmin(d1, d2) <= 250
    } else {
      rep(FALSE, nrow(nb))
    }
    operon_idx <- same_strand_idx & near_bp_idx

    if ("dbCAN_family_id" %in% names(nb)) {
      fams <- as.character(nb$dbCAN_family_id)
      fams <- fams[!is.na(fams) & nzchar(fams)]
      if (length(fams) > 0) {
        for (ff in fams) {
          subs <- fam_to_sub[[ff]]
          if (length(subs) > 0 && any(grepl(rx, tolower(subs), ignore.case = TRUE))) {
            score <- score + 1.4
            hits <- c(hits, paste0("dbCAN:", ff))
          }
          if (length(fam_pref) > 0 && ff %in% fam_pref) {
            score <- score + 1.8
            hits <- c(hits, paste0("prefGH:", ff))
          }
        }
      }

      if (any(operon_idx)) {
        fams_op <- as.character(nb$dbCAN_family_id[operon_idx])
        fams_op <- fams_op[!is.na(fams_op) & nzchar(fams_op)]
        if (length(fams_op) > 0) {
          score <- score + 0.9
          hits <- c(hits, paste0("operon-dbCAN:", fams_op[1]))
          if (length(fam_pref) > 0 && any(fams_op %in% fam_pref)) {
            score <- score + 1.2
            hits <- c(hits, paste0("operon-prefGH:", fams_op[fams_op %in% fam_pref][1]))
          }
        }
      }
    }

    if ("PAZy_substrate_label" %in% names(nb)) {
      pz <- tolower(as.character(nb$PAZy_substrate_label))
      pz <- pz[!is.na(pz) & nzchar(pz)]
      if (length(pz) > 0) {
        match_idx <- grepl(rx, pz, ignore.case = TRUE)
        if (any(match_idx)) {
          score <- score + 1.2 * sum(match_idx)
          hits <- c(hits, paste0("PAZy:", unique(pz[match_idx])[1]))
        }
      }
    }

    prod_txt <- tolower(paste(
      if ("product" %in% names(nb)) nb$product else "",
      if ("gene" %in% names(nb)) nb$gene else "",
      collapse = " | "
    ))
    if (nzchar(prod_txt) && grepl(rx, prod_txt, ignore.case = TRUE)) {
      score <- score + 0.8
      hits <- c(hits, "neighbor-text")
    }

    if (any(operon_idx)) {
      op_txt <- tolower(paste(
        if ("product" %in% names(nb)) nb$product[operon_idx] else "",
        if ("gene" %in% names(nb)) nb$gene[operon_idx] else "",
        collapse = " | "
      ))
      if (nzchar(op_txt) && grepl(rx, op_txt, ignore.case = TRUE)) {
        score <- score + 1.0
        hits <- c(hits, "operon-text")
      }
    }

    if ("dbCAN_dbcan_cgc_id" %in% names(nb)) {
      cgc_vals <- as.character(nb$dbCAN_dbcan_cgc_id)
      cgc_vals <- cgc_vals[!is.na(cgc_vals) & nzchar(cgc_vals)]
      if (length(cgc_vals) > 0) {
        score <- score + 0.6
        hits <- c(hits, paste0("cgc=", unique(cgc_vals)[1]))
      }
    }
    if ("dbCAN_dbcan_cgc_gene_type" %in% names(nb)) {
      cgc_types <- tolower(as.character(nb$dbCAN_dbcan_cgc_gene_type))
      cgc_types <- cgc_types[!is.na(cgc_types) & nzchar(cgc_types)]
      if (length(cgc_types) > 0) {
        if (any(cgc_types %in% c("tc", "tp", "transport", "transporter", "stp"))) {
          score <- score + 0.8
          hits <- c(hits, "cgc-type=TC")
        }
        if (any(cgc_types %in% c("caZyme", "cazyme"))) {
          score <- score + 0.6
          hits <- c(hits, "cgc-type=CAZyme")
        }
      }
    }

    transporters$context_score[i] <- score
    if (length(hits) > 0) {
      transporters$context_hits[i] <- paste(unique(hits), collapse = "; ")
    }
  }

  transporters
}

.dnmb_cct_transporter_glyph_layers <- function(tx, ty, half_span, core_color,
                                               confidence = "medium", is_pts = FALSE,
                                               kind = "GEN") {
  conf_lw <- switch(confidence, high = 1.05, medium = 0.92, low = 0.78, 0.72)
  conf_alpha <- switch(confidence, high = 0.95, medium = 0.82, low = 0.66, 0.55)
  core_lty <- if (isTRUE(is_pts)) "dashed" else "solid"

  shell_df <- data.frame(x = tx - half_span, xend = tx + half_span, y = ty, yend = ty)
  layers <- list(
    ggplot2::geom_segment(
      data = shell_df,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = conf_lw + 0.9, color = "#7F8C8D", alpha = 0.55,
      lineend = "round", inherit.aes = FALSE
    ),
    ggplot2::geom_segment(
      data = shell_df,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = conf_lw + 0.45, color = "#F8F8F8", alpha = 0.95,
      lineend = "round", inherit.aes = FALSE
    ),
    ggplot2::geom_segment(
      data = shell_df,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = conf_lw, color = core_color, alpha = conf_alpha,
      linetype = core_lty, lineend = "round", inherit.aes = FALSE
    )
  )

  if (identical(kind, "PTS")) {
    center_df <- data.frame(x = tx, y = ty)
    layers[[length(layers) + 1L]] <- ggplot2::geom_point(
      data = center_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 23, size = 2.6, fill = core_color, color = "#FFFFFF",
      stroke = 0.35, alpha = conf_alpha, inherit.aes = FALSE
    )
  } else if (identical(kind, "ABC")) {
    cap_df <- data.frame(
      x = c(tx - half_span * 0.72, tx + half_span * 0.72),
      y = c(ty, ty)
    )
    layers[[length(layers) + 1L]] <- ggplot2::geom_point(
      data = cap_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 2.4, fill = core_color, color = "#FFFFFF",
      stroke = 0.28, alpha = conf_alpha, inherit.aes = FALSE
    )
    gate_df <- data.frame(x = tx - half_span * 0.18, xend = tx + half_span * 0.18, y = ty, yend = ty)
    layers[[length(layers) + 1L]] <- ggplot2::geom_segment(
      data = gate_df,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 0.26, color = "#FFFFFF", alpha = 0.9, inherit.aes = FALSE
    )
  } else if (identical(kind, "MFS")) {
    tilt_df <- data.frame(
      x = tx - half_span * 0.22, xend = tx + half_span * 0.22,
      y = ty - 0.05, yend = ty + 0.05
    )
    layers[[length(layers) + 1L]] <- ggplot2::geom_segment(
      data = tilt_df,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 0.32, color = "#FFFFFF", alpha = 0.85, inherit.aes = FALSE
    )
  }

  layers
}

.dnmb_cct_pack_transporters_lane <- function(center_x, half_spans, lane_ranks,
                                            label_widths = NULL,
                                            label_dx = NULL,
                                            desired_x = NULL) {
  n <- length(half_spans)
  if (n == 0) return(data.frame(tx = numeric(0), ty = numeric(0), row_id = integer(0)))
  if (n == 1) return(data.frame(tx = center_x, ty = 8.50, row_id = 1L))

  priorities <- if (missing(lane_ranks) || length(lane_ranks) != n) rep(1, n) else lane_ranks
  priorities <- rank(-priorities, ties.method = "first")
  weight <- rev(priorities) + 1
  if (is.null(label_widths) || length(label_widths) != n) label_widths <- rep(0.4, n)
  if (is.null(label_dx) || length(label_dx) != n) label_dx <- rep(0, n)
  if (is.null(desired_x) || length(desired_x) != n) desired_x <- rep(center_x, n)

  row_patterns <- switch(as.character(n),
    `2` = list(c(2L), c(1L, 1L)),
    `3` = list(c(2L, 1L), c(1L, 2L), c(1L, 1L, 1L)),
    `4` = list(c(2L, 2L), c(1L, 2L, 1L), c(2L, 1L, 1L), c(1L, 1L, 2L), c(1L, 1L, 1L, 1L)),
    list(rep(1L, n))
  )
  perms <- .dnmb_cct_small_permutations(seq_len(n))

  pack_rows <- function(row_groups) {
    y_levels <- c(8.56, 8.44, 8.32)
    tx <- numeric(n)
    ty <- numeric(n)
    row_id <- integer(n)
    gap <- 0.08
    for (rr in seq_along(row_groups)) {
      idx <- row_groups[[rr]]
      spans <- half_spans[idx]
      widths <- 2 * spans
      total_w <- sum(widths) + gap * max(0, length(idx) - 1L)
      left <- center_x - total_w / 2
      cur <- left
      for (j in seq_along(idx)) {
        ii <- idx[j]
        tx[ii] <- cur + widths[j] / 2
        ty[ii] <- y_levels[rr]
        row_id[ii] <- rr
        cur <- tx[ii] + widths[j] / 2 + gap
      }
    }
    data.frame(tx = tx, ty = ty, row_id = row_id)
  }

  best <- NULL
  best_cost <- Inf
  for (pat in row_patterns) {
    for (perm in perms) {
      start <- 1L
      row_groups <- list()
      for (ps in pat) {
        row_groups[[length(row_groups) + 1L]] <- perm[start:(start + ps - 1L)]
        start <- start + ps
      }
      cand <- pack_rows(row_groups)
      span_cost <- (max(cand$tx) - min(cand$tx)) * 2.5
      center_cost <- sum(abs(cand$tx - desired_x) * weight)
      row_cost <- sum((cand$row_id - 1L) * weight * 1.15)
      low_row_penalty <- sum((cand$row_id > 1L) * weight * 0.35)
      label_cost <- 0
      for (rr in unique(cand$row_id)) {
        idx <- which(cand$row_id == rr)
        if (length(idx) >= 2) {
          label_x <- cand$tx[idx] + label_dx[idx]
          ord2 <- idx[order(label_x)]
          for (k in 2:length(ord2)) {
            i1 <- ord2[k - 1L]; i2 <- ord2[k]
            required_gap <- (label_widths[i1] + label_widths[i2]) / 2 + 0.10
            actual_gap <- abs((cand$tx[i2] + label_dx[i2]) - (cand$tx[i1] + label_dx[i1]))
            if (actual_gap < required_gap) {
              label_cost <- label_cost + (required_gap - actual_gap) * 8
            }
          }
        }
      }
      vertical_penalty <- max(cand$row_id) * 0.6
      cost <- span_cost + center_cost + row_cost + low_row_penalty + label_cost + vertical_penalty
      if (cost < best_cost) {
        best <- cand
        best_cost <- cost
      }
    }
  }
  best
}

.dnmb_cct_junction_glyph_layers <- function(x, y, color, size = 2.1, alpha = 0.9) {
  list(
    ggplot2::geom_point(
      data = data.frame(x = x, y = y),
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = size, fill = "#FFFFFF", color = color,
      stroke = 0.45, alpha = alpha, inherit.aes = FALSE
    ),
    ggplot2::geom_point(
      data = data.frame(x = x, y = y),
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = size * 0.34, color = color,
      alpha = alpha, inherit.aes = FALSE
    )
  )
}

.dnmb_cct_place_small_labels <- function(df, x_thresh = 0.7, y_thresh = 0.18) {
  if (is.null(df) || nrow(df) == 0) return(df)
  cand <- data.frame(
    dx = c(0, 0, 0, 0.18, -0.18, 0.26, -0.26, 0.34, -0.34),
    dy = c(0, 0.16, -0.16, 0.10, 0.10, -0.10, -0.10, 0.22, 0.22)
  )
  base_x <- if ("x_lab" %in% names(df)) df$x_lab else df$x
  base_y <- if ("y_lab" %in% names(df)) df$y_lab else df$y
  ord <- order(base_x, base_y)
  out <- df[ord, , drop = FALSE]
  base_x <- base_x[ord]
  base_y <- base_y[ord]
  placed_x <- numeric(0)
  placed_y <- numeric(0)
  out$x_lab <- base_x
  out$y_lab <- base_y
  out$hjust <- 0.5

  for (i in seq_len(nrow(out))) {
    chosen <- NULL
    for (j in seq_len(nrow(cand))) {
      cx <- base_x[i] + cand$dx[j]
      cy <- base_y[i] + cand$dy[j]
      clash <- any(abs(cx - placed_x) < x_thresh & abs(cy - placed_y) < y_thresh)
      if (!clash) {
        chosen <- c(cx, cy)
        out$hjust[i] <- if (cand$dx[j] > 0.05) 0 else if (cand$dx[j] < -0.05) 1 else 0.5
        break
      }
    }
    if (is.null(chosen)) {
      chosen <- c(base_x[i], base_y[i] + 0.28)
      out$hjust[i] <- 0.5
    }
    out$x_lab[i] <- chosen[1]
    out$y_lab[i] <- chosen[2]
    placed_x <- c(placed_x, chosen[1])
    placed_y <- c(placed_y, chosen[2])
  }
  out
}

.dnmb_cct_lane_label_offsets <- function(n, lane_dir = 0) {
  lane_dir <- sign(as.numeric(lane_dir)[1])
  if (!is.finite(lane_dir)) lane_dir <- 0
  if (lane_dir > 0) {
    base <- c(0.14, 0.24, 0.34, 0.46, 0.58, 0.70)
  } else if (lane_dir < 0) {
    base <- c(-0.14, -0.24, -0.34, -0.46, -0.58, -0.70)
  } else {
    base <- c(0, 0.18, -0.18, 0.30, -0.30, 0.42, -0.42)
  }
  if (n <= length(base)) return(base[seq_len(n)])
  c(base, rep(base[length(base)], n - length(base)))
}

.dnmb_cct_thin_labels <- function(df, x_thresh = 0.75, y_thresh = 0.2,
                                  priority_col = "priority") {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!priority_col %in% names(df)) df[[priority_col]] <- seq_len(nrow(df))
  ord <- order(df[[priority_col]], decreasing = TRUE, na.last = TRUE)
  keep <- logical(nrow(df))
  kept_x <- numeric(0)
  kept_y <- numeric(0)
  for (ii in ord) {
    clash <- any(abs(df$x[ii] - kept_x) < x_thresh & abs(df$y[ii] - kept_y) < y_thresh)
    if (!clash) {
      keep[ii] <- TRUE
      kept_x <- c(kept_x, df$x[ii])
      kept_y <- c(kept_y, df$y[ii])
    }
  }
  df[keep, , drop = FALSE]
}

.dnmb_cct_path_tangent_angle <- function(points_df, idx = NULL) {
  if (is.null(points_df) || nrow(points_df) < 2) return(0)
  if (is.null(idx) || !is.finite(idx)) idx <- max(1, round(nrow(points_df) / 2))
  i0 <- max(1, idx - 1L)
  i1 <- min(nrow(points_df), idx + 1L)
  if (i0 == i1) return(0)
  dx <- points_df$x[i1] - points_df$x[i0]
  dy <- points_df$y[i1] - points_df$y[i0]
  ang <- atan2(dy, dx) * 180 / pi
  if (ang > 90) ang <- ang - 180
  if (ang < -90) ang <- ang + 180
  ang
}

.dnmb_cct_hub_entry_path <- function(hub_x, hub_y, target_x, target_y,
                                     lane_rank = 1L,
                                     target_count = 1L,
                                     grid_step = 0.5) {
  dx <- target_x - hub_x
  dy <- target_y - hub_y
  if (abs(dx) <= max(0.30, grid_step * 0.6) || abs(dy) <= 0.02) {
    return(data.frame(x = c(hub_x, target_x), y = c(hub_y, target_y)))
  }

  lane_rank <- max(1L, as.integer(lane_rank)[1])
  side <- sign(hub_x - target_x)
  if (!is.finite(side) || side == 0) side <- ifelse(lane_rank %% 2L == 1L, -1, 1)
  stagger_seq <- c(-0.25, 0.25, -0.50, 0.50, -0.75, 0.75)
  approach_seq <- c(0.5, 1.0, 0.75, 1.25, 0.25, 1.5)
  stagger <- stagger_seq[((lane_rank - 1L) %% length(stagger_seq)) + 1L] * grid_step
  spread_scale <- 1 + 0.22 * max(0, as.integer(target_count)[1] - 1L)
  approach_offset <- approach_seq[((lane_rank - 1L) %% length(approach_seq)) + 1L] * grid_step * spread_scale
  if (abs(dx) < approach_offset + grid_step * 0.45) {
    approach_offset <- max(grid_step * 0.5, abs(dx) - grid_step * 0.25)
  }
  corridor_hi <- hub_y - grid_step * 0.5
  corridor_lo <- target_y + grid_step * 1.1
  corridor_y <- .dnmb_cct_snap_to_grid(
    max(corridor_lo, min(corridor_hi, target_y + grid_step * (1.0 + 0.2 * max(0, target_count - 1L)) + stagger)),
    step = grid_step / 2
  )
  approach_y <- .dnmb_cct_snap_to_grid(
    target_y + grid_step * (0.85 + 0.45 * ((lane_rank - 1L) %% 4L)),
    step = grid_step / 2
  )
  approach_x <- .dnmb_cct_snap_to_grid(target_x + side * approach_offset, step = grid_step / 2)
  if (abs(approach_x - hub_x) < grid_step * 0.4) {
    approach_x <- .dnmb_cct_snap_to_grid(target_x + side * max(grid_step, approach_offset), step = grid_step / 2)
  }
  slot_offset <- grid_step * (0.20 + 0.14 * min(3L, max(0L, target_count - 1L)) + 0.10 * ((lane_rank - 1L) %% 3L))
  slot_x <- .dnmb_cct_snap_to_grid(target_x + side * slot_offset, step = grid_step / 2)
  slot_y <- .dnmb_cct_snap_to_grid(
    target_y + grid_step * (0.30 + 0.22 * ((lane_rank - 1L) %% 3L)),
    step = grid_step / 2
  )
  if (abs(slot_x - target_x) < grid_step * 0.15) {
    slot_x <- .dnmb_cct_snap_to_grid(target_x + side * grid_step * 0.25, step = grid_step / 2)
  }
  if (slot_y > approach_y - grid_step * 0.1) {
    slot_y <- .dnmb_cct_snap_to_grid(target_y + grid_step * 0.25, step = grid_step / 2)
  }

  .dnmb_cct_waypoint_route_points(
    data.frame(
      x = c(hub_x, hub_x, approach_x, approach_x, slot_x, slot_x, target_x),
      y = c(hub_y, corridor_y, corridor_y, approach_y, approach_y, target_y, target_y)
    ),
    grid_step = grid_step
  )
}

.dnmb_cct_waypoint_route_points <- function(points_df, grid_step = 0.5) {
  if (is.null(points_df) || nrow(points_df) < 2) return(points_df)
  out <- points_df[1, , drop = FALSE]
  for (i in 2:nrow(points_df)) {
    seg <- .dnmb_cct_pair_route_points(
      x1 = points_df$x[i - 1L],
      y1 = points_df$y[i - 1L],
      x2 = points_df$x[i],
      y2 = points_df$y[i],
      grid_step = grid_step
    )
    out <- rbind(out, seg[-1, , drop = FALSE])
  }
  keep <- c(TRUE, diff(out$x) != 0 | diff(out$y) != 0)
  out[keep, , drop = FALSE]
}

.dnmb_cct_entry_route_nodes <- function(cs_id, node_ids = NULL) {
  cid <- tolower(as.character(cs_id)[1])
  entry_map <- .dnmb_cct_auto_entry_map()
  inter_map <- .dnmb_cct_auto_intermediates()
  seq_ids <- character(0)
  if (cid %in% names(inter_map)) {
    seq_ids <- vapply(inter_map[[cid]], function(x) as.character(x[1]), character(1))
  }
  final_id <- unname(entry_map[cid])
  seq_ids <- c(seq_ids, final_id)
  seq_ids <- seq_ids[!is.na(seq_ids) & nzchar(seq_ids)]
  seq_ids <- seq_ids[!duplicated(seq_ids)]
  if (!is.null(node_ids)) seq_ids <- seq_ids[seq_ids %in% node_ids]
  seq_ids
}

.dnmb_cct_points_from_start_to_nodes <- function(start_x, start_y, node_ids, node_x, node_y,
                                                 grid_step = 0.5,
                                                 lane_rank = 1L, lane_count = 1L) {
  if (length(node_ids) == 0) return(NULL)
  keep_ids <- node_ids[node_ids %in% names(node_x) & node_ids %in% names(node_y)]
  if (length(keep_ids) == 0) return(NULL)
  # y-offset to separate overlapping routes targeting the same nodes
  y_spread <- if (lane_count > 1L) {
    (lane_rank - (lane_count + 1) / 2) * grid_step * 0.18
  } else 0
  cur_x <- start_x
  cur_y <- start_y
  out <- data.frame(x = cur_x, y = cur_y)
  for (nid in keep_ids) {
    tgt_y <- unname(node_y[nid]) + y_spread
    seg <- .dnmb_cct_pair_route_points(
      x1 = cur_x, y1 = cur_y,
      x2 = unname(node_x[nid]), y2 = tgt_y,
      grid_step = grid_step
    )
    out <- rbind(out, seg[-1, , drop = FALSE])
    cur_x <- unname(node_x[nid])
    cur_y <- tgt_y
  }
  keep <- c(TRUE, diff(out$x) != 0 | diff(out$y) != 0)
  out[keep, , drop = FALSE]
}

.dnmb_cct_route_label_index <- function(points_df, frac = 0.50) {
  if (is.null(points_df) || nrow(points_df) <= 1) return(1L)
  # Use cumulative path length for accurate position along polyline
  dx <- diff(points_df$x)
  dy <- diff(points_df$y)
  seg_len <- sqrt(dx^2 + dy^2)
  cum_len <- c(0, cumsum(seg_len))
  total_len <- cum_len[length(cum_len)]
  if (total_len <= 0) return(1L)
  target_len <- total_len * frac
  idx <- which(cum_len >= target_len)[1]
  as.integer(max(1L, min(nrow(points_df), idx)))
}

.dnmb_cct_edge_points_from_row <- function(ce, edge_idx = 1L, grid_step = 0.5) {
  x1 <- ce$x; y1 <- ce$y; x2 <- ce$xend; y2 <- ce$yend
  if (any(is.na(c(x1, y1, x2, y2)))) return(data.frame(x = numeric(), y = numeric()))
  stagger <- (edge_idx %% 5 - 2) * 0.25
  is_straight <- abs(x2 - x1) < 0.01 || abs(y2 - y1) < 0.01
  if (is_straight) {
    return(data.frame(x = c(x1, x2), y = c(y1, y2)))
  }
  is_mostly_vertical <- abs(y2 - y1) >= abs(x2 - x1)
  if (is_mostly_vertical) {
    mid_y <- round((y2 + stagger) * 4) / 4
    route_pts <- data.frame(
      x = c(x1, x1, x2, x2),
      y = c(y1, mid_y, mid_y, y2)
    )
  } else {
    mid_x <- round((x2 + stagger) * 4) / 4
    route_pts <- data.frame(
      x = c(x1, mid_x, mid_x, x2),
      y = c(y1, y1, y2, y2)
    )
  }
  .dnmb_cct_rounded_route_points(route_pts, radius = grid_step)
}

.dnmb_cct_pathway_hint_nodes <- function(path_ids, matched_steps, transporters, entry_map = NULL) {
  pids <- unique(tolower(as.character(path_ids)))
  pids <- pids[!is.na(pids) & nzchar(pids)]
  if (length(pids) == 0) return(character(0))
  out <- character(0)
  for (pid in pids) {
    out <- union(out, names(.dnmb_cct_pathway_node_bonus(pid, matched_steps = matched_steps)))
    out <- union(out, names(.dnmb_cct_pathway_transport_bonus(pid, transporters = transporters, entry_map = entry_map)))
  }
  out
}

.dnmb_cct_layout_entry_labels <- function(df, center_x = NA_real_) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!is.finite(center_x)) center_x <- stats::median(df$x, na.rm = TRUE)
  out <- df
  out$x_lab <- out$x
  out$y_lab <- out$y + 0.15
  out$hjust <- 0.5

  cluster_key <- paste0(round(out$x / 0.5) * 0.5, "::", round(out$y / 0.5) * 0.5)
  split_idx <- split(seq_len(nrow(out)), cluster_key)
  for (grp in split_idx) {
    if (length(grp) == 1L) {
      ii <- grp[1]
      if (abs(out$x[ii] - center_x) <= 1.0) {
        side <- ifelse(out$x[ii] <= center_x, -1, 1)
        out$x_lab[ii] <- out$x[ii] + side * 0.22
        out$hjust[ii] <- if (side < 0) 1 else 0
      }
      next
    }

    ord <- grp[order(out$priority[grp], decreasing = TRUE, na.last = TRUE)]
    dx_seq <- c(-0.26, 0.26, -0.44, 0.44, -0.62, 0.62)
    dy_seq <- c(0.18, 0.18, 0.32, 0.32, 0.46, 0.46)
    for (k in seq_along(ord)) {
      ii <- ord[k]
      idx <- ((k - 1L) %% length(dx_seq)) + 1L
      shift_x <- dx_seq[idx]
      if (abs(out$x[ii] - center_x) > 1.4) {
        shift_x <- shift_x * 0.75
      }
      out$x_lab[ii] <- out$x[ii] + shift_x
      out$y_lab[ii] <- out$y[ii] + dy_seq[idx]
      out$hjust[ii] <- if (shift_x < -0.05) 1 else if (shift_x > 0.05) 0 else 0.5
    }
  }
  out
}

.dnmb_cct_layout_route_labels <- function(df, center_x = NA_real_,
                                          max_iter = 200L,
                                          label_w = 0.30,
                                          label_h = 0.10,
                                          k_repel = 0.008,
                                          k_attract = 0.06,
                                          damping = 0.85,
                                          temp_start = 0.12,
                                          temp_end = 0.005) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!is.finite(center_x)) center_x <- stats::median(df$x, na.rm = TRUE)
  n <- nrow(df)
  out <- df
  # Initial placement: offset to left or right of anchor depending on side
  side <- ifelse(out$x <= center_x, -1, 1)
  out$x_lab <- out$x + side * 0.22
  out$y_lab <- out$y
  out$hjust <- ifelse(side < 0, 1, 0)

  if (n <= 1L) return(out)

  # Force-directed relaxation loop (simulated annealing)
  anchor_x <- out$x
  anchor_y <- out$y
  vx <- rep(0, n)
  vy <- rep(0, n)

  for (iter in seq_len(max_iter)) {
    temp <- temp_start * ((temp_end / temp_start) ^ ((iter - 1) / max(1, max_iter - 1)))
    fx <- rep(0, n)
    fy <- rep(0, n)

    # Repulsion: push overlapping labels apart
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        dx <- out$x_lab[i] - out$x_lab[j]
        dy <- out$y_lab[i] - out$y_lab[j]
        ox <- label_w - abs(dx)
        oy <- label_h - abs(dy)
        if (ox > 0 && oy > 0) {
          # Overlap detected — push apart proportional to overlap area
          push <- k_repel * ox * oy / max(0.001, sqrt(dx^2 + dy^2))
          if (abs(dx) < 0.001) dx <- (runif(1) - 0.5) * 0.01
          if (abs(dy) < 0.001) dy <- (runif(1) - 0.5) * 0.01
          norm_d <- sqrt(dx^2 + dy^2)
          fx[i] <- fx[i] + push * dx / norm_d
          fy[i] <- fy[i] + push * dy / norm_d
          fx[j] <- fx[j] - push * dx / norm_d
          fy[j] <- fy[j] - push * dy / norm_d
        }
      }
    }
    # Attraction: pull labels toward their anchor point
    for (i in seq_len(n)) {
      dx_a <- anchor_x[i] - out$x_lab[i]
      dy_a <- anchor_y[i] - out$y_lab[i]
      fx[i] <- fx[i] + k_attract * dx_a
      fy[i] <- fy[i] + k_attract * dy_a
    }
    # Update with velocity damping and temperature
    vx <- damping * vx + fx
    vy <- damping * vy + fy
    # Clamp to temperature
    speed <- sqrt(vx^2 + vy^2)
    clamp <- pmin(1, temp / pmax(speed, 1e-6))
    vx <- vx * clamp
    vy <- vy * clamp
    out$x_lab <- out$x_lab + vx
    out$y_lab <- out$y_lab + vy
  }
  # Update hjust based on final label-vs-anchor x
  out$hjust <- ifelse(out$x_lab < out$x - 0.05, 1,
                      ifelse(out$x_lab > out$x + 0.05, 0, 0.5))
  out
}

.dnmb_cct_short_gh_label <- function(gh_family, gene_name = NA_character_) {
  gh <- as.character(gh_family)[1]
  gn <- as.character(gene_name)[1]
  if (is.na(gh) || !nzchar(gh)) gh <- "GH"
  if (is.na(gn) || !nzchar(gn)) return(gh)
  if (.dnmb_cct_is_generic_feature_id(gn)) return(gh)
  gn <- sub("/.*$", "", gn)
  gn <- sub("\\|.*$", "", gn)
  gn <- gsub("_", "-", gn, fixed = TRUE)
  if (nchar(gn) > 10) gn <- substr(gn, 1, 10)
  paste0(gh, " ", gn)
}

.dnmb_cct_short_metabolite_label <- function(label) {
  x <- as.character(label)[1]
  if (is.na(x) || !nzchar(x)) return(x)
  short_map <- c(
    "Glycerol-3-P" = "Gro-3-P",
    "Xylulose" = "Xylu",
    "Ribulose" = "Ribu",
    "Fuculose" = "Fucul"
  )
  out <- unname(short_map[x])
  if (is.na(out) || !nzchar(out)) x else out
}

.dnmb_cct_short_context_summary <- function(context_hits) {
  x <- as.character(context_hits)[1]
  if (is.na(x) || !nzchar(x)) return("")
  parts <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return("")
  pri <- c("operon-prefGH", "prefGH", "operon-dbCAN", "dbCAN", "PAZy", "operon-text", "neighbor-text")
  parts <- c(intersect(pri, parts), setdiff(parts, pri))
  parts <- parts[1:min(2, length(parts))]
  parts <- gsub("^operon-prefGH:", "opGH:", parts)
  parts <- gsub("^prefGH:", "GH:", parts)
  parts <- gsub("^operon-dbCAN:", "op:", parts)
  parts <- gsub("^dbCAN:", "GH:", parts)
  parts <- gsub("^PAZy:", "PZ:", parts)
  parts <- gsub("^operon-text$", "op", parts)
  parts <- gsub("^neighbor-text$", "ctx", parts)
  parts <- gsub("^cgc=", "CGC:", parts)
  parts <- gsub("^cgc-type=", "type:", parts)
  paste(parts, collapse = " / ")
}

.dnmb_cct_estimate_text_width <- function(label) {
  x <- as.character(label)[1]
  if (is.na(x) || !nzchar(x)) return(0.3)
  lines <- unlist(strsplit(x, "\n", fixed = TRUE))
  max(0.3, max(nchar(lines), na.rm = TRUE) * 0.06)
}

.dnmb_cct_lighten_color <- function(color, amount = 0.65) {
  rgb <- grDevices::col2rgb(color) / 255
  mix <- rgb + (1 - rgb) * amount
  grDevices::rgb(mix[1, 1], mix[2, 1], mix[3, 1])
}

.dnmb_cct_gradient_path_layers <- function(points_df, color, linewidth = 0.5,
                                           alpha = 0.7, linetype = "solid",
                                           arrow_last = FALSE,
                                           arrow_length = 0.02,
                                           trim_start = 0,
                                           trim_end = 0) {
  if (is.null(points_df) || nrow(points_df) < 2) return(list())
  points_df <- .dnmb_cct_trim_path_ends(points_df, trim_start = trim_start, trim_end = trim_end)
  if (is.null(points_df) || nrow(points_df) < 2) return(list())
  base_cols <- c(
    .dnmb_cct_lighten_color(color, amount = 0.82),
    .dnmb_cct_lighten_color(color, amount = 0.48),
    color
  )
  base_lw <- linewidth * c(1.55, 1.20, 1.00)
  base_alpha <- c(max(0.06, alpha * 0.18), max(0.10, alpha * 0.40), alpha)

  make_layer <- function(col, lw, al, with_arrow = FALSE) {
    if (nrow(points_df) == 2L) {
      seg_df <- data.frame(
        x = points_df$x[1], y = points_df$y[1],
        xend = points_df$x[2], yend = points_df$y[2]
      )
      return(ggplot2::geom_segment(
        data = seg_df,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        linewidth = lw, color = col, alpha = al,
        linetype = linetype, lineend = "round",
        arrow = if (with_arrow) ggplot2::arrow(length = grid::unit(arrow_length, "inches"), type = "closed") else NULL,
        inherit.aes = FALSE
      ))
    }
    ggplot2::geom_path(
      data = points_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      linewidth = lw, color = col, alpha = al,
      linetype = linetype, lineend = "round", linejoin = "round",
      inherit.aes = FALSE
    )
  }

  layers <- list(
    make_layer(base_cols[1], base_lw[1], base_alpha[1], with_arrow = FALSE),
    make_layer(base_cols[2], base_lw[2], base_alpha[2], with_arrow = FALSE),
    make_layer(base_cols[3], base_lw[3], base_alpha[3], with_arrow = FALSE)
  )
  if (isTRUE(arrow_last)) {
    seg_df <- data.frame(
      x = points_df$x[nrow(points_df) - 1L], y = points_df$y[nrow(points_df) - 1L],
      xend = points_df$x[nrow(points_df)], yend = points_df$y[nrow(points_df)]
    )
    layers[[length(layers) + 1L]] <- ggplot2::geom_segment(
      data = seg_df,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      linewidth = linewidth, color = color, alpha = alpha,
      linetype = linetype, lineend = "round",
      arrow = ggplot2::arrow(length = grid::unit(arrow_length, "inches"), type = "closed"),
      inherit.aes = FALSE
    )
  }
  layers
}

.dnmb_cct_single_path_layer <- function(points_df, color, linewidth = 0.4,
                                        alpha = 0.35, linetype = "solid",
                                        arrow_last = FALSE,
                                        arrow_length = 0.02,
                                        trim_start = 0,
                                        trim_end = 0) {
  if (is.null(points_df) || nrow(points_df) < 2) return(NULL)
  points_df <- .dnmb_cct_trim_path_ends(points_df, trim_start = trim_start, trim_end = trim_end)
  if (is.null(points_df) || nrow(points_df) < 2) return(NULL)
  arr <- if (isTRUE(arrow_last)) ggplot2::arrow(length = grid::unit(arrow_length, "inches"), type = "closed") else NULL
  if (nrow(points_df) == 2) {
    return(ggplot2::geom_segment(
      data = data.frame(x = points_df$x[1], y = points_df$y[1], xend = points_df$x[2], yend = points_df$y[2]),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      linewidth = linewidth, color = color, alpha = alpha,
      linetype = linetype, lineend = "round", arrow = arr,
      inherit.aes = FALSE
    ))
  }
  ggplot2::geom_path(
    data = points_df,
    ggplot2::aes(x = .data$x, y = .data$y),
    linewidth = linewidth, color = color, alpha = alpha,
    linetype = linetype, lineend = "round", linejoin = "round",
    arrow = arr, inherit.aes = FALSE
  )
}

.dnmb_cct_soft_path_points <- function(points_df, n_iter = 2L) {
  if (is.null(points_df) || nrow(points_df) <= 2) return(points_df)
  pts <- as.data.frame(points_df[, c("x", "y"), drop = FALSE])
  for (iter in seq_len(max(1L, as.integer(n_iter)))) {
    if (nrow(pts) <= 2) break
    new_pts <- pts[1, , drop = FALSE]
    for (i in seq_len(nrow(pts) - 1L)) {
      p0 <- pts[i, , drop = FALSE]
      p1 <- pts[i + 1L, , drop = FALSE]
      q <- data.frame(
        x = 0.75 * p0$x + 0.25 * p1$x,
        y = 0.75 * p0$y + 0.25 * p1$y
      )
      r <- data.frame(
        x = 0.25 * p0$x + 0.75 * p1$x,
        y = 0.25 * p0$y + 0.75 * p1$y
      )
      new_pts <- rbind(new_pts, q, r)
    }
    new_pts <- rbind(new_pts, pts[nrow(pts), , drop = FALSE])
    keep <- c(TRUE, diff(new_pts$x) != 0 | diff(new_pts$y) != 0)
    pts <- new_pts[keep, , drop = FALSE]
  }
  pts$x <- .dnmb_cct_snap_to_grid(pts$x, step = 0.25)
  pts$y <- .dnmb_cct_snap_to_grid(pts$y, step = 0.25)
  keep <- c(TRUE, diff(pts$x) != 0 | diff(pts$y) != 0)
  pts <- pts[keep, , drop = FALSE]
  pts
}

.dnmb_cct_trim_path_ends <- function(points_df, trim_start = 0, trim_end = 0) {
  if (is.null(points_df) || nrow(points_df) < 2) return(points_df)
  trim_once <- function(df, amount, from_start = TRUE) {
    if (!is.finite(amount) || amount <= 0 || nrow(df) < 2) return(df)
    idx_seq <- if (from_start) seq_len(nrow(df) - 1L) else rev(seq_len(nrow(df) - 1L))
    remain <- amount
    pts <- df
    for (ii in idx_seq) {
      i0 <- if (from_start) ii else ii + 1L
      i1 <- if (from_start) ii + 1L else ii
      dx <- pts$x[i1] - pts$x[i0]
      dy <- pts$y[i1] - pts$y[i0]
      seg_len <- sqrt(dx^2 + dy^2)
      if (!is.finite(seg_len) || seg_len <= 1e-9) next
      if (remain < seg_len) {
        frac <- remain / seg_len
        if (from_start) {
          pts$x[i0] <- pts$x[i0] + dx * frac
          pts$y[i0] <- pts$y[i0] + dy * frac
        } else {
          # Move endpoint i0 toward i1 by 'frac' of segment length
          pts$x[i0] <- pts$x[i0] + (pts$x[i1] - pts$x[i0]) * frac
          pts$y[i0] <- pts$y[i0] + (pts$y[i1] - pts$y[i0]) * frac
        }
        break
      } else {
        remain <- remain - seg_len
        if (from_start) {
          pts <- pts[-i0, , drop = FALSE]
        } else {
          pts <- pts[-i0, , drop = FALSE]
        }
        if (nrow(pts) < 2) return(pts)
      }
    }
    pts
  }
  out <- trim_once(points_df, trim_start, from_start = TRUE)
  out <- trim_once(out, trim_end, from_start = FALSE)
  keep <- c(TRUE, diff(out$x) != 0 | diff(out$y) != 0)
  out[keep, , drop = FALSE]
}

.dnmb_cct_product_route_weight <- function(product_sugar, primary_product, product_vec = NULL) {
  ps <- tolower(.dnmb_snfg_normalize_sugar_type_v2(product_sugar))[1]
  pp <- tolower(.dnmb_snfg_normalize_sugar_type_v2(primary_product))[1]
  pv <- tolower(.dnmb_snfg_normalize_sugar_type_v2(product_vec))
  if (!is.na(pp) && nzchar(pp) && identical(ps, pp)) return(1.0)
  if (length(pv) > 1 && ps %in% pv) return(0.72)
  0.55
}

.dnmb_cct_path_importance <- function(confidence = "medium", fraction = NA_real_, sink_weight = 1) {
  conf_score <- switch(as.character(confidence),
    high = 1.00, medium = 0.78, low = 0.56, none = 0.34, 0.50
  )
  frac_score <- ifelse(is.na(fraction), 0.50, pmax(0, pmin(1, fraction)))
  sink_score <- pmax(0.35, pmin(1.0, sink_weight))
  0.45 * conf_score + 0.40 * frac_score + 0.15 * sink_score
}

.dnmb_cct_step_target_nodes <- function(step_id, pathway_id = NA_character_) {
  sid <- tolower(as.character(step_id)[1])
  pid <- tolower(as.character(pathway_id)[1])
  if (is.na(sid) || !nzchar(sid)) return(character(0))
  entry_map <- .dnmb_cct_auto_entry_map()
  entry_node <- unname(entry_map[pid])

  is_transport_like <- grepl(
    "pts|transport|permease|porter|abc|binding|mfs|symporter|tm00|ggu|mgl|rbs|xyl[FGH]|lac[EFG]|mal[EFGK]|thu[EFGK]|ara[ETUV]|galP|glc[PTUV]|kguT|nup[ABC]|pot[ABCD]|mtl[EK]|chvE",
    sid, ignore.case = TRUE
  )
  if (is_transport_like && !is.na(entry_node) && nzchar(entry_node)) {
    return(c(entry_node))
  }

  if (sid %in% c("glk", "mgla", "mglb", "mglc")) return(c("Glc-6-P"))
  if (sid %in% c("nagA")) return(c("GlcN-6-P"))
  if (sid %in% c("nagB", "manA")) return(c("Fru-6-P"))
  if (sid %in% c("galK")) return(c("Gal-1-P"))
  if (sid %in% c("galT")) return(c("UDP-Gal", "Glc-1-P", "UDP-Glc"))
  if (sid %in% c("galE")) return(c("UDP-Gal", "UDP-Glc"))
  if (sid %in% c("pgmA")) return(c("Glc-1-P", "Glc-6-P"))
  if (sid %in% c("gntK", "gnd")) return(c("6-PG"))
  if (sid %in% c("kdgK", "eda")) return(c("KDPG", "Pyruvate", "GA3P"))
  if (sid %in% c("xylA")) return(c("Xylulose", "Xu-5-P"))
  if (sid %in% c("xylB")) return(c("Xu-5-P"))
  if (sid %in% c("rbsK")) return(c("R-5-P"))
  if (sid %in% c("araA")) return(c("Ribulose", "Ribulose-5-P", "Xu-5-P"))
  if (sid %in% c("araB")) return(c("Ribulose-5-P", "Xu-5-P"))
  if (sid %in% c("araD", "gguA", "gguB")) return(c("Xu-5-P"))
  if (sid %in% c("glpK", "glpO")) return(c("Glycerol-3-P", "DHAP"))
  if (sid %in% c("fucI")) return(c("Fuculose", "Fuculose-1-P", "DHAP"))
  if (sid %in% c("fucK")) return(c("Fuculose-1-P", "DHAP"))
  if (sid %in% c("fba")) return(c("Fru-1,6-BP", "GA3P"))
  if (sid %in% c("tpi")) return(c("DHAP", "GA3P"))
  if (sid %in% c("thuK", "pstp")) return(c("Trehalose-6-P", "Glc-6-P"))
  if (sid %in% c("lacZ", "lacE", "lacF", "lacG", "lacK")) return(c("Gal-1-P", "Glc-6-P"))
  if (sid %in% c("ggua", "ggub", "chve", "xylg", "xylh", "rbsa", "rbsb", "rbsc", "nupa", "nupb", "nupc'",
                 "male1", "male2", "malf1", "malk1", "malt1", "mtlk", "mtle", "glcv", "glct", "glu",
                 "aglf", "aglk", "tm0027", "tm0028", "tm0029", "tm0030", "tm0031")) {
    if (!is.na(entry_node) && nzchar(entry_node)) return(c(entry_node))
  }
  if (pid %in% c("gluconate")) return(c("6-PG", "KDPG", "Pyruvate"))
  if (pid %in% c("xylose")) return(c("Xylulose", "Xu-5-P", "GA3P"))
  if (pid %in% c("arabinose")) return(c("Ribulose", "Ribulose-5-P", "Xu-5-P", "GA3P"))
  if (pid %in% c("ribose")) return(c("R-5-P", "GA3P"))
  if (pid %in% c("glycerol")) return(c("Glycerol-3-P", "DHAP"))
  if (pid %in% c("fucose")) return(c("Fuculose", "Fuculose-1-P", "DHAP", "GA3P"))
  character(0)
}

.dnmb_cct_sink_priority <- function(pathway_id, sink_id, entry_target = NULL) {
  pid <- tolower(as.character(pathway_id)[1])
  sid <- as.character(sink_id)[1]
  et <- as.character(entry_target)[1]
  if (sid == "Pyruvate") return(1.00)
  if (sid == "Acetyl-CoA") return(if (pid %in% c("maltose","cellobiose","trehalose","galactose","lactose","mannose","fructose","sucrose","mannitol")) 0.88 else 0.72)
  if (sid == "GA3P") {
    if (!is.na(et) && et %in% c("R-5-P", "Xu-5-P", "6-PG", "6-PGL", "KDPG")) return(0.92)
    if (pid %in% c("xylose", "arabinose", "ribose", "gluconate", "glycerol")) return(0.86)
    return(0.52)
  }
  0.5
}

.dnmb_cct_pathway_node_bonus <- function(pathway_id, matched_steps) {
  pid <- tolower(as.character(pathway_id)[1])
  if (nrow(matched_steps) == 0) return(stats::setNames(numeric(0), character(0)))
  ms <- matched_steps[tolower(matched_steps$pathway_id) == pid, , drop = FALSE]
  if (nrow(ms) == 0) return(stats::setNames(numeric(0), character(0)))
  bonus_map <- list()
  for (i in seq_len(nrow(ms))) {
    tgt_nodes <- .dnmb_cct_step_target_nodes(ms$step_id[i], ms$pathway_id[i])
    if (length(tgt_nodes) == 0) next
    imp <- .dnmb_cct_path_importance(ms$confidence[i], sink_weight = 1)
    node_bonus <- 0.06 + 0.14 * imp
    for (nn in tgt_nodes) {
      bonus_map[[nn]] <- max(c(bonus_map[[nn]], node_bonus), na.rm = TRUE)
    }
  }
  stats::setNames(as.numeric(unlist(bonus_map)), names(bonus_map))
}

.dnmb_cct_pathway_transport_bonus <- function(pathway_id, transporters, entry_map = NULL) {
  pid <- tolower(as.character(pathway_id)[1])
  if (is.null(transporters) || !is.data.frame(transporters) || nrow(transporters) == 0) {
    return(stats::setNames(numeric(0), character(0)))
  }
  tr <- transporters[tolower(transporters$pathway) == pid, , drop = FALSE]
  if (nrow(tr) == 0) return(stats::setNames(numeric(0), character(0)))
  if (!"context_score" %in% names(tr)) tr$context_score <- 0
  conf_rank <- c(none = 0L, low = 1L, medium = 2L, high = 3L)
  tr$rank <- conf_rank[tr$confidence]
  tr$rank[is.na(tr$rank)] <- 0L
  tr$imp <- mapply(
    .dnmb_cct_path_importance,
    confidence = tr$confidence,
    fraction = NA_real_,
    sink_weight = 1
  ) + pmin(2, pmax(0, tr$context_score)) * 0.08

  entry_node <- if (!is.null(entry_map)) unname(entry_map[pid]) else NA_character_
  bonus <- list()
  if (!is.na(entry_node) && nzchar(entry_node)) {
    bonus[[entry_node]] <- max(0.08 + 0.14 * max(tr$imp, na.rm = TRUE), na.rm = TRUE)
  }

  for (i in seq_len(nrow(tr))) {
    tgt_nodes <- .dnmb_cct_step_target_nodes(tr$step[i], tr$pathway[i])
    if (length(tgt_nodes) == 0) next
    node_bonus <- 0.04 + 0.08 * tr$imp[i]
    for (nn in tgt_nodes) {
      bonus[[nn]] <- max(c(bonus[[nn]], node_bonus), na.rm = TRUE)
    }
  }
  stats::setNames(as.numeric(unlist(bonus)), names(bonus))
}

.dnmb_cct_continuity_routes <- function() {
  list(
    maltose = c("Glucose", "Glc-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    cellobiose = c("Glucose", "Glc-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    trehalose = c("Trehalose-6-P", "Glc-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    galactose = c("Gal-1-P", "UDP-Gal", "UDP-Glc", "Glc-1-P", "Glc-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    lactose = c("Gal-1-P", "UDP-Gal", "UDP-Glc", "Glc-1-P", "Glc-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    nag = c("GlcNAc-6-P", "GlcN-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    glucosamine = c("GlcN-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    mannose = c("Man-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    mannitol = c("Mannitol-1-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    fructose = c("Fru-1-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    sucrose = c("Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    glycerol = c("Glycerol-3-P", "DHAP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    xylose = c("Xylulose", "Xu-5-P", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    arabinose = c("Ribulose", "Ribulose-5-P", "Xu-5-P", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    ribose = c("R-5-P", "S-7-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    gluconate = c("6-PG", "KDPG", "Pyruvate"),
    fucose = c("Fuculose", "Fuculose-1-P", "DHAP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    rhamnose = c("Rhamnulose", "Rhamnulose-1-P", "DHAP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate")
  )
}

.dnmb_cct_route_group_for_pathway <- function(pathway_id, route_group_members) {
  pid <- tolower(as.character(pathway_id)[1])
  for (nm in names(route_group_members)) {
    members <- tolower(route_group_members[[nm]])
    if (pid %in% members || identical(pid, tolower(nm))) return(nm)
  }
  pid
}

.dnmb_cct_continuity_preferences <- function(pathway_id, route_group_members, entry_target = NULL) {
  grp <- .dnmb_cct_route_group_for_pathway(pathway_id, route_group_members)
  prefs <- c(grp)
  if (!is.null(entry_target) && !is.na(entry_target) && nzchar(entry_target)) {
    if (entry_target %in% c("R-5-P", "Xu-5-P", "6-PG", "6-PGL", "KDPG")) {
      prefs <- c(prefs, "ppp", "ed")
    }
  }
  if (identical(grp, "gluconate")) prefs <- c(prefs, "ppp", "ed")
  unique(c(prefs, "backbone"))
}

.dnmb_cct_continuity_graph <- function(cyto_edges) {
  edges <- cyto_edges[, c("from", "to", "pathway"), drop = FALSE]
  edges$reverse <- FALSE
  rev_edges <- edges[edges$pathway %in% c("ppp"), , drop = FALSE]
  if (nrow(rev_edges) > 0) {
    rev_tmp <- rev_edges$from
    rev_edges$from <- rev_edges$to
    rev_edges$to <- rev_tmp
    rev_edges$reverse <- TRUE
    edges <- rbind(edges, rev_edges)
  }
  edges
}

.dnmb_cct_candidate_sink_ids <- function(pathway_id, entry_target = NULL) {
  pid <- tolower(as.character(pathway_id)[1])
  et <- as.character(entry_target)[1]
  if (!is.na(et) && et %in% c("6-PG", "6-PGL", "Ru-5-P", "R-5-P", "Xu-5-P", "KDPG")) {
    return(c("Pyruvate", "GA3P", "Acetyl-CoA"))
  }
  if (!is.na(et) && et %in% c("Glycerol-3-P", "DHAP", "GA3P")) {
    return(c("Pyruvate", "Acetyl-CoA"))
  }
  if (pid %in% c("gluconate", "xylose", "arabinose", "ribose")) {
    return(c("Pyruvate", "GA3P", "Acetyl-CoA"))
  }
  c("Pyruvate", "Acetyl-CoA")
}

.dnmb_cct_preferred_nodes <- function(pathway_id, entry_target = NULL) {
  pid <- tolower(as.character(pathway_id)[1])
  et <- as.character(entry_target)[1]
  if (!is.na(et) && et %in% c("6-PG", "6-PGL", "Ru-5-P", "R-5-P", "Xu-5-P", "KDPG")) {
    return(c("6-PGL", "6-PG", "Ru-5-P", "R-5-P", "Xu-5-P", "KDPG", "GA3P", "Pyruvate"))
  }
  if (pid %in% c("gluconate")) {
    return(c("6-PG", "KDPG", "Pyruvate"))
  }
  if (pid %in% c("xylose", "arabinose", "ribose")) {
    return(c("Xu-5-P", "R-5-P", "S-7-P", "E-4-P", "Fru-6-P", "GA3P", "Pyruvate"))
  }
  c("Glc-6-P", "Fru-6-P", "Fru-1,6-BP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate", "Acetyl-CoA")
}

.dnmb_cct_continuity_edge_cost <- function(edge_pathway, to_node, preferred_paths,
                                           preferred_nodes, node_bonus = NULL,
                                           reverse = FALSE) {
  pw <- as.character(edge_pathway)[1]
  cost <- if (pw %in% preferred_paths) 0.16 else if (pw %in% c("backbone", "ppp", "ed")) 0.34 else 0.95
  if (pw == "pyruvate_branch") cost <- 1.75
  if (!is.na(to_node) && to_node %in% preferred_nodes) cost <- cost - 0.10
  if (!is.null(node_bonus) && !is.na(to_node) && to_node %in% names(node_bonus)) {
    cost <- cost - unname(node_bonus[to_node])
  }
  if (isTRUE(reverse)) cost <- cost + 0.16
  pmax(0.05, cost)
}

.dnmb_cct_best_continuity_nodes <- function(pathway_id, cyto_nodes, cyto_edges,
                                            route_group_members, entry_map,
                                            matched_steps = NULL,
                                            transporters = NULL) {
  pid <- tolower(as.character(pathway_id)[1])
  graph <- .dnmb_cct_continuity_graph(cyto_edges)
  node_type <- stats::setNames(cyto_nodes$type, cyto_nodes$id)
  node_y <- stats::setNames(cyto_nodes$y, cyto_nodes$id)
  entry_target <- unname(entry_map[pid])
  preferred_paths <- .dnmb_cct_continuity_preferences(pid, route_group_members, entry_target = entry_target)
  preferred_nodes <- .dnmb_cct_preferred_nodes(pid, entry_target = entry_target)
  sink_ids <- .dnmb_cct_candidate_sink_ids(pid, entry_target = entry_target)
  node_bonus <- .dnmb_cct_pathway_node_bonus(pid, matched_steps = matched_steps)
  tr_bonus <- .dnmb_cct_pathway_transport_bonus(pid, transporters = transporters, entry_map = entry_map)
  if (length(tr_bonus) > 0) {
    all_nodes <- union(names(node_bonus), names(tr_bonus))
    merged_bonus <- setNames(rep(0, length(all_nodes)), all_nodes)
    merged_bonus[names(node_bonus)] <- node_bonus
    merged_bonus[names(tr_bonus)] <- pmax(merged_bonus[names(tr_bonus)], tr_bonus)
    node_bonus <- merged_bonus
  }

  candidate_starts <- unique(c(
    entry_target,
    graph$from[graph$pathway %in% preferred_paths]
  ))
  candidate_starts <- candidate_starts[
    candidate_starts %in% names(node_type) &
      node_type[candidate_starts] %in% c("entry_intermediate", "backbone", "ppp", "ed")
  ]
  sink_ids <- sink_ids[sink_ids %in% names(node_type)]
  if (length(candidate_starts) == 0 || length(sink_ids) == 0) return(NULL)

  best_cost <- Inf
  best_path <- NULL
  node_ids <- unique(c(graph$from, graph$to))
  for (st in candidate_starts) {
    dist <- stats::setNames(rep(Inf, length(node_ids)), node_ids)
    prev <- stats::setNames(rep(NA_character_, length(node_ids)), node_ids)
    visited <- stats::setNames(rep(FALSE, length(node_ids)), node_ids)
    dist[st] <- 0

    repeat {
      open_nodes <- names(dist)[!visited]
      if (length(open_nodes) == 0) break
      cur <- open_nodes[which.min(dist[open_nodes])]
      if (!is.finite(dist[cur])) break
      visited[cur] <- TRUE
      if (cur %in% sink_ids) break
      eidx <- which(graph$from == cur)
      if (length(eidx) == 0) next
      for (ii in eidx) {
        ee <- graph[ii, , drop = FALSE]
        nxt <- ee$to[1]
        step_cost <- .dnmb_cct_continuity_edge_cost(
          edge_pathway = ee$pathway[1],
          to_node = nxt,
          preferred_paths = preferred_paths,
          preferred_nodes = preferred_nodes,
          node_bonus = node_bonus,
          reverse = isTRUE(ee$reverse[1])
        )
        new_cost <- dist[cur] + step_cost
        if (new_cost + 1e-9 < dist[nxt]) {
          dist[nxt] <- new_cost
          prev[nxt] <- cur
        }
      }
    }

    sink_costs <- dist[sink_ids]
    sink_costs <- sink_costs[is.finite(sink_costs)]
    if (length(sink_costs) == 0) next
    sink_adj <- vapply(names(sink_costs), function(sk) {
      sink_costs[sk] - 0.22 * .dnmb_cct_sink_priority(pid, sk, entry_target = entry_target)
    }, numeric(1))
    best_sink <- names(sink_adj)[which.min(sink_adj)]
    path <- best_sink
    while (!is.na(prev[path[1]])) {
      path <- c(prev[path[1]], path)
      if (length(path) > length(node_ids) + 2L) break
    }
    if (length(path) < 2) next
    total_cost <- unname(sink_adj[best_sink]) - 0.03 * unname(node_y[st])
    if (total_cost < best_cost) {
      best_cost <- total_cost
      best_path <- path
      attr(best_path, "sink_id") <- best_sink
      attr(best_path, "sink_weight") <- .dnmb_cct_sink_priority(pid, best_sink, entry_target = entry_target)
    }
  }

  best_path
}

.dnmb_cct_pair_route_points <- function(x1, y1, x2, y2, grid_step = 0.5) {
  if (any(is.na(c(x1, y1, x2, y2)))) return(data.frame(x = numeric(), y = numeric()))
  if (abs(x2 - x1) < 0.01 || abs(y2 - y1) < 0.01) {
    return(data.frame(x = c(x1, x2), y = c(y1, y2)))
  }
  # Always vertical-first: go down (or up) first, then horizontal
  pts <- data.frame(x = c(x1, x1, x2), y = c(y1, y2, y2))
  # Radius limited to half the shorter segment to guarantee connectivity
  seg_v <- abs(y2 - y1)
  seg_h <- abs(x2 - x1)
  safe_radius <- min(grid_step * 0.4, seg_v * 0.45, seg_h * 0.45)
  .dnmb_cct_rounded_route_points(pts, radius = max(0.05, safe_radius))
}

.dnmb_cct_route_overlay_points <- function(node_ids, node_x, node_y, grid_step = 0.5) {
  keep_ids <- node_ids[node_ids %in% names(node_x) & node_ids %in% names(node_y)]
  if (length(keep_ids) < 2) return(NULL)
  out <- NULL
  for (i in 2:length(keep_ids)) {
    seg <- .dnmb_cct_pair_route_points(
      x1 = unname(node_x[keep_ids[i - 1]]),
      y1 = unname(node_y[keep_ids[i - 1]]),
      x2 = unname(node_x[keep_ids[i]]),
      y2 = unname(node_y[keep_ids[i]]),
      grid_step = grid_step
    )
    if (is.null(out)) {
      out <- seg
    } else {
      out <- rbind(out, seg[-1, , drop = FALSE])
    }
  }
  keep <- c(TRUE, diff(out$x) != 0 | diff(out$y) != 0)
  out[keep, , drop = FALSE]
}

.dnmb_cct_extra_preferred_lane <- function(substrate_id, extra_chains = NULL) {
  overrides <- c(
    Starch = "glucose",
    Cellulose = "glucose",
    Maltose = "glucose",
    Cellobiose = "glucose",
    Trehalose = "glucose",
    Chitin = "glcnac",
    Xylan = "xylose",
    Pectin = "gala",
    Mannan = "mannose",
    Inulin = "fructose",
    Sucrose = "fructose",
    Lactose = "galactose",
    Raffinose = "galactose"
  )
  if (substrate_id %in% names(overrides)) return(unname(overrides[substrate_id]))

  chain <- extra_chains[[substrate_id]]
  if (!is.null(chain) && length(chain) > 0) return(chain[1])
  "glucose"
}

.dnmb_cct_extra_sugar_profile <- function(substrate_id, extra_chains, extra_branches = NULL,
                                          primary_boost = 0, branch_weight = 0.65) {
  primary_sugar <- tolower(.dnmb_snfg_normalize_sugar_type_v2(
    .dnmb_cct_extra_preferred_lane(substrate_id, extra_chains = extra_chains)
  ))
  chain <- extra_chains[[substrate_id]]
  if (is.null(chain)) chain <- character(0)
  chain_norm <- tolower(.dnmb_snfg_normalize_sugar_type_v2(chain))

  branches <- extra_branches[[substrate_id]]
  branch_norm <- if (is.null(branches) || length(branches) == 0) {
    character(0)
  } else {
    tolower(.dnmb_snfg_normalize_sugar_type_v2(
      vapply(branches, function(br) br$sugar, character(1))
    ))
  }

  sugar_vec <- c(chain_norm, branch_norm)
  if (length(sugar_vec) == 0) return(stats::setNames(numeric(0), character(0)))

  weight_vec <- c(rep(1, length(chain_norm)), rep(branch_weight, length(branch_norm)))
  if (!is.na(primary_sugar) && nzchar(primary_sugar) && primary_sugar %in% sugar_vec) {
    weight_vec[sugar_vec == primary_sugar] <- weight_vec[sugar_vec == primary_sugar] + primary_boost
  }

  prof <- tapply(weight_vec, sugar_vec, sum)
  prof[order(names(prof))]
}

.dnmb_cct_extra_anchor_x <- function(substrate_id, extra_chains, extra_branches,
                                     sugar_lane_x, default_anchor = NULL) {
  chain <- extra_chains[[substrate_id]]
  if (is.null(chain)) chain <- character(0)
  primary_boost <- if (length(chain) >= 4) 2 else if (length(unique(chain)) > 1) 1 else 0
  prof <- .dnmb_cct_extra_sugar_profile(
    substrate_id = substrate_id,
    extra_chains = extra_chains,
    extra_branches = extra_branches,
    primary_boost = primary_boost,
    branch_weight = 0.55
  )
  prof <- prof[names(prof) %in% names(sugar_lane_x)]
  if (length(prof) == 0) {
    if (is.null(default_anchor)) default_anchor <- stats::median(unname(sugar_lane_x))
    return(default_anchor)
  }
  stats::weighted.mean(unname(sugar_lane_x[names(prof)]), w = as.numeric(prof))
}

.dnmb_cct_extra_row_key <- function(n_units) {
  if (is.na(n_units) || n_units <= 2) return("2")
  if (n_units == 3) return("3")
  if (n_units == 4) return("4")
  if (n_units == 5) return("5")
  "6"  # 6+ units (Pullulan) get their own row
}

.dnmb_cct_extra_target_stats <- function(cascade_key, substrate_id, cascades,
                                         sugar_lane_x, default_target = NULL,
                                         extra_chains = NULL, extra_branches = NULL) {
  if (!is.null(cascade_key) && !is.na(cascade_key) && cascade_key %in% names(cascades)) {
    prod <- cascades[[cascade_key]]$products
    if (is.null(prod)) prod <- cascades[[cascade_key]]$product
    prod <- tolower(.dnmb_snfg_normalize_sugar_type_v2(prod))
    prod <- prod[prod %in% names(sugar_lane_x)]
    if (length(prod) > 0) {
      prod_w <- stats::setNames(rep(1, length(prod)), prod)
      prof <- .dnmb_cct_extra_sugar_profile(
        substrate_id = substrate_id,
        extra_chains = extra_chains,
        extra_branches = extra_branches,
        primary_boost = 1,
        branch_weight = 0.45
      )
      if (length(prof) > 0) {
        shared <- intersect(names(prod_w), names(prof))
        if (length(shared) > 0) prod_w[shared] <- prod_w[shared] + prof[shared]
      }
      primary_prod <- cascades[[cascade_key]]$product
      primary_prod <- tolower(.dnmb_snfg_normalize_sugar_type_v2(primary_prod))
      if (length(primary_prod) > 0 && primary_prod[1] %in% names(prod_w)) {
        prod_w[primary_prod[1]] <- prod_w[primary_prod[1]] + 1.5
      }
      xs <- unname(sugar_lane_x[names(prod_w)])
      return(c(
        center = stats::weighted.mean(xs, w = as.numeric(prod_w)),
        left = min(xs),
        right = max(xs)
      ))
    }
  }
  if (is.null(default_target)) default_target <- stats::median(unname(sugar_lane_x))
  c(center = default_target, left = default_target, right = default_target)
}

.dnmb_cct_cascade_cut_x <- function(cascade_key, chain_x, n_mono, step,
                                    matched_families = character(0)) {
  fams <- toupper(as.character(matched_families))
  if (!length(fams)) fams <- character(0)
  # Exo-acting GH families cleave from the non-reducing end (rightmost bond).
  # Endo-acting GH families cleave internal bonds (penultimate bond).
  # If both present, prefer exo (non-reducing end) since that's the released monomer.
  exo_ghs <- list(
    starch    = "^GH13(_20|_31|_39|_45)?$|^GH4$|^GH65$|^GH15$",
    cellulose = "^GH1$|^GH3$|^GH4$",
    xylan     = "^GH43$|^GH3$",
    chitin    = "^GH20$|^GH3$",
    mannan    = "^GH92$|^GH130|^GH5(_7)?$",
    pectin    = "^GH28$|^GH78$",
    betaglucan = "^GH3$|^GH55$",
    arabinoxylan = "^GH43$|^GH3$",
    agarose   = "^GH50$|^GH117$"
  )
  ckey <- tolower(as.character(cascade_key)[1])
  # Di/tri-saccharides: always first bond (only one bond to cut)
  if (ckey %in% c("raffinose","lactose","sucrose","trehalose","maltose","cellobiose")) {
    return(chain_x + step * 0.5)
  }
  # Polymers: check if exo-acting GH is present → cut at non-reducing end (last bond)
  has_exo <- if (ckey %in% names(exo_ghs)) any(grepl(exo_ghs[[ckey]], fams)) else FALSE
  bond_idx <- if (has_exo) {
    max(1, n_mono - 1L)   # last bond (non-reducing end)
  } else {
    max(1, n_mono - 2L)   # penultimate bond (internal endo cleavage)
  }
  chain_x + (bond_idx - 1L) * step + step * 0.5
}

.dnmb_cct_extra_product_exit_x <- function(chain, chain_x, product_sugar,
                                           branches = NULL, sym_size = 0.075, gap = 0.02,
                                           target_x = NA_real_) {
  if (is.null(chain)) chain <- character(0)
  step <- sym_size * 2 + gap
  chain_norm <- tolower(.dnmb_snfg_normalize_sugar_type_v2(chain))
  prod_norm <- tolower(.dnmb_snfg_normalize_sugar_type_v2(product_sugar))
  if (length(prod_norm) == 0) prod_norm <- ""
  prod_norm <- prod_norm[1]
  composite_products <- list(
    sucrose = c("glucose", "fructose"),
    lactose = c("galactose", "glucose"),
    trehalose = c("glucose", "glucose"),
    maltose = c("glucose", "glucose"),
    cellobiose = c("glucose", "glucose"),
    raffinose = c("galactose", "glucose", "fructose")
  )
  composite_anchor <- c(
    sucrose = "fructose",
    lactose = "galactose",
    trehalose = "glucose",
    maltose = "glucose",
    cellobiose = "glucose",
    raffinose = "fructose"
  )

  if (length(chain_norm) > 0) {
    chain_xs <- chain_x + (seq_along(chain_norm) - 1) * step
    if (prod_norm %in% names(composite_products)) {
      target_seq <- composite_products[[prod_norm]]
      n_target <- length(target_seq)
      if (length(chain_norm) >= n_target) {
        centers <- numeric(0)
        anchor_hits <- numeric(0)
        for (i in seq_len(length(chain_norm) - n_target + 1L)) {
          if (identical(chain_norm[i:(i + n_target - 1L)], target_seq)) {
            centers <- c(centers, mean(chain_xs[i:(i + n_target - 1L)]))
            anchor_name <- unname(composite_anchor[prod_norm])
            if (!is.na(anchor_name) && nzchar(anchor_name)) {
              local_idx <- which(target_seq == anchor_name)
              if (length(local_idx) > 0) {
                anchor_hits <- c(anchor_hits, chain_xs[i + local_idx[length(local_idx)] - 1L])
              } else {
                anchor_hits <- c(anchor_hits, mean(chain_xs[i:(i + n_target - 1L)]))
              }
            } else {
              anchor_hits <- c(anchor_hits, mean(chain_xs[i:(i + n_target - 1L)]))
            }
          }
        }
        if (length(anchor_hits) > 0) {
          if (is.finite(target_x)) return(anchor_hits[which.min(abs(anchor_hits - target_x))])
          return(anchor_hits[1])
        }
      }
    }
    hit_x <- chain_xs[chain_norm == prod_norm]
    if (length(hit_x) > 0) {
      if (is.finite(target_x)) return(hit_x[which.min(abs(hit_x - target_x))])
      return(hit_x[1])
    }
  }

  if (!is.null(branches) && length(branches) > 0) {
    br_sugars <- tolower(.dnmb_snfg_normalize_sugar_type_v2(
      vapply(branches, function(br) br$sugar, character(1))
    ))
    br_pos <- vapply(branches, function(br) br$pos, numeric(1))
    br_x <- chain_x + (br_pos - 1) * step
    hit_x <- br_x[br_sugars == prod_norm]
    if (length(hit_x) > 0) {
      if (is.finite(target_x)) return(hit_x[which.min(abs(hit_x - target_x))])
      return(hit_x[1])
    }
  }

  chain_x + max(0, length(chain_norm) - 1) * step / 2
}

.dnmb_cct_small_permutations <- function(x, max_full = 7L) {
  x <- as.integer(x)
  n <- length(x)
  if (n <= 1) return(list(x))
  # Hard cap: for n > max_full the full n! enumeration explodes
  # (factorial(8)=40320, factorial(10)=3.6M, factorial(12)=479M list
  # pointers). The row-layout optimiser that consumes this list only
  # needs a handful of useful orderings, not every permutation, so fall
  # back to a curated shortlist when the row grows beyond max_full.
  if (n > max_full) {
    return(list(
      x,                                      # identity
      rev(x),                                 # reverse
      x[order(x)],                            # ascending
      x[order(-x)]                            # descending
    ))
  }
  out <- vector("list", factorial(n))
  k <- 1L
  for (i in seq_along(x)) {
    tail_perms <- .dnmb_cct_small_permutations(x[-i], max_full = max_full)
    for (tp in tail_perms) {
      out[[k]] <- c(x[i], tp)
      k <- k + 1L
    }
  }
  out
}

.dnmb_cct_extra_pack_row <- function(row_df, min_gap = 0.35, left_bound = 0.5) {
  centers0 <- numeric(nrow(row_df))
  split_idx <- split(seq_len(nrow(row_df)), sprintf("%.3f", row_df$anchor_x))
  for (gidx in split_idx) {
    widths <- row_df$width_eff[gidx]
    anchor <- row_df$anchor_x[gidx[1]]
    span <- sum(widths) + min_gap * max(0, length(gidx) - 1)
    left <- anchor - span / 2
    cur_left <- left
    for (j in seq_along(gidx)) {
      ii <- gidx[j]
      centers0[ii] <- cur_left + widths[j] / 2
      cur_left <- centers0[ii] + widths[j] / 2 + min_gap
    }
  }

  centers <- centers0
  if (nrow(row_df) >= 2) {
    for (i in 2:nrow(row_df)) {
      sep <- (row_df$width_eff[i - 1] + row_df$width_eff[i]) / 2 + min_gap
      centers[i] <- max(centers[i], centers[i - 1] + sep)
    }
  }

  desired_shift <- mean(row_df$anchor_x - centers)
  min_left_edge <- min(centers - row_df$width_eff / 2)
  shift <- max(desired_shift, left_bound - min_left_edge)
  row_df$x_center <- centers + shift
  row_df
}

.dnmb_cct_extra_row_gap <- function(row_key, min_gap = 0.35) {
  rk <- as.character(row_key)[1]
  extra <- switch(rk, `6` = 0.40, `5` = 0.30, `4` = 0.30, 0)
  min_gap + extra
}

.dnmb_cct_extra_chain_step <- function(sym_size = 0.075, gap = 0.18) {
  sym_size * 2 + gap
}

.dnmb_cct_extra_chain_node_x <- function(chain_x, pos, sym_size = 0.075, gap = 0.18) {
  chain_x + (as.numeric(pos)[1] - 1) * .dnmb_cct_extra_chain_step(sym_size = sym_size, gap = gap)
}

.dnmb_cct_extra_route_trim <- function(sym_size = 0.075) {
  max(0.06, sym_size * 1.05)
}

.dnmb_cct_extra_stagger_row <- function(row_df, y_step = 0.26, min_sep = 0.18) {
  if (is.null(row_df) || nrow(row_df) <= 1) return(row_df)
  ord <- order(row_df$x_center, row_df$anchor_x, row_df$id)
  out <- row_df[ord, , drop = FALSE]
  widths <- pmax(out$width_eff, out$chain_width, out$label_width)
  right_edge_by_level <- numeric(0)
  level_assign <- integer(nrow(out))
  for (i in seq_len(nrow(out))) {
    left_edge <- out$x_center[i] - widths[i] / 2
    right_edge <- out$x_center[i] + widths[i] / 2
    placed <- FALSE
    if (length(right_edge_by_level) > 0) {
      for (lev in seq_along(right_edge_by_level)) {
        if (left_edge >= (right_edge_by_level[lev] + min_sep)) {
          level_assign[i] <- lev
          right_edge_by_level[lev] <- right_edge
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      right_edge_by_level <- c(right_edge_by_level, right_edge)
      level_assign[i] <- length(right_edge_by_level)
    }
  }
  offsets <- c(0, y_step, -y_step, 2 * y_step, -2 * y_step, 3 * y_step, -3 * y_step)
  if (max(level_assign, na.rm = TRUE) > length(offsets)) {
    extra_levels <- max(level_assign, na.rm = TRUE) - length(offsets)
    extra_offsets <- y_step * seq_len(extra_levels)
    offsets <- c(offsets, extra_offsets)
  }
  out$y <- out$y + offsets[level_assign]
  out[match(row_df$id, out$id), , drop = FALSE]
}

.dnmb_cct_extra_row_cost <- function(row_df) {
  if (nrow(row_df) == 0) return(0)
  disp_cost <- sum(((row_df$x_center - row_df$anchor_x)^2) / pmax(row_df$width_eff, 0.5))
  interval_dist <- pmax(0, row_df$route_left_x - row_df$x_center, row_df$x_center - row_df$route_right_x)
  lane_cost <- 0.10 * sum(abs(row_df$x_center - row_df$route_target_x)) +
    0.25 * sum(interval_dist)
  cross_cost <- 0
  prox_cost <- 0
  if (nrow(row_df) >= 2) {
    for (i in 1:(nrow(row_df) - 1)) {
      for (j in (i + 1):nrow(row_df)) {
        desired_sep <- (row_df$width_eff[i] + row_df$width_eff[j]) / 2 +
          min(0.60, 0.10 * (row_df$route_span[i] + row_df$route_span[j]))
        actual_sep <- abs(row_df$x_center[i] - row_df$x_center[j])
        if (actual_sep < desired_sep) {
          prox_cost <- prox_cost + (desired_sep - actual_sep)^2 * 3.5
        }
        inv <- (row_df$x_center[i] - row_df$x_center[j]) * (row_df$route_target_x[i] - row_df$route_target_x[j])
        if (inv < 0) {
          overlap <- max(0, min(row_df$route_right_x[i], row_df$route_right_x[j]) -
                            max(row_df$route_left_x[i], row_df$route_left_x[j]))
          span <- max(row_df$route_right_x[i], row_df$route_right_x[j]) -
                  min(row_df$route_left_x[i], row_df$route_left_x[j])
          overlap_frac <- if (span > 0) overlap / span else 0
          cross_cost <- cross_cost + (10 + abs(row_df$route_target_x[i] - row_df$route_target_x[j])) *
            (1 - 0.5 * overlap_frac)
        }
      }
    }
  }
  disp_cost + lane_cost + cross_cost + prox_cost
}

.dnmb_cct_extra_global_cost <- function(extra_df) {
  if (is.null(extra_df) || nrow(extra_df) == 0) return(0)

  cost <- 0
  for (rk in unique(extra_df$row_key)) {
    cost <- cost + .dnmb_cct_extra_row_cost(extra_df[extra_df$row_key == rk, , drop = FALSE])
  }

  if (nrow(extra_df) >= 2) {
    for (i in 1:(nrow(extra_df) - 1)) {
      for (j in (i + 1):nrow(extra_df)) {
        inv <- (extra_df$x_center[i] - extra_df$x_center[j]) *
          (extra_df$route_target_x[i] - extra_df$route_target_x[j])
        if (inv < 0) {
          overlap <- max(0, min(extra_df$route_right_x[i], extra_df$route_right_x[j]) -
                            max(extra_df$route_left_x[i], extra_df$route_left_x[j]))
          span <- max(extra_df$route_right_x[i], extra_df$route_right_x[j]) -
                  min(extra_df$route_left_x[i], extra_df$route_left_x[j])
          overlap_frac <- if (span > 0) overlap / span else 0
          row_dist <- abs(extra_df$y[i] - extra_df$y[j])
          weight <- 14 / (1 + row_dist)
          cost <- cost + (weight + 0.35 * abs(extra_df$route_target_x[i] - extra_df$route_target_x[j])) *
            (1 - 0.5 * overlap_frac)
        }
      }
    }
  }

  cost
}

.dnmb_cct_extra_layout <- function(extra_substrates, extra_chains, sugar_lane_x,
                                   row_y_map, extra_branches = NULL, cascades = NULL,
                                   sym_size = 0.075,
                                   gap = 0.02, min_gap = 0.35) {
  step <- sym_size * 2 + gap
  extra_df <- extra_substrates
  extra_df$n_units <- vapply(extra_df$id, function(id) length(extra_chains[[id]]), integer(1))
  extra_df$row_key <- vapply(extra_df$n_units, .dnmb_cct_extra_row_key, character(1))
  extra_df$y <- unname(row_y_map[extra_df$row_key])
  extra_df$lane_sugar <- vapply(extra_df$id, .dnmb_cct_extra_preferred_lane,
                                character(1), extra_chains = extra_chains)
  default_anchor <- stats::median(unname(sugar_lane_x))
  extra_df$anchor_x <- vapply(extra_df$id, .dnmb_cct_extra_anchor_x, numeric(1),
                              extra_chains = extra_chains,
                              extra_branches = extra_branches,
                              sugar_lane_x = sugar_lane_x,
                              default_anchor = default_anchor)
  target_stats <- vapply(seq_len(nrow(extra_df)), function(i) {
    .dnmb_cct_extra_target_stats(
      cascade_key = extra_df$cascade_key[i],
      substrate_id = extra_df$id[i],
      cascades = cascades,
      sugar_lane_x = sugar_lane_x,
      default_target = extra_df$anchor_x[i],
      extra_chains = extra_chains,
      extra_branches = extra_branches
    )
  }, numeric(3))
  # vapply over an empty sequence returns a 3x0 matrix with no dimnames,
  # which then breaks target_stats["center", ] below. Force the row names
  # so subsetting works for both empty and non-empty inputs.
  rownames(target_stats) <- c("center", "left", "right")
  extra_df$route_target_x <- target_stats["center", ]
  extra_df$route_left_x <- target_stats["left", ]
  extra_df$route_right_x <- target_stats["right", ]
  extra_df$route_span <- pmax(
    abs(extra_df$route_target_x - extra_df$anchor_x),
    extra_df$route_right_x - extra_df$route_left_x
  )
  extra_df$chain_width <- ifelse(
    extra_df$n_units > 0,
    (extra_df$n_units - 1) * step + 2 * sym_size,
    2 * sym_size
  )
  extra_df$label_width <- pmax(0.65, nchar(extra_df$label) * 0.11)
  extra_df$width_eff <- pmax(extra_df$chain_width, extra_df$label_width, 0.60 + 0.20 * extra_df$route_span)
  extra_df$x_center <- extra_df$anchor_x

  for (rk in unique(extra_df$row_key)) {
    ridx <- which(extra_df$row_key == rk)
    if (length(ridx) == 0) next
    row_df <- extra_df[ridx, , drop = FALSE]
    row_gap <- .dnmb_cct_extra_row_gap(rk, min_gap = min_gap)
    perms <- .dnmb_cct_small_permutations(seq_len(nrow(row_df)))
    best_df <- NULL
    best_cost <- Inf
    for (perm in perms) {
      cand_df <- .dnmb_cct_extra_pack_row(row_df[perm, , drop = FALSE], min_gap = row_gap)
      cand_cost <- .dnmb_cct_extra_row_cost(cand_df)
      if (cand_cost < best_cost) {
        best_cost <- cand_cost
        best_df <- cand_df
      }
    }
    extra_df$x_center[match(best_df$id, extra_df$id)] <- best_df$x_center
  }

  # Global coordinate-descent refinement:
  # re-optimize each row while holding others fixed to reduce cross-row crossings.
  row_keys <- unique(extra_df$row_key)
  if (length(row_keys) > 1) {
    max_iter <- 8L
    for (iter in seq_len(max_iter)) {
      improved <- FALSE
      current_cost <- .dnmb_cct_extra_global_cost(extra_df)
      for (rk in row_keys) {
        ridx <- which(extra_df$row_key == rk)
        row_df <- extra_df[ridx, , drop = FALSE]
        row_gap <- .dnmb_cct_extra_row_gap(rk, min_gap = min_gap)
        perms <- .dnmb_cct_small_permutations(seq_len(nrow(row_df)))
        best_all <- extra_df
        best_cost <- current_cost
        for (perm in perms) {
          cand_row <- .dnmb_cct_extra_pack_row(row_df[perm, , drop = FALSE], min_gap = row_gap)
          cand_all <- extra_df
          cand_all$x_center[match(cand_row$id, cand_all$id)] <- cand_row$x_center
          cand_cost <- .dnmb_cct_extra_global_cost(cand_all)
          if (cand_cost + 1e-9 < best_cost) {
            best_cost <- cand_cost
            best_all <- cand_all
          }
        }
        if (best_cost + 1e-9 < current_cost) {
          extra_df <- best_all
          current_cost <- best_cost
          improved <- TRUE
        }
      }
      if (!improved) break
    }
  }

  for (rk in unique(extra_df$row_key)) {
    if (!rk %in% c("4", "5", "6")) next
    ridx <- which(extra_df$row_key == rk)
    if (length(ridx) <= 1) next
    extra_df[ridx, ] <- .dnmb_cct_extra_stagger_row(
      extra_df[ridx, , drop = FALSE],
      y_step = switch(rk, `4` = 0.28, `5` = 0.20, `6` = 0.25, 0.20),
      min_sep = switch(rk, `4` = 0.22, `5` = 0.16, `6` = 0.20, 0.16)
    )
  }

  extra_df$x_start <- extra_df$x_center - (extra_df$n_units - 1) * step / 2
  extra_df$x_left <- extra_df$x_start - sym_size
  extra_df$x_right <- extra_df$x_start + (extra_df$n_units - 1) * step + sym_size
  extra_df
}

#' Draw a connected chain of SNFG sugar symbols
#'
#' Monomers are laid out left-to-right.
#' Bond types control the line style between consecutive units.
#'
#' @param x_start Left-most x position of the first monomer
#' @param y        y position (constant for a horizontal chain)
#' @param monomers Character vector of sugar types
#' @param bonds    Character vector of bond types (length = length(monomers) - 1).
#'   Recognised values: \code{"alpha-1,4"} (solid), \code{"beta-1,4"} (dashed),
#'   \code{"alpha-1,1"} (dotted), \code{"alpha-1,6"} (solid + vertical branch stub).
#' @param size Radius per symbol
#' @param gap  Extra spacing between symbols (plot units)
#' @return List of ggplot2 layers
#' @keywords internal
.dnmb_draw_sugar_chain_v2 <- function(x_start, y, monomers, bonds = NULL,
                                   size = 0.075, gap = 0.18, show_label = TRUE) {

  n_mono <- length(monomers)
  if (n_mono == 0L) return(list())

  step   <- size * 2 + gap
  xs     <- x_start + (seq_len(n_mono) - 1L) * step

  layers <- list()

  # ---- bond segments first (drawn below symbols) ----
  if (!is.null(bonds) && n_mono > 1L) {
    bonds <- rep_len(bonds, n_mono - 1L)

    for (i in seq_along(bonds)) {
      # Parse bond type: "alpha-1,4" -> anomeric = alpha, positions = 1,4
      bond_parts <- strsplit(bonds[i], "-")[[1L]]
      is_alpha   <- identical(bond_parts[1L], "alpha")
      lty        <- if (is_alpha) "solid" else "dashed"
      pos_parts  <- if (length(bond_parts) >= 2L) strsplit(bond_parts[2L], ",")[[1L]] else c("1","4")
      anom_sym   <- if (is_alpha) "\u03b1" else "\u03b2"

      seg_x1 <- xs[i] + size
      seg_x2 <- xs[i + 1L] - size

      seg_df <- data.frame(
        x    = seg_x1,
        xend = seg_x2,
        y    = y,
        yend = y
      )
      layers[[length(layers) + 1L]] <- ggplot2::geom_segment(
        data    = seg_df,
        mapping = ggplot2::aes(x = .data[["x"]], xend = .data[["xend"]],
                               y = .data[["y"]], yend = .data[["yend"]]),
        linetype  = lty,
        linewidth = 0.35,
        color     = "#555555",
        inherit.aes = FALSE
      )

      # Bond notation: compact label at midpoint above the linkage line
      mid_x <- (seg_x1 + seg_x2) / 2
      bond_label <- if (length(pos_parts) >= 2L) {
        paste0(anom_sym, pos_parts[1L], "\u2013", pos_parts[2L])
      } else {
        paste0(anom_sym, pos_parts[1L])
      }
      layers[[length(layers) + 1L]] <- ggplot2::geom_text(
        data    = data.frame(x = mid_x, y = y + size * 2.0,
                             label = bond_label),
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                               label = .data[["label"]]),
        size = 1.2, color = "#888888", hjust = 0.5, inherit.aes = FALSE
      )

      # alpha-1,6 gets a small vertical branch stub
      if (bonds[i] == "alpha-1,6") {
        branch_df <- data.frame(
          x    = (xs[i] + xs[i + 1L]) / 2,
          xend = (xs[i] + xs[i + 1L]) / 2,
          y    = y,
          yend = y + size * 1.2
        )
        layers[[length(layers) + 1L]] <- ggplot2::geom_segment(
          data    = branch_df,
          mapping = ggplot2::aes(x = .data[["x"]], xend = .data[["xend"]],
                                 y = .data[["y"]], yend = .data[["yend"]]),
          linewidth = 0.4,
          color     = "#333333",
          inherit.aes = FALSE
        )
      }
    }
  }

  # ---- monomer symbols on top ----
  for (i in seq_len(n_mono)) {
    if (isTRUE(show_label)) {
      sym_layers <- .dnmb_snfg_symbol_layers_v2(xs[i], y, monomers[i], size = size)
    } else {
      sym_layers <- .dnmb_snfg_symbol_layers_v2(xs[i], y, monomers[i], size = size, label = NULL)
    }
    layers     <- c(layers, sym_layers)
  }

  layers
}

#' Draw a GH cleavage mark between (or at the end of) sugar units
#'
#' @param x,y Position of the cleavage mark (midpoint between two bonds or
#'   at the terminal end)
#' @param gh_family Character label for the GH family, e.g. \code{"GH13"}
#' @param endo_exo  \code{"endo"} draws a red X, \code{"exo"} draws a red
#'   arrow pointing outward
#' @param size Visual size scaling factor
#' @return List of ggplot2 layers
#' @keywords internal
.dnmb_draw_gh_cleavage <- function(x, y, gh_family, endo_exo = "endo",
                                   size = 0.1) {

  layers <- list()

  if (endo_exo == "endo") {
    # Red X: two short crossing line segments
    d <- size * 0.8
    x_df <- data.frame(
      x    = c(x - d, x + d),
      xend = c(x + d, x - d),
      y    = c(y - d, y + d),
      yend = c(y + d, y - d)
    )
    layers[[1L]] <- ggplot2::geom_segment(
      data    = x_df,
      mapping = ggplot2::aes(x = .data[["x"]], xend = .data[["xend"]],
                             y = .data[["y"]], yend = .data[["yend"]]),
      linewidth = 1.2,
      color     = "#ED1C24",
      inherit.aes = FALSE
    )
  } else {
    # Exo: red arrow pointing right (at chain end)
    arrow_df <- data.frame(
      x    = x - size,
      xend = x + size,
      y    = y,
      yend = y
    )
    layers[[1L]] <- ggplot2::geom_segment(
      data    = arrow_df,
      mapping = ggplot2::aes(x = .data[["x"]], xend = .data[["xend"]],
                             y = .data[["y"]], yend = .data[["yend"]]),
      arrow     = ggplot2::arrow(length = ggplot2::unit(size * 40, "pt"),
                                 type = "closed"),
      linewidth = 1.0,
      color     = "#ED1C24",
      inherit.aes = FALSE
    )
  }

  # GH family label below the mark
  label_df <- data.frame(x = x, y = y - size * 2.5, label = gh_family)
  layers[[length(layers) + 1L]] <- ggplot2::geom_text(
    data    = label_df,
    mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                           label = .data[["label"]]),
    size    = 3,
    color   = "#ED1C24",
    fontface = "bold",
    inherit.aes = FALSE
  )

  layers
}

# ---------------------------------------------------------------------------
# Scissors grob for GH cleavage marks
# ---------------------------------------------------------------------------

#' Draw scissors symbol (two crossed blades) for GH cleavage
#' @param x,y Position of the scissors
#' @param size Scale factor
#' @param angle Rotation angle in degrees (0 = horizontal cut)
#' @return List of ggplot2 layers (two curved blade segments + pivot)
#' @keywords internal
.dnmb_scissors_grob_v2 <- function(x, y, size = 0.12, angle = 0) {
  d <- size
  rad <- angle * pi / 180
  # Two blade lines crossing at (x, y)
  # Blade 1: top-left to bottom-right
  b1 <- data.frame(
    x    = x + d * cos(rad + pi * 0.75),
    xend = x + d * cos(rad - pi * 0.25),
    y    = y + d * sin(rad + pi * 0.75),
    yend = y + d * sin(rad - pi * 0.25)
  )
  # Blade 2: bottom-left to top-right
  b2 <- data.frame(
    x    = x + d * cos(rad + pi * 1.25),
    xend = x + d * cos(rad + pi * 0.25),
    y    = y + d * sin(rad + pi * 1.25),
    yend = y + d * sin(rad + pi * 0.25)
  )
  # Handle rings (small circles at blade ends)
  r <- size * 0.3
  ring1 <- .dnmb_snfg_polygon_coords_v2("circle", b1$x, b1$y, r)
  ring2 <- .dnmb_snfg_polygon_coords_v2("circle", b2$x, b2$y, r)
  ring1$group <- "r1"; ring2$group <- "r2"

  list(
    ggplot2::geom_segment(
      data = b1, mapping = ggplot2::aes(x = .data[["x"]], xend = .data[["xend"]],
                                         y = .data[["y"]], yend = .data[["yend"]]),
      linewidth = 0.7, color = "#ED1C24", inherit.aes = FALSE),
    ggplot2::geom_segment(
      data = b2, mapping = ggplot2::aes(x = .data[["x"]], xend = .data[["xend"]],
                                         y = .data[["y"]], yend = .data[["yend"]]),
      linewidth = 0.7, color = "#ED1C24", inherit.aes = FALSE),
    ggplot2::geom_polygon(
      data = ring1, mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                            group = .data[["group"]]),
      fill = NA, color = "#ED1C24", linewidth = 0.4, inherit.aes = FALSE),
    ggplot2::geom_polygon(
      data = ring2, mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                            group = .data[["group"]]),
      fill = NA, color = "#ED1C24", linewidth = 0.8, inherit.aes = FALSE)
  )
}

# ---------------------------------------------------------------------------
# 3-Zone CAZy Carbon Transport Map — data extraction helpers
# ---------------------------------------------------------------------------

#' Extract GH enzymes with locus_tag, gene_name, family from genbank_table
#' @keywords internal
.dnmb_cct_3zone_extract_gh <- function(genbank_table) {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE, row.names = NULL)
  family_col <- .dnmb_pick_column(tbl, c("dbCAN_family_id", "family_id"))
  if (is.null(family_col)) return(NULL)
  tbl <- tbl[!is.na(tbl[[family_col]]) & nzchar(tbl[[family_col]]), , drop = FALSE]
  if (!nrow(tbl)) return(NULL)
  tbl$family <- as.character(tbl[[family_col]])
  # Keep only GH families

  gh <- tbl[grepl("^GH", tbl$family), , drop = FALSE]
  if (!nrow(gh)) return(NULL)
  lt_col <- .dnmb_pick_column(gh, c("locus_tag", "old_locus_tag"))
  gn_col <- .dnmb_pick_column(gh, c("gene", "gene_name"))
  data.frame(
    locus_tag = if (!is.null(lt_col)) gh[[lt_col]] else paste0("gene_", seq_len(nrow(gh))),
    gene_name = if (!is.null(gn_col)) gh[[gn_col]] else NA_character_,
    gh_family = gh$family,
    stringsAsFactors = FALSE
  )
}

#' Extract transporter genes from genbank_table (GapMindCarbon transport steps)
#' @keywords internal
.dnmb_cct_3zone_extract_transporters <- function(genbank_table, output_dir = NULL) {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE, row.names = NULL)
  path_col <- "GapMindCarbon_pathway_id"; step_col <- "GapMindCarbon_step_id"
  result <- NULL

  # Source 0: robust transport reconstruction from rules + steps
  if (!is.null(output_dir)) {
    strong_tr <- .dnmb_cct_transport_extract(output_dir)
    if (!is.null(strong_tr) && nrow(strong_tr) > 0) {
      keep_tr <- mapply(
        .dnmb_cct_is_transport_like_step,
        step_id = strong_tr$step,
        gene_name = if ("gene_name" %in% names(strong_tr)) strong_tr$gene_name else strong_tr$step,
        is_pts = if ("is_pts" %in% names(strong_tr)) strong_tr$is_pts else FALSE
      )
      strong_tr <- strong_tr[keep_tr, , drop = FALSE]
    }
    if (!is.null(strong_tr) && nrow(strong_tr) > 0) {
      result <- data.frame(
        locus_tag = strong_tr$locus_tag,
        gene_name = if ("gene_name" %in% names(strong_tr)) strong_tr$gene_name else strong_tr$step,
        pathway = strong_tr$pathway,
        step = strong_tr$step,
        confidence = strong_tr$confidence,
        step_score = if ("score" %in% names(strong_tr)) strong_tr$score else NA_real_,
        is_pts = if ("is_pts" %in% names(strong_tr)) strong_tr$is_pts else FALSE,
        stringsAsFactors = FALSE
      )
    }
  }

  # Source 1: genbank_table GapMindCarbon columns
  if (path_col %in% names(tbl)) {
    tbl2 <- tbl[!is.na(tbl[[path_col]]) & nzchar(tbl[[path_col]]), , drop = FALSE]
    if (nrow(tbl2) > 0) {
      is_t <- grepl("transport|PTS|permease|ABC|SBP|MFS|EIICB|IIBC|EIIBC|EIIA|crr|porter|uptake", tbl2[[step_col]], ignore.case = TRUE) |
              grepl("transport", tbl2[[path_col]], ignore.case = TRUE) |
              grepl("^pts[A-Z]|^mal[EF]IICB|^fru[A-Z]*-ABC", tbl2[[step_col]], ignore.case = FALSE)
      trans <- tbl2[is_t, , drop = FALSE]
      if (nrow(trans) > 0) {
        lt_col <- .dnmb_pick_column(trans, c("locus_tag","old_locus_tag"))
        gn_col <- .dnmb_pick_column(trans, c("gene","gene_name"))
        cc <- "GapMindCarbon_confidence"
        res1 <- data.frame(
          locus_tag  = if (!is.null(lt_col)) trans[[lt_col]] else paste0("t_", seq_len(nrow(trans))),
          gene_name  = if (!is.null(gn_col)) trans[[gn_col]] else NA_character_,
          pathway = trans[[path_col]], step = trans[[step_col]],
          confidence = if (cc %in% names(trans)) trans[[cc]] else "medium",
          step_score = if ("GapMindCarbon_step_score" %in% names(trans)) as.numeric(trans[["GapMindCarbon_step_score"]]) else NA_real_,
          is_pts = grepl("pts|PTS|EIIC|EIIB|crr", trans[[step_col]], ignore.case = TRUE),
          stringsAsFactors = FALSE, row.names = NULL)
        result <- if (is.null(result)) res1 else rbind(result, res1)
      }
    }
  }
  # Source 2: aa.sum.steps fallback
  if (is.null(result) && !is.null(output_dir)) {
    for (f in c(file.path(output_dir, "dnmb_module_gapmindcarbon", "aa.sum.steps"),
                file.path(output_dir, "dnmb_module_gapmind_carbon", "aa.sum.steps"),
                file.path(output_dir, "aa.sum.steps"))) {
      if (file.exists(f)) {
        st <- utils::read.delim(f, stringsAsFactors = FALSE, row.names = NULL)
        if (all(c("pathway","step","locusId") %in% names(st))) {
          is_t <- grepl("transport|PTS|pts|permease|ABC|SBP|MFS", st$step, ignore.case = TRUE)
          has_l <- !is.na(st$locusId) & nzchar(st$locusId)
          ts <- st[is_t & has_l, , drop = FALSE]
          ts <- ts[!duplicated(ts$locusId), , drop = FALSE]
          if (nrow(ts) > 0) {
            sc <- if ("score" %in% names(ts)) ts$score else 0
            result <- data.frame(
              locus_tag = ts$locusId,
              gene_name = if ("sysName" %in% names(ts)) ts$sysName else NA_character_,
              pathway = ts$pathway, step = ts$step,
              confidence = dplyr::case_when(sc >= 2 ~ "high", sc >= 1 ~ "medium", TRUE ~ "low"),
              step_score = sc,
              is_pts = grepl("pts|PTS|EIIC|EIIB|crr", ts$step, ignore.case = TRUE),
              stringsAsFactors = FALSE, row.names = NULL)
          }
        }
        break
      }
    }
  }
  # Source 3: CGC transporter proteins (gene_type = TP/STP)
  cgc_type_col <- intersect(c("dbCAN_dbcan_cgc_gene_type", "dbcan_cgc_gene_type"), names(tbl))
  cgc_id_col <- intersect(c("dbCAN_dbcan_cgc_id", "dbcan_cgc_id"), names(tbl))
  pul_sub_col <- intersect(c("dbCAN_dbcan_pul_substrate", "dbcan_pul_substrate"), names(tbl))
  if (length(cgc_type_col) > 0) {
    cgc_tp <- tbl[!is.na(tbl[[cgc_type_col[1]]]) &
                  grepl("^TP$|^STP$|^TC$", tbl[[cgc_type_col[1]]], ignore.case = TRUE), , drop = FALSE]
    if (nrow(cgc_tp) > 0) {
      lt_col <- .dnmb_pick_column(cgc_tp, c("locus_tag", "old_locus_tag"))
      gn_col <- .dnmb_pick_column(cgc_tp, c("gene", "gene_name"))
      prod_col <- .dnmb_pick_column(cgc_tp, c("product"))
      # Map CGC substrate to pathway
      pul_map <- .dnmb_cct_pul_substrate_to_cs()
      cgc_pathway <- rep("unknown", nrow(cgc_tp))
      if (length(pul_sub_col) > 0) {
        for (ci in seq_len(nrow(cgc_tp))) {
          sub_val <- tolower(trimws(as.character(cgc_tp[[pul_sub_col[1]]][ci])))
          if (!is.na(sub_val) && nzchar(sub_val)) {
            mapped <- pul_map[strsplit(sub_val, ";\\s*")[[1]][1]]
            if (!is.na(mapped)) cgc_pathway[ci] <- tolower(mapped)
          }
        }
      }
      step_name <- if (!is.null(prod_col)) {
        sub("^hypothetical protein$", "CGC-TP", cgc_tp[[prod_col]], ignore.case = TRUE)
      } else rep("CGC-TP", nrow(cgc_tp))
      res3 <- data.frame(
        locus_tag = if (!is.null(lt_col)) cgc_tp[[lt_col]] else paste0("cgc_", seq_len(nrow(cgc_tp))),
        gene_name = if (!is.null(gn_col)) cgc_tp[[gn_col]] else step_name,
        pathway = cgc_pathway,
        step = step_name,
        confidence = "medium",
        step_score = 1.5,
        is_pts = FALSE,
        stringsAsFactors = FALSE
      )
      res3 <- res3[!is.na(res3$locus_tag) & nzchar(res3$locus_tag), , drop = FALSE]
      result <- if (is.null(result)) res3 else rbind(result, res3)
    }
  }

  if (is.null(result) || nrow(result) == 0) return(result)

  result <- result[!is.na(result$locus_tag) & nzchar(result$locus_tag), , drop = FALSE]
  result$gene_name <- ifelse(is.na(result$gene_name) | !nzchar(result$gene_name), result$step, result$gene_name)
  same_as_locus <- !is.na(result$gene_name) & result$gene_name == result$locus_tag
  result$gene_name[same_as_locus] <- result$step[same_as_locus]
  result$conf_rank <- c(none = 0L, low = 1L, medium = 2L, high = 3L)[result$confidence]
  result$conf_rank[is.na(result$conf_rank)] <- 0L
  result$kind <- mapply(.dnmb_cct_transporter_kind, result$step, result$gene_name, result$is_pts)
  result <- result[order(result$pathway, -result$conf_rank, result$is_pts, result$step, result$locus_tag), , drop = FALSE]
  result <- result[!duplicated(result[, c("pathway", "step", "locus_tag")]), , drop = FALSE]
  result
}

#' Define glycan degradation cascades
#' Uses GH prefix matching (GH13 matches GH13_20, GH13_31 etc.)
#' @keywords internal
.dnmb_cct_3zone_glycan_cascades <- function() {
  list(
    starch = list(
      name = "Starch",
      chains = list(
        c("glucose", "glucose", "glucose", "glucose", "glucose"),
        c("glucose", "glucose", "glucose"),
        c("glucose", "glucose"),
        c("glucose")
      ),
      bonds  = "alpha-1,4",
      gh_prefixes = c("GH13", "GH14", "GH15", "GH57", "GH77", "GH31"),
      product = "Glucose"
    ),
    cellulose = list(
      name = "Cellulose",
      chains = list(
        c("glucose", "glucose", "glucose", "glucose"),
        c("glucose", "glucose"),
        c("glucose")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH5", "GH6", "GH7", "GH9", "GH44", "GH45", "GH48", "GH1"),
      product = "Glucose"
    ),
    xylan = list(
      name = "Xylan",
      chains = list(
        c("xylose", "xylose", "xylose", "xylose"),
        c("xylose", "xylose"),
        c("xylose")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH10", "GH11", "GH30", "GH43", "GH8"),
      product = "Xylose"
    ),
    chitin = list(
      name = "Chitin",
      chains = list(
        c("glcnac", "glcnac", "glcnac", "glcnac"),
        c("glcnac", "glcnac"),
        c("glcnac")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH18", "GH19", "GH20"),
      product = "NAG"
    ),
    pectin = list(
      name = "Pectin",
      chains = list(
        c("gala", "gala", "gala", "gala"),
        c("gala", "gala"),
        c("gala")
      ),
      bonds  = "alpha-1,4",
      gh_prefixes = c("GH28", "GH53", "GH78", "GH88", "GH105"),
      product = "GalA"
    ),
    mannan = list(
      name = "Mannan",
      chains = list(
        c("mannose", "mannose", "mannose", "mannose"),
        c("mannose", "mannose"),
        c("mannose")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH5", "GH26", "GH113", "GH134", "GH130"),
      product = "Mannose"
    ),
    sucrose = list(
      name = "Sucrose",
      chains = list(
        c("glucose", "fructose")
      ),
      bonds  = "alpha-1,2",
      gh_prefixes = c("GH32", "GH68", "GH100"),
      product = "Glucose",
      products = c("glucose", "fructose")
    ),
    maltose = list(
      name = "Maltose",
      chains = list(
        c("glucose", "glucose"),
        c("glucose")
      ),
      bonds  = "alpha-1,4",
      gh_prefixes = c("GH13", "GH4", "GH65"),
      product = "Glucose"
    ),
    lactose = list(
      name = "Lactose",
      chains = list(
        c("galactose", "glucose")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH1", "GH2", "GH35", "GH42"),
      product = "Galactose",
      products = c("galactose", "glucose")
    ),
    cellobiose = list(
      name = "Cellobiose",
      chains = list(
        c("glucose", "glucose"),
        c("glucose")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH1", "GH3", "GH4"),
      product = "Glucose"
    ),
    trehalose = list(
      name = "Trehalose",
      chains = list(
        c("glucose", "glucose"),
        c("glucose")
      ),
      bonds  = "alpha-1,1",
      gh_prefixes = c("GH37", "GH65", "GH13"),
      product = "Glucose"
    ),
    raffinose = list(
      name = "Raffinose",
      chains = list(
        c("galactose", "glucose", "fructose"),
        c("glucose", "fructose")
      ),
      bonds  = c("alpha-1,6", "alpha-1,2"),
      level_bonds = list(
        c("alpha-1,6", "alpha-1,2"),  # Gal-α1,6-Glc-α1,2-Fru
        c("alpha-1,2")                 # Glc-α1,2-Fru (=Sucrose)
      ),
      gh_prefixes = c("GH27", "GH36"),
      product = "Galactose",
      products = c("galactose", "sucrose")
    ),
    inulin = list(
      name = "Inulin",
      chains = list(
        c("fructose", "fructose", "fructose", "fructose"),
        c("fructose", "fructose"),
        c("fructose")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH32", "GH91"),
      product = "Fructose"
    ),
    betaglucan = list(
      name = "beta-Glucan",
      chains = list(
        c("glucose", "glucose", "glucose", "glucose"),
        c("glucose", "glucose"),
        c("glucose")
      ),
      bonds  = "beta-1,3",
      gh_prefixes = c("GH16", "GH17", "GH55", "GH64", "GH81", "GH128", "GH158"),
      product = "Glucose"
    ),
    arabinoxylan = list(
      name = "Arabinoxylan",
      chains = list(
        c("xylose", "xylose", "xylose", "xylose"),
        c("xylose", "xylose"),
        c("xylose")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH10", "GH11", "GH43", "GH51", "GH62"),
      product = "Xylose"
    ),
    agarose = list(
      name = "Agarose",
      chains = list(
        c("galactose", "galactose", "galactose", "galactose"),
        c("galactose", "galactose"),
        c("galactose")
      ),
      bonds  = "beta-1,4",
      gh_prefixes = c("GH16", "GH50", "GH86", "GH96", "GH117", "GH118"),
      product = "Galactose"
    )
  )
}

#' Match GH enzymes to cascade gh_prefixes (prefix matching for subfamilies)
#' @keywords internal
.dnmb_cct_3zone_match_gh <- function(gh_enzymes, gh_prefixes) {
  if (is.null(gh_enzymes) || !is.data.frame(gh_enzymes) || nrow(gh_enzymes) == 0) {
    return(data.frame(gh_family = character(0), locus_tag = character(0),
                      gene_name = character(0), stringsAsFactors = FALSE, row.names = NULL))
  }
  # Build regex: GH13 should match GH13, GH13_20, GH13_31 etc.
  patterns <- paste0("^(", paste(gh_prefixes, collapse = "|"), ")(\\b|_|$)")
  matched <- grepl(patterns, gh_enzymes$gh_family)
  gh_enzymes[matched, , drop = FALSE]
}

#' Classify GH family as endo-acting (cuts mid-chain) or exo-acting (cuts end)
#' @keywords internal
.dnmb_gh_is_endo <- function(gh_family) {
  base <- sub("(_.*|$)", "", gh_family)
  endo <- c("GH5","GH6","GH7","GH9","GH10","GH11","GH13","GH18","GH19",
            "GH23","GH28","GH44","GH45","GH48","GH26","GH8")
  base %in% endo
}

# ---------------------------------------------------------------------------
# Sugar metabolism node positions (gapmind_aa pattern)
# ---------------------------------------------------------------------------
#' @keywords internal
.dnmb_cct_sugar_nodes <- function(cs_xs = NULL, tca_cx_override = NULL) {
  n <- function(id, x, y, label, type, sugar_type = "generic") {
    if (is.na(id) || !nzchar(id)) id <- paste0("node_", x, "_", y)
    data.frame(id = id, x = x, y = y, label = label, type = type,
               sugar_type = sugar_type, stringsAsFactors = FALSE, row.names = NULL)
  }
  # ====================================================================
  # TOP-DOWN LAYOUT (landscape): Extracellular(top) -> Membrane -> Cytoplasm(bottom)
  # y = 9~15  : Extracellular (handled by rendering loop)
  # y = 8     : Membrane (horizontal line)
  # y = 0~7.5 : Cytoplasm
  # ====================================================================
  s  <- 0.5   # uniform vertical step for backbone

  # Carbon sources spread horizontally just below membrane
  # Accept external cs_xs (named vector) or fall back to even spacing
  all_cs <- c("Maltose","Cellobiose","Galactose","Trehalose","Lactose",
              "Mannose","NAG","Glucosamine","Fructose","Sucrose",
              "Mannitol","Glycerol","Xylose","Arabinose","Ribose",
              "Fucose","Rhamnose","Gluconate")
  # Family base positions: compact around backbone (bx=10)
  #   Glc family left of backbone, Fru/Man right, PPP/Xyl further right
  fam_cs_base <- c(Glc=7, GlcNAc=8, Gal=9, Fru=11, Man=12, Xyl=13.5, Other=15)
  dx_cs <- 0.5
  if (is.null(cs_xs)) {
    # Order: Maltose Cellobiose Galactose Trehalose Lactose Mannose NAG Glucosamine
    #        Fructose Sucrose Mannitol Glycerol Xylose Arabinose Ribose Fucose Rhamnose Gluconate
    cs_xs_default <- c(
      fam_cs_base["Glc"], fam_cs_base["Glc"]+dx_cs, fam_cs_base["Gal"], fam_cs_base["Glc"]+2*dx_cs,
      fam_cs_base["Gal"]+dx_cs, fam_cs_base["Man"], fam_cs_base["GlcNAc"], fam_cs_base["GlcNAc"]+dx_cs,
      fam_cs_base["Fru"], fam_cs_base["Fru"]+dx_cs, fam_cs_base["Man"]+dx_cs, fam_cs_base["Other"],
      fam_cs_base["Xyl"], fam_cs_base["Xyl"]+dx_cs, fam_cs_base["Xyl"]+2*dx_cs,
      fam_cs_base["Other"]+dx_cs, fam_cs_base["Other"]+2*dx_cs, fam_cs_base["Other"]+3*dx_cs)
    cs_xs <- stats::setNames(cs_xs_default, all_cs)
  }
  cs_map <- cs_xs

  # Backbone top y -- just below membrane
  by <- 7.5  # Glucose starts here
  cs_y <- by + 0.5  # carbon source row y

  # ---- Backbone centered, entries at their family x ----
  bx <- 10.0                   # Glycolysis backbone at center x
  px <- bx + 3.0               # PPP column (x=13)
  rx <- bx + 4.0               # Right PPP branch (x=14)
  compact_anchor <- function(ids, target_x, fallback, pull_left = 0.70, pull_right = 0.52) {
    xs <- unname(cs_map[ids])
    xs <- xs[is.finite(xs)]
    mean_x <- if (length(xs) > 0) mean(xs) else fallback
    pull <- if (mean_x < target_x) pull_left else pull_right
    .dnmb_cct_snap_to_grid(mean_x * (1 - pull) + target_x * pull, step = 0.5)
  }
  fam_cs <- c(
    Glc = compact_anchor(c("Maltose", "Cellobiose", "Trehalose"), bx, fallback = fam_cs_base["Glc"], pull_left = 0.55, pull_right = 0.45),
    GlcNAc = compact_anchor(c("NAG", "Glucosamine"), bx, fallback = fam_cs_base["GlcNAc"], pull_left = 0.50, pull_right = 0.45),
    Gal = compact_anchor(c("Galactose", "Lactose"), bx, fallback = fam_cs_base["Gal"], pull_left = 0.50, pull_right = 0.40),
    Fru = compact_anchor(c("Fructose", "Sucrose"), bx, fallback = fam_cs_base["Fru"], pull_left = 0.40, pull_right = 0.40),
    Man = compact_anchor(c("Mannose", "Mannitol"), bx, fallback = fam_cs_base["Man"], pull_left = 0.45, pull_right = 0.40),
    Xyl = compact_anchor(c("Xylose", "Arabinose", "Ribose"), px, fallback = fam_cs_base["Xyl"], pull_left = 0.35, pull_right = 0.40),
    Other = compact_anchor(c("Glycerol", "Fucose", "Rhamnose", "Gluconate"), bx + 2.0, fallback = fam_cs_base["Other"], pull_left = 0.35, pull_right = 0.45)
  )
  xyl_lane <- .dnmb_cct_snap_to_grid(fam_cs["Xyl"], step = 0.5)
  arab_lane <- .dnmb_cct_snap_to_grid(fam_cs["Xyl"] + dx_cs, step = 0.5)
  rib_lane <- .dnmb_cct_snap_to_grid(fam_cs["Xyl"] + 2 * dx_cs, step = 0.5)
  ppp_left <- .dnmb_cct_snap_to_grid(max(px, xyl_lane - 0.5), step = 0.5)
  ppp_mid <- .dnmb_cct_snap_to_grid(max(px + 0.5, arab_lane), step = 0.5)






  # ---- TCA cycle: under TCA-entry substrates if override given ----
  if (!is.null(tca_cx_override)) {
    tca_cx <- tca_cx_override
  } else {
    tca_cx <- bx
  }
  # TCA nodes: exact circle positions (NO grid snap — exception for TCA)
  tca_cy <- 0.5;  tca_r <- 1.5
  tca_ang <- pi/2 - (0:7) * (2 * pi / 8)  # 45° intervals, OAA at top
  tx <- tca_cx + tca_r * cos(tca_ang)      # exact circle x
  ty <- tca_cy + tca_r * sin(tca_ang)      # exact circle y





  # Helper: get x of a carbon source, fall back to family x
  csx <- function(cs_name, fallback) {
    if (cs_name %in% names(cs_map) && is.finite(cs_map[cs_name])) {
      .dnmb_cct_snap_to_grid(unname(cs_map[cs_name]), step = 0.5)
    } else {
      .dnmb_cct_snap_to_grid(fallback, step = 0.5)
    }
  }

  rbind(
    # ---- FULL GLYCOLYSIS BACKBONE (vertical column at bx) ----
    n("Glucose",      bx,    by,        "Glucose",       "backbone", "glucose"),
    n("Glc-6-P",      bx,    by-1*s,    "Glc-6-P",       "backbone", "glucose"),
    n("Fru-6-P",      bx,    by-2*s,    "Fru-6-P",       "backbone", "fructose"),
    n("Fru-1,6-BP",   bx,    by-3*s,    "Fru-1,6-BP",    "backbone", "fructose"),
    n("DHAP",         bx-s,  by-4*s,    "DHAP",           "backbone", "intermediate"),
    n("GA3P",         bx,    by-4*s,    "GA3P",           "backbone", "intermediate"),
    n("1,3-BPG",      bx,    by-5*s,    "1,3-BPG",        "backbone", "intermediate"),
    n("3-PG",         bx,    by-6*s,    "3-PG",            "backbone", "intermediate"),
    n("2-PG",         bx,    by-7*s,    "2-PG",            "backbone", "intermediate"),
    n("PEP",          bx,    by-8*s,    "PEP",             "backbone", "intermediate"),
    n("Pyruvate",     bx,    by-9*s,    "Pyruvate",        "backbone", "intermediate"),
    n("Acetyl-CoA",   bx,    by-10*s,   "Acetyl-CoA",      "backbone", "intermediate"),

    # ---- TCA CYCLE (perfect circle, gapmind_aa style) ----
    n("OAA",       tx[1], ty[1], "OAA",       "tca", "intermediate"),
    n("Citrate",   tx[2], ty[2], "Citrate",   "tca", "intermediate"),
    n("Isocit",    tx[3], ty[3], "Isocit",    "tca", "intermediate"),
    n("AKG",       tx[4], ty[4], "\u03b1-KG", "tca", "intermediate"),
    n("SucCoA",    tx[5], ty[5], "Suc-CoA",   "tca", "intermediate"),
    n("Succinate", tx[6], ty[6], "Succ",      "tca", "intermediate"),
    n("Fumarate",  tx[7], ty[7], "Fum",       "tca", "intermediate"),
    n("Malate",    tx[8], ty[8], "Mal",       "tca", "intermediate"),

    # ---- PHOSPHORYLATION ROW: y=7.5, each at its parent carbon source x ----
    # First phosphorylated form — placed directly below parent carbon source
    n("Glc-6-P-entry", csx("Maltose", fam_cs["Glc"]),    7.5, "Glc-6-P",     "entry_intermediate", "glucose"),
    n("Gal-1-P",       csx("Galactose", fam_cs["Gal"]),  7.5, "Gal-1-P",     "entry_intermediate", "galactose"),
    n("GlcNAc-6-P",    csx("NAG", fam_cs["GlcNAc"]),     7.5, "GlcNAc-6-P",  "entry_intermediate", "glcnac"),
    n("Man-6-P",       csx("Mannose", fam_cs["Man"]),     7.5, "Man-6-P",     "entry_intermediate", "mannose"),
    n("Fru-6-P-entry", csx("Fructose", fam_cs["Fru"]),   7.5, "Fru-6-P",     "entry_intermediate", "fructose"),
    n("Xu-5-P-entry",  csx("Xylose", xyl_lane),          7.5, "Xu-5-P",      "entry_intermediate", "xylose"),
    n("R-5-P-entry",   csx("Ribose", rib_lane),          7.5, "R-5-P",       "entry_intermediate", "ribose"),
    n("Glycerol-3-P",  csx("Glycerol", fam_cs["Other"]-dx_cs), 7.5, "Glycerol-3-P", "entry_intermediate", "glycerol"),
    # Second step intermediates (y=7.0) — also aligned to parent carbon source x
    n("Glc-1-P",      csx("Maltose", fam_cs["Glc"]),     7.0, "Glc-1-P",     "entry_intermediate", "glucose"),
    n("Trehalose-6-P", csx("Trehalose", fam_cs["Glc"]+dx_cs), 7.0, "Tre-6-P","entry_intermediate", "glucose"),
    n("UDP-Gal",      csx("Galactose", fam_cs["Gal"])+dx_cs, 7.0, "UDP-Gal", "entry_intermediate", "galactose"),
    n("UDP-Glc",      csx("Galactose", fam_cs["Glc"])+2*dx_cs, 7.0, "UDP-Glc","entry_intermediate", "glucose"),
    n("GlcN-6-P",     csx("NAG", fam_cs["GlcNAc"]),      7.0, "GlcN-6-P",   "entry_intermediate", "glucosamine"),
    n("Fru-1-P",      csx("Fructose", fam_cs["Fru"]),    7.0, "Fru-1-P",     "entry_intermediate", "fructose"),
    n("Mannitol-1-P", csx("Mannitol", fam_cs["Man"]+dx_cs), 7.0, "Mtl-1-P",  "entry_intermediate", "mannose"),
    n("Xylulose",     csx("Xylose", xyl_lane),           7.0, "Xylulose",    "entry_intermediate", "xylose"),
    n("Ribulose",     csx("Arabinose", arab_lane),        7.0, "Ribulose",   "entry_intermediate", "arabinose"),
    n("Ribulose-5-P", csx("Arabinose", arab_lane),        6.5, "Ribu-5-P",   "entry_intermediate", "arabinose"),
    n("Fuculose",     csx("Fucose", fam_cs["Other"]),     7.0, "Fuculose",    "entry_intermediate", "fucose"),
    n("Fuculose-1-P", csx("Fucose", fam_cs["Other"]),     6.5, "Fuc-1-P",    "entry_intermediate", "fucose"),
    n("Rhamnulose",    csx("Rhamnose", fam_cs["Other"]+dx_cs), 7.0, "Rhamnulose","entry_intermediate", "rhamnose"),
    n("Rhamnulose-1-P",csx("Rhamnose", fam_cs["Other"]+dx_cs), 6.5, "Rha-1-P",  "entry_intermediate", "rhamnose"),

    # ---- PPP (right of backbone) — SNFG sugar types for sugar-derived metabolites ----
    n("6-PGL",     ppp_left, 7.0, "6-PGL",     "ppp", "glucose"),     # glucono-δ-lactone
    n("6-PG",      ppp_left, 6.5, "6-PG",      "ppp", "gluconate"),   # 6-phosphogluconate
    n("Ru-5-P",    ppp_mid,  6.0, "Ru-5-P",    "ppp", "ribose"),      # ribulose-5-P
    n("R-5-P",     rib_lane, 5.5, "R-5-P",     "ppp", "ribose"),      # ribose-5-P
    n("Xu-5-P",    xyl_lane, 5.5, "Xu-5-P",    "ppp", "xylose"),      # xylulose-5-P
    n("S-7-P",     ppp_mid,  5.0, "S-7-P",     "ppp", "intermediate"), # sedoheptulose (7C, no SNFG)
    n("E-4-P",     ppp_mid,  4.5, "E-4-P",     "ppp", "intermediate"), # erythrose (4C, no SNFG)

    # ---- ED PATHWAY ----
    n("KDG",      px + 0.8, 6.5, "KDG",  "ed", "gluconate"),     # 2-keto-3-deoxy-gluconate
    n("KDPG",     px+0.5, 6.0, "KDPG", "ed", "gluconate"),     # 2-keto-3-deoxy-6-phosphogluconate

    # ---- PYRUVATE BRANCHES (fermentation products) ----
    n("Formate",   bx-1.5, 3.0, "Formate",  "pyruvate_branch", "organic_acid"),
    n("Lactate",   bx-1.5, 2.5, "Lactate",  "pyruvate_branch", "organic_acid"),
    n("Acetyl-P",  bx+1.5, 2.0, "Acetyl-P", "pyruvate_branch", "intermediate"),
    n("Acetate",   bx+2.0, 2.0, "Acetate",  "pyruvate_branch", "organic_acid"),
    n("Acetaldehyde", bx+1.5, 1.5, "AcAld", "pyruvate_branch", "intermediate"),
    n("Ethanol",   bx+2.0, 1.5, "Ethanol",  "pyruvate_branch", "organic_acid"),

    # ---- GLYOXYLATE SHUNT ----
    n("Glyoxylate", tca_cx, tca_cy, "Glx", "tca", "intermediate"),

    # ---- CARBON SOURCES: dynamic, only present substrates ----
    .dnmb_cct_carbon_source_nodes(cs_xs, cs_y, n)
  )
}

#' Orthogonal (L-shaped) edges for the simplified cytoplasm
#' Each edge goes horizontal then vertical (or vice versa)
#' @keywords internal
.dnmb_cct_3zone_cyto_edges <- function(nodes) {
  node_x <- stats::setNames(nodes$x, nodes$id)
  node_y <- stats::setNames(nodes$y, nodes$id)

  e <- function(from, to, pathway = "backbone", curvature = 0) {
    if (!from %in% names(node_x) || !to %in% names(node_x)) return(NULL)
    data.frame(from = from, to = to,
               x = unname(node_x[from]), y = unname(node_y[from]),
               xend = unname(node_x[to]), yend = unname(node_y[to]),
               pathway = pathway, curvature = curvature,
               stringsAsFactors = FALSE, row.names = NULL)
  }

  # Pathway colors for each substrate family (gapmind_aa-style)
  # Stored as attribute for the rendering function to use

  edges <- list(
    # Full glycolysis backbone (straight vertical)
    e("Glucose", "Glc-6-P", "backbone"),
    e("Glc-6-P", "Fru-6-P", "backbone"),
    e("Fru-6-P", "Fru-1,6-BP", "backbone"),
    e("Fru-1,6-BP", "DHAP", "backbone"),
    e("Fru-1,6-BP", "GA3P", "backbone"),
    e("DHAP", "GA3P", "backbone", 0.15),
    e("GA3P", "1,3-BPG", "backbone"),
    e("1,3-BPG", "3-PG", "backbone"),
    e("3-PG", "2-PG", "backbone"),
    e("2-PG", "PEP", "backbone"),
    e("PEP", "Pyruvate", "backbone"),
    e("Pyruvate", "Acetyl-CoA", "backbone"),
    # AcCoA → Citrate (with OAA); OAA also from Malate
    e("Acetyl-CoA", "Citrate", "backbone", -0.15),
    e("OAA", "Citrate", "tca", 0.15),

    # TCA cycle edges (clockwise, curved to follow circle)
    e("Citrate", "Isocit", "tca", -0.3),
    e("Isocit", "AKG", "tca", -0.3),
    e("AKG", "SucCoA", "tca", -0.3),
    e("SucCoA", "Succinate", "tca", -0.3),
    e("Succinate", "Fumarate", "tca", -0.3),
    e("Fumarate", "Malate", "tca", -0.3),
    e("Malate", "OAA", "tca", -0.3),

    # Full PPP branch
    e("Glc-6-P", "6-PGL", "ppp", 0.15),
    e("6-PGL", "6-PG", "ppp"),
    e("6-PG", "Ru-5-P", "ppp"),
    e("Ru-5-P", "R-5-P", "ppp"),
    e("Ru-5-P", "Xu-5-P", "ppp", -0.15),
    e("Xu-5-P", "S-7-P", "ppp"),
    e("S-7-P", "E-4-P", "ppp"),

    # Entry intermediates → backbone
    e("Glc-1-P", "Glc-6-P", "starch"),
    e("Trehalose-6-P", "Glc-6-P", "starch", -0.1),  # treC / trehalose-6P hydrolase route
    e("Gal-1-P", "UDP-Gal", "galactose", 0.10), # GalT
    e("UDP-Gal", "UDP-Glc", "galactose", 0.10), # GalE
    e("UDP-Glc", "Glc-1-P", "galactose", 0.10), # GalT output
    e("GlcNAc-6-P", "GlcN-6-P", "chitin"),      # nagA: GlcNAc-6-P → GlcN-6-P
    e("GlcN-6-P", "Fru-6-P", "chitin"),          # nagB: GlcN-6-P → Fru-6-P
    e("Man-6-P", "Fru-6-P", "mannan", -0.1),     # manA: Man-6-P → Fru-6-P
    e("Fru-1-P", "Fru-1,6-BP", "fructose"),      # fruK: Fru-1-P → FBP
    e("Glycerol-3-P", "DHAP", "glycerol"),        # glpK/glpD
    e("Mannitol-1-P", "Fru-6-P", "mannan", -0.2), # mtlD: Mtl-1-P → F6P
    e("Xylulose", "Xu-5-P", "xylan", 0.10),      # xylA/xylB
    e("Ribulose", "Ribulose-5-P", "xylan", -0.10), # araA / araB
    e("Ribulose-5-P", "Xu-5-P", "xylan", -0.10),   # araD
    e("Fuculose", "Fuculose-1-P", "fucose", -0.10), # fucI / fucK
    e("Fuculose-1-P", "DHAP", "fucose", -0.15),  # fucI/fucK
    e("Rhamnulose", "Rhamnulose-1-P", "rhamnose"),
    e("Rhamnulose-1-P", "DHAP", "rhamnose", -0.15),  # rhaA/rhaB/rhaD

    # ED pathway: 6-PG → KDPG → Pyruvate + GAP
    e("Gluconate", "KDG", "gluconate", -0.1),     # unphosphorylated Entner-Doudoroff entry
    e("KDG", "KDPG", "gluconate", -0.1),          # kdgK
    e("6-PG", "KDPG", "ppp"),                     # edd
    e("KDPG", "Pyruvate", "ppp", 0.2),            # eda → Pyruvate
    e("KDPG", "GA3P", "ppp", -0.2),               # eda → GAP

    # Pyruvate branches (fermentation)
    e("Pyruvate", "Formate", "pyruvate_branch"),   # pfl
    e("Pyruvate", "Lactate", "pyruvate_branch"),   # ldh
    e("Acetyl-CoA", "Acetyl-P", "pyruvate_branch"), # pta
    e("Acetyl-P", "Acetate", "pyruvate_branch"),    # ackA
    e("Acetyl-CoA", "Acetaldehyde", "pyruvate_branch"), # aldh
    e("Acetaldehyde", "Ethanol", "pyruvate_branch"),    # adh

    # Anaplerotic reactions
    e("PEP", "OAA", "backbone", -0.2),              # ppc: PEP carboxylase (anaplerosis)

    # Glyoxylate shunt
    e("Isocit", "Glyoxylate", "tca", 0.2),        # aceA: isocitrate lyase
    e("Glyoxylate", "Malate", "tca", -0.2),        # aceB: malate synthase

    # Transketolase/Transaldolase cross-links (PPP ↔ Glycolysis)
    e("R-5-P", "S-7-P", "ppp", -0.15),            # tkt: R5P + Xu5P → S7P + GAP
    e("Xu-5-P", "GA3P", "ppp", 0.15),             # tkt: Xu5P + R5P → S7P + GAP
    e("S-7-P", "Fru-6-P", "ppp", 0.2),            # tal: S7P + GAP → F6P + E4P
    e("E-4-P", "Fru-6-P", "ppp", 0.15),           # tal: E4P + Xu5P → F6P + GAP

    # Carbon sources → entry intermediates/backbone
    # Straight lines by default; small curvature only when paths would overlap
    e("Maltose", "Glucose", "starch"),
    e("Cellobiose", "Glucose", "cellulose", -0.15),
    e("Galactose", "Gal-1-P", "galactose"),
    e("Trehalose", "Trehalose-6-P", "starch", -0.25),
    e("Lactose", "Gal-1-P", "galactose", -0.15),
    e("Lactose", "Glucose", "galactose", 0.15),     # β-gal: Lac → Gal + Glc
    e("Mannose", "Man-6-P", "mannan"),
    e("NAG", "GlcNAc-6-P", "chitin"),
    e("Glucosamine", "GlcNAc-6-P", "chitin", -0.15),
    e("Fructose", "Fru-1-P", "fructose"),
    e("Sucrose", "Fru-6-P", "sucrose", -0.15),
    e("Sucrose", "Glucose", "sucrose", 0.15),        # invertase: Suc → Glc + Fru
    e("Mannitol", "Mannitol-1-P", "mannan", -0.15),  # mtlK: Mannitol → Mtl-1-P
    e("Glycerol", "Glycerol-3-P", "glycerol"),
    e("Xylose", "Xylulose", "xylan"),
    e("Arabinose", "Ribulose", "xylan", -0.15),
    e("Ribose", "R-5-P", "pentose"),
    e("Fucose", "Fuculose", "fucose", -0.2),
    e("Rhamnose", "Rhamnulose", "rhamnose", -0.2),
    e("Gluconate", "6-PG", "gluconate", -0.15)
  )

  do.call(rbind, Filter(Negate(is.null), edges))
}

# Merge entry intermediates that are already represented on the main stream
# (e.g. Fru-6-P-entry -> backbone Fru-6-P, Xu-5-P-entry -> PPP Xu-5-P).
.dnmb_cct_merge_entry_into_mainstream <- function(nodes) {
  if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) == 0) return(nodes)

  main_nodes <- nodes[nodes$type %in% c("backbone", "ppp", "tca"), , drop = FALSE]
  entry_nodes <- nodes[nodes$type == "entry_intermediate", , drop = FALSE]
  if (nrow(main_nodes) == 0 || nrow(entry_nodes) == 0) return(nodes)

  main_labels <- unique(main_nodes$label[!is.na(main_nodes$label) & nzchar(main_nodes$label)])
  drop_idx <- nodes$type == "entry_intermediate" &
    !is.na(nodes$label) &
    nodes$label %in% main_labels

  nodes[!drop_idx, , drop = FALSE]
}

.dnmb_cct_snap_to_grid <- function(x, step = 0.5, origin = 0) {
  origin + round((x - origin) / step) * step
}

# Auto-resolve grid overlaps: ensure each 0.5 cell has at most one node
# Displaces conflicting nodes to nearest empty cell (prefer x+0.5, then x-0.5)
.dnmb_cct_resolve_grid_overlaps <- function(nodes, step = 0.5) {
  if (is.null(nodes) || nrow(nodes) < 2) return(nodes)
  # TCA nodes use exact circle positions — exempt from grid snapping
  is_tca <- nodes$type == "tca"
  tca_save_x <- nodes$x[is_tca]
  tca_save_y <- nodes$y[is_tca]
  nodes$x <- .dnmb_cct_snap_to_grid(nodes$x, step = step)
  nodes$y <- .dnmb_cct_snap_to_grid(nodes$y, step = step)
  # Restore TCA exact positions
  nodes$x[is_tca] <- tca_save_x
  nodes$y[is_tca] <- tca_save_y

  # Priority: backbone > ppp > tca > ed > entry_intermediate > carbon_source
  type_priority <- c(backbone = 6, ppp = 5, tca = 4, ed = 3,
                     entry_intermediate = 2, carbon_source = 1, pyruvate_branch = 3)
  nodes$grid_priority <- ifelse(nodes$type %in% names(type_priority),
                                type_priority[nodes$type], 0)
  nodes <- nodes[order(-nodes$grid_priority), , drop = FALSE]

  occupied <- list()  # key = "x,y" -> node index
  for (i in seq_len(nrow(nodes))) {
    key <- paste(nodes$x[i], nodes$y[i], sep = ",")
    if (is.null(occupied[[key]])) {
      occupied[[key]] <- i
    } else {
      # This cell is taken — find nearest empty cell
      found <- FALSE
      for (d in seq_len(8)) {
        for (dy in c(0, -step, step)) {
          for (dx in c(step, -step)) {
            nx <- .dnmb_cct_snap_to_grid(nodes$x[i] + dx * d, step = step)
            ny <- .dnmb_cct_snap_to_grid(nodes$y[i] + dy, step = step)
            nkey <- paste(nx, ny, sep = ",")
            if (is.null(occupied[[nkey]])) {
              nodes$x[i] <- nx; nodes$y[i] <- ny
              occupied[[nkey]] <- i; found <- TRUE; break
            }
          }
          if (found) break
        }
        if (found) break
      }
    }
  }
  nodes$grid_priority <- NULL
  nodes[order(as.numeric(rownames(nodes))), , drop = FALSE]
}

.dnmb_cct_compact_ordered_positions <- function(ids, x, center = NULL, spacing = NULL,
                                                shrink = 0.86, step = 0.25) {
  ids <- as.character(ids)
  x <- stats::setNames(as.numeric(x), ids)
  if (is.null(center) || !is.finite(center)) center <- stats::median(x, na.rm = TRUE)
  out <- center + (x - center) * shrink
  out <- .dnmb_cct_snap_to_grid(out, step = step)

  if (is.null(spacing)) {
    min_gap <- rep(0.4, max(0, length(ids) - 1L))
  } else {
    spacing <- stats::setNames(as.numeric(spacing), ids)
    min_gap <- numeric(max(0, length(ids) - 1L))
    if (length(ids) >= 2) {
      for (i in 2:length(ids)) {
        min_gap[i - 1L] <- max(step, 0.65 * (spacing[ids[i - 1L]] + spacing[ids[i]]) / 2)
      }
    }
  }

  if (length(ids) >= 2) {
    for (i in 2:length(ids)) {
      out[ids[i]] <- max(out[ids[i]], out[ids[i - 1L]] + min_gap[i - 1L])
      out[ids[i]] <- .dnmb_cct_snap_to_grid(out[ids[i]], step = step)
    }
    for (i in (length(ids) - 1L):1L) {
      out[ids[i]] <- min(out[ids[i]], out[ids[i + 1L]] - min_gap[i])
      out[ids[i]] <- .dnmb_cct_snap_to_grid(out[ids[i]], step = step)
    }
  }

  out[ids]
}

.dnmb_cct_cs_lane_group <- function(id) {
  switch(as.character(id)[1],
    Maltose = "glucose",
    Cellobiose = "glucose",
    Trehalose = "glucose",
    Galactose = "galactose",
    Lactose = "galactose",
    NAG = "glcnac",
    Glucosamine = "glcnac",
    Mannose = "mannose",
    Mannitol = "mannose",
    Fructose = "fructose",
    Sucrose = "fructose",
    Glycerol = "glycerol",
    Xylose = "pentose",
    Arabinose = "pentose",
    Ribose = "pentose",
    Fucose = "deoxy",
    Rhamnose = "deoxy",
    Gluconate = "gluconate",
    "other"
  )
}

.dnmb_cct_pack_cs_order <- function(order_ids, init_x, target_x, spacing, evidence,
                                    priority = NULL,
                                    grid_step = 0.5) {
  n <- length(order_ids)
  if (n == 0) {
    return(list(pos = stats::setNames(numeric(0), character(0)), cost = 0))
  }

  gap_vec <- numeric(max(0, n - 1L))
  if (n >= 2) {
    for (i in 2:n) {
      prev_id <- order_ids[i - 1L]
      cur_id <- order_ids[i]
      gap_vec[i - 1L] <- (spacing[prev_id] + spacing[cur_id]) / 2 +
        ifelse(evidence[prev_id] != evidence[cur_id], 0.25, 0)
    }
  }

  offsets <- c(0, if (length(gap_vec) > 0) cumsum(gap_vec) else numeric(0))
  ids_target <- unname(target_x[order_ids])
  ids_init <- unname(init_x[order_ids])
  ids_evidence <- unname(evidence[order_ids])
  if (is.null(priority)) {
    ids_priority <- rep(1, length(order_ids))
  } else {
    ids_priority <- unname(priority[order_ids])
    ids_priority[!is.finite(ids_priority)] <- 1
  }
  pr_norm <- ids_priority / max(1, max(ids_priority, na.rm = TRUE))
  target_w <- 1.25 + 0.36 * abs(ids_target - ids_init) +
    ifelse(ids_evidence == "cazy", 0.55, ifelse(ids_evidence == "transport", 0.25, ifelse(ids_evidence == "enzyme", 0.15, 0)))
  target_w <- target_w + 0.60 * pr_norm
  init_w <- 0.16 + 0.05 * (ids_evidence != "none") + 0.04 * pr_norm

  origin_min <- .dnmb_cct_snap_to_grid(min(c(ids_target, ids_init)) - 1.5, step = grid_step)
  origin_max <- .dnmb_cct_snap_to_grid(max(c(ids_target, ids_init)) + 1.5, step = grid_step)
  candidate_origins <- seq(origin_min, origin_max, by = grid_step)
  if (length(candidate_origins) == 0) candidate_origins <- origin_min
  center_ref <- stats::median(ids_target, na.rm = TRUE)

  best_cost <- Inf
  best_pos <- NULL
  for (origin in candidate_origins) {
    pos <- origin + offsets
    names(pos) <- order_ids
    target_cost <- sum(abs(pos - ids_target) * target_w)
    init_cost <- sum(abs(pos - ids_init) * init_w)
    elbow_cost <- sum(pmax(0, abs(pos - ids_target) - 2.0)) * 0.28
    span_cost <- (max(pos) - min(pos)) * 0.18
    center_penalty <- sum(pmax(0, 1.9 - abs(pos - center_ref)) * pmax(0, 1.20 - pr_norm)) * 0.55
    total_cost <- target_cost + init_cost + elbow_cost + span_cost + center_penalty
    if (total_cost + 1e-9 < best_cost) {
      best_cost <- total_cost
      best_pos <- pos
    }
  }

  list(pos = best_pos, cost = best_cost)
}

.dnmb_cct_optimize_cs_positions <- function(ids, init_x, target_x, spacing, evidence,
                                            priority = NULL,
                                            grid_step = 0.5, n_iter = 10L) {
  ids <- as.character(ids)
  if (length(ids) == 0) return(stats::setNames(numeric(0), character(0)))
  init_x <- stats::setNames(as.numeric(init_x), ids)
  target_x <- stats::setNames(as.numeric(target_x), ids)
  spacing <- stats::setNames(as.numeric(spacing), ids)
  evidence <- stats::setNames(as.character(evidence), ids)
  if (is.null(priority)) priority <- stats::setNames(rep(1, length(ids)), ids)
  priority <- stats::setNames(as.numeric(priority), ids)
  lane_group <- stats::setNames(vapply(ids, .dnmb_cct_cs_lane_group, character(1)), ids)
  group_levels <- unique(unname(lane_group[ids]))
  group_members <- lapply(group_levels, function(grp) {
    mids <- ids[lane_group[ids] == grp]
    mids <- mids[order(target_x[mids], init_x[mids])]
    if (length(mids) <= 1L) return(mids)
    perms <- lapply(.dnmb_cct_small_permutations(seq_along(mids)), function(idx) mids[idx])
    best_perm <- mids
    best_cost <- Inf
    for (perm in perms) {
      target_rank <- seq_along(perm)
      perm_cost <- sum(abs(target_x[perm] - sort(target_x[mids]))) * 0.85 +
        sum(abs(init_x[perm] - sort(init_x[mids]))) * 0.15 +
        sum(abs(match(perm, mids) - target_rank)) * 0.08
      if (perm_cost + 1e-9 < best_cost) {
        best_cost <- perm_cost
        best_perm <- perm
      }
    }
    best_perm
  })
  names(group_members) <- group_levels

  eval_group_order <- function(group_order) {
    order_ids <- unlist(group_members[group_order], use.names = FALSE)
    packed <- .dnmb_cct_pack_cs_order(
      order_ids = order_ids,
      init_x = init_x,
      target_x = target_x,
      spacing = spacing,
      evidence = evidence,
      priority = priority,
      grid_step = grid_step
    )
    order_rank <- match(order_ids, ids)
    group_shift <- match(group_order, group_levels)
    group_centers <- vapply(group_order, function(grp) {
      mids <- group_members[[grp]]
      mean(unname(packed$pos[mids]))
    }, numeric(1))
    group_targets <- vapply(group_order, function(grp) {
      mids <- group_members[[grp]]
      stats::median(unname(target_x[mids]))
    }, numeric(1))
    packed$cost <- packed$cost +
      sum(abs(order_rank - seq_along(order_ids))) * 0.14 +
      sum(abs(group_shift - seq_along(group_order))) * 0.42 +
      sum(abs(group_centers - group_targets)) * 0.52
    packed$group_order <- group_order
    packed
  }

  best <- eval_group_order(group_levels)
  improved <- TRUE
  iter <- 0L
  while (improved && iter < max(3L, as.integer(n_iter))) {
    improved <- FALSE
    iter <- iter + 1L
    cur_order <- best$group_order
    cand_best <- best
    if (length(cur_order) >= 2) {
      for (i in seq_len(length(cur_order) - 1L)) {
        for (j in (i + 1L):length(cur_order)) {
          swap_order <- cur_order
          swap_order[c(i, j)] <- swap_order[c(j, i)]
          cand <- eval_group_order(swap_order)
          if (cand$cost + 1e-9 < cand_best$cost) cand_best <- cand
        }
      }
      for (i in seq_len(length(cur_order))) {
        for (j in seq_len(length(cur_order))) {
          if (i == j) next
          move_order <- cur_order[-i]
          insert_at <- if (j <= 1L) 1L else if (j >= length(cur_order)) length(move_order) + 1L else j
          move_order <- append(move_order, cur_order[i], after = insert_at - 1L)
          cand <- eval_group_order(move_order)
          if (cand$cost + 1e-9 < cand_best$cost) cand_best <- cand
        }
      }
    }
    if (cand_best$cost + 1e-9 < best$cost) {
      best <- cand_best
      improved <- TRUE
    }
  }

  out <- init_x
  out[names(best$pos)] <- best$pos
  out[ids]
}

.dnmb_cct_quarter_arc_points <- function(corner, dir_in, dir_out, radius, n = 12) {
  center <- corner - dir_in * radius + dir_out * radius
  start  <- corner - dir_in * radius
  end    <- corner + dir_out * radius
  a0 <- atan2(start[2] - center[2], start[1] - center[1])
  a1 <- atan2(end[2] - center[2], end[1] - center[1])
  turn <- dir_in[1] * dir_out[2] - dir_in[2] * dir_out[1]

  if (turn >= 0) {
    if (a1 < a0) a1 <- a1 + 2 * pi
    ang <- seq(a0, a1, length.out = n)
  } else {
    if (a1 > a0) a1 <- a1 - 2 * pi
    ang <- seq(a0, a1, length.out = n)
  }

  data.frame(
    x = center[1] + radius * cos(ang),
    y = center[2] + radius * sin(ang)
  )
}

.dnmb_cct_rounded_route_points <- function(points_df, radius = 0.5, n_arc = 12) {
  if (is.null(points_df) || nrow(points_df) < 3) return(points_df)

  # Handle 3-point L-shape (single bend): start → corner → end
  if (nrow(points_df) == 3) {
    pts <- as.matrix(points_df[, c("x", "y"), drop = FALSE])
    d01 <- pts[2, ] - pts[1, ]; l01 <- sum(abs(d01))
    d12 <- pts[3, ] - pts[2, ]; l12 <- sum(abs(d12))
    if (l01 <= 0 || l12 <= 0) return(points_df)
    u01 <- d01 / l01; u12 <- d12 / l12
    r_eff <- min(radius, l01 * 0.40, l12 * 0.40)
    if (!is.finite(r_eff) || r_eff <= 0.02) return(points_df)
    a1 <- pts[2, ] - u01 * r_eff  # straight segment ends here (before arc)
    b1 <- pts[2, ] + u12 * r_eff  # arc ends here (straight segment resumes)
    arc1 <- .dnmb_cct_quarter_arc_points(pts[2, ], u01, u12, r_eff, n = n_arc)
    path <- rbind(
      data.frame(x = pts[1, 1], y = pts[1, 2]),  # start
      data.frame(x = a1[1], y = a1[2]),            # end of first straight
      arc1[-1, , drop = FALSE],                    # curved bend
      data.frame(x = b1[1], y = b1[2]),            # start of second straight
      data.frame(x = pts[3, 1], y = pts[3, 2])    # end (destination)
    )
    # Remove duplicate consecutive points
    keep <- c(TRUE, diff(path$x) != 0 | diff(path$y) != 0)
    return(path[keep, , drop = FALSE])
  }

  pts <- as.matrix(points_df[, c("x", "y"), drop = FALSE])
  d01 <- pts[2, ] - pts[1, ]
  d12 <- pts[3, ] - pts[2, ]
  d23 <- pts[4, ] - pts[3, ]
  l01 <- sum(abs(d01))
  l12 <- sum(abs(d12))
  l23 <- sum(abs(d23))
  if (l01 <= 0 || l12 <= 0 || l23 <= 0) return(points_df)

  u01 <- d01 / l01
  u12 <- d12 / l12
  u23 <- d23 / l23

  r_eff <- min(radius, l01, l23, l12 / 2)
  if (!is.finite(r_eff) || r_eff <= 0.02) return(points_df)

  a1 <- pts[2, ] - u01 * r_eff
  b1 <- pts[2, ] + u12 * r_eff
  a2 <- pts[3, ] - u12 * r_eff
  b2 <- pts[3, ] + u23 * r_eff

  arc1 <- .dnmb_cct_quarter_arc_points(pts[2, ], u01, u12, r_eff, n = n_arc)
  arc2 <- .dnmb_cct_quarter_arc_points(pts[3, ], u12, u23, r_eff, n = n_arc)

  path <- rbind(
    data.frame(x = pts[1, 1], y = pts[1, 2]),
    data.frame(x = a1[1], y = a1[2]),
    arc1[-1, , drop = FALSE],
    data.frame(x = a2[1], y = a2[2]),
    arc2[-1, , drop = FALSE],
    data.frame(x = pts[4, 1], y = pts[4, 2])
  )

  keep <- c(TRUE, diff(path$x) != 0 | diff(path$y) != 0)
  path[keep, , drop = FALSE]
}

# Per-substrate pathway colors (gapmind_aa style)
.dnmb_cct_pathway_colors <- function() {
  c(backbone  = "#444444",
    ppp       = "#2B8CBE",
    tca       = "#E7298A",    # magenta
    starch    = "#F47920",    # orange (Starch/Maltose)
    cellulose = "#00A651",    # green (Cellulose/Cellobiose)
    galactose = "#FFD400",    # yellow (Galactose/Lactose)
    chitin    = "#0072BC",    # blue (Chitin/NAG)
    mannan    = "#00A651",    # green (Mannose)
    fructose  = "#A8D08D",    # light green (Fructose/Sucrose)
    sucrose   = "#A8D08D",    # light green
    xylan     = "#F69EA1",    # pink (Xylose/Arabinose)
    pentose   = "#F69EA1",    # pink
    fucose    = "#ED1C24",    # red
    gluconate = "#6BAED6",    # light blue
    glycerol  = "#3182BD",    # medium blue
    ed        = "#9E5BA3",    # purple (Entner-Doudoroff)
    pyruvate_branch = "#888888",  # grey (fermentation products)
    none      = "#CCCCCC")
}

# ---------------------------------------------------------------------------
# 3-Zone CAZy Carbon Transport Map — main plot function
# ---------------------------------------------------------------------------

#' Plot the full 3-zone CAZy Carbon Transport Map
#'
#' Zone 1 (Extracellular): Glycan chains drawn as SNFG symbols,
#'   GH cleavage scissors with locus_tag + GH_family, degradation cascades
#' Zone 2 (Membrane): Transporters with arrows, PTS cascade
#' Zone 3 (Cytoplasm): Metabolism backbone (glycolysis + PPP + ED) with
#'   enzyme labels + confidence coloring
# PUL/CGC substrate name → carbon source ID mapping
.dnmb_cct_pul_substrate_to_cs <- function() {
  c(
    starch = "Maltose", cellulose = "Cellobiose", xylan = "Xylose",
    chitin = "NAG", pectin = "Gluconate", mannan = "Mannose",
    inulin = "Fructose", fructan = "Fructose", sucrose = "Sucrose",
    "beta-glucan" = "Cellobiose", "beta-mannan" = "Mannose",
    arabinan = "Arabinose", galactan = "Galactose",
    trehalose = "Trehalose", lactose = "Lactose", raffinose = "Sucrose",
    fucose = "Fucose", rhamnose = "Rhamnose", maltose = "Maltose",
    glucan = "Maltose", galactomannan = "Mannose",
    xyloglucan = "Xylose", arabinoxylan = "Xylose",
    laminarin = "Cellobiose", lichenan = "Cellobiose",
    "beta-glucan" = "Cellobiose", betaglucan = "Cellobiose",
    agarose = "Galactose",
    agar = "Galactose", carrageenan = "Galactose",
    cellobiose = "Cellobiose", mannose = "Mannose",
    fructose = "Fructose", glucose = "Maltose",
    galactose = "Galactose", ribose = "Ribose",
    glycerol = "Glycerol", gluconate = "Gluconate"
  )
}

# Aggregate CGC/PUL substrate evidence from genbank_table
.dnmb_cct_aggregate_cgc_pul_evidence <- function(genbank_table) {
  pul_map <- .dnmb_cct_pul_substrate_to_cs()
  evidence_cs <- character(0)

  # Source 1: PUL substrate predictions
  pul_col <- intersect(c("dbCAN_dbcan_pul_substrate", "dbcan_pul_substrate"), names(genbank_table))
  if (length(pul_col) > 0) {
    pul_vals <- unique(as.character(genbank_table[[pul_col[1]]]))
    pul_vals <- pul_vals[!is.na(pul_vals) & nzchar(pul_vals)]
    pul_vals <- unlist(strsplit(pul_vals, ";\\s*"))
    pul_vals <- trimws(tolower(pul_vals))
    pul_vals <- pul_vals[nzchar(pul_vals)]
    mapped <- unname(pul_map[pul_vals])
    evidence_cs <- c(evidence_cs, mapped[!is.na(mapped)])
  }

  # Source 2: dbCAN-sub substrate predictions
  sub_col <- intersect(c("dbCAN_dbcan_sub_substrate", "dbcan_sub_substrate"), names(genbank_table))
  if (length(sub_col) > 0) {
    sub_vals <- unique(as.character(genbank_table[[sub_col[1]]]))
    sub_vals <- sub_vals[!is.na(sub_vals) & nzchar(sub_vals)]
    sub_vals <- unlist(strsplit(sub_vals, ";\\s*"))
    sub_vals <- trimws(tolower(sub_vals))
    sub_vals <- sub_vals[nzchar(sub_vals)]
    mapped <- unname(pul_map[sub_vals])
    evidence_cs <- c(evidence_cs, mapped[!is.na(mapped)])
  }

  unique(evidence_cs)
}

# GH family → carbon source mapping for direct evidence scoring
# Maps each GH family prefix to the carbon sources it acts on
.dnmb_cct_gh_family_to_cs <- function() {
  list(
    # Starch / Maltose
    GH13 = c("Maltose", "Trehalose"), GH14 = "Maltose", GH15 = "Maltose",
    GH31 = "Maltose", GH57 = "Maltose", GH77 = "Maltose",
    GH97 = "Maltose",
    # Cellulose / Cellobiose
    GH5 = "Cellobiose", GH6 = "Cellobiose", GH7 = "Cellobiose",
    GH9 = "Cellobiose", GH44 = "Cellobiose", GH45 = "Cellobiose",
    GH48 = "Cellobiose", GH12 = "Cellobiose",
    GH1 = c("Cellobiose", "Lactose"), GH3 = c("Cellobiose", "Xylose"),
    GH4 = "Cellobiose",
    # Xylan / Arabinose
    GH10 = "Xylose", GH11 = "Xylose", GH30 = "Xylose",
    GH8 = "Xylose", GH67 = "Xylose", GH115 = "Xylose",
    GH43 = c("Xylose", "Arabinose"), GH51 = "Arabinose", GH54 = "Arabinose",
    GH62 = "Arabinose",
    # Chitin / NAG
    GH18 = "NAG", GH19 = "NAG", GH20 = "NAG", GH23 = "NAG",
    GH84 = "NAG", GH85 = "NAG",
    # Pectin / Gluconate
    GH28 = "Gluconate", GH78 = c("Gluconate", "Rhamnose"), GH105 = "Gluconate",
    # Mannan / Mannose
    GH26 = "Mannose", GH113 = "Mannose", GH130 = "Mannose",
    GH134 = "Mannose", GH76 = "Mannose",
    # Fructan / Fructose / Sucrose
    GH32 = c("Fructose", "Sucrose"), GH68 = c("Fructose", "Sucrose"),
    GH91 = "Sucrose",
    # Galactose / Lactose
    GH2 = c("Galactose", "Lactose"), GH35 = c("Galactose", "Lactose"),
    GH42 = c("Galactose", "Lactose"),
    GH27 = "Galactose", GH36 = "Galactose", GH110 = "Galactose",
    # Trehalose
    GH37 = "Trehalose", GH65 = "Trehalose",
    # Fucose
    GH29 = "Fucose", GH95 = "Fucose", GH141 = "Fucose",
    # Rhamnose
    GH106 = "Rhamnose", GH145 = "Rhamnose",
    # Inulin (→ Fructose)
    GH91 = "Fructose"
  )
}

# Metabolic enzyme evidence: scan gene names and product descriptions
# for known pathway-specific enzymes (cytoplasmic, not transport)
.dnmb_cct_metabolic_enzyme_evidence <- function(genbank_table) {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  gene_col <- .dnmb_pick_column(tbl, c("gene", "gene_name"))
  prod_col <- .dnmb_pick_column(tbl, c("product"))
  ec_col   <- .dnmb_pick_column(tbl, c("EC_number", "ec_number"))

  evidence_cs <- character(0)

  # Gene name patterns → carbon sources (cytoplasmic enzymes only)
  gene_to_cs <- list(
    Maltose     = "^(malQ|malZ|amyA|malP|pgm|treC|glgP)",
    Galactose   = "^(galK|galT|galE|galM|pgmA|galU)",
    Lactose     = "^(lacZ|lacA|bgaB|bga|ebgA)",
    Mannose     = "^(manA|manB|manC|pmi|algC)",
    NAG         = "^(nagA|nagB|nagC|nagK|gamA|nag1|nag2)",
    Glucosamine = "^(nagB|gamA|glmS|glmM)",
    Fructose    = "^(fruK|fruA|fruB|scrK|1pfk)",
    Sucrose     = "^(sacA|scrB|sacB|sacC|inv[ABCD]|sucA)",
    Mannitol    = "^(mtlD|mtlK|mtlA)",
    Glycerol    = "^(glpK|glpD|glpA|glpB|glpC|dhaK|dhaL|dhaM)",
    Xylose      = "^(xylA|xylB|xylR)",
    Arabinose   = "^(araA|araB|araD|araR)",
    Ribose      = "^(rbsK|rbsR|rbsD)",
    Fucose      = "^(fucI|fucK|fucA|fucO|fucR)",
    Rhamnose    = "^(rhaA|rhaB|rhaD|rhaR|rhaS)",
    Gluconate   = "^(gntK|gntR|gntU|gntT|gnd|edd|eda|kdgK|kduD)",
    Trehalose   = "^(treA|treC|treR|treS|otsA|otsB)"
  )
  if (!is.null(gene_col)) {
    genes <- tolower(as.character(tbl[[gene_col]]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    for (csid in names(gene_to_cs)) {
      if (any(grepl(gene_to_cs[[csid]], genes, ignore.case = TRUE)))
        evidence_cs <- c(evidence_cs, csid)
    }
  }

  # Product description patterns (fallback when gene names are absent)
  prod_to_cs <- list(
    Maltose     = "alpha-glucosidase|maltase|amylase|pullulanase|maltodextrin",
    Galactose   = "galactokinase|galactose-1-phosphate|UDP-glucose.*epimerase",
    Lactose     = "beta-galactosidase|lactase",
    Mannose     = "mannose-6-phosphate isomerase|phosphomannomutase",
    NAG         = "N-acetylglucosamine|chitobiase|nagA|nagB",
    Fructose    = "1-phosphofructokinase|fructokinase",
    Sucrose     = "invertase|sucrase|levansucrase|sucrose phosphorylase",
    Mannitol    = "mannitol.*(dehydrogenase|kinase)|mannitol-1-phosphate",
    Glycerol    = "glycerol kinase|glycerol-3-phosphate dehydrogenase|glycerol dehydrogenase",
    Xylose      = "xylose isomerase|xylulokinase",
    Arabinose   = "L-arabinose isomerase|ribulokinase|L-ribulose",
    Ribose      = "ribokinase|ribose-5-phosphate",
    Fucose      = "L-fucose isomerase|L-fuculose|fuculokinase",
    Rhamnose    = "L-rhamnose isomerase|rhamnulokinase|L-rhamnulose",
    Gluconate   = "gluconate kinase|gluconokinase|Entner-Doudoroff|KDG|KDPG",
    Trehalose   = "trehalase|trehalose-6-phosphate|trehalose synthase"
  )
  if (!is.null(prod_col)) {
    products <- tolower(as.character(tbl[[prod_col]]))
    products <- products[!is.na(products) & nzchar(products)]
    for (csid in names(prod_to_cs)) {
      if (csid %in% evidence_cs) next  # already found
      if (any(grepl(prod_to_cs[[csid]], products, ignore.case = TRUE)))
        evidence_cs <- c(evidence_cs, csid)
    }
  }

  # EC number patterns (most specific, catches annotations lacking gene names)
  ec_to_cs <- list(
    Maltose     = "^3\\.2\\.1\\.20$|^3\\.2\\.1\\.1$|^2\\.4\\.1\\.1$",   # alpha-glucosidase, alpha-amylase, glycogen phosphorylase
    Galactose   = "^2\\.7\\.1\\.6$|^2\\.7\\.7\\.12$|^5\\.1\\.3\\.2$",   # galK, galT, galE
    Lactose     = "^3\\.2\\.1\\.23$|^3\\.2\\.1\\.108$",                  # beta-galactosidase, lactase
    NAG         = "^3\\.5\\.1\\.25$|^3\\.5\\.99\\.6$",                   # nagA, chitobiase
    Fructose    = "^2\\.7\\.1\\.56$",                                     # 1-phosphofructokinase
    Sucrose     = "^3\\.2\\.1\\.26$|^2\\.4\\.1\\.7$",                    # invertase, sucrose phosphorylase
    Xylose      = "^5\\.3\\.1\\.5$",                                      # xylose isomerase
    Arabinose   = "^5\\.3\\.1\\.4$|^2\\.7\\.1\\.16$",                    # araA, araB
    Fucose      = "^5\\.3\\.1\\.25$|^2\\.7\\.1\\.51$",                   # fucI, fucK
    Rhamnose    = "^5\\.3\\.1\\.14$|^2\\.7\\.1\\.5$",                    # rhaA, rhaB
    Gluconate   = "^2\\.7\\.1\\.12$|^4\\.2\\.1\\.12$|^4\\.1\\.2\\.14$"  # gntK, edd, eda
  )
  if (!is.null(ec_col)) {
    ecs <- as.character(tbl[[ec_col]])
    ecs <- ecs[!is.na(ecs) & nzchar(ecs)]
    for (csid in names(ec_to_cs)) {
      if (csid %in% evidence_cs) next
      if (any(grepl(ec_to_cs[[csid]], ecs)))
        evidence_cs <- c(evidence_cs, csid)
    }
  }

  unique(evidence_cs)
}

#'
#' @param genbank_table The combined genbank/annotation table
#' @param output_dir Output directory
#' @param file_stub PDF filename stub
#' @return List with pdf path
#' @keywords internal
.dnmb_plot_cazy_carbon_3zone_v2 <- function(genbank_table, output_dir,
                                          file_stub = "CAZy_overview") {

  # ====================================================================
  # DATA EXTRACTION
  # ====================================================================
  step_status <- .dnmb_gapmind_carbon_step_status(genbank_table, output_dir = output_dir)
  if (is.null(step_status)) return(NULL)
  valid_locus_tags <- unique(as.character(genbank_table$locus_tag))
  valid_locus_tags <- valid_locus_tags[!is.na(valid_locus_tags) & nzchar(valid_locus_tags)]
  if ("locus_tag" %in% names(step_status)) {
    step_status$locus_tag <- .dnmb_cct_filter_valid_locus_tags(step_status$locus_tag, valid_locus_tags)
  }

  # Define valid carbon source pathway IDs (exclude amino acid, organic acid, etc.)
  valid_carbon_pathways <- tolower(c(
    "glucose", "maltose", "cellobiose", "galactose", "trehalose", "lactose",
    "mannose", "nag", "glucosamine", "fructose", "sucrose", "mannitol",
    "glycerol", "xylose", "arabinose", "ribose", "fucose", "rhamnose",
    "gluconate", "deoxyribose", "deoxyribonate"
  ))
  pstats       <- .dnmb_gapmind_carbon_pathway_stats(step_status, output_dir = output_dir)
  # Filter pstats to carbon source pathways only
  if (!is.null(pstats) && nrow(pstats) > 0) {
    pstats <- pstats[tolower(pstats$pathway_id) %in% valid_carbon_pathways, , drop = FALSE]
  }
  gh_enzymes   <- .dnmb_cct_3zone_extract_gh(genbank_table)
  # Enrich GH enzymes with CGC/PUL annotation from genbank_table
  if (!is.null(gh_enzymes) && nrow(gh_enzymes) > 0) {
    cgc_col <- intersect(c("dbCAN_dbcan_cgc_id", "dbcan_cgc_id"), names(genbank_table))
    pul_col <- intersect(c("dbCAN_dbcan_pul_substrate", "dbcan_pul_substrate"), names(genbank_table))
    sub_col <- intersect(c("dbCAN_dbcan_sub_substrate", "dbcan_sub_substrate"), names(genbank_table))
    if (length(cgc_col) > 0 || length(pul_col) > 0 || length(sub_col) > 0) {
      enrich_cols <- c("locus_tag", cgc_col, pul_col, sub_col)
      enrich_df <- genbank_table[, intersect(enrich_cols, names(genbank_table)), drop = FALSE]
      gh_enzymes <- merge(gh_enzymes, enrich_df, by = "locus_tag", all.x = TRUE)
    }
  }
  transporters <- .dnmb_cct_3zone_extract_transporters(genbank_table, output_dir = output_dir)
  if (!is.null(transporters) && nrow(transporters) > 0 && "locus_tag" %in% names(transporters)) {
    transporters$locus_tag <- .dnmb_cct_filter_valid_locus_tags(transporters$locus_tag, valid_locus_tags)
  }
  cascades     <- .dnmb_cct_3zone_glycan_cascades()
  conf_colors  <- c(high = "#2CA25F", medium = "#FEC44F", low = "#F03B20", none = "#CCCCCC")
  conf_rank    <- c(none = 0L, low = 1L, medium = 2L, high = 3L)
  pathway_to_cs <- c(
    maltose = "Maltose", cellobiose = "Cellobiose", galactose = "Galactose",
    trehalose = "Trehalose", lactose = "Lactose", mannose = "Mannose",
    NAG = "NAG", glucosamine = "Glucosamine", fructose = "Fructose",
    sucrose = "Sucrose", mannitol = "Mannitol", glycerol = "Glycerol",
    xylose = "Xylose", arabinose = "Arabinose", ribose = "Ribose",
    fucose = "Fucose", rhamnose = "Rhamnose", gluconate = "Gluconate",
    glucose = "Maltose", deoxyribose = "Ribose", deoxyribonate = "Gluconate"
  )

  # Filter step_status to ONLY carbon source pathways (exclude amino acid, organic acid, etc.)
  matched_steps <- step_status[
    !is.na(step_status$locus_tag) & nzchar(step_status$locus_tag) &
    tolower(step_status$pathway_id) %in% valid_carbon_pathways,
    , drop = FALSE
  ]
  if (nrow(matched_steps) > 0) {
    matched_steps$rank <- conf_rank[matched_steps$confidence]
    matched_steps$rank[is.na(matched_steps$rank)] <- 0L
  }
  matched_pathway_ids <- unique(tolower(matched_steps$pathway_id))

  route_group_members <- list(
    starch = c("maltose", "trehalose"),
    cellulose = c("cellobiose"),
    galactose = c("galactose", "lactose", "raffinose"),
    chitin = c("nag", "glucosamine"),
    mannan = c("mannose", "mannitol"),
    fructose = c("fructose", "mannitol"),
    sucrose = c("sucrose"),
    xylan = c("xylose", "arabinose"),
    pentose = c("ribose"),
    gluconate = c("gluconate"),
    glycerol = c("glycerol"),
    fucose = c("fucose"),
    rhamnose = c("rhamnose"),
    ppp = c("gluconate", "xylose", "arabinose", "ribose")
  )

  route_best_match <- function(path_ids) {
    if (length(path_ids) == 0 || nrow(matched_steps) == 0) return(NULL)
    idx <- tolower(matched_steps$pathway_id) %in% tolower(path_ids)
    cand <- matched_steps[idx, , drop = FALSE]
    if (nrow(cand) == 0) return(NULL)
    cand <- cand[order(-cand$rank, cand$pathway_id, cand$step_id), , drop = FALSE]
    cand[1, , drop = FALSE]
  }

  route_color_for <- function(path_ids, sugar_type, fallback = "#CFCFCF") {
    hit <- route_best_match(path_ids)
    if (is.null(hit)) return(fallback)
    .dnmb_cct_sugar_route_color(sugar_type, fallback = fallback)
  }

  # Match cascades to genome GH enzymes (prefix matching for subfamilies)
  active_cascades <- list()
  for (cname in names(cascades)) {
    casc <- cascades[[cname]]
    matched <- .dnmb_cct_3zone_match_gh(gh_enzymes, casc$gh_prefixes)
    if (!is.null(matched) && nrow(matched) > 0) {
      casc$matched_enzymes <- matched
      active_cascades[[cname]] <- casc
    }
  }

  # ====================================================================
  # EVIDENCE-BASED SUBSTRATE SPACING
  # ====================================================================
  # Carbon source IDs in canonical order (matches sugar_nodes index)
  cs_ids_ordered <- c("Maltose","Cellobiose","Galactose","Trehalose","Lactose",
                      "Mannose","NAG","Glucosamine","Fructose","Sucrose",
                      "Mannitol","Glycerol","Xylose","Arabinose","Ribose",
                      "Fucose","Rhamnose","Gluconate")
  # Classify each carbon source by evidence level
  cs_evidence <- stats::setNames(rep("none", length(cs_ids_ordered)), cs_ids_ordered)
  # Check pathway stats for transport/gapmind evidence (require >= 25% steps matched)
  if (!is.null(pstats) && nrow(pstats) > 0) {
    for (csid in cs_ids_ordered) {
      pid <- tolower(csid)
      if (pid %in% pstats$pathway_id && any(pstats$fraction[pstats$pathway_id == pid] >= 0.25))
        cs_evidence[csid] <- "transport"
    }
  }
  # Check for CAZy enzyme evidence (active cascades target substrates)
  # Map cascade keys to carbon source IDs
  cascade_to_cs <- c(starch = "Maltose", cellulose = "Cellobiose",
                     xylan = "Xylose", chitin = "NAG",
                     pectin = "Gluconate", mannan = "Mannose",
                     inulin = "Fructose", raffinose = "Sucrose",
                     sucrose = "Sucrose", maltose = "Maltose",
                     lactose = "Lactose", cellobiose = "Cellobiose",
                     trehalose = "Trehalose",
                     betaglucan = "Cellobiose",
                     arabinoxylan = "Xylose", agarose = "Galactose")
  for (cname in names(active_cascades)) {
    cs_target <- cascade_to_cs[cname]
    if (!is.na(cs_target) && cs_target %in% names(cs_evidence))
      cs_evidence[cs_target] <- "cazy"
  }
  # GH family → carbon source direct mapping (replaces broken substrate check)
  if (!is.null(gh_enzymes) && nrow(gh_enzymes) > 0) {
    gh_to_cs <- .dnmb_cct_gh_family_to_cs()
    present_gh <- unique(sub("_.*", "", gh_enzymes$gh_family))  # GH13_20 → GH13
    for (gh in present_gh) {
      cs_targets <- gh_to_cs[[gh]]
      if (!is.null(cs_targets)) {
        for (csid in cs_targets) {
          if (csid %in% names(cs_evidence) && cs_evidence[csid] == "none")
            cs_evidence[csid] <- "cazy"
        }
      }
    }
  }

  # Integrate CGC/PUL substrate predictions as additional evidence
  cgc_pul_cs <- .dnmb_cct_aggregate_cgc_pul_evidence(genbank_table)
  for (csid in cgc_pul_cs) {
    if (csid %in% names(cs_evidence) && cs_evidence[csid] == "none")
      cs_evidence[csid] <- "cazy"
  }

  # Metabolic enzyme evidence: check product/gene annotations for pathway enzymes
  metab_cs <- .dnmb_cct_metabolic_enzyme_evidence(genbank_table)
  for (csid in metab_cs) {
    if (csid %in% names(cs_evidence) && cs_evidence[csid] == "none")
      cs_evidence[csid] <- "enzyme"
  }

  # FILTER: keep only substrates with any evidence
  cs_has_evidence <- cs_evidence != "none"
  if (!any(cs_has_evidence)) {
    # No carbon utilization evidence — return placeholder
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
        label = "No CAZy / transport / CGC-PUL carbon utilization evidence detected",
        size = 5, color = "#888888") +
      ggplot2::theme_void()
    plot_dir <- .dnmb_module_plot_dir(output_dir)
    pdf_path <- file.path(plot_dir, paste0(file_stub, ".pdf"))
    .dnmb_module_plot_save(p, pdf_path, width = 10, height = 6)
    return(list(pdf = pdf_path))
  }
  cs_ids_ordered <- cs_ids_ordered[cs_has_evidence]
  cs_evidence <- cs_evidence[cs_has_evidence]

  cs_priority <- stats::setNames(rep(1, length(cs_ids_ordered)), cs_ids_ordered)
  if (!is.null(pstats) && nrow(pstats) > 0) {
    for (csid in cs_ids_ordered) {
      pid <- tolower(csid)
      if (pid %in% pstats$pathway_id) {
        frac <- pstats$fraction[pstats$pathway_id == pid][1]
        if (!is.na(frac)) cs_priority[csid] <- cs_priority[csid] + 1.10 * frac
      }
    }
  }
  for (cname in names(active_cascades)) {
    cs_target <- cascade_to_cs[cname]
    if (!is.na(cs_target) && cs_target %in% names(cs_priority)) {
      cs_priority[cs_target] <- cs_priority[cs_target] + 0.70
    }
  }
  if (!is.null(transporters) && nrow(transporters) > 0) {
    tr_priority <- transporters
    tr_priority$cs_id <- unname(pathway_to_cs[tr_priority$pathway])
    tr_priority <- tr_priority[!is.na(tr_priority$cs_id) & tr_priority$cs_id %in% cs_ids_ordered, , drop = FALSE]
    if (nrow(tr_priority) > 0) {
      if (!"step_score" %in% names(tr_priority)) tr_priority$step_score <- NA_real_
      tr_priority <- .dnmb_cct_annotate_transport_context(tr_priority, genbank_table)
      tr_priority <- .dnmb_cct_select_transporters(tr_priority, max_per_lane = 4L)
      if (nrow(tr_priority) > 0 && "display_score" %in% names(tr_priority)) {
        tr_score <- stats::aggregate(display_score ~ cs_id, tr_priority, max)
        cs_priority[tr_score$cs_id] <- cs_priority[tr_score$cs_id] + pmin(1.35, 0.10 * tr_score$display_score)
      }
    }
  }

  # Compute cumulative X positions on the glycolysis grid (0.5-unit lattice)
  glyco_grid_step <- 0.5
  # cazy = 1.5 (3 grid cells), transport = 1.0 (2 cells), none = 0.5 (1 cell)
  # group gap = 0.5 (1 cell) at evidence-level transitions
  cs_spacing <- vapply(cs_evidence, function(ev) {
    switch(ev, cazy = 1.5, transport = 1.0, enzyme = 0.8, none = 0.5)
  }, numeric(1))
  cs_x_pos <- stats::setNames(numeric(length(cs_ids_ordered)), cs_ids_ordered)
  cs_x_pos[1] <- .dnmb_cct_snap_to_grid(1.0, step = glyco_grid_step)
  for (i in seq_along(cs_ids_ordered)[-1]) {
    gap <- if (cs_evidence[i] != cs_evidence[i - 1]) glyco_grid_step else 0.0
    next_x <- cs_x_pos[i - 1] + (cs_spacing[i - 1] + cs_spacing[i]) / 2 + gap
    cs_x_pos[i] <- .dnmb_cct_snap_to_grid(next_x, step = glyco_grid_step)
  }

  # Compress carbon-source lanes toward their downstream entry targets.
  provisional_nodes <- .dnmb_cct_sugar_nodes(cs_xs = cs_x_pos, tca_cx_override = NULL)
  cs_entry_map <- .dnmb_cct_auto_entry_map()
  prov_target_nodes <- provisional_nodes[provisional_nodes$type %in% c("backbone", "ppp", "entry_intermediate", "ed"), , drop = FALSE]
  prov_target_x <- stats::setNames(prov_target_nodes$x, prov_target_nodes$id)
  cs_target_x <- vapply(cs_ids_ordered, function(id) {
    tgt <- unname(cs_entry_map[tolower(id)])
    if (!is.na(tgt) && tgt %in% names(prov_target_x)) {
      unname(prov_target_x[tgt])
    } else {
      cs_x_pos[match(id, cs_ids_ordered)]
    }
  }, numeric(1))
  cs_x_pos <- .dnmb_cct_optimize_cs_positions(
    ids = cs_ids_ordered,
    init_x = cs_x_pos,
    target_x = cs_target_x,
    spacing = cs_spacing,
    evidence = cs_evidence,
    priority = cs_priority,
    grid_step = glyco_grid_step
  )
  cs_x_pos <- .dnmb_cct_compact_ordered_positions(
    ids = cs_ids_ordered,
    x = cs_x_pos,
    center = stats::median(cs_target_x, na.rm = TRUE),
    spacing = cs_spacing,
    shrink = 0.82,
    step = glyco_grid_step / 2
  )

  # TCA center: directly below backbone (same x as AcCoA = bx = center of substrates)
  # Pass NULL to let sugar_nodes default TCA to backbone center
  # Simplified cytoplasm nodes + edges (carbs only, grid layout)
  cyto_nodes <- .dnmb_cct_sugar_nodes(cs_xs = cs_x_pos, tca_cx_override = NULL)
  cyto_nodes <- .dnmb_cct_merge_entry_into_mainstream(cyto_nodes)
  cyto_nodes <- .dnmb_cct_resolve_grid_overlaps(cyto_nodes, step = glyco_grid_step)
  cyto_edges <- .dnmb_cct_3zone_cyto_edges(cyto_nodes)

  # Prune entry_intermediate nodes not reachable from active substrates
  active_entry_targets <- unique(unname(cs_entry_map[tolower(cs_ids_ordered)]))
  active_entry_targets <- active_entry_targets[!is.na(active_entry_targets)]
  # Trace reachable intermediates: walk edges from active targets backward
  reachable <- active_entry_targets
  keep_types <- c("backbone", "ppp", "tca", "carbon_source", "ed")
  always_keep <- cyto_nodes$id[cyto_nodes$type %in% keep_types]
  # Iteratively expand reachable set through edges
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (ri in seq_len(nrow(cyto_edges))) {
      e_from <- cyto_edges$from[ri]; e_to <- cyto_edges$to[ri]
      if (e_from %in% reachable && !e_to %in% reachable) {
        reachable <- c(reachable, e_to); changed <- TRUE
      }
      if (e_to %in% reachable && !e_from %in% reachable) {
        reachable <- c(reachable, e_from); changed <- TRUE
      }
    }
  }
  reachable <- unique(c(reachable, always_keep))
  prune_ids <- cyto_nodes$id[cyto_nodes$type == "entry_intermediate" & !cyto_nodes$id %in% reachable]
  if (length(prune_ids) > 0) {
    cyto_nodes <- cyto_nodes[!cyto_nodes$id %in% prune_ids, , drop = FALSE]
    cyto_edges <- cyto_edges[!cyto_edges$from %in% prune_ids & !cyto_edges$to %in% prune_ids, , drop = FALSE]
  }

  node_x <- stats::setNames(cyto_nodes$x, cyto_nodes$id)
  node_y <- stats::setNames(cyto_nodes$y, cyto_nodes$id)

  # Assign confidence to carbon source nodes
  carbon_src <- cyto_nodes[cyto_nodes$type == "carbon_source", , drop = FALSE]
  carbon_src$confidence <- "none"
  carbon_src$fill_color <- conf_colors["none"]
  for (i in seq_len(nrow(carbon_src))) {
    pid <- tolower(carbon_src$id[i])
    if (!is.null(pstats) && pid %in% pstats$pathway_id && pid %in% matched_pathway_ids) {
      frac <- pstats$fraction[pstats$pathway_id == pid][1]
      if (frac >= 0.75) {
        carbon_src$fill_color[i] <- conf_colors["high"]; carbon_src$confidence[i] <- "high"
      } else if (frac >= 0.5) {
        carbon_src$fill_color[i] <- conf_colors["medium"]; carbon_src$confidence[i] <- "medium"
      } else if (frac > 0) {
        carbon_src$fill_color[i] <- conf_colors["low"]; carbon_src$confidence[i] <- "low"
      }
    }
  }
  sugar_lane_x <- vapply(split(carbon_src$x, carbon_src$sugar_type), function(xs) {
    .dnmb_cct_snap_to_grid(mean(xs), step = glyco_grid_step)
  }, numeric(1))

  backbone_nodes <- cyto_nodes[cyto_nodes$type == "backbone", , drop = FALSE]
  entry_inter    <- cyto_nodes[cyto_nodes$type == "entry_intermediate", , drop = FALSE]
  ppp_nodes      <- cyto_nodes[cyto_nodes$type == "ppp", , drop = FALSE]

  # ====================================================================
  # EXTRACELLULAR — substrate row definitions (polysaccharides + disaccharides only)
  # Monosaccharides shown as compact "Free uptake" strip below
  # ====================================================================
  extra_substrates <- data.frame(
    id = c("Starch", "Cellulose", "Xylan", "Chitin", "Pectin", "Mannan",
           "Inulin", "Raffinose", "Sucrose", "Maltose", "Lactose",
           "Cellobiose", "Trehalose",
           "beta-Glucan", "Arabinoxylan", "Agarose"),
    label = c("Starch/Pullulan", "Cellulose", "Xylan", "Chitin", "Pectin", "Mannan",
              "Inulin", "Raffinose", "Sucrose", "Maltose", "Lactose",
              "Cellobiose", "Trehalose",
              "beta-Glucan", "Arabinoxylan", "Agarose"),
    cascade_key = c("starch", "cellulose", "xylan", "chitin", "pectin", "mannan",
                    "inulin", "raffinose", "sucrose", "maltose", "lactose",
                    "cellobiose", "trehalose",
                    "betaglucan", "arabinoxylan", "agarose"),
    stringsAsFactors = FALSE
  )
  # SNFG monomer chain for each substrate
  extra_chains <- list(
    Starch = c("glucose", "glucose", "glucose", "glucose", "glucose"),
    Cellulose = c("glucose", "glucose", "glucose", "glucose"),
    Xylan = c("xylose", "xylose", "xylose", "xylose"),
    Chitin = c("glcnac", "glcnac", "glcnac", "glcnac"),
    Pectin = c("gala", "gala", "gala", "gala"),
    Mannan = c("mannose", "mannose", "mannose", "mannose"),
    Inulin = c("glucose", "fructose", "fructose", "fructose"),
    Raffinose = c("galactose", "glucose", "fructose"),
    Sucrose = c("glucose", "fructose"),
    Maltose = c("glucose", "glucose"),
    Lactose = c("galactose", "glucose"),
    Cellobiose = c("glucose", "glucose"),
    Trehalose = c("glucose", "glucose"),
    `beta-Glucan` = c("glucose", "glucose", "glucose", "glucose"),
    Arabinoxylan = c("xylose", "xylose", "xylose", "xylose"),
    Agarose = c("galactose", "galactose", "galactose", "galactose")
  )
  extra_bonds <- list(
    Starch = "alpha-1,4",                          # amylose: Glc-a1,4-Glc
    Cellulose = "beta-1,4",                         # Glc-b1,4-Glc
    Xylan = "beta-1,4",                             # Xyl-b1,4-Xyl
    Chitin = "beta-1,4",                            # GlcNAc-b1,4-GlcNAc
    Pectin = "alpha-1,4",                           # GalA-a1,4-GalA
    Mannan = "beta-1,4",                            # Man-b1,4-Man
    Inulin = "beta-2,1",                            # Fru-b2,1-Fru
    Raffinose = c("alpha-1,6", "alpha-1,2"),        # Gal-a1,6-Glc-a1,2-Fru
    Sucrose = "alpha-1,2",                          # Glc-a1,b2-Fru
    Maltose = "alpha-1,4",                          # Glc-a1,4-Glc
    Lactose = "beta-1,4",                           # Gal-b1,4-Glc
    Cellobiose = "beta-1,4",                        # Glc-b1,4-Glc
    Trehalose = "alpha-1,1",                         # Glc-a1,1-Glc
    `beta-Glucan` = "beta-1,3",                      # Glc-b1,3-Glc (laminarin)
    Arabinoxylan = "beta-1,4",                        # Xyl-b1,4-Xyl + Ara branches
    Agarose = c("beta-1,3", "alpha-1,4")                # D-Gal-β1,3-anhydro-L-Gal-α1,4 alternating
  )
  # Branches: list of (position_on_backbone, branch_sugar, bond_type)
  # position = 1-based index of backbone monomer where branch attaches
  # Branches: PubChem-accurate side chains on polysaccharide backbones
  extra_branches <- list(
    # Amylopectin: α1,4 backbone + α1,6 branch; debranch = pullulanase/isoamylase (GH13)
    Starch = list(list(pos = 3L, sugar = "glucose", bond = "alpha-1,6",
                       debranch_gh = c("GH13"))),
    # Xylan: β1,4 Xyl + α1,3-Ara; debranch = arabinosidase (GH43, GH51, GH62)
    Xylan  = list(list(pos = 2L, sugar = "arabinose", bond = "alpha-1,3",
                       debranch_gh = c("GH43", "GH51", "GH62"))),
    # Pectin RG-I: GalA + Rha + Gal; debranch = rhamnosidase/galactosidase
    Pectin = list(list(pos = 2L, sugar = "rhamnose",  bond = "alpha-1,2",
                       debranch_gh = c("GH78", "GH106")),
                  list(pos = 3L, sugar = "galactose", bond = "beta-1,4",
                       debranch_gh = c("GH35", "GH42"))),
    # Galactomannan: Man + α1,6-Gal; debranch = α-galactosidase (GH27, GH36)
    Mannan = list(list(pos = 2L, sugar = "galactose", bond = "alpha-1,6",
                       debranch_gh = c("GH27", "GH36"))),
    # Arabinoxylan: Xyl backbone + α1,3-Ara branches
    Arabinoxylan = list(list(pos = 2L, sugar = "arabinose", bond = "alpha-1,3",
                             debranch_gh = c("GH43", "GH51", "GH62")))
  )

  # Monosaccharides for "Free uptake" strip
  free_sugars <- list(
    list(label = "Glc",  sugar = "glucose"),
    list(label = "Gal",  sugar = "galactose"),
    list(label = "Man",  sugar = "mannose"),
    list(label = "Fru",  sugar = "fructose"),
    list(label = "GlcNAc", sugar = "glcnac"),
    list(label = "Xyl",  sugar = "xylose"),
    list(label = "Ara",  sugar = "arabinose"),
    list(label = "Fuc",  sugar = "fucose"),
    list(label = "Rha",  sugar = "rhamnose"),
    list(label = "GlcA", sugar = "glca"),
    list(label = "GalA", sugar = "gala")
  )

  # ====================================================================
  # TOP-DOWN LAYOUT (landscape): Extracellular → Membrane → Cytoplasm
  # y = 9~15  : Extracellular
  # y = 8.5   : Membrane (horizontal double line)
  # y = 0~8   : Cytoplasm
  # ====================================================================
  # Filter extracellular substrates: only show those with CAZy GH evidence
  # (cascade requires actual GH enzyme match, not just transport evidence)
  extra_has_evidence <- vapply(extra_substrates$cascade_key, function(ck) {
    ck %in% names(active_cascades)
  }, logical(1))
  extra_substrates <- extra_substrates[extra_has_evidence, , drop = FALSE]
  extra_chains <- extra_chains[extra_substrates$id]
  extra_bonds <- extra_bonds[extra_substrates$id]
  extra_branches <- extra_branches[intersect(names(extra_branches), extra_substrates$id)]

  n_extra  <- nrow(extra_substrates)
  sym_size <- 0.075
  chain_gap <- 0.18
  chain_step <- .dnmb_cct_extra_chain_step(sym_size = sym_size, gap = chain_gap)
  extra_trim <- .dnmb_cct_extra_route_trim(sym_size = sym_size)
  y_memb   <- 8.5   # horizontal membrane y
  # Count max cascade depth for polysaccharides (for vertical space)
  max_depth <- max(vapply(cascades, function(c) length(c$chains), integer(1)), na.rm = TRUE)
  y_bot    <- min(cyto_nodes$y) - 0.5

  # ====================================================================
  # HIERARCHY GRID: 4 levels (Polymer → Oligomer → Dimer → Monomer)
  # y-levels fixed, x-columns by sugar family
  # ====================================================================
  # ====================================================================
  # GRID LAYOUT: y = backbone monomer count, x = horizontal spread
  # Same count → same y, different x (side by side)
  # ====================================================================
  step_y <- 0.5  # y-spacing per monomer-count level

  # y-positions by backbone monomer count (branch excluded)
  y_6mer <- y_memb + 4.5  # 6-mer polymers (Pullulan)
  y_5mer <- y_memb + 4.0  # 5-mer polymers (Starch)
  y_4mer <- y_memb + 3.5  # 4-mer polymers (Chitin, Xylan, Pectin, Mannan, Cellulose)
  y_3mer <- y_memb + 3.0  # 3-mer (Raffinose, Inulin)
  y_2mer <- y_memb + 2.0  # 2-mer dimers
  y_mono <- y_memb + 1.0  # 1-mer monomers
  y_extra_top <- y_6mer + 0.5

  row_y_map <- c(`6` = y_6mer, `5` = y_5mer, `4` = y_4mer, `3` = y_3mer, `2` = y_2mer)
  extra_layout <- .dnmb_cct_extra_layout(
    extra_substrates = extra_substrates,
    extra_chains = extra_chains,
    sugar_lane_x = sugar_lane_x,
    row_y_map = row_y_map,
    extra_branches = extra_branches,
    cascades = cascades,
    sym_size = sym_size
  )
  extra_x_map <- stats::setNames(extra_layout$x_start, extra_layout$id)
  extra_center_map <- stats::setNames(extra_layout$x_center, extra_layout$id)
  extra_y_map <- stats::setNames(extra_layout$y, extra_layout$id)

  # x_max adapts to both cytoplasm lanes and extracellular chain width
  x_max <- max(max(cs_x_pos) + 0.6, max(extra_layout$x_right) + 0.2)

  # Pathway colors (needed by cascade rendering + edge rendering)
  pw_colors <- .dnmb_cct_pathway_colors()
  node_sugar_type <- stats::setNames(cyto_nodes$sugar_type, cyto_nodes$id)
  mem_xmin <- min(
    c(extra_layout$x_left, extra_layout$x_center, carbon_src$x, cyto_nodes$x),
    na.rm = TRUE
  ) - 0.20
  mem_xmax <- max(
    c(extra_layout$x_right, extra_layout$x_center, carbon_src$x, cyto_nodes$x),
    na.rm = TRUE
  ) + 0.20

  # ====================================================================
  # BUILD THE PLOT — clean canvas
  # ====================================================================
  p <- ggplot2::ggplot()

  # --- Horizontal phospholipid bilayer membrane ---
  memb_xs <- seq(mem_xmin, mem_xmax, by = 0.35)
  head_outer <- data.frame(x = memb_xs, y = y_memb + 0.18)
  head_inner <- data.frame(x = memb_xs, y = y_memb - 0.18)
  tails_outer <- data.frame(x = memb_xs, xend = memb_xs,
                             y = y_memb + 0.14, yend = y_memb + 0.03)
  tails_inner <- data.frame(x = memb_xs, xend = memb_xs,
                             y = y_memb - 0.14, yend = y_memb - 0.03)
  # Light fill between leaflets
  p <- p + ggplot2::geom_rect(
    data = data.frame(xmin = mem_xmin, xmax = mem_xmax,
                      ymin = y_memb - 0.20, ymax = y_memb + 0.20),
    ggplot2::aes(xmin = .data$xmin, ymin = .data$ymin,
                 xmax = .data$xmax, ymax = .data$ymax),
    fill = "#FFF8F0", color = NA, inherit.aes = FALSE)
  # Tails
  p <- p +
    ggplot2::geom_segment(data = tails_outer,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 0.25, color = "#BCAAA4", inherit.aes = FALSE) +
    ggplot2::geom_segment(data = tails_inner,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 0.25, color = "#BCAAA4", inherit.aes = FALSE)
  # Head groups
  p <- p +
    ggplot2::geom_point(data = head_outer,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 1.2, fill = "#D7CCC8", color = "#A1887F",
      stroke = 0.25, inherit.aes = FALSE) +
    ggplot2::geom_point(data = head_inner,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 1.2, fill = "#D7CCC8", color = "#A1887F",
      stroke = 0.25, inherit.aes = FALSE)

  # --- Zone labels ---
  p <- p +
    ggplot2::annotate("text", x = x_max + 0.3, y = y_extra_top + 0.5,
                      label = "Extracellular", hjust = 1, size = 3,
                      fontface = "italic", color = "#BDBDBD") +
    ggplot2::annotate("text", x = x_max + 0.3, y = y_memb,
                      label = "Membrane", hjust = 1, size = 2.5,
                      fontface = "italic", color = "#A1887F") +
    ggplot2::annotate("text", x = x_max + 0.3, y = y_memb - 0.5,
                      label = "Cytoplasm", hjust = 1, size = 3,
                      fontface = "italic", color = "#BDBDBD")

  # ====================================================================
  # MONOMER LEVEL: Draw unique monomers ONCE at y_mono
  # All cascade products converge to these nodes
  # ====================================================================
  all_prods <- unique(.dnmb_cct_canonical_sugar_id(
    unlist(lapply(cascades, function(c) c$products %||% c$product))
  ))
  # Monomer hub x = carbon source hub x (vertical alignment through membrane)
  # Reuse the canonical sugar lanes so outside/inside sugars align vertically
  cs_hub_x <- sugar_lane_x
  composite_product_ids <- c("sucrose", "trehalose", "lactose", "maltose", "cellobiose", "raffinose")
  # Map monomer IDs to their carbon source hub x
  mono_ids <- c("glucose","galactose","mannose","fructose","glcnac",
                "xylose","arabinose","ribose","fucose","rhamnose","glca","gala")
  mono_ids <- unique(c(mono_ids, setdiff(all_prods, composite_product_ids)))
  # Sugar type → monomer ID mapping
  mono_to_sugar <- c(glucose="glucose", galactose="galactose", mannose="mannose",
    fructose="fructose", glcnac="glcnac", xylose="xylose", arabinose="arabinose",
    ribose="ribose", fucose="fucose", rhamnose="fucose", glca="glca", gala="glca",
    glycerol="glycerol", glucosamine="glcnac")
  mono_to_csid <- c(
    sucrose = "Sucrose",
    trehalose = "Trehalose",
    lactose = "Lactose",
    maltose = "Maltose",
    cellobiose = "Cellobiose",
    raffinose = "Raffinose"
  )
  # Assign x from carbon source hub, or spread remaining evenly
  mono_xs <- vapply(mono_ids, function(mid) {
    csid <- unname(mono_to_csid[mid])
    if (!is.na(csid) && csid %in% carbon_src$id) {
      return(unname(carbon_src$x[match(csid, carbon_src$id)]))
    }
    st <- mono_to_sugar[mid]
    if (!is.na(st) && st %in% names(cs_hub_x)) return(unname(cs_hub_x[st]))
    NA_real_
  }, numeric(1))
  # Fill NA with even spacing in remaining gaps
  na_idx <- which(is.na(mono_xs))
  if (length(na_idx) > 0) {
    used <- mono_xs[!is.na(mono_xs)]
    gap_xs <- seq(max(used) + 1.5, max(used) + 1.5 + length(na_idx) * 1.5, length.out = length(na_idx))
    mono_xs[na_idx] <- gap_xs
  }
  mono_x_map <- stats::setNames(mono_xs, mono_ids)
  extra_target_id_map <- stats::setNames(extra_substrates$id, tolower(extra_substrates$id))

  for (mi in seq_along(mono_ids)) {
    mid <- mono_ids[mi]
    mx <- mono_xs[mi]
    mono_col <- .dnmb_cct_sugar_route_color(mid)
    # First render: no label (final overlay at end will add labels above route lines)
    p <- .dnmb_snfg_render_symbol_v2(p, mx, y_mono, mid, r = 0.10, label = NULL)
    # Arrow from monomer to membrane
    p <- p + ggplot2::geom_segment(
      data = data.frame(x = mx, xend = mx,
                        y = y_mono, yend = y_memb + 0.25),
      ggplot2::aes(x=.data$x, xend=.data$xend, y=.data$y, yend=.data$yend),
      arrow = ggplot2::arrow(length = grid::unit(0.02, "inches"), type = "closed"),
      linewidth = 0.2, color = mono_col, alpha = 0.7, inherit.aes = FALSE)
  }

  # ====================================================================
  # ZONE 1: EXTRACELLULAR — SNFG chains + GH scissors on chain
  # ====================================================================
  gh_label_rows <- list()
  for (ri in seq_len(n_extra)) {
    sub   <- extra_substrates[ri, , drop = FALSE]
    ry    <- extra_y_map[sub$id]
    chain_x <- extra_x_map[sub$id]  # per-substrate x position
    chain_cx <- extra_center_map[sub$id]
    chain <- extra_chains[[sub$id]]
    bond_type <- extra_bonds[[sub$id]]

    # Substrate name (above chain)
    p <- p + ggplot2::annotate("text", x = chain_cx, y = ry + 0.25, label = sub$label,
                                hjust = 0.5, size = 1.5, color = "#333333", fontface = "bold")

    # SNFG glycan chain
    n_mono <- length(chain)
    if (n_mono > 0) {
      bonds_vec <- if (n_mono > 1 && !is.null(bond_type)) rep(bond_type, n_mono - 1) else NULL
      chain_layers <- .dnmb_draw_sugar_chain_v2(
        x_start = chain_x, y = ry, monomers = chain,
        bonds = bonds_vec, size = sym_size, gap = chain_gap, show_label = FALSE)
      for (ly in chain_layers) p <- p + ly

      # Draw branches (side chains hanging below the backbone)
      branches <- extra_branches[[sub$id]]
      if (!is.null(branches)) {
        for (bi in seq_along(branches)) {
          br <- branches[[bi]]
          br_x <- .dnmb_cct_extra_chain_node_x(
            chain_x = chain_x,
            pos = br$pos,
            sym_size = sym_size,
            gap = chain_gap
          )  # x of backbone monomer center
          br_y <- ry - sym_size * (2.5 + 1.4 * (bi - 1))  # stagger branch depth
          # Vertical bond line
          is_beta <- grepl("^beta", br$bond)
          lty <- if (is_beta) "dashed" else "solid"
          # Parse bond label: "alpha-1,6" → "α1–6"
          br_label <- sub("^alpha-", "\u03b1", sub("^beta-", "\u03b2", br$bond))
          br_label <- gsub(",", "\u2013", br_label)  # comma → en-dash
          p <- p + ggplot2::geom_segment(
            data = data.frame(x = br_x, xend = br_x, y = ry - sym_size, yend = br_y + sym_size),
            ggplot2::aes(x=.data$x, xend=.data$xend, y=.data$y, yend=.data$yend),
            linewidth = 0.3, color = "#795548", linetype = lty, inherit.aes = FALSE)
          # Branch sugar symbol
          abbr_br <- .dnmb_snfg_abbreviation_v2(br$sugar)
          p <- .dnmb_snfg_render_symbol_v2(p, br_x, br_y, br$sugar, r = sym_size, label = abbr_br)
          # Bond label
          p <- p + ggplot2::annotate("text", x = br_x + 0.05, y = (ry + br_y)/2,
                                      label = br_label, size = 0.9, color = "#795548",
                                      hjust = 0, fontface = "plain")

          # Debranching enzyme: scissors on the branch bond + GH label + arrow to monomer hub
          debranch_ghs <- br$debranch_gh
          if (!is.null(debranch_ghs) && !is.null(gh_enzymes) && nrow(gh_enzymes) > 0) {
            db_matched <- .dnmb_cct_3zone_match_gh(gh_enzymes, debranch_ghs)
            if (nrow(db_matched) > 0) {
              # Small scissors on the branch bond midpoint
              sc_y <- (ry + br_y) / 2
              sc <- .dnmb_scissors_grob_v2(br_x - 0.1, sc_y, size = 0.03)
              for (ly in sc) p <- p + ly
              # Debranching enzyme label
              db_lab <- db_matched$gh_family[1]
              if (!is.na(db_matched$gene_name[1]) && nzchar(db_matched$gene_name[1]))
                db_lab <- .dnmb_cct_short_gh_label(db_matched$gh_family[1], db_matched$gene_name[1])
              db_priority <- 12 + min(4, nrow(db_matched)) +
                ifelse(sub$id %in% c("Pectin", "Mannan"), 2, 0)
              gh_label_rows[[length(gh_label_rows) + 1L]] <- data.frame(
                x = br_x - 0.18, y = sc_y,
                label = db_lab, color = "#1565C0", priority = db_priority,
                stringsAsFactors = FALSE
              )
              # Arrow from detached branch to monomer hub
              if (tolower(br$sugar) %in% names(mono_x_map)) {
                target_mx <- unname(mono_x_map[tolower(br$sugar)])
                mid_y2 <- y_mono + sym_size * 1.5
                branch_col <- .dnmb_cct_sugar_route_color(br$sugar, fallback = "#1565C0")
                branch_imp <- min(1.0, 0.45 + 0.18 * nrow(db_matched))
                branch_path <- .dnmb_cct_rounded_route_points(
                  data.frame(
                    x = c(br_x, br_x, target_mx, target_mx),
                    y = c(br_y, mid_y2, mid_y2, y_mono)
                  ),
                  radius = 0.2
                )
                for (ly in .dnmb_cct_gradient_path_layers(
                  branch_path, color = branch_col,
                  linewidth = 0.18 + 0.18 * branch_imp,
                  alpha = 0.22 + 0.28 * branch_imp,
                  linetype = "dashed", arrow_last = TRUE, arrow_length = 0.015,
                  trim_start = extra_trim, trim_end = 0.11
                )) p <- p + ly
              }
            }
          }
        }
      }
    }

    # GH scissors on the chain + arrow to monomer hub (no sub-level chains)
    ckey <- sub$cascade_key
    if (!is.na(ckey) && ckey %in% names(active_cascades)) {
      casc <- active_cascades[[ckey]]
      matched <- casc$matched_enzymes
      product_sugar <- .dnmb_cct_canonical_sugar_id(casc$product)[1]

      # GH scissors on bond line
      cut_x <- .dnmb_cct_cascade_cut_x(
        cascade_key = ckey,
        chain_x = chain_x,
        n_mono = n_mono,
        step = chain_step,
        matched_families = matched$gh_family
      )
      sc <- .dnmb_scissors_grob_v2(cut_x, ry, size = 0.04)
      for (ly in sc) p <- p + ly

      # GH label below scissors (with CGC/PUL annotation if available)
      show_n <- min(if (n_mono >= 4) 2L else 1L, nrow(matched))
      for (j in seq_len(show_n)) {
        enz <- matched[j, , drop = FALSE]
        gh_lab <- .dnmb_cct_short_gh_label(enz$gh_family, enz$gene_name)
        # Append PUL substrate or CGC ID if available
        pul_sub <- NULL
        for (pc in intersect(c("dbCAN_dbcan_pul_substrate", "dbcan_pul_substrate"), names(enz))) {
          val <- as.character(enz[[pc]])
          if (!is.na(val) && nzchar(val)) { pul_sub <- val; break }
        }
        if (is.null(pul_sub)) {
          for (sc in intersect(c("dbCAN_dbcan_sub_substrate", "dbcan_sub_substrate"), names(enz))) {
            val <- as.character(enz[[sc]])
            if (!is.na(val) && nzchar(val)) { pul_sub <- val; break }
          }
        }
        if (!is.null(pul_sub)) {
          pul_short <- substr(trimws(strsplit(pul_sub, ";")[[1]][1]), 1, 12)
          gh_lab <- paste0(gh_lab, " [", pul_short, "]")
        }
        main_priority <- 10 - j + ifelse(sub$id %in% c("Pectin", "Mannan"), 2, 0) +
          ifelse(grepl("GH28|GH35|GH42|GH78|GH106|GH27|GH36|GH26|GH113", enz$gh_family), 1, 0)
        gh_label_rows[[length(gh_label_rows) + 1L]] <- data.frame(
          x = cut_x, y = ry - 0.15 - (j - 1) * 0.10,
          label = gh_lab, color = "#C62828",
          priority = main_priority, stringsAsFactors = FALSE
        )
      }

      # Arrow from chain → monomer hub: pathway color, 0.5 grid snapped
      prod_sugars <- cascades[[ckey]]$products
      if (is.null(prod_sugars)) prod_sugars <- product_sugar
      prod_sugars_u <- unique(.dnmb_cct_canonical_sugar_id(prod_sugars))
      prod_path_map <- c(
        glucose = "glucose", galactose = "galactose", mannose = "mannose",
        fructose = "fructose", glcnac = "nag", xylose = "xylose",
        arabinose = "arabinose", ribose = "ribose", fucose = "fucose",
        glca = "gluconate", gala = "galacturonate", sucrose = "sucrose"
      )
      supported_minor <- vapply(prod_sugars_u, function(ps) {
        pid <- unname(prod_path_map[ps])
        !is.na(pid) && tolower(pid) %in% matched_pathway_ids
      }, logical(1))
      keep_prod <- prod_sugars_u == product_sugar | supported_minor
      if (any(keep_prod)) {
        prod_sugars_u <- prod_sugars_u[keep_prod]
      }
      # Get pathway color for this cascade
      casc_col <- pw_colors[ckey]
      if (is.na(casc_col)) casc_col <- "#795548"
      for (ps in prod_sugars_u) {
        target_id <- unname(extra_target_id_map[ps])
        target_x <- if (ps %in% names(mono_x_map)) {
          unname(mono_x_map[ps])
        } else if (!is.na(target_id) && target_id %in% names(extra_center_map)) {
          unname(extra_center_map[target_id])
        } else {
          NA_real_
        }
        target_y <- if (ps %in% names(mono_x_map)) {
          y_mono
        } else if (!is.na(target_id) && target_id %in% names(extra_y_map)) {
          unname(extra_y_map[target_id])
        } else {
          NA_real_
        }
        trim_end_val <- if (ps %in% names(mono_x_map)) 0.11 else extra_trim
        if (is.finite(target_x) && is.finite(target_y)) {
          route_col <- .dnmb_cct_sugar_route_color(ps, fallback = casc_col)
          route_wt <- .dnmb_cct_product_route_weight(
            product_sugar = ps,
            primary_product = casc$product,
            product_vec = prod_sugars_u
          )
          # Exit x: use the GH cleavage site (cut_x) as the release point.
          # For exo-acting enzymes, cut is at the non-reducing end (rightmost);
          # for endo-acting enzymes, cut is internal. The monomer exits from
          # the cleavage position, not from an arbitrary chain match.
          exit_x <- cut_x
          # Snap mid_y to 0.5 grid with micro-stagger to avoid overlapping horizontal lines
          ps_idx <- match(ps, prod_sugars_u)
          if (is.na(ps_idx)) ps_idx <- 1L
          micro_stagger <- (ps_idx - 1L) * 0.06 + (ri - 1L) * 0.04
          mid_y3 <- round(((target_y + extra_trim) + 0.5) * 2) / 2 + micro_stagger
          cascade_path <- .dnmb_cct_rounded_route_points(
            data.frame(
              x = c(exit_x, exit_x, target_x, target_x),
              y = c(ry, mid_y3, mid_y3, target_y)
            ),
            radius = 0.5
          )
          # Simple colored path (no gradient, unified width)
          p <- p + .dnmb_cct_single_path_layer(
            cascade_path, color = route_col,
            linewidth = 0.35,
            alpha = 0.45,
            arrow_last = TRUE, arrow_length = 0.02,
            trim_start = extra_trim, trim_end = trim_end_val
          )
        }
      }

    } else {
      # No cascade match — just draw product arrow to membrane if chain exists
      chain_end_x <- chain_x + max(0, n_mono - 1) * chain_step
      p <- p + ggplot2::geom_segment(
        data = data.frame(x = chain_end_x, xend = chain_end_x,
                          y = ry, yend = y_memb + 0.25),
        ggplot2::aes(x=.data$x, xend=.data$xend, y=.data$y, yend=.data$yend),
        arrow = ggplot2::arrow(length = grid::unit(0.02, "inches"), type = "closed"),
        linewidth = 0.2, color = "#999999", inherit.aes = FALSE)
    }
  }

  if (length(gh_label_rows) > 0) {
    gh_label_df <- do.call(rbind, gh_label_rows)
    gh_label_df <- .dnmb_cct_thin_labels(gh_label_df, x_thresh = 0.55, y_thresh = 0.14)
    gh_label_df <- .dnmb_cct_place_small_labels(gh_label_df, x_thresh = 0.55, y_thresh = 0.14)
    p <- p + ggplot2::geom_text(
      data = gh_label_df,
      ggplot2::aes(x = .data$x_lab, y = .data$y_lab, label = .data$label),
      size = 0.82, color = gh_label_df$color,
      fontface = "bold", hjust = gh_label_df$hjust,
      lineheight = 0.82, inherit.aes = FALSE
    )
  }

  # ====================================================================
  # ZONE 2: MEMBRANE — Transporters aligned with carbon source x positions
  # PTS transporters (orange dashed), non-PTS (solid colored arrows)
  # ====================================================================
  # Map carbon source IDs to x positions for alignment
  cs_x_map <- stats::setNames(carbon_src$x, carbon_src$id)
  # Also map pathway names to carbon source IDs
  pathway_to_cs <- c(
    maltose = "Maltose", cellobiose = "Cellobiose", galactose = "Galactose",
    trehalose = "Trehalose", lactose = "Lactose", mannose = "Mannose",
    NAG = "NAG", glucosamine = "Glucosamine", fructose = "Fructose",
    sucrose = "Sucrose", mannitol = "Mannitol", glycerol = "Glycerol",
    xylose = "Xylose", arabinose = "Arabinose", ribose = "Ribose",
    fucose = "Fucose", rhamnose = "Rhamnose", gluconate = "Gluconate",
    glucose = "Maltose", deoxyribose = "Ribose", deoxyribonate = "Gluconate")

  if (!is.null(transporters) && nrow(transporters) > 0) {
    transporters <- transporters[!duplicated(transporters$locus_tag), , drop = FALSE]
    # Determine if PTS
    transporters$is_pts <- grepl("pts|PTS|EIIC|EIIB|crr", transporters$step, ignore.case = TRUE) |
      grepl("pts|PTS", transporters$gene_name, ignore.case = TRUE)
    transporters$cs_id <- unname(pathway_to_cs[transporters$pathway])
    transporters <- transporters[!is.na(transporters$cs_id) & transporters$cs_id %in% names(cs_x_map), , drop = FALSE]
    if (nrow(transporters) > 0) {
      if (!"step_score" %in% names(transporters)) transporters$step_score <- NA_real_
      transporters <- .dnmb_cct_annotate_transport_context(transporters, genbank_table)
      transporters$tx <- unname(cs_x_map[transporters$cs_id])
      tr_entry_map <- .dnmb_cct_auto_entry_map()
      tr_target_nodes <- rbind(
        entry_inter,
        backbone_nodes,
        ppp_nodes,
        cyto_nodes[cyto_nodes$type == "ed", , drop = FALSE]
      )
      tr_target_x <- stats::setNames(tr_target_nodes$x, tr_target_nodes$id)
      transporters$entry_target <- unname(tr_entry_map[tolower(transporters$cs_id)])
      transporters$center_x <- transporters$tx
      transporters$desired_x <- transporters$tx
      has_target <- !is.na(transporters$entry_target) & transporters$entry_target %in% names(tr_target_x)
      if (any(has_target)) {
        tx_target <- unname(tr_target_x[transporters$entry_target[has_target]])
        delta_x <- tx_target - transporters$tx[has_target]
        transporters$center_x[has_target] <- transporters$tx[has_target] +
          pmax(-0.75, pmin(0.75, 0.18 * delta_x))
        transporters$desired_x[has_target] <- transporters$tx[has_target] +
          pmax(-1.00, pmin(1.00, 0.32 * delta_x))
      }
      transporters$matched <- !is.na(transporters$locus_tag) & nzchar(transporters$locus_tag)
      transporters$conf_rank <- conf_rank[transporters$confidence]
      transporters$conf_rank[is.na(transporters$conf_rank)] <- 0L
      if (!"kind" %in% names(transporters)) {
        transporters$kind <- mapply(.dnmb_cct_transporter_kind, transporters$step, transporters$gene_name, is_pts)
      }
      transporters <- .dnmb_cct_select_transporters(transporters, max_per_lane = 8L)
      split_lane <- split(seq_len(nrow(transporters)), paste(transporters$cs_id, round(transporters$tx, 3)))
      transporters$lane_rank <- 1L
      transporters$lane_top_score <- transporters$display_score
      lane_center_ref <- stats::median(cs_x_map, na.rm = TRUE)
      transporters$central_lane <- abs(transporters$tx - lane_center_ref) <= 2.2
      for (grp in split_lane) {
        ord <- grp[order(transporters$display_score[grp], transporters$matched[grp], transporters$conf_rank[grp],
                         transporters$gene_name[grp], decreasing = TRUE, na.last = TRUE)]
        transporters$lane_rank[ord] <- seq_along(ord)
        transporters$lane_top_score[grp] <- max(transporters$display_score[grp], na.rm = TRUE)
      }
      transporters$show_ctx <- !is.na(transporters$context_score) &
        transporters$context_score >= 1.0 &
        (transporters$lane_rank == 1L | transporters$context_score >= 2.2)
      transporters$show_label <- transporters$matched | transporters$lane_rank == 1L
      transporters$ctx_lab <- ifelse(
        transporters$show_ctx,
        vapply(transporters$context_hits, .dnmb_cct_short_context_summary, character(1)),
        ""
      )
      transporters$label_text <- vapply(seq_len(nrow(transporters)), function(i) {
        tr <- transporters[i, , drop = FALSE]
        human_lab <- .dnmb_cct_prefer_human_label(tr$step, tr$gene_name, pathway_id = tr$pathway)
        if (tr$lane_rank == 1L) {
          if (!is.na(human_lab) && nzchar(human_lab)) {
            paste0(
              human_lab,
              ifelse(tr$is_pts, " (PTS)", ""),
              ifelse(tr$matched, paste0("\n", tr$locus_tag), ""),
              ifelse(nzchar(tr$ctx_lab), paste0("\n", tr$ctx_lab), "")
            )
          } else if (tr$matched) {
            paste0(tr$locus_tag, ifelse(nzchar(tr$ctx_lab), paste0("\n", tr$ctx_lab), ""))
          } else {
            paste0(.dnmb_cct_short_step_label(tr$step, tr$pathway),
                   ifelse(nzchar(tr$ctx_lab), paste0("\n", tr$ctx_lab), ""))
          }
        } else {
          if (tr$matched) {
            paste0(human_lab, ifelse(nzchar(human_lab), "\n", ""), tr$locus_tag)
          } else if (!is.na(human_lab) && nzchar(human_lab)) {
            paste0(human_lab, ifelse(tr$is_pts, " (PTS)", ""))
          } else if (tr$matched) {
            tr$locus_tag
          } else {
            .dnmb_cct_short_step_label(tr$step, tr$pathway)
          }
        }
      }, character(1))
      transporters$label_width <- vapply(transporters$label_text, .dnmb_cct_estimate_text_width, numeric(1))
      transporters$seg_half <- mapply(.dnmb_cct_transporter_half_span, transporters$step, transporters$gene_name, transporters$is_pts)
      transporters$label_dx <- 0
      transporters$tx_draw <- transporters$tx
      transporters$ty_draw <- 8.56
      transporters$row_draw <- 1L
      global_cs_center <- stats::median(carbon_src$x, na.rm = TRUE)
      for (grp in split_lane) {
        ord <- grp[order(transporters$lane_rank[grp])]
        lane_dir <- sign(mean(transporters$tx[ord], na.rm = TRUE) - global_cs_center)
        if (!is.finite(lane_dir) || lane_dir == 0) {
          lane_dir <- sign(mean(transporters$desired_x[ord] - transporters$tx[ord], na.rm = TRUE))
        }
        ord_dx <- .dnmb_cct_lane_label_offsets(length(ord), lane_dir = lane_dir)
        transporters$label_dx[ord] <- ord_dx
        packed <- .dnmb_cct_pack_transporters_lane(
          center_x = transporters$center_x[ord[1]],
          half_spans = transporters$seg_half[ord],
          lane_ranks = transporters$display_score[ord],
          label_widths = transporters$label_width[ord],
          label_dx = ord_dx,
          desired_x = transporters$desired_x[ord]
        )
        transporters$tx_draw[ord] <- packed$tx
        transporters$ty_draw[ord] <- packed$ty
        transporters$row_draw[ord] <- packed$row_id
      }
      transporter_label_rows <- list()

      for (i in seq_len(nrow(transporters))) {
        tr <- transporters[i, , drop = FALSE]
        tx <- tr$tx_draw
        conf <- if (!is.na(tr$confidence)) tr$confidence else "medium"
        sugar_type <- tolower(as.character(carbon_src$sugar_type[match(tr$cs_id, carbon_src$id)]))
        frac_val <- if (!is.null(pstats) && tolower(tr$pathway) %in% pstats$pathway_id) {
          pstats$fraction[pstats$pathway_id == tolower(tr$pathway)][1]
        } else {
          NA_real_
        }
        tr_importance <- .dnmb_cct_path_importance(confidence = conf, fraction = frac_val, sink_weight = 1)
        arr_col <- if (isTRUE(tr$matched)) {
          .dnmb_cct_sugar_route_color(sugar_type, fallback = "#666666")
        } else {
          "#BDBDBD"
        }
        seg_alpha <- switch(conf, high = 0.9, medium = 0.75, low = 0.55, 0.4)
        tr_kind <- tr$kind
        seg_half <- tr$seg_half
        ty <- tr$ty_draw
        lbl_dx <- tr$label_dx
        base_lbl_y <- if (tr$row_draw == 1) y_memb + 0.28 else y_memb - 0.28
        lbl_y <- base_lbl_y + if (tr$row_draw == 1) 0.04 * ((tr$lane_rank - 1) %/% 2) else -0.04 * ((tr$lane_rank - 1) %/% 2)

        for (ly in .dnmb_cct_transporter_glyph_layers(
          tx = tx, ty = ty, half_span = seg_half, core_color = arr_col,
          confidence = conf, is_pts = tr$is_pts, kind = tr_kind
        )) p <- p + ly
        if (isTRUE(tr$matched)) {
          for (ly in .dnmb_cct_junction_glyph_layers(
            x = tx, y = if (tr$row_draw == 1) y_memb + 0.18 else y_memb - 0.18, color = arr_col,
            size = 1.4 + 0.8 * tr_importance, alpha = 0.65 + 0.2 * tr_importance
          )) p <- p + ly
        }

        if (isTRUE(tr$show_label)) {
          transporter_label_rows[[length(transporter_label_rows) + 1L]] <- data.frame(
            x = tx + lbl_dx * ifelse(tr$lane_rank == 1L, 0.75, 0.55),
            y = lbl_y,
            label = tr$label_text,
            color = arr_col,
            lane_rank = tr$lane_rank,
            priority = if (isTRUE(tr$matched)) 10 * tr_importance + tr$display_score else tr$display_score,
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(transporter_label_rows) > 0) {
        tlab_df <- do.call(rbind, transporter_label_rows)
        tlab_df <- .dnmb_cct_place_small_labels(tlab_df, x_thresh = 0.34, y_thresh = 0.08)
        p <- p + ggplot2::geom_text(
          data = tlab_df,
          ggplot2::aes(x = .data$x_lab, y = .data$y_lab, label = .data$label),
          size = 0.58, color = tlab_df$color,
          fontface = "bold", lineheight = 0.84,
          inherit.aes = FALSE
        )
      }
    }
  }

  # ====================================================================
  # ZONE 3: CYTOPLASM — Simplified metabolism (carbs only, grid layout)
  # ====================================================================

  # Draw curved edges — gapmind_aa-style glow ribbons for backbone,
  # colored curves for substrate pathways
  pw_colors <- .dnmb_cct_pathway_colors()
  glow_lw    <- c(5, 3.5, 2, 1.0)
  glow_alpha <- c(0.04, 0.08, 0.14, 0.20)
  edge_entry_map <- .dnmb_cct_auto_entry_map()
  route_hint_cache <- new.env(parent = emptyenv())
  route_hint_nodes_for <- function(path_ids) {
    key <- paste(sort(unique(tolower(as.character(path_ids)))), collapse = "|")
    if (!nzchar(key)) return(character(0))
    if (exists(key, envir = route_hint_cache, inherits = FALSE)) {
      return(get(key, envir = route_hint_cache, inherits = FALSE))
    }
    val <- .dnmb_cct_pathway_hint_nodes(
      path_ids = unlist(strsplit(key, "|", fixed = TRUE)),
      matched_steps = matched_steps,
      transporters = transporters,
      entry_map = edge_entry_map
    )
    assign(key, val, envir = route_hint_cache)
    val
  }

  if (!is.null(cyto_edges) && nrow(cyto_edges) > 0) {
    for (i in seq_len(nrow(cyto_edges))) {
      ce <- cyto_edges[i, , drop = FALSE]
      pw  <- ce$pathway

      if (identical(pw, "tca")) next  # TCA circle handles visuals

      # All edges: simple solid gray lines with L-shaped routing
      x1 <- ce$x; y1 <- ce$y; x2 <- ce$xend; y2 <- ce$yend
      if (any(is.na(c(x1, y1, x2, y2)))) next
      is_straight <- abs(x2 - x1) < 0.01 || abs(y2 - y1) < 0.01

      if (identical(pw, "backbone")) {
        # Backbone: solid gray segment, trimmed to avoid passing through node symbols
        bb_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
        bb_trim <- min(0.12, max(0.10, bb_len * 0.15))
        if (bb_len > 2 * bb_trim) {
          bb_frac_s <- bb_trim / bb_len
          bb_frac_e <- 1 - bb_trim / bb_len
          seg_df <- data.frame(
            x = x1 + (x2 - x1) * bb_frac_s, xend = x1 + (x2 - x1) * bb_frac_e,
            y = y1 + (y2 - y1) * bb_frac_s, yend = y1 + (y2 - y1) * bb_frac_e
          )
        } else {
          seg_df <- data.frame(x = x1, xend = x2, y = y1, yend = y2)
        }
        p <- p + ggplot2::geom_segment(data = seg_df,
          ggplot2::aes(x=.data$x, xend=.data$xend, y=.data$y, yend=.data$yend),
          color="#999999", linewidth=0.35, alpha=0.5,
          lineend="round", inherit.aes=FALSE)
      } else {
        # Non-backbone: L-shaped routing
        stagger <- (i %% 5 - 2) * 0.25
        if (is_straight) {
          edge_pts <- data.frame(x = c(x1, x2), y = c(y1, y2))
        } else {
          is_mostly_vertical <- abs(y2 - y1) >= abs(x2 - x1)
          if (is_mostly_vertical) {
            mid_y <- round((y2 + stagger) * 4) / 4
            route_pts <- data.frame(
              x = c(x1, x1, x2, x2),
              y = c(y1, mid_y, mid_y, y2)
            )
          } else {
            mid_x <- round((x2 + stagger) * 4) / 4
            route_pts <- data.frame(
              x = c(x1, mid_x, mid_x, x2),
              y = c(y1, y1, y2, y2)
            )
          }
          edge_pts <- .dnmb_cct_rounded_route_points(route_pts, radius = glyco_grid_step)
        }
        # Adaptive trim: at least node radius (0.10) so lines don't pass through symbols
        edge_total_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
        edge_trim <- min(0.12, max(0.10, edge_total_len * 0.15))
        # Simple gray line with arrow (unified width)
        p <- p + .dnmb_cct_single_path_layer(
          edge_pts, color = "#AAAAAA",
          linewidth = 0.35, alpha = 0.5,
          arrow_last = TRUE, arrow_length = 0.015,
          trim_start = edge_trim, trim_end = edge_trim
        )
      }
    }
  }

  # Continuity overlays for matched carbon pathways:
  # draw a single continuous route from entry through central metabolism.
  continuity_routes <- .dnmb_cct_continuity_routes()
  continuity_entry_map <- .dnmb_cct_auto_entry_map()
  active_continuity <- intersect(tolower(names(continuity_entry_map)), matched_pathway_ids)
  if (length(active_continuity) > 0) {
    for (pid in active_continuity) {
      best_hit <- route_best_match(pid)
      if (is.null(best_hit)) next
      sugar_type <- tolower(as.character(carbon_src$sugar_type[match(tolower(pid), tolower(carbon_src$id))]))
      cont_col <- route_color_for(path_ids = pid, sugar_type = sugar_type)
      generic_nodes <- .dnmb_cct_best_continuity_nodes(
        pathway_id = pid,
        cyto_nodes = cyto_nodes,
        cyto_edges = cyto_edges,
        route_group_members = route_group_members,
        entry_map = continuity_entry_map,
        matched_steps = matched_steps,
        transporters = transporters
      )
      fallback_nodes <- if (pid %in% names(continuity_routes)) continuity_routes[[pid]] else NULL
      cont_node_ids <- if (!is.null(generic_nodes)) generic_nodes else fallback_nodes
      cont_pts <- .dnmb_cct_route_overlay_points(
        node_ids = cont_node_ids,
        node_x = node_x,
        node_y = node_y,
        grid_step = glyco_grid_step
      )
      if (is.null(cont_pts) || nrow(cont_pts) < 2) next
      frac <- if (!is.null(pstats) && pid %in% pstats$pathway_id) {
        pstats$fraction[pstats$pathway_id == pid][1]
      } else {
        NA_real_
      }
      frac <- ifelse(is.na(frac), 0.5, frac)
      cont_imp <- .dnmb_cct_path_importance(
        confidence = best_hit$confidence,
        fraction = frac,
        sink_weight = if (!is.null(generic_nodes) && !is.null(attr(generic_nodes, "sink_weight"))) attr(generic_nodes, "sink_weight") else 1
      )
      cont_alpha <- 0.25 + 0.15 * cont_imp
      cont_lw <- 0.35  # unified with all other edges
      # Simple colored line (no gradient)
      p <- p + .dnmb_cct_single_path_layer(
        cont_pts, color = cont_col, linewidth = cont_lw, alpha = cont_alpha,
        trim_start = 0.12, trim_end = 0.12
      )
    }
  }

  # Draw TCA cycle arc (smooth circle behind TCA nodes)
  tca_nodes <- cyto_nodes[cyto_nodes$type == "tca", , drop = FALSE]
  if (nrow(tca_nodes) >= 3) {
    tca_cx <- mean(tca_nodes$x)
    tca_cy <- mean(tca_nodes$y)
    tca_r  <- max(sqrt((tca_nodes$x - tca_cx)^2 + (tca_nodes$y - tca_cy)^2))
    # Step 1: Draw perfect circle (gapmind_aa style — thick grey dashed)
    theta_seq <- seq(0, 2 * pi, length.out = 200)
    tca_circle <- data.frame(x = tca_cx + tca_r * cos(theta_seq),
                              y = tca_cy + tca_r * sin(theta_seq))
    p <- p + ggplot2::geom_path(
      data = tca_circle,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "#CCCCCC", linewidth = 1.5, linetype = "solid", inherit.aes = FALSE)
  }

  # Draw SNFG nodes for all cytoplasm node types (symbols only, no labels)
  node_r <- 0.10  # unified size
  for (ntype in c("backbone", "ppp", "entry_intermediate", "tca", "ed", "pyruvate_branch")) {
    nset <- cyto_nodes[cyto_nodes$type == ntype, , drop = FALSE]
    for (i in seq_len(nrow(nset))) {
      nd <- nset[i, , drop = FALSE]
      # Symbol only — full name labels are rendered separately below
      p <- .dnmb_snfg_render_symbol_v2(p, nd$x, nd$y, nd$sugar_type, r = node_r, label = NULL)
    }
  }

  # (backbone/ppp/entry_intermediate labels rendered in final overlay section below)

  # Carbon source nodes — CONVERGE: multiple substrates → 1 hub per sugar_type
  # 1) Group by sugar_type, find hub position (mean x of group)
  # 2) Draw converging arrows from each substrate label → hub SNFG symbol
  # 3) Hub connects onward to entry intermediate
  cs_groups <- split(seq_len(nrow(carbon_src)), carbon_src$sugar_type)
  hidden_secondary_cs_labels <- c("Lactose")
  hub_positions <- list()  # sugar_type -> list(x, y)
  route_label_rows <- list()
  for (st in names(cs_groups)) {
    idx <- cs_groups[[st]]
    hub_x <- if (st %in% names(sugar_lane_x)) {
      unname(sugar_lane_x[st])
    } else {
      .dnmb_cct_snap_to_grid(mean(carbon_src$x[idx]), step = glyco_grid_step)
    }
    hub_y <- .dnmb_cct_snap_to_grid(carbon_src$y[idx[1]], step = glyco_grid_step)
    hub_positions[[st]] <- list(x = hub_x, y = hub_y)
    # Draw hub SNFG symbol (no label — final overlay adds label above route lines)
    p <- .dnmb_snfg_render_symbol_v2(p, hub_x, hub_y, st, r = 0.10, label = NULL)
    # Draw converging lines from each substrate to hub
    for (i in idx) {
      nd <- carbon_src[i, , drop = FALSE]
      if (abs(nd$x - hub_x) > 0.2) {
        # Substrate label at its original x, slightly above
        if (!nd$id %in% hidden_secondary_cs_labels) {
          p <- p + ggplot2::annotate("text", x = nd$x, y = nd$y + 0.20,
            label = nd$label, size = 0.8, color = "#111111", hjust = 0.5, angle = 0)
        }
        # Orthogonal guide segments snapped to the glycolysis grid
        guide_y <- .dnmb_cct_snap_to_grid(hub_y + glyco_grid_step, step = glyco_grid_step)
        guide_col <- route_color_for(path_ids = tolower(nd$id), sugar_type = st)
        guide_path <- .dnmb_cct_rounded_route_points(
          data.frame(
            x = c(nd$x, nd$x, hub_x, hub_x),
            y = c(nd$y + 0.15, guide_y, guide_y, hub_y)
          ),
          radius = glyco_grid_step
        )
        p <- p + .dnmb_cct_single_path_layer(
          guide_path, color = .dnmb_cct_reference_gray("light"), linewidth = 0.16, alpha = 0.26
        )
        if (!is.null(route_best_match(tolower(nd$id)))) {
          for (ly in .dnmb_cct_gradient_path_layers(
            guide_path, color = guide_col, linewidth = 0.18, alpha = 0.40
          )) p <- p + ly
        }
      } else {
        # Single substrate = same as hub, just label below
        p <- p + ggplot2::annotate("text", x = nd$x, y = nd$y - 0.20,
          label = nd$label, size = 0.9, fontface = "bold", color = "#111111", hjust = 0.5)
      }
    }
  }
  carbon_src$confidence_f <- factor(carbon_src$confidence,
                                     levels = c("high", "medium", "low", "none"))

  # ---- HUB → ENTRY INTERMEDIATE connections (vertical lines) ----
  # Each hub connects downward to its entry intermediate using auto_entry_map
  entry_map <- .dnmb_cct_auto_entry_map()
  ei_nodes <- cyto_nodes[cyto_nodes$type == "entry_intermediate", , drop = FALSE]
  bb_nodes_all <- cyto_nodes[cyto_nodes$type %in% c("backbone", "ppp"), , drop = FALSE]
  all_targets <- rbind(ei_nodes, bb_nodes_all)
  target_x <- stats::setNames(all_targets$x, all_targets$id)
  target_y <- stats::setNames(all_targets$y, all_targets$id)
  hub_route_specs <- list()
  for (st in names(hub_positions)) {
    hp <- hub_positions[[st]]
    cs_of_type <- carbon_src$id[carbon_src$sugar_type == st]
    route_key <- vapply(cs_of_type, function(csid) {
      paste(.dnmb_cct_entry_route_nodes(tolower(csid), node_ids = names(target_x)), collapse = "|")
    }, character(1))
    route_groups <- split(cs_of_type, route_key)
    for (rk in names(route_groups)) {
      route_csids <- route_groups[[rk]]
      route_nodes <- unlist(strsplit(rk, "|", fixed = TRUE))
      route_nodes <- route_nodes[nzchar(route_nodes)]
      rep_csid <- route_csids[1]
      entry_target <- if (length(route_nodes) > 0) tail(route_nodes, 1) else unname(entry_map[tolower(rep_csid)])
      if (is.na(entry_target) || !entry_target %in% names(target_x)) next
      best_hit <- route_best_match(tolower(route_csids))
      route_col <- route_color_for(path_ids = tolower(route_csids), sugar_type = st)
      frac_val <- NA_real_
      if (!is.null(pstats)) {
        frac_vec <- pstats$fraction[pstats$pathway_id %in% tolower(route_csids)]
        if (length(frac_vec) > 0) frac_val <- max(frac_vec, na.rm = TRUE)
      }
      route_imp <- if (!is.null(best_hit)) {
        .dnmb_cct_path_importance(best_hit$confidence, frac_val, 1)
      } else {
        .dnmb_cct_path_importance("medium", frac_val, 1)
      }
      tx <- .dnmb_cct_snap_to_grid(unname(target_x[entry_target]), step = glyco_grid_step)
      ty <- .dnmb_cct_snap_to_grid(unname(target_y[entry_target]), step = glyco_grid_step)
      hub_route_specs[[length(hub_route_specs) + 1L]] <- data.frame(
        sugar_type = st,
        cs_id = rep_csid,
        path_ids = paste(sort(unique(tolower(route_csids))), collapse = "|"),
        hub_x = hp$x,
        hub_y = hp$y,
        route_nodes = paste(route_nodes, collapse = "|"),
        target_id = entry_target,
        target_x = tx,
        target_y = ty,
        route_col = route_col,
        route_imp = route_imp,
        best_step_id = if (!is.null(best_hit)) best_hit$step_id else "",
        label_step = if (!is.null(best_hit)) .dnmb_cct_short_step_label(best_hit$step_id, best_hit$pathway_id) else "",
        label_locus = if (!is.null(best_hit)) best_hit$locus_tag else "",
        label_rank = if (!is.null(best_hit)) best_hit$rank else 0,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(hub_route_specs) > 0) {
    hub_route_df <- do.call(rbind, hub_route_specs)
    hub_route_df$target_rank <- 1L
    hub_route_df$target_count <- 1L
    split_target <- split(seq_len(nrow(hub_route_df)), hub_route_df$target_id)
    for (grp in split_target) {
      ord <- grp[order(hub_route_df$hub_x[grp], hub_route_df$hub_y[grp])]
      hub_route_df$target_rank[ord] <- seq_along(ord)
      hub_route_df$target_count[grp] <- length(grp)
    }

    for (i in seq_len(nrow(hub_route_df))) {
      rr <- hub_route_df[i, , drop = FALSE]
      route_nodes <- unlist(strsplit(rr$route_nodes, "|", fixed = TRUE))
      route_nodes <- route_nodes[nzchar(route_nodes)]
      main_path <- .dnmb_cct_points_from_start_to_nodes(
        start_x = rr$hub_x,
        start_y = rr$hub_y,
        node_ids = route_nodes,
        node_x = target_x,
        node_y = target_y,
        grid_step = glyco_grid_step,
        lane_rank = rr$target_rank,
        lane_count = rr$target_count
      )
      if (is.null(main_path)) {
        main_path <- .dnmb_cct_hub_entry_path(
          hub_x = rr$hub_x,
          hub_y = rr$hub_y,
          target_x = rr$target_x,
          target_y = rr$target_y,
          lane_rank = rr$target_rank,
          target_count = rr$target_count,
          grid_step = glyco_grid_step
        )
      }
      # Simple colored line from hub to entry intermediate
      hub_col <- if (rr$label_rank > 0) rr$route_col else "#AAAAAA"
      hub_alpha <- if (rr$label_rank > 0) 0.55 else 0.35
      p <- p + .dnmb_cct_single_path_layer(
        main_path, color = hub_col,
        linewidth = 0.35,
        alpha = hub_alpha,
        arrow_last = TRUE, arrow_length = 0.02,
        trim_end = 0.08
      )
      for (ly in .dnmb_cct_junction_glyph_layers(
        x = rr$target_x, y = rr$target_y, color = .dnmb_cct_reference_gray("dark"),
        size = 1.3 + 0.5 * rr$route_imp, alpha = 0.45
      )) p <- p + ly
      if (rr$label_rank > 0) {
        for (ly in .dnmb_cct_junction_glyph_layers(
          x = rr$target_x, y = rr$target_y, color = rr$route_col,
          size = 1.5 + 0.7 * rr$route_imp, alpha = 0.7
        )) p <- p + ly
      }

      hit_ids <- unlist(strsplit(rr$path_ids, "|", fixed = TRUE))
      route_hits <- matched_steps[tolower(matched_steps$pathway_id) %in% hit_ids, , drop = FALSE]
      if (nrow(route_hits) > 0) {
        keep_tr_like <- unlist(mapply(
          .dnmb_cct_is_transport_like_step,
          step_id = route_hits$step_id,
          gene_name = route_hits$step_id,
          is_pts = FALSE,
          SIMPLIFY = TRUE
        ))
        route_hits <- route_hits[keep_tr_like, , drop = FALSE]
      }
      if (nrow(route_hits) == 0 && nzchar(rr$label_step) && nzchar(rr$label_locus)) {
        route_hits <- data.frame(
          step_id = rr$best_step_id,
          locus_tag = rr$label_locus,
          rank = rr$label_rank,
          stringsAsFactors = FALSE
        )
      }
      if (nrow(route_hits) > 0) {
        route_hits <- route_hits[order(-route_hits$rank, route_hits$step_id, route_hits$locus_tag), , drop = FALSE]
        route_hits <- utils::head(route_hits, 3L)  # max 3 labels per route
        for (hi in seq_len(nrow(route_hits))) {
          # Place label at midpoint of the path (between two metabolites)
          path_idx <- .dnmb_cct_route_label_index(
            main_path,
            frac = 0.50 + 0.12 * (hi - 1L)
          )
          path_mid <- main_path[path_idx, , drop = FALSE]
          route_label_rows[[length(route_label_rows) + 1L]] <- data.frame(
            x = path_mid$x + 0.12,
            y = path_mid$y,
            label = paste0(.dnmb_cct_short_step_label(route_hits$step_id[hi], rr$cs_id), "\n", route_hits$locus_tag[hi]),
            step_id = route_hits$step_id[hi],
            locus_tag = route_hits$locus_tag[hi],
            color = rr$route_col,
            target_id = rr$target_id,
            target_rank = rr$target_rank,
            priority = 10 * rr$route_imp + route_hits$rank[hi] + 0.01 * (nrow(route_hits) - hi),
            angle = 0,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  shown_route_keys <- character(0)
  if (length(route_label_rows) > 0) {
    shown_route_keys <- unique(vapply(route_label_rows, function(x) {
      paste(x$step_id[1], x$locus_tag[1], sep = "::")
    }, character(1)))
  }
  edge_label_count <- integer(if (!is.null(cyto_edges)) nrow(cyto_edges) else 0L)
  if (!is.null(cyto_edges) && nrow(cyto_edges) > 0 && nrow(matched_steps) > 0) {
    # Only label high/medium confidence steps to reduce clutter
    ms_ord <- matched_steps[matched_steps$rank >= 2L, , drop = FALSE]
    ms_ord <- ms_ord[order(-ms_ord$rank, ms_ord$pathway_id, ms_ord$step_id), , drop = FALSE]
    for (mi in seq_len(nrow(ms_ord))) {
      ms <- ms_ord[mi, , drop = FALSE]
      key <- paste(ms$step_id, ms$locus_tag, sep = "::")
      if (key %in% shown_route_keys) next
      # Also skip if this locus_tag was already shown (different step, same gene)
      if (ms$locus_tag %in% sub("^.*::", "", shown_route_keys)) next

      tgt_nodes <- .dnmb_cct_step_target_nodes(ms$step_id, ms$pathway_id)
      hint_nodes <- unique(c(
        tgt_nodes,
        .dnmb_cct_pathway_hint_nodes(
          path_ids = ms$pathway_id,
          matched_steps = ms,
          transporters = transporters,
          entry_map = entry_map
        )
      ))
      cand_idx <- integer(0)
      if (length(hint_nodes) > 0) {
        cand_idx <- which(cyto_edges$from %in% hint_nodes | cyto_edges$to %in% hint_nodes)
      }
      if (length(cand_idx) == 0) {
        cand_idx <- which(cyto_edges$pathway %in% c(
          .dnmb_cct_route_group_for_pathway(ms$pathway_id, route_group_members),
          "ppp", "ed", "backbone"
        ))
      }
      if (length(cand_idx) == 0) next

      score_edge <- function(ii) {
        ce <- cyto_edges[ii, , drop = FALSE]
        score <- 0
        if (ce$to %in% tgt_nodes) score <- score + 3
        if (ce$from %in% tgt_nodes) score <- score + 2
        if (ce$to %in% hint_nodes || ce$from %in% hint_nodes) score <- score + 1
        if (ce$pathway == .dnmb_cct_route_group_for_pathway(ms$pathway_id, route_group_members)) score <- score + 0.8
        score - 0.05 * edge_label_count[ii]
      }
      edge_scores <- vapply(cand_idx, score_edge, numeric(1))
      best_idx <- cand_idx[which.max(edge_scores)]
      edge_label_count[best_idx] <- edge_label_count[best_idx] + 1L
      if (edge_label_count[best_idx] > 2L) next  # max 2 labels per edge
      ce <- cyto_edges[best_idx, , drop = FALSE]
      edge_pts <- .dnmb_cct_edge_points_from_row(ce, edge_idx = best_idx, grid_step = glyco_grid_step)
      if (is.null(edge_pts) || nrow(edge_pts) == 0) next
      tgt_node <- if (length(tgt_nodes) > 0) tgt_nodes[1] else ce$to
      edge_sugar <- unname(node_sugar_type[tgt_node])
      if (is.na(edge_sugar) || !nzchar(edge_sugar) || edge_sugar == "intermediate") {
        edge_sugar <- unname(node_sugar_type[ce$to])
      }
      # Place label at midpoint of the edge (between two metabolites)
      path_idx <- .dnmb_cct_route_label_index(
        edge_pts,
        frac = 0.50 + 0.15 * (edge_label_count[best_idx] - 1L)
      )
      path_mid <- edge_pts[path_idx, , drop = FALSE]
      route_label_rows[[length(route_label_rows) + 1L]] <- data.frame(
        x = path_mid$x + 0.12,
        y = path_mid$y,
        label = paste0(.dnmb_cct_short_step_label(ms$step_id, ms$pathway_id), "\n", ms$locus_tag),
        step_id = ms$step_id,
        locus_tag = ms$locus_tag,
        color = route_color_for(path_ids = ms$pathway_id, sugar_type = edge_sugar, fallback = "#8C8C8C"),
        target_id = tgt_node,
        target_rank = edge_label_count[best_idx],
        priority = 8 + ms$rank,
        angle = 0,
        stringsAsFactors = FALSE
      )
      shown_route_keys <- c(shown_route_keys, key)
    }
  }

  # Route labels (gene/locus_tag) are deferred to render LAST — see below after metabolite labels
  if (length(route_label_rows) > 0) {
    route_label_df <- do.call(rbind, route_label_rows)
    route_label_df <- route_label_df[order(-route_label_df$priority), , drop = FALSE]
    route_label_df <- route_label_df[!duplicated(route_label_df$locus_tag), , drop = FALSE]
    center_ref <- stats::median(backbone_nodes$x, na.rm = TRUE)
    route_label_df$priority <- route_label_df$priority -
      2.1 * pmax(0, 1.6 - abs(route_label_df$x - center_ref))
    route_label_df <- .dnmb_cct_layout_route_labels(route_label_df, center_x = center_ref)
  } else {
    route_label_df <- NULL
  }

  # ---- Final extracellular symbol overlay: keep polymer/branch nodes above route lines ----
  for (ri in seq_len(n_extra)) {
    sub <- extra_substrates[ri, , drop = FALSE]
    ry <- extra_y_map[sub$id]
    chain_x <- extra_x_map[sub$id]
    chain <- extra_chains[[sub$id]]
    if (length(chain) > 0) {
      chain_xs <- chain_x + (seq_along(chain) - 1) * chain_step
      for (ci in seq_along(chain)) {
        sug <- chain[[ci]]
        abbr <- .dnmb_snfg_abbreviation_v2(sug)
        p <- .dnmb_snfg_render_symbol_v2(p, chain_xs[ci], ry, sug, r = sym_size, label = abbr)
      }
    }
    branches <- extra_branches[[sub$id]]
    if (!is.null(branches)) {
      for (bi in seq_along(branches)) {
        br <- branches[[bi]]
        br_x <- .dnmb_cct_extra_chain_node_x(
          chain_x = chain_x,
          pos = br$pos,
          sym_size = sym_size,
          gap = chain_gap
        )
        br_y <- ry - sym_size * (2.5 + 1.4 * (bi - 1))
        abbr_br <- .dnmb_snfg_abbreviation_v2(br$sugar)
        p <- .dnmb_snfg_render_symbol_v2(p, br_x, br_y, br$sugar, r = sym_size, label = abbr_br)
      }
    }
  }

  # ---- Final metabolite overlay: keep nodes/labels above every route line ----
  for (mi in seq_along(mono_ids)) {
    mid <- mono_ids[mi]
    mx <- mono_xs[mi]
    abbr <- .dnmb_snfg_abbreviation_v2(mid)
    p <- .dnmb_snfg_render_symbol_v2(p, mx, y_mono, mid, r = 0.10, label = abbr)
  }

  for (st in names(hub_positions)) {
    hp <- hub_positions[[st]]
    abbr <- .dnmb_snfg_abbreviation_v2(st)
    p <- .dnmb_snfg_render_symbol_v2(p, hp$x, hp$y, st, r = 0.10, label = abbr)
  }

  for (ntype in c("backbone", "ppp", "entry_intermediate", "tca", "ed", "pyruvate_branch")) {
    nset <- cyto_nodes[cyto_nodes$type == ntype, , drop = FALSE]
    for (i in seq_len(nrow(nset))) {
      nd <- nset[i, , drop = FALSE]
      # Symbol only — full name labels are rendered below
      p <- .dnmb_snfg_render_symbol_v2(p, nd$x, nd$y, nd$sugar_type, r = node_r, label = NULL)
    }
  }

  # ---- Metabolite labels: unified bold style, rendered LAST (above all symbols) ----
  all_metab_labels <- data.frame(
    x = numeric(0), y = numeric(0), label = character(0),
    nudge_x = numeric(0), nudge_y = numeric(0),
    stringsAsFactors = FALSE
  )
  # All cytoplasmic metabolites: unified style
  for (ntype in c("backbone", "ppp", "entry_intermediate", "ed", "pyruvate_branch")) {
    nset <- cyto_nodes[cyto_nodes$type == ntype, , drop = FALSE]
    if (nrow(nset) > 0) {
      all_metab_labels <- rbind(all_metab_labels, data.frame(
        x = nset$x, y = nset$y, label = nset$label,
        nudge_x = 0, nudge_y = 0.16,
        stringsAsFactors = FALSE
      ))
    }
  }
  # TCA: radial nudge outward from cycle center
  if (nrow(tca_nodes) > 0) {
    tca_cx_lbl <- mean(tca_nodes$x); tca_cy_lbl <- mean(tca_nodes$y)
    tca_nodes$radial_dx <- tca_nodes$x - tca_cx_lbl
    tca_nodes$radial_dy <- tca_nodes$y - tca_cy_lbl
    tca_nodes$radial_dist <- sqrt(tca_nodes$radial_dx^2 + tca_nodes$radial_dy^2)
    tca_ring <- tca_nodes[tca_nodes$radial_dist > 0.25, , drop = FALSE]
    if (nrow(tca_ring) > 0) {
      tca_ring$ux <- tca_ring$radial_dx / tca_ring$radial_dist
      tca_ring$uy <- tca_ring$radial_dy / tca_ring$radial_dist
      all_metab_labels <- rbind(all_metab_labels, data.frame(
        x = tca_ring$x, y = tca_ring$y, label = tca_ring$label,
        nudge_x = tca_ring$ux * 0.20, nudge_y = tca_ring$uy * 0.20,
        stringsAsFactors = FALSE
      ))
    }
  }
  # Single unified ggrepel layer — bold, size 1.3, dark gray
  if (nrow(all_metab_labels) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = all_metab_labels,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 1.3, fontface = "bold", color = "#333333",
      nudge_x = all_metab_labels$nudge_x,
      nudge_y = all_metab_labels$nudge_y,
      segment.size = 0.12, segment.color = "#CCCCCC",
      segment.alpha = 0.3,
      box.padding = 0.06, point.padding = 0.04,
      min.segment.length = 0.12,
      max.overlaps = 60,
      force = 1.2, force_pull = 2.5,
      max.iter = 5000,
      seed = 42,
      inherit.aes = FALSE
    )
  }

  # ---- Route labels (gene/locus_tag): rendered LAST, above all symbols and metabolite labels ----
  if (!is.null(route_label_df) && nrow(route_label_df) > 0) {
    # Thin connector segments from label to anchor point
    conn_df <- route_label_df[abs(route_label_df$x_lab - route_label_df$x) > 0.08 |
                               abs(route_label_df$y_lab - route_label_df$y) > 0.08, , drop = FALSE]
    if (nrow(conn_df) > 0) {
      p <- p + ggplot2::geom_segment(
        data = conn_df,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$x_lab, yend = .data$y_lab),
        linewidth = 0.10, color = "#BBBBBB", alpha = 0.4,
        inherit.aes = FALSE
      )
    }
    p <- p + ggplot2::geom_text(
      data = route_label_df,
      ggplot2::aes(x = .data$x_lab, y = .data$y_lab, label = .data$label),
      size = 0.72, color = "#111111",
      fontface = "bold", hjust = route_label_df$hjust,
      angle = 0,
      lineheight = 0.85, inherit.aes = FALSE
    )
  }

  # ====================================================================
  # LEGENDS
  # ====================================================================
  p <- p + ggplot2::scale_color_manual(
    name = "Pathway\nCompleteness",
    values = c(high = "#2CA25F", medium = "#FEC44F", low = "#F03B20", none = "#CCCCCC"),
    labels = c(high = "High (\u226575%)", medium = "Medium (50-74%)",
               low = "Low (<50%)", none = "Not found"),
    drop = FALSE)

  # SNFG legend — compact, RIGHT side, aligned with TCA circle height
  snfg_leg <- data.frame(
    sugar = c("glucose", "galactose", "mannose", "fructose", "glcnac",
              "fucose", "rhamnose", "xylose", "arabinose", "glca", "gala", "ribose"),
    label = c("Glc", "Gal", "Man", "Fru", "GlcNAc", "Fuc", "Rha", "Xyl",
              "Ara", "GlcA", "GalA", "Rib"),
    stringsAsFactors = FALSE, row.names = NULL)
  # Right side x: just past the rightmost cytoplasm node
  leg_x <- max(cyto_nodes$x) + 1.5
  # y aligned with TCA circle center
  tca_cy_leg <- if (nrow(tca_nodes) > 0) mean(tca_nodes$y) else 0.5
  sp <- 0.20  # tight row spacing
  n_snfg <- nrow(snfg_leg)
  leg_y_top <- tca_cy_leg + (n_snfg / 2) * sp
  p <- p + ggplot2::annotate("text", x = leg_x, y = leg_y_top + 0.3,
    hjust = 0, size = 1.4, fontface = "bold", color = "#333333",
    label = "SNFG Symbols")
  for (i in seq_len(n_snfg)) {
    ly <- leg_y_top - (i - 1) * sp
    p <- .dnmb_snfg_render_symbol_v2(p, leg_x, ly, snfg_leg$sugar[i], r = 0.04)
    p <- p + ggplot2::annotate("text", x = leg_x + 0.16, y = ly,
      label = snfg_leg$label[i], hjust = 0, size = 1.1, color = "#555555")
  }
  # Bond legend below SNFG
  bleg_y <- leg_y_top - n_snfg * sp - 0.15
  p <- p +
    ggplot2::geom_segment(
      data = data.frame(x = leg_x, xend = leg_x + 0.3, y = bleg_y, yend = bleg_y),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 0.5, color = "#444444", linetype = "solid", inherit.aes = FALSE) +
    ggplot2::annotate("text", x = leg_x + 0.4, y = bleg_y,
      label = "\u03b1 bond", hjust = 0, size = 1.0, color = "#555555")
  p <- p +
    ggplot2::geom_segment(
      data = data.frame(x = leg_x, xend = leg_x + 0.3, y = bleg_y - sp, yend = bleg_y - sp),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 0.5, color = "#444444", linetype = "dashed", inherit.aes = FALSE) +
    ggplot2::annotate("text", x = leg_x + 0.4, y = bleg_y - sp,
      label = "\u03b2 bond", hjust = 0, size = 1.0, color = "#555555")
  # Scissors
  sc_ly <- .dnmb_scissors_grob_v2(leg_x + 0.10, bleg_y - 2 * sp, size = 0.05)
  for (l in sc_ly) p <- p + l
  p <- p + ggplot2::annotate("text", x = leg_x + 0.4, y = bleg_y - 2 * sp,
    label = "GH cleavage", hjust = 0, size = 1.0, color = "#555555")
  enz_leg_y <- bleg_y - 3 * sp

  # ====================================================================
  # THEME & OUTPUT
  # ====================================================================
  y_bottom <- y_bot - 0.5
  y_top_plot <- max(y_extra_top + 1.5, max(cyto_nodes$y) + 1.5)
  x_left <- mem_xmin - 0.05
  x_right <- max(mem_xmax + 0.10, leg_x + 2.0)  # include right-side legend
  x_span <- max(1, x_right - x_left)
  y_span <- max(1, y_top_plot - y_bottom)
  n_active <- length(cs_ids_ordered)
  # Derive dimensions from coordinate span to match coord_fixed(ratio=1)
  scale_factor <- 0.78
  plot_height <- max(15.5, min(40, y_span * scale_factor))
  min_width <- max(8, 5 + n_active * 0.55)
  plot_width <- max(min_width, min(40, x_span * scale_factor + 0.9))

  p <- p +
    ggplot2::coord_fixed(ratio = 1,
      xlim = c(x_left, x_right), ylim = c(y_bottom, y_top_plot), clip = "off") +
    ggplot2::labs(
      title = "CAZy Carbon Transport Map (3-Zone)",
      subtitle = paste0(
        "Extracellular: SNFG glycan chains + GH cleavage | ",
        "Membrane: transporters | Cytoplasm: grid metabolism")) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0,
                                          margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = 7.5, hjust = 0, color = "#888888",
                                             margin = ggplot2::margin(b = 8)),
      legend.position = c(0.85, 0.12),
      legend.justification = c(0, 0),
      legend.direction = "vertical",
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 7, face = "bold"),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.width = grid::unit(12, "pt"),
      legend.key.height = grid::unit(5, "pt"),
      plot.margin = ggplot2::margin(14, 12, 16, 12))

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(file_stub, ".pdf"))
  .dnmb_module_plot_save(p, pdf_path, width = plot_width, height = plot_height)
  list(pdf = pdf_path)
}

 # ---- Comprehensive metabolism node position table ----
 # Glycolysis backbone + PPP branch + all sugar/acid entry points + SNFG annotation
# Layout: compact vertical backbone, sugars fan out left, PPP/pentoses right, acids below

# ====================================================================
# AUTO-LAYOUT ALGORITHM: Compute cytoplasm node positions from pathway data
# This replaces hardcoded positions with a data-driven layout
# ====================================================================
.dnmb_cct_auto_entry_map <- function() {
  # Each carbon source → which backbone node it enters
  # This is universal biochemistry, not genome-specific
  out <- c(
    # → Glucose → Glc-6-P
    glucose = "Glc-6-P", maltose = "Glc-6-P", cellobiose = "Glc-6-P",
    trehalose = "Glc-6-P",
    # → Galactose → Gal-1-P → Glc-1-P → Glc-6-P
    galactose = "Glc-6-P", lactose = "Glc-6-P",
    # → Fructose → Fru-6-P or Fru-1-P → Fru-1,6-BP
    fructose = "Fru-6-P", sucrose = "Fru-6-P",
    # → Mannose → Man-6-P → Fru-6-P
    mannose = "Fru-6-P", mannitol = "Fru-6-P",
    # → GlcNAc → GlcNAc-6-P → GlcN-6-P → Fru-6-P
    NAG = "Fru-6-P", glucosamine = "Fru-6-P",
    # → Xylose → Xu-5-P (PPP)
    xylose = "Xu-5-P", arabinose = "Xu-5-P",
    # → Ribose → R-5-P (PPP)
    ribose = "R-5-P",
    # → Gluconate → 6-PG (PPP)
    gluconate = "6-PG",
    # → Glycerol → DHAP
    glycerol = "DHAP",
    # → Fucose → L-fuculose-1-P → DHAP + L-lactaldehyde
    fucose = "DHAP",
    # → Rhamnose → L-rhamnulose-1-P → DHAP + L-lactaldehyde
    rhamnose = "DHAP"
  )
  stats::setNames(unname(out), tolower(names(out)))
}

# Intermediate steps between carbon source and backbone entry point
# Returns list of (intermediate_id, sugar_type, enzyme_gene) tuples
.dnmb_cct_auto_intermediates <- function() {
  out <- list(
    glucose     = list(),  # direct: Glc → glk → G6P
    maltose     = list(c("Glucose", "glucose", "maltase/GH1")),
    cellobiose  = list(c("Glucose", "glucose", "GH1/GH4")),
    trehalose   = list(c("Trehalose-6-P", "glucose", "treP"), c("Glucose", "glucose", "treC")),
    galactose   = list(c("Gal-1-P", "galactose", "galK"), c("UDP-Gal", "galactose", "galT"), c("UDP-Glc", "glucose", "galE"), c("Glc-1-P", "glucose", "pgmA")),
    lactose     = list(c("Glucose", "glucose", "lacZ"), c("Gal-1-P", "galactose", "galK"), c("UDP-Gal", "galactose", "galT"), c("UDP-Glc", "glucose", "galE"), c("Glc-1-P", "glucose", "pgmA")),
    fructose    = list(c("Fru-1-P", "fructose", "fruK")),
    sucrose     = list(c("Glucose+Fructose", "fructose", "GH32/scrK")),
    mannose     = list(c("Man-6-P", "mannose", "manA")),
    mannitol    = list(c("Mannitol-1-P", "mannose", "mtlA"), c("Fru-6-P", "fructose", "mtlD")),
    NAG         = list(c("GlcNAc-6-P", "glcnac", "nagE"), c("GlcN-6-P", "glucosamine", "nagA")),
    glucosamine = list(c("GlcN-6-P", "glucosamine", "nagB")),
    raffinose   = list(c("Galactose", "galactose", "GH27/GH36"), c("Sucrose", "sucrose", "sucrose")),
    xylose      = list(c("Xylulose", "xylose", "xylA"), c("Xu-5-P", "xylose", "xylB")),
    arabinose   = list(c("Ribulose", "arabinose", "araA"), c("Ribulose-5-P", "arabinose", "araB"), c("Xu-5-P", "xylose", "araD")),
    ribose      = list(c("R-5-P", "ribose", "rbsK")),
    gluconate   = list(c("6-PG", "gluconate", "gntK")),
    glycerol    = list(c("Glycerol-3-P", "glycerol", "glpK")),
    fucose      = list(c("Fuculose", "fucose", "fucI"), c("Fuculose-1-P", "fucose", "fucK")),
    rhamnose    = list(c("Rhamnulose", "rhamnose", "rhaA"), c("Rhamnulose-1-P", "rhamnose", "rhaB"))
  )
  names(out) <- tolower(names(out))
  out
}

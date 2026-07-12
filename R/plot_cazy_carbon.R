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
    Deoxyribose = "deoxyribose", Fucose = "fucose",
    Rhamnose = "rhamnose", Gluconate = "glca"
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
# SNFG symbol rendering (native circles + polygons, no ggstar dependency)
# ---------------------------------------------------------------------------

# One residue has one map footprint in every zone. Composite carbohydrates
# grow horizontally by residue count; their individual units are not shrunk.
.dnmb_cct_sugar_icon_radius <- function() 0.10

#' Return polygon vertex coordinates for an SNFG shape
#' @param shape_type One of "circle", "square", "diamond", "triangle", "star"
#' @param x,y Centre coordinates
#' @param size Radius / half-width of the symbol (plot units)
#' @return data.frame with columns \code{x}, \code{y}
#' @keywords internal
.dnmb_snfg_polygon_coords_v2 <- function(
    shape_type, x = 0, y = 0,
    size = .dnmb_cct_sugar_icon_radius()) {
  normalize_bbox <- function(unit_x, unit_y) {
    x_range <- range(unit_x)
    y_range <- range(unit_y)
    scale <- 2 * size / max(diff(x_range), diff(y_range))
    data.frame(
      x = x + (unit_x - mean(x_range)) * scale,
      y = y + (unit_y - mean(y_range)) * scale
    )
  }

  switch(shape_type,
    circle = {
      # 20-gon approximation
      theta <- seq(0, 2 * pi, length.out = 21L)[-21L]
      normalize_bbox(cos(theta), sin(theta))
    },
    square = {
      # axis-aligned square
      normalize_bbox(c(-1, 1, 1, -1), c(-1, -1, 1, 1))
    },
    diamond = {
      # 45-degree rotated square
      normalize_bbox(c(0, 1, 0, -1), c(1, 0, -1, 0))
    },
    triangle = {
      # Equilateral, pointing up. Bbox normalization preserves the shape while
      # giving it the same maximum footprint as circles and squares.
      angles <- c(pi / 2, pi / 2 + 2 * pi / 3, pi / 2 + 4 * pi / 3)
      normalize_bbox(cos(angles), sin(angles))
    },
    star = {
      # 5-pointed star: 10 vertices, alternating outer / inner radius
      angles  <- seq(pi / 2, pi / 2 + 2 * pi, length.out = 11L)[-11L]
      radii   <- rep(c(1, 0.38), 5L)
      normalize_bbox(radii * cos(angles), radii * sin(angles))
    },
    pentagon = {
      # regular pentagon, point up
      angles <- seq(pi / 2, pi / 2 + 2 * pi, length.out = 6L)[-6L]
      normalize_bbox(cos(angles), sin(angles))
    },
    # fallback: small circle
    {
      theta <- seq(0, 2 * pi, length.out = 21L)[-21L]
      normalize_bbox(cos(theta), sin(theta))
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
    deoxyribose  = list(shape = "star",     fill = "#D7A1C4", color = "#8A5879"),
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

.dnmb_native_circle_layer_v2 <- function(x, y, radius, fill, color = NA,
                                         linewidth = 0, alpha = 1) {
  radius <- max(0, as.numeric(radius)[1])
  ggplot2::annotation_custom(
    grob = grid::circleGrob(
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.5, "npc"),
      r = grid::unit(0.5, "npc"),
      gp = grid::gpar(
        fill = fill,
        col = color,
        lwd = as.numeric(linewidth)[1] * 72.27 / 25.4,
        alpha = as.numeric(alpha)[1]
      )
    ),
    xmin = x - radius,
    xmax = x + radius,
    ymin = y - radius,
    ymax = y + radius
  )
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
.dnmb_snfg_symbol_layers_v2 <- function(
    x, y, sugar_type, size = .dnmb_cct_sugar_icon_radius(),
    label = NULL, alpha = 1) {

  sp <- .dnmb_snfg_sugar_spec_v2(sugar_type)

  if (identical(sp$shape, "circle")) {
    layers <- list(
      .dnmb_native_circle_layer_v2(
        x, y, radius = size * 1.18,
        fill = "#FFFFFF", color = NA, linewidth = 0, alpha = 0.92 * alpha
      ),
      .dnmb_native_circle_layer_v2(
        x, y, radius = size,
        fill = sp$fill, color = sp$color, linewidth = 0.30, alpha = alpha
      )
    )
  } else {
    verts <- .dnmb_snfg_polygon_coords_v2(sp$shape, x, y, size)
    halo_verts <- .dnmb_snfg_polygon_coords_v2(sp$shape, x, y, size * 1.18)
    layers <- list(
      ggplot2::geom_polygon(
        data = halo_verts,
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
        fill = "#FFFFFF", color = NA, alpha = 0.92 * alpha,
        inherit.aes = FALSE
      )
    )
    shape_layers <- if (isTRUE(sp$half_fill) && identical(sp$shape, "diamond")) {
      .dnmb_snfg_split_diamond_layers(
        cx = x, cy = y, r = size,
        fill_color = sp$fill, border_color = sp$color, border_lw = 0.3,
        alpha = alpha
      )
    } else {
      list(
        ggplot2::geom_polygon(
          data = verts,
          mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
          fill = sp$fill, color = sp$color, linewidth = 0.3, alpha = alpha,
          inherit.aes = FALSE
        )
      )
    }
    layers <- c(layers, shape_layers)
  }

  # Abbreviation label inside the symbol
  if (!is.null(label) && nzchar(label)) {
    # White text for dark fills, black for light fills
    txt_col <- if (mean(grDevices::col2rgb(sp$fill)) < 180) "white" else "#333333"
    label_lines <- strsplit(as.character(label), "\n", fixed = TRUE)[[1]]
    max_chars <- max(nchar(label_lines))
    text_scale <- if (max_chars <= 3L) {
      11.0
    } else if (max_chars <= 5L) {
      9.5
    } else if (max_chars <= 7L) {
      7.8
    } else {
      6.6
    }
    lbl_df <- data.frame(x = x, y = y, label = label)
    layers[[length(layers) + 1L]] <- ggplot2::geom_text(
      data = lbl_df,
      mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]),
      size = size * text_scale, fontface = "bold", color = txt_col, alpha = alpha,
      lineheight = 0.82,
      inherit.aes = FALSE)
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
.dnmb_snfg_render_symbol_v2 <- function(
    p, cx, cy, sugar_name, r = .dnmb_cct_sugar_icon_radius(),
    label = NULL, alpha = 1) {
  layers <- .dnmb_snfg_symbol_layers_v2(
    cx, cy, sugar_name, size = r, label = label, alpha = alpha
  )
  for (ly in layers) p <- p + ly
  p
}

# Cytoplasmic carbon sources need their own visible glyphs even when several
# substrates converge on the same downstream sugar-family hub. Disaccharides
# are represented by their two constituent SNFG symbols; the remaining carbon
# sources use one substrate-specific symbol.
.dnmb_cct_carbon_source_render_catalog <- function() {
  list(
    Maltose = list(components = c("glucose", "glucose"), bond = "alpha"),
    Cellobiose = list(components = c("glucose", "glucose"), bond = "beta"),
    Galactose = list(components = "galactose", bond = NA_character_),
    Trehalose = list(components = c("glucose", "glucose"), bond = "alpha"),
    Lactose = list(components = c("galactose", "glucose"), bond = "beta"),
    Mannose = list(components = "mannose", bond = NA_character_),
    NAG = list(components = "glcnac", bond = NA_character_),
    Glucosamine = list(components = "glucosamine", bond = NA_character_),
    Fructose = list(components = "fructose", bond = NA_character_),
    Sucrose = list(components = c("glucose", "fructose"), bond = "alpha"),
    Mannitol = list(components = "mannitol", bond = NA_character_),
    Glycerol = list(components = "glycerol", bond = NA_character_),
    Xylose = list(components = "xylose", bond = NA_character_),
    Arabinose = list(components = "arabinose", bond = NA_character_),
    Ribose = list(components = "ribose", bond = NA_character_),
    Deoxyribose = list(components = "deoxyribose", bond = NA_character_),
    Fucose = list(components = "fucose", bond = NA_character_),
    Rhamnose = list(components = "rhamnose", bond = NA_character_),
    Gluconate = list(components = "gluconate", bond = NA_character_)
  )
}

.dnmb_cct_carbon_source_render_spec <- function(source_id) {
  source_id <- trimws(as.character(source_id)[1])
  catalog <- .dnmb_cct_carbon_source_render_catalog()
  catalog_names <- names(catalog)
  catalog_idx <- match(tolower(source_id), tolower(catalog_names))
  registered <- !is.na(catalog_idx)
  if (registered) {
    canonical_id <- catalog_names[catalog_idx]
    item <- catalog[[catalog_idx]]
  } else {
    canonical_id <- source_id
    item <- list(components = tolower(source_id), bond = NA_character_)
  }
  n_comp <- length(item$components)
  data.frame(
    source_id = rep(canonical_id, n_comp),
    component_index = seq_len(n_comp),
    sugar_type = as.character(item$components),
    x_offset_r = seq(
      -(n_comp - 1) / 2, (n_comp - 1) / 2,
      length.out = n_comp
    ) * 2.18,
    bond = rep(as.character(item$bond)[1], n_comp),
    registered = rep(registered, n_comp),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.dnmb_cct_render_carbon_source_node <- function(p, x, y, source_id,
                                                 r = .dnmb_cct_sugar_icon_radius(),
                                                 alpha = 1) {
  spec <- .dnmb_cct_carbon_source_render_spec(source_id)
  component_r <- r
  component_x <- x + spec$x_offset_r * r
  source_abbr <- .dnmb_cct_carbon_source_abbreviation(source_id)

  if (nrow(spec) > 1L) {
    bond_style <- if (identical(spec$bond[1], "beta")) "dashed" else "solid"
    bond_df <- data.frame(
      x = head(component_x, -1L), xend = tail(component_x, -1L),
      y = rep(y, nrow(spec) - 1L), yend = rep(y, nrow(spec) - 1L)
    )
    p <- p +
      ggplot2::geom_segment(
        data = bond_df,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        color = "white", linewidth = 1.25, alpha = 0.92 * alpha,
        lineend = "round", inherit.aes = FALSE
      ) +
      ggplot2::geom_segment(
        data = bond_df,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        color = "#4F5963", linewidth = 0.42, linetype = bond_style,
        alpha = alpha, lineend = "round", inherit.aes = FALSE
      )
  }

  for (i in seq_len(nrow(spec))) {
    component_label <- if (nrow(spec) == 1L) {
      .dnmb_cct_compact_inside_label(source_abbr, max_chars = 5L)
    } else {
      NULL
    }
    p <- .dnmb_snfg_render_symbol_v2(
      p, component_x[i], y, spec$sugar_type[i], r = component_r,
      label = component_label, alpha = alpha
    )
  }

  if (nrow(spec) > 1L && nzchar(source_abbr)) {
    label_df <- data.frame(x = x, y = y, label = source_abbr)
    p <- p +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = r * 12.0, fontface = "bold", color = "white",
        alpha = alpha, lineheight = 0.82, inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = r * 9.2, fontface = "bold", color = "#263238",
        alpha = alpha, lineheight = 0.82, inherit.aes = FALSE
      )
  }
  p
}

# Sugar abbreviation lookup
.dnmb_snfg_abbreviation_v2 <- function(sugar_type) {
  abbr <- c(
    glucose = "Glc", galactose = "Gal", mannose = "Man", fructose = "Fru",
    glcnac = "GlcNAc", galnac = "GalNAc", fucose = "Fuc", rhamnose = "Rha",
    xylose = "Xyl", arabinose = "Ara", glca = "GlcA", gala = "GalA",
    ribose = "Rib", deoxyribose = "dRib", glycerol = "Gro",
    glucosamine = "GlcN", nag = "GlcNAc",
    trehalose = "Tre", sucrose = "Suc", maltose = "Mal",
    cellobiose = "Cel", lactose = "Lac", mannitol = "Mtl",
    gluconate = "Gnt", intermediate = "", organic_acid = "",
    phospho_sugar = "", generic = "")
  key <- tolower(sugar_type)
  if (key %in% names(abbr)) unname(abbr[key]) else ""
}

.dnmb_cct_carbon_source_abbreviation <- function(source_id) {
  source_id <- tolower(trimws(as.character(source_id)[1]))
  abbr <- c(
    maltose = "Mal", cellobiose = "Cel", galactose = "Gal",
    trehalose = "Tre", lactose = "Lac", mannose = "Man",
    nag = "GlcNAc", glucosamine = "GlcN", fructose = "Fru",
    sucrose = "Suc", mannitol = "Mtl", glycerol = "Gro",
    xylose = "Xyl", arabinose = "Ara", ribose = "Rib",
    deoxyribose = "dRib", fucose = "Fuc", rhamnose = "Rha",
    gluconate = "Gnt"
  )
  if (source_id %in% names(abbr)) unname(abbr[source_id]) else ""
}

.dnmb_cct_compact_inside_label <- function(label, max_chars = 7L) {
  label <- trimws(as.character(label)[1])
  if (is.na(label) || !nzchar(label)) return("")
  if (max(nchar(strsplit(label, "\n", fixed = TRUE)[[1]])) <= max_chars) {
    return(label)
  }
  if (grepl("-", label, fixed = TRUE)) {
    return(sub("^([^-]+)-", "\\1-\n", label))
  }
  if (identical(label, "GlcNAc")) return("Glc\nNAc")
  if (identical(label, "GalNAc")) return("Gal\nNAc")
  label
}

# Cytoplasmic sugar and sugar-phosphate nodes carry their abbreviation inside
# the glyph. Non-sugar endpoints retain their external labels.
.dnmb_cct_cytoplasmic_inside_label <- function(id, label, type, sugar_type) {
  id <- as.character(id)[1]
  type <- as.character(type)[1]
  if (!type %in% c("backbone", "ppp", "entry_intermediate", "ed")) return("")
  if (id %in% c("Pyruvate", "Acetyl-CoA")) return("")

  short_by_id <- c(
    Glucose = "Glc",
    Xylulose = "Xul",
    Ribulose = "Rul",
    `Deoxyribose-5-P` = "dRib-5-P",
    Fuculose = "Fucol",
    Rhamnulose = "Rhul",
    `Glycerol-3-P` = "Gro-3-P"
  )
  short <- if (id %in% names(short_by_id)) {
    unname(short_by_id[id])
  } else {
    as.character(label)[1]
  }
  .dnmb_cct_compact_inside_label(short, max_chars = 7L)
}

.dnmb_cct_snfg_legend_catalog <- function() {
  data.frame(
    sugar = c(
      "glucose", "galactose", "mannose", "fructose", "glcnac", "galnac",
      "glucosamine", "fucose", "rhamnose", "xylose", "arabinose", "ribose",
      "deoxyribose", "glca", "gala", "glycerol", "mannitol", "gluconate",
      "trehalose", "sucrose", "maltose", "cellobiose", "lactose"
    ),
    label = c(
      "D-Glucose (Glc)", "D-Galactose (Gal)", "D-Mannose (Man)",
      "D-Fructose (Fru)", "N-Acetyl-D-glucosamine (GlcNAc)",
      "N-Acetyl-D-galactosamine (GalNAc)", "D-Glucosamine (GlcN)",
      "L-Fucose (Fuc)", "L-Rhamnose (Rha)", "D-Xylose (Xyl)",
      "L-Arabinose (Ara)", "D-Ribose (Rib)",
      "2-Deoxy-D-ribose (dRib)", "D-Glucuronic acid (GlcA)",
      "D-Galacturonic acid (GalA)", "Glycerol (Gro)",
      "D-Mannitol (Mtl)", "D-Gluconate", "Trehalose (Tre)",
      "Sucrose (Suc)", "Maltose (Mal)", "Cellobiose (Cel)", "Lactose (Lac)"
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.dnmb_cct_used_snfg_types <- function(cyto_nodes = NULL, source_ids = NULL,
                                      monomer_ids = NULL, extra_chains = NULL,
                                      extra_branches = NULL) {
  used <- character(0)
  if (is.data.frame(cyto_nodes) && "sugar_type" %in% names(cyto_nodes)) {
    used <- c(used, as.character(cyto_nodes$sugar_type))
  }
  if (length(source_ids)) {
    source_components <- unlist(lapply(source_ids, function(source_id) {
      .dnmb_cct_carbon_source_render_spec(source_id)$sugar_type
    }), use.names = FALSE)
    used <- c(used, source_components)
  }
  if (length(monomer_ids)) used <- c(used, as.character(monomer_ids))
  if (length(extra_chains)) used <- c(used, unlist(extra_chains, use.names = FALSE))
  if (length(extra_branches)) {
    branch_types <- unlist(lapply(extra_branches, function(branches) {
      if (is.null(branches) || !length(branches)) return(character(0))
      vapply(branches, function(branch) as.character(branch$sugar)[1], character(1))
    }), use.names = FALSE)
    used <- c(used, branch_types)
  }
  used <- unique(.dnmb_cct_canonical_sugar_id(used))
  used <- used[!is.na(used) & nzchar(used)]
  setdiff(used, c("generic", "intermediate", "organic_acid", "phospho_sugar"))
}

.dnmb_cct_snfg_legend_data <- function(sugar_types) {
  used <- unique(.dnmb_cct_canonical_sugar_id(sugar_types))
  used <- used[!is.na(used) & nzchar(used)]
  used <- setdiff(used, c("generic", "intermediate", "organic_acid", "phospho_sugar"))
  catalog <- .dnmb_cct_snfg_legend_catalog()
  known <- catalog[catalog$sugar %in% used, , drop = FALSE]
  unknown <- setdiff(used, catalog$sugar)
  if (length(unknown)) {
    unknown_df <- data.frame(
      sugar = unknown,
      label = tools::toTitleCase(gsub("_", " ", unknown, fixed = TRUE)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    known <- rbind(known, unknown_df)
  }
  rownames(known) <- NULL
  known
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

.dnmb_cct_monomer_lane_map <- function() {
  c(
    glucose = "glucose", galactose = "galactose", mannose = "mannose",
    fructose = "fructose", glcnac = "glcnac", xylose = "xylose",
    arabinose = "arabinose", ribose = "ribose", fucose = "fucose",
    rhamnose = "rhamnose", glca = "glca", gala = "glca",
    glycerol = "glycerol", glucosamine = "glcnac",
    deoxyribose = "deoxyribose"
  )
}

.dnmb_cct_extracellular_monomer_ids <- function() {
  c(
    "glucose", "galactose", "mannose", "fructose", "glcnac",
    "xylose", "arabinose", "ribose", "deoxyribose", "fucose",
    "rhamnose", "glca", "gala"
  )
}

.dnmb_cct_transporter_pathway_map <- function() {
  c(
    glucose = "Glucose", maltose = "Maltose",
    cellobiose = "Cellobiose", galactose = "Galactose",
    trehalose = "Trehalose", lactose = "Lactose",
    mannose = "Mannose", nag = "NAG",
    glucosamine = "Glucosamine", fructose = "Fructose",
    sucrose = "Sucrose", mannitol = "Mannitol",
    glycerol = "Glycerol", xylose = "Xylose",
    arabinose = "Arabinose", ribose = "Ribose",
    fucose = "Fucose", rhamnose = "Rhamnose",
    gluconate = "Gluconate", deoxyribose = "Deoxyribose"
  )
}

.dnmb_cct_transporter_route_id <- function(cs_id, target_key) {
  paste(
    "transport", trimws(as.character(cs_id)[1L]),
    trimws(as.character(target_key)[1L]), sep = ":"
  )
}

# Resolve the visible extracellular source for each selected transporter lane.
# Existing chains and monomer glyphs are preferred; a source with no visible
# extracellular representation receives an explicit substrate glyph at y_mono.
.dnmb_cct_transporter_source_anchors <- function(cs_ids, carbon_src,
                                                 mono_x_map, y_mono,
                                                 extra_center_map = numeric(),
                                                 extra_y_map = numeric()) {
  empty <- data.frame(
    cs_id = character(), source_kind = character(), source_id = character(),
    source_x = numeric(), source_y = numeric(), create_glyph = logical(),
    glyph_source_id = character(), sugar_type = character(),
    cyto_x = numeric(), cyto_y = numeric(), stringsAsFactors = FALSE
  )
  cs_ids <- unique(trimws(as.character(cs_ids)))
  cs_ids <- cs_ids[!is.na(cs_ids) & nzchar(cs_ids)]
  if (!length(cs_ids) || is.null(carbon_src) || !nrow(carbon_src)) return(empty)

  source_sugar <- c(
    Glucose = "glucose", Maltose = "maltose", Cellobiose = "cellobiose",
    Galactose = "galactose", Trehalose = "trehalose",
    Lactose = "lactose", Mannose = "mannose", NAG = "glcnac",
    Glucosamine = "glucosamine", Fructose = "fructose",
    Sucrose = "sucrose", Mannitol = "mannitol", Glycerol = "glycerol",
    Xylose = "xylose", Arabinose = "arabinose", Ribose = "ribose",
    Deoxyribose = "deoxyribose", Fucose = "fucose",
    Rhamnose = "rhamnose", Gluconate = "gluconate"
  )
  chain_names <- names(extra_center_map)
  mono_names <- names(mono_x_map)
  rows <- lapply(cs_ids, function(cs_id) {
    cyto_idx <- match(tolower(cs_id), tolower(as.character(carbon_src$id)))
    if (is.na(cyto_idx)) return(NULL)
    cyto_x <- suppressWarnings(as.numeric(carbon_src$x[cyto_idx]))
    cyto_y <- suppressWarnings(as.numeric(carbon_src$y[cyto_idx]))
    if (!is.finite(cyto_x) || !is.finite(cyto_y)) return(NULL)

    chain_idx <- match(tolower(cs_id), tolower(chain_names))
    if (!is.na(chain_idx)) {
      source_x <- suppressWarnings(as.numeric(extra_center_map[chain_idx]))
      source_y <- suppressWarnings(as.numeric(extra_y_map[chain_idx]))
      if (is.finite(source_x) && is.finite(source_y)) {
        return(data.frame(
          cs_id = cs_id, source_kind = "complex_chain",
          source_id = chain_names[chain_idx], source_x = source_x,
          source_y = source_y, create_glyph = FALSE,
          glyph_source_id = cs_id,
          sugar_type = unname(source_sugar[cs_id]) %||% tolower(cs_id),
          cyto_x = cyto_x, cyto_y = cyto_y,
          stringsAsFactors = FALSE
        ))
      }
    }

    sugar_id <- unname(source_sugar[cs_id])
    if (is.na(sugar_id) || !nzchar(sugar_id)) sugar_id <- tolower(cs_id)
    mono_idx <- match(tolower(sugar_id), tolower(mono_names))
    if (!is.na(mono_idx)) {
      source_x <- suppressWarnings(as.numeric(mono_x_map[mono_idx]))
      if (is.finite(source_x) && is.finite(y_mono)) {
        return(data.frame(
          cs_id = cs_id, source_kind = "monosaccharide",
          source_id = mono_names[mono_idx], source_x = source_x,
          source_y = as.numeric(y_mono), create_glyph = FALSE,
          glyph_source_id = cs_id, sugar_type = sugar_id,
          cyto_x = cyto_x, cyto_y = cyto_y,
          stringsAsFactors = FALSE
        ))
      }
    }

    data.frame(
      cs_id = cs_id, source_kind = "explicit_source",
      source_id = paste0("uptake:", tolower(cs_id)), source_x = cyto_x,
      source_y = as.numeric(y_mono), create_glyph = TRUE,
      glyph_source_id = cs_id, sugar_type = sugar_id,
      cyto_x = cyto_x, cyto_y = cyto_y,
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(empty)
  explicit_idx <- which(out$create_glyph)
  if (length(explicit_idx)) {
    occupied_x <- suppressWarnings(as.numeric(mono_x_map))
    occupied_x <- occupied_x[is.finite(occupied_x)]
    explicit_idx <- explicit_idx[order(out$source_x[explicit_idx], out$cs_id[explicit_idx])]
    candidate_offsets <- c(0, as.vector(rbind(seq(0.5, 4, by = 0.5), -seq(0.5, 4, by = 0.5))))
    for (idx in explicit_idx) {
      desired <- .dnmb_cct_snap_to_grid(out$source_x[idx], step = 0.25)
      candidates <- desired + candidate_offsets
      free <- vapply(candidates, function(candidate) {
        !length(occupied_x) || all(abs(candidate - occupied_x) >= 0.50 - 1e-10)
      }, logical(1))
      chosen <- if (any(free)) candidates[which(free)[1L]] else desired
      out$source_x[idx] <- chosen
      occupied_x <- c(occupied_x, chosen)
    }
  }
  out[order(out$source_x, out$cs_id), , drop = FALSE]
}

.dnmb_cct_separate_lanes <- function(x, min_gap = 0.5, snap_step = 0.25) {
  original_names <- names(x)
  x <- as.numeric(x)
  names(x) <- if (!is.null(original_names)) original_names else as.character(seq_along(x))
  finite <- is.finite(x)
  if (sum(finite) < 2L) return(x)

  idx <- which(finite)
  ord <- idx[order(x[idx], idx)]
  target <- .dnmb_cct_snap_to_grid(x[ord], step = snap_step)
  placed <- target
  for (i in seq_along(placed)[-1L]) {
    placed[i] <- max(target[i], placed[i - 1L] + min_gap)
  }
  placed <- placed + .dnmb_cct_snap_to_grid(
    mean(target) - mean(placed), step = snap_step
  )
  x[ord] <- placed
  x
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
  generic_re <- "transport|permease|porter|uptake|symporter|antiporter|binding.protein|sodium.?solute|\\bsbp\\b|\\babc\\b|\\bmfs\\b|\\bsss\\b|\\bpts\\b|\\beiic|\\beiib|\\beiia"
  # Specific gene names (use word boundaries, check each field separately)
  gene_re <- "\\b(mgl[ABC]|msiK|malE|malF|malG|malK|thu[EFGK]|rbs[ABC]|xyl[FGHT]|ggu[AB]|glc[TUVP]|galP|lac[PEFG]|chvE|gts[ABCD]|kguT|ara[ESTUV]|mtl[AEK]|fucP|deoP|nagE|nagF|nagPcb|fruA|fruB|fruD|scrT|sacP|cscB|treB|manP|bglF|celB|ascB|ptsG|ptsH|ptsI|crr)\\b"
  sid <- as.character(step_id)[1]
  gnm <- as.character(gene_name)[1]
  if (grepl(generic_re, sid, ignore.case = TRUE)) return(TRUE)
  if (!is.na(gnm) && grepl(generic_re, gnm, ignore.case = TRUE)) return(TRUE)
  if (grepl(gene_re, sid, ignore.case = TRUE)) return(TRUE)
  if (!is.na(gnm) && grepl(gene_re, gnm, ignore.case = TRUE)) return(TRUE)
  FALSE
}

# Build one physical transporter entity per locus while retaining every
# supported locus/pathway/model relationship in a separate membership table.
.dnmb_cct_transporter_entities <- function(transporters) {
  empty_entities <- data.frame(
    entity_id = character(), locus_tag = character(),
    model_names = character(), substrate_names = character(),
    pathway_names = character(), n_models = integer(),
    n_substrates = integer(), n_memberships = integer(),
    shared = logical(), label = character(),
    stringsAsFactors = FALSE
  )
  empty_memberships <- data.frame(
    entity_id = character(), locus_tag = character(),
    pathway = character(), step = character(), cs_id = character(),
    confidence = character(), step_score = numeric(),
    shared = logical(), stringsAsFactors = FALSE
  )
  if (is.null(transporters) || !is.data.frame(transporters) || !nrow(transporters)) {
    return(list(entities = empty_entities, memberships = empty_memberships))
  }

  required <- c("locus_tag", "pathway", "step")
  missing_required <- setdiff(required, names(transporters))
  if (length(missing_required)) {
    stop(
      "transporter entity construction requires columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }

  rows <- as.data.frame(transporters, stringsAsFactors = FALSE, row.names = NULL)
  rows$.source_order <- seq_len(nrow(rows))
  rows$locus_tag <- trimws(.dnmb_cct_normalize_locus_tag(rows$locus_tag))
  rows$pathway <- trimws(as.character(rows$pathway))
  rows$step <- trimws(as.character(rows$step))
  if (!"cs_id" %in% names(rows)) rows$cs_id <- rows$pathway
  rows$cs_id <- trimws(as.character(rows$cs_id))
  missing_cs <- is.na(rows$cs_id) | !nzchar(rows$cs_id)
  rows$cs_id[missing_cs] <- rows$pathway[missing_cs]

  if (!"confidence" %in% names(rows)) rows$confidence <- NA_character_
  rows$confidence <- tolower(trimws(as.character(rows$confidence)))
  if (!"step_score" %in% names(rows)) rows$step_score <- NA_real_
  rows$step_score <- suppressWarnings(as.numeric(rows$step_score))
  rows$.confidence_rank <- unname(c(none = 0L, low = 1L, medium = 2L, high = 3L)[rows$confidence])
  rows$.confidence_rank[is.na(rows$.confidence_rank)] <- 0L

  valid_locus <- !is.na(rows$locus_tag) & nzchar(rows$locus_tag)
  valid_membership <- !is.na(rows$pathway) & nzchar(rows$pathway) &
    !is.na(rows$step) & nzchar(rows$step) &
    !is.na(rows$cs_id) & nzchar(rows$cs_id)
  foreground <- rows$confidence %in% c("high", "medium") |
    (!is.na(rows$step_score) & rows$step_score >= 1)
  rows <- rows[valid_locus & valid_membership & foreground, , drop = FALSE]
  if (!nrow(rows)) {
    return(list(entities = empty_entities, memberships = empty_memberships))
  }

  score_order <- ifelse(is.na(rows$step_score), -Inf, rows$step_score)
  rows <- rows[order(
    rows$locus_tag, rows$pathway, rows$step,
    -rows$.confidence_rank, -score_order, rows$.source_order
  ), , drop = FALSE]
  membership_key <- paste(rows$locus_tag, rows$pathway, rows$step, sep = "\r")
  rows <- rows[!duplicated(membership_key), , drop = FALSE]

  entity_rows <- lapply(split(seq_len(nrow(rows)), rows$locus_tag), function(idx) {
    member <- rows[idx, , drop = FALSE]
    models <- unique(member$step)
    substrates <- unique(member$cs_id)
    pathways <- unique(member$pathway)
    is_shared <- length(substrates) > 1L
    locus <- member$locus_tag[1]
    model_text <- paste(models, collapse = " | ")
    substrate_text <- paste(substrates, collapse = " | ")
    data.frame(
      entity_id = locus,
      locus_tag = locus,
      model_names = model_text,
      substrate_names = substrate_text,
      pathway_names = paste(pathways, collapse = " | "),
      n_models = length(models),
      n_substrates = length(substrates),
      n_memberships = nrow(member),
      shared = is_shared,
      label = paste0(
        locus,
        "\nModels: ", model_text,
        "\nSubstrates: ", substrate_text,
        "\nShared: ", if (is_shared) "yes" else "no"
      ),
      stringsAsFactors = FALSE
    )
  })
  entities <- do.call(rbind, entity_rows)
  rownames(entities) <- NULL
  entities <- entities[order(entities$locus_tag), , drop = FALSE]

  rows$entity_id <- rows$locus_tag
  rows$shared <- entities$shared[match(rows$locus_tag, entities$locus_tag)]
  private_columns <- c(".source_order", ".confidence_rank")
  key_columns <- c(
    "entity_id", "locus_tag", "pathway", "step", "cs_id",
    "confidence", "step_score", "shared"
  )
  memberships <- rows[, c(key_columns, setdiff(names(rows), c(key_columns, private_columns))), drop = FALSE]
  rownames(memberships) <- NULL

  list(entities = entities, memberships = memberships)
}

.dnmb_cct_transporter_half_span <- function(step_id, gene_name = NA_character_, is_pts = FALSE) {
  txt <- paste(step_id, gene_name)
  if (isTRUE(is_pts)) return(0.34)
  if (grepl("abc|sbp|binding|malE|mgl[ABC]|rbs[ABC]|xyl[FGH]|thu[EFG]", txt, ignore.case = TRUE)) return(0.30)
  if (grepl("mfs|sss|sodium.?solute|permease|symporter|transporter|porter|xylT|gnt|lac|gal|man|srl|iolT|kgtP", txt, ignore.case = TRUE)) return(0.26)
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
  if (grepl("mfs|sss|sodium.?solute|permease|symporter|transporter|porter|xylT|gnt|lac|gal|man|srl|iolT|kgtP|thuK|glpO", txt, ignore.case = TRUE)) {
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
                                             desired_x = NULL,
                                             stable_ids = NULL,
                                             y_memb = 8.5,
                                             xlim = NULL,
                                             min_gap = 0.08,
                                             coincident_tolerance = 1e-6,
                                             coincident_fan_gap = 0.22,
                                             row_step = 0.12,
                                             max_rows = 3L) {
  n <- length(half_spans)
  if (n == 0) return(data.frame(tx = numeric(0), ty = numeric(0), row_id = integer(0)))
  half_spans <- pmax(0.04, as.numeric(half_spans))
  desired_x <- if (is.null(desired_x) || length(desired_x) != n) {
    rep(as.numeric(center_x)[1], n)
  } else {
    as.numeric(desired_x)
  }
  priorities <- if (missing(lane_ranks) || length(lane_ranks) != n) {
    rep(1, n)
  } else {
    as.numeric(lane_ranks)
  }
  priorities[!is.finite(priorities)] <- 0
  stable_ids <- if (is.null(stable_ids) || length(stable_ids) != n) {
    sprintf("%08d", seq_len(n))
  } else {
    as.character(stable_ids)
  }
  missing_ids <- is.na(stable_ids) | !nzchar(stable_ids)
  stable_ids[missing_ids] <- sprintf("%08d", which(missing_ids))
  min_gap <- max(0, as.numeric(min_gap)[1])
  coincident_tolerance <- max(0, as.numeric(coincident_tolerance)[1])
  coincident_fan_gap <- max(0, as.numeric(coincident_fan_gap)[1])
  row_step <- max(0.10, as.numeric(row_step)[1])
  max_rows <- max(1L, min(3L, as.integer(max_rows)[1]))

  # Distinct loci can share the same substrate anchor. Fan those centers out
  # before row packing so their glyphs and route trunks cannot read as one
  # serial transporter. Stable locus IDs make the layout input-order invariant.
  original_x <- desired_x
  finite_ord <- order(original_x, stable_ids, na.last = NA)
  if (length(finite_ord) > 1L && coincident_fan_gap > 0) {
    group_start <- 1L
    for (jj in seq.int(2L, length(finite_ord) + 1L)) {
      group_done <- jj > length(finite_ord) ||
        abs(original_x[finite_ord[jj]] - original_x[finite_ord[jj - 1L]]) >
          coincident_tolerance
      if (!group_done) next
      idx <- finite_ord[group_start:(jj - 1L)]
      if (length(idx) > 1L) {
        idx <- idx[order(-priorities[idx], stable_ids[idx])]
        fan_gap <- max(coincident_fan_gap, max(half_spans[idx]) + min_gap)
        fan_step <- integer(length(idx))
        if (length(idx) > 1L) {
          kk <- seq_len(length(idx) - 1L)
          fan_step[-1L] <- ceiling(kk / 2) * ifelse(kk %% 2L == 1L, 1L, -1L)
        }
        candidate_x <- original_x[idx[1L]] + fan_step * fan_gap

        if (!is.null(xlim) && length(xlim) >= 2L && all(is.finite(xlim[1:2]))) {
          bounds <- sort(as.numeric(xlim[1:2]))
          left <- min(candidate_x - half_spans[idx])
          right <- max(candidate_x + half_spans[idx])
          if (right > bounds[2]) candidate_x <- candidate_x - (right - bounds[2])
          left <- min(candidate_x - half_spans[idx])
          if (left < bounds[1]) candidate_x <- candidate_x + (bounds[1] - left)
        }
        desired_x[idx] <- candidate_x
      }
      group_start <- jj
    }
  }

  # Greedy interval coloring keeps non-colliding glyphs on the membrane
  # centerline. Only genuinely overlapping intervals use the two auxiliary
  # membrane rows.
  row_right <- rep(-Inf, max_rows)
  row_id <- integer(n)
  ord <- order(desired_x, -priorities, stable_ids, na.last = TRUE)
  for (ii in ord) {
    free <- which(desired_x[ii] - half_spans[ii] >= row_right + min_gap)
    row <- if (length(free)) free[1] else which.min(row_right)
    row_id[ii] <- row
    placed_x <- max(desired_x[ii], row_right[row] + min_gap + half_spans[ii])
    row_right[row] <- placed_x + half_spans[ii]
  }

  tx <- desired_x
  for (row in seq_len(max_rows)) {
    idx <- which(row_id == row)
    if (!length(idx)) next
    idx <- idx[order(desired_x[idx], -priorities[idx], stable_ids[idx])]
    if (length(idx) > 1L) {
      offsets <- numeric(length(idx))
      for (jj in 2:length(idx)) {
        offsets[jj] <- offsets[jj - 1L] + half_spans[idx[jj - 1L]] +
          half_spans[idx[jj]] + min_gap
      }
      target <- desired_x[idx] - offsets
      tx[idx] <- stats::isoreg(seq_along(target), target)$yf + offsets
    }

    if (!is.null(xlim) && length(xlim) >= 2L && all(is.finite(xlim[1:2]))) {
      bounds <- sort(as.numeric(xlim[1:2]))
      left <- min(tx[idx] - half_spans[idx])
      right <- max(tx[idx] + half_spans[idx])
      if (right > bounds[2]) tx[idx] <- tx[idx] - (right - bounds[2])
      left <- min(tx[idx] - half_spans[idx])
      if (left < bounds[1]) tx[idx] <- tx[idx] + (bounds[1] - left)
    }
  }

  row_offsets <- c(0, row_step, -row_step)[seq_len(max_rows)]
  data.frame(
    tx = tx,
    ty = y_memb + row_offsets[row_id],
    row_id = row_id
  )
}

# Route transporter memberships through one compact bus per substrate lane.
# At shared transporter loci, each substrate keeps its own coloured landing
# track so no membership disappears into a neutral common trunk.
.dnmb_cct_transporter_bus_layout <- function(memberships, y_memb = 8.5,
                                             bus_offset = 0.70,
                                             tier_step = 0.065,
                                             max_bus_offset = 1.12,
                                             source_offset = 0.76,
                                             glyph_clearance = 0.11,
                                             membrane_half_height = 0.20,
                                             route_gap = 0.10,
                                             interval_gap = 0.05,
                                             corner_radius = 0.08,
                                             suppress_redundant_medium = TRUE) {
  empty_segment <- data.frame(
    group = character(), x = numeric(), xend = numeric(),
    y = numeric(), yend = numeric(), color = character(),
    confidence = character(), linewidth = numeric(), alpha = numeric(),
    linetype = character(), stringsAsFactors = FALSE
  )
  empty_group <- data.frame(
    group = character(), cs_id = character(), lane_x = numeric(),
    source_x = numeric(), source_y = numeric(),
    x_left = numeric(), x_right = numeric(), tier = integer(),
    bus_y = numeric(), confidence = character(), color = character(),
    direct = logical(), stringsAsFactors = FALSE
  )
  empty_membership <- data.frame(
    group = character(), target_key = character(), cs_id = character(),
    lane_x = numeric(), tx_draw = numeric(), ty_draw = numeric(),
    bus_y = numeric(), direct = logical(), draw = logical(),
    stringsAsFactors = FALSE
  )
  empty <- list(
    groups = empty_group, memberships = empty_membership,
    stems = empty_segment, buses = empty_segment,
    trunks = empty_segment, direct = empty_segment,
    exterior_routes = list(), exterior_source_routes = list(),
    exterior_trunk_routes = list(), exterior_render_routes = list()
  )
  if (is.null(memberships) || !is.data.frame(memberships) || !nrow(memberships)) {
    return(empty)
  }

  required <- c("cs_id", "lane_x", "tx_draw", "ty_draw")
  missing_required <- setdiff(required, names(memberships))
  if (length(missing_required)) {
    stop(
      "transporter bus layout requires columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  rows <- as.data.frame(memberships, stringsAsFactors = FALSE, row.names = NULL)
  if (!"confidence" %in% names(rows)) rows$confidence <- "medium"
  if (!"color" %in% names(rows)) rows$color <- "#78909C"
  rows$cs_id <- trimws(as.character(rows$cs_id))
  rows$confidence <- tolower(trimws(as.character(rows$confidence)))
  rows$confidence[!rows$confidence %in% c("high", "medium")] <- "medium"
  rows$color <- as.character(rows$color)
  rows$color[is.na(rows$color) | !nzchar(rows$color)] <- "#78909C"
  rows$lane_x <- suppressWarnings(as.numeric(rows$lane_x))
  if (!"source_x" %in% names(rows)) rows$source_x <- rows$lane_x
  if (!"source_y" %in% names(rows)) rows$source_y <- y_memb + source_offset
  if (!"source_trim" %in% names(rows)) rows$source_trim <- 0.10
  rows$source_x <- suppressWarnings(as.numeric(rows$source_x))
  rows$source_y <- suppressWarnings(as.numeric(rows$source_y))
  rows$source_trim <- suppressWarnings(as.numeric(rows$source_trim))
  rows$source_trim[!is.finite(rows$source_trim) | rows$source_trim < 0] <- 0.10
  missing_source_x <- !is.finite(rows$source_x)
  rows$source_x[missing_source_x] <- rows$lane_x[missing_source_x]
  rows$lane_x <- rows$source_x
  rows$tx_draw <- suppressWarnings(as.numeric(rows$tx_draw))
  rows$ty_draw <- suppressWarnings(as.numeric(rows$ty_draw))
  rows <- rows[
    !is.na(rows$cs_id) & nzchar(rows$cs_id) &
      is.finite(rows$lane_x) & is.finite(rows$source_y) &
      is.finite(rows$tx_draw) & is.finite(rows$ty_draw),
    , drop = FALSE
  ]
  if (!nrow(rows)) return(empty)
  membrane_half_height <- max(0, as.numeric(membrane_half_height)[1L])
  route_gap <- max(0, as.numeric(route_gap)[1L])
  min_bus_offset <- max(
    membrane_half_height,
    max(rows$ty_draw + glyph_clearance - y_memb, na.rm = TRUE)
  ) + route_gap
  bus_offset <- max(as.numeric(bus_offset)[1L], min_bus_offset)
  max_bus_offset <- max(as.numeric(max_bus_offset)[1L], bus_offset)

  target_id <- if ("locus_tag" %in% names(rows)) {
    as.character(rows$locus_tag)
  } else if ("anchor_id" %in% names(rows)) {
    as.character(rows$anchor_id)
  } else {
    rep(NA_character_, nrow(rows))
  }
  missing_target <- is.na(target_id) | !nzchar(target_id)
  target_id[missing_target] <- paste0(
    "xy:", formatC(rows$tx_draw[missing_target], digits = 6, format = "f"),
    ":", formatC(rows$ty_draw[missing_target], digits = 6, format = "f")
  )
  rows$target_key <- target_id
  rows <- rows[order(
    rows$lane_x, rows$cs_id, rows$target_key,
    -as.integer(rows$confidence == "high"), rows$tx_draw
  ), , drop = FALSE]
  rows <- rows[!duplicated(rows[, c("cs_id", "lane_x", "target_key")]), , drop = FALSE]

  # A medium connector adds no visual information when the same substrate lane
  # already has high-confidence uptake. Its locus remains in the evidence ledger.
  rows$draw <- TRUE
  if (isTRUE(suppress_redundant_medium)) {
    lane_has_high <- ave(
      rows$confidence == "high",
      interaction(rows$cs_id, rows$lane_x, drop = TRUE),
      FUN = any
    )
    rows$draw <- rows$confidence == "high" | !lane_has_high
  }
  draw_rows <- rows[rows$draw, , drop = FALSE]
  if (!nrow(draw_rows)) return(empty)

  lane_key <- paste0(
    draw_rows$cs_id, "@",
    formatC(draw_rows$lane_x, digits = 6, format = "f")
  )
  draw_rows$group <- lane_key
  lane_split <- split(seq_len(nrow(draw_rows)), draw_rows$group)
  group_rows <- lapply(names(lane_split), function(group_id) {
    idx <- lane_split[[group_id]]
    lane_x <- stats::median(draw_rows$source_x[idx])
    source_y <- stats::median(draw_rows$source_y[idx])
    target_x <- unique(draw_rows$tx_draw[idx])
    direct <- length(target_x) == 1L && abs(target_x - lane_x) <= 1e-8
    data.frame(
      group = group_id,
      cs_id = draw_rows$cs_id[idx[1L]],
      lane_x = lane_x,
      source_x = lane_x,
      source_y = source_y,
      x_left = min(c(lane_x, target_x)),
      x_right = max(c(lane_x, target_x)),
      tier = 0L,
      bus_y = NA_real_,
      confidence = if (any(draw_rows$confidence[idx] == "high")) "high" else "medium",
      color = draw_rows$color[idx[1L]],
      direct = direct,
      stringsAsFactors = FALSE
    )
  })
  groups <- dplyr::bind_rows(group_rows)
  groups <- groups[order(groups$x_left, groups$x_right, groups$group), , drop = FALSE]

  # Greedy interval coloring guarantees that horizontal buses on a tier do not
  # overlap. The tier spacing contracts only when needed to remain below the
  # extracellular monomer row.
  bus_idx <- which(!groups$direct)
  if (length(bus_idx)) {
    tier_right <- numeric(0)
    for (ii in bus_idx) {
      free <- which(groups$x_left[ii] >= tier_right + interval_gap)
      tier <- if (length(free)) free[1L] else length(tier_right) + 1L
      if (tier > length(tier_right)) tier_right <- c(tier_right, -Inf)
      groups$tier[ii] <- tier
      tier_right[tier] <- groups$x_right[ii]
    }
    n_tiers <- max(groups$tier[bus_idx])
    step_used <- if (n_tiers <= 1L) 0 else min(
      tier_step,
      max(0, (max_bus_offset - bus_offset) / (n_tiers - 1L))
    )
    groups$bus_y[bus_idx] <- y_memb + bus_offset +
      (groups$tier[bus_idx] - 1L) * step_used
  }
  draw_rows$bus_y <- groups$bus_y[match(draw_rows$group, groups$group)]
  draw_rows$direct <- groups$direct[match(draw_rows$group, groups$group)]

  style_for <- function(confidence, bus = FALSE) {
    high <- confidence == "high"
    data.frame(
      linewidth = ifelse(high, if (bus) 0.28 else 0.24, if (bus) 0.18 else 0.16),
      alpha = ifelse(high, 0.68, 0.28),
      linetype = ifelse(high, "solid", "dashed"),
      stringsAsFactors = FALSE
    )
  }
  segment_frame <- function(group, x, xend, y, yend, color, confidence,
                            bus = FALSE) {
    style <- style_for(confidence, bus = bus)
    data.frame(
      group = group, x = x, xend = xend, y = y, yend = yend,
      color = color, confidence = confidence,
      linewidth = style$linewidth, alpha = style$alpha,
      linetype = style$linetype, stringsAsFactors = FALSE
    )
  }

  buses <- if (length(bus_idx)) {
    bus_groups <- groups[bus_idx, , drop = FALSE]
    segment_frame(
      bus_groups$group, bus_groups$x_left, bus_groups$x_right,
      bus_groups$bus_y, bus_groups$bus_y,
      bus_groups$color, bus_groups$confidence, bus = TRUE
    )
  } else empty_segment

  # If a transporter trunk is already on the source x, it also serves as the
  # substrate stem. The trunk is extended to source_y below so this omission
  # never leaves a gap when the same group also has off-axis targets.
  stem_groups <- groups[bus_idx, , drop = FALSE]
  has_lane_trunk <- vapply(stem_groups$group, function(group_id) {
    idx <- which(draw_rows$group == group_id)
    any(abs(draw_rows$tx_draw[idx] - draw_rows$lane_x[idx]) <= 0.015)
  }, logical(1))
  stem_groups <- stem_groups[!has_lane_trunk, , drop = FALSE]
  stems <- if (nrow(stem_groups)) {
    segment_frame(
      stem_groups$group, stem_groups$lane_x, stem_groups$lane_x,
      stem_groups$source_y, stem_groups$bus_y,
      stem_groups$color, stem_groups$confidence
    )
  } else empty_segment

  # One trunk per physical transporter. Every substrate bus intersects this
  # trunk at its own tier, preserving shared-locus traceability without fans.
  trunk_rows <- draw_rows[!draw_rows$direct, , drop = FALSE]
  trunk_split <- split(seq_len(nrow(trunk_rows)), trunk_rows$target_key)
  trunks <- if (length(trunk_split)) {
    dplyr::bind_rows(lapply(names(trunk_split), function(target_key) {
      idx <- trunk_split[[target_key]]
      target_conf <- if (any(trunk_rows$confidence[idx] == "high")) "high" else "medium"
      target_colors <- unique(trunk_rows$color[idx])
      trunk_color <- if (length(unique(trunk_rows$group[idx])) > 1L) "#52616B" else target_colors[1L]
      aligned_source <- any(
        draw_rows$target_key == target_key &
          abs(draw_rows$lane_x - draw_rows$tx_draw) <= 1e-8
      )
      trunk_top <- max(
        trunk_rows$bus_y[idx],
        if (aligned_source) trunk_rows$source_y[idx] else -Inf
      )
      segment_frame(
        paste0("target:", target_key),
        trunk_rows$tx_draw[idx[1L]], trunk_rows$tx_draw[idx[1L]],
        trunk_top,
        trunk_rows$ty_draw[idx[1L]] + glyph_clearance,
        trunk_color, target_conf
      )
    }))
  } else empty_segment

  direct_rows <- draw_rows[draw_rows$direct, , drop = FALSE]
  # A shared trunk from another bus already passes through an aligned source.
  trunk_targets <- names(trunk_split)
  direct_rows <- direct_rows[!direct_rows$target_key %in% trunk_targets, , drop = FALSE]
  direct <- if (nrow(direct_rows)) {
    segment_frame(
      direct_rows$group, direct_rows$lane_x, direct_rows$tx_draw,
      direct_rows$source_y,
      direct_rows$ty_draw + glyph_clearance,
      direct_rows$color, direct_rows$confidence
    )
  } else empty_segment

  rows$group <- paste0(
    rows$cs_id, "@", formatC(rows$lane_x, digits = 6, format = "f")
  )
  rows$bus_y <- groups$bus_y[match(rows$group, groups$group)]
  rows$direct <- groups$direct[match(rows$group, groups$group)]
  membership_out <- rows[, c(
    "group", "target_key", "cs_id", "lane_x", "source_x", "source_y",
    "source_trim", "tx_draw", "ty_draw", "bus_y", "direct", "draw"
  ), drop = FALSE]

  # The plotted exterior connectors are complete raw orthogonal paths. Apply
  # the corner radius only once after all bus coordinates are fixed, avoiding
  # curve-then-step artifacts at the horizontal lane.
  exterior_routes <- lapply(seq_len(nrow(draw_rows)), function(i) {
    row <- draw_rows[i, , drop = FALSE]
    target_y <- row$ty_draw + glyph_clearance
    raw <- if (isTRUE(row$direct)) {
      data.frame(x = c(row$source_x, row$tx_draw),
                 y = c(row$source_y, target_y))
    } else {
      data.frame(
        x = c(row$source_x, row$source_x, row$tx_draw, row$tx_draw),
        y = c(row$source_y, row$bus_y, row$bus_y, target_y)
      )
    }
    raw <- raw[c(TRUE, diff(raw$x) != 0 | diff(raw$y) != 0), , drop = FALSE]
    points <- .dnmb_cct_round_orthogonal_route(
      raw, radius = max(0.08, as.numeric(corner_radius)[1L])
    )
    style <- style_for(row$confidence)
    list(
      route_id = .dnmb_cct_transporter_route_id(row$cs_id, row$target_key),
      cs_id = row$cs_id, target_key = row$target_key,
      source_x = row$source_x, source_y = row$source_y,
      target_x = row$tx_draw, target_y = target_y,
      points = points, color = row$color,
      confidence = row$confidence, linewidth = style$linewidth,
      alpha = style$alpha, linetype = style$linetype,
      trim_start = row$source_trim, arrow_last = TRUE
    )
  })

  # Parallel landing tracks remain inside the transporter glyph envelope. One
  # route stays on the physical centre; additional routes alternate right and
  # left in deterministic substrate order so every colour remains visible.
  render_target_x <- draw_rows$tx_draw
  target_split <- split(seq_len(nrow(draw_rows)), draw_rows$target_key)
  for (target_key in names(target_split)) {
    idx <- target_split[[target_key]]
    if (length(idx) < 2L) next
    aligned <- abs(draw_rows$source_x[idx] - draw_rows$tx_draw[idx]) <= 1e-8
    ord <- idx[order(
      -as.integer(aligned),
      -as.integer(draw_rows$confidence[idx] == "high"),
      draw_rows$cs_id[idx], draw_rows$group[idx], draw_rows$source_x[idx]
    )]
    kk <- seq_len(length(ord) - 1L)
    track_step <- integer(length(ord))
    track_step[-1L] <- ceiling(kk / 2) * ifelse(kk %% 2L == 1L, 1L, -1L)
    track_offset <- track_step * 0.035
    if (max(abs(track_offset)) > 0.07) {
      track_offset <- track_offset * 0.07 / max(abs(track_offset))
    }
    render_target_x[ord] <- stats::median(draw_rows$tx_draw[idx]) + track_offset
  }

  exterior_render_routes <- lapply(seq_len(nrow(draw_rows)), function(i) {
    row <- draw_rows[i, , drop = FALSE]
    landing_x <- render_target_x[i]
    target_y <- row$ty_draw + glyph_clearance
    has_fan_in <- abs(landing_x - row$tx_draw) > 1e-8
    landing_y <- target_y + if (has_fan_in) 0.075 else 0
    raw <- if (isTRUE(row$direct) && abs(row$source_x - landing_x) <= 1e-8) {
      data.frame(x = c(row$source_x, landing_x), y = c(row$source_y, landing_y))
    } else {
      join_y <- if (isTRUE(row$direct)) y_memb + bus_offset else row$bus_y
      data.frame(
        x = c(row$source_x, row$source_x, landing_x, landing_x),
        y = c(row$source_y, join_y, join_y, landing_y)
      )
    }
    raw <- raw[c(TRUE, diff(raw$x) != 0 | diff(raw$y) != 0), , drop = FALSE]
    points <- .dnmb_cct_round_orthogonal_route(
      raw, radius = max(0.08, as.numeric(corner_radius)[1L])
    )
    if (has_fan_in) {
      points <- rbind(
        points,
        data.frame(x = row$tx_draw, y = target_y)
      )
    }
    style <- style_for(row$confidence)
    list(
      route_id = .dnmb_cct_transporter_route_id(row$cs_id, row$target_key),
      cs_id = row$cs_id, target_key = row$target_key,
      source_x = row$source_x, source_y = row$source_y,
      transporter_x = row$tx_draw, landing_x = landing_x,
      landing_y = landing_y,
      target_x = row$tx_draw, target_y = target_y,
      points = points,
      color = row$color, confidence = row$confidence,
      linewidth = style$linewidth, alpha = style$alpha,
      linetype = style$linetype, trim_start = row$source_trim,
      arrow_last = TRUE
    )
  })
  exterior_source_routes <- exterior_render_routes
  exterior_trunk_routes <- list()
  list(
    groups = groups, memberships = membership_out,
    stems = stems, buses = buses, trunks = trunks, direct = direct,
    exterior_routes = exterior_routes,
    exterior_source_routes = exterior_source_routes,
    exterior_trunk_routes = exterior_trunk_routes,
    exterior_render_routes = exterior_render_routes
  )
}

.dnmb_cct_transporter_interior_routes <- function(memberships, y_memb = 8.5,
                                                  glyph_clearance = 0.11,
                                                  source_clearance = 0.11,
                                                  membrane_half_height = 0.20,
                                                  route_gap = 0.10,
                                                  lane_step = 0.045,
                                                  max_lane_offset = 0.18,
                                                  corner_radius = 0.08) {
  if (is.null(memberships) || !is.data.frame(memberships) || !nrow(memberships)) {
    return(list())
  }
  required <- c(
    "cs_id", "tx_draw", "ty_draw", "cyto_x", "cyto_y",
    "confidence", "color"
  )
  if (length(missing <- setdiff(required, names(memberships)))) {
    stop(
      "transporter interior routing requires columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  rows <- as.data.frame(memberships, stringsAsFactors = FALSE, row.names = NULL)
  rows <- rows[
    is.finite(rows$tx_draw) & is.finite(rows$ty_draw) &
      is.finite(rows$cyto_x) & is.finite(rows$cyto_y),
    , drop = FALSE
  ]
  if (!nrow(rows)) return(list())
  if (!"locus_tag" %in% names(rows)) rows$locus_tag <- seq_len(nrow(rows))
  rows <- rows[order(rows$cyto_x, rows$cs_id, rows$tx_draw, rows$locus_tag), , drop = FALSE]
  rows <- rows[!duplicated(rows[, c("cs_id", "locus_tag", "tx_draw", "cyto_x")]), , drop = FALSE]
  membrane_half_height <- max(0, as.numeric(membrane_half_height)[1L])
  route_gap <- max(0, as.numeric(route_gap)[1L])
  max_lane_offset <- max(0, as.numeric(max_lane_offset)[1L])

  route_specs <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, , drop = FALSE]
    start_y <- row$ty_draw - glyph_clearance
    end_y <- row$cyto_y + source_clearance
    lower_envelope <- min(y_memb - membrane_half_height, start_y)
    lane_y <- lower_envelope - route_gap - max_lane_offset
    lane_y <- max(end_y + 0.06, lane_y)
    raw <- if (abs(row$tx_draw - row$cyto_x) <= 1e-8) {
      data.frame(x = c(row$tx_draw, row$cyto_x), y = c(start_y, end_y))
    } else {
      data.frame(
        x = c(row$tx_draw, row$tx_draw, row$cyto_x, row$cyto_x),
        y = c(start_y, lane_y, lane_y, end_y)
      )
    }
    list(
      route_id = .dnmb_cct_transporter_route_id(row$cs_id, row$locus_tag),
      group = paste0("inside:", row$cs_id),
      priority = if (tolower(row$confidence) == "high") 2 else 1,
      raw = raw, row = row,
      min_lane_y = end_y + 0.06,
      max_lane_y = lower_envelope - route_gap
    )
  })
  metadata <- data.frame(
    route_id = vapply(route_specs, `[[`, character(1), "route_id"),
    group = vapply(route_specs, `[[`, character(1), "group"),
    priority = vapply(route_specs, `[[`, numeric(1), "priority"),
    stringsAsFactors = FALSE
  )
  unrounded_paths <- .dnmb_cct_separate_route_lanes(
    routes = lapply(route_specs, `[[`, "raw"), metadata = metadata,
    lane_step = lane_step, max_offset = max_lane_offset,
    horizontal_tolerance = 0.01, y_tolerance = 0.05,
    min_overlap = 0.10, connection_mode = "move_internal",
    round_radius = 0
  )
  paths <- lapply(seq_along(route_specs), function(i) {
    spec <- route_specs[[i]]
    pts <- unrounded_paths[[i]]
    if (nrow(spec$raw) <= 2L) return(pts)
    dx <- diff(pts$x)
    dy <- diff(pts$y)
    horizontal_y <- pts$y[-nrow(pts)][abs(dx) > 1e-8 & abs(dy) <= 1e-8]
    lane_y <- if (length(horizontal_y)) stats::median(horizontal_y) else {
      stats::median(spec$raw$y[2:(nrow(spec$raw) - 1L)])
    }
    lane_y <- min(spec$max_lane_y, max(spec$min_lane_y, lane_y))
    raw <- data.frame(
      x = c(spec$raw$x[1L], spec$raw$x[1L],
            spec$raw$x[nrow(spec$raw)], spec$raw$x[nrow(spec$raw)]),
      y = c(spec$raw$y[1L], lane_y, lane_y, spec$raw$y[nrow(spec$raw)])
    )
    raw <- raw[c(TRUE, diff(raw$x) != 0 | diff(raw$y) != 0), , drop = FALSE]
    .dnmb_cct_round_orthogonal_route(
      raw, radius = max(0.08, as.numeric(corner_radius)[1L])
    )
  })
  lapply(seq_along(route_specs), function(i) {
    row <- route_specs[[i]]$row
    high <- tolower(row$confidence) == "high"
    c(route_specs[[i]][c("route_id", "group", "priority")], list(
      cs_id = row$cs_id, target_key = as.character(row$locus_tag),
      points = paths[[i]], color = row$color,
      confidence = row$confidence,
      linewidth = if (high) 0.24 else 0.16,
      alpha = if (high) 0.68 else 0.28,
      linetype = if (high) "solid" else "dashed"
    ))
  })
}

.dnmb_cct_drawn_transporter_entities <- function(entities, memberships) {
  if (is.null(entities) || !is.data.frame(entities) || !nrow(entities)) {
    return(entities)
  }
  if (!"locus_tag" %in% names(entities) ||
      is.null(memberships) || !is.data.frame(memberships) ||
      !nrow(memberships) ||
      !all(c("target_key", "draw") %in% names(memberships))) {
    return(entities[0, , drop = FALSE])
  }
  drawn_targets <- unique(as.character(
    memberships$target_key[memberships$draw]
  ))
  entities[
    as.character(entities$locus_tag) %in% drawn_targets,
    , drop = FALSE
  ]
}

.dnmb_cct_layout_transporter_labels <- function(labels, xlim, y_memb = 8.5,
                                                base_offset = 0.40,
                                                track_step = 0.16,
                                                max_aux_tracks = 2L,
                                                min_gap = 0.08) {
  if (is.null(labels) || !nrow(labels)) return(labels)
  out <- as.data.frame(labels, stringsAsFactors = FALSE, row.names = NULL)
  if (!"anchor_x" %in% names(out)) out$anchor_x <- out$x
  if (!"anchor_y" %in% names(out)) out$anchor_y <- out$y
  if (!"priority" %in% names(out)) out$priority <- 0
  if (!"hjust" %in% names(out)) out$hjust <- 0.5

  xlim <- sort(as.numeric(xlim)[seq_len(2L)])
  if (any(!is.finite(xlim)) || xlim[1] >= xlim[2]) {
    stop("transporter label layout requires a finite increasing xlim", call. = FALSE)
  }
  base_offset <- max(0.20, as.numeric(base_offset)[1])
  track_step <- max(0.08, as.numeric(track_step)[1])
  max_aux_tracks <- max(0L, as.integer(max_aux_tracks)[1])
  min_gap <- max(0, as.numeric(min_gap)[1])

  out$label_width <- pmax(
    0.20,
    vapply(strsplit(as.character(out$label), "\n", fixed = TRUE), function(parts) {
      max(nchar(parts), na.rm = TRUE) * 0.06
    }, numeric(1))
  )
  out$x_lab <- pmin(
    pmax(out$anchor_x, xlim[1] + out$label_width / 2),
    xlim[2] - out$label_width / 2
  )
  out$track_level <- 0L

  # Assign the primary track first. A label is raised only when its anchored
  # interval collides with a label already occupying a lower track.
  occupied <- replicate(max_aux_tracks + 1L, integer(0), simplify = FALSE)
  ord <- order(
    out$anchor_x, -out$priority, as.character(out$label),
    na.last = TRUE
  )
  for (ii in ord) {
    collision_count <- vapply(seq_along(occupied), function(track) {
      idx <- occupied[[track]]
      if (!length(idx)) return(0L)
      required <- (out$label_width[idx] + out$label_width[ii]) / 2 + min_gap
      sum(abs(out$x_lab[idx] - out$x_lab[ii]) < required)
    }, integer(1))
    free_track <- which(collision_count == 0L)[1]
    if (is.na(free_track)) free_track <- which.min(collision_count)
    out$track_level[ii] <- free_track - 1L
    occupied[[free_track]] <- c(occupied[[free_track]], ii)
  }

  # Pack each occupied track left-to-right while preserving anchor order.
  # This handles dense sites after all three tracks have been used.
  for (track in seq.int(0L, max_aux_tracks)) {
    idx <- which(out$track_level == track)
    if (!length(idx)) next
    idx <- idx[order(out$anchor_x[idx], out$label[idx])]
    widths <- out$label_width[idx]
    pos <- out$x_lab[idx]
    pos[1] <- max(pos[1], xlim[1] + widths[1] / 2)
    if (length(idx) > 1L) {
      for (jj in 2:length(idx)) {
        required <- (widths[jj - 1L] + widths[jj]) / 2 + min_gap
        pos[jj] <- max(pos[jj], pos[jj - 1L] + required)
      }
      overflow <- pos[length(pos)] + widths[length(widths)] / 2 - xlim[2]
      if (overflow > 0) pos <- pos - overflow
      for (jj in rev(seq_len(length(idx) - 1L))) {
        required <- (widths[jj] + widths[jj + 1L]) / 2 + min_gap
        pos[jj] <- min(pos[jj], pos[jj + 1L] - required)
      }
    }
    underflow <- xlim[1] - (pos[1] - widths[1] / 2)
    if (underflow > 0) pos <- pos + underflow
    out$x_lab[idx] <- pos
  }

  out$y_lab <- y_memb + base_offset + out$track_level * track_step
  out$hjust <- 0.5
  out$draw_leader <- abs(out$x_lab - out$anchor_x) > 0.04 |
    out$track_level > 0L
  out
}

.dnmb_cct_transporter_label_leaders <- function(labels, clearance = 0.05,
                                                glyph_clearance = 0.11) {
  empty <- data.frame(
    group = character(), x = numeric(), y = numeric(), color = character(),
    stringsAsFactors = FALSE
  )
  required <- c("anchor_x", "anchor_y", "x_lab", "y_lab")
  if (is.null(labels) || !nrow(labels) || !all(required %in% names(labels))) {
    return(empty)
  }
  rows <- labels[
    is.finite(labels$anchor_x) & is.finite(labels$anchor_y) &
      is.finite(labels$x_lab) & is.finite(labels$y_lab),
    , drop = FALSE
  ]
  if ("draw_leader" %in% names(rows)) {
    keep <- !is.na(rows$draw_leader) & rows$draw_leader
    rows <- rows[keep, , drop = FALSE]
  }
  if (!nrow(rows)) return(empty)
  if (!"color" %in% names(rows)) rows$color <- "#52616B"

  dplyr::bind_rows(lapply(seq_len(nrow(rows)), function(i) {
    y_start <- rows$anchor_y[i] + glyph_clearance
    y_elbow <- max(y_start, rows$y_lab[i] - clearance)
    data.frame(
      group = paste0("transporter-label-", i),
      x = c(rows$anchor_x[i], rows$anchor_x[i], rows$x_lab[i]),
      y = c(y_start, y_elbow, y_elbow),
      color = rep(as.character(rows$color[i]), 3L),
      stringsAsFactors = FALSE
    )
  }))
}

.dnmb_cct_transporter_text_layer <- function(labels) {
  if (is.null(labels) || !nrow(labels)) return(NULL)
  ggplot2::geom_text(
    data = labels,
    ggplot2::aes(x = .data$x_lab, y = .data$y_lab, label = .data$label),
    color = labels$color,
    size = labels$size,
    fontface = labels$fontface,
    hjust = labels$hjust,
    lineheight = 0.9,
    check_overlap = FALSE,
    inherit.aes = FALSE
  )
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

.dnmb_cct_layout_external_labels <- function(labels, xlim,
                                              row_step = 0.15,
                                              max_aux_rows = 2L,
                                              min_gap = 0.08) {
  if (is.null(labels) || !nrow(labels)) return(labels)
  out <- as.data.frame(labels, stringsAsFactors = FALSE, row.names = NULL)
  if (!"label_kind" %in% names(out)) out$label_kind <- "fixed"
  if (!"anchor_x" %in% names(out)) out$anchor_x <- out$x
  if (!"anchor_y" %in% names(out)) out$anchor_y <- out$y
  if (!"priority" %in% names(out)) out$priority <- 0
  if (!"hjust" %in% names(out)) out$hjust <- 0.5
  if (!"color" %in% names(out)) out$color <- "#455A64"

  xlim <- sort(as.numeric(xlim)[seq_len(2L)])
  if (any(!is.finite(xlim)) || xlim[1] >= xlim[2]) {
    stop("external label layout requires a finite increasing xlim", call. = FALSE)
  }
  row_step <- max(0.05, as.numeric(row_step)[1])
  max_aux_rows <- max(0L, as.integer(max_aux_rows)[1])
  min_gap <- max(0, as.numeric(min_gap)[1])

  out$x_lab <- out$x
  out$y_lab <- out$y
  out$row_level <- 0L
  out$draw_leader <- FALSE
  out$hjust <- ifelse(is.na(out$hjust), 0.5, out$hjust)

  gh_idx <- which(out$label_kind == "gh" &
                    is.finite(out$anchor_x) & is.finite(out$anchor_y) &
                    is.finite(out$x) & is.finite(out$y))
  if (!length(gh_idx)) return(out)

  gh <- out[gh_idx, , drop = FALSE]
  gh$label_width <- pmax(
    0.20,
    vapply(strsplit(as.character(gh$label), "\n", fixed = TRUE), function(parts) {
      max(nchar(parts), na.rm = TRUE) * 0.055
    }, numeric(1))
  )

  # Use a single global row grid. Labels stay on their primary row whenever
  # their horizontal intervals fit; only a collision opens one of two lower
  # auxiliary rows.
  grid_origin <- max(gh$y, na.rm = TRUE)
  gh$primary_y <- grid_origin - round((grid_origin - gh$y) / row_step) * row_step
  gh$x_lab <- pmin(
    pmax(gh$anchor_x, xlim[1] + gh$label_width / 2),
    xlim[2] - gh$label_width / 2
  )
  gh$y_lab <- gh$primary_y
  gh$row_level <- 0L

  placed <- data.frame(
    x = numeric(), y = numeric(), width = numeric(),
    stringsAsFactors = FALSE
  )
  horizontal_step <- max(0.24, max(gh$label_width, na.rm = TRUE) + min_gap)
  horizontal_offsets <- c(0, horizontal_step, -horizontal_step)
  ord <- order(
    -gh$priority, -gh$primary_y, gh$anchor_x,
    as.character(gh$label), na.last = TRUE
  )
  for (ii in ord) {
    candidate_level <- rep(
      seq.int(0L, max_aux_rows),
      each = length(horizontal_offsets)
    )
    candidate_x <- gh$anchor_x[ii] + rep(
      horizontal_offsets,
      times = max_aux_rows + 1L
    )
    candidate_y <- gh$primary_y[ii] - candidate_level * row_step
    within_bounds <-
      candidate_x - gh$label_width[ii] / 2 >= xlim[1] &
      candidate_x + gh$label_width[ii] / 2 <= xlim[2]
    collision_count <- vapply(seq_along(candidate_x), function(candidate_i) {
      if (!within_bounds[candidate_i]) return(.Machine$integer.max)
      cx <- candidate_x[candidate_i]
      cy <- candidate_y[candidate_i]
      same_row <- abs(placed$y - cy) < row_step * 0.45
      if (!any(same_row)) return(0L)
      required <- (placed$width[same_row] + gh$label_width[ii]) / 2 + min_gap
      sum(abs(placed$x[same_row] - cx) < required)
    }, integer(1))
    chosen <- which(collision_count == 0L)[1]
    if (is.na(chosen)) {
      valid <- which(within_bounds)
      if (!length(valid)) valid <- seq_along(candidate_x)
      penalty <- collision_count[valid] * 100L +
        candidate_level[valid] * 10L +
        rep(c(0L, 1L, 1L), times = max_aux_rows + 1L)[valid]
      chosen <- valid[which.min(penalty)]
    }
    gh$x_lab[ii] <- pmin(
      pmax(candidate_x[chosen], xlim[1] + gh$label_width[ii] / 2),
      xlim[2] - gh$label_width[ii] / 2
    )
    gh$row_level[ii] <- candidate_level[chosen]
    gh$y_lab[ii] <- candidate_y[chosen]
    placed <- rbind(placed, data.frame(
      x = gh$x_lab[ii], y = gh$y_lab[ii], width = gh$label_width[ii],
      stringsAsFactors = FALSE
    ))
  }

  # Pack each occupied row left-to-right without changing anchor order. This
  # resolves the rare case where more than three labels share one cut site.
  row_key <- sprintf("%.5f", gh$y_lab)
  for (rk in unique(row_key)) {
    idx <- which(row_key == rk)
    idx <- idx[order(gh$x_lab[idx], gh$anchor_x[idx], gh$label[idx])]
    if (length(idx) <= 1L) next
    pos <- gh$x_lab[idx]
    widths <- gh$label_width[idx]
    pos[1] <- max(pos[1], xlim[1] + widths[1] / 2)
    for (jj in 2:length(idx)) {
      required <- (widths[jj - 1L] + widths[jj]) / 2 + min_gap
      pos[jj] <- max(pos[jj], pos[jj - 1L] + required)
    }
    overflow <- pos[length(pos)] + widths[length(widths)] / 2 - xlim[2]
    if (overflow > 0) pos <- pos - overflow
    for (jj in rev(seq_len(length(idx) - 1L))) {
      required <- (widths[jj] + widths[jj + 1L]) / 2 + min_gap
      pos[jj] <- min(pos[jj], pos[jj + 1L] - required)
    }
    underflow <- xlim[1] - (pos[1] - widths[1] / 2)
    if (underflow > 0) pos <- pos + underflow
    gh$x_lab[idx] <- pos
  }

  gh$draw_leader <- gh$row_level > 0L | abs(gh$x_lab - gh$anchor_x) > 0.04
  gh$hjust <- 0.5
  out[gh_idx, c("x_lab", "y_lab", "row_level", "draw_leader", "hjust")] <-
    gh[, c("x_lab", "y_lab", "row_level", "draw_leader", "hjust")]
  out
}

.dnmb_cct_external_label_leaders <- function(labels, label_clearance = 0.045) {
  empty <- data.frame(
    group = character(), x = numeric(), y = numeric(), color = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(labels) || !nrow(labels) || !"draw_leader" %in% names(labels)) {
    return(empty)
  }
  rows <- labels[
    labels$label_kind == "gh" & labels$draw_leader &
      is.finite(labels$anchor_x) & is.finite(labels$anchor_y) &
      is.finite(labels$x_lab) & is.finite(labels$y_lab),
    , drop = FALSE
  ]
  if (!nrow(rows)) return(empty)

  dplyr::bind_rows(lapply(seq_len(nrow(rows)), function(i) {
    direction <- sign(rows$anchor_y[i] - rows$y_lab[i])
    if (!is.finite(direction) || direction == 0) direction <- 1
    y_start <- rows$anchor_y[i] - direction * label_clearance
    y_elbow <- rows$y_lab[i] + direction * label_clearance
    data.frame(
      group = paste0("external-label-", i),
      x = c(rows$anchor_x[i], rows$anchor_x[i], rows$x_lab[i]),
      y = c(y_start, y_elbow, y_elbow),
      color = rep(as.character(rows$color[i]), 3L),
      stringsAsFactors = FALSE
    )
  }))
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
                                     grid_step = 0.5,
                                     rounded = TRUE) {
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

  waypoints <- data.frame(
    x = c(hub_x, hub_x, approach_x, approach_x, slot_x, slot_x, target_x),
    y = c(hub_y, corridor_y, corridor_y, approach_y, approach_y, target_y, target_y)
  )
  waypoints <- waypoints[
    c(TRUE, diff(waypoints$x) != 0 | diff(waypoints$y) != 0),
    , drop = FALSE
  ]
  if (!isTRUE(rounded)) return(waypoints)
  .dnmb_cct_waypoint_route_points(waypoints, grid_step = grid_step)
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

.dnmb_cct_initial_entry_targets <- function(cs_id, node_ids = NULL) {
  cid <- tolower(trimws(as.character(cs_id)[1L]))
  branched_targets <- list(
    lactose = c("Glucose", "Gal-1-P"),
    sucrose = c("Glucose", "Fru-6-P")
  )
  targets <- branched_targets[[cid]]
  if (is.null(targets)) {
    route_nodes <- .dnmb_cct_entry_route_nodes(cid, node_ids = node_ids)
    targets <- if (length(route_nodes)) route_nodes[1L] else character(0)
  }
  targets <- unique(as.character(targets))
  targets <- targets[!is.na(targets) & nzchar(targets)]
  if (!is.null(node_ids)) targets <- targets[targets %in% node_ids]
  targets
}

.dnmb_cct_points_from_start_to_nodes <- function(start_x, start_y, node_ids, node_x, node_y,
                                                 grid_step = 0.5,
                                                 lane_rank = 1L, lane_count = 1L,
                                                 obstacle_x = node_x,
                                                 obstacle_y = node_y) {
  if (length(node_ids) == 0) return(NULL)
  keep_ids <- node_ids[node_ids %in% names(node_x) & node_ids %in% names(node_y)]
  if (length(keep_ids) == 0) return(NULL)
  # Keep every biochemical node coordinate exact. Shared horizontal portions
  # are separated while these routes are still orthogonal waypoints; corners
  # are rounded only after the lane heights have been finalized.
  y_spread <- 0
  cur_x <- start_x
  cur_y <- start_y
  out <- data.frame(x = cur_x, y = cur_y)
  for (nid in keep_ids) {
    tgt_y <- unname(node_y[nid]) + y_spread
    seg <- .dnmb_cct_safe_pair_route_points(
      x1 = cur_x, y1 = cur_y,
      x2 = unname(node_x[nid]), y2 = tgt_y,
      obstacle_x = obstacle_x, obstacle_y = obstacle_y,
      grid_step = grid_step,
      rounded = FALSE
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

.dnmb_cct_route_label_point <- function(points_df, frac = 0.50) {
  if (is.null(points_df) || !is.data.frame(points_df) || !nrow(points_df) ||
      !all(c("x", "y") %in% names(points_df))) {
    return(data.frame(x = NA_real_, y = NA_real_))
  }
  pts <- points_df[is.finite(points_df$x) & is.finite(points_df$y), c("x", "y"), drop = FALSE]
  if (!nrow(pts)) return(data.frame(x = NA_real_, y = NA_real_))
  if (nrow(pts) == 1L) return(pts)
  keep <- c(TRUE, diff(pts$x) != 0 | diff(pts$y) != 0)
  pts <- pts[keep, , drop = FALSE]
  if (nrow(pts) == 1L) return(pts)

  seg_len <- sqrt(diff(pts$x)^2 + diff(pts$y)^2)
  total_len <- sum(seg_len)
  if (!is.finite(total_len) || total_len <= 0) return(pts[1, , drop = FALSE])
  frac <- min(1, max(0, as.numeric(frac)[1]))
  if (!is.finite(frac)) frac <- 0.5
  target_len <- total_len * frac
  cumulative <- c(0, cumsum(seg_len))
  seg_idx <- which(cumulative[-1L] >= target_len)[1]
  if (is.na(seg_idx)) return(pts[nrow(pts), , drop = FALSE])
  local_frac <- (target_len - cumulative[seg_idx]) / seg_len[seg_idx]
  data.frame(
    x = pts$x[seg_idx] + local_frac * (pts$x[seg_idx + 1L] - pts$x[seg_idx]),
    y = pts$y[seg_idx] + local_frac * (pts$y[seg_idx + 1L] - pts$y[seg_idx])
  )
}

.dnmb_cct_edge_points_from_row <- function(ce, edge_idx = 1L, grid_step = 0.5) {
  x1 <- ce$x; y1 <- ce$y; x2 <- ce$xend; y2 <- ce$yend
  if (any(is.na(c(x1, y1, x2, y2)))) return(data.frame(x = numeric(), y = numeric()))
  .dnmb_cct_pair_route_points(x1, y1, x2, y2, grid_step = grid_step)
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
                                          label_w = 0.38,
                                          label_h = 0.14,
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
  max(0.3, max(nchar(lines), na.rm = TRUE) * 0.08)
}

.dnmb_cct_final_text_layer <- function(labels, obstacles = NULL,
                                       xlim = c(NA, NA), ylim = c(NA, NA),
                                       seed = 42L, force = 1.1,
                                       force_pull = 2.2) {
  if (is.null(labels) || !nrow(labels)) return(NULL)
  labels <- as.data.frame(labels, stringsAsFactors = FALSE, row.names = NULL)
  defaults <- list(
    color = "#222222", size = 1, fontface = "plain", hjust = 0.5,
    priority = 0, nudge_x = 0, nudge_y = 0, alpha = 1
  )
  for (column in names(defaults)) {
    if (!column %in% names(labels)) labels[[column]] <- defaults[[column]]
  }
  labels <- labels[!is.na(labels$label) & nzchar(labels$label), , drop = FALSE]
  if (!nrow(labels)) return(NULL)
  labels <- labels[order(labels$priority, decreasing = TRUE, na.last = TRUE), , drop = FALSE]

  plot_df <- labels[, c(
    "x", "y", "label", "color", "size", "fontface", "hjust", "alpha",
    "nudge_x", "nudge_y"
  ), drop = FALSE]
  if (!is.null(obstacles) && nrow(obstacles)) {
    obstacles <- unique(as.data.frame(obstacles[, c("x", "y"), drop = FALSE]))
    obstacles <- obstacles[is.finite(obstacles$x) & is.finite(obstacles$y), , drop = FALSE]
    if (nrow(obstacles)) {
      obstacle_rows <- data.frame(
        x = obstacles$x, y = obstacles$y, label = "",
        color = "#00000000", size = 0.01, fontface = "plain", hjust = 0.5,
        nudge_x = 0, nudge_y = 0, alpha = 0,
        stringsAsFactors = FALSE
      )
      plot_df <- rbind(plot_df, obstacle_rows)
    }
  }

  ggrepel::geom_text_repel(
    data = plot_df,
    ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
    color = plot_df$color,
    size = plot_df$size,
    alpha = plot_df$alpha,
    fontface = plot_df$fontface,
    hjust = plot_df$hjust,
    nudge_x = plot_df$nudge_x,
    nudge_y = plot_df$nudge_y,
    bg.color = "#FFFFFF",
    bg.r = 0.045,
    box.padding = grid::unit(0.08, "lines"),
    point.padding = grid::unit(0.10, "lines"),
    min.segment.length = 0.04,
    segment.size = 0.12,
    segment.color = "#B0BEC5",
    segment.alpha = 0.65,
    force = force,
    force_pull = force_pull,
    max.overlaps = Inf,
    max.time = 3,
    max.iter = 15000,
    xlim = xlim,
    ylim = ylim,
    seed = seed,
    inherit.aes = FALSE
  )
}

.dnmb_cct_obstacle_perimeter <- function(x, y, rx = 0.1, ry = rx,
                                         points = 8L, include_center = TRUE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  rx <- rep_len(pmax(0, as.numeric(rx)), length(x))
  ry <- rep_len(pmax(0, as.numeric(ry)), length(x))
  keep <- is.finite(x) & is.finite(y) & is.finite(rx) & is.finite(ry)
  x <- x[keep]
  y <- y[keep]
  rx <- rx[keep]
  ry <- ry[keep]
  if (!length(x)) return(data.frame(x = numeric(), y = numeric()))

  theta <- seq(0, 2 * pi, length.out = max(4L, as.integer(points)[1]) + 1L)
  theta <- theta[-length(theta)]
  perimeter <- do.call(rbind, lapply(seq_along(x), function(i) {
    data.frame(
      x = x[i] + rx[i] * cos(theta),
      y = y[i] + ry[i] * sin(theta)
    )
  }))
  if (isTRUE(include_center)) {
    perimeter <- rbind(data.frame(x = x, y = y), perimeter)
  }
  rownames(perimeter) <- NULL
  perimeter
}

.dnmb_cct_collapse_route_labels <- function(route_labels) {
  if (is.null(route_labels) || !nrow(route_labels)) return(route_labels)
  route_labels <- as.data.frame(route_labels, stringsAsFactors = FALSE, row.names = NULL)
  required <- c("label", "locus_tag", "member_keys", "target_id", "priority")
  if (!all(required %in% names(route_labels))) {
    stop("route label deduplication requires membership and target columns", call. = FALSE)
  }
  exact_key <- paste(
    route_labels$target_id,
    route_labels$locus_tag,
    route_labels$member_keys,
    route_labels$label,
    sep = "\r"
  )
  ord <- order(-route_labels$priority)
  route_labels <- route_labels[ord, , drop = FALSE]
  exact_key <- exact_key[ord]
  route_labels[!duplicated(exact_key), , drop = FALSE]
}

.dnmb_cct_pathway_presence <- function(step_status, pstats = NULL,
                                       valid_pathways = NULL) {
  empty <- data.frame(
    pathway_id = character(), fraction = numeric(),
    n_transport = integer(), n_intracellular = integer(),
    cytoplasm_status = character(), stringsAsFactors = FALSE
  )
  if (is.null(step_status) || !nrow(step_status)) return(empty)

  steps <- as.data.frame(step_status, stringsAsFactors = FALSE, row.names = NULL)
  required <- c("pathway_id", "step_id", "confidence", "locus_tag")
  if (!all(required %in% names(steps))) return(empty)
  steps$pathway_id <- tolower(as.character(steps$pathway_id))
  if (!is.null(valid_pathways)) {
    steps <- steps[steps$pathway_id %in% tolower(valid_pathways), , drop = FALSE]
  }
  if (!nrow(steps)) return(empty)

  conf_rank <- c(none = 0L, low = 1L, medium = 2L, high = 3L)
  steps$rank <- unname(conf_rank[tolower(as.character(steps$confidence))])
  steps$rank[is.na(steps$rank)] <- 0L
  steps <- steps[
    steps$rank >= 2L & !is.na(steps$locus_tag) & nzchar(steps$locus_tag),
    , drop = FALSE
  ]

  path_ids <- unique(c(
    if (!is.null(valid_pathways)) tolower(valid_pathways) else character(),
    steps$pathway_id,
    if (!is.null(pstats) && "pathway_id" %in% names(pstats)) {
      tolower(as.character(pstats$pathway_id))
    } else character()
  ))
  path_ids <- path_ids[!is.na(path_ids) & nzchar(path_ids)]
  if (!length(path_ids)) return(empty)

  out <- data.frame(
    pathway_id = path_ids,
    fraction = 0,
    n_transport = 0L,
    n_intracellular = 0L,
    stringsAsFactors = FALSE
  )
  if (!is.null(pstats) && nrow(pstats) &&
      all(c("pathway_id", "fraction") %in% names(pstats))) {
    ps <- as.data.frame(pstats, stringsAsFactors = FALSE, row.names = NULL)
    ps$pathway_id <- tolower(as.character(ps$pathway_id))
    idx <- match(out$pathway_id, ps$pathway_id)
    found <- !is.na(idx)
    out$fraction[found] <- suppressWarnings(as.numeric(ps$fraction[idx[found]]))
    out$fraction[!is.finite(out$fraction)] <- 0
  }

  if (nrow(steps)) {
    steps$is_transport <- mapply(
      .dnmb_cct_is_transport_like_step,
      step_id = steps$step_id,
      gene_name = steps$step_id,
      MoreArgs = list(is_pts = FALSE)
    )
    count_distinct <- function(df, keep_transport) {
      sub <- df[df$is_transport == keep_transport, , drop = FALSE]
      if (!nrow(sub)) return(integer())
      vapply(split(sub$step_id, sub$pathway_id), function(x) {
        length(unique(as.character(x)))
      }, integer(1))
    }
    tr_count <- count_distinct(steps, TRUE)
    in_count <- count_distinct(steps, FALSE)
    if (length(tr_count)) {
      out$n_transport[match(names(tr_count), out$pathway_id)] <- unname(tr_count)
    }
    if (length(in_count)) {
      out$n_intracellular[match(names(in_count), out$pathway_id)] <- unname(in_count)
    }
  }

  active <- (out$fraction >= 0.50 & out$n_intracellular >= 2L) |
    (out$fraction >= 0.25 & out$n_transport >= 1L & out$n_intracellular >= 1L)
  partial <- !active & (out$n_transport > 0L | out$n_intracellular > 0L |
                          out$fraction > 0)
  out$cytoplasm_status <- ifelse(active, "active", ifelse(partial, "partial", "reference"))
  out
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

.dnmb_cct_path_step_key <- function(pathway_id, step_id) {
  pathway_id <- tolower(trimws(as.character(pathway_id)))
  step_id <- tolower(trimws(as.character(step_id)))
  paste(pathway_id, step_id, sep = "\r")
}

# Only reactions with an explicit directed biochemical edge are eligible for a
# map label. Evidence without one of these mappings remains in the ledger.
.dnmb_cct_exact_step_edge_map <- function() {
  data.frame(
    pathway_id = rep("deoxyribose", 8L),
    step_id = c("deoK", "deoC", "deoC", "adh", "ackA", "pta", "acs", "ald-dh-CoA"),
    from = c(
      "Deoxyribose", "Deoxyribose-5-P", "Deoxyribose-5-P",
      "Acetaldehyde", "Acetate", "Acetyl-P", "Acetate", "Acetaldehyde"
    ),
    to = c(
      "Deoxyribose-5-P", "GA3P", "Acetaldehyde", "Acetate",
      "Acetyl-P", "Acetyl-CoA", "Acetyl-CoA", "Acetyl-CoA"
    ),
    label_anchor = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

.dnmb_cct_exact_step_edge_matches <- function(step_evidence, cyto_edges,
                                               edge_map = .dnmb_cct_exact_step_edge_map()) {
  empty <- data.frame()
  if (is.null(step_evidence) || !is.data.frame(step_evidence) || !nrow(step_evidence) ||
      is.null(cyto_edges) || !is.data.frame(cyto_edges) || !nrow(cyto_edges) ||
      is.null(edge_map) || !is.data.frame(edge_map) || !nrow(edge_map)) {
    return(empty)
  }
  if (!all(c("pathway_id", "step_id") %in% names(step_evidence)) ||
      !all(c("from", "to") %in% names(cyto_edges)) ||
      !all(c("pathway_id", "step_id", "from", "to", "label_anchor") %in% names(edge_map))) {
    return(empty)
  }

  evidence <- as.data.frame(step_evidence, stringsAsFactors = FALSE, row.names = NULL)
  evidence$.evidence_order <- seq_len(nrow(evidence))
  evidence$.path_step_key <- .dnmb_cct_path_step_key(
    evidence$pathway_id, evidence$step_id
  )
  mapping <- as.data.frame(edge_map, stringsAsFactors = FALSE, row.names = NULL)
  mapping$.map_order <- seq_len(nrow(mapping))
  mapping$.path_step_key <- .dnmb_cct_path_step_key(
    mapping$pathway_id, mapping$step_id
  )
  names(mapping)[names(mapping) %in% c("pathway_id", "step_id", "from", "to")] <-
    c("mapped_pathway_id", "mapped_step_id", "mapped_from", "mapped_to")

  matched <- merge(
    evidence, mapping,
    by = ".path_step_key", all = FALSE, sort = FALSE
  )
  if (!nrow(matched)) return(empty)

  edge_key <- paste(cyto_edges$from, cyto_edges$to, sep = "\r")
  matched$edge_index <- match(
    paste(matched$mapped_from, matched$mapped_to, sep = "\r"),
    edge_key
  )
  matched <- matched[!is.na(matched$edge_index), , drop = FALSE]
  if (!nrow(matched)) return(empty)
  matched <- matched[order(matched$.evidence_order, matched$.map_order), , drop = FALSE]
  rownames(matched) <- NULL
  matched
}

.dnmb_cct_exact_step_labels <- function(step_evidence, cyto_edges,
                                        edge_map = .dnmb_cct_exact_step_edge_map()) {
  empty <- data.frame(
    x = numeric(), y = numeric(), label = character(), base_label = character(),
    step_id = character(), pathway_id = character(), confidence = character(),
    rank = integer(), locus_tag = character(), member_loci = character(),
    member_keys = character(), color = character(), target_id = character(),
    target_rank = integer(), priority = numeric(), angle = numeric(),
    nudge_x = numeric(), nudge_y = numeric(), size = numeric(),
    fontface = character(), hjust = numeric(), alpha = numeric(),
    stringsAsFactors = FALSE
  )
  matches <- .dnmb_cct_exact_step_edge_matches(
    step_evidence = step_evidence,
    cyto_edges = cyto_edges,
    edge_map = edge_map
  )
  if (!nrow(matches)) return(empty)

  if (!"rank" %in% names(matches)) {
    confidence_rank <- c(none = 0L, low = 1L, medium = 2L, high = 3L)
    matches$rank <- unname(confidence_rank[tolower(as.character(matches$confidence))])
    matches$rank[is.na(matches$rank)] <- 0L
  }
  if (!"confidence" %in% names(matches)) {
    matches$confidence <- ifelse(matches$rank >= 3L, "high", "medium")
  }
  if (!"locus_tag" %in% names(matches)) matches$locus_tag <- ""

  split_rows <- split(
    seq_len(nrow(matches)),
    .dnmb_cct_path_step_key(matches$mapped_pathway_id, matches$mapped_step_id)
  )
  rows <- lapply(split_rows, function(idx) {
    group <- matches[idx, , drop = FALSE]
    anchor <- group[group$label_anchor %in% TRUE, , drop = FALSE]
    if (!nrow(anchor)) anchor <- group[1, , drop = FALSE]
    anchor <- anchor[1, , drop = FALSE]
    edge <- cyto_edges[anchor$edge_index, , drop = FALSE]

    loci <- sort(unique(trimws(as.character(group$locus_tag))))
    loci <- loci[!is.na(loci) & nzchar(loci)]
    if (!length(loci)) return(NULL)
    rank <- max(suppressWarnings(as.integer(group$rank)), na.rm = TRUE)
    if (!is.finite(rank) || rank < 2L) return(NULL)
    confidence <- if (rank >= 3L) "high" else "medium"
    step_label <- as.character(anchor$mapped_step_id)
    label <- paste0(step_label, "\n", paste(loci, collapse = " | "))
    edge_points <- .dnmb_cct_edge_points_from_row(
      edge, edge_idx = anchor$edge_index, grid_step = 0.5
    )
    label_point <- .dnmb_cct_route_label_point(edge_points, frac = 0.5)
    x <- label_point$x[1]
    y <- label_point$y[1]
    is_split <- identical(tolower(step_label), "deoc")

    data.frame(
      x = x, y = y, label = label, base_label = label,
      step_id = step_label,
      pathway_id = as.character(anchor$mapped_pathway_id),
      confidence = confidence, rank = rank,
      locus_tag = paste(loci, collapse = ","),
      member_loci = paste(loci, collapse = ";"),
      member_keys = paste(
        paste0(tolower(step_label), "::", loci), collapse = ";"
      ),
      color = if (rank >= 3L) "#007C83" else "#B7791F",
      target_id = paste(unique(group$mapped_to), collapse = ";"),
      target_rank = 1L, priority = 24 + rank, angle = 0,
      nudge_x = if (is_split) 0.24 else 0.12,
      nudge_y = if (is_split) -0.08 else 0,
      size = 1.14, fontface = "bold", hjust = 0,
      alpha = 1,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(empty)
  out <- dplyr::bind_rows(rows)
  out[order(out$pathway_id, out$step_id), , drop = FALSE]
}

.dnmb_cct_step_target_nodes <- function(step_id, pathway_id = NA_character_) {
  sid <- tolower(trimws(as.character(step_id)[1]))
  pid <- tolower(trimws(as.character(pathway_id)[1]))
  if (is.na(sid) || !nzchar(sid)) return(character(0))
  entry_map <- .dnmb_cct_auto_entry_map()
  entry_node <- unname(entry_map[pid])

  is_transport_like <- grepl(
    "pts|transport|permease|porter|abc|binding|mfs|symporter|tm00|ggu|mgl|rbs|xyl[fgh]|lac[efg]|mal[efgk]|thu[efgk]|ara[etuv]|galp|glc[ptuv]|kgut|nup[abc]|pot[abcd]|mtl[ek]|chve",
    sid, ignore.case = TRUE
  )
  if (is_transport_like && !is.na(entry_node) && nzchar(entry_node)) {
    return(c(entry_node))
  }

  if (sid %in% c("glk", "mgla", "mglb", "mglc")) return(c("Glc-6-P"))
  if (sid %in% c("naga")) return(c("GlcN-6-P"))
  if (sid %in% c("nagb", "mana")) return(c("Fru-6-P"))
  if (sid %in% c("galk")) return(c("Gal-1-P"))
  if (sid %in% c("galt")) return(c("UDP-Gal", "Glc-1-P", "UDP-Glc"))
  if (sid %in% c("gale")) return(c("UDP-Gal", "UDP-Glc"))
  if (sid %in% c("pgma")) return(c("Glc-1-P", "Glc-6-P"))
  if (sid %in% c("gntk", "gnd")) return(c("6-PG"))
  if (sid %in% c("kdgk", "eda")) return(c("KDPG", "Pyruvate", "GA3P"))
  if (sid %in% c("xyla")) return(c("Xylulose", "Xu-5-P"))
  if (sid %in% c("xylb")) return(c("Xu-5-P"))
  if (sid %in% c("rbsk")) return(c("R-5-P"))
  if (sid %in% c("araa")) return(c("Ribulose", "Ribulose-5-P", "Xu-5-P"))
  if (sid %in% c("arab")) return(c("Ribulose-5-P", "Xu-5-P"))
  if (sid %in% c("arad", "ggua", "ggub")) return(c("Xu-5-P"))
  if (sid %in% c("glpk", "glpo")) return(c("Glycerol-3-P", "DHAP"))
  if (sid %in% c("fuci")) return(c("Fuculose", "Fuculose-1-P", "DHAP"))
  if (sid %in% c("fuck")) return(c("Fuculose-1-P", "DHAP"))
  if (sid %in% c("fba")) return(c("Fru-1,6-BP", "GA3P"))
  if (sid %in% c("tpi")) return(c("DHAP", "GA3P"))
  if (pid == "deoxyribose" && sid == "deok") return(c("Deoxyribose-5-P"))
  if (pid == "deoxyribose" && sid == "deoc") {
    return(c("Deoxyribose-5-P", "GA3P", "Acetaldehyde"))
  }
  if (pid == "deoxyribose" && sid == "adh") return(c("Acetaldehyde", "Acetate"))
  if (pid == "deoxyribose" && sid == "acka") return(c("Acetate", "Acetyl-P"))
  if (pid == "deoxyribose" && sid == "pta") return(c("Acetyl-P", "Acetyl-CoA"))
  if (pid == "deoxyribose" && sid == "acs") return(c("Acetate", "Acetyl-CoA"))
  if (pid == "deoxyribose" && sid == "ald-dh-coa") return(c("Acetaldehyde", "Acetyl-CoA"))
  if (sid %in% c("thuk", "pstp")) return(c("Trehalose-6-P", "Glc-6-P"))
  if (sid %in% c("lacz", "lace", "lacf", "lacg", "lack")) return(c("Gal-1-P", "Glc-6-P"))
  if (sid %in% c("ggua", "ggub", "chve", "xylg", "xylh", "rbsa", "rbsb", "rbsc", "nupa", "nupb", "nupc",
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
  if (pid == "deoxyribose") {
    return(c("Deoxyribose-5-P", "GA3P", "Acetaldehyde", "Acetate", "Acetyl-CoA"))
  }
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
    deoxyribose = c("Deoxyribose-5-P", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    gluconate = c("6-PG", "KDPG", "Pyruvate"),
    fucose = c("Fuculose", "Fuculose-1-P", "DHAP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate"),
    rhamnose = c("Rhamnulose", "Rhamnulose-1-P", "DHAP", "GA3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate")
  )
}

.dnmb_cct_core_merge_nodes <- function() {
  c(
    "Glc-6-P", "Fru-6-P", "Fru-1,6-BP", "DHAP", "GA3P",
    "1,3-BPG", "3-PG", "2-PG", "PEP", "Pyruvate", "Acetyl-CoA",
    "6-PGL", "6-PG", "Ru-5-P", "R-5-P", "Xu-5-P", "S-7-P", "E-4-P",
    "KDPG"
  )
}

.dnmb_cct_route_to_first_core_merge <- function(node_ids,
                                                merge_nodes = .dnmb_cct_core_merge_nodes()) {
  node_ids <- as.character(node_ids)
  node_ids <- node_ids[!is.na(node_ids) & nzchar(node_ids)]
  if (length(node_ids) < 2L || node_ids[1L] %in% merge_nodes) return(NULL)
  merge_index <- which(node_ids[-1L] %in% merge_nodes)
  if (!length(merge_index)) return(NULL)
  node_ids[seq_len(merge_index[1L] + 1L)]
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
  if (pid %in% c("gluconate", "xylose", "arabinose", "ribose", "deoxyribose")) {
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
  if (pid == "deoxyribose") {
    return(c("Deoxyribose-5-P", "GA3P", "Acetaldehyde", "Acetate", "Acetyl-CoA", "Pyruvate"))
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

  # The overlay must retain the pathway-specific entry chemistry. Starting at
  # any low-cost backbone/PPP node can otherwise skip named intermediates (for
  # example Glucose, Xylulose, or Deoxyribose-5-P).
  seed_route <- .dnmb_cct_continuity_routes()[[pid]]
  explicit_start <- if (length(seed_route)) seed_route[1L] else NA_character_
  candidate_starts <- unique(c(explicit_start, entry_target))
  candidate_starts <- candidate_starts[
    candidate_starts %in% names(node_type) &
      node_type[candidate_starts] %in% c("entry_intermediate", "backbone", "ppp", "ed")
  ]
  if (!is.na(explicit_start) && explicit_start %in% candidate_starts) {
    candidate_starts <- explicit_start
  }
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

.dnmb_cct_path_obstacle_score <- function(points_df, obstacle_x, obstacle_y,
                                          clearance = 0.14) {
  if (is.null(points_df) || nrow(points_df) < 2) return(Inf)
  obstacles <- data.frame(x = as.numeric(obstacle_x), y = as.numeric(obstacle_y))
  obstacles <- obstacles[is.finite(obstacles$x) & is.finite(obstacles$y), , drop = FALSE]
  if (!nrow(obstacles)) {
    return(sum(sqrt(diff(points_df$x)^2 + diff(points_df$y)^2)))
  }

  start_dist <- sqrt((obstacles$x - points_df$x[1])^2 +
                     (obstacles$y - points_df$y[1])^2)
  end_dist <- sqrt((obstacles$x - points_df$x[nrow(points_df)])^2 +
                   (obstacles$y - points_df$y[nrow(points_df)])^2)
  obstacles <- obstacles[
    start_dist > clearance * 0.85 & end_dist > clearance * 0.85,
    , drop = FALSE
  ]

  min_dist <- rep(Inf, nrow(obstacles))
  if (nrow(obstacles)) {
    for (i in seq_len(nrow(points_df) - 1L)) {
      ax <- points_df$x[i]
      ay <- points_df$y[i]
      bx <- points_df$x[i + 1L]
      by <- points_df$y[i + 1L]
      vx <- bx - ax
      vy <- by - ay
      denom <- vx^2 + vy^2
      if (denom <= 1e-12) next
      projection <- ((obstacles$x - ax) * vx + (obstacles$y - ay) * vy) / denom
      projection <- pmin(1, pmax(0, projection))
      px <- ax + projection * vx
      py <- ay + projection * vy
      distance <- sqrt((obstacles$x - px)^2 + (obstacles$y - py)^2)
      min_dist <- pmin(min_dist, distance)
    }
  }

  collision <- min_dist < clearance
  collision_penalty <- sum(ifelse(
    collision, 1000 + 250 * (clearance - min_dist), 0
  ))
  path_length <- sum(sqrt(diff(points_df$x)^2 + diff(points_df$y)^2))
  bend_penalty <- 0.04 * max(0, nrow(points_df) - 2L)
  collision_penalty + path_length + bend_penalty
}

# Choose an orthogonal hub route that avoids unrelated metabolite nodes.
.dnmb_cct_safe_pair_route_points <- function(x1, y1, x2, y2, obstacle_x,
                                             obstacle_y, grid_step = 0.5,
                                             clearance = 0.14,
                                             rounded = TRUE) {
  if (any(!is.finite(c(x1, y1, x2, y2)))) {
    return(data.frame(x = numeric(), y = numeric()))
  }
  clean_candidate <- function(points) {
    keep <- c(TRUE, diff(points$x) != 0 | diff(points$y) != 0)
    points[keep, , drop = FALSE]
  }
  candidates <- list(
    clean_candidate(data.frame(x = c(x1, x1, x2), y = c(y1, y2, y2))),
    clean_candidate(data.frame(x = c(x1, x2, x2), y = c(y1, y1, y2)))
  )
  if (abs(x2 - x1) < 0.01 || abs(y2 - y1) < 0.01) {
    candidates <- c(list(data.frame(x = c(x1, x2), y = c(y1, y2))), candidates)
  }

  y_direction <- if (y2 < y1) -1 else 1
  x_direction <- if (x2 < x1) -1 else 1
  corridor_y <- unique(c(
    y1 + y_direction * grid_step / 2,
    y2 - y_direction * grid_step / 2,
    y1 - y_direction * grid_step / 2,
    y2 + y_direction * grid_step / 2
  ))
  for (cy in corridor_y) {
    candidates[[length(candidates) + 1L]] <- clean_candidate(data.frame(
      x = c(x1, x1, x2, x2), y = c(y1, cy, cy, y2)
    ))
  }
  corridor_x <- unique(c(
    x1 + x_direction * grid_step / 2,
    x2 - x_direction * grid_step / 2,
    x1 - x_direction * grid_step / 2,
    x2 + x_direction * grid_step / 2
  ))
  for (cx in corridor_x) {
    candidates[[length(candidates) + 1L]] <- clean_candidate(data.frame(
      x = c(x1, cx, cx, x2), y = c(y1, y1, y2, y2)
    ))
  }

  scores <- vapply(
    candidates, .dnmb_cct_path_obstacle_score, numeric(1),
    obstacle_x = obstacle_x, obstacle_y = obstacle_y, clearance = clearance
  )
  best <- candidates[[which.min(scores)]]
  if (!isTRUE(rounded)) return(best)
  .dnmb_cct_rounded_route_points(best, radius = min(grid_step * 0.35, clearance * 0.9))
}

.dnmb_cct_raw_pair_route_points <- function(x1, y1, x2, y2,
                                            grid_step = 0.5) {
  if (any(!is.finite(c(x1, y1, x2, y2)))) {
    return(data.frame(x = numeric(), y = numeric()))
  }
  if (abs(x2 - x1) < 0.01 || abs(y2 - y1) < 0.01) {
    return(data.frame(x = c(x1, x2), y = c(y1, y2)))
  }

  # Keep a horizontal corridor between two vertical endpoint stems. Lane
  # offsets can then move the corridor without moving either metabolite node.
  corridor_y <- .dnmb_cct_snap_to_grid(
    (y1 + y2) / 2,
    step = max(0.01, grid_step / 2)
  )
  if (abs(corridor_y - y1) < 0.01 || abs(corridor_y - y2) < 0.01) {
    corridor_y <- (y1 + y2) / 2
  }
  data.frame(
    x = c(x1, x1, x2, x2),
    y = c(y1, corridor_y, corridor_y, y2)
  )
}

.dnmb_cct_route_overlay_edges <- function(node_ids, node_x, node_y,
                                          grid_step = 0.5) {
  keep_ids <- node_ids[node_ids %in% names(node_x) & node_ids %in% names(node_y)]
  if (length(keep_ids) < 2L) return(list())
  lapply(seq_len(length(keep_ids) - 1L), function(edge_index) {
    from_id <- keep_ids[edge_index]
    to_id <- keep_ids[edge_index + 1L]
    pts <- .dnmb_cct_raw_pair_route_points(
      x1 = unname(node_x[from_id]), y1 = unname(node_y[from_id]),
      x2 = unname(node_x[to_id]), y2 = unname(node_y[to_id]),
      grid_step = grid_step
    )
    attr(pts, "from_id") <- from_id
    attr(pts, "to_id") <- to_id
    pts
  })
}

.dnmb_cct_validate_extracellular_routes <- function(specs, anchors,
                                                     tolerance = 1e-6) {
  if (is.null(specs) || !length(specs)) return(list())
  if (!is.list(specs)) {
    stop("extracellular route specs must be a list", call. = FALSE)
  }
  required_anchor_cols <- c("kind", "id", "x", "y")
  if (is.null(anchors) || !is.data.frame(anchors) ||
      !all(required_anchor_cols %in% names(anchors))) {
    stop(
      "extracellular anchors require kind, id, x, and y columns",
      call. = FALSE
    )
  }
  anchors <- as.data.frame(anchors, stringsAsFactors = FALSE, row.names = NULL)
  anchors$kind <- trimws(as.character(anchors$kind))
  anchors$id <- trimws(as.character(anchors$id))
  anchors$x <- suppressWarnings(as.numeric(anchors$x))
  anchors$y <- suppressWarnings(as.numeric(anchors$y))
  anchors <- anchors[
    !is.na(anchors$kind) & nzchar(anchors$kind) &
      !is.na(anchors$id) & nzchar(anchors$id) &
      is.finite(anchors$x) & is.finite(anchors$y),
    , drop = FALSE
  ]
  tolerance <- max(0, suppressWarnings(as.numeric(tolerance)[1L]))
  if (!is.finite(tolerance)) tolerance <- 1e-6
  axis_tolerance <- max(1e-8, tolerance)

  anchor_for <- function(kind, id) {
    idx <- which(anchors$kind == kind & anchors$id == id)
    if (length(idx) != 1L) return(NULL)
    anchors[idx, , drop = FALSE]
  }
  validated <- list()
  for (spec in specs) {
    if (!is.list(spec)) next
    required_spec <- c(
      "source_kind", "source_id", "target_kind", "target_id", "points"
    )
    if (!all(required_spec %in% names(spec))) next
    source_kind <- trimws(as.character(spec$source_kind)[1L])
    source_id <- trimws(as.character(spec$source_id)[1L])
    target_kind <- trimws(as.character(spec$target_kind)[1L])
    target_id <- trimws(as.character(spec$target_id)[1L])
    if (any(is.na(c(source_kind, source_id, target_kind, target_id))) ||
        any(!nzchar(c(source_kind, source_id, target_kind, target_id)))) next
    source <- anchor_for(source_kind, source_id)
    target <- anchor_for(target_kind, target_id)
    if (is.null(source) || is.null(target)) next

    points <- spec$points
    if (is.null(points) || !is.data.frame(points) || nrow(points) < 2L ||
        !all(c("x", "y") %in% names(points))) next
    points <- as.data.frame(points[, c("x", "y"), drop = FALSE])
    points$x <- suppressWarnings(as.numeric(points$x))
    points$y <- suppressWarnings(as.numeric(points$y))
    if (any(!is.finite(points$x)) || any(!is.finite(points$y))) next
    source_distance <- sqrt(
      (points$x[1L] - source$x)^2 + (points$y[1L] - source$y)^2
    )
    target_distance <- sqrt(
      (points$x[nrow(points)] - target$x)^2 +
        (points$y[nrow(points)] - target$y)^2
    )
    if (source_distance > tolerance || target_distance > tolerance) next

    points$x[1L] <- source$x
    points$y[1L] <- source$y
    points$x[nrow(points)] <- target$x
    points$y[nrow(points)] <- target$y
    points <- points[
      c(TRUE, diff(points$x) != 0 | diff(points$y) != 0),
      , drop = FALSE
    ]
    if (nrow(points) < 2L) next
    dx <- diff(points$x)
    dy <- diff(points$y)
    if (any(abs(dx) > axis_tolerance & abs(dy) > axis_tolerance)) next
    spec$points <- points
    validated[[length(validated) + 1L]] <- spec
  }
  validated
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

.dnmb_cct_round_orthogonal_route <- function(points_df, radius = 0.04,
                                             n_arc = 5L,
                                             axis_tolerance = 1e-8) {
  if (is.null(points_df) || !is.data.frame(points_df) || nrow(points_df) < 3L) {
    return(points_df)
  }
  pts <- as.data.frame(points_df[, c("x", "y"), drop = FALSE])
  keep <- is.finite(pts$x) & is.finite(pts$y)
  pts <- pts[keep, , drop = FALSE]
  if (nrow(pts) < 3L) return(pts)
  pts <- pts[c(TRUE, diff(pts$x) != 0 | diff(pts$y) != 0), , drop = FALSE]
  if (nrow(pts) < 3L || !is.finite(radius) || radius <= 0) return(pts)

  out <- pts[1L, , drop = FALSE]
  for (i in 2:(nrow(pts) - 1L)) {
    prev <- c(pts$x[i - 1L], pts$y[i - 1L])
    corner <- c(pts$x[i], pts$y[i])
    next_pt <- c(pts$x[i + 1L], pts$y[i + 1L])
    vin <- corner - prev
    vout <- next_pt - corner
    len_in <- sqrt(sum(vin^2))
    len_out <- sqrt(sum(vout^2))
    if (len_in <= axis_tolerance || len_out <= axis_tolerance) next
    uin <- vin / len_in
    uout <- vout / len_out
    orthogonal <- abs(sum(uin * uout)) <= axis_tolerance &&
      ((abs(uin[1]) <= axis_tolerance || abs(uin[2]) <= axis_tolerance) &&
       (abs(uout[1]) <= axis_tolerance || abs(uout[2]) <= axis_tolerance))
    if (!orthogonal) {
      out <- rbind(out, pts[i, , drop = FALSE])
      next
    }
    r_eff <- min(radius, len_in * 0.30, len_out * 0.30)
    if (!is.finite(r_eff) || r_eff <= axis_tolerance) {
      out <- rbind(out, pts[i, , drop = FALSE])
      next
    }
    tangent_in <- corner - uin * r_eff
    arc <- .dnmb_cct_quarter_arc_points(
      corner = corner, dir_in = uin, dir_out = uout,
      radius = r_eff, n = max(3L, as.integer(n_arc)[1L])
    )
    out <- rbind(
      out,
      data.frame(x = tangent_in[1L], y = tangent_in[2L]),
      arc[-1L, , drop = FALSE]
    )
  }
  out <- rbind(out, pts[nrow(pts), , drop = FALSE])
  out[c(TRUE, diff(out$x) != 0 | diff(out$y) != 0), , drop = FALSE]
}

.dnmb_cct_route_lane_assignments <- function(routes, metadata = NULL,
                                             lane_step = 0.05,
                                             max_offset = 0.18,
                                             horizontal_tolerance = 0.01,
                                             y_tolerance = 0.04,
                                             min_overlap = 0.18) {
  if (is.null(routes)) routes <- list()
  if (!is.list(routes)) stop("route lanes require a list of polylines", call. = FALSE)
  n_routes <- length(routes)
  empty <- data.frame(
    route_index = integer(), route_id = character(), segment_index = integer(),
    group = character(), priority = numeric(), x_min = numeric(), x_max = numeric(),
    original_y = numeric(), lane_y = numeric(), offset = numeric(),
    component = integer(), stringsAsFactors = FALSE
  )
  if (!n_routes) return(empty)

  if (is.null(metadata)) metadata <- data.frame(route_id = as.character(seq_len(n_routes)))
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, row.names = NULL)
  if (nrow(metadata) != n_routes) {
    stop("route metadata must have one row per polyline", call. = FALSE)
  }
  if (!"route_id" %in% names(metadata)) metadata$route_id <- as.character(seq_len(n_routes))
  if (!"priority" %in% names(metadata)) metadata$priority <- 0
  if (!"group" %in% names(metadata)) metadata$group <- "route"
  metadata$route_id <- as.character(metadata$route_id)
  metadata$group <- as.character(metadata$group)
  metadata$priority <- suppressWarnings(as.numeric(metadata$priority))
  metadata$priority[!is.finite(metadata$priority)] <- 0
  if (anyNA(metadata$route_id) || any(!nzchar(metadata$route_id)) || anyDuplicated(metadata$route_id)) {
    stop("route_id values must be unique and non-empty", call. = FALSE)
  }

  segment_rows <- list()
  for (route_index in seq_len(n_routes)) {
    pts <- routes[[route_index]]
    if (is.null(pts) || !is.data.frame(pts) || nrow(pts) < 2L ||
        !all(c("x", "y") %in% names(pts))) next
    for (segment_index in seq_len(nrow(pts) - 1L)) {
      x1 <- suppressWarnings(as.numeric(pts$x[segment_index]))
      x2 <- suppressWarnings(as.numeric(pts$x[segment_index + 1L]))
      y1 <- suppressWarnings(as.numeric(pts$y[segment_index]))
      y2 <- suppressWarnings(as.numeric(pts$y[segment_index + 1L]))
      if (any(!is.finite(c(x1, x2, y1, y2)))) next
      if (abs(y2 - y1) > horizontal_tolerance || abs(x2 - x1) < min_overlap) next
      segment_rows[[length(segment_rows) + 1L]] <- data.frame(
        route_index = route_index,
        route_id = metadata$route_id[route_index],
        segment_index = segment_index,
        group = metadata$group[route_index],
        priority = metadata$priority[route_index],
        x_min = min(x1, x2), x_max = max(x1, x2),
        original_y = mean(c(y1, y2)),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(segment_rows)) return(empty)
  segments <- do.call(rbind, segment_rows)
  n_segments <- nrow(segments)
  adjacency <- vector("list", n_segments)
  for (i in seq_len(n_segments)) adjacency[[i]] <- integer()
  if (n_segments > 1L) {
    for (i in seq_len(n_segments - 1L)) {
      for (j in (i + 1L):n_segments) {
        if (segments$route_index[i] == segments$route_index[j]) {
          # Collinear pieces of one uninterrupted run must share a component.
          # Otherwise different competitors on its left and right can assign
          # opposite lane ranks and create an artificial H-V-H notch midway.
          interval_gap <- max(
            segments$x_min[i], segments$x_min[j]
          ) - min(segments$x_max[i], segments$x_max[j])
          contiguous <- interval_gap <= horizontal_tolerance + 1e-12 &&
            abs(segments$original_y[i] - segments$original_y[j]) <= y_tolerance
          if (contiguous) {
            adjacency[[i]] <- c(adjacency[[i]], j)
            adjacency[[j]] <- c(adjacency[[j]], i)
          }
          next
        }
        overlap <- min(segments$x_max[i], segments$x_max[j]) -
          max(segments$x_min[i], segments$x_min[j])
        if (overlap + 1e-12 < min_overlap ||
            abs(segments$original_y[i] - segments$original_y[j]) > y_tolerance) next
        adjacency[[i]] <- c(adjacency[[i]], j)
        adjacency[[j]] <- c(adjacency[[j]], i)
      }
    }
  }

  component <- integer(n_segments)
  next_component <- 0L
  for (seed in seq_len(n_segments)) {
    if (component[seed] != 0L) next
    next_component <- next_component + 1L
    queue <- seed
    component[seed] <- next_component
    while (length(queue)) {
      current <- queue[1L]
      queue <- queue[-1L]
      unseen <- adjacency[[current]][component[adjacency[[current]]] == 0L]
      if (length(unseen)) {
        component[unseen] <- next_component
        queue <- c(queue, unseen)
      }
    }
  }
  segments$component <- component
  segments$lane_y <- segments$original_y
  segments$offset <- 0

  lane_step <- max(0, as.numeric(lane_step)[1L])
  max_offset <- max(0, as.numeric(max_offset)[1L])
  for (component_id in unique(component)) {
    idx <- which(component == component_id)
    route_ids <- unique(segments$route_id[idx])
    if (length(route_ids) < 2L) next
    route_rows <- lapply(route_ids, function(route_id) {
      ridx <- idx[segments$route_id[idx] == route_id]
      data.frame(
        route_id = route_id,
        original_y = stats::median(segments$original_y[ridx]),
        priority = max(segments$priority[ridx]),
        group = sort(unique(segments$group[ridx]))[1L],
        stringsAsFactors = FALSE
      )
    })
    route_rows <- do.call(rbind, route_rows)
    route_rows <- route_rows[order(
      route_rows$original_y, -route_rows$priority,
      route_rows$group, route_rows$route_id
    ), , drop = FALSE]
    center_y <- stats::median(route_rows$original_y)
    half_span <- lane_step * (nrow(route_rows) - 1L) / 2
    original_span <- max(abs(route_rows$original_y - center_y))
    allowed_span <- max(0, max_offset - original_span)
    if (half_span > allowed_span && half_span > 0) {
      effective_step <- lane_step * allowed_span / half_span
    } else {
      effective_step <- lane_step
    }
    lane_y <- center_y +
      (seq_len(nrow(route_rows)) - (nrow(route_rows) + 1) / 2) * effective_step
    lane_map <- stats::setNames(lane_y, route_rows$route_id)
    for (route_id in route_rows$route_id) {
      ridx <- idx[segments$route_id[idx] == route_id]
      target_y <- unname(lane_map[route_id])
      shift <- target_y - segments$original_y[ridx]
      shift <- pmax(-max_offset, pmin(max_offset, shift))
      segments$offset[ridx] <- shift
      segments$lane_y[ridx] <- segments$original_y[ridx] + shift
    }
  }
  segments
}

.dnmb_cct_separate_route_lanes <- function(routes, metadata = NULL,
                                           lane_step = 0.05,
                                           max_offset = 0.18,
                                           horizontal_tolerance = 0.01,
                                           y_tolerance = 0.04,
                                           min_overlap = 0.18,
                                           transition_span = 0.04,
                                           connection_mode = c("move_internal", "transition"),
                                           round_radius = 0) {
  connection_mode <- match.arg(connection_mode)
  assignments <- .dnmb_cct_route_lane_assignments(
    routes = routes, metadata = metadata,
    lane_step = lane_step, max_offset = max_offset,
    horizontal_tolerance = horizontal_tolerance,
    y_tolerance = y_tolerance, min_overlap = min_overlap
  )
  if (!length(routes) || !nrow(assignments)) {
    if (length(routes) && round_radius > 0) {
      routes <- lapply(routes, function(points) {
        .dnmb_cct_round_orthogonal_route(points, radius = round_radius)
      })
    }
    attr(routes, "lane_assignments") <- assignments
    return(routes)
  }

  out <- vector("list", length(routes))
  for (route_index in seq_along(routes)) {
    pts <- routes[[route_index]]
    route_assignments <- assignments[assignments$route_index == route_index, , drop = FALSE]
    if (!nrow(route_assignments) || is.null(pts) || nrow(pts) < 2L) {
      out[[route_index]] <- if (round_radius > 0) {
        .dnmb_cct_round_orthogonal_route(pts, radius = round_radius)
      } else pts
      next
    }
    offset_map <- stats::setNames(route_assignments$offset, route_assignments$segment_index)
    lane_y_map <- stats::setNames(route_assignments$lane_y, route_assignments$segment_index)

    if (connection_mode == "move_internal") {
      adjusted <- as.data.frame(pts[, c("x", "y"), drop = FALSE])
      proposals <- vector("list", nrow(adjusted))
      for (segment_index in as.integer(names(lane_y_map))) {
        if (abs(offset_map[as.character(segment_index)]) <= 1e-12) next
        target_y <- unname(lane_y_map[as.character(segment_index)])
        if (segment_index > 1L) {
          proposals[[segment_index]] <- c(proposals[[segment_index]], target_y)
        }
        if (segment_index + 1L < nrow(adjusted)) {
          proposals[[segment_index + 1L]] <- c(proposals[[segment_index + 1L]], target_y)
        }
      }
      for (vertex in seq_along(proposals)) {
        if (length(proposals[[vertex]])) adjusted$y[vertex] <- stats::median(proposals[[vertex]])
      }
      adjusted$x[1L] <- pts$x[1L]
      adjusted$y[1L] <- pts$y[1L]
      adjusted$x[nrow(adjusted)] <- pts$x[nrow(pts)]
      adjusted$y[nrow(adjusted)] <- pts$y[nrow(pts)]

      # If the separated horizontal corridor starts or ends at a protected
      # route endpoint, add a vertical stem at that endpoint. This keeps the
      # endpoint exact without inserting a height-changing notch inside the
      # horizontal run. Internal elbow vertices are moved directly above.
      expanded <- adjusted[1L, c("x", "y"), drop = FALSE]
      for (segment_index in seq_len(nrow(adjusted) - 1L)) {
        key <- as.character(segment_index)
        has_lane <- key %in% names(lane_y_map) &&
          abs(offset_map[key]) > 1e-12 &&
          abs(pts$y[segment_index + 1L] - pts$y[segment_index]) <= horizontal_tolerance
        if (has_lane && segment_index == 1L) {
          expanded <- rbind(expanded, data.frame(
            x = pts$x[1L], y = unname(lane_y_map[key])
          ))
        }
        if (has_lane && segment_index + 1L == nrow(adjusted)) {
          expanded <- rbind(expanded, data.frame(
            x = pts$x[nrow(pts)], y = unname(lane_y_map[key])
          ))
        }
        expanded <- rbind(expanded, adjusted[segment_index + 1L, c("x", "y"), drop = FALSE])
      }
      adjusted <- expanded
      adjusted <- adjusted[c(TRUE, diff(adjusted$x) != 0 | diff(adjusted$y) != 0), , drop = FALSE]
      out[[route_index]] <- if (round_radius > 0) {
        .dnmb_cct_round_orthogonal_route(adjusted, radius = round_radius)
      } else adjusted
      next
    }

    adjusted <- pts[1L, c("x", "y"), drop = FALSE]
    for (segment_index in seq_len(nrow(pts) - 1L)) {
      p0 <- pts[segment_index, c("x", "y"), drop = FALSE]
      p1 <- pts[segment_index + 1L, c("x", "y"), drop = FALSE]
      key <- as.character(segment_index)
      offset <- if (key %in% names(offset_map)) unname(offset_map[key]) else 0
      if (!is.finite(offset) || abs(offset) <= 1e-12) {
        adjusted <- rbind(adjusted, p1)
        next
      }
      dx <- p1$x - p0$x
      inset <- min(max(0, transition_span), abs(dx) * 0.25)
      direction <- sign(dx)
      before <- data.frame(x = p0$x + direction * inset, y = p0$y)
      lane_before <- data.frame(x = before$x, y = p0$y + offset)
      lane_after <- data.frame(x = p1$x - direction * inset, y = p1$y + offset)
      after <- data.frame(x = lane_after$x, y = p1$y)
      adjusted <- rbind(adjusted, before, lane_before, lane_after, after, p1)
    }
    adjusted <- adjusted[c(TRUE, diff(adjusted$x) != 0 | diff(adjusted$y) != 0), , drop = FALSE]
    out[[route_index]] <- adjusted
  }
  attr(out, "lane_assignments") <- assignments
  out
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
                                           branches = NULL,
                                           sym_size = .dnmb_cct_sugar_icon_radius(),
                                           gap = 0.02,
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

.dnmb_cct_extra_chain_step <- function(
    sym_size = .dnmb_cct_sugar_icon_radius(), gap = 0.18) {
  sym_size * 2 + gap
}

.dnmb_cct_extra_chain_node_x <- function(
    chain_x, pos, sym_size = .dnmb_cct_sugar_icon_radius(), gap = 0.18) {
  chain_x + (as.numeric(pos)[1] - 1) * .dnmb_cct_extra_chain_step(sym_size = sym_size, gap = gap)
}

.dnmb_cct_extra_route_trim <- function(
    sym_size = .dnmb_cct_sugar_icon_radius()) {
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
                                   sym_size = .dnmb_cct_sugar_icon_radius(),
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
                                   size = .dnmb_cct_sugar_icon_radius(),
                                   gap = 0.18, show_label = TRUE,
                                   show_bond_label = TRUE) {

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
        paste0(anom_sym, pos_parts[1L], "-", pos_parts[2L])
      } else {
        paste0(anom_sym, pos_parts[1L])
      }
      if (isTRUE(show_bond_label)) {
        layers[[length(layers) + 1L]] <- ggplot2::geom_text(
          data    = data.frame(x = mid_x, y = y + size * 2.0,
                               label = bond_label),
          mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                 label = .data[["label"]]),
          size = 1.2, color = "#888888", hjust = 0.5, inherit.aes = FALSE
        )
      }

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

#' Compute the geometry for the GH-cleavage scissors glyph
#' @keywords internal
.dnmb_scissors_geometry_v2 <- function(x, y, size = 0.12, angle = 0) {
  d <- max(0.001, as.numeric(size)[1])
  rad <- as.numeric(angle)[1] * pi / 180
  rotate <- function(dx, dy) {
    data.frame(
      x = x + dx * cos(rad) - dy * sin(rad),
      y = y + dx * sin(rad) + dy * cos(rad)
    )
  }

  handles <- rotate(c(-0.78, -0.78) * d, c(0.42, -0.42) * d)
  handles$group <- c("upper", "lower")
  pivot <- rotate(0, 0)
  tips <- rotate(c(1.12, 1.12) * d, c(-0.42, 0.42) * d)
  tips$group <- c("lower", "upper")

  shanks <- data.frame(
    x = handles$x, y = handles$y,
    xend = pivot$x, yend = pivot$y,
    group = handles$group
  )
  blades <- data.frame(
    x = pivot$x, y = pivot$y,
    xend = tips$x, yend = tips$y,
    group = tips$group
  )
  list(
    handles = handles,
    pivot = pivot,
    tips = tips,
    shanks = shanks,
    blades = blades,
    handle_radius = 0.27 * d,
    pivot_radius = 0.12 * d
  )
}

#' Draw a foreground scissors symbol for GH cleavage
#' @param x,y Position of the pivot
#' @param size Scale factor in plot coordinates
#' @param angle Rotation angle in degrees
#' @return List of ggplot2 layers
#' @keywords internal
.dnmb_scissors_grob_v2 <- function(x, y, size = 0.12, angle = 0) {
  geom <- .dnmb_scissors_geometry_v2(x, y, size = size, angle = angle)
  all_segments <- rbind(geom$shanks, geom$blades)

  layers <- list(
    ggplot2::geom_segment(
      data = all_segments,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      linewidth = 1.25, color = "#FFFFFF", lineend = "round", inherit.aes = FALSE),
    ggplot2::geom_segment(
      data = geom$blades,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      linewidth = 0.62, color = "#687078", lineend = "round", inherit.aes = FALSE),
    ggplot2::geom_segment(
      data = geom$shanks,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      linewidth = 0.62, color = "#C62828", lineend = "round", inherit.aes = FALSE)
  )
  for (i in seq_len(nrow(geom$handles))) {
    layers[[length(layers) + 1L]] <- .dnmb_native_circle_layer_v2(
      geom$handles$x[i], geom$handles$y[i], radius = geom$handle_radius,
      fill = "#FFFFFF", color = "#C62828", linewidth = 0.42
    )
  }
  layers[[length(layers) + 1L]] <- .dnmb_native_circle_layer_v2(
    geom$pivot$x, geom$pivot$y, radius = geom$pivot_radius,
    fill = "#FFFFFF", color = "#7F1D1D", linewidth = 0.38
  )
  layers
}

# ---------------------------------------------------------------------------
# 3-Zone CAZy Carbon Transport Map — data extraction helpers
# ---------------------------------------------------------------------------

#' Extract GH enzymes with locus_tag, gene_name, family from genbank_table
#' @keywords internal
.dnmb_cct_3zone_extract_gh <- function(genbank_table) {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE, row.names = NULL)
  family_cols <- intersect(
    c("dbCAN_dbcan_all_families", "dbCAN_dbcan_hit", "dbCAN_family_id", "family_id"),
    names(tbl)
  )
  if (!length(family_cols) || !nrow(tbl)) return(NULL)

  gh_by_row <- lapply(seq_len(nrow(tbl)), function(i) {
    families <- character()
    for (column in family_cols) {
      families <- .dnmb_dbcan_family_tokens(tbl[[column]][i])
      if (length(families)) break
    }
    families[grepl("^GH", families)]
  })
  keep <- lengths(gh_by_row) > 0L
  if (!any(keep)) return(NULL)

  row_index <- rep(which(keep), lengths(gh_by_row[keep]))
  gh <- tbl[row_index, , drop = FALSE]
  gh$family <- unlist(gh_by_row[keep], use.names = FALSE)
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
      is_t <- grepl("transport|PTS|permease|ABC|SBP|MFS|EIICB|IIBC|EIIBC|EIIA|crr|porter|uptake|\\bdeoP\\b", tbl2[[step_col]], ignore.case = TRUE) |
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
          is_t <- grepl("transport|PTS|pts|permease|ABC|SBP|MFS|\\bdeoP\\b", st$step, ignore.case = TRUE)
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
              "Deoxyribose","Fucose","Rhamnose","Gluconate")
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
      fam_cs_base["Other"]+dx_cs, fam_cs_base["Other"]+2*dx_cs,
      fam_cs_base["Other"]+3*dx_cs, fam_cs_base["Other"]+4*dx_cs)
    cs_xs <- stats::setNames(cs_xs_default, all_cs)
  }
  cs_map <- cs_xs

  # Backbone top y -- just below membrane
  by <- 7.2  # shift the upper cytoplasmic grid below the membrane gutter
  cs_y <- by + 0.45

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
    Other = compact_anchor(c("Glycerol", "Deoxyribose", "Fucose", "Rhamnose", "Gluconate"), bx + 2.0, fallback = fam_cs_base["Other"], pull_left = 0.35, pull_right = 0.45)
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
    n("Isocit",    tx[3], ty[3], "Isocitrate", "tca", "intermediate"),
    n("AKG",       tx[4], ty[4], "a-KG",       "tca", "intermediate"),
    n("SucCoA",    tx[5], ty[5], "Suc-CoA",    "tca", "intermediate"),
    n("Succinate", tx[6], ty[6], "Succinate",  "tca", "intermediate"),
    n("Fumarate",  tx[7], ty[7], "Fumarate",   "tca", "intermediate"),
    n("Malate",    tx[8], ty[8], "Malate",     "tca", "intermediate"),

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
    n("Deoxyribose-5-P", csx("Deoxyribose", fam_cs["Other"]), 7.0,
      "dRibose-5-P", "entry_intermediate", "deoxyribose"),
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
    n("Glyoxylate", tca_cx + 0.34, tca_cy + 0.22, "Glyoxylate", "tca_shunt", "intermediate"),

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
    e("Deoxyribose-5-P", "GA3P", "deoxyribose", -0.12), # deoC product 1
    e("Deoxyribose-5-P", "Acetaldehyde", "deoxyribose", 0.12), # deoC product 2
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
    e("Acetaldehyde", "Acetate", "deoxyribose"),        # deoxyribose adh
    e("Acetate", "Acetyl-P", "deoxyribose"),            # deoxyribose ackA
    e("Acetyl-P", "Acetyl-CoA", "deoxyribose"),         # deoxyribose pta
    e("Acetate", "Acetyl-CoA", "deoxyribose"),          # deoxyribose acs
    e("Acetaldehyde", "Acetyl-CoA", "deoxyribose"),     # deoxyribose ald-dh-CoA

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
    e("Deoxyribose", "Deoxyribose-5-P", "deoxyribose"),
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
  is_carbon_source <- nodes$type == "carbon_source"
  tca_save_x <- nodes$x[is_tca]
  tca_save_y <- nodes$y[is_tca]
  carbon_source_x <- nodes$x[is_carbon_source]
  nodes$x <- .dnmb_cct_snap_to_grid(nodes$x, step = step)
  nodes$y <- .dnmb_cct_snap_to_grid(nodes$y, step = step)
  # Restore TCA exact positions
  nodes$x[is_tca] <- tca_save_x
  nodes$y[is_tca] <- tca_save_y
  # Carbon-source lanes are optimized on a 0.25-unit lattice.  Preserve that
  # spacing so a source is not snapped onto a neighboring pathway's chain.
  nodes$x[is_carbon_source] <- carbon_source_x

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
    Deoxyribose = "deoxy",
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
    deoxyribose = "#8C6BB1",  # muted violet
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
    chitin = "NAG", mannan = "Mannose",
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
    fructose = "Fructose", glucose = "Glucose",
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
    Sucrose     = "^(sacA|scrB|sacB|sacC|inv[ABCD])",
    Mannitol    = "^(mtlD|mtlK|mtlA)",
    Glycerol    = "^(glpK|glpD|glpA|glpB|glpC|dhaK|dhaL|dhaM)",
    Xylose      = "^(xylA|xylB|xylR)",
    Arabinose   = "^(araA|araB|araD|araR)",
    Ribose      = "^(rbsK|rbsR|rbsD)",
    Deoxyribose = "^(deoK|deoC)",
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
    Deoxyribose = "deoxyribose kinase|deoxyribose-5-phosphate|deoxyribose-phosphate aldolase",
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
    Deoxyribose = "^2\\.7\\.1\\.15$|^4\\.1\\.2\\.4$",                    # deoK, deoC
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

# Conservative annotation scan for the central-carbon reactions that are not
# assessed by GapMindCarbon's substrate-utilization paths.  A strong hit needs
# agreement between at least two independent annotation axes (gene name,
# product, EC, or KO); single exact gene/EC/KO matches remain partial evidence.
.dnmb_cct_core_metabolism_evidence <- function(genbank_table) {
  if (is.list(genbank_table) && !is.data.frame(genbank_table) &&
      "features" %in% names(genbank_table)) {
    genbank_table <- genbank_table$features
  }
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE, row.names = NULL)

  empty_hits <- data.frame(
    reaction_id = character(), pathway = character(), component = character(),
    component_required = logical(), locus_tag = character(), gene_name = character(),
    product = character(), evidence_level = character(), evidence_axes = character(),
    stringsAsFactors = FALSE
  )

  component <- function(name, required = TRUE, genes = character(),
                        products = character(), ecs = character(), kos = character(),
                        exclude_genes = character(), exclude_products = character(),
                        exclude_kos = character(), allow_cs_ambiguity = FALSE) {
    list(
      name = name, required = required, genes = genes, products = products,
      ecs = ecs, kos = kos, exclude_genes = exclude_genes,
      exclude_products = exclude_products, exclude_kos = exclude_kos,
      allow_cs_ambiguity = allow_cs_ambiguity
    )
  }

  catalog <- list(
    list(id = "C01", pathway = "TCA cycle", from = "OAA", to = "Citrate",
         label = "Citrate synthase", mode = "all", components = list(
           component(
             "CS", genes = c("gltA", "citZ"),
             products = c("citrate synthase"), ecs = "2.3.3.1", kos = "K01647",
             exclude_genes = c("prpC", "mmgD", "cimA", "leuA"),
             exclude_products = c("methylcitrate", "2-methylcitrate", "citramalate",
                                  "homocitrate", "isopropylmalate"),
             exclude_kos = "K01659", allow_cs_ambiguity = TRUE
           )
         )),
    list(id = "C02", pathway = "TCA cycle", from = "Citrate", to = "Isocit",
         label = "Aconitate hydratase", mode = "all", components = list(
           component(
             "ACN", genes = c("acn", "acnA", "acnB", "citB"),
             products = c("aconitate hydratase", "aconitase"),
             ecs = "4.2.1.3", kos = c("K01681", "K01682"),
             exclude_genes = c("prpD", "leuC", "leuD"),
             exclude_products = c("methylcitrate", "isopropylmalate", "homoaconitate")
           )
         )),
    list(id = "C03", pathway = "TCA cycle", from = "Isocit", to = "AKG",
         label = "Isocitrate dehydrogenase", mode = "all", components = list(
           component(
             "ICD", genes = c("icd", "icdA", "icdB", "idh", "idhA", "idp"),
             products = c("isocitrate dehydrogenase"),
             ecs = c("1.1.1.41", "1.1.1.42"), kos = c("K00030", "K00031"),
             exclude_genes = c("leuB"),
             exclude_products = c("isopropylmalate", "isocitrate/isopropylmalate")
           )
         )),
    list(id = "C04", pathway = "TCA cycle", from = "AKG", to = "SucCoA",
         label = "2-Oxoglutarate DH E1/E2", mode = "all", components = list(
           component(
             "E1", genes = c("sucA", "odhA", "ogdh"),
             products = c("2-oxoglutarate dehydrogenase.*E1",
                          "alpha-ketoglutarate dehydrogenase.*E1"),
             ecs = "1.2.4.2", kos = "K00164",
             exclude_genes = c("pdhA", "bkdA"),
             exclude_products = c("pyruvate dehydrogenase", "branched-chain.*dehydrogenase")
           ),
           component(
             "E2", genes = c("sucB", "odhB"),
             products = c("dihydrolipo.*succinyltransferase",
                          "2-oxoglutarate dehydrogenase.*E2"),
             ecs = "2.3.1.61", kos = "K00658",
             exclude_genes = c("pdhC", "bkdB"),
             exclude_products = c("pyruvate dehydrogenase", "branched-chain.*dehydrogenase")
           )
         )),
    list(id = "C05", pathway = "TCA cycle", from = "SucCoA", to = "Succinate",
         label = "Succinate-CoA ligase alpha/beta", mode = "all", components = list(
           component(
             "alpha", genes = "sucD", products = c("succinate--CoA ligase.*alpha",
                                                    "succinyl-CoA synthetase.*alpha"),
             ecs = c("6.2.1.4", "6.2.1.5"), kos = "K01902"
           ),
           component(
             "beta", genes = "sucC", products = c("succinate--CoA ligase.*beta",
                                                   "succinyl-CoA synthetase.*beta"),
             ecs = c("6.2.1.4", "6.2.1.5"), kos = "K01903"
           )
         )),
    list(id = "C06", pathway = "TCA cycle", from = "Succinate", to = "Fumarate",
         label = "Succinate dehydrogenase A/B", mode = "all", components = list(
           component(
             "A", genes = "sdhA", products = "succinate dehydrogenase.*flavoprotein",
             ecs = "1.3.5.1", kos = "K00239",
             exclude_genes = "frdA", exclude_products = "fumarate reductase"
           ),
           component(
             "B", genes = "sdhB", products = "succinate dehydrogenase.*iron-sulfur",
             ecs = "1.3.5.1", kos = "K00240",
             exclude_genes = "frdB", exclude_products = "fumarate reductase"
           ),
           component(
             "membrane", genes = c("sdhC", "sdhD"),
             products = c("succinate dehydrogenase.*cytochrome",
                          "succinate dehydrogenase.*membrane"),
             ecs = "1.3.5.1", kos = c("K00241", "K00242"),
             exclude_genes = c("frdC", "frdD"), exclude_products = "fumarate reductase"
           ),
           component(
             "FRD alternative", required = FALSE,
             genes = c("frdA", "frdB", "frdC", "frdD"),
             products = "fumarate reductase", ecs = "1.3.5.4",
             kos = c("K00244", "K00245", "K00246", "K00247")
           )
         )),
    list(id = "C07", pathway = "TCA cycle", from = "Fumarate", to = "Malate",
         label = "Fumarate hydratase", mode = "all", components = list(
           component(
             "FUM", genes = c("fum", "fumA", "fumB", "fumC"),
             products = c("fumarate hydratase", "fumarase"),
             ecs = "4.2.1.2", kos = c("K01676", "K01677", "K01678"),
             exclude_products = c("fumarylacetoacetate hydrolase", "fumarate reductase",
                                  "tartrate hydratase")
           )
         )),
    list(id = "C08", pathway = "TCA cycle", from = "Malate", to = "OAA",
         label = "Malate DH / Mqo", mode = "any", components = list(
           component(
             "MDH", genes = "mdh", products = "malate dehydrogenase",
             ecs = "1.1.1.37", kos = "K00024",
             exclude_products = c("lactate dehydrogenase", "malic enzyme",
                                  "isopropylmalate", "methylmalonate-semialdehyde")
           ),
           component(
             "Mqo", genes = c("mqo", "mqoA"),
             products = c("malate dehydrogenase \\(quinone\\)",
                          "malate:quinone oxidoreductase", "malate quinone oxidoreductase"),
             ecs = "1.1.5.4", kos = "K00116",
             exclude_products = c("lactate dehydrogenase", "malic enzyme")
           )
         )),
    list(id = "C09", pathway = "Glyoxylate shunt", from = "Isocit", to = "Glyoxylate",
         label = "Isocitrate lyase", mode = "all", components = list(
           component(
             "ICL", genes = "aceA", products = "isocitrate lyase",
             ecs = "4.1.3.1", kos = "K01637",
             exclude_genes = "prpB", exclude_products = "methylisocitrate lyase",
             exclude_kos = "K03417"
           )
         )),
    list(id = "C10", pathway = "Glyoxylate shunt", from = "Glyoxylate", to = "Malate",
         label = "Malate synthase", mode = "all", components = list(
           component(
             "MS", genes = c("aceB", "glcB"), products = c("malate synthase A?", "malate synthase"),
             ecs = "2.3.3.9", kos = "K01638",
             exclude_genes = "bshA",
             exclude_products = c("N-acetylglucosaminyl.*malate synthase", "glucosaminyl.*malate synthase")
           )
         )),
    list(id = "C11", pathway = "Glycolysis", from = "Glucose", to = "Glc-6-P",
         label = "Glucose phosphorylation", mode = "any", components = list(
           component(
             "Glk", genes = c("glk", "glcK"), products = c("glucokinase", "glucose kinase"),
             ecs = "2.7.1.2", kos = "K00845",
             exclude_products = c("regulatory protein", "transcriptional regulator")
           ),
           component(
             "Hxk", genes = c("hxk", "hexK"), products = "hexokinase",
             ecs = "2.7.1.1", kos = "K00844",
             exclude_products = "regulatory protein"
           )
         )),
    list(id = "C12", pathway = "Glycolysis", from = "Glc-6-P", to = "Fru-6-P",
         label = "Glucose-6-P isomerase", mode = "all", components = list(
           component(
             "Pgi", genes = c("pgi", "pgiA"),
             products = "glucose-6-phosphate isomerase",
             ecs = "5.3.1.9", kos = "K01810",
             exclude_genes = "manA", exclude_products = "mannose-6-phosphate isomerase"
           )
         )),
    list(id = "C13", pathway = "Glycolysis", from = "Fru-6-P", to = "Fru-1,6-BP",
         label = "6-Phosphofructokinase", mode = "all", components = list(
           component(
             "Pfk", genes = c("pfk", "pfkA"), products = "6-phosphofructokinase",
             ecs = "2.7.1.11", kos = c("K00850", "K21071"),
             exclude_genes = c("fruK", "pfkB"),
             exclude_products = c("1-phosphofructokinase", "ribokinase family")
           )
         )),
    list(id = "C14", pathway = "Glycolysis", from = "Fru-1,6-BP", to = "GA3P",
         label = "Fructose-bisphosphate aldolase", mode = "all", components = list(
           component(
             "Fba", genes = c("fba", "fbaA", "fbaB"),
             products = "fructose.*bisphosphate aldolase",
             ecs = "4.1.2.13", kos = c("K01623", "K01624", "K11645"),
             exclude_products = c("2-dehydro-3-deoxy", "tagatose-bisphosphate aldolase")
           )
         )),
    list(id = "C15", pathway = "Glycolysis", from = "DHAP", to = "GA3P",
         label = "Triose-phosphate isomerase", mode = "all", components = list(
           component(
             "Tpi", genes = c("tpi", "tpiA"), products = "triose-phosphate isomerase",
             ecs = "5.3.1.1", kos = "K01803"
           )
         )),
    list(id = "C16", pathway = "Glycolysis", from = "GA3P", to = "1,3-BPG",
         label = "Glyceraldehyde-3-P DH", mode = "all", components = list(
           component(
             "Gap", genes = c("gap", "gapA", "gapB"),
             products = "glyceraldehyde-3-phosphate dehydrogenase",
             ecs = "1.2.1.12", kos = "K00134",
             exclude_genes = "gapN", exclude_products = "non-phosphorylating"
           )
         )),
    list(id = "C17", pathway = "Glycolysis", from = "1,3-BPG", to = "3-PG",
         label = "Phosphoglycerate kinase", mode = "all", components = list(
           component(
             "Pgk", genes = "pgk", products = "phosphoglycerate kinase",
             ecs = "2.7.2.3", kos = "K00927"
           )
         )),
    list(id = "C18", pathway = "Glycolysis", from = "3-PG", to = "2-PG",
         label = "Phosphoglycerate mutase", mode = "all", components = list(
           component(
             "Gpm", genes = c("gpm", "gpmA", "gpmB", "gpmI"),
             products = c("phosphoglycerate mutase", "2,3-bisphosphoglycerate-independent phosphoglycerate mutase"),
             ecs = c("5.4.2.11", "5.4.2.12"),
             kos = c("K01834", "K15633", "K15634"),
             exclude_products = c("histidine phosphatase family", "cofactor-dependent phosphoglycerate mutase regulator")
           )
         )),
    list(id = "C19", pathway = "Glycolysis", from = "2-PG", to = "PEP",
         label = "Enolase", mode = "all", components = list(
           component(
             "Eno", genes = "eno", products = c("enolase", "phosphopyruvate hydratase"),
             ecs = "4.2.1.11", kos = "K01689",
             exclude_products = c("methylthiopentyl", "mandelate racemase")
           )
         )),
    list(id = "C20", pathway = "Glycolysis", from = "PEP", to = "Pyruvate",
         label = "Pyruvate kinase", mode = "all", components = list(
           component(
             "Pyk", genes = c("pyk", "pykA", "pykF"), products = "pyruvate kinase",
             ecs = "2.7.1.40", kos = "K00873",
             exclude_products = c("pyruvate phosphate dikinase", "phosphoenolpyruvate synthase")
           )
         ))
  )

  pick <- function(candidates, default = "") {
    hit <- intersect(candidates, names(tbl))
    if (!length(hit)) return(rep(default, nrow(tbl)))
    out <- as.character(tbl[[hit[1]]])
    out[is.na(out)] <- default
    trimws(out)
  }
  locus <- pick(c("locus_tag", "old_locus_tag"))
  gene <- pick(c("gene", "gene_name"))
  preferred <- pick(c("EggNOG_Preferred_name", "eggnog_preferred_name"))
  product <- pick(c("product", "description"))
  egg_desc <- pick(c("EggNOG_Description", "eggnog_description"))
  ec_text <- paste(pick(c("EC_number", "ec_number")), pick(c("EggNOG_EC", "eggnog_ec")))
  ko_text <- paste(pick(c("EggNOG_KEGG_ko", "eggnog_kegg_ko")))

  has_token <- function(values, tokens) {
    if (!length(tokens)) return(rep(FALSE, length(values)))
    tokens <- toupper(tokens)
    vapply(toupper(values), function(value) {
      observed <- unlist(strsplit(value, "[^A-Z0-9.]+"))
      any(observed[nzchar(observed)] %in% tokens)
    }, logical(1))
  }
  has_gene <- function(values, aliases) {
    if (!length(aliases)) return(rep(FALSE, length(values)))
    aliases <- tolower(aliases)
    vapply(tolower(values), function(value) {
      observed <- unlist(strsplit(value, "[^a-z0-9]+"))
      any(observed[nzchar(observed)] %in% aliases)
    }, logical(1))
  }
  has_pattern <- function(values, patterns) {
    if (!length(patterns)) return(rep(FALSE, length(values)))
    grepl(paste(patterns, collapse = "|"), values, ignore.case = TRUE, perl = TRUE)
  }

  hit_rows <- list()
  step_rows <- list()
  for (spec in catalog) {
    component_levels <- character()
    for (rule in spec$components) {
      trusted_gene <- has_gene(gene, rule$genes)
      preferred_gene <- has_gene(preferred, rule$genes)
      gene_axis <- trusted_gene | preferred_gene
      product_axis <- has_pattern(product, rule$products)
      ec_axis <- has_token(ec_text, rule$ecs)
      ko_axis <- has_token(ko_text, rule$kos)

      negative_gene <- has_gene(gene, rule$exclude_genes) |
        has_gene(preferred, rule$exclude_genes)
      negative_product <- has_pattern(product, rule$exclude_products)
      negative_ko <- has_token(ko_text, rule$exclude_kos)
      cs_exception <- isTRUE(rule$allow_cs_ambiguity) & gene_axis & ko_axis &
        !negative_gene & !negative_ko
      rejected <- negative_gene | negative_ko | (negative_product & !cs_exception)

      strong <- !rejected & (
        (gene_axis & (product_axis | ec_axis | ko_axis)) |
          (product_axis & (ec_axis | ko_axis))
      )
      candidate <- !rejected & !strong & (trusted_gene | ec_axis | ko_axis)
      keep <- (strong | candidate) & nzchar(locus)
      level <- ifelse(strong, "strong", ifelse(candidate, "candidate", "none"))
      component_levels[rule$name] <- if (any(strong & keep)) {
        "strong"
      } else if (any(candidate & keep)) {
        "candidate"
      } else {
        "none"
      }

      if (any(keep)) {
        axes <- vapply(which(keep), function(i) {
          paste(c(
            if (gene_axis[i]) "gene" else NULL,
            if (product_axis[i]) "product" else NULL,
            if (ec_axis[i]) "EC" else NULL,
            if (ko_axis[i]) "KO" else NULL
          ), collapse = "+")
        }, character(1))
        display_gene <- gene[keep]
        missing_gene <- !nzchar(display_gene)
        display_gene[missing_gene] <- preferred[keep][missing_gene]
        display_gene[!nzchar(display_gene)] <- "unassigned"
        display_product <- product[keep]
        display_product[!nzchar(display_product)] <- egg_desc[keep][!nzchar(display_product)]
        display_product[!nzchar(display_product)] <- "not annotated"
        hit_rows[[length(hit_rows) + 1L]] <- data.frame(
          reaction_id = spec$id,
          pathway = spec$pathway,
          component = rule$name,
          component_required = isTRUE(rule$required),
          locus_tag = locus[keep],
          gene_name = display_gene,
          product = display_product,
          evidence_level = level[keep],
          evidence_axes = axes,
          stringsAsFactors = FALSE
        )
      }
    }

    required_names <- vapply(
      Filter(function(rule) isTRUE(rule$required), spec$components),
      function(rule) rule$name, character(1)
    )
    required_levels <- component_levels[required_names]
    strong_complete <- if (identical(spec$mode, "any")) {
      any(component_levels == "strong")
    } else {
      length(required_levels) > 0L && all(required_levels == "strong")
    }
    any_evidence <- any(component_levels != "none")
    status <- if (strong_complete) "active" else if (any_evidence) "partial" else "reference"
    step_rows[[length(step_rows) + 1L]] <- data.frame(
      reaction_id = spec$id,
      pathway = spec$pathway,
      from = spec$from,
      to = spec$to,
      reaction_label = spec$label,
      status = status,
      stringsAsFactors = FALSE
    )
  }

  hits <- if (length(hit_rows)) dplyr::bind_rows(hit_rows) else empty_hits
  if (nrow(hits)) {
    hits <- hits[!duplicated(hits[, c("reaction_id", "component", "locus_tag")]), , drop = FALSE]
    hits <- hits[order(hits$reaction_id, hits$component, hits$locus_tag), , drop = FALSE]
  }
  steps <- dplyr::bind_rows(step_rows)
  summary_hits <- hits
  if (nrow(summary_hits)) {
    keep_summary <- unlist(lapply(
      split(seq_len(nrow(summary_hits)), paste(summary_hits$reaction_id, summary_hits$component)),
      function(idx) {
        strong_idx <- idx[summary_hits$evidence_level[idx] == "strong"]
        if (length(strong_idx)) strong_idx else idx
      }
    ), use.names = FALSE)
    summary_hits <- summary_hits[sort(unique(keep_summary)), , drop = FALSE]
    active_ids <- steps$reaction_id[steps$status == "active"]
    summary_hits <- summary_hits[
      !(summary_hits$reaction_id %in% active_ids & !summary_hits$component_required),
      , drop = FALSE
    ]
  }
  steps$locus_tags <- vapply(steps$reaction_id, function(id) {
    vals <- unique(summary_hits$locus_tag[summary_hits$reaction_id == id])
    paste(vals[nzchar(vals)], collapse = "; ")
  }, character(1))
  steps$gene_names <- vapply(steps$reaction_id, function(id) {
    vals <- unique(summary_hits$gene_name[summary_hits$reaction_id == id])
    paste(vals[nzchar(vals)], collapse = "; ")
  }, character(1))
  steps$product_evidence <- vapply(steps$reaction_id, function(id) {
    sub <- summary_hits[summary_hits$reaction_id == id, , drop = FALSE]
    if (!nrow(sub)) return("")
    paste(paste0(sub$locus_tag, ": ", sub$product), collapse = "; ")
  }, character(1))
  steps$components <- vapply(steps$reaction_id, function(id) {
    sub <- summary_hits[summary_hits$reaction_id == id, , drop = FALSE]
    if (!nrow(sub)) return("")
    paste(unique(paste0(sub$component, "=", sub$locus_tag)), collapse = "; ")
  }, character(1))
  list(steps = steps, hits = hits, display_hits = summary_hits)
}

# Build compact gene/locus labels beside conservatively supported central-
# carbon reactions.  Positions are deterministic so the labels remain tied to
# the same reaction rather than drifting into unrelated pathway lanes.
.dnmb_cct_core_direct_labels <- function(core_steps, display_hits, nodes) {
  empty <- data.frame(
    reaction_id = character(), pathway = character(), x = numeric(), y = numeric(),
    label = character(), hjust = numeric(), vjust = numeric(), size = numeric(),
    status = character(), stringsAsFactors = FALSE
  )
  if (is.null(core_steps) || !is.data.frame(core_steps) || !nrow(core_steps) ||
      is.null(display_hits) || !is.data.frame(display_hits) || !nrow(display_hits) ||
      is.null(nodes) || !is.data.frame(nodes) || !nrow(nodes)) {
    return(empty)
  }

  node_x <- stats::setNames(nodes$x, nodes$id)
  node_y <- stats::setNames(nodes$y, nodes$id)
  detected <- core_steps[
    core_steps$status %in% c("active", "partial") &
      core_steps$from %in% names(node_x) & core_steps$to %in% names(node_x),
    , drop = FALSE
  ]
  if (!nrow(detected)) return(empty)

  label_for <- function(reaction_id) {
    sub <- display_hits[display_hits$reaction_id == reaction_id, , drop = FALSE]
    if (!nrow(sub)) return("")
    sub <- sub[order(sub$component, sub$locus_tag), , drop = FALSE]
    sub <- sub[!duplicated(sub$locus_tag), , drop = FALSE]
    genes <- trimws(as.character(sub$gene_name))
    missing_gene <- is.na(genes) | !nzchar(genes) | genes == "unassigned"
    genes[missing_gene] <- trimws(as.character(sub$component[missing_gene]))
    paste(paste(genes, sub$locus_tag, sep = "  "), collapse = "\n")
  }

  rows <- list()
  tca_ids <- c("OAA", "Citrate", "Isocit", "AKG", "SucCoA",
               "Succinate", "Fumarate", "Malate")
  tca_nodes <- nodes[nodes$id %in% tca_ids, , drop = FALSE]
  tca_cx <- if (nrow(tca_nodes)) mean(tca_nodes$x) else NA_real_
  tca_cy <- if (nrow(tca_nodes)) mean(tca_nodes$y) else NA_real_
  tca_r <- if (nrow(tca_nodes)) {
    max(sqrt((tca_nodes$x - tca_cx)^2 + (tca_nodes$y - tca_cy)^2))
  } else {
    NA_real_
  }

  for (i in seq_len(nrow(detected))) {
    step <- detected[i, , drop = FALSE]
    label <- label_for(step$reaction_id)
    if (!nzchar(label)) next
    x_from <- unname(node_x[step$from]); y_from <- unname(node_y[step$from])
    x_to <- unname(node_x[step$to]); y_to <- unname(node_y[step$to])
    x <- mean(c(x_from, x_to)); y <- mean(c(y_from, y_to))
    hjust <- 0.5; vjust <- 0.5

    if (identical(step$pathway, "TCA cycle") && is.finite(tca_r)) {
      a_from <- atan2(y_from - tca_cy, x_from - tca_cx)
      a_to <- atan2(y_to - tca_cy, x_to - tca_cx)
      delta <- ((a_to - a_from + pi) %% (2 * pi)) - pi
      mid_angle <- a_from + delta / 2
      label_r <- if (step$reaction_id %in% c("C01", "C08")) {
        max(0.65, tca_r - 0.40)
      } else {
        tca_r + 0.33
      }
      x <- tca_cx + label_r * cos(mid_angle)
      y <- tca_cy + label_r * sin(mid_angle)
      hjust <- if (step$reaction_id %in% c("C01", "C02", "C03", "C04")) 0 else 1
      if (step$reaction_id %in% c("C01", "C08")) y <- y - 0.28
      if (step$reaction_id %in% c("C04", "C05")) y <- y - 0.06
    } else if (identical(step$pathway, "Glyoxylate shunt")) {
      y <- y - if (identical(step$reaction_id, "C09")) 0.30 else 0.22
    } else if (identical(step$pathway, "Glycolysis")) {
      if (step$reaction_id %in% c("C11", "C12", "C13")) {
        x <- min(x_from, x_to) - 0.28
        hjust <- 1
      } else if (identical(step$reaction_id, "C14")) {
        x <- max(x_from, x_to) + 0.24
        hjust <- 0
      } else if (identical(step$reaction_id, "C15")) {
        y <- y - 0.20
      } else {
        x <- max(x_from, x_to) + 0.26
        hjust <- 0
      }
    }

    label_lines <- strsplit(label, "\n", fixed = TRUE)[[1L]]
    longest <- max(nchar(label_lines), na.rm = TRUE)
    size <- max(0.88, min(1.12, 1.16 - 0.055 * (length(label_lines) - 1L) -
      0.010 * max(0, longest - 24L)))
    rows[[length(rows) + 1L]] <- data.frame(
      reaction_id = step$reaction_id, pathway = step$pathway,
      x = x, y = y, label = label, hjust = hjust, vjust = vjust,
      size = size, status = step$status,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows)) dplyr::bind_rows(rows) else empty
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
  core_metabolism <- .dnmb_cct_core_metabolism_evidence(genbank_table)
  core_steps <- core_metabolism$steps
  core_display_hits <- core_metabolism$display_hits
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
  pathway_to_cs <- .dnmb_cct_transporter_pathway_map()

  # Filter step_status to ONLY carbon source pathways (exclude amino acid, organic acid, etc.)
  matched_steps <- step_status[
    !is.na(step_status$locus_tag) & nzchar(step_status$locus_tag) &
    tolower(step_status$pathway_id) %in% valid_carbon_pathways,
    , drop = FALSE
  ]
  if (nrow(matched_steps) > 0) {
    matched_steps$rank <- conf_rank[matched_steps$confidence]
    matched_steps$rank[is.na(matched_steps$rank)] <- 0L
    # score-0 GapMind assignments are candidate references, not foreground
    # evidence.  They remain available in the module tables but do not turn an
    # entire utilization route on in this map.
    matched_steps <- matched_steps[matched_steps$rank >= 2L, , drop = FALSE]
  }
  matched_pathway_ids <- unique(tolower(matched_steps$pathway_id))
  pathway_presence <- .dnmb_cct_pathway_presence(
    step_status = step_status,
    pstats = pstats,
    valid_pathways = valid_carbon_pathways
  )
  cytoplasm_status <- stats::setNames(
    pathway_presence$cytoplasm_status,
    pathway_presence$pathway_id
  )
  pathway_state_for <- function(path_ids) {
    states <- unname(cytoplasm_status[tolower(as.character(path_ids))])
    states <- states[!is.na(states)]
    if (any(states == "active")) return("active")
    if (any(states == "partial")) return("partial")
    "reference"
  }
  pathway_alpha_for <- function(path_ids, zone = c("route", "node", "label")) {
    zone <- match.arg(zone)
    state <- pathway_state_for(path_ids)
    values <- switch(
      zone,
      route = c(reference = 0.12, partial = 0.34, active = 0.82),
      node = c(reference = 0.22, partial = 0.58, active = 1.00),
      label = c(reference = 0.28, partial = 0.65, active = 1.00)
    )
    unname(values[state])
  }

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
    deoxyribose = c("deoxyribose"),
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
                      "Deoxyribose","Fucose","Rhamnose","Gluconate")
  # Classify each carbon source by evidence level
  cs_evidence <- stats::setNames(rep("none", length(cs_ids_ordered)), cs_ids_ordered)
  # Cytoplasmic support is deliberately separate from transport and GH
  # support.  A membrane hit alone must not foreground a metabolic route.
  for (csid in cs_ids_ordered) {
    pid <- tolower(csid)
    state <- unname(cytoplasm_status[pid])
    if (is.na(state)) next
    if (state == "active") cs_evidence[csid] <- "transport"
    if (state == "partial") cs_evidence[csid] <- "enzyme"
  }
  # Check for CAZy enzyme evidence (active cascades target substrates)
  # Map cascade keys to carbon source IDs
  cascade_to_cs <- c(starch = "Maltose", cellulose = "Cellobiose",
                     xylan = "Xylose", chitin = "NAG",
                     pectin = "Galacturonate", mannan = "Mannose",
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

  # Retain the full biochemical reference frame. Unsupported substrates are
  # rendered at low opacity instead of being deleted, so presence/absence is
  # visible without implying that a missing route was never assessed.

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
    tr_priority$cs_id <- unname(pathway_to_cs[tolower(tr_priority$pathway)])
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
  keep_types <- c("backbone", "ppp", "tca", "tca_shunt", "carbon_source", "ed")
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
  cyto_node_status <- stats::setNames(rep("reference", nrow(cyto_nodes)), cyto_nodes$id)
  continuity_seed_routes <- .dnmb_cct_continuity_routes()
  status_rank <- c(reference = 0L, partial = 1L, active = 2L)
  for (pid in names(cytoplasm_status)) {
    state <- unname(cytoplasm_status[pid])
    if (is.na(state) || state == "reference") next
    path_nodes <- unique(c(
      .dnmb_cct_entry_route_nodes(pid, node_ids = cyto_nodes$id),
      if (pid %in% names(continuity_seed_routes)) continuity_seed_routes[[pid]] else character()
    ))
    path_nodes <- intersect(path_nodes, names(cyto_node_status))
    if (!length(path_nodes)) next
    upgrade <- status_rank[state] > status_rank[cyto_node_status[path_nodes]]
    cyto_node_status[path_nodes[upgrade]] <- state
  }
  if (!is.null(core_steps) && nrow(core_steps)) {
    for (i in seq_len(nrow(core_steps))) {
      step_state <- core_steps$status[i]
      step_nodes <- intersect(
        c(core_steps$from[i], core_steps$to[i]),
        names(cyto_node_status)
      )
      if (!length(step_nodes)) next
      upgrade <- status_rank[step_state] > status_rank[cyto_node_status[step_nodes]]
      cyto_node_status[step_nodes[upgrade]] <- step_state
    }
  }
  exact_edge_hits <- .dnmb_cct_exact_step_edge_matches(
    step_evidence = matched_steps,
    cyto_edges = cyto_edges
  )
  if (nrow(exact_edge_hits)) {
    for (i in seq_len(nrow(exact_edge_hits))) {
      hit_rank <- suppressWarnings(as.integer(exact_edge_hits$rank[i]))
      if (!is.finite(hit_rank) || hit_rank < 2L) next
      hit_state <- if (hit_rank >= 3L) "active" else "partial"
      hit_nodes <- intersect(
        c(exact_edge_hits$mapped_from[i], exact_edge_hits$mapped_to[i]),
        names(cyto_node_status)
      )
      if (!length(hit_nodes)) next
      upgrade <- status_rank[hit_state] > status_rank[cyto_node_status[hit_nodes]]
      cyto_node_status[hit_nodes[upgrade]] <- hit_state
    }
  }
  cyto_node_alpha <- vapply(cyto_node_status, function(state) {
    switch(state, active = 1, partial = 0.58, reference = 0.22)
  }, numeric(1))

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
  sugar_lane_ids <- unique(carbon_src$sugar_type[order(carbon_src$x)])
  sugar_lane_x <- stats::setNames(vapply(sugar_lane_ids, function(st) {
    .dnmb_cct_snap_to_grid(
      mean(carbon_src$x[carbon_src$sugar_type == st]),
      step = glyco_grid_step / 2
    )
  }, numeric(1)), sugar_lane_ids)
  sugar_lane_x <- .dnmb_cct_separate_lanes(
    sugar_lane_x, min_gap = glyco_grid_step, snap_step = glyco_grid_step / 2
  )

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
  sugar_icon_r <- .dnmb_cct_sugar_icon_radius()
  sym_size <- sugar_icon_r
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
  y_mono <- y_memb + 1.4  # reserve a clear gutter above transporter buses
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
  extra_top_actual <- if (nrow(extra_layout) > 0) max(extra_layout$y, na.rm = TRUE) + 0.35 else (y_memb + 1.0)

  # Monomer aliases can share a preferred sugar lane (for example GlcA/GalA).
  # Reserve the positions after deterministic lane separation when sizing the
  # membrane, otherwise the outer alias can sit just beyond the bilayer.
  monomer_lane_ids <- .dnmb_cct_extracellular_monomer_ids()
  monomer_lane_keys <- unname(.dnmb_cct_monomer_lane_map()[monomer_lane_ids])
  monomer_bound_x <- unname(sugar_lane_x[monomer_lane_keys])
  names(monomer_bound_x) <- monomer_lane_ids
  monomer_bound_x <- .dnmb_cct_separate_lanes(
    monomer_bound_x,
    min_gap = glyco_grid_step,
    snap_step = glyco_grid_step / 2
  )

  # x_max adapts to both cytoplasm lanes and extracellular chain width
  x_max <- max(max(cs_x_pos) + 0.6, max(extra_layout$x_right) + 0.2)

  # Pathway colors (needed by cascade rendering + edge rendering)
  pw_colors <- .dnmb_cct_pathway_colors()
  node_sugar_type <- stats::setNames(cyto_nodes$sugar_type, cyto_nodes$id)
  mem_xmin <- min(
    c(
      extra_layout$x_left, extra_layout$x_center,
      carbon_src$x, cyto_nodes$x, monomer_bound_x
    ),
    na.rm = TRUE
  ) - 0.20
  mem_xmax <- max(
    c(
      extra_layout$x_right, extra_layout$x_center,
      carbon_src$x, cyto_nodes$x, monomer_bound_x
    ),
    na.rm = TRUE
  ) + 0.45

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
  zone_label_x <- mem_xmax + 0.14
  p <- p +
    ggplot2::annotate("text", x = zone_label_x, y = extra_top_actual + 0.18,
                      label = "Extracellular", hjust = 0, size = 3,
                      fontface = "italic", color = "#BDBDBD") +
    ggplot2::annotate("text", x = zone_label_x, y = y_memb,
                      label = "Membrane", hjust = 0, size = 2.5,
                      fontface = "italic", color = "#A1887F") +
    ggplot2::annotate("text", x = zone_label_x, y = y_memb - 0.5,
                      label = "Cytoplasm", hjust = 0, size = 3,
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
  mono_ids <- monomer_lane_ids
  mono_ids <- unique(c(mono_ids, setdiff(all_prods, composite_product_ids)))
  # Sugar type → monomer ID mapping
  mono_to_sugar <- .dnmb_cct_monomer_lane_map()
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
  mono_xs <- .dnmb_cct_separate_lanes(
    mono_xs, min_gap = glyco_grid_step, snap_step = glyco_grid_step / 2
  )
  mono_x_map <- stats::setNames(mono_xs, mono_ids)
  extra_target_id_map <- stats::setNames(extra_substrates$id, tolower(extra_substrates$id))

  for (mi in seq_along(mono_ids)) {
    mid <- mono_ids[mi]
    mx <- mono_xs[mi]
    mono_st <- unname(mono_to_sugar[mid])
    mono_paths <- if (!is.na(mono_st)) {
      tolower(carbon_src$id[carbon_src$sugar_type == mono_st])
    } else {
      tolower(unname(mono_to_csid[mid]))
    }
    mono_alpha <- pathway_alpha_for(mono_paths, zone = "node")
    # First render: no label (final overlay at end will add labels above route lines)
    p <- .dnmb_snfg_render_symbol_v2(
      p, mx, y_mono, mid, r = sugar_icon_r,
      label = NULL, alpha = mono_alpha
    )
  }

  # ====================================================================
  # ZONE 1: EXTRACELLULAR — SNFG chains + GH scissors on chain
  # ====================================================================
  gh_label_rows <- list()
  extra_label_rows <- list()
  scissor_specs <- list()
  gh_label_df <- NULL
  gh_ledger_df <- NULL
  extra_release_specs <- list()
  extra_route_anchor_rows <- list(
    data.frame(
      kind = "monosaccharide", id = names(mono_x_map),
      x = unname(mono_x_map), y = rep(y_mono, length(mono_x_map)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      kind = "complex_chain", id = names(extra_center_map),
      x = unname(extra_center_map), y = unname(extra_y_map[names(extra_center_map)]),
      stringsAsFactors = FALSE
    )
  )
  for (ri in seq_len(n_extra)) {
    sub   <- extra_substrates[ri, , drop = FALSE]
    ry    <- extra_y_map[sub$id]
    chain_x <- extra_x_map[sub$id]  # per-substrate x position
    chain_cx <- extra_center_map[sub$id]
    chain <- extra_chains[[sub$id]]
    bond_type <- extra_bonds[[sub$id]]
    branches <- extra_branches[[sub$id]]

    extra_label_rows[[length(extra_label_rows) + 1L]] <- data.frame(
      x = chain_cx, y = ry + 0.30, label = sub$label,
      color = "#333333", size = 1.7, fontface = "bold", hjust = 0.5,
      priority = 30, nudge_x = 0, nudge_y = 0,
      label_kind = "title", anchor_x = chain_cx, anchor_y = ry,
      stringsAsFactors = FALSE
    )

    # SNFG glycan chain
    n_mono <- length(chain)
    if (n_mono > 0) {
      bonds_vec <- if (n_mono > 1 && !is.null(bond_type)) rep(bond_type, n_mono - 1) else NULL
      chain_layers <- .dnmb_draw_sugar_chain_v2(
        x_start = chain_x, y = ry, monomers = chain,
        bonds = bonds_vec, size = sym_size, gap = chain_gap, show_label = FALSE,
        show_bond_label = FALSE)
      for (ly in chain_layers) p <- p + ly
      if (!is.null(bonds_vec) && length(bonds_vec)) {
        for (bond_i in seq_along(bonds_vec)) {
          parts <- strsplit(bonds_vec[bond_i], "-")[[1L]]
          anomer <- if (identical(parts[1L], "alpha")) "\u03b1" else "\u03b2"
          positions <- if (length(parts) >= 2L) strsplit(parts[2L], ",")[[1L]] else c("1", "4")
          bond_label <- paste0(anomer, paste(positions, collapse = "-"))
          extra_label_rows[[length(extra_label_rows) + 1L]] <- data.frame(
            x = chain_x + (bond_i - 0.5) * chain_step,
            y = ry + 0.13,
            label = bond_label,
            color = "#777777", size = 1.30, fontface = "plain", hjust = 0.5,
            priority = 8, nudge_x = 0, nudge_y = 0,
            label_kind = "bond", anchor_x = chain_x + (bond_i - 0.5) * chain_step,
            anchor_y = ry,
            stringsAsFactors = FALSE
          )
        }
      }

      # Draw branches (side chains hanging below the backbone)
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
          br_label <- gsub(",", "-", br_label)
          p <- p + ggplot2::geom_segment(
            data = data.frame(x = br_x, xend = br_x, y = ry - sym_size, yend = br_y + sym_size),
            ggplot2::aes(x=.data$x, xend=.data$xend, y=.data$y, yend=.data$yend),
            linewidth = 0.3, color = "#795548", linetype = lty, inherit.aes = FALSE)
          # Branch sugar symbol
          abbr_br <- .dnmb_snfg_abbreviation_v2(br$sugar)
          p <- .dnmb_snfg_render_symbol_v2(p, br_x, br_y, br$sugar, r = sym_size, label = abbr_br)
          extra_label_rows[[length(extra_label_rows) + 1L]] <- data.frame(
            x = br_x + 0.05, y = (ry + br_y) / 2,
            label = br_label, color = "#795548", size = 0.9,
            fontface = "plain", hjust = 0,
            priority = 10, nudge_x = 0, nudge_y = 0,
            label_kind = "branch_bond", anchor_x = br_x,
            anchor_y = (ry + br_y) / 2,
            stringsAsFactors = FALSE
          )

          # Debranching enzyme: scissors on the branch bond + GH label + arrow to monomer hub
          debranch_ghs <- br$debranch_gh
          if (!is.null(debranch_ghs) && !is.null(gh_enzymes) && nrow(gh_enzymes) > 0) {
            db_matched <- .dnmb_cct_3zone_match_gh(gh_enzymes, debranch_ghs)
            if (nrow(db_matched) > 0) {
              # Small scissors on the branch bond midpoint
              sc_y <- (ry + br_y) / 2
              scissor_specs[[length(scissor_specs) + 1L]] <- data.frame(
                x = br_x - 0.055, y = sc_y, size = 0.05, angle = 0,
                stringsAsFactors = FALSE
              )
              # Debranching enzyme label
              for (db_i in seq_len(nrow(db_matched))) {
                db_lab <- .dnmb_cct_short_gh_label(
                  db_matched$gh_family[db_i], db_matched$gene_name[db_i]
                )
                db_lab <- paste0(db_lab, "\n", db_matched$locus_tag[db_i])
                db_priority <- 12 + min(4, nrow(db_matched)) +
                  ifelse(sub$id %in% c("Pectin", "Mannan"), 2, 0)
                gh_label_rows[[length(gh_label_rows) + 1L]] <- data.frame(
                  x = br_x - 0.18, y = br_y - 0.13,
                  label = db_lab, color = "#1565C0", priority = db_priority,
                  label_kind = "gh", anchor_x = br_x - 0.055, anchor_y = sc_y,
                  locus_tag = db_matched$locus_tag[db_i],
                  gh_family = db_matched$gh_family[db_i],
                  gene_name = db_matched$gene_name[db_i],
                  substrate = sub$label,
                  pul_substrate = "",
                  stringsAsFactors = FALSE
                )
              }
              # Arrow from detached branch to monomer hub
              if (tolower(br$sugar) %in% names(mono_x_map)) {
                target_mx <- unname(mono_x_map[tolower(br$sugar)])
                branch_anchor_id <- paste("branch", sub$id, bi, sep = ":")
                extra_route_anchor_rows[[length(extra_route_anchor_rows) + 1L]] <- data.frame(
                  kind = "branch_sugar", id = branch_anchor_id,
                  x = br_x, y = br_y, stringsAsFactors = FALSE
                )
                # Leave the branch locally before taking a horizontal lane.
                # The shared lane packer below separates same-row products.
                mid_y2 <- br_y - max(0.14, sym_size * 1.6)
                branch_col <- .dnmb_cct_sugar_route_color(br$sugar, fallback = "#1565C0")
                branch_imp <- min(1.0, 0.45 + 0.18 * nrow(db_matched))
                extra_release_specs[[length(extra_release_specs) + 1L]] <- list(
                  route_id = paste("branch", sub$id, bi, tolower(br$sugar), sep = ":"),
                  source_kind = "branch_sugar",
                  source_id = branch_anchor_id,
                  target_kind = "monosaccharide",
                  target_id = tolower(br$sugar),
                  group = paste0("release:", tolower(br$sugar)),
                  priority = branch_imp,
                  points = data.frame(
                    x = c(br_x, br_x, target_mx, target_mx),
                    y = c(br_y, mid_y2, mid_y2, y_mono)
                  ),
                  color = branch_col,
                  linewidth = 0.18 + 0.18 * branch_imp,
                  alpha = 0.22 + 0.28 * branch_imp,
                  linetype = "dashed",
                  gradient = TRUE,
                  trim_start = extra_trim,
                  trim_end = 0.11
                )
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
      scissor_specs[[length(scissor_specs) + 1L]] <- data.frame(
        x = cut_x, y = ry, size = 0.067, angle = 0,
        stringsAsFactors = FALSE
      )
      scissor_anchor_id <- paste("scissor", sub$id, ckey, sep = ":")
      extra_route_anchor_rows[[length(extra_route_anchor_rows) + 1L]] <- data.frame(
        kind = "scissor", id = scissor_anchor_id,
        x = cut_x, y = ry, stringsAsFactors = FALSE
      )

      # GH label below scissors (with CGC/PUL annotation if available)
      gh_primary_y <- ry - 0.22
      if (!is.null(branches) && length(branches)) {
        branch_ys <- vapply(seq_along(branches), function(branch_i) {
          ry - sym_size * (2.5 + 1.4 * (branch_i - 1L))
        }, numeric(1))
        gh_primary_y <- min(gh_primary_y, min(branch_ys) - 0.13)
      }
      for (j in seq_len(nrow(matched))) {
        enz <- matched[j, , drop = FALSE]
        gh_lab <- .dnmb_cct_short_gh_label(enz$gh_family, enz$gene_name)
        gh_lab <- paste0(gh_lab, "\n", enz$locus_tag)
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
          gh_lab <- paste0(gh_lab, "\nPUL substrate: ", trimws(pul_sub))
        }
        main_priority <- 10 - j + ifelse(sub$id %in% c("Pectin", "Mannan"), 2, 0) +
          ifelse(grepl("GH28|GH35|GH42|GH78|GH106|GH27|GH36|GH26|GH113", enz$gh_family), 1, 0)
        gh_label_rows[[length(gh_label_rows) + 1L]] <- data.frame(
          x = cut_x, y = gh_primary_y,
          label = gh_lab, color = "#C62828",
          priority = main_priority,
          label_kind = "gh", anchor_x = cut_x, anchor_y = ry,
          locus_tag = enz$locus_tag,
          gh_family = enz$gh_family,
          gene_name = enz$gene_name,
          substrate = sub$label,
          pul_substrate = if (is.null(pul_sub)) "" else trimws(pul_sub),
          stringsAsFactors = FALSE
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
          # Exit just below the cleavage site. Same-row release routes receive
          # deterministic horizontal lanes after every cascade is collected.
          mid_y3 <- ry - max(0.14, sym_size * 1.8)
          extra_release_specs[[length(extra_release_specs) + 1L]] <- list(
            route_id = paste("cascade", sub$id, ckey, ps, sep = ":"),
            source_kind = "scissor",
            source_id = scissor_anchor_id,
            target_kind = if (ps %in% names(mono_x_map)) {
              "monosaccharide"
            } else {
              "complex_chain"
            },
            target_id = if (ps %in% names(mono_x_map)) ps else target_id,
            group = paste0("release:", ps),
            priority = route_wt,
            points = data.frame(
              x = c(exit_x, exit_x, target_x, target_x),
              y = c(ry, mid_y3, mid_y3, target_y)
            ),
            color = route_col,
            linewidth = 0.35,
            alpha = 0.45,
            linetype = "solid",
            gradient = FALSE,
            trim_start = extra_trim,
            trim_end = trim_end_val
          )
        }
      }

    }
  }

  if (length(extra_release_specs)) {
    extra_release_specs <- .dnmb_cct_validate_extracellular_routes(
      specs = extra_release_specs,
      anchors = dplyr::bind_rows(extra_route_anchor_rows),
      tolerance = 1e-6
    )
  }
  if (length(extra_release_specs)) {
    release_metadata <- data.frame(
      route_id = vapply(extra_release_specs, `[[`, character(1), "route_id"),
      priority = vapply(extra_release_specs, `[[`, numeric(1), "priority"),
      group = vapply(extra_release_specs, `[[`, character(1), "group"),
      stringsAsFactors = FALSE
    )
    release_paths <- .dnmb_cct_separate_route_lanes(
      routes = lapply(extra_release_specs, `[[`, "points"),
      metadata = release_metadata,
      lane_step = 0.065,
      max_offset = 0.22,
      horizontal_tolerance = 0.01,
      y_tolerance = 0.18,
      min_overlap = 0.14,
      connection_mode = "move_internal",
      round_radius = 0.08
    )
    for (route_index in seq_along(extra_release_specs)) {
      spec <- extra_release_specs[[route_index]]
      release_path <- release_paths[[route_index]]
      if (isTRUE(spec$gradient)) {
        for (ly in .dnmb_cct_gradient_path_layers(
          release_path, color = spec$color,
          linewidth = spec$linewidth, alpha = spec$alpha,
          linetype = spec$linetype,
          arrow_last = TRUE, arrow_length = 0.015,
          trim_start = spec$trim_start, trim_end = spec$trim_end
        )) p <- p + ly
      } else {
        p <- p + .dnmb_cct_single_path_layer(
          release_path, color = spec$color,
          linewidth = spec$linewidth, alpha = spec$alpha,
          linetype = spec$linetype,
          arrow_last = TRUE, arrow_length = 0.02,
          trim_start = spec$trim_start, trim_end = spec$trim_end
        )
      }
    }
  }

  if (length(gh_label_rows) > 0) {
    gh_label_df <- do.call(rbind, gh_label_rows)
    gh_label_df$locus_tag <- as.character(gh_label_df$locus_tag)
    gh_label_df <- gh_label_df[
      !is.na(gh_label_df$locus_tag) & nzchar(gh_label_df$locus_tag),
      , drop = FALSE
    ]
    gh_loci <- sort(unique(gh_label_df$locus_tag))
    gh_anchor_map <- stats::setNames(sprintf("G%02d", seq_along(gh_loci)), gh_loci)
    gh_label_df$label <- unname(gh_anchor_map[gh_label_df$locus_tag])
    gh_ledger_rows <- lapply(gh_loci, function(locus) {
      rows <- gh_label_df[gh_label_df$locus_tag == locus, , drop = FALSE]
      families <- unique(as.character(rows$gh_family))
      genes <- unique(as.character(rows$gene_name))
      genes <- genes[!is.na(genes) & nzchar(genes)]
      substrates <- unique(as.character(rows$substrate))
      substrates <- substrates[!is.na(substrates) & nzchar(substrates)]
      pul <- unique(as.character(rows$pul_substrate))
      pul <- pul[!is.na(pul) & nzchar(pul)]
      data.frame(
        anchor_id = unname(gh_anchor_map[locus]),
        locus_tag = locus,
        gh_family = paste(families, collapse = " | "),
        gene_name = if (length(genes)) paste(genes, collapse = " | ") else "Gene symbol not annotated",
        substrate_names = paste(substrates, collapse = " | "),
        pul_substrate = if (length(pul)) paste(pul, collapse = " | ") else "No direct PUL substrate call",
        stringsAsFactors = FALSE
      )
    })
    gh_ledger_df <- dplyr::bind_rows(gh_ledger_rows)
    gh_ledger_df$ledger_text <- paste0(
      gh_ledger_df$anchor_id, "  ", gh_ledger_df$locus_tag,
      " | dbCAN family: ", gh_ledger_df$gh_family,
      " | Gene: ", gh_ledger_df$gene_name,
      " | Cleavage context: ", gh_ledger_df$substrate_names,
      " | PUL substrate: ", gh_ledger_df$pul_substrate
    )
    gh_label_df$size <- 1.55
    gh_label_df$fontface <- "bold"
    gh_label_df$hjust <- 0.5
    gh_label_df$nudge_x <- 0
    gh_label_df$nudge_y <- 0
  }

  # ====================================================================
  # ZONE 2: MEMBRANE — Transporters aligned with carbon source x positions
  # PTS transporters (orange dashed), non-PTS (solid colored arrows)
  # ====================================================================
  # Map carbon source IDs to x positions for alignment
  transport_source_nodes <- carbon_src
  glucose_node <- cyto_nodes[
    cyto_nodes$id == "Glucose",
    c("id", "x", "y", "label", "type", "sugar_type"),
    drop = FALSE
  ]
  if (nrow(glucose_node) && !"Glucose" %in% transport_source_nodes$id) {
    transport_source_nodes <- dplyr::bind_rows(
      transport_source_nodes, glucose_node[1L, , drop = FALSE]
    )
  }
  cs_x_map <- stats::setNames(
    transport_source_nodes$x, transport_source_nodes$id
  )
  # Also map pathway names to carbon source IDs
  pathway_to_cs <- .dnmb_cct_transporter_pathway_map()

  transporter_label_df <- NULL
  transporter_ledger_df <- NULL
  transporter_entities_df <- NULL
  transporter_bus_layout <- NULL
  transporter_source_anchors <- NULL
  transporter_interior_routes <- list()
  if (!is.null(transporters) && nrow(transporters) > 0) {
    transporters$is_pts <- grepl(
      "pts|EIIC|EIIB|crr", transporters$step, ignore.case = TRUE
    ) | grepl("pts", transporters$gene_name, ignore.case = TRUE)
    transporters$cs_id <- unname(pathway_to_cs[tolower(transporters$pathway)])
    transporters <- transporters[
      !is.na(transporters$cs_id) & transporters$cs_id %in% names(cs_x_map),
      , drop = FALSE
    ]
    if (nrow(transporters) > 0) {
      if (!"step_score" %in% names(transporters)) transporters$step_score <- NA_real_
      transporters <- .dnmb_cct_annotate_transport_context(transporters, genbank_table)
      transporters$conf_rank <- unname(conf_rank[tolower(transporters$confidence)])
      transporters$conf_rank[is.na(transporters$conf_rank)] <- 0L
      transporters$kind <- mapply(
        .dnmb_cct_transporter_kind,
        transporters$step, transporters$gene_name, transporters$is_pts
      )

      positive <- transporters$conf_rank >= 2L |
        (!is.na(transporters$step_score) & transporters$step_score >= 1)
      transporters <- transporters[positive, , drop = FALSE]
    }
    if (nrow(transporters) > 0) {
      # Keep every high-confidence locus. Medium calls are ranked within each
      # substrate lane and contribute only the strongest additional locus.
      high_loci <- unique(transporters$locus_tag[
        transporters$conf_rank >= 3L |
          (!is.na(transporters$step_score) & transporters$step_score >= 2)
      ])
      medium_pool <- transporters[!transporters$locus_tag %in% high_loci, , drop = FALSE]
      if (nrow(medium_pool)) {
        medium_pool <- .dnmb_cct_select_transporters(medium_pool, max_per_lane = 1L)
      }
      selected_loci <- unique(c(high_loci, medium_pool$locus_tag))
      transporters <- transporters[transporters$locus_tag %in% selected_loci, , drop = FALSE]

      entity_result <- .dnmb_cct_transporter_entities(transporters)
      transporter_entities_df <- entity_result$entities
      transporters <- entity_result$memberships
      if (nrow(transporter_entities_df)) {
        transporter_entities_df$anchor_id <- sprintf(
          "T%02d", seq_len(nrow(transporter_entities_df))
        )
        transporters$anchor_id <- transporter_entities_df$anchor_id[
          match(transporters$locus_tag, transporter_entities_df$locus_tag)
        ]
        transporter_source_anchors <- .dnmb_cct_transporter_source_anchors(
          cs_ids = unique(transporters$cs_id), carbon_src = transport_source_nodes,
          mono_x_map = mono_x_map, y_mono = y_mono,
          extra_center_map = extra_center_map, extra_y_map = extra_y_map
        )
        source_idx <- match(
          transporters$cs_id, transporter_source_anchors$cs_id
        )
        transporters$source_x <- transporter_source_anchors$source_x[source_idx]
        transporters$source_y <- transporter_source_anchors$source_y[source_idx]
        transporters$source_kind <- transporter_source_anchors$source_kind[source_idx]
        transporters$source_trim <- ifelse(
          transporters$source_kind == "complex_chain", extra_trim, sugar_icon_r
        )
        transporters$cyto_x <- transporter_source_anchors$cyto_x[source_idx]
        transporters$cyto_y <- transporter_source_anchors$cyto_y[source_idx]
        transporters$lane_x <- transporters$source_x
        transporters <- transporters[
          is.finite(transporters$source_x) &
            is.finite(transporters$source_y) &
            is.finite(transporters$cyto_x) &
            is.finite(transporters$cyto_y),
          , drop = FALSE
        ]

        product_col <- .dnmb_pick_column(genbank_table, c("product", "description"))
        product_map <- if (!is.null(product_col) && "locus_tag" %in% names(genbank_table)) {
          stats::setNames(as.character(genbank_table[[product_col]]), genbank_table$locus_tag)
        } else {
          character()
        }
        transporter_entities_df$product <- unname(product_map[transporter_entities_df$locus_tag])
        missing_product <- is.na(transporter_entities_df$product) |
          !nzchar(transporter_entities_df$product)
        transporter_entities_df$product[missing_product] <- "Product not annotated"

        entity_members <- split(seq_len(nrow(transporters)), transporters$locus_tag)
        transporter_entities_df$anchor_x <- vapply(
          transporter_entities_df$locus_tag,
          function(locus) stats::median(unique(transporters$lane_x[entity_members[[locus]]]), na.rm = TRUE),
          numeric(1)
        )
        transporter_entities_df$confidence <- vapply(
          transporter_entities_df$locus_tag,
          function(locus) {
            rows <- transporters[entity_members[[locus]], , drop = FALSE]
            if (any(rows$confidence == "high")) "high" else "medium"
          }, character(1)
        )
        transporter_entities_df$is_pts <- vapply(
          transporter_entities_df$locus_tag,
          function(locus) any(transporters$is_pts[entity_members[[locus]]]),
          logical(1)
        )
        transporter_entities_df$kind <- vapply(
          transporter_entities_df$locus_tag,
          function(locus) {
            kinds <- transporters$kind[entity_members[[locus]]]
            kinds[order(match(kinds, c("PTS", "ABC", "MFS", "GEN")))][1]
          }, character(1)
        )
        transporter_entities_df$seg_half <- vapply(
          transporter_entities_df$locus_tag,
          function(locus) {
            rows <- transporters[entity_members[[locus]], , drop = FALSE]
            max(mapply(
              .dnmb_cct_transporter_half_span,
              rows$step, rows$gene_name, rows$is_pts
            ))
          }, numeric(1)
        )
        transporter_entities_df$anchor_y <- y_memb
        transporter_entities_df$row_draw <- 1L
        packed <- .dnmb_cct_pack_transporters_lane(
          center_x = mean(transporter_entities_df$anchor_x),
          half_spans = transporter_entities_df$seg_half,
          lane_ranks = ifelse(
            transporter_entities_df$confidence == "high", 2, 1
          ),
          desired_x = transporter_entities_df$anchor_x,
          stable_ids = transporter_entities_df$locus_tag,
          y_memb = y_memb,
          xlim = c(mem_xmin, mem_xmax),
          min_gap = 0.08,
          row_step = 0.12,
          max_rows = 3L
        )
        transporter_entities_df$anchor_x <- packed$tx
        transporter_entities_df$anchor_y <- packed$ty
        transporter_entities_df$row_draw <- packed$row_id

        transporters$tx_draw <- transporter_entities_df$anchor_x[
          match(transporters$locus_tag, transporter_entities_df$locus_tag)
        ]
        transporters$ty_draw <- transporter_entities_df$anchor_y[
          match(transporters$locus_tag, transporter_entities_df$locus_tag)
        ]

        # Many substrate lanes may converge on the same physical transporter.
        # Consolidate them into one bus per substrate and one trunk per locus;
        # the evidence ledger still retains every selected membership.
        connector_rows <- unique(transporters[, c(
          "anchor_id", "locus_tag", "cs_id", "lane_x", "source_x", "source_y",
          "source_kind", "source_trim",
          "cyto_x", "cyto_y", "tx_draw", "ty_draw", "confidence"
        ), drop = FALSE])
        connector_rows$color <- vapply(seq_len(nrow(connector_rows)), function(i) {
          member <- connector_rows[i, , drop = FALSE]
          sugar_type <- transporter_source_anchors$sugar_type[
            match(member$cs_id, transporter_source_anchors$cs_id)
          ]
          .dnmb_cct_sugar_route_color(sugar_type, fallback = "#78909C")
        }, character(1))
        transporter_bus_layout <- .dnmb_cct_transporter_bus_layout(
          connector_rows,
          y_memb = y_memb,
          bus_offset = 0.70,
          tier_step = 0.065,
          max_bus_offset = 1.12,
          source_offset = 0.76,
          glyph_clearance = 0.11,
          membrane_half_height = 0.20,
          route_gap = 0.10,
          interval_gap = 0.05,
          corner_radius = 0.08,
          suppress_redundant_medium = TRUE
        )
        for (route in transporter_bus_layout$exterior_render_routes) {
          p <- p + .dnmb_cct_single_path_layer(
            route$points, color = route$color,
            linewidth = route$linewidth, alpha = route$alpha,
            linetype = route$linetype,
            arrow_last = isTRUE(route$arrow_last),
            arrow_length = 0.015, trim_start = route$trim_start
          )
        }

        draw_members <- transporter_bus_layout$memberships[
          transporter_bus_layout$memberships$draw, , drop = FALSE
        ]
        draw_keys <- paste(draw_members$cs_id, draw_members$target_key, sep = "\r")
        connector_keys <- paste(connector_rows$cs_id, connector_rows$locus_tag, sep = "\r")
        drawn_connector_rows <- connector_rows[connector_keys %in% draw_keys, , drop = FALSE]
        transporter_interior_routes <- .dnmb_cct_transporter_interior_routes(
          drawn_connector_rows, y_memb = y_memb,
          glyph_clearance = 0.11, source_clearance = 0.11,
          membrane_half_height = 0.20, route_gap = 0.10,
          lane_step = 0.045, max_lane_offset = 0.18,
          corner_radius = 0.08
        )
        for (route in transporter_interior_routes) {
          p <- p + .dnmb_cct_single_path_layer(
            route$points, color = route$color,
            linewidth = route$linewidth, alpha = route$alpha,
            linetype = route$linetype, arrow_last = TRUE,
            arrow_length = 0.015
          )
        }

        drawn_source_ids <- unique(drawn_connector_rows$cs_id)
        explicit_sources <- transporter_source_anchors[
          transporter_source_anchors$create_glyph &
            transporter_source_anchors$cs_id %in% drawn_source_ids,
          , drop = FALSE
        ]
        for (i in seq_len(nrow(explicit_sources))) {
          src <- explicit_sources[i, , drop = FALSE]
          p <- .dnmb_cct_render_carbon_source_node(
            p, src$source_x, src$source_y, src$glyph_source_id,
            r = sugar_icon_r
          )
        }

        transporter_ledger_df <- transporter_entities_df
        transporter_ledger_df$evidence_label <- ifelse(
          transporter_ledger_df$confidence == "high",
          "high confidence",
          "medium candidate"
        )
        transporter_ledger_df$shared_label <- ifelse(
          transporter_ledger_df$shared,
          "shared predicted transporter",
          "predicted transporter"
        )
        transporter_ledger_df$ledger_text <- paste0(
          transporter_ledger_df$anchor_id, "  ",
          transporter_ledger_df$locus_tag,
          " | Product: ", transporter_ledger_df$product,
          " | GapMind models: ", transporter_ledger_df$model_names,
          " | Supported substrates: ", transporter_ledger_df$substrate_names,
          " | Evidence: ", transporter_ledger_df$evidence_label,
          " | ", transporter_ledger_df$shared_label
        )
        transporter_entities_df <- .dnmb_cct_drawn_transporter_entities(
          transporter_entities_df, transporter_bus_layout$memberships
        )

        for (i in seq_len(nrow(transporter_entities_df))) {
          entity <- transporter_entities_df[i, , drop = FALSE]
          core_color <- if (isTRUE(entity$shared)) "#455A64" else {
            first_member <- transporters[transporters$locus_tag == entity$locus_tag, , drop = FALSE][1, ]
            sugar_type <- transporter_source_anchors$sugar_type[
              match(first_member$cs_id, transporter_source_anchors$cs_id)
            ]
            .dnmb_cct_sugar_route_color(sugar_type, fallback = "#607D8B")
          }
          for (ly in .dnmb_cct_transporter_glyph_layers(
            tx = entity$anchor_x, ty = entity$anchor_y,
            half_span = entity$seg_half, core_color = core_color,
            confidence = entity$confidence, is_pts = entity$is_pts,
            kind = entity$kind
          )) p <- p + ly
          for (ly in .dnmb_cct_junction_glyph_layers(
            x = entity$anchor_x, y = entity$anchor_y,
            color = core_color,
            size = if (isTRUE(entity$shared)) 2.4 else 1.8,
            alpha = if (entity$confidence == "high") 0.9 else 0.65
          )) p <- p + ly
        }

        transporter_label_df <- data.frame(
          x = transporter_entities_df$anchor_x,
          y = transporter_entities_df$anchor_y,
          anchor_x = transporter_entities_df$anchor_x,
          anchor_y = transporter_entities_df$anchor_y,
          label = transporter_entities_df$anchor_id,
          color = ifelse(transporter_entities_df$confidence == "high", "#102A43", "#52616B"),
          size = 1.60,
          fontface = "bold",
          hjust = 0.5,
          priority = ifelse(transporter_entities_df$confidence == "high", 20, 10),
          nudge_x = 0,
          nudge_y = ifelse(transporter_entities_df$anchor_y >= y_memb, 0.26, -0.26),
          stringsAsFactors = FALSE
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
          color="#B8B8B8", linewidth=0.24, alpha=0.15,
          lineend="round", inherit.aes=FALSE)
      } else {
        # Keep orthogonal routes inside their endpoint bounding box.  The old
        # target-based stagger could overshoot short PPP edges and fold back
        # into pale circular loops around Xu-5-P/R-5-P/S-7-P/E-4-P.
        edge_pts <- .dnmb_cct_edge_points_from_row(
          ce, edge_idx = i, grid_step = glyco_grid_step
        )
        # Adaptive trim: at least node radius (0.10) so lines don't pass through symbols
        edge_total_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
        edge_trim <- min(0.12, max(0.10, edge_total_len * 0.15))
        # Simple gray line with arrow (unified width)
        p <- p + .dnmb_cct_single_path_layer(
          edge_pts, color = "#B8B8B8",
          linewidth = 0.24, alpha = 0.12,
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
  active_continuity <- intersect(
    tolower(names(continuity_entry_map)),
    names(cytoplasm_status)[cytoplasm_status %in% c("active", "partial")]
  )
  continuity_specs <- list()
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
      # Substrate colour ends at the first shared central-carbon merge. The
      # shared glycolysis/PPP reaction remains a single core/reference edge,
      # rather than becoming one parallel coloured line per carbon source.
      cont_node_ids <- .dnmb_cct_route_to_first_core_merge(cont_node_ids)
      if (is.null(cont_node_ids)) next
      cont_edges <- .dnmb_cct_route_overlay_edges(
        node_ids = cont_node_ids,
        node_x = node_x,
        node_y = node_y,
        grid_step = glyco_grid_step
      )
      if (!length(cont_edges)) next
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
      cont_state <- pathway_state_for(pid)
      cont_alpha <- pathway_alpha_for(pid, zone = "route")
      cont_lw <- if (cont_state == "active") 0.50 else 0.34
      cont_lty <- if (cont_state == "active") "solid" else "dashed"
      for (edge_index in seq_along(cont_edges)) {
        edge_points <- cont_edges[[edge_index]]
        continuity_specs[[length(continuity_specs) + 1L]] <- list(
          route_id = paste0("continuity:", pid, ":", edge_index),
          group = .dnmb_cct_route_group_for_pathway(pid, route_group_members),
          priority = cont_imp,
          points = edge_points,
          from_id = attr(edge_points, "from_id"),
          to_id = attr(edge_points, "to_id"),
          color = cont_col,
          linewidth = cont_lw,
          alpha = cont_alpha,
          linetype = cont_lty
        )
      }
    }
  }
  if (length(continuity_specs)) {
    continuity_metadata <- data.frame(
      route_id = vapply(continuity_specs, `[[`, character(1), "route_id"),
      priority = vapply(continuity_specs, `[[`, numeric(1), "priority"),
      group = vapply(continuity_specs, function(spec) {
        group <- as.character(spec$group)[1L]
        if (is.na(group) || !nzchar(group)) "ungrouped" else group
      }, character(1)),
      stringsAsFactors = FALSE
    )
    continuity_paths <- .dnmb_cct_separate_route_lanes(
      routes = lapply(continuity_specs, `[[`, "points"),
      metadata = continuity_metadata,
      lane_step = 0.04,
      max_offset = 0.16,
      horizontal_tolerance = 0.01,
      y_tolerance = 0.035,
      min_overlap = 0.20,
      connection_mode = "move_internal",
      round_radius = 0.04
    )
    for (route_index in seq_along(continuity_specs)) {
      spec <- continuity_specs[[route_index]]
      p <- p + .dnmb_cct_single_path_layer(
        continuity_paths[[route_index]],
        color = spec$color, linewidth = spec$linewidth,
        alpha = spec$alpha, linetype = spec$linetype,
        trim_start = 0.12, trim_end = 0.12
      )
    }
  }

  # Foreground every reaction that has an exact directed step mapping. This
  # includes branches that a single continuity path cannot represent, such as
  # the two products of deoxyribose-phosphate aldolase.
  if (nrow(exact_edge_hits)) {
    exact_by_edge <- split(seq_len(nrow(exact_edge_hits)), exact_edge_hits$edge_index)
    for (idx in exact_by_edge) {
      edge_index <- exact_edge_hits$edge_index[idx[1L]]
      if (!is.finite(edge_index) || edge_index < 1L || edge_index > nrow(cyto_edges)) next
      ce <- cyto_edges[edge_index, , drop = FALSE]
      edge_pts <- .dnmb_cct_edge_points_from_row(
        ce, edge_idx = edge_index, grid_step = glyco_grid_step
      )
      if (nrow(edge_pts) < 2L) next
      ranks <- suppressWarnings(as.integer(exact_edge_hits$rank[idx]))
      ranks <- ranks[is.finite(ranks)]
      if (!length(ranks) || max(ranks) < 2L) next
      edge_rank <- max(ranks)
      edge_color <- unname(pw_colors[as.character(ce$pathway[1L])])
      if (is.na(edge_color) || !nzchar(edge_color)) edge_color <- "#52616B"
      p <- p + .dnmb_cct_single_path_layer(
        edge_pts,
        color = edge_color,
        linewidth = if (edge_rank >= 3L) 0.56 else 0.38,
        alpha = if (edge_rank >= 3L) 0.88 else 0.64,
        linetype = if (edge_rank >= 3L) "solid" else "dashed",
        arrow_last = TRUE,
        arrow_length = 0.016,
        trim_start = 0.11,
        trim_end = 0.11
      )
    }
  }

  # Draw a faint biochemical reference ring, then foreground each TCA
  # reaction independently from the conservative genome-annotation scan.
  tca_nodes <- cyto_nodes[cyto_nodes$type == "tca", , drop = FALSE]
  tca_shunt_nodes <- cyto_nodes[cyto_nodes$type == "tca_shunt", , drop = FALSE]
  core_state_color <- function(state) {
    switch(state, active = "#2B6F77", partial = "#C78932", "#B8B8B8")
  }
  core_state_alpha <- function(state) {
    switch(state, active = 0.90, partial = 0.66, reference = 0.12, 0.12)
  }

  glycolysis_steps <- core_steps[core_steps$pathway == "Glycolysis", , drop = FALSE]
  if (nrow(glycolysis_steps)) {
    for (i in seq_len(nrow(glycolysis_steps))) {
      step <- glycolysis_steps[i, , drop = FALSE]
      if (identical(step$status, "reference") ||
          !all(c(step$from, step$to) %in% names(node_x))) next
      endpoint_pairs <- list(c(step$from, step$to))
      if (identical(step$reaction_id, "C14") && "DHAP" %in% names(node_x)) {
        endpoint_pairs[[2L]] <- c(step$from, "DHAP")
      }
      for (pair in endpoint_pairs) {
        x1 <- unname(node_x[pair[1]]); y1 <- unname(node_y[pair[1]])
        x2 <- unname(node_x[pair[2]]); y2 <- unname(node_y[pair[2]])
        edge_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
        trim <- min(0.12, edge_len * 0.22)
        frac <- if (edge_len > 0) trim / edge_len else 0
        seg <- data.frame(
          x = x1 + (x2 - x1) * frac,
          y = y1 + (y2 - y1) * frac,
          xend = x2 - (x2 - x1) * frac,
          yend = y2 - (y2 - y1) * frac
        )
        p <- p + ggplot2::geom_segment(
          data = seg,
          ggplot2::aes(x = .data$x, y = .data$y,
                       xend = .data$xend, yend = .data$yend),
          color = core_state_color(step$status),
          linewidth = if (step$status == "active") 0.60 else 0.46,
          alpha = core_state_alpha(step$status),
          linetype = if (step$status == "partial") "dashed" else "solid",
          lineend = "round", inherit.aes = FALSE
        )
      }
    }
  }

  if (nrow(tca_nodes) >= 3) {
    tca_cx <- mean(tca_nodes$x)
    tca_cy <- mean(tca_nodes$y)
    tca_r  <- max(sqrt((tca_nodes$x - tca_cx)^2 + (tca_nodes$y - tca_cy)^2))
    theta_seq <- seq(0, 2 * pi, length.out = 200)
    tca_circle <- data.frame(x = tca_cx + tca_r * cos(theta_seq),
                              y = tca_cy + tca_r * sin(theta_seq))
    p <- p + ggplot2::geom_path(
      data = tca_circle,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "#B8B8B8", linewidth = 0.70, alpha = 0.16,
      linetype = "solid", inherit.aes = FALSE)

    tca_step_rows <- core_steps[core_steps$pathway == "TCA cycle", , drop = FALSE]
    for (i in seq_len(nrow(tca_step_rows))) {
      step <- tca_step_rows[i, , drop = FALSE]
      if (!all(c(step$from, step$to) %in% names(node_x))) next
      a_from <- atan2(unname(node_y[step$from]) - tca_cy,
                      unname(node_x[step$from]) - tca_cx)
      a_to <- atan2(unname(node_y[step$to]) - tca_cy,
                    unname(node_x[step$to]) - tca_cx)
      delta <- ((a_to - a_from + pi) %% (2 * pi)) - pi
      trim_angle <- min(0.055, abs(delta) * 0.18)
      theta <- seq(
        a_from + sign(delta) * trim_angle,
        a_from + delta - sign(delta) * trim_angle,
        length.out = 28
      )
      arc <- data.frame(x = tca_cx + tca_r * cos(theta),
                        y = tca_cy + tca_r * sin(theta))
      state <- step$status[1]
      p <- p + ggplot2::geom_path(
        data = arc,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = core_state_color(state),
        linewidth = if (state == "active") 0.82 else if (state == "partial") 0.64 else 0.34,
        alpha = core_state_alpha(state),
        linetype = if (state == "partial") "dashed" else "solid",
        lineend = "round", inherit.aes = FALSE
      )
    }
  }

  shunt_steps <- core_steps[core_steps$pathway == "Glyoxylate shunt", , drop = FALSE]
  if (nrow(tca_shunt_nodes) > 0 && nrow(shunt_steps)) {
    for (i in seq_len(nrow(shunt_steps))) {
      step <- shunt_steps[i, , drop = FALSE]
      if (!all(c(step$from, step$to) %in% names(node_x))) next
      state <- step$status[1]
      p <- p + ggplot2::geom_curve(
        data = data.frame(
          x = unname(node_x[step$from]), y = unname(node_y[step$from]),
          xend = unname(node_x[step$to]), yend = unname(node_y[step$to])
        ),
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        color = core_state_color(state),
        linewidth = if (state == "active") 0.62 else 0.46,
        linetype = if (state == "active") "solid" else "dashed",
        alpha = core_state_alpha(state),
        curvature = if (step$reaction_id == "C09") 0.12 else -0.12,
        lineend = "round", inherit.aes = FALSE
      )
    }
  }
  core_direct_label_df <- .dnmb_cct_core_direct_labels(
    core_steps = core_steps,
    display_hits = core_display_hits,
    nodes = cyto_nodes
  )

  # Draw SNFG nodes for all cytoplasm node types (symbols only, no labels)
  node_r <- sugar_icon_r
  for (ntype in c("backbone", "ppp", "entry_intermediate", "tca", "tca_shunt", "ed", "pyruvate_branch")) {
    nset <- cyto_nodes[cyto_nodes$type == ntype, , drop = FALSE]
    for (i in seq_len(nrow(nset))) {
      nd <- nset[i, , drop = FALSE]
      # Symbol only — full name labels are rendered separately below
      node_alpha_i <- unname(cyto_node_alpha[nd$id])
      if (!is.finite(node_alpha_i)) node_alpha_i <- 0.22
      p <- .dnmb_snfg_render_symbol_v2(
        p, nd$x, nd$y, nd$sugar_type, r = node_r, label = NULL,
        alpha = node_alpha_i
      )
    }
  }

  # (backbone/ppp/entry_intermediate labels rendered in final overlay section below)

  # Carbon source nodes — CONVERGE: multiple substrates → 1 hub per sugar_type
  # 1) Group by sugar_type, find hub position (mean x of group)
  # 2) Draw converging arrows from each substrate label → hub SNFG symbol
  # 3) Hub connects onward to entry intermediate
  cs_groups <- split(seq_len(nrow(carbon_src)), carbon_src$sugar_type)
  hub_positions <- list()  # sugar_type -> list(x, y)
  route_label_rows <- list()
  metabolic_ledger_df <- NULL
  hub_guide_specs <- list()
  for (st in names(cs_groups)) {
    idx <- cs_groups[[st]]
    hub_x <- if (st %in% names(sugar_lane_x)) {
      unname(sugar_lane_x[st])
    } else {
      .dnmb_cct_snap_to_grid(mean(carbon_src$x[idx]), step = glyco_grid_step)
    }
    hub_y <- .dnmb_cct_snap_to_grid(carbon_src$y[idx[1]], step = glyco_grid_step)
    hub_path_ids <- tolower(carbon_src$id[idx])
    hub_alpha <- pathway_alpha_for(hub_path_ids, zone = "node")
    source_at_hub <- any(
      abs(carbon_src$x[idx] - hub_x) < 0.02 &
        abs(carbon_src$y[idx] - hub_y) < 0.02
    )
    hub_positions[[st]] <- list(
      x = hub_x, y = hub_y, alpha = hub_alpha, path_ids = hub_path_ids,
      source_at_hub = source_at_hub
    )
    # A source exactly on the convergence coordinate doubles as the hub. Its
    # substrate-specific glyph is drawn in the final overlay instead of
    # stacking a generic family glyph underneath it.
    if (!source_at_hub) {
      p <- .dnmb_snfg_render_symbol_v2(
        p, hub_x, hub_y, st, r = sugar_icon_r,
        label = NULL, alpha = hub_alpha
      )
    }
    # Draw converging lines from each substrate to hub
    for (i in idx) {
      nd <- carbon_src[i, , drop = FALSE]
      if (abs(nd$x - hub_x) > 0.2) {
        # Keep source-to-hub convergence in its own cytoplasmic band instead
        # of letting it merge with the membrane centerline and T-ID glyphs.
        guide_y <- .dnmb_cct_snap_to_grid(
          hub_y + glyco_grid_step / 2,
          step = glyco_grid_step / 2
        )
        guide_col <- route_color_for(path_ids = tolower(nd$id), sugar_type = st)
        guide_path <- data.frame(
          x = c(nd$x, nd$x, hub_x, hub_x),
          y = c(nd$y + 0.15, guide_y, guide_y, hub_y)
        )
        guide_state <- pathway_state_for(tolower(nd$id))
        hub_guide_specs[[length(hub_guide_specs) + 1L]] <- list(
          route_id = paste("source", st, nd$id, sep = ":"),
          group = paste0("hub:", st),
          priority = if (guide_state == "active") 2 else if (guide_state == "partial") 1 else 0,
          points = guide_path,
          color = guide_col,
          state = guide_state,
          supported = !is.null(route_best_match(tolower(nd$id))) && guide_state != "reference",
          path_id = nd$id
        )
      }
    }
  }
  if (length(hub_guide_specs)) {
    guide_metadata <- data.frame(
      route_id = vapply(hub_guide_specs, `[[`, character(1), "route_id"),
      priority = vapply(hub_guide_specs, `[[`, numeric(1), "priority"),
      group = vapply(hub_guide_specs, `[[`, character(1), "group"),
      stringsAsFactors = FALSE
    )
    hub_guide_paths <- .dnmb_cct_separate_route_lanes(
      routes = lapply(hub_guide_specs, `[[`, "points"),
      metadata = guide_metadata,
      lane_step = 0.04,
      max_offset = 0.12,
      horizontal_tolerance = 0.01,
      y_tolerance = 0.04,
      min_overlap = 0.16,
      connection_mode = "move_internal",
      round_radius = 0.04
    )
    for (route_index in seq_along(hub_guide_specs)) {
      spec <- hub_guide_specs[[route_index]]
      guide_path <- hub_guide_paths[[route_index]]
      p <- p + .dnmb_cct_single_path_layer(
        guide_path, color = .dnmb_cct_reference_gray("light"),
        linewidth = 0.16, alpha = 0.12
      )
      if (isTRUE(spec$supported)) {
        for (ly in .dnmb_cct_gradient_path_layers(
          guide_path, color = spec$color,
          linewidth = if (spec$state == "active") 0.34 else 0.22,
          alpha = pathway_alpha_for(spec$path_id, zone = "route"),
          linetype = if (spec$state == "active") "solid" else "dashed"
        )) p <- p + ly
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
    target_rows <- dplyr::bind_rows(lapply(cs_of_type, function(csid) {
      targets <- .dnmb_cct_initial_entry_targets(
        tolower(csid), node_ids = names(target_x)
      )
      if (!length(targets)) {
        fallback <- unname(entry_map[tolower(csid)])
        targets <- fallback[!is.na(fallback) & fallback %in% names(target_x)]
      }
      if (!length(targets)) return(NULL)
      data.frame(
        cs_id = rep(csid, length(targets)),
        target_id = targets,
        stringsAsFactors = FALSE
      )
    }))
    if (!nrow(target_rows)) next
    target_groups <- split(target_rows$cs_id, target_rows$target_id)
    for (entry_target in names(target_groups)) {
      route_csids <- sort(unique(target_groups[[entry_target]]))
      route_nodes <- entry_target
      rep_csid <- route_csids[1]
      # The hub owns only the connection to the first intracellular node.
      # Substrate-specific continuity then carries that node to the first core
      # merge, avoiding a duplicate foreground path over the same reactions.
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
    route_obstacle_x <- c(target_x, carbon_src$x, hub_route_df$hub_x)
    route_obstacle_y <- c(target_y, carbon_src$y, hub_route_df$hub_y)
    hub_paths <- vector("list", nrow(hub_route_df))
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
        lane_count = rr$target_count,
        obstacle_x = route_obstacle_x,
        obstacle_y = route_obstacle_y
      )
      if (is.null(main_path)) {
        main_path <- .dnmb_cct_hub_entry_path(
          hub_x = rr$hub_x,
          hub_y = rr$hub_y,
          target_x = rr$target_x,
          target_y = rr$target_y,
          lane_rank = rr$target_rank,
          target_count = rr$target_count,
          grid_step = glyco_grid_step,
          rounded = FALSE
        )
      }
      hub_paths[[i]] <- main_path
    }
    hub_metadata <- data.frame(
      route_id = paste0(
        "hub:", seq_len(nrow(hub_route_df)), ":",
        hub_route_df$path_ids, ":", hub_route_df$target_id
      ),
      priority = hub_route_df$route_imp,
      group = paste0("target:", hub_route_df$target_id),
      stringsAsFactors = FALSE
    )
    hub_paths <- .dnmb_cct_separate_route_lanes(
      routes = hub_paths,
      metadata = hub_metadata,
      lane_step = 0.045,
      max_offset = 0.14,
      horizontal_tolerance = 0.01,
      y_tolerance = 0.04,
      min_overlap = 0.18,
      connection_mode = "move_internal",
      round_radius = 0.04
    )

    for (i in seq_len(nrow(hub_route_df))) {
      rr <- hub_route_df[i, , drop = FALSE]
      main_path <- hub_paths[[i]]
      # Keep the biochemical reference route faint, then foreground only
      # pathways with positive intracellular evidence.
      hub_state <- pathway_state_for(hit_ids <- unlist(strsplit(rr$path_ids, "|", fixed = TRUE)))
      p <- p + .dnmb_cct_single_path_layer(
        main_path, color = "#B8B8B8",
        linewidth = 0.22,
        alpha = 0.12,
        arrow_last = TRUE, arrow_length = 0.02,
        trim_end = 0.08
      )
      if (hub_state != "reference") {
        p <- p + .dnmb_cct_single_path_layer(
          main_path, color = rr$route_col,
          linewidth = if (hub_state == "active") 0.48 else 0.32,
          alpha = pathway_alpha_for(hit_ids, zone = "route"),
          linetype = if (hub_state == "active") "solid" else "dashed",
          arrow_last = TRUE, arrow_length = 0.02,
          trim_end = 0.08
        )
      }
      # The SNFG metabolite at the route target already marks the junction.
      # A second circular junction glyph creates a false ring behind stars such
      # as Xu-5-P and R-5-P, so metabolic targets use the native symbol only.
      # Step labels are added only through the exact directed-edge mapping
      # below. A pathway-level best hit is not a biochemical reaction anchor.
    }
  }

  # Map labels require an exact, directed (pathway, step) -> reaction-edge
  # mapping. There is intentionally no nearest-node or pathway-level fallback.
  route_label_df <- .dnmb_cct_exact_step_labels(
    step_evidence = matched_steps,
    cyto_edges = cyto_edges
  )
  if (!nrow(route_label_df)) route_label_df <- NULL

  # The ledger is deliberately broader than the map: every high/medium
  # intracellular hit remains reviewable even when no displayed edge can
  # localize it without ambiguity.
  ledger_hits <- matched_steps[
    !is.na(matched_steps$locus_tag) & nzchar(matched_steps$locus_tag) &
      matched_steps$rank >= 2L,
    , drop = FALSE
  ]
  if (nrow(ledger_hits)) {
    localized_matches <- .dnmb_cct_exact_step_edge_matches(
      step_evidence = ledger_hits,
      cyto_edges = cyto_edges
    )
    localized_keys <- if (nrow(localized_matches)) {
      unique(paste(
        localized_matches$.path_step_key,
        localized_matches$locus_tag,
        sep = "\r"
      ))
    } else {
      character()
    }
    ledger_hits$.path_step_key <- .dnmb_cct_path_step_key(
      ledger_hits$pathway_id, ledger_hits$step_id
    )
    ledger_hits$.localized <- paste(
      ledger_hits$.path_step_key, ledger_hits$locus_tag, sep = "\r"
    ) %in% localized_keys
    metabolic_loci <- sort(unique(as.character(ledger_hits$locus_tag)))
    product_col <- .dnmb_pick_column(genbank_table, c("product", "description"))
    metabolic_product_map <- if (!is.null(product_col) && "locus_tag" %in% names(genbank_table)) {
      stats::setNames(as.character(genbank_table[[product_col]]), genbank_table$locus_tag)
    } else {
      character()
    }
    metabolic_ledger_df <- dplyr::bind_rows(lapply(metabolic_loci, function(locus) {
      rows <- ledger_hits[ledger_hits$locus_tag == locus, , drop = FALSE]
      steps <- unique(as.character(rows$step_id))
      pathways <- unique(as.character(rows$pathway_id))
      steps <- steps[!is.na(steps) & nzchar(steps)]
      pathways <- pathways[!is.na(pathways) & nzchar(pathways)]
      product <- unname(metabolic_product_map[locus])
      if (is.na(product) || !nzchar(product)) product <- "Product not annotated"
      evidence <- if (any(rows$confidence == "high" | rows$rank >= 3L)) {
        "high confidence"
      } else {
        "medium confidence"
      }
      data.frame(
        locus_tag = locus,
        evidence_label = evidence,
        ledger_text = paste0(
          locus, " | Product: ", product,
          " | GapMind intracellular steps: ", paste(steps, collapse = " | "),
          " | Supported pathways: ", paste(pathways, collapse = " | "),
          " | Evidence: ", evidence,
          " | Map localization: ",
          if (any(rows$.localized)) "exact directed reaction edge" else "ledger only"
        ),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    metabolic_ledger_df <- NULL
  }

  # ---- Final extracellular symbol overlay: keep polymer/branch nodes above route lines ----
  extra_obstacle_rows <- list()
  for (ri in seq_len(n_extra)) {
    sub <- extra_substrates[ri, , drop = FALSE]
    ry <- extra_y_map[sub$id]
    chain_x <- extra_x_map[sub$id]
    chain <- extra_chains[[sub$id]]
    if (length(chain) > 0) {
      chain_xs <- chain_x + (seq_along(chain) - 1) * chain_step
      extra_obstacle_rows[[length(extra_obstacle_rows) + 1L]] <- .dnmb_cct_obstacle_perimeter(
        x = chain_xs, y = rep(ry, length(chain_xs)),
        rx = sym_size * 1.15, ry = sym_size * 1.15
      )
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
        extra_obstacle_rows[[length(extra_obstacle_rows) + 1L]] <- .dnmb_cct_obstacle_perimeter(
          x = br_x, y = br_y, rx = sym_size * 1.15, ry = sym_size * 1.15
        )
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
    mono_st <- unname(mono_to_sugar[mid])
    mono_paths <- if (!is.na(mono_st)) {
      tolower(carbon_src$id[carbon_src$sugar_type == mono_st])
    } else {
      tolower(unname(mono_to_csid[mid]))
    }
    p <- .dnmb_snfg_render_symbol_v2(
      p, mx, y_mono, mid, r = sugar_icon_r, label = abbr,
      alpha = pathway_alpha_for(mono_paths, zone = "node")
    )
  }

  for (st in names(hub_positions)) {
    hp <- hub_positions[[st]]
    if (!isTRUE(hp$source_at_hub)) {
      abbr <- .dnmb_snfg_abbreviation_v2(st)
      p <- .dnmb_snfg_render_symbol_v2(
        p, hp$x, hp$y, st, r = sugar_icon_r, label = abbr,
        alpha = hp$alpha %||% 0.22
      )
    }
  }

  # Draw every cytoplasmic source after all route lines. Source abbreviations
  # are embedded in their glyphs, including composite disaccharide symbols.
  for (i in seq_len(nrow(carbon_src))) {
    nd <- carbon_src[i, , drop = FALSE]
    p <- .dnmb_cct_render_carbon_source_node(
      p, nd$x, nd$y, nd$id, r = sugar_icon_r,
      alpha = pathway_alpha_for(tolower(nd$id), zone = "node")
    )
  }

  for (ntype in c("backbone", "ppp", "entry_intermediate", "tca", "tca_shunt", "ed", "pyruvate_branch")) {
    nset <- cyto_nodes[cyto_nodes$type == ntype, , drop = FALSE]
    for (i in seq_len(nrow(nset))) {
      nd <- nset[i, , drop = FALSE]
      inside_label <- .dnmb_cct_cytoplasmic_inside_label(
        nd$id, nd$label, nd$type, nd$sugar_type
      )
      node_alpha_i <- unname(cyto_node_alpha[nd$id])
      if (!is.finite(node_alpha_i)) node_alpha_i <- 0.22
      p <- .dnmb_snfg_render_symbol_v2(
        p, nd$x, nd$y, nd$sugar_type, r = node_r,
        label = if (nzchar(inside_label)) inside_label else NULL,
        alpha = node_alpha_i
      )
    }
  }

  # Cleavage marks are annotations, so draw them only after every route and
  # extracellular symbol. This keeps the blades, handles, and pivot visible.
  scissor_df <- if (length(scissor_specs)) dplyr::bind_rows(scissor_specs) else data.frame()
  if (nrow(scissor_df)) {
    for (i in seq_len(nrow(scissor_df))) {
      scissor_layers <- .dnmb_scissors_grob_v2(
        scissor_df$x[i], scissor_df$y[i],
        size = scissor_df$size[i], angle = scissor_df$angle[i]
      )
      for (layer in scissor_layers) p <- p + layer
    }
  }

  # ---- External labels for non-sugar metabolites only ----
  all_metab_labels <- data.frame(
    x = numeric(0), y = numeric(0), label = character(0),
    nudge_x = numeric(0), nudge_y = numeric(0),
    alpha = numeric(0),
    stringsAsFactors = FALSE
  )
  # Sugar and sugar-phosphate abbreviations are already inside their glyphs.
  for (ntype in c("backbone", "ppp", "entry_intermediate", "ed", "pyruvate_branch")) {
    nset <- cyto_nodes[cyto_nodes$type == ntype, , drop = FALSE]
    if (nrow(nset) > 0) {
      inside_labels <- mapply(
        .dnmb_cct_cytoplasmic_inside_label,
        nset$id, nset$label, nset$type, nset$sugar_type,
        USE.NAMES = FALSE
      )
      nset <- nset[!nzchar(inside_labels), , drop = FALSE]
    }
    if (nrow(nset) > 0) {
      all_metab_labels <- rbind(all_metab_labels, data.frame(
        x = nset$x, y = nset$y, label = nset$label,
        nudge_x = 0, nudge_y = 0.16,
        alpha = unname(cyto_node_alpha[nset$id]),
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
        alpha = unname(cyto_node_alpha[tca_ring$id]),
        stringsAsFactors = FALSE
      ))
    }
  }
  if (nrow(tca_shunt_nodes) > 0) {
    all_metab_labels <- rbind(all_metab_labels, data.frame(
      x = tca_shunt_nodes$x,
      y = tca_shunt_nodes$y,
      label = tca_shunt_nodes$label,
      nudge_x = 0,
      nudge_y = 0,
      alpha = unname(cyto_node_alpha[tca_shunt_nodes$id]),
      stringsAsFactors = FALSE
    ))
  }
  if (nrow(all_metab_labels)) {
    all_metab_labels$color <- "#333333"
    all_metab_labels$size <- 1.50
    all_metab_labels$fontface <- "bold"
    all_metab_labels$hjust <- 0.5
    all_metab_labels$priority <- 20
  }

  extra_label_df <- dplyr::bind_rows(
    extra_label_rows,
    if (!is.null(gh_label_df)) list(gh_label_df) else list()
  )
  extra_label_df <- .dnmb_cct_layout_external_labels(
    extra_label_df,
    xlim = c(mem_xmin, mem_xmax),
    row_step = 0.15,
    max_aux_rows = 2L,
    min_gap = 0.08
  )
  cyto_label_df <- dplyr::bind_rows(all_metab_labels, route_label_df)

  extra_obstacles <- dplyr::bind_rows(
    extra_obstacle_rows,
    .dnmb_cct_obstacle_perimeter(
      x = mono_xs, y = rep(y_mono, length(mono_xs)), rx = 0.12, ry = 0.12
    ),
    if (nrow(scissor_df)) {
      .dnmb_cct_obstacle_perimeter(
        x = scissor_df$x, y = scissor_df$y,
        rx = scissor_df$size * 1.25, ry = scissor_df$size * 0.75
      )
    } else NULL
  )
  membrane_obstacles <- if (!is.null(transporter_entities_df) &&
      nrow(transporter_entities_df) &&
      all(c("anchor_x", "anchor_y", "seg_half") %in% names(transporter_entities_df))) {
    .dnmb_cct_obstacle_perimeter(
      x = transporter_entities_df$anchor_x,
      y = transporter_entities_df$anchor_y,
      rx = pmax(0.12, transporter_entities_df$seg_half + 0.03), ry = 0.10
    )
  } else {
    data.frame(x = numeric(), y = numeric())
  }
  hub_obstacles <- dplyr::bind_rows(lapply(hub_positions, function(hp) {
    .dnmb_cct_obstacle_perimeter(hp$x, hp$y, rx = 0.12, ry = 0.12)
  }))
  cyto_obstacles <- dplyr::bind_rows(
    .dnmb_cct_obstacle_perimeter(cyto_nodes$x, cyto_nodes$y, rx = 0.12, ry = 0.12),
    .dnmb_cct_obstacle_perimeter(carbon_src$x, carbon_src$y, rx = 0.12, ry = 0.12),
    hub_obstacles
  )

  extra_leader_df <- .dnmb_cct_external_label_leaders(extra_label_df)
  if (nrow(extra_leader_df)) {
    p <- p + ggplot2::geom_path(
      data = extra_leader_df,
      ggplot2::aes(
        x = .data$x, y = .data$y, group = .data$group
      ),
      color = extra_leader_df$color,
      linewidth = 0.12, alpha = 0.52, lineend = "round",
      show.legend = FALSE, inherit.aes = FALSE
    )
  }
  if (nrow(extra_label_df)) {
    p <- p + ggplot2::geom_text(
      data = extra_label_df,
      ggplot2::aes(x = .data$x_lab, y = .data$y_lab, label = .data$label),
      color = extra_label_df$color,
      size = extra_label_df$size,
      fontface = extra_label_df$fontface,
      hjust = extra_label_df$hjust,
      alpha = if ("alpha" %in% names(extra_label_df)) extra_label_df$alpha else 1,
      lineheight = 0.9,
      check_overlap = FALSE,
      inherit.aes = FALSE
    )
  }

  transporter_label_df <- .dnmb_cct_layout_transporter_labels(
    transporter_label_df,
    xlim = c(mem_xmin, mem_xmax),
    y_memb = y_memb,
    base_offset = 0.26,
    track_step = 0.11,
    max_aux_tracks = 2L,
    min_gap = 0.08
  )
  transporter_leader_df <- .dnmb_cct_transporter_label_leaders(
    transporter_label_df
  )
  if (nrow(transporter_leader_df)) {
    p <- p + ggplot2::geom_path(
      data = transporter_leader_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      color = transporter_leader_df$color,
      linewidth = 0.12, alpha = 0.58, lineend = "round",
      show.legend = FALSE, inherit.aes = FALSE
    )
  }
  membrane_text_layer <- .dnmb_cct_transporter_text_layer(transporter_label_df)
  if (!is.null(membrane_text_layer)) p <- p + membrane_text_layer

  cyto_text_layer <- .dnmb_cct_final_text_layer(
    cyto_label_df,
    obstacles = cyto_obstacles,
    xlim = c(mem_xmin, mem_xmax),
    ylim = c(y_bot - 0.40, y_memb - 0.28),
    seed = 44L,
    force = 1.35,
    force_pull = 2.5
  )
  if (!is.null(cyto_text_layer)) p <- p + cyto_text_layer

  # ====================================================================
  # LEGENDS
  # ====================================================================
  # Build the legend from every sugar type actually rendered in this map.
  used_snfg_types <- .dnmb_cct_used_snfg_types(
    cyto_nodes = cyto_nodes,
    source_ids = carbon_src$id,
    monomer_ids = mono_ids,
    extra_chains = extra_chains,
    extra_branches = extra_branches
  )
  snfg_leg <- .dnmb_cct_snfg_legend_data(used_snfg_types)
  snfg_legend_width <- max(1.02, 0.052 * max(nchar(snfg_leg$label)) + 0.36)
  # Pull the legend into the unused lower-right interior space.
  leg_x <- max(cyto_nodes$x) - 0.28
  # y aligned slightly above the TCA circle center so the block sits tighter.
  tca_cy_leg <- if (nrow(tca_nodes) > 0) mean(tca_nodes$y) + 0.20 else 0.7
  sp <- max(
    if (nrow(snfg_leg) > 14L) 0.165 else 0.19,
    sugar_icon_r * 2.35
  )
  n_snfg <- nrow(snfg_leg)
  leg_y_top <- tca_cy_leg + (n_snfg / 2) * sp
  p <- p + ggplot2::annotate("text", x = leg_x, y = leg_y_top + 0.22,
    hjust = 0, size = 1.9, fontface = "bold", color = "#333333",
    label = "SNFG / sugar symbols")
  for (i in seq_len(n_snfg)) {
    ly <- leg_y_top - (i - 1) * sp
    p <- .dnmb_snfg_render_symbol_v2(
      p, leg_x, ly, snfg_leg$sugar[i], r = sugar_icon_r
    )
    p <- p + ggplot2::annotate("text", x = leg_x + 0.19, y = ly,
      label = snfg_leg$label[i], hjust = 0, size = 1.35, color = "#555555")
  }
  # Bond legend below SNFG
  bleg_y <- leg_y_top - n_snfg * sp - 0.15
  p <- p +
    ggplot2::geom_segment(
      data = data.frame(x = leg_x, xend = leg_x + 0.24, y = bleg_y, yend = bleg_y),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 0.75, color = "#444444", linetype = "solid", inherit.aes = FALSE) +
    ggplot2::annotate("text", x = leg_x + 0.33, y = bleg_y,
      label = "\u03b1 bond", hjust = 0, size = 1.45, color = "#555555")
  p <- p +
    ggplot2::geom_segment(
      data = data.frame(x = leg_x, xend = leg_x + 0.24, y = bleg_y - sp, yend = bleg_y - sp),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      linewidth = 0.75, color = "#444444", linetype = "dashed", inherit.aes = FALSE) +
    ggplot2::annotate("text", x = leg_x + 0.33, y = bleg_y - sp,
      label = "\u03b2 bond", hjust = 0, size = 1.45, color = "#555555")
  # Scissors
  sc_ly <- .dnmb_scissors_grob_v2(leg_x + 0.09, bleg_y - 2 * sp, size = 0.10)
  for (l in sc_ly) p <- p + l
  p <- p + ggplot2::annotate("text", x = leg_x + 0.33, y = bleg_y - 2 * sp,
    label = "GH cleavage", hjust = 0, size = 1.45, color = "#555555")
  enz_leg_y <- bleg_y - 3 * sp
  p <- p + ggplot2::annotate(
    "text", x = leg_x, y = enz_leg_y,
    label = "Pathway evidence", hjust = 0, size = 1.75,
    fontface = "bold", color = "#333333"
  )
  evidence_legend <- data.frame(
    y = enz_leg_y - sp * seq_len(3),
    color = c("#2A7F62", "#7F8C8D", "#B8B8B8"),
    alpha = c(0.85, 0.55, 0.25),
    linewidth = c(0.75, 0.55, 0.35),
    linetype = c("solid", "dashed", "solid"),
    label = c("Supported pathway", "Partial evidence", "Reference / not supported"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(evidence_legend))) {
    ev <- evidence_legend[i, , drop = FALSE]
    p <- p + ggplot2::geom_segment(
      data = data.frame(
        x = leg_x, xend = leg_x + 0.24,
        y = ev$y, yend = ev$y
      ),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      color = ev$color, alpha = ev$alpha, linewidth = ev$linewidth,
      linetype = ev$linetype, inherit.aes = FALSE
    ) + ggplot2::annotate(
      "text", x = leg_x + 0.33, y = ev$y,
      label = ev$label, hjust = 0, size = 1.35, color = "#555555"
    )
  }

  # Central-carbon loci are written beside the supported reaction itself.
  # This final overlay keeps them readable without a separate lookup panel.
  if (nrow(core_direct_label_df)) {
    direct_colors <- vapply(
      core_direct_label_df$status, core_state_color, character(1)
    )
    direct_alpha <- pmax(
      0.55,
      vapply(core_direct_label_df$status, core_state_alpha, numeric(1))
    )
    p <- p + ggplot2::geom_text(
      data = core_direct_label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = core_direct_label_df$hjust,
      vjust = core_direct_label_df$vjust,
      size = core_direct_label_df$size,
      color = direct_colors,
      alpha = direct_alpha,
      lineheight = 0.90,
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  # ====================================================================
  # THEME & OUTPUT
  # ====================================================================
  y_bottom <- y_bot - 0.5
  y_top_plot <- max(extra_top_actual + 0.55, max(cyto_nodes$y) + 1.15)
  x_left <- mem_xmin - 0.05
  x_right <- max(zone_label_x + 1.30, leg_x + snfg_legend_width)

  core_ledger_df <- core_steps[core_steps$status != "reference", , drop = FALSE]
  if (nrow(core_ledger_df)) {
    core_ledger_df$ledger_text <- paste0(
      core_ledger_df$reaction_label,
      " | Components/loci: ", core_ledger_df$components,
      " | Products: ", core_ledger_df$product_evidence,
      " | Status: ", core_ledger_df$status,
      " | Evidence: conservative gene/product/EC/KO annotation agreement"
    )
  }
  evidence_ledger_df <- dplyr::bind_rows(
    if (!is.null(transporter_ledger_df) && nrow(transporter_ledger_df)) {
      data.frame(
        ledger_text = transporter_ledger_df$ledger_text,
        ledger_type = "Transporter",
        ledger_priority = ifelse(
          transporter_ledger_df$confidence == "high",
          ifelse(transporter_ledger_df$shared, 4.0, 3.5),
          ifelse(transporter_ledger_df$shared, 3.0, 2.5)
        ),
        ledger_color = ifelse(
          transporter_ledger_df$confidence == "high", "#102A43", "#52616B"
        ),
        ledger_fontface = ifelse(
          transporter_ledger_df$confidence == "high", "bold", "plain"
        ),
        stringsAsFactors = FALSE
      )
    } else NULL,
    if (!is.null(gh_ledger_df) && nrow(gh_ledger_df)) {
      data.frame(
        ledger_text = gh_ledger_df$ledger_text,
        ledger_type = "GH cleavage",
        ledger_priority = 1.5,
        ledger_color = "#9B1C1C",
        ledger_fontface = "plain",
        stringsAsFactors = FALSE
      )
    } else NULL,
    if (!is.null(metabolic_ledger_df) && nrow(metabolic_ledger_df)) {
      data.frame(
        ledger_text = metabolic_ledger_df$ledger_text,
        ledger_type = "Intracellular",
        ledger_priority = ifelse(
          metabolic_ledger_df$evidence_label == "high confidence", 2.0, 1.0
        ),
        ledger_color = "#214E3A",
        ledger_fontface = ifelse(
          metabolic_ledger_df$evidence_label == "high confidence", "bold", "plain"
        ),
        stringsAsFactors = FALSE
      )
    } else NULL,
    if (nrow(core_ledger_df)) {
      data.frame(
        ledger_text = core_ledger_df$ledger_text,
        ledger_type = "Core metabolism",
        ledger_priority = ifelse(core_ledger_df$status == "active", 2.75, 1.75),
        ledger_color = vapply(core_ledger_df$status, core_state_color, character(1)),
        ledger_fontface = ifelse(core_ledger_df$status == "active", "bold", "plain"),
        stringsAsFactors = FALSE
      )
    } else NULL
  )
  if (nrow(evidence_ledger_df)) {
    map_y_bottom <- y_bottom
    max_rows_per_block <- 24L
    n_blocks <- ceiling(nrow(evidence_ledger_df) / max_rows_per_block)
    rows_per_block <- ceiling(nrow(evidence_ledger_df) / n_blocks)
    evidence_ledger_df$block <- rep(
      seq_len(n_blocks), each = rows_per_block,
      length.out = nrow(evidence_ledger_df)
    )
    evidence_ledger_df$row_in_block <- ave(
      seq_len(nrow(evidence_ledger_df)),
      evidence_ledger_df$block,
      FUN = seq_along
    )
    ledger_x0 <- x_left + 0.30
    block_width <- (x_right - x_left - 0.60) / n_blocks
    wrap_chars <- max(60L, floor(block_width / 0.075))
    evidence_ledger_df$wrapped_text <- vapply(
      evidence_ledger_df$ledger_text,
      function(value) paste(strwrap(value, width = wrap_chars), collapse = "\n"),
      character(1)
    )
    evidence_ledger_df$n_lines <- pmax(
      1L,
      lengths(strsplit(evidence_ledger_df$wrapped_text, "\n", fixed = TRUE))
    )
    evidence_ledger_df$x <- ledger_x0 +
      (evidence_ledger_df$block - 1L) * block_width
    ledger_y_title <- map_y_bottom - 0.55
    ledger_y_top <- ledger_y_title - 0.72
    evidence_ledger_df$y <- NA_real_
    ledger_bottom <- ledger_y_top
    for (block_id in seq_len(n_blocks)) {
      idx <- which(evidence_ledger_df$block == block_id)
      row_height <- 0.25 * evidence_ledger_df$n_lines[idx] + 0.20
      row_bottom <- ledger_y_top - cumsum(row_height)
      evidence_ledger_df$y[idx] <- row_bottom + row_height / 2
      ledger_bottom <- min(ledger_bottom, min(row_bottom))
    }
    y_bottom <- ledger_bottom - 0.45
    p <- p + ggplot2::annotate(
      "text",
      x = ledger_x0, y = ledger_y_title,
      label = "Full functional evidence",
      hjust = 0, size = 2.5,
      fontface = "bold", color = "#243B53"
    ) + ggplot2::annotate(
      "text",
      x = ledger_x0, y = ledger_y_title - 0.28,
      label = paste0(
        "T-ID = membrane transporter locus; G-ID = extracellular GH locus; ",
        "intracellular and central-carbon genes are labeled beside reactions; ",
        "all retained evidence is listed below"
      ),
      hjust = 0, size = 1.55, color = "#607D8B"
    ) + ggplot2::geom_text(
      data = evidence_ledger_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$wrapped_text),
      hjust = 0, vjust = 0.5, size = 1.55,
      color = evidence_ledger_df$ledger_color,
      fontface = evidence_ledger_df$ledger_fontface,
      lineheight = 0.95, inherit.aes = FALSE
    )
  }
  x_span <- max(1, x_right - x_left)
  y_span <- max(1, y_top_plot - y_bottom)
  n_active <- length(cs_ids_ordered)
  # Match the device to the fixed-coordinate content so the one-page map does
  # not retain the large blank band created by the former square minimum.
  unit_in <- 0.82
  plot_height <- max(10, min(24, y_span * unit_in + 0.9))
  min_width <- max(8, 5 + n_active * 0.55)
  plot_width <- max(min_width, min(30, x_span * unit_in + 0.9))

  p <- p +
    ggplot2::coord_fixed(ratio = 1,
      xlim = c(x_left, x_right), ylim = c(y_bottom, y_top_plot),
      expand = FALSE, clip = "off") +
    ggplot2::labs(
      title = "CAZy Carbon Transport Map",
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
      plot.margin = ggplot2::margin(8, 8, 10, 8))

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(file_stub, ".pdf"))
  .dnmb_module_plot_save(
    p, pdf_path,
    width = plot_width, height = plot_height,
    device = grDevices::cairo_pdf
  )
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
    # → 2-deoxyribose-5-P → GA3P + acetaldehyde
    deoxyribose = "GA3P",
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
    deoxyribose = list(c("Deoxyribose-5-P", "deoxyribose", "deoK")),
    gluconate   = list(c("6-PG", "gluconate", "gntK")),
    glycerol    = list(c("Glycerol-3-P", "glycerol", "glpK")),
    fucose      = list(c("Fuculose", "fucose", "fucI"), c("Fuculose-1-P", "fucose", "fucK")),
    rhamnose    = list(c("Rhamnulose", "rhamnose", "rhaA"), c("Rhamnulose-1-P", "rhamnose", "rhaB"))
  )
  names(out) <- tolower(names(out))
  out
}

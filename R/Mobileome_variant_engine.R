# ---------------------------------------------------------------------------
# Built-in table of ~130 universally essential bacterial marker gene families
# Source categories: DEG (Database of Essential Genes), BUSCO bacterial markers,
# CheckM universal markers.
# Each entry: gene_name (canonical), category, product_keywords, marker_score
# ---------------------------------------------------------------------------
.dnmb_essential_marker_genes <- function() {
  dplyr::bind_rows(
    # --- Ribosomal proteins 30S (rpsA-rpsU) --------------------------------
    tibble::tibble(
      gene_name = c("rpsa", "rpsb", "rpsc", "rpsd", "rpse", "rpsf", "rpsg",
                     "rpsh", "rpsi", "rpsj", "rpsk", "rpsl", "rpsm", "rpsn",
                     "rpso", "rpsp", "rpsq", "rpsr", "rpss", "rpst", "rpsu"),
      category = "ribosomal_protein_30S",
      product_keywords = "30s ribosomal protein",
      marker_score = 0.93
    ),
    # --- Ribosomal proteins 50S (rplA-rplY) --------------------------------
    tibble::tibble(
      gene_name = c("rpla", "rplb", "rplc", "rpld", "rple", "rplf", "rplg",  # nolint
                     "rplh", "rpli", "rplj", "rplk", "rpll", "rplm", "rpln",
                     "rplo", "rplp", "rplq", "rplr", "rpls", "rplt", "rplu",
                     "rplv", "rplw", "rplx", "rply"),
      category = "ribosomal_protein_50S",
      product_keywords = "50s ribosomal protein",
      marker_score = 0.93
    ),
    # --- DNA replication ---------------------------------------------------
    tibble::tibble(
      gene_name = c("dnaa", "dnab", "dnac", "dnae", "dnag", "dnan", "dnaq",
                     "dnax", "gyra", "gyrb", "liga", "pola", "ssb"),
      category = "dna_replication",
      product_keywords = "dna polymerase|dna gyrase|dna ligase|dna helicase|replicative helicase|single-strand|chromosomal replication initiator",
      marker_score = 0.92
    ),
    # --- Transcription -----------------------------------------------------
    tibble::tibble(
      gene_name = c("rpoa", "rpob", "rpoc", "rpod", "nusa", "nusb", "nusg",
                     "rho"),
      category = "transcription",
      product_keywords = "rna polymerase|transcription termination|transcription antitermination",
      marker_score = 0.92
    ),
    # --- Translation factors -----------------------------------------------
    tibble::tibble(
      gene_name = c("infa", "infb", "infc", "fusa", "tsf", "tufa", "tufb",
                     "efp", "frr", "prfa", "prfb"),
      category = "translation_factor",
      product_keywords = "translation initiation factor|elongation factor|peptide chain release factor|ribosome recycling factor",
      marker_score = 0.85
    ),
    # --- tRNA synthetases --------------------------------------------------
    tibble::tibble(
      gene_name = c("alas", "args", "asns", "asps", "cyss", "glns", "glus",
                     "glys", "hiss", "iles", "leus", "lyss", "mets", "phes",
                     "pros", "sers", "thrs", "trps", "tyrs", "vals"),
      category = "trna_synthetase",
      product_keywords = "aminoacyl-trna synthetase|trna ligase|tRNA synthetase|--tRNA ligase",
      marker_score = 0.83
    ),
    # --- Cell division -----------------------------------------------------
    tibble::tibble(
      gene_name = c("ftsa", "ftsb", "ftsi", "ftsk", "ftsl", "ftsn", "ftsq",
                     "ftsw", "ftsz", "minc", "mind", "mine", "mreb",
                     "mura", "murb", "murc", "murd", "mure", "murf", "murg"),
      category = "cell_division",
      product_keywords = "cell division|septum|divisome|peptidoglycan|murein|shape-determining",
      marker_score = 0.73
    ),
    # --- Protein export / secretion ----------------------------------------
    tibble::tibble(
      gene_name = c("seca", "secd", "sece", "secf", "secg", "secy", "ffh",
                     "ftsy", "yidc", "lepb"),
      category = "protein_export",
      product_keywords = "preprotein translocase|signal recognition particle|sec translocase|leader peptidase|membrane protein insertase",
      marker_score = 0.73
    ),
    # --- Chaperones --------------------------------------------------------
    tibble::tibble(
      gene_name = c("groel", "groes", "dnak", "dnaj", "grpe"),
      category = "chaperone",
      product_keywords = "chaperonin|co-chaperonin|molecular chaperone|heat shock",
      marker_score = 0.63
    ),
    # --- Lipid / fatty acid biosynthesis -----------------------------------
    tibble::tibble(
      gene_name = c("acca", "accb", "accc", "accd", "fabd", "fabg", "fabh",
                     "fabi", "fabz", "plsb", "plsc"),
      category = "lipid_biosynthesis",
      product_keywords = "acetyl-coa carboxylase|acyl carrier|3-oxoacyl|enoyl|malonyl|glycerol-3-phosphate acyltransferase",
      marker_score = 0.73
    ),
    # --- Cell envelope (Gram-negative LPS) ---------------------------------
    tibble::tibble(
      gene_name = c("lpxa", "lpxb", "lpxc", "lpxd", "lpxh", "lpxk",
                     "kdsa", "kdsb"),
      category = "cell_envelope",
      product_keywords = "lipid a|kdo|3-deoxy-d-manno|udp-3-o-acyl",
      marker_score = 0.70
    ),
    # --- Energy: ATP synthase subunits -------------------------------------
    tibble::tibble(
      gene_name = c("atpa", "atpb", "atpc", "atpd", "atpe", "atpf", "atpg",
                     "atph"),
      category = "energy_atp_synthase",
      product_keywords = "atp synthase|f0f1 atp synthase|f1f0 atp synthase",
      marker_score = 0.63
    ),
    # --- Other essential ---------------------------------------------------
    tibble::tibble(
      gene_name = c("map", "def", "coaa", "coad", "coae",
                     "hema", "hemb", "hemc", "hemd", "heme"),
      category = "other_essential",
      product_keywords = "methionine aminopeptidase|peptide deformylase|pantothenate kinase|phosphopantetheine adenylyltransferase|dephospho-coa kinase|glutamyl-trna reductase|porphobilinogen|hydroxymethylbilane|uroporphyrinogen",
      marker_score = 0.65
    )
  )
}

# ---------------------------------------------------------------------------
# Gene essentiality prediction (two-tier system)
# Tier 1 (default ON): Completeness marker gene matching
# Tier 2 (default OFF): BPGAconverter core gene frequency integration
# ---------------------------------------------------------------------------
.dnmb_predict_gene_essentiality <- function(genes, core_gene_table = NULL) {
  if (!nrow(genes)) {
    return(genes %>%
      dplyr::mutate(
        redundancy_key = character(),
        redundancy_count = integer(),
        essentiality_score = numeric(),
        essentiality_class = character(),
        essentiality_reasons = character(),
        marker_gene_match = character(),
        essentiality_tier = character()
      ))
  }

  # --- Prepare annotation text and redundancy ------------------------------
  genes <- genes %>%
    dplyr::mutate(
      annotation_lc = stringr::str_to_lower(stringr::str_squish(paste(
        .data$feature_types,
        .data$gene,
        .data$product,
        .data$note
      ))),
      redundancy_key = dplyr::case_when(
        !is.na(.data$gene) & nzchar(.data$gene) ~ stringr::str_to_lower(.data$gene),
        !is.na(.data$product) & nzchar(.data$product) ~ .dnmb_normalize_product(.data$product),
        TRUE ~ .data$gene_id
      )
    )

  redundancy <- genes %>%
    dplyr::count(.data$redundancy_key, name = "redundancy_count")

  genes <- genes %>%
    dplyr::left_join(redundancy, by = "redundancy_key")

  # =========================================================================
  # TIER 1: Marker gene matching
  # =========================================================================
  markers <- .dnmb_essential_marker_genes()

  gene_name_lc <- stringr::str_to_lower(dplyr::coalesce(genes$gene, ""))
  product_lc   <- stringr::str_to_lower(dplyr::coalesce(genes$product, ""))

  # --- Exact gene name match (highest priority) ----------------------------
  marker_lookup <- stats::setNames(seq_len(nrow(markers)), markers$gene_name)
  name_match_idx <- marker_lookup[gene_name_lc]

  marker_gene_match <- rep(NA_character_, nrow(genes))
  marker_score      <- rep(NA_real_, nrow(genes))
  marker_category   <- rep(NA_character_, nrow(genes))
  marker_tier       <- rep(NA_character_, nrow(genes))

  has_name_match <- !is.na(name_match_idx)
  if (any(has_name_match)) {
    matched_rows <- name_match_idx[has_name_match]
    marker_gene_match[has_name_match] <- markers$gene_name[matched_rows]
    marker_score[has_name_match]      <- markers$marker_score[matched_rows]
    marker_category[has_name_match]   <- markers$category[matched_rows]
    marker_tier[has_name_match]       <- "tier1_name_match"
  }


  # --- Product keyword match (medium priority, for unmatched genes) --------
  unmatched <- is.na(marker_gene_match)
  if (any(unmatched)) {
    unique_categories <- unique(markers[, c("category", "product_keywords", "marker_score")])
    for (i in seq_len(nrow(unique_categories))) {
      kw_pattern <- unique_categories$product_keywords[i]
      cat_name   <- unique_categories$category[i]
      cat_score  <- unique_categories$marker_score[i]
      hits <- unmatched & grepl(kw_pattern, product_lc, perl = TRUE)
      if (any(hits)) {
        # Product match gets a small penalty vs exact name match
        marker_gene_match[hits] <- paste0(cat_name, "_product")
        marker_score[hits]      <- cat_score - 0.05
        marker_category[hits]   <- cat_name
        marker_tier[hits]       <- "tier1_product_match"
        unmatched[hits]         <- FALSE
      }
    }
  }

  # =========================================================================
  # TIER 2 (stub): BPGAconverter core gene frequency
  # When core_gene_table is provided, override scores for matched locus_tags
  # =========================================================================
  if (!is.null(core_gene_table)) {
    # Expected columns: locus_tag, pangenome_class (core/soft_core/shell/cloud)
    if (all(c("locus_tag", "pangenome_class") %in% colnames(core_gene_table))) {
      tier2_scores <- dplyr::case_when(
        core_gene_table$pangenome_class == "core"      ~ 0.95,
        core_gene_table$pangenome_class == "soft_core" ~ 0.75,
        core_gene_table$pangenome_class == "shell"     ~ 0.45,
        core_gene_table$pangenome_class == "cloud"     ~ 0.15,
        TRUE ~ NA_real_
      )
      tier2_lookup <- stats::setNames(tier2_scores, core_gene_table$locus_tag)
      tier2_match  <- tier2_lookup[genes$locus_tag]
      has_tier2 <- !is.na(tier2_match)
      if (any(has_tier2)) {
        # Tier 2 takes precedence when available
        marker_score[has_tier2]    <- tier2_match[has_tier2]
        marker_tier[has_tier2]     <- "tier2_core_gene"
        marker_category[has_tier2] <- paste0("pangenome_", core_gene_table$pangenome_class[
          match(genes$locus_tag[has_tier2], core_gene_table$locus_tag)
        ])
      }
    }
  }

  # =========================================================================
  # FALLBACK: Original keyword heuristics for genes not matched by Tier 1/2
  # =========================================================================
  text <- genes$annotation_lc
  feature_types <- dplyr::coalesce(genes$feature_types, "")
  is_plasmid <- grepl("plasmid", stringr::str_to_lower(dplyr::coalesce(genes$definition, "")), fixed = TRUE)

  rrna <- grepl("rRNA", feature_types, fixed = TRUE)
  structural_rna <- grepl("tRNA|tmRNA", feature_types)
  translation_core <- grepl("ribosomal protein|ribosome maturation|translation initiation factor|translation elongation factor|release factor|aminoacyl-tRNA synthetase|peptide chain release factor|transfer-messenger rna", text)
  transcription_core <- grepl("rna polymerase|\\brpo[abczn]\\b|transcription termination factor rho|nus[abg]", text)
  replication_core <- grepl("\\bdnaa\\b|\\bdnab\\b|\\bdnag\\b|\\bdnae\\b|\\bdnan\\b|dna polymerase iii|dna ligase|\\bliga\\b|\\bgyr[ab]\\b|\\bpar[ce]\\b|topoisomerase iv|replicative helicase|chromosome partition", text)
  division_wall <- grepl("\\bfts[azqiwklnh]\\b|cell division|septum|divisome|\\bmreb\\b|murein|peptidoglycan|\\bmur[abcdefgijmnp]\\b|shape-determining protein", text)
  export_core <- grepl("\\bsec[adfyg]\\b|signal recognition particle|\\bffh\\b|\\bftsy\\b|leader peptidase|translocase", text)
  homeostasis_core <- grepl("\\bgroel\\b|\\bdnak\\b|\\bgroes\\b|chaperon", text)
  energy_core <- grepl("\\batp synthase\\b", text)

  mobile_accessory <- grepl("transposase|insertion sequence|integrase|recombinase|resolvase|site-specific recombinase|mobile element", text)
  phage_accessory <- grepl("phage|capsid|tail protein|prophage|portal protein|terminase", text)
  plasmid_transfer <- grepl("plasmid|conjug|tra[[:digit:]]|type iv secretion|relaxase", text)
  condition_defense <- grepl("resistan|efflux|beta-lactamase|toxin|antitoxin|restriction endonuclease|crispr|\\bcas[[:digit:]]", text)
  niche_adaptation <- grepl("flagell|chemotaxis|biofilm|sporulation|adhesin|pilus|fimbr", text)
  annotation_uncertain <- grepl("hypothetical protein|uncharacterized protein|domain of unknown function|duf", text)
  redundancy_detected <- genes$redundancy_count > 1

  fallback_score <- 0.05 +
    0.95 * rrna +
    0.9 * structural_rna +
    0.65 * translation_core +
    0.65 * transcription_core +
    0.7 * replication_core +
    0.55 * division_wall +
    0.45 * export_core +
    0.2 * homeostasis_core +
    0.2 * energy_core -
    0.8 * mobile_accessory -
    0.7 * phage_accessory -
    0.55 * plasmid_transfer -
    0.45 * condition_defense -
    0.35 * niche_adaptation -
    0.25 * annotation_uncertain -
    0.15 * is_plasmid -
    0.15 * redundancy_detected

  fallback_score <- pmax(0.01, pmin(0.99, fallback_score))

  # =========================================================================
  # Combine: marker score takes priority; fallback for unmatched
  # =========================================================================
  has_marker <- !is.na(marker_score)

  # For marker-matched genes, apply negative modifiers from mobile/phage etc.
  negative_mod <- -(0.8 * mobile_accessory +
                    0.7 * phage_accessory +
                    0.55 * plasmid_transfer +
                    0.15 * is_plasmid +
                    0.15 * redundancy_detected)

  score <- dplyr::if_else(
    has_marker,
    pmax(0.01, pmin(0.99, marker_score + negative_mod)),
    fallback_score
  )

  # Fill tier for fallback genes
  marker_tier <- dplyr::if_else(is.na(marker_tier), "keyword_heuristic", marker_tier)
  marker_gene_match <- dplyr::if_else(is.na(marker_gene_match), NA_character_, marker_gene_match)

  essentiality_class <- dplyr::case_when(
    score >= 0.8 ~ "high",
    score >= 0.45 ~ "medium",
    TRUE ~ "low"
  )

  # --- Build reason strings ------------------------------------------------
  reason_matrix <- data.frame(
    rRNA = ifelse(rrna, "rRNA", NA_character_),
    structural_RNA = ifelse(structural_rna, "structural_RNA", NA_character_),
    translation_core = ifelse(translation_core, "translation_core", NA_character_),
    transcription_core = ifelse(transcription_core, "transcription_core", NA_character_),
    replication_chromosome = ifelse(replication_core, "replication_chromosome", NA_character_),
    cell_division_wall = ifelse(division_wall, "cell_division_wall", NA_character_),
    protein_export_core = ifelse(export_core, "protein_export_core", NA_character_),
    protein_homeostasis = ifelse(homeostasis_core, "protein_homeostasis", NA_character_),
    energy_core = ifelse(energy_core, "energy_core", NA_character_),
    mobile_accessory = ifelse(mobile_accessory, "mobile_accessory", NA_character_),
    phage_accessory = ifelse(phage_accessory, "phage_accessory", NA_character_),
    plasmid_transfer = ifelse(plasmid_transfer, "plasmid_transfer", NA_character_),
    condition_specific_defense = ifelse(condition_defense, "condition_specific_defense", NA_character_),
    niche_adaptation = ifelse(niche_adaptation, "niche_adaptation", NA_character_),
    annotation_uncertain = ifelse(annotation_uncertain, "annotation_uncertain", NA_character_),
    plasmid_location = ifelse(is_plasmid, "plasmid_location", NA_character_),
    redundancy_detected = ifelse(redundancy_detected, "redundancy_detected", NA_character_),
    marker_gene = ifelse(has_marker, paste0("marker:", marker_gene_match), NA_character_)
  )
  essentiality_reasons <- vapply(
    seq_len(nrow(reason_matrix)),
    function(i) {
      vals <- stats::na.omit(as.character(reason_matrix[i, ]))
      if (!length(vals)) NA_character_ else paste(vals, collapse = "; ")
    },
    character(1)
  )
  phenotype_strings <- .dnmb_predict_phenotype_strings_vector(genes$annotation_text)

  genes %>%
    dplyr::mutate(
      essentiality_score = score,
      essentiality_class = essentiality_class,
      essentiality_reasons = essentiality_reasons,
      marker_gene_match = marker_gene_match,
      essentiality_tier = marker_tier,
      gene_phenotype_categories = phenotype_strings$categories,
      gene_phenotype_prediction = phenotype_strings$predictions
    )
}

.dnmb_normalize_product <- function(product) {
  text <- stringr::str_to_lower(dplyr::coalesce(product, ""))
  text <- stringr::str_replace_all(text, "[^a-z0-9]+", " ")
  text <- stringr::str_remove_all(text, "\\b(protein|subunit|putative|probable|family|type|chain|domain)\\b")
  stringr::str_squish(text)
}

.dnmb_family_variant_rules <- function() {
  tibble::tribble(
    ~family, ~fallback_model_type, ~fallback_motif, ~transposition_modes, ~dr_expected_bp, ~tir_expected, ~rule_source, ~recognition_mode, ~supplementary_strategy,
    "IS1", "at_rich", "AT", "copy_and_paste; cointegrate", 8L, TRUE, "ISfinder_general_features", "dna_tsd", "TSD_plus_TIR_plus_local_context",
    "IS3", "motif", "TG", "copy_paste", 4L, TRUE, "ISfinder_general_features", "dna_tsd", "TSD_plus_TIR_plus_local_context",
    "IS4", "motif", "CT", "cut_and_paste", 9L, TRUE, "ISfinder_general_features", "dna_tsd", "short_TSD_plus_context_enrichment",
    "IS30", "none", NA_character_, "copy_and_paste", 2L, TRUE, "ISfinder_general_features", "dna_context", "short_TSD_plus_TIR_architecture",
    "IS66", "motif", "GTAA", "unknown_dde", 8L, TRUE, "ISfinder_general_features", "dna_tsd", "TSD_set_plus_context_enrichment",
    "IS91", "motif_set", "GAAC; CAAG", "rolling_circle", 0L, FALSE, "ISfinder_general_features", "rolling_circle_context", "ori_like_end_signal_plus_context",
    "IS110", "none", NA_character_, "unknown_dedd", 0L, FALSE, "ISfinder_general_features", "rna_guided_candidate", "bridgeRNA_candidate_plus_context",
    "IS200/IS605", "none", NA_character_, "peel_and_paste", 0L, FALSE, "ISfinder_general_features", "hairpin_or_rna_assisted", "hairpin_or_omegaRNA_candidate_plus_context",
    "IS256", "none", NA_character_, "copy_paste", 8L, TRUE, "ISfinder_general_features", "dna_context", "TIR_plus_context_enrichment",
    "IS630", "motif", "TA", "cut_and_paste", 2L, TRUE, "ISfinder_general_features", "dna_context", "short_TSD_plus_context_enrichment",
    "IS481", "motif", "TGT", "copy_paste_possible", 4L, TRUE, "ISfinder_general_features", "dna_context", "short_TSD_plus_context_enrichment",
    "IS701", "none", NA_character_, "unknown_dde", 4L, TRUE, "ISfinder_table1", "dna_tsd", "TSD_plus_context_enrichment",
    "ISH3", "motif", "CGT", "unknown_dde", 4L, TRUE, "ISfinder_table1", "dna_context", "short_TSD_plus_context_enrichment",
    "IS1634", "motif", "C", "unknown_dde", 5L, TRUE, "ISfinder_table1", "dna_context", "very_short_TSD_plus_context_only",
    "IS5", "motif", "GA", "unknown_dde", 4L, TRUE, "ISfinder_table1", "dna_context", "short_TSD_plus_context_enrichment",
    "new", "none", NA_character_, "unknown", NA_integer_, NA, "none", "unresolved", "cluster_specific_followup"
  )
}

.dnmb_build_background_kmer_index <- function(metadata, genes = NULL) {
  if (is.null(metadata) || !nrow(metadata) || !"sequence" %in% colnames(metadata)) {
    return(list(sequences = character(), total_bp = 0L))
  }

  seqs <- toupper(dplyr::coalesce(metadata$sequence, ""))
  seqs <- seqs[nzchar(seqs)]
  list(
    sequences = seqs,
    total_bp = sum(nchar(seqs), na.rm = TRUE)
  )
}

.dnmb_count_exact_motif_occurrences <- function(sequences, motif) {
  if (!length(sequences) || is.na(motif) || !nzchar(motif)) {
    return(0L)
  }
  motif <- toupper(motif)
  if (!grepl("^[ACGT]+$", motif)) {
    return(0L)
  }

  sum(vapply(sequences, function(seq) {
    hits <- gregexpr(motif, seq, fixed = TRUE)[[1]]
    if (length(hits) == 1L && identical(hits[[1]], -1L)) {
      0L
    } else {
      length(hits)
    }
  }, integer(1)))
}

.dnmb_score_tsd_enrichment <- function(tsd_rows, background_index) {
  if (is.null(tsd_rows) || !nrow(tsd_rows)) {
    return(tibble::tibble(
      tsd_seq = character(),
      tsd_len_bp = integer(),
      obs_count = integer(),
      obs_prop = numeric(),
      bg_count = integer(),
      bg_prop = numeric(),
      log2_enrichment = numeric()
    ))
  }

  observed <- tsd_rows %>%
    dplyr::count(.data$tsd_seq, .data$tsd_len_bp, name = "obs_count") %>%
    dplyr::group_by(.data$tsd_len_bp) %>%
    dplyr::mutate(obs_prop = .data$obs_count / sum(.data$obs_count)) %>%
    dplyr::ungroup()

  seqs <- background_index$sequences
  if (is.null(seqs)) {
    seqs <- character()
  }
  if (!length(seqs)) {
    return(observed %>%
      dplyr::mutate(
        bg_count = 0L,
        bg_prop = 0,
        log2_enrichment = NA_real_
      ) %>%
      dplyr::arrange(dplyr::desc(.data$obs_count), .data$tsd_len_bp, .data$tsd_seq))
  }

  possible_sites <- observed %>%
    dplyr::distinct(.data$tsd_len_bp) %>%
    dplyr::mutate(
      total_positions = vapply(
        .data$tsd_len_bp,
        function(k) sum(pmax(0L, nchar(seqs) - as.integer(k) + 1L)),
        numeric(1)
      )
    )

  observed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      bg_count = .dnmb_count_exact_motif_occurrences(seqs, .data$tsd_seq)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(possible_sites, by = "tsd_len_bp") %>%
    dplyr::mutate(
      bg_prop = dplyr::if_else(.data$total_positions > 0, .data$bg_count / .data$total_positions, 0),
      log2_enrichment = log2((.data$obs_prop + 1e-6) / (.data$bg_prop + 1e-6))
    ) %>%
    dplyr::select(
      .data$tsd_seq,
      .data$tsd_len_bp,
      .data$obs_count,
      .data$obs_prop,
      .data$bg_count,
      .data$bg_prop,
      .data$log2_enrichment
    ) %>%
    dplyr::arrange(dplyr::desc(.data$log2_enrichment), dplyr::desc(.data$obs_count), .data$tsd_len_bp, .data$tsd_seq)
}

.dnmb_format_enriched_examples <- function(enrichment_tbl, top_n = 5L) {
  if (is.null(enrichment_tbl) || !nrow(enrichment_tbl)) {
    return(NA_character_)
  }

  top_tbl <- enrichment_tbl %>%
    dplyr::filter(!is.na(.data$log2_enrichment)) %>%
    dplyr::arrange(dplyr::desc(.data$log2_enrichment), dplyr::desc(.data$obs_count), .data$tsd_seq) %>%
    dplyr::slice_head(n = as.integer(top_n))

  if (!nrow(top_tbl)) {
    return(NA_character_)
  }

  paste0(
    top_tbl$tsd_seq,
    "(",
    top_tbl$obs_count,
    ",log2=",
    sprintf("%.2f", top_tbl$log2_enrichment),
    ")",
    collapse = "; "
  )
}

.dnmb_build_target_models <- function(elements, metadata, genes) {
  if (!nrow(elements)) {
    return(tibble::tibble())
  }

  rules <- .dnmb_family_variant_rules()
  seq_map <- stats::setNames(metadata$sequence, metadata$contig)
  background_index <- .dnmb_build_background_kmer_index(metadata = metadata, genes = genes)

  element_context <- lapply(seq_len(nrow(elements)), function(i) {
    context <- .dnmb_classify_interval_context(
      contig = elements$contig[[i]],
      start = elements$start[[i]],
      end = elements$end[[i]],
      genes = genes
    )
    tibble::tibble(
      element_id = elements$element_id[[i]],
      context_class = context$context_class
    )
  }) %>% dplyr::bind_rows()

  elements2 <- elements %>%
    dplyr::left_join(element_context, by = "element_id")

  families <- unique(elements2$element_family)
  families <- families[!is.na(families) & nzchar(families)]

  models <- lapply(seq_along(families), function(i) {
    family_name <- families[[i]]
    family_elements <- elements2 %>% dplyr::filter(.data$element_family == .env$family_name)
    rule <- rules %>% dplyr::filter(.data$family == .env$family_name)
    if (!nrow(rule)) {
      rule <- rules %>% dplyr::filter(.data$family == "new")
    }

    tsd_rows <- .dnmb_collect_empirical_tsd_signals(
      family_elements = family_elements,
      seq_map = seq_map,
      expected_len = rule$dr_expected_bp[[1]]
    )
    enrichment_tbl <- .dnmb_score_tsd_enrichment(tsd_rows = tsd_rows, background_index = background_index)

    motif_source <- "none"
    model_type <- "none"
    model_motif <- NA_character_
    motif_specificity <- NA_real_
    dominant_tsd_len <- NA_integer_
    observed_tsd_count <- nrow(tsd_rows)
    enriched_tsd_examples <- .dnmb_format_enriched_examples(enrichment_tbl)
    top_enrichment_log2 <- if (nrow(enrichment_tbl)) max(enrichment_tbl$log2_enrichment, na.rm = TRUE) else NA_real_

    if (nrow(tsd_rows)) {
      dominant_tsd_len <- as.integer(names(which.max(table(tsd_rows$tsd_len_bp)))[1])
      dominant_seq_values <- tsd_rows$tsd_seq[tsd_rows$tsd_len_bp == dominant_tsd_len]
      dominant_seqs <- unique(dominant_seq_values)
      if (length(dominant_seqs)) {
        seq_freq <- sort(table(dominant_seq_values), decreasing = TRUE)
        consensus_motif <- .dnmb_iupac_consensus(dominant_seqs)
        consensus_specificity <- .dnmb_sequence_specificity(dominant_seqs)
        motif_source <- "observed_genome"
        if (length(dominant_seqs) == 1L || (consensus_specificity >= 0.8 && length(dominant_seqs) <= 3L)) {
          model_motif <- consensus_motif
          motif_specificity <- consensus_specificity
          model_type <- "empirical_tsd"
        } else if (sum(head(seq_freq, min(3L, length(seq_freq)))) / sum(seq_freq) >= 0.6) {
          model_motif <- paste(names(head(seq_freq, min(3L, length(seq_freq)))), collapse = "; ")
          motif_specificity <- as.numeric(max(seq_freq) / sum(seq_freq))
          model_type <- "motif_set"
        } else {
          model_type <- "context_only"
          model_motif <- NA_character_
          motif_specificity <- as.numeric(max(seq_freq) / sum(seq_freq))
        }
      }
    }

    empirical_motif_len <- if (!is.na(model_motif) && nzchar(model_motif)) nchar(gsub("[^A-Z]", "", model_motif)) else 0L
    if (
      identical(model_type, "empirical_tsd") &&
      (
        empirical_motif_len <= 2L ||
        (empirical_motif_len <= 3L && (observed_tsd_count < 5L || dplyr::coalesce(motif_specificity, 0) < 0.75)) ||
        dplyr::coalesce(motif_specificity, 0) < 0.55
      )
    ) {
      model_type <- "context_only"
      model_motif <- NA_character_
    }

    if (
      model_type == "context_only" &&
      nrow(enrichment_tbl)
    ) {
      top_exact <- enrichment_tbl %>%
        dplyr::filter(.data$tsd_len_bp == dominant_tsd_len) %>%
        dplyr::arrange(dplyr::desc(.data$log2_enrichment), dplyr::desc(.data$obs_count), .data$tsd_seq)
      if (nrow(top_exact)) {
        keep_exact <- top_exact %>%
          dplyr::filter(
            (.data$tsd_len_bp <= 2 & .data$obs_count >= 3 & .data$log2_enrichment >= 2.5) |
              (.data$tsd_len_bp == 3 & .data$obs_count >= 3 & .data$log2_enrichment >= 2.0) |
              (.data$tsd_len_bp >= 4 & .data$obs_count >= 2 & .data$log2_enrichment >= 1.2)
          ) %>%
          dplyr::slice_head(n = 3)
        if (nrow(keep_exact)) {
          model_type <- "motif_set"
          model_motif <- paste(keep_exact$tsd_seq, collapse = "; ")
          motif_specificity <- dplyr::coalesce(max(keep_exact$obs_prop, na.rm = TRUE), motif_specificity)
          motif_source <- "observed_genome_enriched"
        }
      }
    }

    if (identical(model_type, "none") && nrow(rule) && !is.na(rule$fallback_model_type[[1]]) && rule$fallback_model_type[[1]] != "none") {
      model_type <- rule$fallback_model_type[[1]]
      model_motif <- rule$fallback_motif[[1]]
      motif_specificity <- 0.35
      motif_source <- "fallback_rule"
    }

    if (
      identical(motif_source, "fallback_rule") &&
      !is.na(model_motif) &&
      nchar(gsub("[^A-Z]", "", model_motif)) < 4L
    ) {
      model_type <- "none"
      model_motif <- NA_character_
    }

    context_props <- family_elements %>%
      dplyr::count(.data$context_class, name = "n") %>%
      dplyr::mutate(prop = .data$n / sum(.data$n))
    intergenic_pref <- context_props$prop[match("intergenic", context_props$context_class)]
    promoter_pref <- context_props$prop[match("promoter", context_props$context_class)]
    gene_pref <- context_props$prop[match("gene_body", context_props$context_class)]
    intergenic_pref <- ifelse(length(intergenic_pref) && !is.na(intergenic_pref), intergenic_pref, 0)
    promoter_pref <- ifelse(length(promoter_pref) && !is.na(promoter_pref), promoter_pref, 0)
    gene_pref <- ifelse(length(gene_pref) && !is.na(gene_pref), gene_pref, 0)

    tibble::tibble(
      target_model_id = paste0("TM_", sprintf("%03d", i)),
      family = family_name,
      n_elements = nrow(family_elements),
      observed_tsd_count = observed_tsd_count,
      dominant_tsd_len = dominant_tsd_len,
      observed_tsd_examples = .dnmb_top_sequence_examples(tsd_rows, top_n = 5L),
      enriched_tsd_examples = enriched_tsd_examples,
      top_enrichment_log2 = top_enrichment_log2,
      model_type = model_type,
      model_motif = model_motif,
      motif_specificity = motif_specificity,
      motif_source = motif_source,
      recognition_strategy = dplyr::case_when(
        motif_source == "observed_genome" & model_type == "empirical_tsd" ~ "exact_or_consensus_tsd_from_input_genome",
        motif_source == "observed_genome" & model_type == "motif_set" ~ "top_exact_tsd_set_from_input_genome",
        motif_source == "observed_genome_enriched" & model_type == "motif_set" ~ "enriched_exact_tsd_set_from_input_genome",
        model_type == "context_only" ~ "context_only_from_input_genome_copy_distribution",
        motif_source == "fallback_rule" ~ "weak_family_rule_fallback",
        TRUE ~ "no_reliable_recognition_model"
      ),
      recognition_mode = rule$recognition_mode[[1]],
      supplementary_strategy = rule$supplementary_strategy[[1]],
      transposition_modes = rule$transposition_modes[[1]],
      dr_expected_bp = rule$dr_expected_bp[[1]],
      tir_expected = rule$tir_expected[[1]],
      preferred_intergenic = intergenic_pref,
      preferred_promoter = promoter_pref,
      preferred_gene_body = gene_pref,
      model_confidence = dplyr::case_when(
        motif_source == "observed_genome" && observed_tsd_count >= 3 && motif_specificity >= 0.7 ~ "high",
        motif_source == "observed_genome_enriched" && !is.na(top_enrichment_log2) && top_enrichment_log2 >= 2 ~ "medium",
        motif_source == "observed_genome" && observed_tsd_count >= 1 ~ "medium",
        motif_source == "fallback_rule" ~ "low",
        TRUE ~ "low"
      ),
      rule_source = rule$rule_source[[1]]
    )
  })

  dplyr::bind_rows(models)
}

.dnmb_top_sequence_examples <- function(tsd_rows, top_n = 5L) {
  if (is.null(tsd_rows) || !nrow(tsd_rows)) {
    return(NA_character_)
  }
  seq_counts <- tsd_rows %>%
    dplyr::count(.data$tsd_seq, sort = TRUE, name = "n") %>%
    dplyr::slice_head(n = as.integer(top_n))
  paste0(seq_counts$tsd_seq, "(", seq_counts$n, ")", collapse = "; ")
}

.dnmb_collect_empirical_tsd_signals <- function(family_elements, seq_map, expected_len = NA_integer_) {
  if (!nrow(family_elements)) {
    return(tibble::tibble(tsd_seq = character(), tsd_len_bp = integer(), source = character()))
  }

  rows <- list()
  for (i in seq_len(nrow(family_elements))) {
    seq <- seq_map[[family_elements$contig[[i]]]]
    if (is.null(seq) || !nzchar(seq)) {
      next
    }

    candidate_lens <- unique(c(
      if (!is.na(expected_len) && expected_len > 0) expected_len else integer(),
      2:12
    ))
    candidate_lens <- candidate_lens[candidate_lens > 0]
    if (!length(candidate_lens)) {
      next
    }

    for (k in candidate_lens) {
      if (family_elements$start[[i]] <= k || (family_elements$end[[i]] + k) > nchar(seq)) {
        next
      }
      left_seq <- substring(seq, family_elements$start[[i]] - k, family_elements$start[[i]] - 1L)
      right_seq <- substring(seq, family_elements$end[[i]] + 1L, family_elements$end[[i]] + k)
      if (identical(left_seq, right_seq)) {
        rows[[length(rows) + 1L]] <- tibble::tibble(
          tsd_seq = left_seq,
          tsd_len_bp = k,
          source = "boundary_duplication"
        )
      }
    }

    if (isTRUE(family_elements$tsd_found[[i]]) && !is.na(family_elements$tsd_seq[[i]]) && nzchar(family_elements$tsd_seq[[i]])) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        tsd_seq = family_elements$tsd_seq[[i]],
        tsd_len_bp = family_elements$tsd_len_bp[[i]],
        source = "sequence_engine"
      )
    }
  }

  if (!length(rows)) {
    return(tibble::tibble(tsd_seq = character(), tsd_len_bp = integer(), source = character()))
  }

  dplyr::bind_rows(rows) %>%
    dplyr::filter(!is.na(.data$tsd_seq), nzchar(.data$tsd_seq))
}

.dnmb_iupac_consensus <- function(seqs) {
  seqs <- unique(toupper(seqs[!is.na(seqs) & nzchar(seqs)]))
  if (!length(seqs)) {
    return(NA_character_)
  }
  n <- unique(nchar(seqs))
  if (length(n) != 1L) {
    seqs <- seqs[nchar(seqs) == min(n)]
    n <- min(nchar(seqs))
  }
  nuc_map <- c(
    "A" = "A", "C" = "C", "G" = "G", "T" = "T",
    "AG" = "R", "CT" = "Y", "CG" = "S", "AT" = "W", "GT" = "K", "AC" = "M",
    "CGT" = "B", "AGT" = "D", "ACT" = "H", "ACG" = "V", "ACGT" = "N"
  )
  chars <- strsplit(seqs, "", fixed = TRUE)
  consensus <- character(n)
  for (i in seq_len(n)) {
    bases <- sort(unique(vapply(chars, `[[`, character(1), i)))
    key <- paste(bases, collapse = "")
    consensus[[i]] <- nuc_map[[key]]
    if (is.null(consensus[[i]])) {
      consensus[[i]] <- "N"
    }
  }
  paste(consensus, collapse = "")
}

.dnmb_sequence_specificity <- function(seqs) {
  seqs <- unique(toupper(seqs[!is.na(seqs) & nzchar(seqs)]))
  if (!length(seqs)) {
    return(NA_real_)
  }
  n <- unique(nchar(seqs))
  if (length(n) != 1L) {
    seqs <- seqs[nchar(seqs) == min(n)]
    n <- min(nchar(seqs))
  }
  chars <- strsplit(seqs, "", fixed = TRUE)
  pos_scores <- vapply(seq_len(n), function(i) {
    bases <- vapply(chars, `[[`, character(1), i)
    max(table(bases)) / length(bases)
  }, numeric(1))
  mean(pos_scores)
}

.dnmb_predict_target_sites <- function(target_models, metadata, genes, max_target_sites_per_family = 2000L, site_scan_strategy = c("balanced", "full")) {
  if (!nrow(target_models)) {
    return(tibble::tibble())
  }
  max_target_sites_per_family <- as.integer(max_target_sites_per_family)
  site_scan_strategy <- match.arg(site_scan_strategy)
  if (site_scan_strategy == "balanced") {
    target_models <- target_models %>%
      dplyr::mutate(
        active_model_type = dplyr::case_when(
          .data$model_type == "context_only" & .data$model_confidence != "low" & .data$n_elements >= 3 ~ "context_only",
          .data$model_type %in% c("empirical_tsd", "motif_set") &
            .data$motif_source == "observed_genome" &
            .data$model_confidence == "high" &
            dplyr::coalesce(.data$dominant_tsd_len, 0L) >= 4L ~ .data$model_type,
          .data$model_type %in% c("empirical_tsd", "motif_set") &
            .data$motif_source == "observed_genome" &
            .data$model_confidence == "medium" &
            dplyr::coalesce(.data$dominant_tsd_len, 0L) >= 6L &
            .data$observed_tsd_count >= 5L ~ .data$model_type,
          TRUE ~ "none"
        )
      ) %>%
      dplyr::filter(.data$active_model_type != "none")
  } else {
    target_models <- target_models %>%
      dplyr::mutate(active_model_type = dplyr::if_else(.data$model_type == "none", "none", .data$model_type)) %>%
      dplyr::filter(.data$active_model_type != "none")
  }
  if (!nrow(target_models)) {
    return(tibble::tibble())
  }
  target_models <- target_models %>%
    dplyr::mutate(
      model_type = .data$active_model_type
    ) %>%
    dplyr::select(
      -active_model_type
    )
  seq_map <- stats::setNames(metadata$sequence, metadata$contig)
  context_cache <- .dnmb_prepare_context_cache(metadata = metadata, genes = genes)

  rows <- list()
  for (i in seq_len(nrow(target_models))) {
    model <- target_models[i, , drop = FALSE]
    family_sites <- list()
    for (contig in metadata$contig) {
      seq <- seq_map[[contig]]
      if (is.null(seq) || !nzchar(seq)) {
        next
      }
      if (site_scan_strategy == "balanced") {
        hits <- .dnmb_score_context_candidates_one_contig(
          seq = seq,
          contig = contig,
          model = model,
          context_sites = context_cache$context_sites_by_contig[[contig]]
        )
      } else {
        hits <- .dnmb_scan_target_sites_one_contig(
          seq = seq,
          contig = contig,
          model = model,
          genes_contig = context_cache$genes_by_contig[[contig]],
          context_sites = context_cache$context_sites_by_contig[[contig]],
          permissive_regions = context_cache$permissive_regions_by_contig[[contig]]
        )
      }
      if (nrow(hits)) {
        family_sites[[length(family_sites) + 1L]] <- hits
      }
    }
    family_sites <- dplyr::bind_rows(family_sites)
    if (!nrow(family_sites)) {
      next
    }
    family_sites <- family_sites %>%
      dplyr::arrange(dplyr::desc(.data$target_site_score), dplyr::desc(.data$motif_match_score), .data$contig, .data$site_start) %>%
      dplyr::slice_head(n = max_target_sites_per_family)
    rows[[length(rows) + 1L]] <- family_sites
  }
  dplyr::bind_rows(rows)
}

.dnmb_empty_targetable_regions <- function() {
  tibble::tibble(
    target_model_id = character(),
    target_family = character(),
    contig = character(),
    region_start = integer(),
    region_end = integer(),
    region_size_bp = integer(),
    region_type = character(),
    recognition_strategy = character(),
    model_type = character(),
    model_motif = character(),
    motif_source = character(),
    transposition_modes = character(),
    target_site_count = integer(),
    context_classes = character(),
    affected_genes = character(),
    affected_gene_count = integer(),
    highest_essentiality_class = character(),
    max_target_site_score = numeric(),
    mean_target_site_score = numeric(),
    max_essentiality_penalty = numeric(),
    representative_motifs = character()
  )
}

.dnmb_build_targetable_regions <- function(target_models, target_sites, metadata, genes, site_scan_strategy = c("balanced", "full"), merge_gap_bp = 50L) {
  site_scan_strategy <- match.arg(site_scan_strategy)
  if (!nrow(target_sites)) {
    return(.dnmb_empty_targetable_regions())
  }

  model_tbl <- target_models %>%
    dplyr::transmute(
      target_model_id = .data$target_model_id,
      target_family = .data$family,
      recognition_strategy = .data$recognition_strategy,
      model_type = .data$model_type,
      model_motif = .data$model_motif,
      motif_source = .data$motif_source,
      transposition_modes = .data$transposition_modes
    )

  ts <- target_sites %>%
    dplyr::left_join(model_tbl, by = c("target_model_id", "target_family")) %>%
    dplyr::arrange(.data$target_family, .data$contig, .data$site_start, .data$site_end)

  region_rows <- list()
  region_idx <- 1L
  family_contig_groups <- split(ts, interaction(ts$target_family, ts$contig, drop = TRUE))

  for (group in family_contig_groups) {
    if (!nrow(group)) {
      next
    }
    cur_start <- group$site_start[[1]]
    cur_end <- group$site_end[[1]]
    cur_rows <- 1L

    emit_region <- function(idx_vec) {
      block <- group[idx_vec, , drop = FALSE]
      max_target_site_score <- if (all(is.na(block$target_site_score))) NA_real_ else max(block$target_site_score, na.rm = TRUE)
      mean_target_site_score <- if (all(is.na(block$target_site_score))) NA_real_ else mean(block$target_site_score, na.rm = TRUE)
      max_essentiality_penalty <- if (all(is.na(block$essentiality_penalty))) NA_real_ else max(block$essentiality_penalty, na.rm = TRUE)
      highest_essentiality_class <- if (any(block$affected_gene_essentiality_class == "high", na.rm = TRUE)) {
        "high"
      } else if (any(block$affected_gene_essentiality_class == "medium", na.rm = TRUE)) {
        "medium"
      } else if (any(block$affected_gene_essentiality_class == "low", na.rm = TRUE)) {
        "low"
      } else {
        NA_character_
      }
      tibble::tibble(
        target_model_id = block$target_model_id[[1]],
        target_family = block$target_family[[1]],
        contig = block$contig[[1]],
        region_start = min(block$site_start, na.rm = TRUE),
        region_end = max(block$site_end, na.rm = TRUE),
        region_size_bp = max(block$site_end, na.rm = TRUE) - min(block$site_start, na.rm = TRUE) + 1L,
        region_type = dplyr::case_when(
          block$model_type[[1]] == "context_only" ~ "context_targetable_region",
          TRUE ~ "motif_supported_targetable_region"
        ),
        recognition_strategy = block$recognition_strategy[[1]],
        model_type = block$model_type[[1]],
        model_motif = block$model_motif[[1]],
        motif_source = block$motif_source[[1]],
        transposition_modes = block$transposition_modes[[1]],
        target_site_count = nrow(block),
        context_classes = paste(unique(stats::na.omit(block$context_class)), collapse = "; "),
        affected_genes = paste(unique(stats::na.omit(block$affected_gene_id)), collapse = "; "),
        affected_gene_count = dplyr::n_distinct(stats::na.omit(block$affected_gene_id)),
        highest_essentiality_class = highest_essentiality_class,
        max_target_site_score = max_target_site_score,
        mean_target_site_score = mean_target_site_score,
        max_essentiality_penalty = max_essentiality_penalty,
        representative_motifs = paste(unique(stats::na.omit(block$motif_seq)), collapse = "; ")
      )
    }

    start_idx <- 1L
    if (nrow(group) > 1L) {
      for (j in 2:nrow(group)) {
        if (group$site_start[[j]] <= (cur_end + as.integer(merge_gap_bp))) {
          cur_end <- max(cur_end, group$site_end[[j]])
        } else {
          region_rows[[region_idx]] <- emit_region(start_idx:(j - 1L))
          region_idx <- region_idx + 1L
          start_idx <- j
          cur_start <- group$site_start[[j]]
          cur_end <- group$site_end[[j]]
        }
      }
    }
    region_rows[[region_idx]] <- emit_region(start_idx:nrow(group))
    region_idx <- region_idx + 1L
  }

  dplyr::bind_rows(region_rows) %>%
    dplyr::arrange(.data$target_family, .data$contig, .data$region_start, .data$region_end)
}

.dnmb_classify_element_evidence <- function(elements) {
  if (is.null(elements) || !nrow(elements)) {
    return(elements)
  }

  if (!"source_call_types" %in% colnames(elements)) {
    elements$source_call_types <- NA_character_
  }
  if (!"tir_found" %in% colnames(elements)) {
    elements$tir_found <- NA
  }
  if (!"tsd_found" %in% colnames(elements)) {
    elements$tsd_found <- NA
  }
  if (!"reference_hit_id" %in% colnames(elements)) {
    elements$reference_hit_id <- NA_character_
  }

  elements %>%
    dplyr::mutate(
      source_call_types = dplyr::coalesce(as.character(.data$source_call_types), ""),
      annotation_supported = grepl("annotation", .data$source_call_types, fixed = TRUE),
      sequence_supported = grepl("sequence", .data$source_call_types, fixed = TRUE),
      native_tir_supported = dplyr::coalesce(as.logical(.data$tir_found), FALSE),
      native_tsd_supported = dplyr::coalesce(as.logical(.data$tsd_found), FALSE),
      reference_supported = !is.na(.data$reference_hit_id) & nzchar(.data$reference_hit_id),
      validation_label = dplyr::case_when(
        .data$native_tir_supported & .data$native_tsd_supported ~ "TIR+TSD",
        .data$native_tir_supported ~ "TIR",
        .data$native_tsd_supported ~ "TSD",
        .data$reference_supported ~ "REF",
        TRUE ~ "none"
      ),
      evidence_overlap_class = dplyr::case_when(
        .data$annotation_supported & .data$sequence_supported ~ "annotation_sequence_overlap",
        .data$annotation_supported ~ "annotation_only",
        .data$sequence_supported ~ "sequence_only",
        TRUE ~ "weak"
      ),
      evidence_class = dplyr::case_when(
        .data$annotation_supported & .data$sequence_supported & (.data$native_tir_supported | .data$native_tsd_supported) ~ "annotation_sequence_native",
        .data$annotation_supported & .data$sequence_supported ~ "annotation_sequence_overlap",
        !.data$annotation_supported & .data$sequence_supported & (.data$native_tir_supported | .data$native_tsd_supported) ~ "sequence_native",
        !.data$annotation_supported & .data$sequence_supported ~ "sequence_only",
        .data$annotation_supported ~ "annotation_only",
        TRUE ~ "weak"
      ),
      primary_support_group = dplyr::case_when(
        .data$evidence_class %in% c("annotation_sequence_native", "sequence_native") ~ "native_supported",
        .data$evidence_class == "annotation_sequence_overlap" ~ "overlap_supported",
        .data$evidence_class == "annotation_only" ~ "annotation_only",
        .data$evidence_class == "sequence_only" ~ "sequence_only",
        TRUE ~ "weak"
      ),
      integrated_evidence_label = dplyr::case_when(
        .data$evidence_class == "annotation_sequence_native" ~ "ANN+SEQ+NAT",
        .data$evidence_class == "annotation_sequence_overlap" ~ "ANN+SEQ",
        .data$evidence_class == "sequence_native" ~ "SEQ+NAT",
        .data$evidence_class == "sequence_only" ~ "SEQ",
        .data$evidence_class == "annotation_only" ~ "ANN",
        TRUE ~ "WEAK"
      )
    )
}

.dnmb_summarize_evidence_labels <- function(labels) {
  labels <- stats::na.omit(as.character(labels))
  labels <- labels[nzchar(labels)]
  if (!length(labels)) {
    return(NA_character_)
  }

  level_order <- c("ANN", "ANN+SEQ", "ANN+SEQ+NAT", "SEQ", "SEQ+NAT", "WEAK")
  label_counts <- table(factor(labels, levels = level_order))
  label_counts <- label_counts[label_counts > 0]
  if (!length(label_counts)) {
    return(NA_character_)
  }
  paste0(names(label_counts), "=", as.integer(label_counts), collapse = "; ")
}

.dnmb_summarize_support_groups <- function(groups) {
  groups <- stats::na.omit(as.character(groups))
  groups <- groups[nzchar(groups)]
  if (!length(groups)) {
    return(NA_character_)
  }

  level_order <- c("native_supported", "overlap_supported", "annotation_only", "sequence_only", "weak")
  label_map <- c(
    "native_supported" = "NATIVE",
    "overlap_supported" = "OVERLAP",
    "annotation_only" = "ANN",
    "sequence_only" = "SEQ",
    "weak" = "WEAK"
  )
  group_counts <- table(factor(groups, levels = level_order))
  group_counts <- group_counts[group_counts > 0]
  if (!length(group_counts)) {
    return(NA_character_)
  }
  paste0(unname(label_map[names(group_counts)]), "=", as.integer(group_counts), collapse = "; ")
}

.dnmb_build_master_table <- function(elements, sequence_elements, target_models, targetable_regions, target_sites, variant_catalog) {
  elements <- .dnmb_classify_element_evidence(elements)
  model_summary <- target_models %>%
    dplyr::transmute(
      target_family = .data$family,
      target_model_id = .data$target_model_id,
      recognition_strategy = .data$recognition_strategy,
      model_type = .data$model_type,
      model_motif = .data$model_motif,
      observed_tsd_count = .data$observed_tsd_count,
      dominant_tsd_len = .data$dominant_tsd_len,
      observed_tsd_examples = .data$observed_tsd_examples,
      model_confidence = .data$model_confidence,
      transposition_modes = .data$transposition_modes,
      recognition_rule_source = .data$rule_source
    )

  element_summary <- if (nrow(elements)) {
    elements %>%
      dplyr::group_by(.data$element_family) %>%
      dplyr::summarise(
        source_element_count = dplyr::n(),
        source_elements = paste(unique(stats::na.omit(.data$element_id)), collapse = "; "),
        source_element_ranges = paste(unique(paste0(.data$contig, ":", .data$start, "-", .data$end)), collapse = "; "),
        source_confidences = paste(unique(stats::na.omit(.data$confidence)), collapse = "; "),
        source_call_types = paste(unique(stats::na.omit(.data$source_call_types)), collapse = "; "),
        merged_tsd_examples = paste(unique(stats::na.omit(.data$tsd_seq)), collapse = "; "),
        merged_reference_hits = paste(unique(stats::na.omit(.data$reference_hit_id)), collapse = "; "),
        annotation_support_count = sum(.data$annotation_supported, na.rm = TRUE),
        merged_sequence_support_count = sum(.data$sequence_supported, na.rm = TRUE),
        merged_native_support_count = sum(.data$native_tir_supported | .data$native_tsd_supported, na.rm = TRUE),
        annotation_only_count = sum(.data$evidence_class == "annotation_only", na.rm = TRUE),
        annotation_sequence_overlap_count = sum(.data$evidence_class == "annotation_sequence_overlap", na.rm = TRUE),
        annotation_sequence_native_count = sum(.data$evidence_class == "annotation_sequence_native", na.rm = TRUE),
        sequence_only_count = sum(.data$evidence_class == "sequence_only", na.rm = TRUE),
        sequence_native_count = sum(.data$evidence_class == "sequence_native", na.rm = TRUE),
        family_evidence_labels = paste(unique(stats::na.omit(.data$integrated_evidence_label)), collapse = "; "),
        family_evidence_summary = .dnmb_summarize_evidence_labels(.data$integrated_evidence_label),
        family_support_group_summary = .dnmb_summarize_support_groups(.data$primary_support_group),
        .groups = "drop"
      ) %>%
      dplyr::rename(target_family = "element_family")
  } else {
    tibble::tibble(
      target_family = character(),
      source_element_count = integer(),
      source_elements = character(),
      source_element_ranges = character(),
      source_confidences = character(),
      source_call_types = character(),
      merged_tsd_examples = character(),
      merged_reference_hits = character(),
      annotation_support_count = integer(),
      merged_sequence_support_count = integer(),
      merged_native_support_count = integer(),
      annotation_only_count = integer(),
      annotation_sequence_overlap_count = integer(),
      annotation_sequence_native_count = integer(),
      sequence_only_count = integer(),
      sequence_native_count = integer(),
      family_evidence_labels = character(),
      family_evidence_summary = character(),
      family_support_group_summary = character()
    )
  }

  sequence_summary <- if (nrow(sequence_elements)) {
    sequence_elements %>%
      dplyr::group_by(.data$element_family) %>%
      dplyr::summarise(
        sequence_supported_elements = paste(unique(stats::na.omit(.data$sequence_element_id)), collapse = "; "),
        tir_examples = paste(
          unique(stats::na.omit(paste0(
            .data$sequence_element_id, "@",
            .data$tir_start1, "-", .data$tir_end2,
            "(len=", .data$tir_len_bp, ",id=", sprintf("%.2f", dplyr::coalesce(.data$tir_identity, 0)), ")"
          ))),
          collapse = "; "
        ),
        sequence_tsd_examples = paste(unique(stats::na.omit(.data$tsd_seq)), collapse = "; "),
        .groups = "drop"
      ) %>%
      dplyr::rename(target_family = "element_family")
  } else {
    tibble::tibble(
      target_family = character(),
      sequence_supported_elements = character(),
      tir_examples = character(),
      sequence_tsd_examples = character()
    )
  }

  build_region_row <- function(region_row) {
    sites_in_region <- target_sites %>%
      dplyr::filter(
        .data$target_model_id == region_row$target_model_id[[1]],
        .data$target_family == region_row$target_family[[1]],
        .data$contig == region_row$contig[[1]],
        .data$site_start >= region_row$region_start[[1]],
        .data$site_end <= region_row$region_end[[1]]
      )
    variants_in_region <- variant_catalog %>%
      dplyr::filter(
        .data$variant_type == "transposition_insertion",
        .data$source_family == region_row$target_family[[1]],
        .data$contig == region_row$contig[[1]],
        .data$site_start >= region_row$region_start[[1]],
        .data$site_end <= region_row$region_end[[1]]
      )

    tibble::tibble(
      master_row_type = "targetable_region",
      target_family = region_row$target_family[[1]],
      target_model_id = region_row$target_model_id[[1]],
      contig = region_row$contig[[1]],
      region_start = region_row$region_start[[1]],
      region_end = region_row$region_end[[1]],
      region_size_bp = region_row$region_size_bp[[1]],
      region_type = region_row$region_type[[1]],
      recognition_strategy = region_row$recognition_strategy[[1]],
      model_type = region_row$model_type[[1]],
      model_motif = region_row$model_motif[[1]],
      transposition_modes = region_row$transposition_modes[[1]],
      target_site_count = region_row$target_site_count[[1]],
      representative_target_sites = if (nrow(sites_in_region)) paste(head(paste0(sites_in_region$site_start, "-", sites_in_region$site_end), 10L), collapse = "; ") else NA_character_,
      representative_target_motifs = region_row$representative_motifs[[1]],
      context_classes = region_row$context_classes[[1]],
      affected_genes = region_row$affected_genes[[1]],
      affected_gene_count = region_row$affected_gene_count[[1]],
      highest_essentiality_class = region_row$highest_essentiality_class[[1]],
      max_target_site_score = region_row$max_target_site_score[[1]],
      mean_target_site_score = region_row$mean_target_site_score[[1]],
      max_essentiality_penalty = region_row$max_essentiality_penalty[[1]],
      variant_rows_supported = nrow(variants_in_region),
      top_variant_ids = if (nrow(variants_in_region)) paste(head(variants_in_region$variant_id[order(-variants_in_region$viability_score)], 10L), collapse = "; ") else NA_character_,
      top_variant_viability = if (nrow(variants_in_region)) max(variants_in_region$viability_score, na.rm = TRUE) else NA_real_
    )
  }

  region_master <- if (nrow(targetable_regions)) {
    dplyr::bind_rows(lapply(seq_len(nrow(targetable_regions)), function(i) build_region_row(targetable_regions[i, , drop = FALSE])))
  } else {
    tibble::tibble()
  }

  if (nrow(region_master)) {
    region_master <- region_master %>%
      dplyr::left_join(model_summary, by = c("target_family", "target_model_id"), suffix = c("", "_model")) %>%
      dplyr::left_join(element_summary, by = "target_family") %>%
      dplyr::left_join(sequence_summary, by = "target_family")
  }

  families_without_regions <- setdiff(model_summary$target_family, unique(targetable_regions$target_family))
  no_region_rows <- if (length(families_without_regions)) {
    model_summary %>%
      dplyr::filter(.data$target_family %in% families_without_regions) %>%
      dplyr::left_join(element_summary, by = "target_family") %>%
      dplyr::left_join(sequence_summary, by = "target_family") %>%
      dplyr::mutate(
        master_row_type = "no_reliable_targetable_region",
        contig = NA_character_,
        region_start = NA_integer_,
        region_end = NA_integer_,
        region_size_bp = NA_integer_,
        region_type = NA_character_,
        target_site_count = 0L,
        representative_target_sites = NA_character_,
        representative_target_motifs = NA_character_,
        context_classes = NA_character_,
        affected_genes = NA_character_,
        affected_gene_count = 0L,
        highest_essentiality_class = NA_character_,
        max_target_site_score = NA_real_,
        mean_target_site_score = NA_real_,
        max_essentiality_penalty = NA_real_,
        variant_rows_supported = 0L,
        top_variant_ids = NA_character_,
        top_variant_viability = NA_real_
      )
  } else {
    tibble::tibble()
  }

  dplyr::bind_rows(region_master, no_region_rows) %>%
    dplyr::select(
      .data$master_row_type,
      .data$target_family,
      .data$target_model_id,
      .data$contig,
      .data$region_start,
      .data$region_end,
      .data$region_size_bp,
      .data$region_type,
      .data$recognition_strategy,
      .data$model_type,
      .data$model_motif,
      .data$observed_tsd_count,
      .data$dominant_tsd_len,
      .data$observed_tsd_examples,
      .data$model_confidence,
      .data$recognition_rule_source,
      .data$transposition_modes,
      .data$source_element_count,
      .data$source_elements,
      .data$source_element_ranges,
      .data$source_confidences,
      .data$source_call_types,
      .data$annotation_support_count,
      .data$merged_sequence_support_count,
      .data$merged_native_support_count,
      .data$annotation_only_count,
      .data$annotation_sequence_overlap_count,
      .data$annotation_sequence_native_count,
      .data$sequence_only_count,
      .data$sequence_native_count,
      .data$family_evidence_labels,
      .data$family_evidence_summary,
      .data$family_support_group_summary,
      .data$merged_tsd_examples,
      .data$merged_reference_hits,
      .data$sequence_supported_elements,
      .data$tir_examples,
      .data$sequence_tsd_examples,
      .data$target_site_count,
      .data$representative_target_sites,
      .data$representative_target_motifs,
      .data$context_classes,
      .data$affected_genes,
      .data$affected_gene_count,
      .data$highest_essentiality_class,
      .data$max_target_site_score,
      .data$mean_target_site_score,
      .data$max_essentiality_penalty,
      .data$variant_rows_supported,
      .data$top_variant_ids,
      .data$top_variant_viability
    ) %>%
    dplyr::arrange(.data$target_family, .data$contig, .data$region_start)
}

.dnmb_build_compact_summary_table <- function(
  elements,
  sequence_elements,
  target_models,
  targetable_regions,
  variant_catalog,
  comparative = NULL
) {
  families <- sort(unique(c(
    elements$element_family,
    sequence_elements$element_family,
    target_models$family,
    targetable_regions$target_family
  )))
  families <- families[!is.na(families) & nzchar(families)]
  if (!length(families)) {
    return(tibble::tibble())
  }

  elements <- .dnmb_classify_element_evidence(elements)

  ann_tbl <- if (nrow(elements)) {
    elements %>%
      dplyr::group_by(.data$element_family) %>%
      dplyr::summarise(
        merged_element_count = dplyr::n(),
        annotation_support_count = sum(.data$annotation_supported, na.rm = TRUE),
        merged_sequence_support_count = sum(.data$sequence_supported, na.rm = TRUE),
        merged_native_support_count = sum(.data$native_tir_supported | .data$native_tsd_supported, na.rm = TRUE),
        annotation_only_count = sum(.data$evidence_class == "annotation_only", na.rm = TRUE),
        annotation_sequence_overlap_count = sum(.data$evidence_class == "annotation_sequence_overlap", na.rm = TRUE),
        annotation_sequence_native_count = sum(.data$evidence_class == "annotation_sequence_native", na.rm = TRUE),
        sequence_only_count = sum(.data$evidence_class == "sequence_only", na.rm = TRUE),
        sequence_native_count = sum(.data$evidence_class == "sequence_native", na.rm = TRUE),
        integrated_evidence_labels = paste(unique(stats::na.omit(.data$integrated_evidence_label)), collapse = "; "),
        integrated_evidence_summary = .dnmb_summarize_evidence_labels(.data$integrated_evidence_label),
        support_group_summary = .dnmb_summarize_support_groups(.data$primary_support_group),
        primary_support_group = dplyr::case_when(
          sum(.data$primary_support_group == "native_supported", na.rm = TRUE) > 0 ~ "native_supported",
          sum(.data$primary_support_group == "overlap_supported", na.rm = TRUE) > 0 ~ "overlap_supported",
          sum(.data$primary_support_group == "annotation_only", na.rm = TRUE) > 0 ~ "annotation_only",
          sum(.data$primary_support_group == "sequence_only", na.rm = TRUE) > 0 ~ "sequence_only",
          TRUE ~ "weak"
        ),
        merged_confidence = paste(unique(stats::na.omit(.data$confidence)), collapse = "; "),
        .groups = "drop"
      ) %>%
      dplyr::rename(family = .data$element_family)
  } else {
    tibble::tibble(
      family = character(),
      merged_element_count = integer(),
      annotation_support_count = integer(),
      merged_sequence_support_count = integer(),
      merged_native_support_count = integer(),
      annotation_only_count = integer(),
      annotation_sequence_overlap_count = integer(),
      annotation_sequence_native_count = integer(),
      sequence_only_count = integer(),
      sequence_native_count = integer(),
      integrated_evidence_labels = character(),
      integrated_evidence_summary = character(),
      primary_support_group = character(),
      support_group_summary = character(),
      merged_confidence = character()
    )
  }

  seq_tbl <- if (nrow(sequence_elements)) {
    sequence_elements %>%
      dplyr::group_by(.data$element_family) %>%
      dplyr::summarise(
        sequence_support_count = dplyr::n(),
        tir_support_count = sum(isTRUE(.data$tir_found), na.rm = TRUE),
        tsd_support_count = sum(isTRUE(.data$tsd_found), na.rm = TRUE),
        top_tir_len = suppressWarnings(max(.data$tir_len_bp, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(top_tir_len = ifelse(is.finite(.data$top_tir_len), .data$top_tir_len, NA)) %>%
      dplyr::rename(family = .data$element_family)
  } else {
    tibble::tibble(family = character(), sequence_support_count = integer(), tir_support_count = integer(), tsd_support_count = integer(), top_tir_len = numeric())
  }

  model_tbl <- if (nrow(target_models)) {
    target_models %>%
      dplyr::group_by(.data$family) %>%
      dplyr::summarise(
        recognition_mode = paste(unique(stats::na.omit(.data$recognition_mode)), collapse = "; "),
        model_types = paste(unique(stats::na.omit(.data$model_type)), collapse = "; "),
        model_confidence = paste(unique(stats::na.omit(.data$model_confidence)), collapse = "; "),
        model_motif = paste(unique(stats::na.omit(.data$model_motif)), collapse = "; "),
        observed_tsd_count = max(.data$observed_tsd_count, na.rm = TRUE),
        enriched_tsd_examples = paste(unique(stats::na.omit(.data$enriched_tsd_examples)), collapse = "; "),
        supplementary_strategy = paste(unique(stats::na.omit(.data$supplementary_strategy)), collapse = "; "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(observed_tsd_count = ifelse(is.finite(.data$observed_tsd_count), .data$observed_tsd_count, 0)) 
  } else {
    tibble::tibble(family = character(), recognition_mode = character(), model_types = character(), model_confidence = character(), model_motif = character(), observed_tsd_count = integer(), enriched_tsd_examples = character(), supplementary_strategy = character())
  }

  region_tbl <- if (nrow(targetable_regions)) {
    targetable_regions %>%
      dplyr::group_by(.data$target_family) %>%
      dplyr::summarise(
        targetable_region_count = dplyr::n(),
        max_target_score = max(.data$max_target_site_score, na.rm = TRUE),
        highest_region_essentiality = paste(unique(stats::na.omit(.data$highest_essentiality_class)), collapse = "; "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(max_target_score = ifelse(is.finite(.data$max_target_score), .data$max_target_score, NA)) %>%
      dplyr::rename(family = .data$target_family)
  } else {
    tibble::tibble(family = character(), targetable_region_count = integer(), max_target_score = numeric(), highest_region_essentiality = character())
  }

  variant_tbl <- if (nrow(variant_catalog)) {
    variant_catalog %>%
      dplyr::group_by(.data$source_family) %>%
      dplyr::summarise(
        variant_count = dplyr::n(),
        max_variant_viability = max(.data$viability_score, na.rm = TRUE),
        chronology_hint = paste(unique(stats::na.omit(.data$variant_type)), collapse = "; "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(max_variant_viability = ifelse(is.finite(.data$max_variant_viability), .data$max_variant_viability, NA)) %>%
      dplyr::rename(family = .data$source_family)
  } else {
    tibble::tibble(family = character(), variant_count = integer(), max_variant_viability = numeric(), chronology_hint = character())
  }

  comp_tbl <- if (!is.null(comparative) && nrow(comparative$hotspots)) {
    comparative$hotspots %>%
      dplyr::group_by(.data$family) %>%
      dplyr::summarise(
        comparative_hotspot_count = dplyr::n(),
        comparative_support_score = max(.data$comparative_support_score, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        comparative$chronology %>%
          dplyr::group_by(.data$family) %>%
          dplyr::summarise(
            chronology_class = paste(unique(stats::na.omit(.data$chronology_class)), collapse = "; "),
            proposed_order = paste(unique(stats::na.omit(.data$proposed_order)), collapse = " | "),
            .groups = "drop"
          ),
        by = "family"
      )
  } else {
    tibble::tibble(
      family = character(),
      comparative_hotspot_count = integer(),
      comparative_support_score = numeric(),
      chronology_class = character(),
      proposed_order = character()
    )
  }

  tibble::tibble(family = families) %>%
    dplyr::left_join(ann_tbl, by = "family") %>%
    dplyr::left_join(seq_tbl, by = "family") %>%
    dplyr::left_join(model_tbl, by = "family") %>%
    dplyr::left_join(region_tbl, by = "family") %>%
    dplyr::left_join(variant_tbl, by = "family") %>%
    dplyr::left_join(comp_tbl, by = "family") %>%
    dplyr::mutate(
      annotation_support_count = dplyr::coalesce(.data$annotation_support_count, 0L),
      merged_sequence_support_count = dplyr::coalesce(.data$merged_sequence_support_count, 0L),
      merged_native_support_count = dplyr::coalesce(.data$merged_native_support_count, 0L),
      annotation_only_count = dplyr::coalesce(.data$annotation_only_count, 0L),
      annotation_sequence_overlap_count = dplyr::coalesce(.data$annotation_sequence_overlap_count, 0L),
      annotation_sequence_native_count = dplyr::coalesce(.data$annotation_sequence_native_count, 0L),
      sequence_only_count = dplyr::coalesce(.data$sequence_only_count, 0L),
      sequence_native_count = dplyr::coalesce(.data$sequence_native_count, 0L),
      sequence_support_count = dplyr::coalesce(.data$sequence_support_count, 0L),
      tir_support_count = dplyr::coalesce(.data$tir_support_count, 0L),
      tsd_support_count = dplyr::coalesce(.data$tsd_support_count, 0L),
      observed_tsd_count = dplyr::coalesce(.data$observed_tsd_count, 0L),
      targetable_region_count = dplyr::coalesce(.data$targetable_region_count, 0L),
      variant_count = dplyr::coalesce(.data$variant_count, 0L),
      comparative_hotspot_count = dplyr::coalesce(.data$comparative_hotspot_count, 0L),
      evidence_badges = vapply(
        seq_len(dplyr::n()),
        function(i) {
          badges <- character()
          if (.data$annotation_support_count[[i]] > 0) badges <- c(badges, "ANN")
          if (.data$sequence_support_count[[i]] > 0) badges <- c(badges, "SEQ")
          if (.data$tsd_support_count[[i]] > 0 || .data$observed_tsd_count[[i]] > 0) badges <- c(badges, "TSD")
          if (.data$tir_support_count[[i]] > 0) badges <- c(badges, "TIR")
          if (!is.na(.data$enriched_tsd_examples[[i]]) && nzchar(.data$enriched_tsd_examples[[i]])) badges <- c(badges, "ENR")
          if (.data$comparative_hotspot_count[[i]] > 0) badges <- c(badges, "CMP")
          if (!length(badges)) "-" else paste(badges, collapse = "|")
        },
        character(1)
      ),
      evidence_tier = dplyr::case_when(
        (.data$annotation_sequence_native_count > 0 | .data$sequence_native_count > 0) & .data$comparative_hotspot_count > 0 ~ "high_confidence_multilayer",
        .data$annotation_sequence_native_count > 0 | .data$sequence_native_count > 0 ~ "sequence_native_supported",
        .data$annotation_sequence_overlap_count > 0 ~ "annotation_sequence_supported",
        .data$annotation_support_count > 0 ~ "annotation_only",
        .data$sequence_only_count > 0 ~ "sequence_only",
        TRUE ~ "weak_or_unresolved"
      )
    ) %>%
    dplyr::arrange(.data$family)
}

.dnmb_score_context_candidates_one_contig <- function(seq, contig, model, context_sites = NULL, flank_window = 12L) {
  if (is.null(context_sites) || !nrow(context_sites)) {
    return(tibble::tibble())
  }

  candidates <- context_sites
  if (identical(model$model_type[[1]], "context_only")) {
    motif_hits <- candidates %>%
      dplyr::mutate(
        motif_seq = "context_only",
        motif_match_score = pmax(.data$motif_match_score, 0.25)
      )
  } else {
    motifs <- unique(trimws(unlist(strsplit(dplyr::coalesce(model$model_motif[[1]], ""), ";", fixed = TRUE))))
    motifs <- motifs[nzchar(motifs)]
    if (!length(motifs)) {
      return(tibble::tibble())
    }
    regions <- .dnmb_candidate_windows_from_sites(candidates$site_start, nchar(seq), flank_window = flank_window)
    motif_hits <- .dnmb_find_model_hits(seq = seq, model = model, regions = regions)
    if (!nrow(motif_hits)) {
      return(tibble::tibble())
    }
    motif_hits <- .dnmb_map_hits_to_context_sites(
      motif_hits = motif_hits,
      context_sites = candidates,
      max_distance = flank_window
    )
  }

  if (!nrow(motif_hits)) {
    return(tibble::tibble())
  }

  motif_hits %>%
    dplyr::mutate(
      target_model_id = model$target_model_id[[1]],
      target_family = model$family[[1]],
      contig = contig,
      motif_source = model$motif_source[[1]],
      transposition_modes = model$transposition_modes[[1]],
      context_pref = dplyr::case_when(
        .data$context_class == "intergenic" ~ dplyr::coalesce(model$preferred_intergenic[[1]], 0),
        .data$context_class == "promoter" ~ dplyr::coalesce(model$preferred_promoter[[1]], 0),
        TRUE ~ dplyr::coalesce(model$preferred_gene_body[[1]], 0)
      ),
      essentiality_penalty = vapply(
        seq_len(dplyr::n()),
        function(i) .dnmb_context_essentiality_penalty(.data$context_class[[i]], .data$affected_gene_essentiality_class[[i]]),
        numeric(1)
      ),
      permissive_bonus = dplyr::case_when(
        .data$context_class == "intergenic" ~ 0.35,
        .data$context_class == "promoter" ~ 0.2,
        TRUE ~ 0.05
      ),
      target_site_score = pmax(0, pmin(1,
        0.55 * .data$motif_match_score +
          0.25 * .data$context_pref +
          0.20 * .data$permissive_bonus -
          0.60 * .data$essentiality_penalty
      )),
      target_site_confidence = dplyr::case_when(
        model$model_confidence[[1]] == "high" & .data$target_site_score >= 0.5 ~ "high",
        .data$target_site_score >= 0.25 ~ "medium",
        TRUE ~ "low"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$target_site_score), dplyr::desc(.data$motif_match_score), .data$site_start)
}

.dnmb_candidate_windows_from_sites <- function(site_starts, seq_len, flank_window = 12L) {
  site_starts <- sort(unique(as.integer(site_starts[!is.na(site_starts)])))
  if (!length(site_starts)) {
    return(tibble::tibble(region_start = integer(), region_end = integer()))
  }
  starts <- pmax(1L, site_starts - as.integer(flank_window))
  ends <- pmin(as.integer(seq_len), site_starts + as.integer(flank_window))
  regions <- tibble::tibble(region_start = starts, region_end = ends) %>%
    dplyr::arrange(.data$region_start, .data$region_end)

  merged <- list()
  cur_start <- regions$region_start[[1]]
  cur_end <- regions$region_end[[1]]
  if (nrow(regions) > 1L) {
    for (i in 2:nrow(regions)) {
      if (regions$region_start[[i]] <= (cur_end + 1L)) {
        cur_end <- max(cur_end, regions$region_end[[i]])
      } else {
        merged[[length(merged) + 1L]] <- tibble::tibble(region_start = cur_start, region_end = cur_end)
        cur_start <- regions$region_start[[i]]
        cur_end <- regions$region_end[[i]]
      }
    }
  }
  merged[[length(merged) + 1L]] <- tibble::tibble(region_start = cur_start, region_end = cur_end)
  dplyr::bind_rows(merged)
}

.dnmb_map_hits_to_context_sites <- function(motif_hits, context_sites, max_distance = 12L) {
  if (!nrow(motif_hits) || !nrow(context_sites)) {
    return(tibble::tibble())
  }

  context_sites <- context_sites %>%
    dplyr::arrange(.data$site_start)

  cs <- context_sites$site_start
  hs <- motif_hits$site_start
  ncs <- length(cs)
  left_idx <- findInterval(hs, cs)
  left_idx[left_idx < 1L] <- 1L
  right_idx <- pmin(left_idx + 1L, ncs)
  dist_left <- abs(hs - cs[left_idx])
  dist_right <- abs(hs - cs[right_idx])
  chosen_idx <- ifelse(dist_left <= dist_right, left_idx, right_idx)
  chosen_dist <- pmin(dist_left, dist_right)
  keep <- which(chosen_dist <= as.integer(max_distance))

  if (!length(keep)) {
    return(tibble::tibble())
  }

  context_sel <- context_sites[chosen_idx[keep], c(
    "context_class",
    "affected_gene_id",
    "affected_gene",
    "affected_product",
    "affected_gene_essentiality_score",
    "affected_gene_essentiality_class",
    "phenotype_categories",
    "phenotype_prediction"
  )]

  dplyr::bind_cols(
    motif_hits[keep, , drop = FALSE],
    context_sel
  )
}

.dnmb_scan_target_sites_one_contig <- function(seq, contig, model, genes_contig = NULL, context_sites = NULL, permissive_regions = NULL) {
  if (is.null(genes_contig)) {
    genes_contig <- tibble::tibble()
  }
  if (identical(model$model_type[[1]], "none")) {
    return(tibble::tibble())
  }
  if (identical(model$model_type[[1]], "context_only")) {
    motif_hits <- if (is.null(context_sites)) tibble::tibble() else context_sites
  } else {
    use_permissive <- (
      (
        identical(model$motif_source[[1]], "fallback_rule") ||
          (!is.na(model$dominant_tsd_len[[1]]) && model$dominant_tsd_len[[1]] <= 3L)
      ) &&
        (
          identical(model$model_type[[1]], "at_rich") ||
            (identical(model$model_type[[1]], "motif") && !is.na(model$model_motif[[1]]) && nchar(gsub("[^A-Z]", "", model$model_motif[[1]])) < 4L) ||
            identical(model$model_type[[1]], "empirical_tsd")
        )
    )

    motif_hits <- .dnmb_find_model_hits(
      seq = seq,
      model = model,
      regions = if (isTRUE(use_permissive)) permissive_regions else NULL
    )
  }
  if (!nrow(motif_hits)) {
    return(tibble::tibble())
  }
  motif_hits <- .dnmb_prune_target_hits(motif_hits, model)

  if (all(c("context_class", "affected_gene_id", "affected_gene_essentiality_class") %in% colnames(motif_hits))) {
    classified <- motif_hits
  } else {
    classified <- .dnmb_annotate_hits_with_context(
      motif_hits = motif_hits,
      genes_contig = genes_contig,
      context_sites = context_sites
    )
  }

  classified %>%
    dplyr::mutate(
      target_model_id = model$target_model_id[[1]],
      target_family = model$family[[1]],
      contig = contig,
      motif_source = model$motif_source[[1]],
      transposition_modes = model$transposition_modes[[1]],
      context_pref = dplyr::case_when(
        .data$context_class == "intergenic" ~ dplyr::coalesce(model$preferred_intergenic[[1]], 0),
        .data$context_class == "promoter" ~ dplyr::coalesce(model$preferred_promoter[[1]], 0),
        TRUE ~ dplyr::coalesce(model$preferred_gene_body[[1]], 0)
      ),
      essentiality_penalty = vapply(
        seq_len(dplyr::n()),
        function(i) .dnmb_context_essentiality_penalty(.data$context_class[[i]], .data$affected_gene_essentiality_class[[i]]),
        numeric(1)
      ),
      permissive_bonus = dplyr::case_when(
        .data$context_class == "intergenic" ~ 0.35,
        .data$context_class == "promoter" ~ 0.2,
        TRUE ~ 0.05
      ),
      target_site_score = pmax(0, pmin(1,
        0.55 * .data$motif_match_score +
          0.25 * .data$context_pref +
          0.20 * .data$permissive_bonus -
          0.60 * .data$essentiality_penalty
      )),
      target_site_confidence = dplyr::case_when(
        model$model_confidence[[1]] == "high" & .data$target_site_score >= 0.5 ~ "high",
        .data$target_site_score >= 0.25 ~ "medium",
        TRUE ~ "low"
      )
    ) %>%
    dplyr::select(
      target_model_id,
      target_family,
      contig,
      site_start,
      site_end,
      motif_seq,
      motif_match_score,
      motif_source,
      transposition_modes,
      context_class,
      affected_gene_id,
      affected_gene,
      affected_product,
      affected_gene_essentiality_score,
      affected_gene_essentiality_class,
      phenotype_categories,
      phenotype_prediction,
      essentiality_penalty,
      target_site_score,
      target_site_confidence
    )
}

.dnmb_annotate_hits_with_context <- function(motif_hits, genes_contig, context_sites = NULL) {
  if (!nrow(motif_hits)) {
    return(motif_hits)
  }

  hdt <- data.table::as.data.table(motif_hits)
  hdt[, hit_id := .I]

  if (!nrow(genes_contig)) {
    hdt[, `:=`(
      context_class = "intergenic",
      affected_gene_id = NA_character_,
      affected_gene = NA_character_,
      affected_product = NA_character_,
      affected_gene_essentiality_score = NA_real_,
      affected_gene_essentiality_class = NA_character_,
      phenotype_categories = NA_character_,
      phenotype_prediction = NA_character_
    )]
    return(tibble::as_tibble(hdt[, !c("hit_id"), with = FALSE]))
  }

  gdt <- data.table::as.data.table(genes_contig)[, .(
    start,
    end,
    gene_id,
    gene,
    product,
    essentiality_score,
    essentiality_class,
    gene_phenotype_categories,
    gene_phenotype_prediction
  )]
  hdt[, start_hit := site_start]
  hdt[, end_hit := site_end]
  data.table::setkey(hdt, start_hit, end_hit)
  data.table::setkey(gdt, start, end)

  overlaps <- data.table::foverlaps(
    hdt,
    gdt,
    by.x = c("start_hit", "end_hit"),
    by.y = c("start", "end"),
    type = "any",
    nomatch = 0L
  )

  if (nrow(overlaps)) {
    overlaps <- overlaps[order(hit_id, -essentiality_score, start)]
    overlaps_best <- overlaps[!duplicated(hit_id), .(
      hit_id,
      context_class = "gene_body",
      affected_gene_id = gene_id,
      affected_gene = gene,
      affected_product = product,
      affected_gene_essentiality_score = essentiality_score,
      affected_gene_essentiality_class = essentiality_class,
      phenotype_categories = gene_phenotype_categories,
      phenotype_prediction = gene_phenotype_prediction
    )]
    hdt <- overlaps_best[hdt, on = "hit_id"]
  } else {
    hdt[, `:=`(
      context_class = NA_character_,
      affected_gene_id = NA_character_,
      affected_gene = NA_character_,
      affected_product = NA_character_,
      affected_gene_essentiality_score = NA_real_,
      affected_gene_essentiality_class = NA_character_,
      phenotype_categories = NA_character_,
      phenotype_prediction = NA_character_
    )]
  }

  remaining <- hdt[is.na(context_class)]
  if (nrow(remaining)) {
    if (!is.null(context_sites) && nrow(context_sites)) {
      cdt <- data.table::as.data.table(context_sites)[, .(
        site_start,
        context_class,
        affected_gene_id,
        affected_gene,
        affected_product,
        affected_gene_essentiality_score,
        affected_gene_essentiality_class,
        phenotype_categories,
        phenotype_prediction
      )]
      data.table::setkey(cdt, site_start)
      nearest <- cdt[remaining[, .(hit_id, site_start)], on = "site_start", roll = "nearest"]
      nearest <- nearest[, .(
        hit_id,
        context_class,
        affected_gene_id,
        affected_gene,
        affected_product,
        affected_gene_essentiality_score,
        affected_gene_essentiality_class,
        phenotype_categories,
        phenotype_prediction
      )]
      hdt <- nearest[hdt, on = "hit_id"]
      for (col in c("context_class", "affected_gene_id", "affected_gene", "affected_product",
                    "affected_gene_essentiality_score", "affected_gene_essentiality_class",
                    "phenotype_categories", "phenotype_prediction")) {
        data.table::set(
          hdt,
          j = col,
          value = dplyr::coalesce(hdt[[col]], hdt[[paste0("i.", col)]])
        )
      }
      drop_cols <- names(hdt)[startsWith(names(hdt), "i.")]
      if (length(drop_cols)) {
        hdt[, (drop_cols) := NULL]
      }
    }

    hdt[is.na(context_class), `:=`(
      context_class = "intergenic",
      affected_gene_id = NA_character_,
      affected_gene = NA_character_,
      affected_product = NA_character_,
      affected_gene_essentiality_score = NA_real_,
      affected_gene_essentiality_class = NA_character_,
      phenotype_categories = NA_character_,
      phenotype_prediction = NA_character_
    )]
  }

  hdt[, c("hit_id", "start_hit", "end_hit") := NULL]
  tibble::as_tibble(hdt)
}

.dnmb_prepare_context_cache <- function(metadata, genes, promoter_bp = 150L) {
  contigs <- metadata$contig
  seq_lens <- stats::setNames(metadata$sequence_length_bp, metadata$contig)
  genes_by_contig <- stats::setNames(vector("list", length(contigs)), contigs)
  context_sites_by_contig <- stats::setNames(vector("list", length(contigs)), contigs)
  permissive_regions_by_contig <- stats::setNames(vector("list", length(contigs)), contigs)

  for (contig in contigs) {
    genes_contig <- genes %>%
      dplyr::filter(.data$contig == .env$contig) %>%
      dplyr::arrange(.data$start, .data$end)
    genes_by_contig[[contig]] <- genes_contig
    context_sites_by_contig[[contig]] <- .dnmb_generate_context_only_sites_from_genes(
      contig_genes = genes_contig,
      seq_len = seq_lens[[contig]],
      promoter_bp = promoter_bp
    )
    permissive_regions_by_contig[[contig]] <- .dnmb_build_permissive_regions_from_genes(
      contig_genes = genes_contig,
      seq_len = seq_lens[[contig]],
      promoter_bp = promoter_bp
    )
  }

  list(
    genes_by_contig = genes_by_contig,
    context_sites_by_contig = context_sites_by_contig,
    permissive_regions_by_contig = permissive_regions_by_contig
  )
}

.dnmb_find_model_hits <- function(seq, model, regions = NULL) {
  model_type <- model$model_type[[1]]
  motif <- model$model_motif[[1]]
  if (is.na(model_type) || model_type == "none") {
    return(tibble::tibble())
  }

  if (model_type %in% c("empirical_tsd", "motif", "motif_set")) {
    motifs <- unique(trimws(unlist(strsplit(dplyr::coalesce(motif, ""), ";", fixed = TRUE))))
    motifs <- motifs[nzchar(motifs)]
    rows <- lapply(motifs, function(one_motif) {
      if (grepl("^[ACGT]+$", one_motif)) {
        return(.dnmb_scan_exact_regions(
          seq = seq,
          motif = one_motif,
          match_score = dplyr::coalesce(model$motif_specificity[[1]], 0.4),
          regions = regions
        ))
      }
      regex <- .dnmb_iupac_to_regex(one_motif)
      if (!nzchar(regex)) {
        return(NULL)
      }
      .dnmb_scan_regex_regions(
        seq = seq,
        regex = regex,
        motif_len = nchar(one_motif),
        match_score = dplyr::coalesce(model$motif_specificity[[1]], 0.4),
        regions = regions
      )
    })
    return(dplyr::bind_rows(rows))
  }

  if (model_type == "at_rich") {
    width <- 8L
    if (nchar(seq) < width) {
      return(tibble::tibble())
    }
    if (is.null(regions) || !nrow(regions)) {
      regions <- tibble::tibble(region_start = 1L, region_end = nchar(seq))
    }
    rows <- lapply(seq_len(nrow(regions)), function(i) {
      rstart <- regions$region_start[[i]]
      rend <- regions$region_end[[i]]
      if ((rend - rstart + 1L) < width) {
        return(NULL)
      }
      starts <- rstart:(rend - width + 1L)
      windows <- substring(seq, starts, starts + width - 1L)
      at_frac <- stringr::str_count(windows, "[AT]") / width
      keep <- which(at_frac >= 0.75)
      if (!length(keep)) {
        return(NULL)
      }
      tibble::tibble(
        site_start = starts[keep],
        site_end = starts[keep] + width - 1L,
        motif_seq = windows[keep],
        motif_match_score = at_frac[keep]
      )
    })
    return(dplyr::bind_rows(rows))
  }

  tibble::tibble()
}

.dnmb_scan_regex_regions <- function(seq, regex, motif_len, match_score, regions = NULL) {
  if (is.null(regions) || !nrow(regions)) {
    regions <- tibble::tibble(region_start = 1L, region_end = nchar(seq))
  }
  rows <- lapply(seq_len(nrow(regions)), function(i) {
    rstart <- regions$region_start[[i]]
    rend <- regions$region_end[[i]]
    if ((rend - rstart + 1L) < motif_len) {
      return(NULL)
    }
    subseq <- substring(seq, rstart, rend)
    starts_local <- gregexpr(paste0("(?=", regex, ")"), subseq, perl = TRUE)[[1]]
    if (length(starts_local) == 1L && starts_local[1] == -1L) {
      return(NULL)
    }
    starts <- as.integer(starts_local) + rstart - 1L
    ends <- starts + motif_len - 1L
    tibble::tibble(
      site_start = starts,
      site_end = ends,
      motif_seq = substring(seq, starts, ends),
      motif_match_score = match_score
    )
  })
  dplyr::bind_rows(rows)
}

.dnmb_scan_exact_regions <- function(seq, motif, match_score, regions = NULL) {
  if (is.null(regions) || !nrow(regions)) {
    regions <- tibble::tibble(region_start = 1L, region_end = nchar(seq))
  }
  dna <- Biostrings::DNAString(seq)
  patt <- Biostrings::DNAString(motif)
  rows <- lapply(seq_len(nrow(regions)), function(i) {
    rstart <- regions$region_start[[i]]
    rend <- regions$region_end[[i]]
    if ((rend - rstart + 1L) < nchar(motif)) {
      return(NULL)
    }
    view <- Biostrings::subseq(dna, start = rstart, end = rend)
    matches <- Biostrings::matchPattern(patt, view, fixed = TRUE)
    if (length(matches) == 0L) {
      return(NULL)
    }
    starts <- Biostrings::start(matches) + rstart - 1L
    ends <- starts + nchar(motif) - 1L
    tibble::tibble(
      site_start = starts,
      site_end = ends,
      motif_seq = substring(seq, starts, ends),
      motif_match_score = match_score
    )
  })
  dplyr::bind_rows(rows)
}

.dnmb_prune_target_hits <- function(motif_hits, model) {
  if (!nrow(motif_hits)) {
    return(motif_hits)
  }

  motif_len <- median(nchar(dplyr::coalesce(motif_hits$motif_seq, "")))
  if ((is.finite(motif_len) && motif_len <= 3) || nrow(motif_hits) > 5000L) {
    bin_size <- if (nrow(motif_hits) > 20000L) 50L else 25L
    motif_hits <- motif_hits %>%
      dplyr::mutate(site_bin = (.data$site_start - 1L) %/% bin_size) %>%
      dplyr::group_by(.data$site_bin) %>%
      dplyr::slice_max(order_by = .data$motif_match_score, n = 1L, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(-site_bin)
  }

  if (nrow(motif_hits) > 1000L) {
    motif_hits <- motif_hits %>%
      dplyr::arrange(dplyr::desc(.data$motif_match_score), .data$site_start) %>%
      dplyr::slice_head(n = 1000L)
  }

  motif_hits
}

.dnmb_build_permissive_regions <- function(contig, genes, seq_len, promoter_bp = 150L) {
  contig_genes <- genes %>%
    dplyr::filter(.data$contig == contig) %>%
    dplyr::arrange(.data$start, .data$end)
  .dnmb_build_permissive_regions_from_genes(contig_genes = contig_genes, seq_len = seq_len, promoter_bp = promoter_bp)
}

.dnmb_build_permissive_regions_from_genes <- function(contig_genes, seq_len, promoter_bp = 150L) {
  if (!nrow(contig_genes)) {
    return(tibble::tibble(region_start = 1L, region_end = seq_len))
  }

  regions <- list()

  promoter_regions <- contig_genes %>%
    dplyr::transmute(
      region_start = dplyr::if_else(.data$strand == "+", pmax(1L, .data$start - promoter_bp), .data$end + 1L),
      region_end = dplyr::if_else(.data$strand == "+", .data$start - 1L, pmin(seq_len, .data$end + promoter_bp))
    ) %>%
    dplyr::filter(.data$region_end >= .data$region_start)
  regions[[length(regions) + 1L]] <- promoter_regions

  if (nrow(contig_genes) >= 2L) {
    gap_regions <- tibble::tibble(
      region_start = contig_genes$end[-nrow(contig_genes)] + 1L,
      region_end = contig_genes$start[-1L] - 1L
    ) %>%
      dplyr::filter(.data$region_end >= .data$region_start)
    if (nrow(gap_regions)) {
      regions[[length(regions) + 1L]] <- gap_regions
    }
  }

  dplyr::bind_rows(regions) %>%
    dplyr::distinct(region_start, region_end)
}

.dnmb_generate_context_only_sites <- function(contig, genes, promoter_bp = 150L) {
  contig_genes <- genes %>%
    dplyr::filter(.data$contig == contig) %>%
    dplyr::arrange(.data$start, .data$end)
  .dnmb_generate_context_only_sites_from_genes(contig_genes = contig_genes, promoter_bp = promoter_bp)
}

.dnmb_generate_context_only_sites_from_genes <- function(contig_genes, seq_len = NA_integer_, promoter_bp = 150L) {
  if (!nrow(contig_genes)) {
    return(tibble::tibble())
  }

  promoter_sites <- contig_genes %>%
    dplyr::transmute(
      site_start = dplyr::if_else(.data$strand == "+", pmax(1L, .data$start - 1L), .data$end + 1L),
      site_end = .data$site_start,
      motif_seq = "context_only",
      motif_match_score = 0.3,
      context_class = "promoter",
      affected_gene_id = .data$gene_id,
      affected_gene = .data$gene,
      affected_product = .data$product,
      affected_gene_essentiality_score = .data$essentiality_score,
      affected_gene_essentiality_class = .data$essentiality_class,
      phenotype_categories = vapply(.data$annotation_text, .dnmb_collapse_categories, character(1)),
      phenotype_prediction = vapply(.data$annotation_text, .dnmb_collapse_effects, character(1))
    )

  intergenic_sites <- tibble::tibble()
  if (nrow(contig_genes) >= 2L) {
    gap_starts <- contig_genes$end[-nrow(contig_genes)] + 1L
    gap_ends <- contig_genes$start[-1L] - 1L
    gap_keep <- which(gap_ends >= gap_starts)
    if (length(gap_keep)) {
      left_genes <- contig_genes[-nrow(contig_genes), , drop = FALSE][gap_keep, , drop = FALSE]
      right_genes <- contig_genes[-1L, , drop = FALSE][gap_keep, , drop = FALSE]
      chosen_gene <- ifelse(
        dplyr::coalesce(left_genes$essentiality_score, 0) >= dplyr::coalesce(right_genes$essentiality_score, 0),
        left_genes$gene_id,
        right_genes$gene_id
      )
      chosen_idx <- match(chosen_gene, contig_genes$gene_id)
      intergenic_sites <- tibble::tibble(
        site_start = floor((gap_starts[gap_keep] + gap_ends[gap_keep]) / 2),
        site_end = floor((gap_starts[gap_keep] + gap_ends[gap_keep]) / 2),
        motif_seq = "context_only",
        motif_match_score = 0.2,
        context_class = "intergenic",
        affected_gene_id = contig_genes$gene_id[chosen_idx],
        affected_gene = contig_genes$gene[chosen_idx],
        affected_product = contig_genes$product[chosen_idx],
        affected_gene_essentiality_score = contig_genes$essentiality_score[chosen_idx],
        affected_gene_essentiality_class = contig_genes$essentiality_class[chosen_idx],
        phenotype_categories = vapply(contig_genes$annotation_text[chosen_idx], .dnmb_collapse_categories, character(1)),
        phenotype_prediction = vapply(contig_genes$annotation_text[chosen_idx], .dnmb_collapse_effects, character(1))
      )
    }
  }

  dplyr::bind_rows(promoter_sites, intergenic_sites) %>%
    dplyr::distinct(site_start, .keep_all = TRUE) %>%
    dplyr::mutate(
      essentiality_penalty = vapply(
        seq_len(dplyr::n()),
        function(i) .dnmb_context_essentiality_penalty(.data$context_class[[i]], .data$affected_gene_essentiality_class[[i]]),
        numeric(1)
      ),
      pre_score = pmax(0, pmin(1, .data$motif_match_score + 0.25 - .data$essentiality_penalty))
    ) %>%
    dplyr::arrange(dplyr::desc(.data$pre_score), .data$site_start) %>%
    dplyr::select(-pre_score)
}

.dnmb_iupac_to_regex <- function(motif) {
  if (is.na(motif) || !nzchar(motif)) {
    return("")
  }
  chars <- strsplit(toupper(motif), "", fixed = TRUE)[[1]]
  map <- c(
    A = "A", C = "C", G = "G", T = "T",
    R = "[AG]", Y = "[CT]", S = "[CG]", W = "[AT]", K = "[GT]", M = "[AC]",
    B = "[CGT]", D = "[AGT]", H = "[ACT]", V = "[ACG]", N = "[ACGT]"
  )
  paste(vapply(chars, function(x) dplyr::coalesce(map[[x]], "[ACGT]"), character(1)), collapse = "")
}

.dnmb_classify_interval_context <- function(contig, start, end, genes, promoter_bp = 150L) {
  contig_name <- contig
  contig_genes <- genes %>%
    dplyr::filter(.data$contig == .env$contig_name) %>%
    dplyr::arrange(.data$start, .data$end)
  .dnmb_classify_interval_context_local(
    start = start,
    end = end,
    genes_contig = contig_genes,
    promoter_bp = promoter_bp
  )
}

.dnmb_classify_interval_context_local <- function(start, end, genes_contig, promoter_bp = 150L) {
  gene_hits <- genes_contig %>%
    dplyr::filter(.data$start <= end, .data$end >= start) %>%
    dplyr::arrange(dplyr::desc(.data$essentiality_score), .data$start)
  if (nrow(gene_hits)) {
    top <- gene_hits[1, , drop = FALSE]
    return(list(
      context_class = "gene_body",
      affected_gene_id = top$gene_id[[1]],
      affected_gene = top$gene[[1]],
      affected_product = top$product[[1]],
      affected_gene_essentiality_score = top$essentiality_score[[1]],
      affected_gene_essentiality_class = top$essentiality_class[[1]],
      phenotype_categories = .dnmb_collapse_categories(top$annotation_text[[1]]),
      phenotype_prediction = .dnmb_collapse_effects(top$annotation_text[[1]])
    ))
  }

  promoter_hits <- genes_contig %>%
    dplyr::mutate(
      promoter_start = dplyr::if_else(.data$strand == "+", pmax(1L, .data$start - promoter_bp), .data$end + 1L),
      promoter_end = dplyr::if_else(.data$strand == "+", .data$start - 1L, .data$end + promoter_bp)
    ) %>%
    dplyr::filter(.data$promoter_start <= end, .data$promoter_end >= start) %>%
    dplyr::arrange(dplyr::desc(.data$essentiality_score), .data$start)
  if (nrow(promoter_hits)) {
    top <- promoter_hits[1, , drop = FALSE]
    return(list(
      context_class = "promoter",
      affected_gene_id = top$gene_id[[1]],
      affected_gene = top$gene[[1]],
      affected_product = top$product[[1]],
      affected_gene_essentiality_score = top$essentiality_score[[1]],
      affected_gene_essentiality_class = top$essentiality_class[[1]],
      phenotype_categories = .dnmb_collapse_categories(top$annotation_text[[1]]),
      phenotype_prediction = .dnmb_collapse_effects(top$annotation_text[[1]])
    ))
  }

  nearest <- genes_contig %>%
    dplyr::mutate(
      dist = dplyr::case_when(
        .data$end < start ~ start - .data$end,
        .data$start > end ~ .data$start - end,
        TRUE ~ 0L
      )
    ) %>%
    dplyr::arrange(.data$dist, dplyr::desc(.data$essentiality_score))
  if (nrow(nearest)) {
    top <- nearest[1, , drop = FALSE]
    return(list(
      context_class = "intergenic",
      affected_gene_id = top$gene_id[[1]],
      affected_gene = top$gene[[1]],
      affected_product = top$product[[1]],
      affected_gene_essentiality_score = top$essentiality_score[[1]],
      affected_gene_essentiality_class = top$essentiality_class[[1]],
      phenotype_categories = .dnmb_collapse_categories(top$annotation_text[[1]]),
      phenotype_prediction = .dnmb_collapse_effects(top$annotation_text[[1]])
    ))
  }

  list(
    context_class = "intergenic",
    affected_gene_id = NA_character_,
    affected_gene = NA_character_,
    affected_product = NA_character_,
    affected_gene_essentiality_score = NA_real_,
    affected_gene_essentiality_class = NA_character_,
    phenotype_categories = NA_character_,
    phenotype_prediction = NA_character_
  )
}

.dnmb_context_essentiality_penalty <- function(context_class, essentiality_class) {
  if (is.na(essentiality_class) || !nzchar(essentiality_class)) {
    return(ifelse(context_class == "intergenic", 0.05, 0.15))
  }
  if (context_class == "gene_body") {
    return(dplyr::case_when(
      essentiality_class == "high" ~ 0.95,
      essentiality_class == "medium" ~ 0.55,
      TRUE ~ 0.2
    ))
  }
  if (context_class == "promoter") {
    return(dplyr::case_when(
      essentiality_class == "high" ~ 0.7,
      essentiality_class == "medium" ~ 0.4,
      TRUE ~ 0.15
    ))
  }
  dplyr::case_when(
    essentiality_class == "high" ~ 0.12,
    essentiality_class == "medium" ~ 0.07,
    TRUE ~ 0.03
  )
}

.dnmb_build_variant_catalog <- function(elements, target_sites, scenarios, genes, max_variant_targets_per_element = 50L) {
  max_variant_targets_per_element <- as.integer(max_variant_targets_per_element)
  insertion_rows <- list()
  if (nrow(elements) && nrow(target_sites)) {
    for (i in seq_len(nrow(elements))) {
      family_sites <- target_sites %>%
        dplyr::filter(
          .data$target_family == elements$element_family[[i]],
          !(.data$contig == elements$contig[[i]] & .data$site_end >= elements$start[[i]] & .data$site_start <= elements$end[[i]])
        ) %>%
        dplyr::arrange(dplyr::desc(.data$target_site_score)) %>%
        dplyr::slice_head(n = max_variant_targets_per_element)
      if (!nrow(family_sites)) {
        next
      }

      donor_context <- .dnmb_classify_interval_context(
        contig = elements$contig[[i]],
        start = elements$start[[i]],
        end = elements$end[[i]],
        genes = genes
      )
      source_activity <- max(0.1, min(1,
        0.15 * elements$confidence_score[[i]] +
          0.2 * as.integer(isTRUE(elements$tir_found[[i]])) +
          0.15 * as.integer(!is.na(elements$sequence_element_ids[[i]]) && nzchar(elements$sequence_element_ids[[i]])) +
          0.1 * as.integer(!is.na(elements$reference_hit_id[[i]]) && nzchar(elements$reference_hit_id[[i]]))
      ))
      donor_move_penalty <- if (grepl("cut_and_paste|rolling_circle|peel_and_paste", dplyr::coalesce(family_sites$transposition_modes[[1]], ""))) {
        .dnmb_context_essentiality_penalty(donor_context$context_class, donor_context$affected_gene_essentiality_class) * 0.25
      } else {
        0.02
      }

      rows_i <- family_sites %>%
        dplyr::mutate(
          variant_id = paste0("VAR_INS_", sprintf("%04d", seq_len(dplyr::n())), "_", elements$element_id[[i]]),
          variant_type = "transposition_insertion",
          source_element_id = elements$element_id[[i]],
          source_contig = elements$contig[[i]],
          source_start = elements$start[[i]],
          source_end = elements$end[[i]],
          source_family = elements$element_family[[i]],
          source_activity_score = source_activity,
          donor_context_class = donor_context$context_class,
          donor_essentiality_class = donor_context$affected_gene_essentiality_class,
          likelihood_score = pmax(0, pmin(1, 0.45 * .data$target_site_score + 0.55 * source_activity)),
          viability_score = pmax(0, pmin(1, 0.6 * .data$target_site_score + 0.4 * source_activity - donor_move_penalty)),
          viability_class = dplyr::case_when(
            .data$viability_score >= 0.7 ~ "high",
            .data$viability_score >= 0.35 ~ "medium",
            TRUE ~ "low"
          ),
          affected_gene_count = dplyr::if_else(!is.na(.data$affected_gene_id), 1L, 0L),
          affected_genes = .data$affected_gene_id,
          affected_products = .data$affected_product,
          rationale = stringr::str_squish(paste(
            "Genome-native target-site model:",
            .data$motif_source,
            "context=", .data$context_class,
            "modes=", .data$transposition_modes
          ))
        ) %>%
        dplyr::select(
          variant_id,
          variant_type,
          source_element_id,
          source_contig,
          source_start,
          source_end,
          source_family,
          contig,
          site_start,
          site_end,
          transposition_modes,
          context_class,
          affected_gene_count,
          affected_genes,
          affected_products,
          affected_gene_essentiality_class,
          phenotype_categories,
          phenotype_prediction,
          source_activity_score,
          motif_match_score,
          target_site_score,
          likelihood_score,
          viability_score,
          viability_class,
          rationale
        )
      insertion_rows[[length(insertion_rows) + 1L]] <- rows_i
    }
  }

  excision_rows <- if (nrow(elements)) {
    lapply(seq_len(nrow(elements)), function(i) {
      donor_context <- .dnmb_classify_interval_context(
        contig = elements$contig[[i]],
        start = elements$start[[i]],
        end = elements$end[[i]],
        genes = genes
      )
      source_activity <- max(0.1, min(1,
        0.15 * elements$confidence_score[[i]] +
          0.2 * as.integer(isTRUE(elements$tir_found[[i]]))
      ))
      viability_score <- pmax(0, pmin(1,
        0.45 * source_activity + 0.35 * (1 - .dnmb_context_essentiality_penalty(donor_context$context_class, donor_context$affected_gene_essentiality_class))
      ))
      tibble::tibble(
        variant_id = paste0("VAR_EXC_", elements$element_id[[i]]),
        variant_type = "donor_excision_footprint",
        source_element_id = elements$element_id[[i]],
        source_contig = elements$contig[[i]],
        source_start = elements$start[[i]],
        source_end = elements$end[[i]],
        source_family = elements$element_family[[i]],
        contig = elements$contig[[i]],
        site_start = elements$start[[i]],
        site_end = elements$end[[i]],
        transposition_modes = "donor_site_change",
        context_class = donor_context$context_class,
        affected_gene_count = ifelse(is.na(donor_context$affected_gene_id), 0L, 1L),
        affected_genes = donor_context$affected_gene_id,
        affected_products = donor_context$affected_product,
        affected_gene_essentiality_class = donor_context$affected_gene_essentiality_class,
        phenotype_categories = donor_context$phenotype_categories,
        phenotype_prediction = donor_context$phenotype_prediction,
        source_activity_score = source_activity,
        motif_match_score = NA_real_,
        target_site_score = NA_real_,
        likelihood_score = source_activity,
        viability_score = viability_score,
        viability_class = dplyr::case_when(
          viability_score >= 0.7 ~ "high",
          viability_score >= 0.35 ~ "medium",
          TRUE ~ "low"
        ),
        rationale = "Excision can leave a donor-site footprint or local scar even when transposition target is elsewhere."
      )
    }) %>% dplyr::bind_rows()
  } else {
    tibble::tibble()
  }

  scenario_rows <- if (nrow(scenarios)) {
    scenarios %>%
      dplyr::transmute(
        variant_id = paste0("VAR_", .data$scenario_id),
        variant_type = .data$scenario_type,
        source_element_id = .data$source_elements,
        source_contig = .data$contig,
        source_start = .data$interval_start,
        source_end = .data$interval_end,
        source_family = .data$source_element_families,
        contig = .data$contig,
        site_start = .data$interval_start,
        site_end = .data$interval_end,
        transposition_modes = .data$orientation_class,
        context_class = "interval",
        affected_gene_count = .data$affected_gene_count,
        affected_genes = .data$affected_genes,
        affected_products = .data$affected_products,
        affected_gene_essentiality_class = .data$highest_essentiality_class,
        phenotype_categories = .data$phenotype_categories,
        phenotype_prediction = .data$phenotype_prediction,
        source_activity_score = NA_real_,
        motif_match_score = NA_real_,
        target_site_score = .data$confidence_score,
        likelihood_score = .data$confidence_score,
        viability_score = .data$viability_score,
        viability_class = .data$viability_class,
        rationale = .data$rationale
      )
  } else {
    tibble::tibble()
  }

  dplyr::bind_rows(
    dplyr::bind_rows(insertion_rows),
    excision_rows,
    scenario_rows
  ) %>%
    dplyr::arrange(dplyr::desc(.data$viability_score), dplyr::desc(.data$likelihood_score))
}

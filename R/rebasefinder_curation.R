.dnmb_rebasefinder_flag <- function(text, pattern) {
  text <- base::tolower(base::as.character(text))
  !base::is.na(text) & base::grepl(pattern, text, perl = TRUE)
}

.dnmb_rebasefinder_annotation_flags <- function(text) {
  text <- base::tolower(base::as.character(text))
  typeiv_gmrsd <- (
    .dnmb_rebasefinder_flag(text, "pf03235|duf262|gmrsd[_ -]?n") &
      .dnmb_rebasefinder_flag(text, "pf07510|gmrsd[_ -]?c|\\bgmrsd\\b")
  ) | .dnmb_rebasefinder_flag(text, "\\bgmrsd\\b")
  typeiv_mrr <- .dnmb_rebasefinder_flag(text, "pf04471|mrr[_ -]?cat|\\bmrr family\\b")
  typeiv_resiii_pld <-
    .dnmb_rebasefinder_flag(text, "pf04851|resiii|dexhc_re_i_iii_res") &
    .dnmb_rebasefinder_flag(text, "pf13091|pldc_2|pld[- ]like|phospholipase d/nuclease")
  typeiv_specific <- typeiv_gmrsd | typeiv_mrr | typeiv_resiii_pld |
    .dnmb_rebasefinder_flag(
      text,
      "pvruts1i|pvurts1i|sauusi|mspji|modification-dependent restriction endonuclease"
    )
  dna_mtase <- .dnmb_rebasefinder_flag(
    text,
    "dna[^|;]{0,45}(methylase|methyltransferase)|(?:methylase|methyltransferase)[^|;]{0,45}dna|adenine-specific dna-methyl|cytosine-specific dna-methyl|(?:type[ _-]?(?:i|iii|1|3)\\b|restriction[- ]modification system)[^|;]{0,45}(subunit m|mod subunit|modification methylase)|type[ _-]?ii[gls]?[^|;]{0,55}(restriction enzyme/methyltransferase|restriction-modification enzyme|methyltransferase)|restriction enzyme/methyltransferase|\\bhsdm\\b|hsdm_n|rmtype1_m|rmtype3_mod|n6_n4_mtase|methyltransf[ _-]?d12|pf02086|eco57i|c5_dna_methyl"
  )
  rm_rease <- .dnmb_rebasefinder_flag(
    text,
    "restriction[- ]modification[^|;]{0,45}(endonuclease|restriction enzyme|subunit r|r subunit)|restriction (endonuclease|enzyme)|type[ _-]?(i|ii|iii|iv)[^|;]{0,40}(restriction endonuclease|restriction enzyme|res subunit|subunit r|r subunit)|type[ _-]?ii[gls]?[^|;]{0,40}restriction-modification enzyme|\\bmmei\\b|\\bhsdr\\b|hsdr_n|ecor124i|resiii|rmtype1_r|rmtype3_res|modification-dependent restriction|\\bmrr\\b|\\bmcr[abc]\\b"
  ) | typeiv_specific
  rm_specificity <- .dnmb_rebasefinder_flag(
    text,
    "specificity subunit|target recognition domain|(?:type[ _-]?i|restriction[- ]modification system)[^|;]{0,45}subunit s|\\bhsds\\b|methylase_s|restriction endonuclease subunit s|rmtype1_s"
  )
  non_dna_mtase <- .dnmb_rebasefinder_flag(
    text,
    "(?:rrna|trna|ribosomal rna|rna)[^|;]{0,55}(?:methylase|methyltransferase)|(?:methylase|methyltransferase)[^|;]{0,55}(?:rrna|trna|ribosomal rna)|protein(?:(?!dna)[^|;]){0,25}(?:methylase|methyltransferase)|(?:methylase|methyltransferase)(?:[, -]+protein(?:[- ][a-z]+)?)\\b|methylated-dna[^|;]{0,45}protein-cysteine|o-?6-alkylguanine[^|;]{0,45}(transferase|methyltransferase)|\\b(?:ogt|mgmt)\\b|serine hydroxymethyltransferase|aminomethyltransferase|glycine hydroxymethyltransferase|chemotaxis[^|;]{0,35}methyl|ubiquinone[^|;]{0,35}methyl|menaquinone[^|;]{0,35}methyl|methylthiotransferase|spou(?:_methylase)?|ftsj(?:-like)?|\\b(?:rsm|rlm|trm|prm|erm)[a-z0-9]*\\b|\\b(?:ksga|cfr|ruma|rumb|fmu|rimo|hemk)\\b|lysine methyltransferase|glutamine methyltransferase|\\bche[rb]\\b|\\bubi[eg]\\b"
  )
  direct_repair_mtase <- .dnmb_rebasefinder_flag(
    text,
    "methylated-dna[^|;]{0,45}protein-cysteine|o-?6-alkylguanine[^|;]{0,45}(transferase|methyltransferase)|\\b(?:ogt|mgmt)\\b"
  )
  other_defense <- .dnmb_rebasefinder_flag(
    text,
    "\\bbrex\\b|brex__|\\bpglx\\b|\\bpglxi\\b|\\bbrx[a-z0-9]*\\b|\\bdisarm\\b|\\bdnd[a-z0-9]*\\b|phosphorothioat|\\bdrm[a-z0-9]*\\b|crispr|\\bcas[0-9][a-z0-9]*\\b|abortive infection|anti[- ]phage|\\bgao[_ -]?rl\\b|toxin[- ]antitoxin|\\babieii?\\b|type iv ta system|\\bmqs[ar]\\b|\\brele\\b|\\bhig[ab]\\b|\\bmaz[ef]\\b|\\bvap[bc]\\b"
  )
  repair_nuclease <- .dnmb_rebasefinder_flag(
    text,
    "\\bvsr\\b|very short patch|mismatch repair|endonuclease (iii|iv|v|viii)|ap endonuclease|exonuclease|\\bmut[hlms]\\b|\\buvr[abcd]\\b|\\brecbcd?\\b|\\badd[ab]\\b|dna repair nuclease|homing endonuclease|\\bribonuclease\\b"
  )
  unrelated_helicase <- .dnmb_rebasefinder_flag(
    text,
    "atp-dependent rna helicase|rna helicase|\\bdbpa\\b|\\bhrpa\\b|\\bsrmb\\b|dna helicase recq|\\brecq\\b|\\buvrd\\b|\\bpcr[ab]\\b|\\brep helicase|\\bding\\b|\\bpria\\b|\\bseca\\b|\\bmfd\\b|\\brecg\\b|molecular chaperone dnak|fe-s protein assembly chaperone hsca|\\bdnak\\b|\\bhsca\\b|hsp70"
  )
  toxin_atpase <- .dnmb_rebasefinder_flag(
    text,
    "pf13304|aaa_21|cog1106|k06926|putative abieii toxin|type iv ta system|pf07693|sll1717|kap p-loop"
  )
  unrelated_metabolic <- .dnmb_rebasefinder_flag(
    text,
    "phosphonate[^|;]{0,35}phnh|\\bphnh\\b|serine hydroxymethyltransferase|aminomethyltransferase|glycine hydroxymethyltransferase"
  )
  mobile_nuclease <- .dnmb_rebasefinder_flag(
    text,
    "transposase|insertion sequence|\\bins[a-z0-9]+\\b|integrase|recombinase|reverse gyrase|terminase|\\biscb\\b|\\btnpb\\b|homing endonuclease|anti-restriction|\\barda\\b"
  )
  non_rm_restriction_fold <- .dnmb_rebasefinder_flag(
    text,
    "\\byaeq\\b|pf07152|\\byran\\b|upf0102"
  )
  base::data.frame(
    dna_mtase = dna_mtase,
    rm_rease = rm_rease,
    typeiv_gmrsd = typeiv_gmrsd,
    typeiv_mrr = typeiv_mrr,
    typeiv_specific = typeiv_specific,
    rm_specificity = rm_specificity,
    rm_specific = dna_mtase | rm_rease | rm_specificity,
    explicit_non_dna_mtase = non_dna_mtase | direct_repair_mtase,
    non_dna_mtase = (non_dna_mtase & !dna_mtase) | direct_repair_mtase,
    other_defense = other_defense,
    repair_nuclease = repair_nuclease,
    unrelated_helicase = unrelated_helicase,
    toxin_atpase = toxin_atpase,
    unrelated_metabolic = unrelated_metabolic,
    mobile_nuclease = mobile_nuclease,
    non_rm_restriction_fold = non_rm_restriction_fold,
    stringsAsFactors = FALSE
  )
}

.dnmb_rebasefinder_non_rm_mtase_mask <- function(text) {
  flags <- .dnmb_rebasefinder_annotation_flags(text)
  flags$explicit_non_dna_mtase | flags$other_defense | flags$unrelated_metabolic
}

.dnmb_rebasefinder_non_rm_rease_mask <- function(text) {
  flags <- .dnmb_rebasefinder_annotation_flags(text)
  flags$other_defense | flags$repair_nuclease | flags$unrelated_helicase |
    flags$toxin_atpase | flags$unrelated_metabolic | flags$mobile_nuclease
}

.dnmb_rebasefinder_reference_rows <- function(subject_ids,
                                              rebase_data,
                                              alignment_length = NA_real_,
                                              preferred_family = NA_character_,
                                              preferred_role = NA_character_) {
  subject_ids <- base::as.character(subject_ids)
  alignment_length <- base::rep_len(suppressWarnings(base::as.numeric(alignment_length)), base::length(subject_ids))
  preferred_family <- base::gsub(
    "_", " ",
    base::rep_len(base::as.character(preferred_family), base::length(subject_ids)),
    fixed = TRUE
  )
  preferred_role <- base::rep_len(base::as.character(preferred_role), base::length(subject_ids))
  n <- base::length(subject_ids)
  out <- base::data.frame(
    subject_id = subject_ids,
    reference_index = base::rep(NA_integer_, n),
    reference_name = .dnmb_rebasefinder_clean_rebase_subject(subject_ids),
    reference_aa_len = base::rep(NA_real_, n),
    reference_family = base::rep(NA_character_, n),
    reference_role = base::rep(NA_character_, n),
    reference_enz_type = base::rep(NA_character_, n),
    reference_is_gold_standard = base::rep(FALSE, n),
    reference_is_putative = base::rep(FALSE, n),
    reference_resolved_subject_id = subject_ids,
    stringsAsFactors = FALSE
  )
  if (!base::is.data.frame(rebase_data) || !base::nrow(rebase_data) ||
      !"enzyme_name" %in% base::names(rebase_data)) return(out)

  enzyme_names <- base::as.character(rebase_data$enzyme_name)
  requested_names <- base::unique(out$reference_name[
    !base::is.na(out$reference_name) & base::nzchar(out$reference_name)
  ])
  requested_positions <- base::which(enzyme_names %in% requested_names)
  candidate_map <- if (base::length(requested_positions)) {
    base::split(requested_positions, enzyme_names[requested_positions])
  } else {
    list()
  }

  n_references <- base::nrow(rebase_data)
  seq_len <- base::rep(NA_real_, n_references)
  if ("seq_length" %in% base::names(rebase_data)) {
    seq_len[requested_positions] <- suppressWarnings(
      base::as.numeric(rebase_data$seq_length[requested_positions])
    )
  } else if ("sequence" %in% base::names(rebase_data)) {
    seq_len[requested_positions] <- base::nchar(
      base::as.character(rebase_data$sequence[requested_positions])
    )
  }
  gold <- base::rep(FALSE, n_references)
  if ("is_gold_standard" %in% base::names(rebase_data)) {
    gold[requested_positions] <- rebase_data$is_gold_standard[requested_positions] %in% TRUE
  }
  enz_types <- base::rep(NA_character_, n_references)
  if ("enz_type" %in% base::names(rebase_data)) {
    enz_types[requested_positions] <- base::as.character(rebase_data$enz_type[requested_positions])
  }
  reference_families <- base::rep(NA_character_, n_references)
  if ("rm_type" %in% base::names(rebase_data)) {
    reference_families[requested_positions] <- base::gsub(
      "_", " ",
      base::as.character(rebase_data$rm_type[requested_positions]),
      fixed = TRUE
    )
  }
  reference_roles <- base::rep(NA_character_, n_references)
  if ("subunit" %in% base::names(rebase_data)) {
    reference_roles[requested_positions] <- base::as.character(
      rebase_data$subunit[requested_positions]
    )
  }
  invalid_role <- requested_positions[
    base::is.na(reference_roles[requested_positions]) |
      !reference_roles[requested_positions] %in% c("M", "R", "S", "RM")
  ]
  reference_roles[invalid_role] <- .dnmb_rebasefinder_role_from_hit(enzyme_names[invalid_role])
  infer_role <- requested_positions[
    base::is.na(reference_roles[requested_positions]) &
      !base::is.na(enz_types[requested_positions])
  ]
  if (base::length(infer_role)) {
    reference_roles[infer_role] <- base::vapply(infer_role, function(j) {
      enz_type <- enz_types[[j]]
      if (base::grepl("methyltransferase", enz_type, ignore.case = TRUE)) {
        if (base::grepl("restriction enzyme/methyltransferase|type iig", enz_type, ignore.case = TRUE)) "RM" else "M"
      } else if (base::grepl("specificity", enz_type, ignore.case = TRUE)) {
        "S"
      } else if (base::grepl("restriction", enz_type, ignore.case = TRUE)) {
        "R"
      } else {
        NA_character_
      }
    }, character(1))
  }

  for (i in base::seq_along(subject_ids)) {
    subject <- subject_ids[[i]]
    clean <- out$reference_name[[i]]
    if (base::is.na(subject) || !base::nzchar(subject) || base::is.na(clean)) next
    candidates <- candidate_map[[clean]]
    if (base::is.null(candidates)) candidates <- integer()
    if (!base::length(candidates)) next
    suffix <- if (base::grepl("_[0-9]+$", subject, perl = TRUE)) {
      suppressWarnings(base::as.integer(base::sub("^.*_([0-9]+)$", "\\1", subject, perl = TRUE)))
    } else {
      NA_integer_
    }
    if (!base::is.na(suffix) && suffix %in% candidates) {
      chosen <- suffix
    } else if (base::length(candidates) == 1L) {
      chosen <- candidates[[1]]
    } else {
      family_match <- !base::is.na(preferred_family[[i]]) &
        reference_families[candidates] == preferred_family[[i]]
      family_match[base::is.na(family_match)] <- FALSE
      role_match <- !base::is.na(preferred_role[[i]]) & (
        reference_roles[candidates] == preferred_role[[i]] |
          (preferred_role[[i]] == "RM" & reference_roles[candidates] %in% c("M", "R", "RM")) |
          (preferred_role[[i]] %in% c("M", "R") & reference_roles[candidates] == "RM")
      )
      role_match[base::is.na(role_match)] <- FALSE
      distance <- base::abs(seq_len[candidates] - alignment_length[[i]])
      distance[base::is.na(distance)] <- Inf
      ord <- base::order(!(family_match & role_match), !family_match, !role_match,
                         distance, !gold[candidates], candidates)
      chosen <- candidates[ord[[1]]]
    }
    out$reference_index[[i]] <- chosen
    if (base::length(candidates) > 1L) {
      out$reference_resolved_subject_id[[i]] <- base::paste0(clean, "_", chosen)
    }
    out$reference_aa_len[[i]] <- seq_len[[chosen]]
    out$reference_family[[i]] <- reference_families[[chosen]]
    enz_type <- enz_types[[chosen]]
    out$reference_enz_type[[i]] <- enz_type
    out$reference_role[[i]] <- reference_roles[[chosen]]
    out$reference_is_gold_standard[[i]] <- gold[[chosen]]
    out$reference_is_putative[[i]] <- !base::is.na(enz_type) & base::grepl("putative", enz_type, ignore.case = TRUE)
  }
  out
}

.dnmb_rebasefinder_alignment_quality <- function(identity,
                                                 query_coverage,
                                                 reference_coverage,
                                                 evalue = NA_real_) {
  identity <- suppressWarnings(base::as.numeric(identity))
  query_coverage <- suppressWarnings(base::as.numeric(query_coverage))
  reference_coverage <- suppressWarnings(base::as.numeric(reference_coverage))
  evalue <- suppressWarnings(base::as.numeric(evalue))
  strong <- !base::is.na(identity) & identity >= 0.50 &
    !base::is.na(query_coverage) & query_coverage >= 0.80 &
    !base::is.na(reference_coverage) & reference_coverage >= 0.80 &
    (base::is.na(evalue) | evalue <= 1e-5)
  moderate <- !strong & !base::is.na(identity) & identity >= 0.30 &
    !base::is.na(query_coverage) & query_coverage >= 0.70 &
    !base::is.na(reference_coverage) & reference_coverage >= 0.70 &
    (base::is.na(evalue) | evalue <= 1e-3)
  base::ifelse(strong, "strong", base::ifelse(moderate, "moderate", "weak"))
}

.dnmb_rebasefinder_strict_blast_support <- function(hits) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits) ||
      !base::all(c("blast_alignment_quality", "blast_role_compatible") %in% base::names(hits))) {
    return(base::rep(FALSE, base::nrow(hits)))
  }
  supported <- hits$blast_alignment_quality %in% c("strong", "moderate") &
    hits$blast_role_compatible %in% TRUE
  if ("blast_family_compatible" %in% base::names(hits)) {
    supported <- supported &
      (base::is.na(hits$blast_family_compatible) | hits$blast_family_compatible %in% TRUE)
  }
  supported[base::is.na(supported)] <- FALSE
  supported
}

.dnmb_rebasefinder_consistent_structure_support <- function(hits) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  n <- base::nrow(hits)
  if (!n) return(logical())
  supported <- if ("structure_supported" %in% base::names(hits)) {
    hits$structure_supported %in% TRUE
  } else if ("structure_pass" %in% base::names(hits)) {
    hits$structure_pass %in% TRUE
  } else {
    base::rep(FALSE, n)
  }
  if ("structure_candidate_consistent" %in% base::names(hits)) {
    supported <- supported & hits$structure_candidate_consistent %in% TRUE
  } else if ("structure_pass" %in% base::names(hits)) {
    supported <- base::rep(FALSE, n)
  }
  supported[base::is.na(supported)] <- FALSE
  supported
}

.dnmb_rebasefinder_select_primary_blast <- function(hits,
                                                    genes,
                                                    blast_tbl,
                                                    rebase_data,
                                                    min_identity = 0.10,
                                                    min_length = 50,
                                                    max_evalue = 0.001) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(hits)

  raw_fields <- c(
    hit_label = "raw_hit_label",
    family_id = "raw_family_id",
    enzyme_role = "raw_enzyme_role",
    evidence_mode = "raw_evidence_mode",
    typing_eligible = "raw_typing_eligible",
    blast_identity = "raw_blast_identity",
    blast_evalue = "raw_blast_evalue",
    blast_bitscore = "raw_blast_bitscore",
    blast_length = "raw_blast_length",
    operon_id = "raw_operon_id"
  )
  for (field in base::names(raw_fields)) {
    raw_name <- raw_fields[[field]]
    if (!raw_name %in% base::names(hits)) {
      hits[[raw_name]] <- if (field %in% base::names(hits)) hits[[field]] else NA
    }
  }

  blast_tbl <- .dnmb_rebasefinder_standardize_blast(blast_tbl)
  if (!base::nrow(blast_tbl) || !base::is.data.frame(rebase_data) || !base::nrow(rebase_data)) {
    return(hits)
  }
  keep <- !base::is.na(blast_tbl$pct_identity) & blast_tbl$pct_identity >= min_identity
  if ("length" %in% base::names(blast_tbl)) {
    keep <- keep & !base::is.na(blast_tbl$length) & blast_tbl$length >= min_length
  }
  if ("evalue" %in% base::names(blast_tbl)) {
    keep <- keep & !base::is.na(blast_tbl$evalue) & blast_tbl$evalue <= max_evalue
  }
  blast_tbl <- blast_tbl[keep, , drop = FALSE]
  if (!base::nrow(blast_tbl)) return(hits)

  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  query_id <- .dnmb_module_clean_annotation_key(blast_tbl$query_id)
  gene_idx <- base::match(query_id, genes$locus_tag)
  preferred_hit_idx <- base::match(query_id, .dnmb_module_clean_annotation_key(hits$query))
  preferred_family <- base::rep(NA_character_, base::nrow(blast_tbl))
  preferred_role <- base::rep(NA_character_, base::nrow(blast_tbl))
  preferred_subject <- base::rep(NA_character_, base::nrow(blast_tbl))
  have_preferred_hit <- !base::is.na(preferred_hit_idx)
  if (base::any(have_preferred_hit)) {
    family_col <- if ("raw_family_id" %in% base::names(hits)) "raw_family_id" else "family_id"
    role_col <- if ("raw_enzyme_role" %in% base::names(hits)) "raw_enzyme_role" else "enzyme_role"
    preferred_family[have_preferred_hit] <- base::as.character(
      hits[[family_col]][preferred_hit_idx[have_preferred_hit]]
    )
    preferred_role[have_preferred_hit] <- base::as.character(
      hits[[role_col]][preferred_hit_idx[have_preferred_hit]]
    )
    subject_col <- if ("raw_hit_label" %in% base::names(hits)) "raw_hit_label" else "hit_label"
    preferred_subject[have_preferred_hit] <- base::as.character(
      hits[[subject_col]][preferred_hit_idx[have_preferred_hit]]
    )
  }
  preferred_family <- base::gsub("_", " ", preferred_family, fixed = TRUE)
  query_len <- if ("translation" %in% base::names(genes)) {
    base::nchar(.dnmb_normalize_translation(genes$translation[gene_idx]))
  } else {
    base::rep(NA_real_, base::nrow(blast_tbl))
  }
  aligned_query <- if ("dnmb_aligned_query_length" %in% base::names(blast_tbl)) {
    suppressWarnings(base::as.numeric(blast_tbl$dnmb_aligned_query_length))
  } else if (base::all(c("qstart", "qend") %in% base::names(blast_tbl))) {
    span <- base::abs(suppressWarnings(base::as.numeric(blast_tbl$qend)) -
                       suppressWarnings(base::as.numeric(blast_tbl$qstart))) + 1
    base::ifelse(base::is.na(span), suppressWarnings(base::as.numeric(blast_tbl$length)), span)
  } else {
    suppressWarnings(base::as.numeric(blast_tbl$length))
  }
  aligned_subject <- if ("dnmb_aligned_subject_length" %in% base::names(blast_tbl)) {
    suppressWarnings(base::as.numeric(blast_tbl$dnmb_aligned_subject_length))
  } else if (base::all(c("sstart", "send") %in% base::names(blast_tbl))) {
    span <- base::abs(suppressWarnings(base::as.numeric(blast_tbl$send)) -
                       suppressWarnings(base::as.numeric(blast_tbl$sstart))) + 1
    base::ifelse(base::is.na(span), suppressWarnings(base::as.numeric(blast_tbl$length)), span)
  } else {
    suppressWarnings(base::as.numeric(blast_tbl$length))
  }
  refs <- .dnmb_rebasefinder_reference_rows(
    blast_tbl$rebase_enzyme,
    rebase_data,
    alignment_length = aligned_subject,
    preferred_family = preferred_family,
    preferred_role = preferred_role
  )
  qcov <- base::pmin(1, aligned_query / query_len)
  scov <- base::pmin(1, aligned_subject / refs$reference_aa_len)
  quality <- .dnmb_rebasefinder_alignment_quality(
    blast_tbl$pct_identity,
    qcov,
    scov,
    if ("evalue" %in% base::names(blast_tbl)) blast_tbl$evalue else NA_real_
  )
  reported_reference_family <- refs$reference_family
  same_raw_subject <- !base::is.na(preferred_subject) & base::nzchar(preferred_subject) &
    .dnmb_rebasefinder_clean_rebase_subject(blast_tbl$rebase_enzyme) ==
      .dnmb_rebasefinder_clean_rebase_subject(preferred_subject)
  preferred_role_compatible <- !base::is.na(preferred_role) & (
    refs$reference_role == preferred_role |
      (preferred_role == "RM" & refs$reference_role %in% c("M", "R", "RM")) |
      (preferred_role %in% c("M", "R") & refs$reference_role == "RM")
  )
  preferred_role_compatible[base::is.na(preferred_role_compatible)] <- FALSE
  explicit_reference_suffix <- base::grepl(
    "_[0-9]+$",
    base::as.character(blast_tbl$rebase_enzyme),
    perl = TRUE
  )
  raw_family_override <- quality == "strong" & same_raw_subject & preferred_role_compatible &
    !explicit_reference_suffix &
    preferred_family %in% c("Type I", "Type II", "Type III") &
    (base::is.na(refs$reference_family) | refs$reference_family != preferred_family)
  raw_family_override[base::is.na(raw_family_override)] <- FALSE
  refs$reference_family[raw_family_override] <- preferred_family[raw_family_override]

  annotation_text <- .dnmb_rebasefinder_gene_annotation_text(genes)
  query_annotation <- annotation_text[gene_idx]
  flags <- .dnmb_rebasefinder_annotation_flags(query_annotation)
  role_hint <- base::rep(NA_character_, base::nrow(blast_tbl))
  role_count <- base::as.integer(flags$dna_mtase) + base::as.integer(flags$rm_rease) +
    base::as.integer(flags$rm_specificity)
  role_hint[role_count == 1L & flags$dna_mtase] <- "M"
  role_hint[role_count == 1L & flags$rm_rease] <- "R"
  role_hint[role_count == 1L & flags$rm_specificity] <- "S"
  hit_idx <- base::match(query_id, .dnmb_module_clean_annotation_key(hits$query))
  context_hit <- base::rep(FALSE, base::nrow(blast_tbl))
  if ("evidence_mode" %in% base::names(hits)) {
    context_hit <- context_hit | (!base::is.na(hit_idx) &
      base::grepl("operon_context", base::as.character(hits$evidence_mode[hit_idx]), fixed = TRUE))
  }
  if ("hit_label" %in% base::names(hits)) {
    context_hit <- context_hit | (!base::is.na(hit_idx) &
      base::grepl("_context:", base::as.character(hits$hit_label[hit_idx]), fixed = TRUE))
  }
  for (col in c("typei_operon_supported", "typeii_operon_supported", "typeiii_operon_supported")) {
    if (col %in% base::names(hits)) {
      context_hit <- context_hit | (!base::is.na(hit_idx) & hits[[col]][hit_idx] %in% TRUE)
    }
  }
  context_role <- base::rep(NA_character_, base::nrow(blast_tbl))
  family_hint <- base::rep(NA_character_, base::nrow(blast_tbl))
  family_tokens <- base::cbind(
    `Type I` = .dnmb_rebasefinder_flag(query_annotation, "\\bhsd[mrs]\\b|type[ _-]?i\\b|rmtype1|tigr00348|tigr00497|ecor124i"),
    `Type II` = .dnmb_rebasefinder_flag(query_annotation, "type[ _-]?ii(?:g)?\\b|type[ _-]?2\\b"),
    `Type III` = .dnmb_rebasefinder_flag(query_annotation, "type[ _-]?iii\\b|type[ _-]?3\\b|resiii|modiii|pf04851|ipr006935"),
    `Type IV` = flags$typeiv_specific | .dnmb_rebasefinder_flag(query_annotation, "type[ _-]?iv\\b|type[ _-]?4\\b")
  )
  family_count <- base::rowSums(family_tokens, na.rm = TRUE)
  for (known_family in base::colnames(family_tokens)) {
    family_hint[family_count == 1L & family_tokens[, known_family]] <- known_family
  }
  hit_context_family <- base::rep(NA_character_, base::nrow(hits))
  context_family_cols <- c(
    `Type I` = "typei_operon_supported",
    `Type II` = "typeii_operon_supported",
    `Type III` = "typeiii_operon_supported"
  )
  context_family_matrix <- base::matrix(
    FALSE,
    nrow = base::nrow(hits),
    ncol = base::length(context_family_cols),
    dimnames = list(NULL, base::names(context_family_cols))
  )
  for (known_family in base::names(context_family_cols)) {
    col <- context_family_cols[[known_family]]
    if (col %in% base::names(hits)) context_family_matrix[, known_family] <- hits[[col]] %in% TRUE
  }
  one_supported_family <- base::rowSums(context_family_matrix) == 1L
  if (base::any(one_supported_family)) {
    hit_context_family[one_supported_family] <- base::colnames(context_family_matrix)[
      base::max.col(context_family_matrix[one_supported_family, , drop = FALSE], ties.method = "first")
    ]
  }
  if ("hit_label" %in% base::names(hits)) {
    label <- base::as.character(hits$hit_label)
    label_family <- base::ifelse(
      base::grepl("^typeI_context:", label), "Type I",
      base::ifelse(base::grepl("^typeII_context:", label), "Type II",
                   base::ifelse(base::grepl("^typeIII_context:", label), "Type III", NA_character_))
    )
    use_label_family <- base::is.na(hit_context_family) & !base::is.na(label_family)
    hit_context_family[use_label_family] <- label_family[use_label_family]
  }
  context_family <- base::rep(NA_character_, base::nrow(blast_tbl))
  context_role[context_hit] <- base::as.character(hits$enzyme_role[hit_idx[context_hit]])
  context_family[context_hit] <- hit_context_family[hit_idx[context_hit]]
  missing_context_family <- context_hit & base::is.na(context_family)
  context_family[missing_context_family] <- base::as.character(hits$family_id[hit_idx[missing_context_family]])
  use_context_role <- context_hit & context_role %in% c("M", "R", "S", "RM")
  role_hint[use_context_role] <- context_role[use_context_role]
  role_compatible <- base::is.na(role_hint) | base::is.na(refs$reference_role) |
    refs$reference_role == role_hint |
    (role_hint == "RM" & refs$reference_role %in% c("M", "R", "RM"))
  canonical_family <- c("Type I", "Type II", "Type III", "Type IV")
  use_context_family <- context_hit & context_family %in% canonical_family
  family_hint[use_context_family] <- context_family[use_context_family]
  use_family_hint <- family_hint %in% canonical_family
  family_compatible <- !use_family_hint | refs$reference_family == family_hint
  family_compatible[base::is.na(family_compatible)] <- !use_family_hint[base::is.na(family_compatible)]
  existing_subject <- base::rep(NA_character_, base::nrow(blast_tbl))
  existing_subject[!base::is.na(hit_idx)] <- base::as.character(hits$hit_label[hit_idx[!base::is.na(hit_idx)]])
  matches_existing_subject <- !base::is.na(existing_subject) & base::nzchar(existing_subject) & (
    base::as.character(blast_tbl$rebase_enzyme) == existing_subject |
      refs$reference_name == .dnmb_rebasefinder_clean_rebase_subject(existing_subject)
  )
  exact_existing_subject <- !base::is.na(existing_subject) & base::nzchar(existing_subject) &
    (base::as.character(blast_tbl$rebase_enzyme) == existing_subject |
       refs$reference_resolved_subject_id == existing_subject)
  quality_rank <- base::match(quality, c("weak", "moderate", "strong"))
  min_coverage <- base::pmin(qcov, scov)
  bitscore <- if ("bitscore" %in% base::names(blast_tbl)) {
    suppressWarnings(base::as.numeric(blast_tbl$bitscore))
  } else base::rep(NA_real_, base::nrow(blast_tbl))
  evalue <- if ("evalue" %in% base::names(blast_tbl)) {
    suppressWarnings(base::as.numeric(blast_tbl$evalue))
  } else base::rep(NA_real_, base::nrow(blast_tbl))

  candidates <- base::data.frame(
    query = query_id,
    curated_blast_subject_id = base::as.character(blast_tbl$rebase_enzyme),
    curated_blast_resolved_subject_id = refs$reference_resolved_subject_id,
    curated_blast_reference_index = refs$reference_index,
    curated_blast_match = refs$reference_name,
    curated_blast_family = refs$reference_family,
    curated_blast_role = refs$reference_role,
    curated_blast_identity = suppressWarnings(base::as.numeric(blast_tbl$pct_identity)),
    curated_blast_evalue = evalue,
    curated_blast_bitscore = bitscore,
    curated_blast_length = suppressWarnings(base::as.numeric(blast_tbl$length)),
    curated_blast_hsp_count = if ("dnmb_hsp_count" %in% base::names(blast_tbl)) suppressWarnings(base::as.integer(blast_tbl$dnmb_hsp_count)) else 1L,
    curated_blast_query_coverage = qcov,
    curated_blast_reference_coverage = scov,
    curated_blast_min_coverage = min_coverage,
    curated_blast_alignment_quality = quality,
    curated_blast_role_compatible = role_compatible,
    curated_blast_family_compatible = family_compatible,
    query_aa_len = query_len,
    rebase_reference_aa_len = refs$reference_aa_len,
    blast_query_coverage = qcov,
    blast_reference_coverage = scov,
    blast_min_coverage = min_coverage,
    blast_alignment_quality = quality,
    blast_role_compatible = role_compatible,
    blast_family_compatible = family_compatible,
    blast_matches_existing_subject = matches_existing_subject,
    rebase_reference_enz_type = refs$reference_enz_type,
    rebase_reference_reported_family = reported_reference_family,
    rebase_reference_family_raw_override = raw_family_override,
    rebase_reference_is_gold_standard = refs$reference_is_gold_standard,
    rebase_reference_is_putative = refs$reference_is_putative,
    .exact_existing_subject = exact_existing_subject,
    .quality_rank = quality_rank,
    stringsAsFactors = FALSE
  )
  supported_family_roles <- candidates$blast_alignment_quality %in% c("strong", "moderate") &
    candidates$blast_role_compatible %in% TRUE &
    candidates$curated_blast_family %in% canonical_family &
    candidates$curated_blast_role %in% c("M", "R", "S", "RM")
  family_role_slugs <- c(
    `Type I` = "typei",
    `Type II` = "typeii",
    `Type III` = "typeiii",
    `Type IV` = "typeiv"
  )
  hit_queries <- .dnmb_module_clean_annotation_key(hits$query)
  for (known_family in base::names(family_role_slugs)) {
    col <- base::paste0("blast_supported_", family_role_slugs[[known_family]], "_roles")
    hits[[col]] <- base::rep(NA_character_, base::nrow(hits))
    idx_family <- base::which(supported_family_roles & candidates$curated_blast_family == known_family)
    if (!base::length(idx_family)) next
    role_map <- base::tapply(
      candidates$curated_blast_role[idx_family],
      candidates$query[idx_family],
      function(x) base::paste(base::sort(base::unique(x)), collapse = "/")
    )
    map_idx <- base::match(hit_queries, base::names(role_map))
    found_role <- !base::is.na(map_idx)
    hits[[col]][found_role] <- base::as.character(role_map[map_idx[found_role]])
  }
  raw_candidates <- candidates[candidates$blast_matches_existing_subject %in% TRUE, , drop = FALSE]
  if (base::nrow(raw_candidates)) {
    raw_ord <- base::order(
      raw_candidates$query,
      !raw_candidates$.exact_existing_subject,
      -raw_candidates$.quality_rank,
      -base::ifelse(base::is.na(raw_candidates$blast_min_coverage), -Inf, raw_candidates$blast_min_coverage),
      -base::ifelse(base::is.na(raw_candidates$curated_blast_bitscore), -Inf, raw_candidates$curated_blast_bitscore),
      -base::ifelse(base::is.na(raw_candidates$curated_blast_identity), -Inf, raw_candidates$curated_blast_identity),
      base::ifelse(base::is.na(raw_candidates$curated_blast_evalue), Inf, raw_candidates$curated_blast_evalue)
    )
    raw_candidates <- raw_candidates[raw_ord, , drop = FALSE]
    raw_candidates <- raw_candidates[!base::duplicated(raw_candidates$query), , drop = FALSE]
  }

  ord <- base::order(
    candidates$query,
    !candidates$blast_family_compatible,
    !candidates$blast_role_compatible,
    -candidates$.quality_rank,
    !candidates$blast_matches_existing_subject,
    -base::ifelse(base::is.na(candidates$blast_min_coverage), -Inf, candidates$blast_min_coverage),
    -base::ifelse(base::is.na(candidates$curated_blast_bitscore), -Inf, candidates$curated_blast_bitscore),
    -base::ifelse(base::is.na(candidates$curated_blast_identity), -Inf, candidates$curated_blast_identity),
    base::ifelse(base::is.na(candidates$curated_blast_evalue), Inf, candidates$curated_blast_evalue)
  )
  candidates <- candidates[ord, , drop = FALSE]
  candidates <- candidates[!base::duplicated(candidates$query), , drop = FALSE]
  candidates$.exact_existing_subject <- NULL
  candidates$.quality_rank <- NULL

  idx <- base::match(.dnmb_module_clean_annotation_key(hits$query), candidates$query)
  found <- !base::is.na(idx)
  for (col in base::setdiff(base::names(candidates), "query")) {
    if (!col %in% base::names(hits)) hits[[col]] <- .dnmb_na_vector_like(candidates[[col]], base::nrow(hits))
    hits[[col]][found] <- candidates[[col]][idx[found]]
  }
  promote <- found & hits$curated_blast_role_compatible %in% TRUE &
    hits$curated_blast_family_compatible %in% TRUE &
    hits$curated_blast_alignment_quality %in% c("strong", "moderate")
  hits$curated_blast_promoted <- promote

  # When a weak alternate candidate is not promoted, keep the generic BLAST
  # fields tied to the still-active raw hit. The alternate remains available
  # in the curated_blast_* provenance columns for review.
  raw_idx <- if (base::nrow(raw_candidates)) {
    base::match(.dnmb_module_clean_annotation_key(hits$query), raw_candidates$query)
  } else {
    base::rep(NA_integer_, base::nrow(hits))
  }
  restore_raw <- !promote & !base::is.na(raw_idx)
  raw_field_map <- c(
    blast_identity = "curated_blast_identity",
    blast_evalue = "curated_blast_evalue",
    blast_bitscore = "curated_blast_bitscore",
    blast_length = "curated_blast_length",
    blast_query_coverage = "blast_query_coverage",
    blast_reference_coverage = "blast_reference_coverage",
    blast_min_coverage = "blast_min_coverage",
    blast_alignment_quality = "blast_alignment_quality",
    blast_role_compatible = "blast_role_compatible",
    blast_family_compatible = "blast_family_compatible",
    blast_matches_existing_subject = "blast_matches_existing_subject",
    rebase_reference_aa_len = "rebase_reference_aa_len",
    rebase_reference_enz_type = "rebase_reference_enz_type",
    rebase_reference_reported_family = "rebase_reference_reported_family",
    rebase_reference_family_raw_override = "rebase_reference_family_raw_override",
    rebase_reference_is_gold_standard = "rebase_reference_is_gold_standard",
    rebase_reference_is_putative = "rebase_reference_is_putative"
  )
  if (base::any(restore_raw)) {
    for (target in base::names(raw_field_map)) {
      source <- raw_field_map[[target]]
      hits[[target]][restore_raw] <- raw_candidates[[source]][raw_idx[restore_raw]]
    }
  }
  missing_active <- found & !promote & base::is.na(raw_idx)
  if (base::any(missing_active)) {
    for (field in c("identity", "evalue", "bitscore", "length")) {
      target <- base::paste0("blast_", field)
      source <- base::paste0("raw_blast_", field)
      hits[[target]][missing_active] <- hits[[source]][missing_active]
    }
    for (field in c(
      "blast_query_coverage", "blast_reference_coverage", "blast_min_coverage",
      "blast_alignment_quality", "blast_role_compatible", "blast_family_compatible",
      "rebase_reference_aa_len", "rebase_reference_enz_type",
      "rebase_reference_reported_family"
    )) {
      hits[[field]][missing_active] <- NA
    }
    hits$blast_matches_existing_subject[missing_active] <- FALSE
    hits$rebase_reference_family_raw_override[missing_active] <- FALSE
    hits$rebase_reference_is_gold_standard[missing_active] <- FALSE
    hits$rebase_reference_is_putative[missing_active] <- FALSE
  }
  if (base::any(promote)) {
    hits$hit_label[promote] <- hits$curated_blast_subject_id[promote]
    family_known <- promote & !base::is.na(hits$curated_blast_family) &
      hits$curated_blast_family %in% canonical_family
    hits$family_id[family_known] <- hits$curated_blast_family[family_known]
    role_known <- promote & !base::is.na(hits$curated_blast_role) &
      hits$curated_blast_role %in% c("M", "R", "S", "RM")
    hits$enzyme_role[role_known] <- hits$curated_blast_role[role_known]
    hits$blast_identity[promote] <- hits$curated_blast_identity[promote]
    hits$blast_evalue[promote] <- hits$curated_blast_evalue[promote]
    hits$blast_bitscore[promote] <- hits$curated_blast_bitscore[promote]
    hits$blast_length[promote] <- hits$curated_blast_length[promote]
    hits$evidence_mode[promote] <- base::ifelse(
      hits$blast_alignment_quality[promote] == "strong",
      "full_length_rebase_homology",
      "moderate_full_length_rebase_homology"
    )
  }
  hits
}

.dnmb_rebasefinder_curation_annotation_text <- function(genes) {
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  cols <- base::grep(
    "^(product|gene|gene_synonym|note|Description|Preferred_name|PFAMs|KEGG_ko)$|^EggNOG_(Description|Preferred_name|PFAMs|KEGG_ko)$|^Signature[ ._]?(accession|description)[ ._](CDD|Pfam|NCBIfam|SUPERFAMILY|InterPro|TIGRFAM|Hamap|SMART|PIRSF|Gene3D|FunFam|PANTHER|PRINTS|ProSiteProfiles|ProSitePatterns)$|^DefenseFinder_(system_type|system_subtype|gene_name|profiles_in_system)$|^PADLOC_(system|protein_name|target_description|hmm_accession|hmm_name|all_domains|best_hits)$",
    base::names(genes),
    ignore.case = TRUE,
    value = TRUE,
    perl = TRUE
  )
  if (!base::length(cols)) return(base::rep(NA_character_, base::nrow(genes)))
  base::vapply(base::seq_len(base::nrow(genes)), function(i) {
    values <- base::as.character(genes[i, cols, drop = TRUE])
    values <- values[!base::is.na(values) & base::nzchar(values) & values != "NA"]
    if (!base::length(values)) NA_character_ else base::tolower(base::paste(values, collapse = " | "))
  }, character(1))
}

.dnmb_rebasefinder_join_reasons <- function(...) {
  values <- list(...)
  n <- base::max(base::vapply(values, base::length, integer(1)))
  values <- base::lapply(values, base::rep_len, length.out = n)
  base::vapply(base::seq_len(n), function(i) {
    one <- base::vapply(values, function(x) base::as.character(x[[i]]), character(1))
    one <- base::unique(one[!base::is.na(one) & base::nzchar(one)])
    if (!base::length(one)) NA_character_ else base::paste(one, collapse = "; ")
  }, character(1))
}

.dnmb_rebasefinder_curate_hits <- function(hits, genes) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(hits)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  idx <- base::match(.dnmb_module_clean_annotation_key(hits$query), genes$locus_tag)
  text <- .dnmb_rebasefinder_curation_annotation_text(genes)
  hit_text <- text[idx]
  flags <- .dnmb_rebasefinder_annotation_flags(hit_text)
  product <- if ("product" %in% base::names(genes)) base::as.character(genes$product[idx]) else NA_character_
  product_flags <- .dnmb_rebasefinder_annotation_flags(product)
  defense_cols <- base::grep("^(DefenseFinder|PADLOC)_", base::names(genes), value = TRUE)
  defense_text <- if (base::length(defense_cols)) {
    base::vapply(base::seq_len(base::nrow(genes)), function(i) {
      values <- base::as.character(genes[i, defense_cols, drop = TRUE])
      values <- values[!base::is.na(values) & base::nzchar(values) & values != "NA"]
      if (!base::length(values)) NA_character_ else base::tolower(base::paste(values, collapse = " | "))
    }, character(1))[idx]
  } else {
    base::rep(NA_character_, base::nrow(hits))
  }
  defense_flags <- .dnmb_rebasefinder_annotation_flags(defense_text)
  signature_cols <- base::grep("^Signature[ ._]", base::names(genes), value = TRUE)
  signature_text <- if (base::length(signature_cols)) {
    base::vapply(base::seq_len(base::nrow(genes)), function(i) {
      values <- base::as.character(genes[i, signature_cols, drop = TRUE])
      values <- values[!base::is.na(values) & base::nzchar(values) & values != "NA"]
      if (!base::length(values)) NA_character_ else base::tolower(base::paste(values, collapse = " | "))
    }, character(1))[idx]
  } else {
    base::rep(NA_character_, base::nrow(hits))
  }
  signature_flags <- .dnmb_rebasefinder_annotation_flags(signature_text)
  typeiv_specific_supported <- signature_flags$typeiv_specific | product_flags$typeiv_specific
  typeiv_gmrsd_supported <- signature_flags$typeiv_gmrsd | product_flags$typeiv_gmrsd
  typeiv_mrr_supported <- signature_flags$typeiv_mrr | product_flags$typeiv_mrr

  role <- base::as.character(hits$enzyme_role)
  family <- base::as.character(hits$family_id)
  evidence_role <- base::rep(NA_character_, base::nrow(hits))
  for (column in c("curated_blast_role", "raw_enzyme_role", "supplemental_blast_role")) {
    if (!column %in% base::names(hits)) next
    value <- base::as.character(hits[[column]])
    use <- base::is.na(evidence_role) & value %in% c("R", "RM")
    evidence_role[use] <- value[use]
  }
  lock_restriction_role <- evidence_role %in% "R" & flags$rm_rease
  lock_fused_role <- evidence_role %in% "RM" & flags$rm_rease & flags$dna_mtase
  role[lock_restriction_role] <- "R"
  role[lock_fused_role] <- "RM"
  role_evidence_locked <- lock_restriction_role | lock_fused_role
  product_inferred_role <- base::ifelse(product_flags$dna_mtase & product_flags$rm_rease, "RM",
    base::ifelse(product_flags$rm_specificity, "S",
      base::ifelse(product_flags$dna_mtase, "M",
        base::ifelse(product_flags$rm_rease, "R", NA_character_))))
  blast_promoted <- if ("curated_blast_promoted" %in% base::names(hits)) {
    hits$curated_blast_promoted %in% TRUE
  } else {
    base::rep(FALSE, base::nrow(hits))
  }
  preserve_restriction_role <- role %in% c("R", "RM") &
    product_inferred_role %in% "M" & flags$rm_rease
  override_from_product <- !blast_promoted & !base::is.na(product_inferred_role) &
    !preserve_restriction_role & !role_evidence_locked
  role[override_from_product] <- product_inferred_role[override_from_product]
  missing_role <- base::is.na(role) | !base::nzchar(role)
  inferred_role <- base::ifelse(flags$dna_mtase & flags$rm_rease, "RM",
    base::ifelse(flags$rm_specificity, "S",
      base::ifelse(flags$dna_mtase, "M",
        base::ifelse(flags$rm_rease, "R", NA_character_))))
  role[missing_role & !base::is.na(inferred_role)] <- inferred_role[missing_role & !base::is.na(inferred_role)]
  hits$enzyme_role <- role
  if (!"blast_min_coverage" %in% base::names(hits)) hits$blast_min_coverage <- NA_real_
  if (!"blast_query_coverage" %in% base::names(hits)) hits$blast_query_coverage <- NA_real_
  if (!"blast_reference_coverage" %in% base::names(hits)) hits$blast_reference_coverage <- NA_real_
  if (!"blast_role_compatible" %in% base::names(hits)) hits$blast_role_compatible <- NA
  if (!"blast_family_compatible" %in% base::names(hits)) hits$blast_family_compatible <- NA
  blast_compatible <- hits$blast_role_compatible %in% TRUE &
    (base::is.na(hits$blast_family_compatible) | hits$blast_family_compatible %in% TRUE)
  strong_blast <- if ("blast_alignment_quality" %in% base::names(hits)) {
    hits$blast_alignment_quality == "strong" & blast_compatible
  } else base::rep(FALSE, base::nrow(hits))
  moderate_blast <- if ("blast_alignment_quality" %in% base::names(hits)) {
    hits$blast_alignment_quality == "moderate" & blast_compatible
  } else base::rep(FALSE, base::nrow(hits))
  strong_blast[base::is.na(strong_blast)] <- FALSE
  moderate_blast[base::is.na(moderate_blast)] <- FALSE
  structure_supported <- .dnmb_rebasefinder_consistent_structure_support(hits)
  homology_model_supported <- if ("homology_model_supported" %in% base::names(hits)) {
    hits$homology_model_supported %in% TRUE
  } else {
    base::rep(FALSE, base::nrow(hits))
  }
  operon_supported <- base::rep(FALSE, base::nrow(hits))
  for (col in c("typei_operon_supported", "typeii_operon_supported", "typeiii_operon_supported")) {
    if (col %in% base::names(hits)) operon_supported <- operon_supported | hits[[col]] %in% TRUE
  }
  operon_candidate <- if ("typeii_operon_candidate" %in% base::names(hits)) {
    hits$typeii_operon_candidate %in% TRUE
  } else {
    base::rep(FALSE, base::nrow(hits))
  }
  mtase_motif <- if ("mtase_motif_verified" %in% base::names(hits)) hits$mtase_motif_verified %in% TRUE else FALSE
  raw_rease_component <- if ("rease_operon_component_raw" %in% base::names(hits)) {
    base::as.character(hits$rease_operon_component_raw)
  } else if ("rease_operon_component" %in% base::names(hits)) {
    base::as.character(hits$rease_operon_component)
  } else base::rep(NA_character_, base::nrow(hits))
  context_component <- if ("operon_component" %in% base::names(hits)) {
    base::as.character(hits$operon_component)
  } else base::rep(NA_character_, base::nrow(hits))

  role_annotation <- (role == "M" & flags$dna_mtase) |
    (role == "R" & flags$rm_rease) |
    (role == "S" & flags$rm_specificity) |
    (role == "RM" & flags$dna_mtase & flags$rm_rease)
  role_annotation[base::is.na(role_annotation)] <- FALSE
  explicit_role_annotation <- (role == "M" & product_flags$dna_mtase) |
    (role == "R" & product_flags$rm_rease) |
    (role == "S" & product_flags$rm_specificity)
  explicit_role_annotation[base::is.na(explicit_role_annotation)] <- FALSE
  domain_support <- role_annotation & !explicit_role_annotation
  motif_support <- (role == "M" & mtase_motif & flags$dna_mtase) |
    (role == "R" & raw_rease_component == "R_motor_nuclease" & flags$rm_rease) |
    (role == "RM" & mtase_motif & flags$dna_mtase & flags$rm_rease &
       raw_rease_component %in% c("R_nuclease", "R_motor_nuclease"))
  motif_support[base::is.na(motif_support)] <- FALSE

  annotation_role_conflict <- (role %in% c("M", "RM") & flags$rm_rease & !flags$dna_mtase) |
    (role == "R" & flags$dna_mtase & !flags$rm_rease) |
    (role == "S" & (flags$dna_mtase | flags$rm_rease) & !flags$rm_specificity)
  positive_rm_evidence <- strong_blast | moderate_blast | role_annotation |
    operon_supported | motif_support | structure_supported
  explicit_other_defense <- product_flags$other_defense | defense_flags$other_defense
  other_defense_conflict <- flags$other_defense & !explicit_other_defense & positive_rm_evidence
  other_defense_call <- flags$other_defense & !other_defense_conflict
  explicit_non_dna_product <- product_flags$explicit_non_dna_mtase
  mixed_dna_mtase_support <- product_flags$dna_mtase & signature_flags$dna_mtase
  mixed_dna_mtase_support[base::is.na(mixed_dna_mtase_support)] <- FALSE
  non_dna_annotation_conflict <- flags$explicit_non_dna_mtase &
    mixed_dna_mtase_support & !explicit_non_dna_product
  blast_role_conflict <- "blast_role_compatible" %in% base::names(hits) &
    !base::is.na(hits$blast_role_compatible) & !hits$blast_role_compatible
  family_conflict <- "blast_family_compatible" %in% base::names(hits) &
    !base::is.na(hits$blast_family_compatible) & !hits$blast_family_compatible
  role_conflict <- annotation_role_conflict | blast_role_conflict
  context_resolves_blast_conflict <- operon_supported & role_annotation
  unresolved_role_conflict <- annotation_role_conflict |
    (blast_role_conflict & !context_resolves_blast_conflict)
  evidence_conflict <- annotation_role_conflict | non_dna_annotation_conflict | other_defense_conflict |
    ((blast_role_conflict | family_conflict) & !context_resolves_blast_conflict)
  raw_partial <- if ("partial_status" %in% base::names(hits)) {
    hits$partial_status == "partial_or_short"
  } else base::rep(FALSE, base::nrow(hits))
  raw_partial[base::is.na(raw_partial)] <- FALSE
  explicit_partial <- if ("partial_reason" %in% base::names(hits)) {
    base::grepl("annotation suggests|pseudogene|frameshift|internal stop|contig edge", hits$partial_reason, ignore.case = TRUE)
  } else base::rep(FALSE, base::nrow(hits))
  reference_complete <- strong_blast &
    !base::is.na(hits$blast_query_coverage) & hits$blast_query_coverage >= 0.90 &
    !base::is.na(hits$blast_reference_coverage) & hits$blast_reference_coverage >= 0.90
  partial_rescued_by_reference <- raw_partial & !explicit_partial & reference_complete
  partial <- raw_partial & !partial_rescued_by_reference

  homology_rescue <- strong_blast & role_annotation
  r_candidate <- role %in% c("R", "RM")
  unrelated_helicase_noise <- r_candidate & flags$unrelated_helicase &
    !(homology_rescue | (operon_supported & flags$rm_rease) | structure_supported)
  explicit_repair_product <- product_flags$repair_nuclease
  repair_noise <- flags$repair_nuclease &
    (r_candidate | (!flags$dna_mtase & !role_annotation)) &
    (explicit_repair_product | !(homology_rescue & flags$rm_rease))
  toxin_noise <- r_candidate & flags$toxin_atpase & !(homology_rescue & flags$rm_rease)
  mobile_noise <- r_candidate & flags$mobile_nuclease & !(homology_rescue & flags$rm_rease)
  non_dna_mtase_noise <- explicit_non_dna_product |
    (flags$explicit_non_dna_mtase & !mixed_dna_mtase_support)
  rm_specificity_override <- role %in% "S" &
    product_flags$rm_specificity & flags$rm_specificity
  non_dna_mtase_noise <- non_dna_mtase_noise & !rm_specificity_override
  non_rm_fold_noise <- r_candidate & flags$non_rm_restriction_fold &
    !(homology_rescue | structure_supported)
  generic_no_positive <- !(strong_blast | moderate_blast | role_annotation |
    operon_supported | operon_candidate | motif_support | structure_supported)
  reference_rec_known <- if ("reference_rec_seq" %in% base::names(hits)) {
    .dnmb_rebasefinder_valid_recognition(hits$reference_rec_seq)
  } else base::rep(FALSE, base::nrow(hits))
  gold_standard <- if ("rebase_reference_is_gold_standard" %in% base::names(hits)) {
    hits$rebase_reference_is_gold_standard %in% TRUE
  } else base::rep(FALSE, base::nrow(hits))

  exclusion_reason <- .dnmb_rebasefinder_join_reasons(
    base::ifelse(non_dna_mtase_noise, "non_DNA_methyltransferase", NA_character_),
    base::ifelse(other_defense_call, "non_RM_defense_system", NA_character_),
    base::ifelse(repair_noise, "DNA_repair_or_non_RM_nuclease", NA_character_),
    base::ifelse(unrelated_helicase_noise, "unrelated_helicase_or_chaperone", NA_character_),
    base::ifelse(toxin_noise, "toxin_or_non_RM_ATPase", NA_character_),
    base::ifelse(flags$unrelated_metabolic, "unrelated_metabolic_enzyme", NA_character_),
    base::ifelse(mobile_noise, "mobile_or_homing_nuclease", NA_character_),
    base::ifelse(non_rm_fold_noise, "non_RM_restriction_like_fold", NA_character_),
    base::ifelse(generic_no_positive, "no_positive_RM_evidence", NA_character_)
  )
  hard_noise <- !base::is.na(exclusion_reason)

  axes <- base::vapply(base::seq_len(base::nrow(hits)), function(i) {
    values <- c(
      if (strong_blast[[i]]) "full_length_REBASE" else if (moderate_blast[[i]]) "moderate_REBASE" else NULL,
      if (role_annotation[[i]]) "RM_annotation_or_domain" else NULL,
      if (operon_supported[[i]]) "coherent_operon" else NULL,
      if (motif_support[[i]]) "role_specific_motif" else NULL,
      if (structure_supported[[i]]) "structure" else NULL
    )
    if (!base::length(values)) NA_character_ else base::paste(values, collapse = "+")
  }, character(1))
  independent_axes <- base::vapply(base::strsplit(base::ifelse(base::is.na(axes), "", axes), "+", fixed = TRUE), function(x) {
    base::sum(base::nzchar(x))
  }, integer(1))

  score <- 45 * strong_blast + 30 * moderate_blast + 20 * role_annotation +
    20 * operon_supported + 15 * motif_support + 25 * structure_supported +
    20 * typeiv_gmrsd_supported + 10 * typeiv_mrr_supported +
    5 * reference_rec_known + 5 * gold_standard -
    25 * partial - 35 * evidence_conflict
  score <- base::pmax(0, base::pmin(100, score))
  score[hard_noise] <- 0

  high <- !hard_noise & !evidence_conflict & !partial & (
    (strong_blast & role_annotation) |
      (structure_supported & role_annotation) |
      (role_annotation & operon_supported & motif_support) |
      (role == "R" & family %in% "Type IV" & typeiv_gmrsd_supported)
  )
  high[base::is.na(high)] <- FALSE
  medium <- !high & !hard_noise & !evidence_conflict & !partial & (
    strong_blast |
      (moderate_blast & independent_axes >= 2L) |
      (role_annotation & operon_supported) |
      (role == "M" & flags$dna_mtase & mtase_motif) |
      (role == "R" & family %in% "Type IV" & typeiv_specific_supported)
  )
  medium[base::is.na(medium)] <- FALSE
  tier <- base::ifelse(other_defense_call, "other_defense",
    base::ifelse(hard_noise, "excluded_noise",
    base::ifelse(high, "high", base::ifelse(medium, "medium", "review")))
  )
  keep <- tier %in% c("high", "medium")
  decision <- base::ifelse(tier == "high", "retain_high_confidence",
    base::ifelse(tier == "medium", "retain_probable",
      base::ifelse(tier == "review", "manual_review",
        base::ifelse(tier == "other_defense", "separate_other_defense", "exclude_non_RM_noise"))))

  annotation_class <- base::ifelse(other_defense_call, "other_defense_system",
    base::ifelse(non_dna_mtase_noise, "non_DNA_methyltransferase",
      base::ifelse(repair_noise, "repair_or_non_RM_nuclease",
          base::ifelse(unrelated_helicase_noise, "unrelated_helicase_or_chaperone",
          base::ifelse(toxin_noise, "toxin_or_non_RM_ATPase",
            base::ifelse(mobile_noise, "mobile_or_homing_nuclease",
              base::ifelse(non_rm_fold_noise, "non_RM_restriction_like_fold",
                base::ifelse(role_annotation, "RM_specific", "unresolved"))))))))
  association <- base::ifelse(other_defense_call, "other_defense_system",
    base::ifelse(hard_noise, "non_RM_excluded",
    base::ifelse(role == "M" & flags$dna_mtase & !operon_supported, "DNA_MTase_RM_association_unproven",
      base::ifelse(role == "RM" & flags$dna_mtase & flags$rm_rease &
          (strong_blast | motif_support | structure_supported),
        "fused_RM_supported",
        base::ifelse(family %in% "Type IV" & role == "R" & typeiv_specific_supported,
          "modification_dependent_restriction_supported",
          base::ifelse(operon_supported, "classic_RM_supported", "RM_association_unproven"))))))

  final_role <- role
  final_role[unresolved_role_conflict] <- "ambiguous"
  final_role[hard_noise] <- "non_RM"
  final_role[other_defense_call] <- "other_defense"
  final_component <- base::ifelse(
    final_role == "R" & !base::is.na(context_component),
    context_component,
    base::ifelse(final_role == "R", raw_rease_component, final_role)
  )
  reason <- .dnmb_rebasefinder_join_reasons(
    base::ifelse(strong_blast, "near_full_length_REBASE_match", NA_character_),
    base::ifelse(moderate_blast, "moderate_full_length_REBASE_match", NA_character_),
    base::ifelse(role_annotation, "role_specific_annotation_or_domain", NA_character_),
    base::ifelse(operon_supported, "operon_context", NA_character_),
    base::ifelse(operon_candidate & !operon_supported, "candidate_operon_context", NA_character_),
    base::ifelse(motif_support, "role_specific_motif_support", NA_character_),
    base::ifelse(structure_supported, "role_consistent_structure_support", NA_character_),
    base::ifelse(homology_model_supported, "homology_model_motif_geometry_support", NA_character_),
    base::ifelse(partial, "partial_or_short", NA_character_),
    base::ifelse(partial_rescued_by_reference, "length_warning_resolved_by_full_length_reference", NA_character_),
    base::ifelse(unresolved_role_conflict, "role_conflict", NA_character_),
    base::ifelse(context_resolves_blast_conflict & (blast_role_conflict | family_conflict), "blast_conflict_resolved_by_coherent_operon", NA_character_),
    base::ifelse(family_conflict & !context_resolves_blast_conflict, "family_conflict", NA_character_),
    base::ifelse(other_defense_conflict, "conflicting_other_defense_annotation", NA_character_),
    base::ifelse(non_dna_annotation_conflict, "conflicting_non_DNA_methyltransferase_annotation", NA_character_),
    exclusion_reason
  )

  hits$pre_curation_typing_eligible <- hits$typing_eligible
  if (!"blast_typing_eligible" %in% base::names(hits)) {
    hits$blast_typing_eligible <- if ("raw_typing_eligible" %in% base::names(hits)) {
      hits$raw_typing_eligible
    } else {
      hits$typing_eligible
    }
  }
  hits$annotation_class <- annotation_class
  hits$rm_association_class <- association
  hits$role_conflict <- role_conflict
  hits$unresolved_role_conflict <- unresolved_role_conflict
  hits$family_conflict <- family_conflict
  hits$evidence_conflict <- evidence_conflict
  hits$context_resolves_blast_conflict <- context_resolves_blast_conflict
  hits$other_defense_conflict <- other_defense_conflict
  hits$other_defense_confirmed <- other_defense_call
  hits$non_dna_mtase_conflict <- non_dna_annotation_conflict
  hits$typeiv_specific_supported <- typeiv_specific_supported
  hits$typeiv_gmrsd_supported <- typeiv_gmrsd_supported
  hits$typeiv_mrr_supported <- typeiv_mrr_supported
  hits$partial_rescued_by_reference <- partial_rescued_by_reference
  hits$reference_recognition_metadata <- reference_rec_known
  hits$evidence_axes <- axes
  hits$independent_evidence_axes <- independent_axes
  hits$curation_score <- score
  hits$curation_tier <- tier
  hits$curation_decision <- decision
  hits$curation_reason <- reason
  hits$exclusion_reason <- exclusion_reason
  hits$curation_keep <- keep
  final_family <- family
  final_family[hard_noise] <- "non_RM"
  final_family[other_defense_call] <- "other_defense"
  unclassified_mtase <- keep & final_role %in% c("M", "RM") & flags$dna_mtase &
    (base::is.na(final_family) | !base::nzchar(final_family) | final_family == "unknown")
  final_family[unclassified_mtase] <- "Unclassified DNA MTase"
  hits$final_family <- final_family
  hits$family_id[unclassified_mtase] <- final_family[unclassified_mtase]
  hits$final_role <- final_role
  hits$final_component <- final_component
  hits$typing_eligible <- keep
  tier_rank <- base::match(tier, c("high", "medium", "review", "other_defense", "excluded_noise"))
  hits <- hits[base::order(
    tier_rank,
    -hits$curation_score,
    !homology_model_supported,
    base::ifelse(base::is.na(hits$blast_min_coverage), Inf, -hits$blast_min_coverage),
    base::ifelse(base::is.na(hits$blast_identity), Inf, -hits$blast_identity),
    hits$query
  ), , drop = FALSE]
  base::rownames(hits) <- NULL
  hits
}

.dnmb_rebasefinder_curation_table <- function(hits, genes) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(base::data.frame())
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  gene_cols <- base::intersect(
    c("locus_tag", "product", "gene", "protein_id", "contig", "start", "end", "direction"),
    base::names(genes)
  )
  gene_idx <- base::match(.dnmb_module_clean_annotation_key(hits$query), genes$locus_tag)
  gene_part <- genes[gene_idx, gene_cols, drop = FALSE]
  hit_part <- hits
  base::names(hit_part)[base::names(hit_part) == "query"] <- "locus_tag"
  hit_part <- hit_part[, base::setdiff(base::names(hit_part), gene_cols), drop = FALSE]
  out <- base::cbind(gene_part, hit_part, stringsAsFactors = FALSE)
  out$passed_blast_filter <- out$curation_keep
  out$rm_type <- out$final_family
  out$subunit <- out$final_role
  out$blast_match <- out$hit_label
  if (!"rec_seq" %in% base::names(out)) out$rec_seq <- NA_character_
  front <- c(
    "locus_tag", "product", "gene", "protein_id", "contig", "start", "end", "direction",
    "curation_tier", "curation_score", "curation_keep", "curation_decision", "curation_reason", "exclusion_reason",
    "final_family", "final_role", "final_component", "rm_association_class", "annotation_class",
    "passed_blast_filter", "rm_type", "subunit", "blast_match", "rec_seq", "reference_rec_seq",
    "blast_alignment_quality", "blast_role_compatible", "blast_identity", "blast_query_coverage",
    "blast_reference_coverage", "blast_min_coverage", "rebase_reference_aa_len",
    "evidence_axes", "independent_evidence_axes", "role_conflict", "family_conflict",
    "evidence_conflict", "partial_rescued_by_reference", "reference_recognition_metadata"
  )
  front <- base::intersect(front, base::names(out))
  out[, c(front, base::setdiff(base::names(out), front)), drop = FALSE]
}

.dnmb_rebasefinder_style_workbook_sheet <- function(wb, sheet, tbl) {
  if (!requireNamespace("openxlsx", quietly = TRUE) || !base::ncol(tbl)) {
    return(invisible(NULL))
  }
  columns <- base::names(tbl)
  widths <- base::rep(16, base::length(columns))
  widths[base::grepl("locus_tag$|^query$|protein_id|hit_label|blast_match", columns, ignore.case = TRUE)] <- 22
  widths[base::grepl("^product$|annotation_class", columns, ignore.case = TRUE)] <- 42
  widths[base::grepl("^contig$|organism", columns, ignore.case = TRUE)] <- 28
  widths[base::grepl("support|reason|definition|context_partners|context_sources|recognition_source", columns, ignore.case = TRUE)] <- 44
  widths[base::grepl("_link$|_url$", columns, ignore.case = TRUE)] <- 15
  widths[base::grepl("^(start|end|direction)$|identity|coverage|evalue|bitscore|score$|_score$|_position$|_length$", columns, ignore.case = TRUE)] <- 13
  widths <- base::pmin(48, base::pmax(10, widths))
  openxlsx::setColWidths(wb, sheet, cols = base::seq_along(columns), widths = widths)

  header_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
    halign = "center", valign = "center", wrapText = TRUE,
    border = "Bottom", borderColour = "#B4C6E7"
  )
  openxlsx::addStyle(
    wb, sheet, header_style,
    rows = 1L, cols = base::seq_along(columns), gridExpand = TRUE, stack = TRUE
  )
  openxlsx::setRowHeights(wb, sheet, rows = 1L, heights = 30)

  wrap_cols <- base::which(base::grepl(
    "^product$|support|reason|definition|context_partners|context_sources|recognition_source",
    columns, ignore.case = TRUE
  ))
  if (base::nrow(tbl) && base::length(wrap_cols)) {
    wrap_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
    openxlsx::addStyle(
      wb, sheet, wrap_style,
      rows = 2L:(base::nrow(tbl) + 1L), cols = wrap_cols,
      gridExpand = TRUE, stack = TRUE
    )
  }
  invisible(NULL)
}

.dnmb_rebasefinder_write_curation_workbook <- function(output_dir,
                                                       hits,
                                                       genes,
                                                       raw_defenseviz = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) return(NULL)
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- base::file.path(output_dir, "R-M_REBASE_analysis.xlsx")
  raw_backup <- base::file.path(output_dir, "R-M_REBASE_analysis_raw_DefenseViz.xlsx")
  source_sheets <- if (base::file.exists(path)) {
    tryCatch(openxlsx::getSheetNames(path), error = function(e) character())
  } else {
    character()
  }
  source_is_raw <- base::length(source_sheets) > 0L &&
    !base::any(c("All_Augmented_Evidence", "Scoring_Criteria") %in% source_sheets)
  if (source_is_raw) {
    base::file.copy(path, raw_backup, overwrite = TRUE)
  }
  all_tbl <- .dnmb_rebasefinder_curation_table(hits, genes)
  curated <- all_tbl[all_tbl$curation_keep %in% TRUE, , drop = FALSE]
  high_identity <- curated[curated$blast_alignment_quality %in% "strong", , drop = FALSE]
  standalone_mtase <- curated[curated$rm_association_class %in% "DNA_MTase_RM_association_unproven", , drop = FALSE]
  review <- all_tbl[all_tbl$curation_tier %in% "review", , drop = FALSE]
  other_defense <- all_tbl[all_tbl$curation_tier %in% "other_defense", , drop = FALSE]
  excluded <- all_tbl[all_tbl$curation_tier %in% "excluded_noise", , drop = FALSE]

  summary <- if (base::nrow(curated)) {
    count <- stats::aggregate(
      base::rep(1L, base::nrow(curated)),
      by = list(
        rm_type = curated$final_family,
        subunit_clean = curated$final_role,
        curation_tier = curated$curation_tier,
        association_class = curated$rm_association_class
      ),
      FUN = base::sum
    )
    base::names(count)[[5]] <- "count"
    count[base::order(count$rm_type, count$subunit_clean, count$curation_tier, count$association_class), , drop = FALSE]
  } else {
    base::data.frame(
      rm_type = character(), subunit_clean = character(),
      curation_tier = character(), association_class = character(),
      count = integer(), stringsAsFactors = FALSE
    )
  }
  criteria <- base::data.frame(
    rule = c(
      "strong_REBASE", "moderate_REBASE", "high_tier", "medium_tier",
      "review_tier", "other_defense", "excluded_noise", "recognition_sequence"
    ),
    definition = c(
      "identity >= 50%; query and reference coverage >= 80%; evalue <= 1e-5; compatible role",
      "identity >= 30%; query and reference coverage >= 70%; evalue <= 1e-3; compatible role",
      "near-full REBASE plus role-specific evidence, or equivalent structure/context evidence",
      "near-full homology or at least two independent compatible evidence axes",
      "partial, conflicting, standalone, or context-only candidate requiring validation",
      "BREX, DISARM, Dnd, CRISPR, and toxin-antitoxin evidence retained outside classic R-M calls",
      "explicit non-DNA MTase, repair/mobile nuclease, unrelated helicase/ATPase, or non-RM fold",
      "reference_rec_seq is matched-reference metadata and is never auto-transferred to query rec_seq"
    ),
    stringsAsFactors = FALSE
  )

  wb <- openxlsx::createWorkbook()
  sheets <- list(
    RM_Comprehensive = curated,
    High_Identity_50pct = high_identity,
    Standalone_DNA_MTase = standalone_mtase,
    Review_Candidates = review,
    Other_Defense = other_defense,
    Excluded_Noise = excluded,
    All_Augmented_Evidence = all_tbl,
    Raw_DefenseViz = if (base::is.data.frame(raw_defenseviz)) raw_defenseviz else base::data.frame(),
    Type_Summary = summary,
    Scoring_Criteria = criteria
  )
  for (sheet in base::names(sheets)) {
    openxlsx::addWorksheet(wb, sheet)
    tbl <- sheets[[sheet]]
    openxlsx::writeData(wb, sheet, tbl, withFilter = base::nrow(tbl) > 0L)
    if (base::ncol(tbl)) {
      openxlsx::freezePane(wb, sheet, firstRow = TRUE)
      dnmb_write_hyperlink_columns(wb, sheet, tbl, startRow = 1L, startCol = 1L)
      .dnmb_rebasefinder_style_workbook_sheet(wb, sheet, tbl)
    }
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  list(
    xlsx = base::normalizePath(path, winslash = "/", mustWork = FALSE),
    raw_xlsx = if (base::file.exists(raw_backup)) base::normalizePath(raw_backup, winslash = "/", mustWork = FALSE) else NA_character_,
    curated = curated,
    standalone_mtase = standalone_mtase,
    review = review,
    other_defense = other_defense,
    excluded = excluded,
    all = all_tbl
  )
}

.curation_hit_table <- function(queries, roles, families = "Type II") {
  n <- length(queries)
  data.frame(
    query = queries,
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = rep_len(families, n),
    hit_label = paste0(rep_len(roles, n), ".Test"),
    enzyme_role = roles,
    evidence_mode = "high_confidence",
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = TRUE,
    blast_identity = 0.8,
    blast_evalue = 1e-30,
    blast_bitscore = 250,
    blast_length = 250,
    blast_alignment_quality = "strong",
    blast_role_compatible = TRUE,
    blast_min_coverage = 0.9,
    blast_query_coverage = 0.9,
    blast_reference_coverage = 0.9,
    partial_status = "full_length",
    mtase_motif_verified = FALSE,
    rease_operon_component_raw = NA_character_,
    typei_operon_supported = FALSE,
    typeii_operon_supported = FALSE,
    typeiii_operon_supported = FALSE,
    reference_rec_seq = NA_character_,
    stringsAsFactors = FALSE
  )
}

test_that("curation hard filters non-RM methylases and generic helicases", {
  ids <- c(
    "true_m", "rrna_m", "glyA", "rna_hel", "true_r",
    "hypothetical_r", "fragment_m", "brex_m", "context_r"
  )
  genes <- data.frame(
    locus_tag = ids,
    product = c(
      "DNA adenine methylase",
      "16S rRNA methyltransferase RsmD",
      "serine hydroxymethyltransferase",
      "ATP-dependent RNA helicase DbpA",
      "type I restriction endonuclease subunit R",
      "hypothetical protein",
      "Vsr endonuclease domain-containing protein",
      "BREX PglX adenine-specific DNA methyltransferase",
      "P-loop NTPase"
    ),
    EggNOG_PFAMs = c(NA, "N6_N4_Mtase", rep(NA, 7)),
    translation = rep(strrep("A", 300), length(ids)),
    stringsAsFactors = FALSE
  )
  roles <- c("M", "M", "M", "R", "R", "R", "M", "M", "R")
  hits <- .curation_hit_table(ids, roles, families = c(
    "Type II", "Type II", "Type II", "Type III", "Type I",
    "Type II", "Type I", "Type II", "Type II"
  ))
  hits$mtase_motif_verified[hits$query == "true_m"] <- TRUE
  hits$blast_alignment_quality[hits$query %in% c("fragment_m", "context_r")] <- "weak"
  hits$blast_min_coverage[hits$query == "fragment_m"] <- 0.19
  hits$partial_status[hits$query == "fragment_m"] <- "partial_or_short"
  hits$typeii_operon_supported[hits$query == "context_r"] <- TRUE
  hits$rease_operon_component_raw[hits$query == "context_r"] <- "R_motor"

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)
  rownames(out) <- out$query

  expect_identical(out["true_m", "curation_tier"], "high")
  expect_identical(out["true_r", "curation_tier"], "high")
  expect_identical(out["hypothetical_r", "curation_tier"], "medium")
  expect_identical(out["context_r", "curation_tier"], "review")
  expect_identical(out["rrna_m", "exclusion_reason"], "non_DNA_methyltransferase")
  expect_match(out["glyA", "exclusion_reason"], "unrelated_metabolic_enzyme")
  expect_match(out["rna_hel", "exclusion_reason"], "unrelated_helicase")
  expect_match(out["fragment_m", "exclusion_reason"], "DNA_repair")
  expect_match(out["brex_m", "exclusion_reason"], "non_RM_defense_system")
  expect_identical(out["fragment_m", "final_family"], "non_RM")
  expect_identical(out["fragment_m", "rm_association_class"], "non_RM_excluded")
  expect_identical(out["brex_m", "final_family"], "other_defense")
  expect_false(any(out[c("rrna_m", "glyA", "rna_hel", "fragment_m", "brex_m"), "curation_keep"]))
})

test_that("annotation-only R candidates without a family are assigned to review", {
  genes <- data.frame(
    locus_tag = "family_unknown_r",
    product = "restriction endonuclease-like protein",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("family_unknown_r", "R", NA_character_)
  hits$blast_alignment_quality <- NA_character_
  hits$blast_role_compatible <- NA
  hits$blast_family_compatible <- NA

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$curation_tier, "review")
  expect_identical(out$curation_decision, "manual_review")
  expect_identical(out$rm_association_class, "RM_association_unproven")
  expect_false(out$curation_keep)
})

test_that("an explicit Vsr repair nuclease is not rescued by REBASE homology", {
  genes <- data.frame(
    locus_tag = "vsr_repair",
    product = "very short patch repair endonuclease",
    translation = strrep("A", 180),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("vsr_repair", "R", "Type II")

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$curation_tier, "excluded_noise")
  expect_match(out$exclusion_reason, "DNA_repair_or_non_RM_nuclease")
  expect_identical(out$final_family, "non_RM")
  expect_false(out$curation_keep)
})

test_that("an explicit RNA MTase family conflict cannot remain in a retained tier", {
  genes <- data.frame(
    locus_tag = "mixed_c5_mtase",
    product = "DNA cytosine methyltransferase",
    `Signature description_Pfam` = "C-5 cytosine-specific DNA methylase",
    `Signature description_PANTHER` = "tRNA (cytosine(38)-C(5))-methyltransferase",
    translation = strrep("A", 420),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- .curation_hit_table("mixed_c5_mtase", "M", "Type II")

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_true(out$evidence_conflict)
  expect_identical(out$curation_tier, "review")
  expect_false(out$curation_keep)
  expect_match(out$curation_reason, "non_DNA_methyltransferase")
})

test_that("a generic HemK or PrmA call is excluded despite a stray DNA MTase motif annotation", {
  genes <- data.frame(
    locus_tag = "hemk_prma",
    product = "methyltransferase",
    EggNOG_Preferred_name = "hemK1",
    EggNOG_PFAMs = "MTS,PrmA",
    `Signature description_Pfam` = "Ribosomal protein L11 methyltransferase (PrmA)",
    `Signature description_PANTHER` = "HEMK METHYLTRANSFERASE",
    `Signature description_ProSitePatterns` = "N-6 Adenine-specific DNA methylases signature",
    translation = strrep("A", 380),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- .curation_hit_table("hemk_prma", "M", "Type II")
  hits$blast_alignment_quality <- NA_character_
  hits$blast_role_compatible <- NA

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$curation_tier, "excluded_noise")
  expect_false(out$curation_keep)
  expect_match(out$exclusion_reason, "non_DNA_methyltransferase")
})

test_that("an explicit RNA MTase product is excluded even when the selected role is R", {
  genes <- data.frame(
    locus_tag = "rsm_wrong_role",
    product = "16S rRNA methyltransferase RsmD",
    translation = strrep("A", 220),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("rsm_wrong_role", "R", "Type II")

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$curation_tier, "excluded_noise")
  expect_false(out$curation_keep)
  expect_match(out$exclusion_reason, "non_DNA_methyltransferase")
  expect_identical(out$final_family, "non_RM")
})

test_that("signature-level RNA MTase evidence cannot bypass filtering through an R role", {
  genes <- data.frame(
    locus_tag = c("mixed_r", "rsm_signature_r"),
    product = c("DNA cytosine methyltransferase", "hypothetical protein"),
    `Signature description_Pfam` = c(
      "C-5 cytosine-specific DNA methylase; restriction endonuclease domain",
      "restriction endonuclease domain"
    ),
    `Signature description_PANTHER` = c(
      "tRNA (cytosine(38)-C(5))-methyltransferase",
      "16S rRNA methyltransferase RsmD"
    ),
    translation = strrep("A", 300),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- .curation_hit_table(c("mixed_r", "rsm_signature_r"), c("R", "R"), "Type II")

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)
  rownames(out) <- out$query

  expect_true(out["mixed_r", "evidence_conflict"])
  expect_identical(out["mixed_r", "curation_tier"], "review")
  expect_false(out["mixed_r", "curation_keep"])
  expect_identical(out["rsm_signature_r", "curation_tier"], "excluded_noise")
  expect_match(out["rsm_signature_r", "exclusion_reason"], "non_DNA_methyltransferase")
})

test_that("coverage-aware primary selection beats a short high-identity fragment", {
  genes <- data.frame(
    locus_tag = "q1",
    product = "hypothetical protein",
    translation = strrep("A", 125),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("q1", "M", "Type I")
  hits$hit_label <- "M.Fragment"
  rebase <- data.frame(
    enzyme_name = c("M.Fragment", "R.Full"),
    enz_type = c("putative Type I methyltransferase", "Type II restriction enzyme"),
    rec_seq = "",
    org_name = "",
    sequence = c(strrep("A", 657), strrep("A", 125)),
    seq_length = c(657, 125),
    is_gold_standard = c(FALSE, TRUE),
    rm_type = c("Type_I", "Type_II"),
    subunit = c("M", "R"),
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = c("q1", "q1"),
    rebase_enzyme = c("M.Fragment", "R.Full"),
    pct_identity = c(0.69, 0.55),
    length = c(123, 120),
    evalue = c(1e-40, 1e-35),
    bitscore = c(180, 170),
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_identical(out$raw_hit_label, "M.Fragment")
  expect_identical(out$hit_label, "R.Full")
  expect_identical(out$enzyme_role, "R")
  expect_identical(out$blast_alignment_quality, "strong")
  expect_equal(out$blast_reference_coverage, 0.96, tolerance = 0.01)
})

test_that("an unpromoted alternate keeps generic BLAST fields on the active raw hit", {
  genes <- data.frame(
    locus_tag = "short_conflict",
    product = "type I restriction endonuclease subunit R",
    translation = strrep("A", 125),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("short_conflict", "M", "Type I")
  hits$hit_label <- "M.RawFragment"
  hits$blast_identity <- 0.69
  hits$blast_evalue <- 1e-40
  hits$blast_bitscore <- 180
  hits$blast_length <- 123
  rebase <- data.frame(
    enzyme_name = c("M.RawFragment", "R.AlternateFragment"),
    enz_type = c("putative Type I methyltransferase", "putative Type I restriction enzyme"),
    rec_seq = "",
    org_name = "",
    sequence = c(strrep("A", 657), strrep("A", 1231)),
    seq_length = c(657, 1231),
    is_gold_standard = FALSE,
    rm_type = "Type_I",
    subunit = c("M", "R"),
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = c("short_conflict", "short_conflict"),
    rebase_enzyme = c("M.RawFragment", "R.AlternateFragment"),
    pct_identity = c(0.69, 0.66),
    length = c(123, 122),
    evalue = c(1e-40, 1e-38),
    bitscore = c(180, 175),
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_false(out$curated_blast_promoted)
  expect_identical(out$hit_label, "M.RawFragment")
  expect_identical(out$curated_blast_match, "R.AlternateFragment")
  expect_equal(out$curated_blast_identity, 0.66)
  expect_equal(out$curated_blast_reference_coverage, 122 / 1231)
  expect_equal(out$blast_identity, 0.69)
  expect_equal(out$blast_reference_coverage, 123 / 657)
  expect_identical(out$rebase_reference_aa_len, 657)
  expect_false(out$blast_role_compatible)
})

test_that("an unpromoted alternate cannot populate generic fields without an active-subject row", {
  genes <- data.frame(
    locus_tag = "missing_raw_subject",
    product = "type I restriction endonuclease subunit R",
    translation = strrep("A", 125),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("missing_raw_subject", "M", "Type I")
  hits$hit_label <- "M.FilteredRaw"
  hits$blast_identity <- 0.69
  hits$blast_evalue <- 1e-40
  hits$blast_bitscore <- 180
  hits$blast_length <- 40
  rebase <- data.frame(
    enzyme_name = c("M.FilteredRaw", "R.AlternateFragment"),
    enz_type = c("putative Type I methyltransferase", "putative Type I restriction enzyme"),
    rec_seq = "",
    org_name = "",
    sequence = c(strrep("A", 657), strrep("A", 1231)),
    seq_length = c(657, 1231),
    is_gold_standard = FALSE,
    rm_type = "Type_I",
    subunit = c("M", "R"),
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = "missing_raw_subject",
    rebase_enzyme = "R.AlternateFragment",
    pct_identity = 0.66,
    length = 122,
    evalue = 1e-38,
    bitscore = 175,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_false(out$curated_blast_promoted)
  expect_identical(out$hit_label, "M.FilteredRaw")
  expect_identical(out$curated_blast_match, "R.AlternateFragment")
  expect_equal(out$blast_identity, 0.69)
  expect_true(is.na(out$blast_reference_coverage))
  expect_true(is.na(out$blast_alignment_quality))
  expect_true(is.na(out$blast_role_compatible))
  expect_false(out$rebase_reference_is_gold_standard)
})

test_that("a single supported Type III context overrides the previous family during BLAST selection", {
  genes <- data.frame(
    locus_tag = "mod1",
    product = "restriction-modification system Mod subunit",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("mod1", "M", "Type II")
  hits$typeiii_operon_supported <- TRUE
  rebase <- data.frame(
    enzyme_name = c("M.TypeII", "M.TypeIII"),
    enz_type = c("Type II methyltransferase", "Type III methyltransferase"),
    rec_seq = "",
    org_name = "",
    sequence = c(strrep("A", 300), strrep("A", 300)),
    seq_length = 300,
    is_gold_standard = TRUE,
    rm_type = c("Type_II", "Type_III"),
    subunit = "M",
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = "mod1",
    rebase_enzyme = c("M.TypeII", "M.TypeIII"),
    pct_identity = c(0.90, 0.80),
    length = 300,
    evalue = c(1e-80, 1e-70),
    bitscore = c(500, 450),
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_identical(out$family_id, "Type III")
  expect_identical(out$hit_label, "M.TypeIII")
  expect_identical(out$blast_supported_typei_roles, NA_character_)
  expect_identical(out$blast_supported_typeii_roles, "M")
  expect_identical(out$blast_supported_typeiii_roles, "M")
})

test_that("unsuffixed duplicate REBASE names use the raw family and role as a tie-breaker", {
  genes <- data.frame(
    locus_tag = "dup_mod",
    product = "site-specific DNA methyltransferase",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("dup_mod", "M", "Type III")
  hits$hit_label <- "M.Duplicate"
  rebase <- data.frame(
    enzyme_name = c("M.Duplicate", "M.Duplicate"),
    enz_type = c("Type II methyltransferase", "Type III methyltransferase"),
    rec_seq = c("GATC", "CCGG"),
    org_name = "",
    sequence = c(strrep("A", 300), strrep("A", 300)),
    seq_length = 300,
    is_gold_standard = TRUE,
    rm_type = c("Type_II", "Type_III"),
    subunit = "M",
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = "dup_mod",
    rebase_enzyme = "M.Duplicate",
    pct_identity = 0.85,
    length = 300,
    evalue = 1e-80,
    bitscore = 500,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_identical(out$family_id, "Type III")
  expect_identical(out$curated_blast_reference_index, 2L)
  expect_identical(out$curated_blast_resolved_subject_id, "M.Duplicate_2")
  expect_identical(out$blast_supported_typeiii_roles, "M")
  expect_true(is.na(out$blast_supported_typeii_roles))
})

test_that("strong exact-subject homology preserves a compatible raw Type III family", {
  genes <- data.frame(
    locus_tag = "mod_exact",
    product = "site-specific DNA methyltransferase",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("mod_exact", "M", "Type III")
  hits$hit_label <- "M.ReferenceConflict"
  rebase <- data.frame(
    enzyme_name = "M.ReferenceConflict",
    enz_type = "Type II methyltransferase",
    rec_seq = "GATC",
    org_name = "",
    sequence = strrep("A", 300),
    seq_length = 300,
    is_gold_standard = TRUE,
    rm_type = "Type_II",
    subunit = "M",
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = "mod_exact",
    rebase_enzyme = "M.ReferenceConflict",
    pct_identity = 1,
    length = 300,
    evalue = 0,
    bitscore = 600,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_identical(out$family_id, "Type III")
  expect_identical(out$rebase_reference_reported_family, "Type II")
  expect_true(out$rebase_reference_family_raw_override)
  expect_identical(out$blast_supported_typeiii_roles, "M")
})

test_that("an explicitly suffixed REBASE subject keeps its exact reported family", {
  genes <- data.frame(
    locus_tag = "mod_exact_duplicate",
    product = "site-specific DNA methyltransferase",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("mod_exact_duplicate", "M", "Type II")
  hits$hit_label <- "M.Duplicate_2"
  rebase <- data.frame(
    enzyme_name = c("M.Duplicate", "M.Duplicate"),
    enz_type = c("Type II methyltransferase", "Type III methyltransferase"),
    rec_seq = c("GATC", "CCGG"),
    org_name = "",
    sequence = c(strrep("A", 300), strrep("A", 300)),
    seq_length = 300,
    is_gold_standard = TRUE,
    rm_type = c("Type_II", "Type_III"),
    subunit = "M",
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = "mod_exact_duplicate",
    rebase_enzyme = "M.Duplicate_2",
    pct_identity = 0.95,
    length = 300,
    evalue = 0,
    bitscore = 600,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_identical(out$curated_blast_resolved_subject_id, "M.Duplicate_2")
  expect_identical(out$curated_blast_family, "Type III")
  expect_identical(out$family_id, "Type III")
  expect_false(out$rebase_reference_family_raw_override)
})

test_that("protein wording does not hide Eco57I MTases and R-subunit fragments", {
  flags <- DNMB:::.dnmb_rebasefinder_annotation_flags(c(
    "Eco57I restriction-modification methylase domain-containing protein",
    "protein arginine methyltransferase",
    "type I restriction-modification enzyme R subunit C-terminal domain-containing protein"
  ))

  expect_true(flags$dna_mtase[[1]])
  expect_false(flags$explicit_non_dna_mtase[[1]])
  expect_true(flags$explicit_non_dna_mtase[[2]])
  expect_true(flags$rm_rease[[3]])
})

test_that("protein followed by an explicit DNA methylase is not non-DNA evidence", {
  flags <- DNMB:::.dnmb_rebasefinder_annotation_flags(c(
    "Vaccinia Virus protein VP39, DNA methylase, subunit A, domain 2",
    "protein arginine methyltransferase"
  ))

  expect_true(flags$dna_mtase[[1]])
  expect_false(flags$explicit_non_dna_mtase[[1]])
  expect_true(flags$explicit_non_dna_mtase[[2]])
})

test_that("HsdS methylase specificity wording is not treated as a protein methylase", {
  genes <- data.frame(
    locus_tag = "hsdS",
    product = "restriction endonuclease subunit S",
    `Signature description_Gene3D` = "Bipartite methylase S protein, DNA methylase specificity domains",
    `Signature description_Pfam` = "Type I restriction modification DNA specificity domain",
    EggNOG_Preferred_name = "hsdS",
    EggNOG_PFAMs = "Methylase_S",
    translation = strrep("A", 420),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- .curation_hit_table("hsdS", "S", "Type I")

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_false(out$curation_tier == "excluded_noise")
  expect_false(grepl("non_DNA_methyltransferase", out$exclusion_reason))
  expect_identical(out$final_role, "S")
})

test_that("ribonuclease filtering does not match deoxyribonuclease wording", {
  flags <- DNMB:::.dnmb_rebasefinder_annotation_flags(c(
    "type I restriction subunit R deoxyribonuclease",
    "ribonuclease R"
  ))

  expect_false(flags$repair_nuclease[[1]])
  expect_true(flags$repair_nuclease[[2]])
})

test_that("a FokI nuclease domain alone is not Type IV-specific evidence", {
  flags <- DNMB:::.dnmb_rebasefinder_annotation_flags(c(
    "FokI_C family restriction endonuclease domain-containing protein",
    "MspJI family modification-dependent restriction endonuclease"
  ))

  expect_false(flags$typeiv_specific[[1]])
  expect_true(flags$typeiv_specific[[2]])
})

test_that("ResIII plus PLD identifies a Type IV restriction helicase candidate", {
  genes <- data.frame(
    locus_tag = c("existing", "resiii_pld", "rna_helicase"),
    contig = "ctg1",
    start = c(1, 1001, 5001),
    end = c(900, 4200, 6200),
    direction = "+",
    product = c("DNA methyltransferase", "DEAD/DEAH box helicase", "ATP-dependent RNA helicase"),
    `Signature accession_Pfam` = c("PF00145", "PF04851;PF13091", "PF00270"),
    `Signature description_Pfam` = c(
      "C-5 cytosine-specific DNA methylase",
      "Type III restriction enzyme, res subunit; PLD-like domain",
      "DEAD box helicase"
    ),
    translation = strrep("A", 500),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- .curation_hit_table("existing", "M", "Type II")

  out <- DNMB:::.dnmb_rebasefinder_add_typeiv_candidates(hits, genes)

  expect_true("resiii_pld" %in% out$query)
  expect_identical(out$family_id[out$query == "resiii_pld"], "Type IV")
  expect_false("rna_helicase" %in% out$query)
})

test_that("non-overlapping HSPs are combined for bidirectional coverage", {
  genes <- data.frame(
    locus_tag = "multi",
    product = "hypothetical protein",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("multi", "R", "Type II")
  hits$hit_label <- "R.Multi"
  rebase <- data.frame(
    enzyme_name = "R.Multi",
    enz_type = "Type II restriction enzyme",
    rec_seq = "",
    org_name = "",
    sequence = strrep("A", 300),
    seq_length = 300,
    is_gold_standard = TRUE,
    rm_type = "Type_II",
    subunit = "R",
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = c("multi", "multi"),
    rebase_enzyme = c("R.Multi", "R.Multi"),
    pct_identity = c(0.60, 0.60),
    length = c(135, 135),
    qstart = c(1, 166),
    qend = c(135, 300),
    sstart = c(1, 166),
    send = c(135, 300),
    evalue = c(1e-20, 1e-18),
    bitscore = c(120, 115),
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_identical(out$curated_blast_hsp_count, 2L)
  expect_equal(out$blast_query_coverage, 0.90, tolerance = 0.001)
  expect_equal(out$blast_reference_coverage, 0.90, tolerance = 0.001)
  expect_identical(out$blast_alignment_quality, "strong")
})

test_that("numbered REBASE prefixes retain M R and S roles", {
  expect_identical(
    DNMB:::.dnmb_rebasefinder_role_from_hit(c("M12.Test", "R3.Test", "S20.Test", "NoPrefix")),
    c("M", "R", "S", NA_character_)
  )
})

test_that("Type I subunit M annotations are not misclassified as restriction subunits", {
  flags <- DNMB:::.dnmb_rebasefinder_annotation_flags(
    "type I restriction-modification system subunit M HsdM"
  )

  expect_true(flags$dna_mtase)
  expect_false(flags$rm_rease)

  genes <- data.frame(
    locus_tag = "hsdM",
    product = "type I restriction-modification system subunit M HsdM",
    translation = strrep("A", 500),
    stringsAsFactors = FALSE
  )
  out <- DNMB:::.dnmb_rebasefinder_curate_hits(
    .curation_hit_table("hsdM", "M", "Type I"),
    genes
  )
  expect_identical(out$curation_tier, "high")
})

test_that("Type IIG fused restriction methyltransferases retain the RM role", {
  genes <- data.frame(
    locus_tag = "iig",
    product = "Type IIG restriction enzyme/methyltransferase",
    translation = strrep("A", 700),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("iig", "RM", "Type II")

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$final_role, "RM")
  expect_identical(out$curation_tier, "high")
  expect_identical(out$rm_association_class, "fused_RM_supported")
})

test_that("MmeI domain evidence supports a moderate Type IIL RM fusion", {
  genes <- data.frame(
    locus_tag = "mme_iil",
    product = "class I SAM-dependent DNA methyltransferase",
    `Signature description_Pfam` = paste(
      "MmeI, target recognition domain",
      "MmeI, DNA-methyltransferase domain",
      "MmeI, helicase spacer domain",
      sep = "; "
    ),
    translation = strrep("A", 735),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- .curation_hit_table("mme_iil", "RM", "Type II")
  hits$blast_alignment_quality <- "moderate"
  hits$blast_min_coverage <- 0.73
  hits$blast_query_coverage <- 0.76
  hits$blast_reference_coverage <- 0.73

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$final_role, "RM")
  expect_identical(out$curation_tier, "medium")
  expect_identical(out$annotation_class, "RM_specific")
})

test_that("consistent structure support confirms a fused Type IIL RM association", {
  genes <- data.frame(
    locus_tag = "mme_structure",
    product = "class I SAM-dependent DNA methyltransferase",
    `Signature description_Pfam` = paste(
      "MmeI, target recognition domain",
      "MmeI, DNA-methyltransferase domain",
      "MmeI, helicase spacer domain",
      sep = "; "
    ),
    translation = strrep("A", 735),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- .curation_hit_table("mme_structure", "RM", "Type II")
  hits$blast_alignment_quality <- "moderate"
  hits$blast_min_coverage <- 0.73
  hits$blast_query_coverage <- 0.76
  hits$blast_reference_coverage <- 0.73
  hits$structure_supported <- TRUE
  hits$structure_candidate_consistent <- TRUE

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$final_role, "RM")
  expect_identical(out$curation_tier, "high")
  expect_identical(out$rm_association_class, "fused_RM_supported")
})

test_that("an R-specific domain prevents a broad MTase product from overwriting the R role", {
  genes <- data.frame(
    locus_tag = "res_with_mtase_product",
    product = "DNA methyltransferase",
    `Signature description_Pfam` = "Type III restriction enzyme, res subunit",
    EggNOG_Preferred_name = "resA",
    translation = strrep("A", 900),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- .curation_hit_table("res_with_mtase_product", "M", "Type III")
  hits$raw_enzyme_role <- "R"
  hits$curated_blast_role <- "R"
  hits$blast_alignment_quality <- "weak"
  hits$blast_min_coverage <- 0.45

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$final_role, "R")
  expect_false(out$role_conflict)
})

test_that("Type IV family domains rescue GmrSD and Mrr but not YaeQ or YraN", {
  genes <- data.frame(
    locus_tag = c("gmrsd", "mrr", "yaeq", "yran"),
    product = c("hypothetical protein", "restriction endonuclease", "YaeQ family protein", "YraN family protein"),
    translation = c(strrep("A", 600), strrep("A", 159), strrep("A", 180), strrep("A", 160)),
    `Signature accession_Pfam` = c("PF03235;PF07510", "PF04471", "PF07152", "UPF0102"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table(genes$locus_tag, rep("R", 4), rep("Type IV", 4))
  hits$blast_alignment_quality <- "weak"
  hits$blast_min_coverage <- 0.2
  hits$partial_status <- "full_length"

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)
  rownames(out) <- out$query

  expect_identical(out["gmrsd", "curation_tier"], "high")
  expect_identical(out["mrr", "curation_tier"], "medium")
  expect_identical(out["yaeq", "curation_tier"], "excluded_noise")
  expect_identical(out["yran", "curation_tier"], "excluded_noise")
})

test_that("reference recognition metadata is not an independent retention axis", {
  genes <- data.frame(
    locus_tag = "q1",
    product = "hypothetical protein",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("q1", "R", "Type II")
  hits$blast_alignment_quality <- "moderate"
  hits$reference_rec_seq <- "GATC"

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_identical(out$independent_evidence_axes, 1L)
  expect_identical(out$curation_tier, "review")
  expect_true(out$reference_recognition_metadata)
})

test_that("compact full-reference matches override length-only partial warnings", {
  genes <- data.frame(
    locus_tag = "compact_r",
    product = "hypothetical protein",
    translation = strrep("A", 125),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("compact_r", "R", "Type II")
  hits$partial_status <- "partial_or_short"
  hits$partial_reason <- "aa_len=125 below expected_min=150"
  hits$blast_query_coverage <- 0.96
  hits$blast_reference_coverage <- 0.95

  out <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)

  expect_true(out$partial_rescued_by_reference)
  expect_identical(out$curation_tier, "medium")
})

test_that("family-specific annotations constrain primary BLAST family selection", {
  genes <- data.frame(
    locus_tag = "mod1",
    product = "Type III modification methylase Mod subunit",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("mod1", "M", "Type III")
  rebase <- data.frame(
    enzyme_name = c("M.TypeII", "M.TypeIII"),
    enz_type = c("Type II methyltransferase", "Type III methyltransferase"),
    rec_seq = "",
    org_name = "",
    sequence = c(strrep("A", 300), strrep("A", 300)),
    seq_length = 300,
    is_gold_standard = FALSE,
    rm_type = c("Type_II", "Type_III"),
    subunit = "M",
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = "mod1",
    rebase_enzyme = c("M.TypeII", "M.TypeIII"),
    pct_identity = c(0.85, 0.60),
    length = c(295, 290),
    evalue = c(1e-50, 1e-40),
    bitscore = c(300, 240),
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_identical(out$curated_blast_family, "Type III")
  expect_true(out$blast_family_compatible)
})

test_that("context-derived roles cannot be overwritten by a different BLAST role", {
  genes <- data.frame(
    locus_tag = "context_r",
    product = "P-loop NTPase",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("context_r", "R", "Type II")
  hits$evidence_mode <- "operon_context_split_R"
  hits$typeii_operon_supported <- TRUE
  rebase <- data.frame(
    enzyme_name = c("M.Wrong", "RM.Fusion", "R.Right"),
    enz_type = c("Type II methyltransferase", "Type IIG restriction enzyme/methyltransferase", "Type II restriction enzyme"),
    rec_seq = "",
    org_name = "",
    sequence = rep(strrep("A", 300), 3),
    seq_length = 300,
    is_gold_standard = FALSE,
    rm_type = "Type_II",
    subunit = c("M", "RM", "R"),
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = "context_r",
    rebase_enzyme = c("M.Wrong", "RM.Fusion", "R.Right"),
    pct_identity = c(0.90, 0.80, 0.55),
    length = c(295, 290, 285),
    evalue = c(1e-60, 1e-50, 1e-30),
    bitscore = c(350, 300, 220),
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_select_primary_blast(hits, genes, blast, rebase)

  expect_identical(out$enzyme_role, "R")
  expect_identical(out$curated_blast_match, "R.Right")
  expect_true(out$blast_role_compatible)
})

test_that("Type II context rejects AbiEii and generic TA ATPases", {
  mtase <- paste0(strrep("A", 25), "YAGGA", strrep("A", 150), "DPPY", strrep("A", 80))
  genes <- data.frame(
    locus_tag = c("motor", "toxin_atpase", "mtase"),
    contig = "ctg1",
    start = c(100, 900, 1800),
    end = c(800, 1700, 2700),
    direction = "+",
    product = c("Sll1717 KAP P-loop ATPase PF07693", "AbiEii toxin AAA_21 PF13304", "DNA adenine methylase"),
    translation = c(
      paste0(strrep("A", 40), "AGGGGGKS", strrep("A", 200), "SAT", strrep("A", 200)),
      paste0(strrep("A", 200), "DEAD", strrep("A", 200)),
      mtase
    ),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table("mtase", "M", "Type II")

  out <- DNMB:::.dnmb_rebasefinder_add_typeii_context(hits, genes)

  expect_false(any(c("motor", "toxin_atpase") %in% out$query))
})

test_that("curation workbook separates retained review noise and raw evidence", {
  genes <- data.frame(
    locus_tag = c("kept", "review", "noise"),
    product = c("DNA methyltransferase", "restriction endonuclease-like protein", "rRNA methyltransferase"),
    stringsAsFactors = FALSE
  )
  hits <- .curation_hit_table(c("kept", "review", "noise"), c("M", "R", "M"))
  hits$curation_tier <- c("high", "review", "excluded_noise")
  hits$curation_score <- c(90, 35, 0)
  hits$curation_keep <- c(TRUE, FALSE, FALSE)
  hits$curation_reason <- c("full", "review", "noise")
  hits$exclusion_reason <- c(NA, NA, "non_DNA_methyltransferase")
  hits$final_family <- hits$family_id
  hits$final_role <- hits$enzyme_role
  hits$final_component <- hits$enzyme_role
  hits$rm_association_class <- c("classic_RM_supported", "RM_association_unproven", "RM_association_unproven")
  hits$annotation_class <- c("RM_specific", "unresolved", "non_DNA_methyltransferase")
  hits$evidence_axes <- c("full_length_REBASE", "RM_annotation_or_domain", NA)
  hits$independent_evidence_axes <- c(1L, 1L, 0L)
  hits$role_conflict <- FALSE
  out_dir <- tempfile("rebase-curation-")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  openxlsx::write.xlsx(data.frame(raw = 1), file.path(out_dir, "R-M_REBASE_analysis.xlsx"))

  written <- DNMB:::.dnmb_rebasefinder_write_curation_workbook(
    out_dir, hits, genes, raw_defenseviz = data.frame(locus_tag = c("raw1", "raw2"))
  )

  expect_true(file.exists(written$raw_xlsx))
  expect_identical(
    openxlsx::getSheetNames(written$xlsx),
    c("RM_Comprehensive", "High_Identity_50pct", "Standalone_DNA_MTase",
      "Review_Candidates", "Other_Defense", "Excluded_Noise",
      "All_Augmented_Evidence", "Raw_DefenseViz", "Type_Summary", "Scoring_Criteria")
  )
  expect_identical(openxlsx::read.xlsx(written$xlsx, "RM_Comprehensive")$locus_tag, "kept")
  expect_identical(openxlsx::read.xlsx(written$xlsx, "Review_Candidates")$locus_tag, "review")
  expect_identical(openxlsx::read.xlsx(written$xlsx, "Excluded_Noise")$locus_tag, "noise")
  expect_equal(nrow(openxlsx::read.xlsx(written$xlsx, "Raw_DefenseViz")), 2L)

  raw_before <- openxlsx::read.xlsx(written$raw_xlsx)
  written_again <- DNMB:::.dnmb_rebasefinder_write_curation_workbook(
    out_dir, hits, genes, raw_defenseviz = data.frame(locus_tag = "new_raw")
  )
  expect_identical(openxlsx::read.xlsx(written_again$raw_xlsx), raw_before)
})

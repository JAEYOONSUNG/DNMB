.rebasefinder_test_hits <- function(query = "mtase1",
                                    family_id = "Type II",
                                    enzyme_role = "M") {
  data.frame(
    query = query,
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = family_id,
    hit_label = "M.TestSystem",
    enzyme_role = enzyme_role,
    evidence_mode = "high_confidence",
    substrate_label = NA_character_,
    support = "test anchor",
    typing_eligible = TRUE,
    blast_identity = 0.8,
    blast_evalue = 1e-40,
    blast_bitscore = 250,
    blast_length = 280,
    rec_seq = NA_character_,
    stringsAsFactors = FALSE
  )
}

.rebasefinder_test_normalize_ws <- function(x) {
  trimws(gsub("[[:space:]]+", " ", as.character(x)))
}

test_that("REBASEfinder gene order separates records sharing a contig description", {
  genes <- data.frame(
    locus_tag = c("record1_a", "record2_a", "record1_b", "record2_b"),
    contig = "shared WGS chromosome description",
    contig_number = c(1, 2, 1, 2),
    start = c(100, 100, 800, 800),
    end = c(700, 700, 1400, 1400),
    direction = "+",
    stringsAsFactors = FALSE
  )

  ordered <- DNMB:::.dnmb_rebasefinder_gene_order(genes)

  expect_identical(
    ordered$locus_tag,
    c("record1_a", "record1_b", "record2_a", "record2_b")
  )
  expect_equal(length(unique(ordered$contig)), 2L)
  expect_true(all(grepl("::record_[12]$", ordered$contig)))
})

test_that("REBASEfinder input preserves rows and fields containing tabs and newlines", {
  genes <- data.frame(
    locus_tag = c("gene1", "gene2"),
    start = c(1, 401),
    end = c(300, 700),
    direction = c("+", "-"),
    translation = c(strrep("A", 100), strrep("G", 100)),
    product = c("DNA methyltransferase\nsubunit M", "ATP-dependent\tendonuclease"),
    note = c(
      "first line\nsecond line\ttabbed annotation",
      "single-line annotation"
    ),
    stringsAsFactors = FALSE
  )
  out_dir <- tempfile("rebasefinder-input-")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  path <- DNMB:::.dnmb_rebasefinder_prepare_input(genes, out_dir)
  roundtrip <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\""
  )

  expect_equal(nrow(roundtrip), nrow(genes))
  expect_identical(roundtrip$locus_tag, genes$locus_tag)
  expect_equal(roundtrip$start, genes$start)
  expect_equal(roundtrip$end, genes$end)
  expect_identical(roundtrip$direction, genes$direction)
  expect_identical(roundtrip$translation, genes$translation)
  expect_identical(
    .rebasefinder_test_normalize_ws(roundtrip$product),
    .rebasefinder_test_normalize_ws(genes$product)
  )
  expect_identical(
    .rebasefinder_test_normalize_ws(roundtrip$note),
    .rebasefinder_test_normalize_ws(genes$note)
  )
})

test_that("operon paths do not bridge an opposite-strand intervening gene", {
  genes <- data.frame(
    locus_tag = c("left", "blocker", "right"),
    contig = "ctg1",
    start = c(100, 500, 900),
    end = c(400, 800, 1200),
    direction = c("+", "-", "+"),
    stringsAsFactors = FALSE
  )
  ordered <- DNMB:::.dnmb_rebasefinder_gene_order(genes)

  expect_true(DNMB:::.dnmb_rebasefinder_coherent_gene_path(ordered, 1L, 1L))
  expect_false(DNMB:::.dnmb_rebasefinder_coherent_gene_path(ordered, 1L, 3L))

  ordered$direction[[2]] <- "+"
  expect_true(DNMB:::.dnmb_rebasefinder_coherent_gene_path(ordered, 1L, 3L))
})

test_that("Type II context rescues a split same-strand R unit and excludes the opposite strand", {
  p_loop_sat <- paste0(
    strrep("A", 45),
    "GRTGSGKT",
    strrep("A", 80),
    "DEAH",
    strrep("A", 46),
    "SAT",
    strrep("A", 328)
  )
  dead_endonuclease <- paste0(strrep("A", 257), "DEAD", strrep("A", 306))
  mtase <- paste0(
    strrep("A", 39),
    "YAGGA",
    strrep("A", 152),
    "DPPY",
    strrep("A", 86)
  )
  genes <- data.frame(
    locus_tag = c("r_atpase", "r_endonuclease", "mtase1", "opposite_r"),
    contig = "ctg1",
    start = c(100, 1885, 3629, 4530),
    end = c(1644, 3585, 4489, 6230),
    direction = c("+", "+", "+", "-"),
    product = c(
      "P-loop NTPase fold protein",
      "ATP-dependent endonuclease",
      "DNA adenine methylase",
      "ATP-dependent endonuclease"
    ),
    translation = c(p_loop_sat, dead_endonuclease, mtase, dead_endonuclease),
    stringsAsFactors = FALSE
  )

  annotation_only_r <- .rebasefinder_test_hits(
    query = "r_endonuclease",
    family_id = NA_character_,
    enzyme_role = NA_character_
  )
  annotation_only_r$hit_label <- NA_character_
  annotation_only_r$evidence_mode <- "annotation_only"
  annotation_only_r$typing_eligible <- FALSE
  out <- DNMB:::.dnmb_rebasefinder_add_typeii_context(
    rbind(.rebasefinder_test_hits(), annotation_only_r),
    genes,
    max_operon_gap = 5000,
    max_intervening = 1L,
    max_neighbors = 3L
  )

  expect_true(all(c("r_atpase", "r_endonuclease", "mtase1") %in% out$query))
  expect_false("opposite_r" %in% out$query)
  expect_true(all(out$family_id[match(c("r_atpase", "r_endonuclease"), out$query)] == "Type II"))
  expect_true(all(out$enzyme_role[match(c("r_atpase", "r_endonuclease"), out$query)] == "R"))
  expect_match(out$hit_label[match("r_endonuclease", out$query)], "^typeII_context:R:")
  expect_true(all(grepl(
    "operon_context",
    out$evidence_mode[match(c("r_atpase", "r_endonuclease"), out$query)],
    fixed = TRUE
  )))
  expect_identical(
    out$operon_component[match(c("r_atpase", "r_endonuclease"), out$query)],
    c("R_motor", "R_annotation")
  )
  motor_partners <- out$typeii_context_partners[match("r_atpase", out$query)]
  nuclease_partners <- out$typeii_context_partners[match("r_endonuclease", out$query)]
  expect_match(motor_partners, "mtase1:M")
  expect_match(motor_partners, "r_endonuclease:R_annotation")
  expect_false(grepl("r_atpase:", motor_partners, fixed = TRUE))
  expect_match(nuclease_partners, "mtase1:M")
  expect_match(nuclease_partners, "r_atpase:R_motor")
  expect_false(grepl("r_endonuclease:", nuclease_partners, fixed = TRUE))
  expect_false(any(out$typing_eligible[match(c("r_atpase", "r_endonuclease"), out$query)]))
  expect_true(isTRUE(out$typeii_operon_supported[match("mtase1", out$query)]))
  expect_true(isTRUE(out$mtase_motif_verified[match("mtase1", out$query)]))
  expect_identical(out$mtase_sam_motif[match("mtase1", out$query)], "YAGGA")
  expect_identical(out$mtase_catalytic_motif[match("mtase1", out$query)], "DPPY")
})

test_that("Type II motif context does not relabel established Type I or III systems", {
  mtase <- paste0(strrep("A", 25), "YAGGA", strrep("A", 120), "DPPY", strrep("A", 80))
  rease <- paste0(strrep("A", 180), "DEAD", strrep("A", 220))
  genes <- data.frame(
    locus_tag = c("known_m", "known_r"),
    contig = "ctg1",
    start = c(100, 1000),
    end = c(900, 2200),
    direction = "+",
    product = c("DNA methyltransferase", "ATP-dependent restriction endonuclease"),
    translation = c(mtase, rease),
    stringsAsFactors = FALSE
  )

  for (family in c("Type I", "Type III")) {
    hits <- rbind(
      .rebasefinder_test_hits("known_m", family, "M"),
      .rebasefinder_test_hits("known_r", family, "R")
    )
    out <- DNMB:::.dnmb_rebasefinder_add_typeii_context(hits, genes)
    expect_identical(out$family_id, rep(family, 2L), info = family)
    if ("typeii_operon_supported" %in% names(out)) {
      expect_false(any(out$typeii_operon_supported %in% TRUE), info = family)
    }
  }
})

test_that("a single unannotated nuclease-like pattern does not create a Type II operon", {
  mtase <- paste0(strrep("A", 25), "YAGGA", strrep("A", 120), "DPPY", strrep("A", 80))
  weak_r <- paste0(strrep("A", 80), "PDAAAAAELK", strrep("A", 100))
  genes <- data.frame(
    locus_tag = c("mtase", "weak_r"),
    contig = "ctg1",
    start = c(100, 1000),
    end = c(900, 1800),
    direction = "+",
    product = c("DNA adenine methylase", "hypothetical protein"),
    translation = c(mtase, weak_r),
    stringsAsFactors = FALSE
  )
  hits <- .rebasefinder_test_hits("mtase", "Type II", "M")
  hits$typing_eligible <- FALSE

  context <- DNMB:::.dnmb_rebasefinder_add_typeii_context(hits, genes)
  rownames(context) <- context$query

  expect_false(context["mtase", "typeii_operon_candidate"] %in% TRUE)
  expect_false("weak_r" %in% context$query)

  context$partial_status <- "full_length"
  curated <- DNMB:::.dnmb_rebasefinder_curate_hits(context, genes)
  rownames(curated) <- curated$query
  expect_identical(curated["mtase", "curation_tier"], "medium")
  expect_identical(curated["mtase", "rm_association_class"], "DNA_MTase_RM_association_unproven")
})

test_that("full-length Type II homology links hypothetical M and R without inventing a nuclease label", {
  genes <- data.frame(
    locus_tag = c("m_homology", "r_homology"),
    contig = "ctg1",
    start = c(100, 1100),
    end = c(900, 2100),
    direction = "+",
    product = "hypothetical protein",
    translation = c(strrep("A", 267), strrep("A", 334)),
    stringsAsFactors = FALSE
  )
  hits <- rbind(
    .rebasefinder_test_hits("m_homology", "Type II", "M"),
    .rebasefinder_test_hits("r_homology", "Type II", "R")
  )
  hits$blast_alignment_quality <- "strong"
  hits$blast_role_compatible <- TRUE
  hits$blast_family_compatible <- TRUE

  out <- DNMB:::.dnmb_rebasefinder_add_typeii_context(hits, genes)
  rownames(out) <- out$query

  expect_true(out["m_homology", "typeii_operon_supported"])
  expect_true(out["r_homology", "typeii_operon_supported"])
  expect_identical(out["r_homology", "operon_component"], "R_homology")
  expect_identical(out["r_homology", "typeii_context_sources"], "rebase_homology")
  expect_match(out["m_homology", "typeii_context_partners"], "r_homology:R_homology")
})

test_that("Type I context respects max_intervening even on a coherent strand", {
  genes <- data.frame(
    locus_tag = c("hsdM", "spacer1", "spacer2", "hsdR", "hsdS"),
    contig = "ctg1",
    start = seq(100, 4100, by = 1000),
    end = seq(900, 4900, by = 1000),
    direction = "+",
    product = "hypothetical protein",
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- rbind(
    .rebasefinder_test_hits("hsdM", "Type I", "M"),
    .rebasefinder_test_hits("hsdR", "Type I", "R"),
    .rebasefinder_test_hits("hsdS", "Type I", "S")
  )
  hits$blast_alignment_quality <- "strong"
  hits$blast_role_compatible <- TRUE
  hits$blast_family_compatible <- TRUE

  out <- DNMB:::.dnmb_rebasefinder_add_typei_context(
    hits, genes, max_operon_gap = 5000, max_intervening = 1L, max_neighbors = 4L
  )

  expect_false(any(out$typei_operon_supported %in% TRUE))
  expect_false(any(out$typei_context_status == "complete_mrs", na.rm = TRUE))
})

test_that("trusted Type I endpoints bridge a short same-strand defense-island insertion", {
  genes <- data.frame(
    locus_tag = c("hsdM", "hsdS", "kwaA", "kwaB", "insert", "hsdR"),
    contig = "ctg1",
    start = seq(100, 5100, by = 1000),
    end = seq(900, 5900, by = 1000),
    direction = "+",
    product = c(
      "type I restriction-modification system subunit M",
      "restriction endonuclease subunit S",
      "anti-phage protein KwaA",
      "anti-phage protein KwaB",
      "hypothetical protein",
      "type I restriction-modification system endonuclease subunit R"
    ),
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- rbind(
    .rebasefinder_test_hits("hsdM", "Type I", "M"),
    .rebasefinder_test_hits("hsdS", "Type I", "S"),
    .rebasefinder_test_hits("hsdR", "Type I", "R")
  )
  hits$blast_alignment_quality <- "strong"
  hits$blast_role_compatible <- TRUE
  hits$blast_family_compatible <- TRUE

  out <- DNMB:::.dnmb_rebasefinder_add_typei_context(
    hits, genes, max_operon_gap = 6000, max_intervening = 1L, max_neighbors = 4L
  )
  rownames(out) <- out$query

  expect_true(out["hsdS", "typei_operon_supported"])
  expect_identical(out["hsdS", "typei_context_status"], "complete_mrs")
  expect_match(out["hsdS", "typei_context_partners"], "hsdR:R")
  expect_length(unique(out[c("hsdM", "hsdS", "hsdR"), "operon_id"]), 1L)
  expect_match(out["hsdS", "operon_id"], "^DNMB_TypeI_")
})

test_that("Type I seeds stage a distant annotated HsdR for supplemental BLAST", {
  genes <- data.frame(
    locus_tag = c("hsdM", "hsdS", "insert1", "insert2", "insert3", "hsdR_candidate"),
    contig = "ctg1",
    start = seq(100, 5100, by = 1000),
    end = seq(900, 5900, by = 1000),
    direction = "+",
    product = c(
      "type I restriction-modification system subunit M",
      "type I restriction specificity subunit HsdS",
      "hypothetical protein",
      "hypothetical protein",
      "hypothetical protein",
      "type I restriction-modification system endonuclease subunit HsdR"
    ),
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- rbind(
    .rebasefinder_test_hits("hsdM", "Type I", "M"),
    .rebasefinder_test_hits("hsdS", "Type I", "S")
  )
  hits$blast_alignment_quality <- "strong"
  hits$blast_role_compatible <- TRUE
  hits$blast_family_compatible <- TRUE

  out <- DNMB:::.dnmb_rebasefinder_add_typei_context(
    hits, genes, max_operon_gap = 6000, max_intervening = 1L, max_neighbors = 4L
  )
  rownames(out) <- out$query

  expect_true("hsdR_candidate" %in% out$query)
  expect_identical(out["hsdR_candidate", "family_id"], "Type I")
  expect_identical(out["hsdR_candidate", "enzyme_role"], "R")
  expect_identical(out["hsdR_candidate", "evidence_mode"], "operon_context")
  expect_false(out["hsdR_candidate", "typing_eligible"])
  expect_true(out["hsdR_candidate", "typei_operon_supported"])
  expect_true(DNMB:::.dnmb_rebasefinder_supplemental_query_mask(out)[match("hsdR_candidate", out$query)])
})

test_that("a distant strong Type I R sequence is staged unless annotated as helicase noise", {
  strong_r <- paste0(
    strrep("A", 40), "PDAAAAADEK", strrep("A", 80), "GSGKT",
    strrep("A", 120), "DEAH", strrep("A", 70), "SAT", strrep("A", 300)
  )
  genes <- data.frame(
    locus_tag = c("hsdM", "hsdS", "insert1", "insert2", "insert3", "motif_r"),
    contig = "ctg1",
    start = seq(100, 5100, by = 1000),
    end = seq(900, 5900, by = 1000),
    direction = "+",
    product = c(
      "type I restriction-modification system subunit M",
      "type I restriction specificity subunit HsdS",
      rep("hypothetical protein", 4)
    ),
    translation = c(rep(strrep("A", 300), 5), strong_r),
    stringsAsFactors = FALSE
  )
  hits <- rbind(
    .rebasefinder_test_hits("hsdM", "Type I", "M"),
    .rebasefinder_test_hits("hsdS", "Type I", "S")
  )
  hits$blast_alignment_quality <- "strong"
  hits$blast_role_compatible <- TRUE
  hits$blast_family_compatible <- TRUE

  staged <- DNMB:::.dnmb_rebasefinder_add_typei_context(
    hits, genes, max_operon_gap = 6000, max_intervening = 1L, max_neighbors = 4L
  )
  rownames(staged) <- staged$query
  expect_true("motif_r" %in% staged$query)
  expect_match(staged["motif_r", "hit_label"], ":R:sequence_motif$")
  expect_true(DNMB:::.dnmb_rebasefinder_supplemental_query_mask(staged)[match("motif_r", staged$query)])

  genes$product[genes$locus_tag == "motif_r"] <- "ATP-dependent RNA helicase DbpA"
  filtered <- DNMB:::.dnmb_rebasefinder_add_typei_context(
    hits, genes, max_operon_gap = 6000, max_intervening = 1L, max_neighbors = 4L
  )
  expect_false("motif_r" %in% filtered$query)
})

test_that("complementary adjacent Type I R fragments are grouped for review", {
  genes <- data.frame(
    locus_tag = c("split_r_left", "split_r_right", "far_r"),
    contig = "ctg1",
    start = c(100, 1250, 10000),
    end = c(1308, 2550, 11200),
    direction = "-",
    product = "type I restriction endonuclease subunit R",
    translation = strrep("A", 430),
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query = genes$locus_tag,
    family_id = "Type I",
    enzyme_role = "R",
    hit_label = "R.Shared",
    curated_blast_subject_id = "R.Shared",
    blast_reference_coverage = c(0.40, 0.53, 0.45),
    blast_query_coverage = c(0.99, 0.88, 0.95),
    support = NA_character_,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_add_split_r_evidence(hits, genes)
  rownames(out) <- out$query

  expect_true(out["split_r_left", "split_r_pair_candidate"])
  expect_true(out["split_r_right", "split_r_pair_candidate"])
  expect_false(out["far_r", "split_r_pair_candidate"])
  expect_identical(out["split_r_left", "split_r_partner"], "split_r_right")
  expect_equal(out["split_r_left", "split_r_combined_reference_coverage"], 0.93)
  expect_identical(
    out["split_r_left", "operon_id"],
    out["split_r_right", "operon_id"]
  )
  expect_match(out["split_r_left", "support"], "split_R_pair=")
})

test_that("Type I context rejects a different family unless alternate Type I homology supports it", {
  genes <- data.frame(
    locus_tag = c("hsdR", "other_m", "hsdS"),
    contig = "ctg1",
    start = c(100, 1000, 1900),
    end = c(900, 1800, 2700),
    direction = "+",
    product = c(
      "type I restriction endonuclease subunit R",
      "DNA methyltransferase",
      "type I restriction specificity subunit HsdS"
    ),
    translation = strrep("A", 300),
    stringsAsFactors = FALSE
  )
  hits <- rbind(
    .rebasefinder_test_hits("hsdR", "Type I", "R"),
    .rebasefinder_test_hits("other_m", "Type II", "M")
  )
  hits$blast_alignment_quality <- "strong"
  hits$blast_role_compatible <- TRUE
  hits$blast_family_compatible <- TRUE

  blocked <- DNMB:::.dnmb_rebasefinder_add_typei_context(hits, genes)
  expect_false(any(blocked$typei_operon_supported %in% TRUE))

  hits$blast_supported_typei_roles <- c(NA_character_, "M")
  reassigned <- DNMB:::.dnmb_rebasefinder_add_typei_context(hits, genes)
  rownames(reassigned) <- reassigned$query
  expect_true(reassigned["other_m", "typei_operon_supported"])
  expect_identical(reassigned["other_m", "typei_context_status"], "complete_mrs")
  expect_match(reassigned["hsdR", "typei_context_sources"], "M=rebase_alternate")
})

test_that("context refresh clears stale supported calls when the seed disappears", {
  genes <- data.frame(
    locus_tag = "not_type_i",
    contig = "ctg1",
    start = 100,
    end = 900,
    direction = "+",
    product = "hypothetical protein",
    translation = strrep("A", 267),
    stringsAsFactors = FALSE
  )
  hits <- .rebasefinder_test_hits("not_type_i", "Type II", "M")
  hits$typei_context_status <- "complete_mrs"
  hits$typei_context_roles <- "M/R/S"
  hits$typei_operon_supported <- TRUE
  hits$typei_context_anchor <- "old_anchor"
  hits$operon_id <- "DNMB_typeI_old_anchor"
  hits$support <- "test anchor; typeI_context=complete_mrs; partners=old"
  hits$raw_typing_eligible <- FALSE

  out <- DNMB:::.dnmb_rebasefinder_add_typei_context(hits, genes)

  expect_false(any(out$typei_operon_supported %in% TRUE))
  expect_true(is.na(out$typei_context_status))
  expect_false(grepl("typeI_context=", out$support, fixed = TRUE))
  expect_false(out$typing_eligible)
  expect_true(is.na(out$operon_id))
})

test_that("empty REBASE normalization retains the supplemental BLAST schema", {
  out <- DNMB:::.dnmb_rebasefinder_normalize_hits(data.frame())
  expect_equal(nrow(out), 0L)
  expect_true(all(c(
    "blast_identity", "blast_evalue", "blast_bitscore", "blast_length",
    "rec_seq", "operon_id"
  ) %in% names(out)))
  expect_type(out$blast_identity, "double")
  expect_type(out$rec_seq, "character")
})

test_that("DefenseViz recognition sites normalize as reference-only metadata", {
  raw <- data.frame(
    locus_tag = "low_identity_r",
    rm_type = "Type II",
    subunit = "R",
    blast_match = "R.LowIdentity",
    blast_identity = 0.20,
    blast_evalue = 1e-6,
    blast_bitscore = 80,
    blast_length = 120,
    rec_seq = "GATC",
    operon_id = "1",
    partner_locus_tag = NA_character_,
    stringsAsFactors = FALSE
  )
  out <- DNMB:::.dnmb_rebasefinder_normalize_hits(raw)

  expect_identical(out$reference_rec_seq, "GATC")
  expect_true(is.na(out$rec_seq))
  expect_true(is.na(out$substrate_label))
  expect_identical(out$recognition_match, "R.LowIdentity")
  expect_identical(out$recognition_match_subject_id, "R.LowIdentity")
  expect_false(out$typing_eligible)
  expect_match(out$support, "reference_rec_seq=GATC")
})

test_that("methyltransferase calls require SAM and catalytic motifs in the same protein", {
  paired <- paste0(strrep("A", 25), "YAGGA", strrep("A", 90), "DPPY", strrep("A", 80))
  sam_only <- paste0(strrep("A", 25), "YAGGA", strrep("A", 174))
  catalytic_only <- paste0(strrep("A", 120), "DPPY", strrep("A", 80))

  type_i <- lapply(
    list(paired = paired, sam_only = sam_only, catalytic_only = catalytic_only),
    DNMB:::.dnmb_rebasefinder_typei_sequence_signals
  )
  type_iii <- lapply(
    list(paired = paired, sam_only = sam_only, catalytic_only = catalytic_only),
    DNMB:::.dnmb_rebasefinder_typeiii_sequence_signals
  )

  expect_identical(type_i$paired$role, "M")
  expect_true(is.na(type_i$sam_only$role))
  expect_true(is.na(type_i$catalytic_only$role))
  expect_identical(type_iii$paired$role, "M")
  expect_true(is.na(type_iii$sam_only$role))
  expect_true(is.na(type_iii$catalytic_only$role))
})

test_that("Bairoch parsing retains restriction entries without methylation fields", {
  bairoch <- tempfile("rebase-bairoch-", fileext = ".txt")
  on.exit(unlink(bairoch, force = TRUE), add = TRUE)
  writeLines(
    c(
      "ID   M.TestSystem",
      "RS   GATC, ?;",
      "MS   2(m6A),-2(m6A);",
      "//",
      "ID   R.TestSystem",
      "RS   GATC, 1;",
      "//"
    ),
    bairoch
  )

  parsed <- DNMB:::.dnmb_rebasefinder_parse_bairoch(bairoch)
  rownames(parsed) <- parsed$enzyme_name

  expect_setequal(parsed$enzyme_name, c("M.TestSystem", "R.TestSystem"))
  expect_identical(parsed["M.TestSystem", "rec_seq"], "GATC")
  expect_identical(parsed["M.TestSystem", "meth_type"], "m6A")
  expect_identical(parsed["R.TestSystem", "rec_seq"], "GATC")
  expect_true(is.na(parsed["R.TestSystem", "meth_type"]))
  expect_true(is.na(parsed["R.TestSystem", "meth_pos"]))
})

test_that("REBASE recognition metadata uses exact hits then cognate system partners", {
  hits <- rbind(
    .rebasefinder_test_hits(query = "exact", enzyme_role = "R"),
    .rebasefinder_test_hits(query = "partner", enzyme_role = "R")
  )
  hits$hit_label <- c("R.ExactSystem", "R.PartnerSystem")
  rebase_data <- data.frame(
    enzyme_name = c("R.ExactSystem", "R.PartnerSystem", "M.PartnerSystem"),
    enz_type = c("Type II restriction enzyme", "Type II restriction enzyme", "Type II methyltransferase"),
    rec_seq = c("CCGG", "", "GATC"),
    rm_type = "Type_II",
    subunit = c("R", "R", "M"),
    stringsAsFactors = FALSE
  )

  enriched <- DNMB:::.dnmb_rebasefinder_enrich_recognition(hits, rebase_data)
  rownames(enriched) <- enriched$query

  expect_true(is.na(enriched["exact", "rec_seq"]))
  expect_identical(enriched["exact", "reference_rec_seq"], "CCGG")
  expect_identical(enriched["exact", "recognition_match"], "R.ExactSystem")
  expect_identical(enriched["exact", "recognition_donor"], "R.ExactSystem")
  expect_match(enriched["exact", "recognition_source"], "exact", ignore.case = TRUE)
  expect_true(is.na(enriched["partner", "rec_seq"]))
  expect_identical(enriched["partner", "reference_rec_seq"], "GATC")
  expect_identical(enriched["partner", "recognition_match"], "R.PartnerSystem")
  expect_identical(enriched["partner", "recognition_donor"], "M.PartnerSystem")
  expect_match(
    enriched["partner", "recognition_source"],
    "system|partner|cognate",
    ignore.case = TRUE
  )
})

test_that("recognition enrichment preserves an existing exact reference site", {
  hits <- .rebasefinder_test_hits("query1", "Type II", "R")
  hits$hit_label <- "R.Duplicate"
  hits$reference_rec_seq <- "CCGG"
  lookup <- data.frame(
    enzyme_name = "R.Duplicate",
    rec_seq = "GATC",
    is_gold_standard = TRUE,
    stringsAsFactors = FALSE
  )

  enriched <- DNMB:::.dnmb_rebasefinder_enrich_recognition(hits, lookup)

  expect_identical(enriched$reference_rec_seq, "CCGG")
  expect_identical(enriched$recognition_source, "pipeline_reference_exact")
  expect_identical(enriched$recognition_match_subject_id, "R.Duplicate")
})

test_that("REBASE PacBio index parsing creates direct strain links", {
  html <- tempfile("rebase-pacbio-index-", fileext = ".html")
  on.exit(unlink(html, force = TRUE), add = TRUE)
  writeLines(
    c(
      "<html><body><table>",
      "<tr><td><a href='/cgi-bin/pacbioget?42846'>Weeksella virosa NCTC 11635</a></td><td>UHIU01000004.1</td></tr>",
      "</table></body></html>"
    ),
    html
  )

  index <- DNMB:::.dnmb_rebasefinder_parse_rebase_pacbio_index(html)
  mapped <- DNMB:::.dnmb_rebasefinder_match_rebase_pacbio(
    c("UHIU01000004", "CP999999"),
    index
  )

  expect_equal(nrow(index), 1L)
  expect_identical(index$nucleotide_accession, "UHIU01000004.1")
  expect_identical(index$organism, "Weeksella virosa NCTC 11635")
  expect_identical(
    index$pacbio_url,
    "https://rebase.neb.com/cgi-bin/pacbioget?42846"
  )
  expect_identical(mapped$pacbio_available, c(TRUE, FALSE))
  expect_identical(mapped$pacbio_url[[1]], index$pacbio_url)
  expect_true(is.na(mapped$pacbio_url[[2]]))
})

test_that("secondary BLAST recognition metadata preserves the primary assignment", {
  cache_root <- tempfile("rebasefinder-cache-")
  cache_dir <- DNMB:::.dnmb_rebasefinder_cache_dir(cache_root, create = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines(
    c(
      paste0(
        ">REBASE:M.Primary",
        "\tEnzType:Type II methyltransferase",
        "\tRecSeq:",
        "\tGenBank:CP000001",
        "\tLocus:primary_locus",
        "\tProteinId:PRIMARY_1"
      ),
      "AAAA",
      paste0(
        ">REBASE:Wvi11635IP",
        "\tEnzType:Type II restriction enzyme",
        "\tRecSeq:GAGNNNNNTAC",
        "\tGenBank:UHIU01000004",
        "\tLocus:recognition_locus",
        "\tProteinId:RECOGNITION_1"
      ),
      "AAAA"
    ),
    file.path(cache_dir, "REBASE_protein_seqs.txt")
  )

  bairoch <- data.frame(
    enzyme_name = "Unrelated",
    rec_seq = "GATC",
    meth_type = NA_character_,
    meth_pos = NA_character_,
    meth_all = NA_character_,
    stringsAsFactors = FALSE
  )
  attr(bairoch, "dnmb_bairoch_schema") <- 2L
  saveRDS(bairoch, file.path(cache_dir, "rebase_bairoch_lookup.rds"))

  pacbio_index <- data.frame(
    nucleotide_accession = "UHIU01000004.1",
    organism = "Weeksella virosa NCTC 11635",
    rebase_org_id = "42846",
    pacbio_url = "https://rebase.neb.com/cgi-bin/pacbioget?42846",
    .accession_key = "UHIU01000004",
    stringsAsFactors = FALSE
  )
  saveRDS(pacbio_index, file.path(cache_dir, "rebase_pacbio_index.rds"))

  hits <- .rebasefinder_test_hits(query = "query1")
  hits$hit_label <- "M.Primary"
  blast <- data.frame(
    query_id = c("query1", "query1"),
    rebase_enzyme = c("M.Primary", "Wvi11635IP"),
    pct_identity = c(0.70, 0.61),
    evalue = c(1e-60, 1e-45),
    bitscore = c(310, 280),
    length = c(285, 280),
    stringsAsFactors = FALSE
  )
  rebase_data <- data.frame(
    enzyme_name = c("M.Primary", "Wvi11635IP"),
    rec_seq = c("", "GAGNNNNNTAC"),
    is_gold_standard = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  enriched <- DNMB:::.dnmb_rebasefinder_enrich_blast_metadata(
    hits,
    blast_tbl = blast,
    rebase_data = rebase_data,
    cache_root = cache_root
  )

  expect_identical(enriched$hit_label, "M.Primary")
  expect_identical(enriched$recognition_match, "Wvi11635IP")
  expect_identical(enriched$reference_rec_seq, "GAGNNNNNTAC")
  expect_false(enriched$recognition_transfer_eligible)
  expect_true(is.na(enriched$rec_seq))
  expect_identical(enriched$recognition_reference_accession, "UHIU01000004")
  expect_true(enriched$recognition_match_pacbio_available)
  expect_identical(
    enriched$recognition_match_pacbio_link,
    "https://rebase.neb.com/cgi-bin/pacbioget?42846"
  )
})

test_that("conflicting supplemental BLAST keeps recognition as reference-only metadata", {
  cache_root <- tempfile("rebasefinder-conflict-cache-")
  cache_dir <- DNMB:::.dnmb_rebasefinder_cache_dir(cache_root, create = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines(
    c(
      paste0(
        ">REBASE:M.Aac101CI",
        "\tEnzType:Type II methyltransferase",
        "\tRecSeq:GATC",
        "\tGenBank:CP000003",
        "\tLocus:aac_locus",
        "\tProteinId:AAC_1"
      ),
      "AAAA"
    ),
    file.path(cache_dir, "REBASE_protein_seqs.txt")
  )
  bairoch <- data.frame(
    enzyme_name = "M.Aac101CI",
    rec_seq = "GATC",
    meth_type = "m6A",
    meth_pos = "2",
    meth_all = "2(m6A)",
    stringsAsFactors = FALSE
  )
  attr(bairoch, "dnmb_bairoch_schema") <- 2L
  saveRDS(bairoch, file.path(cache_dir, "rebase_bairoch_lookup.rds"))
  saveRDS(
    data.frame(
      nucleotide_accession = "CP999999",
      organism = "Unrelated strain",
      rebase_org_id = "1",
      pacbio_url = "https://rebase.neb.com/cgi-bin/pacbioget?1",
      .accession_key = "CP999999",
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "rebase_pacbio_index.rds")
  )

  hits <- .rebasefinder_test_hits(query = "context_r", enzyme_role = "R")
  hits$hit_label <- "typeII_context:R:nuclease_motif"
  hits$typing_eligible <- FALSE
  hits$supplemental_blast_match <- "M.Aac101CI"
  hits$supplemental_blast_context_conflict <- TRUE
  blast <- data.frame(
    query_id = "context_r",
    rebase_enzyme = "M.Aac101CI",
    pct_identity = 0.72,
    evalue = 1e-50,
    bitscore = 300,
    length = 280,
    stringsAsFactors = FALSE
  )
  rebase_data <- data.frame(
    enzyme_name = "M.Aac101CI",
    rec_seq = "GATC",
    is_gold_standard = TRUE,
    stringsAsFactors = FALSE
  )

  enriched <- DNMB:::.dnmb_rebasefinder_enrich_blast_metadata(
    hits,
    blast_tbl = blast,
    rebase_data = rebase_data,
    cache_root = cache_root
  )

  expect_identical(enriched$hit_label, "typeII_context:R:nuclease_motif")
  expect_identical(enriched$recognition_match, "M.Aac101CI")
  expect_identical(enriched$reference_rec_seq, "GATC")
  expect_false(enriched$recognition_transfer_eligible)
  expect_true(is.na(enriched$rec_seq))
  expect_true(is.na(enriched$rebase_reference_accession))
  expect_identical(enriched$recognition_reference_accession, "CP000003")
})

test_that("recognition donor provenance is separate and short hits are ignored", {
  cache_root <- tempfile("rebasefinder-donor-cache-")
  cache_dir <- DNMB:::.dnmb_rebasefinder_cache_dir(cache_root, create = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines(
    c(
      paste0(
        ">REBASE:R.TestSystem\tEnzType:Type II restriction enzyme\tRecSeq:",
        "\tGenBank:CP000010\tLocus:r_locus\tProteinId:R_1"
      ),
      "AAAA",
      paste0(
        ">REBASE:M.TestSystem\tEnzType:Type II methyltransferase\tRecSeq:GATC",
        "\tGenBank:CP000011\tLocus:m_locus\tProteinId:M_1"
      ),
      "AAAA"
    ),
    file.path(cache_dir, "REBASE_protein_seqs.txt")
  )
  bairoch <- data.frame(
    enzyme_name = "M.TestSystem",
    rec_seq = "GATC",
    meth_type = "m6A",
    meth_pos = "2",
    meth_all = "2(m6A)",
    stringsAsFactors = FALSE
  )
  attr(bairoch, "dnmb_bairoch_schema") <- 2L
  saveRDS(bairoch, file.path(cache_dir, "rebase_bairoch_lookup.rds"))
  saveRDS(
    data.frame(
      nucleotide_accession = c("CP000010", "CP000011"),
      organism = c("Matched strain", "Donor strain"),
      rebase_org_id = c("10", "11"),
      pacbio_url = paste0("https://rebase.neb.com/cgi-bin/pacbioget?", c("10", "11")),
      .accession_key = c("CP000010", "CP000011"),
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "rebase_pacbio_index.rds")
  )
  rebase_data <- data.frame(
    enzyme_name = c("R.TestSystem", "M.TestSystem"),
    rec_seq = c("", "GATC"),
    is_gold_standard = TRUE,
    stringsAsFactors = FALSE
  )
  blast <- data.frame(
    query_id = "query_r",
    rebase_enzyme = "R.TestSystem",
    pct_identity = 0.80,
    evalue = 1e-40,
    bitscore = 250,
    length = 220,
    stringsAsFactors = FALSE
  )
  hits <- .rebasefinder_test_hits("query_r", "Type II", "R")
  hits$hit_label <- "R.TestSystem"
  hits$rec_seq <- NA_character_

  enriched <- DNMB:::.dnmb_rebasefinder_enrich_blast_metadata(
    hits,
    blast_tbl = blast,
    rebase_data = rebase_data,
    cache_root = cache_root
  )

  expect_identical(enriched$recognition_match, "R.TestSystem")
  expect_identical(enriched$recognition_donor, "M.TestSystem")
  expect_identical(enriched$reference_rec_seq, "GATC")
  expect_true(is.na(enriched$rec_seq))
  expect_false(enriched$recognition_transfer_eligible)
  expect_identical(enriched$recognition_match_reference_accession, "CP000010")
  expect_identical(enriched$recognition_reference_accession, "CP000011")
  expect_identical(
    enriched$recognition_match_pacbio_link,
    "https://rebase.neb.com/cgi-bin/pacbioget?10"
  )
  expect_identical(
    enriched$recognition_donor_pacbio_link,
    "https://rebase.neb.com/cgi-bin/pacbioget?11"
  )

  short_hits <- .rebasefinder_test_hits("short_query", "Type II", "R")
  short_hits$hit_label <- "typeII_context:R:nuclease_motif"
  short_blast <- blast
  short_blast$query_id <- "short_query"
  short_blast$rebase_enzyme <- "M.TestSystem"
  short_blast$length <- 20
  short <- DNMB:::.dnmb_rebasefinder_enrich_blast_metadata(
    short_hits,
    blast_tbl = short_blast,
    rebase_data = rebase_data,
    cache_root = cache_root,
    min_length = 50
  )
  expect_true(is.na(short$reference_rec_seq))
  expect_true(is.na(short$recognition_match))
  expect_true(is.na(short$rec_seq))
})

test_that("duplicate REBASE subject suffix preserves the exact accession and PacBio strain", {
  cache_root <- tempfile("rebasefinder-duplicate-cache-")
  cache_dir <- DNMB:::.dnmb_rebasefinder_cache_dir(cache_root, create = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines(
    c(
      ">REBASE:Eco47I\tEnzType:Type II restriction enzyme\tRecSeq:GATC\tGenBank:GQ231553\tLocus:eco_a\tProteinId:ECO_A",
      "AAAA",
      ">REBASE:Eco47I\tEnzType:Type II restriction enzyme\tRecSeq:CCGG\tGenBank:X82105\tLocus:eco_b\tProteinId:ECO_B",
      "AAAA"
    ),
    file.path(cache_dir, "REBASE_protein_seqs.txt")
  )
  rebase_data <- data.frame(
    enzyme_name = c("Eco47I", "Eco47I"),
    rec_seq = c("GATC", "CCGG"),
    is_gold_standard = TRUE,
    stringsAsFactors = FALSE
  )
  saveRDS(rebase_data, file.path(cache_dir, "rebase_data.rds"))
  bairoch <- data.frame(
    enzyme_name = "Unrelated",
    rec_seq = "AAGCTT",
    meth_type = NA_character_,
    meth_pos = NA_character_,
    meth_all = NA_character_,
    stringsAsFactors = FALSE
  )
  attr(bairoch, "dnmb_bairoch_schema") <- 2L
  saveRDS(bairoch, file.path(cache_dir, "rebase_bairoch_lookup.rds"))
  saveRDS(
    data.frame(
      nucleotide_accession = c("GQ231553", "X82105"),
      organism = c("A strain", "B strain"),
      rebase_org_id = c("101", "102"),
      pacbio_url = paste0("https://rebase.neb.com/cgi-bin/pacbioget?", c("101", "102")),
      .accession_key = c("GQ231553", "X82105"),
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "rebase_pacbio_index.rds")
  )
  hits <- .rebasefinder_test_hits("duplicate_query", "Type II", "R")
  hits$hit_label <- "Eco47I_2"
  blast <- data.frame(
    query_id = "duplicate_query",
    rebase_enzyme = "Eco47I_2",
    pct_identity = 0.75,
    evalue = 1e-50,
    bitscore = 280,
    length = 250,
    stringsAsFactors = FALSE
  )

  enriched <- DNMB:::.dnmb_rebasefinder_enrich_blast_metadata(
    hits,
    blast_tbl = blast,
    rebase_data = rebase_data,
    cache_root = cache_root
  )

  expect_identical(enriched$reference_rec_seq, "CCGG")
  expect_identical(enriched$recognition_match, "Eco47I")
  expect_identical(enriched$recognition_match_subject_id, "Eco47I_2")
  expect_identical(enriched$recognition_match_reference_accession, "X82105")
  expect_identical(enriched$recognition_reference_accession, "X82105")
  expect_identical(enriched$recognition_donor_accession, "X82105")
  expect_identical(enriched$recognition_match_pacbio_organism, "B strain")
  expect_identical(
    enriched$recognition_donor_pacbio_link,
    "https://rebase.neb.com/cgi-bin/pacbioget?102"
  )
})

test_that("duplicate metadata follows the sorted REBASE record rather than raw FASTA order", {
  cache_root <- tempfile("rebasefinder-reordered-duplicate-")
  cache_dir <- DNMB:::.dnmb_rebasefinder_cache_dir(cache_root, create = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)

  raw_a <- strrep("A", 80)
  raw_b <- strrep("C", 80)
  writeLines(
    c(
      ">REBASE:Eco47I\tEnzType:Type II restriction enzyme\tRecSeq:GATC\tGenBank:CP100001\tLocus:eco_a\tProteinId:ECO_A\tOrgName:Raw A strain",
      raw_a,
      ">REBASE:Eco47I\tEnzType:Type II restriction enzyme\tRecSeq:CCGG\tGenBank:CP100002\tLocus:eco_b\tProteinId:ECO_B\tOrgName:Raw B strain",
      raw_b
    ),
    file.path(cache_dir, "REBASE_protein_seqs.txt")
  )
  # DefenseViz sorts parsed records before writing its BLAST FASTA. Thus _2
  # identifies row 2 below (raw record A), not raw FASTA occurrence 2.
  rebase_data <- data.frame(
    enzyme_name = c("Eco47I", "Eco47I"),
    enz_type = "Type II restriction enzyme",
    rec_seq = c("CCGG", "GATC"),
    org_name = c("Raw B strain", "Raw A strain"),
    sequence = c(raw_b, raw_a),
    seq_length = 80,
    is_gold_standard = TRUE,
    stringsAsFactors = FALSE
  )
  saveRDS(rebase_data, file.path(cache_dir, "rebase_data.rds"))
  bairoch <- data.frame(
    enzyme_name = character(), rec_seq = character(), meth_type = character(),
    meth_pos = character(), meth_all = character(), stringsAsFactors = FALSE
  )
  attr(bairoch, "dnmb_bairoch_schema") <- 2L
  saveRDS(bairoch, file.path(cache_dir, "rebase_bairoch_lookup.rds"))
  saveRDS(
    data.frame(
      nucleotide_accession = "CP100001",
      organism = "Raw A strain",
      rebase_org_id = "1001",
      pacbio_url = "https://rebase.neb.com/cgi-bin/pacbioget?1001",
      .accession_key = "CP100001",
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "rebase_pacbio_index.rds")
  )

  metadata <- DNMB:::.dnmb_rebasefinder_reference_header_metadata(
    c("Eco47I_1", "Eco47I_2"),
    cache_root = cache_root,
    rebase_data = rebase_data
  )
  expect_identical(metadata$genbank_accession, c("CP100002", "CP100001"))
  metadata_reordered <- DNMB:::.dnmb_rebasefinder_reference_header_metadata(
    c("Eco47I_1", "Eco47I_2"),
    cache_root = cache_root,
    rebase_data = rebase_data[c(2, 1), , drop = FALSE]
  )
  expect_identical(metadata_reordered$genbank_accession, c("CP100001", "CP100002"))

  hits <- .rebasefinder_test_hits("reordered_duplicate", "Type II", "R")
  hits$hit_label <- "Eco47I_2"
  blast <- data.frame(
    query_id = "reordered_duplicate",
    rebase_enzyme = "Eco47I_2",
    pct_identity = 0.90,
    evalue = 1e-40,
    bitscore = 180,
    length = 80,
    stringsAsFactors = FALSE
  )
  enriched <- DNMB:::.dnmb_rebasefinder_enrich_blast_metadata(
    hits,
    blast_tbl = blast,
    rebase_data = rebase_data,
    cache_root = cache_root
  )

  expect_identical(enriched$reference_rec_seq, "GATC")
  expect_identical(enriched$rebase_reference_accession, "CP100001")
  expect_identical(enriched$rebase_reference_locus, "eco_a")
  expect_identical(enriched$recognition_match_reference_protein_id, "ECO_A")
  expect_identical(enriched$rebase_match_pacbio_organism, "Raw A strain")
})

test_that("missing REBASE header names are negatively cached", {
  cache_root <- tempfile("rebasefinder-negative-header-cache-")
  cache_dir <- DNMB:::.dnmb_rebasefinder_cache_dir(cache_root, create = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(
    c(
      ">REBASE:M.Present\tEnzType:Type II methyltransferase\tRecSeq:GATC",
      strrep("A", 80)
    ),
    file.path(cache_dir, "REBASE_protein_seqs.txt")
  )

  metadata <- DNMB:::.dnmb_rebasefinder_reference_header_metadata(
    "M.Absent",
    cache_root = cache_root,
    rebase_data = data.frame(
      enzyme_name = "M.Absent",
      sequence = strrep("C", 80),
      stringsAsFactors = FALSE
    )
  )

  expect_identical(metadata$enzyme_name, "M.Absent")
  expect_true(is.na(metadata$genbank_accession))
  cache <- readRDS(file.path(cache_dir, "rebase_reference_header_metadata.rds"))$data
  absent <- cache[cache$enzyme_name == "M.Absent", , drop = FALSE]
  expect_equal(nrow(absent), 1L)
  expect_identical(as.integer(absent$.raw_occurrence), 0L)
})

test_that("resolved unsuffixed duplicate metadata follows the selected REBASE record", {
  cache_root <- tempfile("rebasefinder-resolved-duplicate-")
  cache_dir <- DNMB:::.dnmb_rebasefinder_cache_dir(cache_root, create = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(
    c(
      ">REBASE:M.Duplicate\tEnzType:Type II methyltransferase\tRecSeq:GATC\tGenBank:CP000101\tLocus:dup_a\tProteinId:DUP_A",
      "AAAA",
      ">REBASE:M.Duplicate\tEnzType:Type III methyltransferase\tRecSeq:CCGG\tGenBank:CP000102\tLocus:dup_b\tProteinId:DUP_B",
      "AAAA"
    ),
    file.path(cache_dir, "REBASE_protein_seqs.txt")
  )
  rebase_data <- data.frame(
    enzyme_name = c("M.Duplicate", "M.Duplicate"),
    enz_type = c("Type II methyltransferase", "Type III methyltransferase"),
    rec_seq = c("GATC", "CCGG"),
    sequence = c(strrep("A", 300), strrep("A", 300)),
    seq_length = 300,
    is_gold_standard = TRUE,
    rm_type = c("Type_II", "Type_III"),
    subunit = "M",
    stringsAsFactors = FALSE
  )
  saveRDS(rebase_data, file.path(cache_dir, "rebase_data.rds"))
  bairoch <- data.frame(
    enzyme_name = character(), rec_seq = character(), meth_type = character(),
    meth_pos = character(), meth_all = character(), stringsAsFactors = FALSE
  )
  attr(bairoch, "dnmb_bairoch_schema") <- 2L
  saveRDS(bairoch, file.path(cache_dir, "rebase_bairoch_lookup.rds"))
  saveRDS(
    data.frame(
      nucleotide_accession = "CP000102",
      organism = "Duplicate B strain",
      rebase_org_id = "102",
      pacbio_url = "https://rebase.neb.com/cgi-bin/pacbioget?102",
      .accession_key = "CP000102",
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "rebase_pacbio_index.rds")
  )

  hits <- .rebasefinder_test_hits("dup_query", "Type III", "M")
  hits$hit_label <- "M.Duplicate"
  hits$curated_blast_subject_id <- "M.Duplicate"
  hits$curated_blast_resolved_subject_id <- "M.Duplicate_2"
  hits$curated_blast_reference_index <- 2L
  hits$curated_blast_promoted <- TRUE
  hits$query_aa_len <- 300
  blast <- data.frame(
    query_id = "dup_query",
    rebase_enzyme = "M.Duplicate",
    pct_identity = 0.85,
    length = 300,
    evalue = 1e-80,
    bitscore = 500,
    stringsAsFactors = FALSE
  )

  enriched <- DNMB:::.dnmb_rebasefinder_enrich_blast_metadata(
    hits,
    blast_tbl = blast,
    rebase_data = rebase_data,
    cache_root = cache_root
  )

  expect_identical(enriched$reference_rec_seq, "CCGG")
  expect_identical(enriched$recognition_match_subject_id, "M.Duplicate_2")
  expect_identical(enriched$rebase_reference_accession, "CP000102")
  expect_identical(enriched$recognition_match_reference_accession, "CP000102")
  expect_identical(enriched$rebase_match_pacbio_organism, "Duplicate B strain")
})

test_that("query GenBank accession links the analyzed strain to PacBio metadata", {
  cache_root <- tempfile("rebasefinder-query-pacbio-")
  cache_dir <- DNMB:::.dnmb_rebasefinder_cache_dir(cache_root, create = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)
  genbank <- tempfile(fileext = ".gbff")
  writeLines(c(
    "LOCUS       TEST 100 bp DNA",
    "ACCESSION   CP123456",
    "VERSION     CP123456.2",
    "//"
  ), genbank)
  saveRDS(
    data.frame(
      nucleotide_accession = "CP123456",
      organism = "Query strain",
      rebase_org_id = "777",
      pacbio_url = "https://rebase.neb.com/cgi-bin/pacbioget?777",
      .accession_key = "CP123456",
      stringsAsFactors = FALSE
    ),
    file.path(cache_dir, "rebase_pacbio_index.rds")
  )
  bairoch <- data.frame(
    enzyme_name = character(), rec_seq = character(), meth_type = character(),
    meth_pos = character(), meth_all = character(), stringsAsFactors = FALSE
  )
  attr(bairoch, "dnmb_bairoch_schema") <- 2L
  saveRDS(bairoch, file.path(cache_dir, "rebase_bairoch_lookup.rds"))

  enriched <- DNMB:::.dnmb_rebasefinder_enrich_blast_metadata(
    .rebasefinder_test_hits("query1"),
    blast_tbl = data.frame(),
    rebase_data = data.frame(),
    cache_root = cache_root,
    query_genbank = genbank
  )

  expect_identical(enriched$query_strain_accessions, "CP123456.2")
  expect_true(enriched$query_strain_pacbio_available)
  expect_identical(enriched$query_strain_pacbio_organism, "Query strain")
  expect_identical(enriched$query_strain_pacbio_link, "https://rebase.neb.com/cgi-bin/pacbioget?777")
  expect_identical(enriched$query_strain_pacbio_index_link, "https://rebase.neb.com/cgi-bin/pblist")
})

test_that("blank GenBank ACCESSION and VERSION headers do not become accessions", {
  genbank <- tempfile(fileext = ".gbk")
  on.exit(unlink(genbank, force = TRUE), add = TRUE)
  writeLines(c(
    "LOCUS       CONTIG_1 100 bp DNA",
    "ACCESSION   ",
    "VERSION     ",
    "//",
    "LOCUS       CONTIG_2 100 bp DNA",
    "ACCESSION   ",
    "//"
  ), genbank)

  accessions <- DNMB:::.dnmb_rebasefinder_genbank_accessions(genbank)

  expect_identical(accessions, character())
})

test_that("REBASE module details omit literal NA values and retain clickable links", {
  tbl <- data.frame(
    locus_tag = "gene1",
    product = "restriction enzyme",
    REBASEfinder_family_id = "Type II",
    REBASEfinder_hit_label = "R.TestSystem",
    REBASEfinder_rec_seq = NA_character_,
    REBASEfinder_reference_rec_seq = "GATC",
    REBASEfinder_recognition_match = "R.TestSystem",
    REBASEfinder_recognition_donor = "M.TestSystem",
    REBASEfinder_rebase_match_link = "https://rebase.neb.com/rebase/enz/R.TestSystem.html",
    REBASEfinder_recognition_match_rebase_link = "https://rebase.neb.com/rebase/enz/R.TestSystem.html",
    REBASEfinder_recognition_donor_rebase_link = "https://rebase.neb.com/rebase/enz/M.TestSystem.html",
    REBASEfinder_rebase_reference_ncbi_link = "https://www.ncbi.nlm.nih.gov/nuccore/CP000010",
    REBASEfinder_recognition_match_reference_ncbi_link = "https://www.ncbi.nlm.nih.gov/nuccore/CP000010",
    REBASEfinder_recognition_reference_ncbi_link = "https://www.ncbi.nlm.nih.gov/nuccore/CP000011",
    REBASEfinder_recognition_donor_ncbi_link = "https://www.ncbi.nlm.nih.gov/nuccore/CP000011",
    REBASEfinder_recognition_match_pacbio_link = "https://rebase.neb.com/cgi-bin/pacbioget?10",
    REBASEfinder_recognition_donor_pacbio_link = "https://rebase.neb.com/cgi-bin/pacbioget?11",
    REBASEfinder_query_strain_pacbio_link = "https://rebase.neb.com/cgi-bin/pacbioget?12",
    REBASEfinder_query_strain_pacbio_index_link = "https://rebase.neb.com/cgi-bin/pblist",
    stringsAsFactors = FALSE
  )
  details <- DNMB:::dnmb_build_module_details_table(tbl)

  expect_false(any(grepl("=NA(?:;|$)", details$reference_id, perl = TRUE)))
  expect_false(any(grepl("=NA(?:;|$)", details$resource_links, perl = TRUE)))
  expect_true(all(c(
    "rebase_match_link", "recognition_match_rebase_link",
    "recognition_donor_rebase_link", "rebase_reference_ncbi_link",
    "recognition_match_reference_ncbi_link", "recognition_reference_ncbi_link",
    "recognition_donor_ncbi_link", "recognition_match_pacbio_link",
    "recognition_donor_pacbio_link", "query_strain_pacbio_link",
    "query_strain_pacbio_index_link"
  ) %in% names(details)))
  expect_identical(
    details$recognition_match_pacbio_link,
    "https://rebase.neb.com/cgi-bin/pacbioget?10"
  )
  expect_identical(
    details$recognition_donor_pacbio_link,
    "https://rebase.neb.com/cgi-bin/pacbioget?11"
  )
  expect_identical(
    details$query_strain_pacbio_link,
    "https://rebase.neb.com/cgi-bin/pacbioget?12"
  )
  expect_identical(
    details$query_strain_pacbio_index_link,
    "https://rebase.neb.com/cgi-bin/pblist"
  )
  expect_identical(
    DNMB:::dnmb_excel_hyperlink_label("REBASEfinder_rebase_reference_ncbi_link"),
    "NCBI"
  )
})

test_that("Excel hyperlink formulas use openxlsx formula syntax", {
  table <- data.frame(
    locus_tag = "gene1",
    pacbio_link = "https://rebase.neb.com/cgi-bin/pacbioget?10",
    stringsAsFactors = FALSE
  )
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "links")
  openxlsx::writeData(wb, "links", table)
  DNMB:::dnmb_write_hyperlink_columns(wb, "links", table)
  xlsx <- tempfile(fileext = ".xlsx")
  extract_dir <- tempfile("xlsx-xml-")
  on.exit(unlink(c(xlsx, extract_dir), recursive = TRUE, force = TRUE), add = TRUE)
  openxlsx::saveWorkbook(wb, xlsx, overwrite = TRUE)
  utils::unzip(xlsx, files = "xl/worksheets/sheet1.xml", exdir = extract_dir)
  xml <- paste(readLines(file.path(extract_dir, "xl/worksheets/sheet1.xml"), warn = FALSE), collapse = "")
  expect_match(xml, "<f>HYPERLINK\\(")
  expect_false(grepl("<f>=HYPERLINK", xml, fixed = TRUE))
})

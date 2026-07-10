test_that("REBASEfinder structure parser uses curated chain roles", {
  tsv <- tempfile(fileext = ".tsv")
  writeLines(c(
    "query\ttarget\tprob\tevalue\tbits\talntmscore\trmsd\tlddt",
    "Qrease\tEcoRI_TypeIIP_1ERI__1ERI_A\t0.9\t1e-20\t200\t0.7\t2.0\t0.8",
    "Qfok\tFokI_TypeIIS_1FOK__1FOK_A\t0.9\t1e-20\t200\t0.7\t2.0\t0.8",
    "Qmtase\tHhaI_C5MTase_2HMY__2HMY_B\t0.9\t1e-20\t200\t0.7\t2.0\t0.8",
    "Qdna\tEcoRI_TypeIIP_1ERI__1ERI_B\t0.9\t1e-20\t200\t0.7\t2.0\t0.8"
  ), tsv)

  parsed <- DNMB:::.dnmb_rebasefinder_read_structure_validation(tsv)
  rownames(parsed) <- parsed$query

  expect_equal(parsed["Qrease", "structure_family"], "Type II")
  expect_equal(parsed["Qrease", "structure_role"], "R")
  expect_equal(parsed["Qfok", "structure_family"], "Type II")
  expect_equal(parsed["Qfok", "structure_family_raw"], "Type IIS")
  expect_equal(parsed["Qfok", "structure_role"], "R")
  expect_equal(parsed["Qmtase", "structure_role"], "M")
  expect_equal(parsed["Qdna", "structure_status"], "structure_excluded")
  expect_false(parsed["Qdna", "structure_pass"])
})

test_that("structure family and fused-role aliases are canonicalized", {
  expect_true(DNMB:::.dnmb_rebasefinder_structure_family_compatible("Type I", "Type ISP"))
  expect_true(DNMB:::.dnmb_rebasefinder_structure_family_compatible("Type II", "Type IIG"))
  expect_true(DNMB:::.dnmb_rebasefinder_structure_family_compatible("Type II", "Type IIL"))
  expect_true(DNMB:::.dnmb_rebasefinder_structure_role_compatible("RM", "R/M"))
  expect_false(DNMB:::.dnmb_rebasefinder_structure_role_compatible("R", "R/M"))
  expect_identical(DNMB:::.dnmb_rebasefinder_canonical_structure_role(c("R/M", "M/R")), c("RM", "RM"))
  partial_row <- data.frame(
    REBASEfinder_structure_status = "structure_supported_partial_model",
    REBASEfinder_enzyme_role = "RM",
    REBASEfinder_structure_role = "RM",
    stringsAsFactors = FALSE
  )
  expect_identical(
    DNMB:::.dnmb_rebasefinder_structure_motif_status(partial_row),
    "structure_supported_role_consistent"
  )
})

test_that("Foldseek validation retains alignment mapping and rejects tiny local matches", {
  tsv <- tempfile(fileext = ".tsv")
  writeLines(c(
    paste(
      c("query", "target", "prob", "evalue", "bits", "alntmscore", "qcov", "tcov",
        "qstart", "qend", "tstart", "tend", "qaln", "taln"),
      collapse = "\t"
    ),
    paste(
      c("Qglobal", "HhaI_C5MTase_2HMY__2HMY_B", "0.9", "1e-20", "200", "0.7",
        "0.8", "0.9", "1", "300", "5", "304", "ACD-EF", "ACDGEF"),
      collapse = "\t"
    ),
    paste(
      c("Qtiny", "HhaI_C5MTase_2HMY__2HMY_B", "0.99", "1e-30", "250", "0.8",
        "0.1", "0.9", "50", "80", "5", "35", "ACDEF", "ACDEF"),
      collapse = "\t"
    )
  ), tsv)

  parsed <- DNMB:::.dnmb_rebasefinder_read_structure_validation(tsv)
  rownames(parsed) <- parsed$query

  expect_true(parsed["Qglobal", "structure_pass"])
  expect_equal(parsed["Qglobal", "structure_query_coverage"], 0.8)
  expect_identical(parsed["Qglobal", "structure_query_alignment"], "ACD-EF")
  expect_false(parsed["Qtiny", "structure_pass"])
  expect_identical(parsed["Qtiny", "structure_status"], "structure_low_coverage")
})

test_that("REBASEfinder structure support promotes weak hits and adds structure-only hits", {
  hits <- data.frame(
    query = "Q1",
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = NA_character_,
    hit_label = "weak_hit",
    enzyme_role = NA_character_,
    evidence_mode = "annotation_only",
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  genes <- data.frame(
    locus_tag = c("Q1", "Qnew"),
    translation = c("MAAAAA", "MBBBBB"),
    stringsAsFactors = FALSE
  )
  structures <- data.frame(
    query = c("Q1", "Qnew"),
    structure_hit = c("HhaI_C5MTase_2HMY__2HMY_B", "EcoRI_TypeIIP_1ERI__1ERI_A"),
    structure_evalue = c(1e-20, 1e-30),
    structure_bitscore = c(200, 300),
    structure_probability = c(0.9, 0.95),
    structure_tmscore = c(0.7, 0.8),
    structure_rmsd = c(2, 1),
    structure_reference_id = c("HhaI_C5MTase_2HMY", "EcoRI_TypeIIP_1ERI"),
    structure_chain = c("B", "A"),
    structure_family = c("Type II", "Type II"),
    structure_role = c("M", "R"),
    structure_class = c("Type II C5 cytosine methyltransferase", "Type IIP restriction endonuclease-DNA complex"),
    structure_chain_role = c("M", "R"),
    structure_pass = c(TRUE, TRUE),
    structure_status = c("structure_supported", "structure_supported"),
    stringsAsFactors = FALSE
  )

  merged <- DNMB:::.dnmb_rebasefinder_merge_structure_validation(hits, genes, structures)
  rownames(merged) <- merged$query

  expect_equal(merged["Q1", "evidence_mode"], "structure_supported")
  expect_true(merged["Q1", "typing_eligible"])
  expect_equal(merged["Q1", "family_id"], "Type II")
  expect_equal(merged["Q1", "enzyme_role"], "M")
  expect_equal(merged["Qnew", "evidence_mode"], "structure_only")
  expect_equal(merged["Qnew", "enzyme_role"], "R")
})

test_that("REBASEfinder structure query FASTA includes every candidate with a translation", {
  out_dir <- tempfile("rebase_queries_")
  dir.create(out_dir)
  genes <- data.frame(
    locus_tag = c("blast1", "anno1", "typeiii1", "struct1"),
    translation = rep("MPEPTIDE", 4),
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query = c("blast1", "anno1", "typeiii1", "struct1"),
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = c("Type II", NA, "Type III", "Type II"),
    hit_label = c("M.Test", NA, "M.TypeIII", "R.Struct"),
    enzyme_role = c("M", NA, "M", "R"),
    evidence_mode = c("low_confidence", "annotation_only", "annotation_only", "structure_only"),
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = c(FALSE, FALSE, FALSE, TRUE),
    blast_identity = c(0.31, NA, NA, NA),
    stringsAsFactors = FALSE
  )

  fasta <- DNMB:::.dnmb_rebasefinder_write_structure_queries(out_dir, genes, hits)
  headers <- readLines(fasta, warn = FALSE)
  headers <- sub("^>", "", headers[grepl("^>", headers)])

  expect_true(any(grepl("^blast1\\b", headers)))
  expect_true(any(grepl("^anno1\\b", headers)))
  expect_true(any(grepl("^typeiii1\\b", headers)))
  expect_true(any(grepl("^struct1\\b", headers)))
  expect_true(any(grepl("^anno1\\b.*priority=annotation_candidate", headers)))
})

test_that("REBASEfinder structure coverage reports missing and unchecked Foldseek queries", {
  out_dir <- tempfile("rebase_coverage_")
  dir.create(out_dir)
  dir.create(file.path(out_dir, "query_structures"))
  genes <- data.frame(
    locus_tag = c("q_supported", "q_unchecked", "q_missing"),
    start = c(1, 101, 201),
    end = c(80, 180, 280),
    translation = rep("MPEPTIDE", 3),
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query = genes$locus_tag,
    family_id = c("Type III", "Type I", "Type IV"),
    enzyme_role = c("M", "S", "R"),
    evidence_mode = c("structure_supported", "operon_context", "operon_context"),
    hit_label = c("M.Test", "S.Test", "Mrr.Test"),
    blast_identity = c(NA, NA, NA),
    stringsAsFactors = FALSE
  )
  fasta <- DNMB:::.dnmb_rebasefinder_write_structure_queries(out_dir, genes, hits)
  writeLines("MODEL\nEND", file.path(out_dir, "query_structures", "q_supported.pdb"))
  writeLines("MODEL\nEND", file.path(out_dir, "query_structures", "q_unchecked.pdb"))
  foldseek <- file.path(out_dir, "foldseek_results.tsv")
  writeLines(c(
    "query\ttarget\tprob\tevalue\tbits\talntmscore\trmsd\tlddt",
    "q_supported\tEcoP15I_ModRes_4ZCF__4ZCF_A\t1.0\t1e-20\t200\t0.7\t2\t0.8"
  ), foldseek)
  structure_tbl <- data.frame(
    query = "q_supported",
    structure_hit = "EcoP15I_ModRes_4ZCF__4ZCF_A",
    structure_pass = TRUE,
    structure_status = "structure_supported",
    stringsAsFactors = FALSE
  )

  cov <- DNMB:::.dnmb_rebasefinder_write_structure_coverage(
    out_dir, genes, hits,
    fasta_path = fasta,
    structure_validation_path = foldseek,
    structure_tbl = structure_tbl
  )
  tbl <- read.delim(cov$tsv, stringsAsFactors = FALSE, check.names = FALSE)

  expect_equal(cov$n_queries, 3)
  expect_equal(cov$n_structure_files, 2)
  expect_equal(cov$n_foldseek_hits, 1)
  expect_equal(tbl$coverage_status[match("q_supported", tbl$query)], "foldseek_supported")
  expect_equal(tbl$coverage_status[match("q_unchecked", tbl$query)], "structure_available_no_foldseek_hit")
  expect_equal(tbl$coverage_status[match("q_missing", tbl$query)], "structure_missing")
  expect_true(file.exists(cov$missing_faa))
  expect_true(any(grepl("^>q_missing\\b", readLines(cov$missing_faa, warn = FALSE))))
})

test_that("partial structure coverage maps safely onto a larger hit table", {
  hits <- data.frame(
    query = c("covered", "not_covered", "also_not_covered"),
    structure_status = NA_character_,
    stringsAsFactors = FALSE
  )
  coverage <- data.frame(
    query = "covered",
    coverage_status = "structure_missing",
    structure_file = NA_character_,
    structure_file_exists = FALSE,
    foldseek_hit_present = FALSE,
    structure_family_consistent = NA,
    structure_role_consistent = NA,
    structure_candidate_consistent = NA,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_apply_structure_coverage(hits, coverage)

  expect_equal(nrow(out), 3L)
  expect_identical(out$structure_status[[1]], "structure_missing")
  expect_true(all(c(
    "structure_family_consistent", "structure_role_consistent",
    "structure_candidate_consistent"
  ) %in% names(out)))
  expect_true(all(is.na(out$structure_candidate_consistent)))
})

test_that("partial query models are labeled separately from Foldseek alignment coverage", {
  hits <- data.frame(
    query = "partial_model",
    support = "fold=positive; qcov=0.99",
    structure_status = "structure_supported",
    stringsAsFactors = FALSE
  )
  coverage <- data.frame(
    query = "partial_model",
    coverage_status = "foldseek_supported_partial_model",
    structure_file = "partial_model.pdb",
    structure_file_exists = TRUE,
    structure_model_coverage = 0.54,
    structure_model_scope = "partial_model",
    foldseek_hit_present = TRUE,
    structure_family = "Type II",
    structure_role = "RM",
    structure_chain_role = "RM",
    structure_family_consistent = TRUE,
    structure_role_consistent = TRUE,
    structure_candidate_consistent = TRUE,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_apply_structure_coverage(hits, coverage)

  expect_true(out$structure_supported)
  expect_true(out$structure_pass)
  expect_identical(out$structure_status, "structure_supported_partial_model")
  expect_equal(out$structure_model_coverage, 0.54)
  expect_identical(out$structure_chain_role, "RM")
  expect_match(out$support, "model_cov_fullseq=0.540")
  expect_match(out$support, "structure_scope=partial_model")
})

test_that("Type III context rescues nearby missing Res partner", {
  genes <- data.frame(
    locus_tag = c("mod1", "res1", "far1"),
    contig = "ctg1",
    start = c(100, 900, 10000),
    end = c(700, 1800, 11000),
    direction = "+",
    product = c("Type III modification methylase Mod subunit",
                "Type III restriction enzyme Res subunit helicase",
                "hypothetical protein"),
    translation = c("MAAAAA", "MBBBBB", "MCCCCC"),
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query = "mod1",
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = "Type III",
    hit_label = "M.EcoP15I",
    enzyme_role = "M",
    evidence_mode = "high_confidence",
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = TRUE,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_add_typeiii_context(
    hits,
    genes,
    max_operon_gap = 5000,
    max_intervening = 1L,
    max_neighbors = 2L
  )
  rownames(out) <- out$query

  expect_equal(out["mod1", "typeiii_context_status"], "complete_mod_res")
  expect_true(out["mod1", "typeiii_operon_supported"])
  expect_true("res1" %in% out$query)
  expect_equal(out["res1", "evidence_mode"], "operon_context")
  expect_equal(out["res1", "enzyme_role"], "R")
  expect_true(out["res1", "typeiii_operon_supported"])
})

test_that("a weak cross-family Mod hit is reclassified by a coherent Type III Mod-Res operon", {
  genes <- data.frame(
    locus_tag = c("mod_weak_typeii", "res_typeiii"),
    contig = "ctg1",
    start = c(100, 820),
    end = c(780, 2600),
    direction = "+",
    product = c("DNA methyltransferase", "DNA methyltransferase"),
    `Signature accession_Pfam` = c("PF01555", "PF04851"),
    `Signature description_Pfam` = c(
      "DNA methylase",
      "Type III restriction enzyme, res subunit"
    ),
    translation = c(strrep("A", 220), strrep("A", 800)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  hits <- data.frame(
    query = genes$locus_tag,
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = c("Type II", "Type III"),
    hit_label = c("M.Weak", "R.ResIII"),
    enzyme_role = c("M", "R"),
    evidence_mode = "high_confidence",
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = FALSE,
    blast_identity = c(0.59, 0.54),
    blast_alignment_quality = "weak",
    blast_role_compatible = TRUE,
    blast_family_compatible = TRUE,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_add_typeiii_context(hits, genes)
  rownames(out) <- out$query

  expect_true(out["mod_weak_typeii", "typeiii_operon_supported"])
  expect_true(out["res_typeiii", "typeiii_operon_supported"])
  expect_identical(out["mod_weak_typeii", "family_id"], "Type III")
  expect_identical(out["mod_weak_typeii", "typeiii_context_previous_family"], "Type II")
  expect_true(out["mod_weak_typeii", "typeiii_context_reclassified"])
})

test_that("Type III context rescues motif-inferred nearby missing Res partner", {
  res_seq <- paste0(
    paste(rep("A", 20), collapse = ""),
    "AGGGGGKS",
    paste(rep("A", 70), collapse = ""),
    "DEAH",
    paste(rep("A", 80), collapse = ""),
    "SAT",
    paste(rep("A", 417), collapse = ""),
    "PDAAAAADEK",
    paste(rep("A", 40), collapse = "")
  )
  genes <- data.frame(
    locus_tag = c("mod1", "res1"),
    contig = "ctg1",
    start = c(100, 900),
    end = c(700, 1800),
    direction = "+",
    product = c("Type III modification methylase Mod subunit",
                "hypothetical protein"),
    translation = c("MFDGTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADPPYAAAA", res_seq),
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query = "mod1",
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = "Type III",
    hit_label = "M.EcoP15I",
    enzyme_role = "M",
    evidence_mode = "high_confidence",
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = TRUE,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_add_typeiii_context(
    hits,
    genes,
    max_operon_gap = 5000,
    max_intervening = 1L,
    max_neighbors = 2L
  )
  rownames(out) <- out$query

  expect_equal(out["mod1", "typeiii_context_status"], "complete_mod_res")
  expect_true(out["mod1", "typeiii_operon_supported"])
  expect_true("res1" %in% out$query)
  expect_equal(out["res1", "evidence_mode"], "operon_context")
  expect_equal(out["res1", "enzyme_role"], "R")
  expect_true(out["res1", "typeiii_operon_supported"])
  expect_match(out["res1", "support"], "ResIII-WA\\+ResIII-WB\\+ResIII-MIII\\+ResIII-PD")
})

test_that("Type III context does not rescue a DEAD-only helicase", {
  helicase_seq <- paste0(
    paste(rep("A", 20), collapse = ""),
    "AGGGGGKS",
    paste(rep("A", 70), collapse = ""),
    "DEAH",
    paste(rep("A", 500), collapse = ""),
    "PDAAAAADEK",
    paste(rep("A", 40), collapse = "")
  )
  genes <- data.frame(
    locus_tag = c("res1", "mod1"),
    contig = "ctg1",
    start = c(100, 1200),
    end = c(1000, 2400),
    direction = "+",
    product = c("ATP-dependent RNA helicase DbpA", "Type III modification methylase Mod subunit"),
    translation = c(helicase_seq, "MFDGTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADPPYAAAA"),
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query = "mod1",
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = "Type III",
    hit_label = "M.TypeIII",
    enzyme_role = "M",
    evidence_mode = "high_confidence",
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = TRUE,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_add_typeiii_context(
    hits,
    genes,
    max_operon_gap = 5000,
    max_intervening = 1L,
    max_neighbors = 2L
  )
  rownames(out) <- out$query

  expect_equal(out["mod1", "typeiii_context_status"], "missing_res_check_neighbors")
  expect_false("res1" %in% out$query)
})

test_that("Type I context rescues HsdR and HsdS around a Type I methylase", {
  genes <- data.frame(
    locus_tag = c("hsdR1", "hsdM1", "hsdS1"),
    contig = "ctg1",
    start = c(100, 1200, 2400),
    end = c(1000, 2200, 3400),
    direction = "+",
    gene = c(NA, NA, NA),
    product = c("type I restriction endonuclease subunit R",
                "class I SAM-dependent DNA methyltransferase",
                "restriction endonuclease subunit S"),
    translation = c(
      paste0(paste(rep("A", 300), collapse = ""), "DEAH", paste(rep("A", 80), collapse = ""), "AGGGGGKS"),
      "MFDGTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADPPYAAAA",
      "MAAAAA"
    ),
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query = "hsdM1",
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = "Type I",
    hit_label = "M.TypeI",
    enzyme_role = "M",
    evidence_mode = "high_confidence",
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = TRUE,
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_add_typei_context(
    hits,
    genes,
    max_operon_gap = 5000,
    max_intervening = 1L,
    max_neighbors = 3L
  )
  rownames(out) <- out$query

  expect_equal(out["hsdM1", "typei_context_status"], "complete_mrs")
  expect_true(out["hsdM1", "typei_operon_supported"])
  expect_equal(out["hsdR1", "enzyme_role"], "R")
  expect_equal(out["hsdS1", "enzyme_role"], "S")
  expect_true(out["hsdR1", "typei_operon_supported"])
  expect_true(out["hsdS1", "typei_operon_supported"])
})

test_that("Type IV candidates are rescued from Mrr-like restriction annotations", {
  seq <- paste0("M", "D", paste(rep("A", 8), collapse = ""), "EAK", paste(rep("A", 30), collapse = ""), "DK")
  genes <- data.frame(
    locus_tag = "typeiv1",
    contig = "ctg1",
    start = 100,
    end = 700,
    direction = "+",
    product = "Mrr_cat PF04471 modification-dependent restriction endonuclease",
    translation = seq,
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query = character(),
    source = character(),
    family_system = character(),
    family_id = character(),
    hit_label = character(),
    enzyme_role = character(),
    evidence_mode = character(),
    substrate_label = character(),
    support = character(),
    typing_eligible = logical(),
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_add_typeiv_candidates(hits, genes)

  expect_true("typeiv1" %in% out$query)
  expect_equal(out$family_id[match("typeiv1", out$query)], "Type IV")
  expect_equal(out$enzyme_role[match("typeiv1", out$query)], "R")
})

test_that("REBASE overview motif scan preserves Type III Res motif columns", {
  seq <- paste0(
    paste(rep("A", 20), collapse = ""),
    "AGGGGGKS",
	    paste(rep("A", 70), collapse = ""),
	    "DEAH",
	    paste(rep("A", 30), collapse = ""),
	    "SAT",
	    paste(rep("A", 87), collapse = ""),
	    "PDAAAAADEK",
    paste(rep("A", 40), collapse = "")
  )
  tbl <- data.frame(locus_tag = "res1", translation = seq, stringsAsFactors = FALSE)
  motifs <- DNMB:::.dnmb_rebasefinder_scan_motifs(tbl)

  expect_true(all(c("ResIII-WA", "ResIII-WB", "ResIII-MIII", "ResIII-PD") %in% names(motifs)))
  expect_match(motifs[["ResIII-WA"]][[1]], "^present")
  expect_match(motifs[["ResIII-WB"]][[1]], "^present")
  expect_match(motifs[["ResIII-MIII"]][[1]], "^present")
  expect_match(motifs[["ResIII-PD"]][[1]], "^present")
})

test_that("an isolated Type I HsdR motif III remains raw and does not verify", {
  seq <- paste0(strrep("A", 260), "TAT", strrep("A", 220))
  tbl <- data.frame(
    locus_tag = "hsdr1",
    REBASEfinder_family_id = "Type I",
    REBASEfinder_enzyme_role = "R",
    translation = seq,
    stringsAsFactors = FALSE
  )
  motifs <- DNMB:::.dnmb_rebasefinder_scan_motifs(tbl)
  hits <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl)
  raw <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl, role_relevant_only = FALSE)

  expect_match(motifs[["HsdR-MIII"]][[1]], "^present")
  expect_equal(nrow(hits), 0L)
  expect_true(any(raw$motif == "HsdR-MIII" & raw$evidence_level == "raw_unqualified"))
  expect_false(DNMB:::.dnmb_rebasefinder_motif_verified(tbl)[[1]])
})

test_that("REBASE overview motif scan separates canonical SAM and catalytic MTase signals", {
  seq <- paste0(
    paste(rep("A", 20), collapse = ""),
    "FDGTG",
    paste(rep("A", 90), collapse = ""),
    "DPPY",
    paste(rep("A", 180), collapse = "")
  )
  tbl <- data.frame(
    locus_tag = "mod1",
    REBASEfinder_enzyme_role = "M",
    translation = seq,
    stringsAsFactors = FALSE
  )
  motifs <- DNMB:::.dnmb_rebasefinder_scan_motifs(tbl)

  expect_match(motifs[["SAM"]][[1]], "^present")
  expect_match(motifs[["Amino-IV"]][[1]], "^present")
  expect_equal(DNMB:::.dnmb_rebasefinder_typeiii_sequence_role(seq), "M")
  loose <- paste0(
    paste(rep("A", 40), collapse = ""),
    "KDGNG",
    paste(rep("A", 70), collapse = ""),
    "DPPY"
  )
  loose_tbl <- data.frame(
    locus_tag = "mod2",
    REBASEfinder_enzyme_role = "M",
    translation = loose,
    stringsAsFactors = FALSE
  )
  loose_motifs <- DNMB:::.dnmb_rebasefinder_scan_motifs(loose_tbl)
  expect_match(loose_motifs[["SAM"]][[1]], "^absent")
  expect_match(loose_motifs[["Amino-IV"]][[1]], "^present")

  manual_like <- paste0(strrep("A", 26), "FAGSA", strrep("A", 85), "DPPY", strrep("A", 40))
  manual_tbl <- data.frame(
    locus_tag = "mod3",
    REBASEfinder_enzyme_role = "M",
    translation = manual_like,
    stringsAsFactors = FALSE
  )
  manual_motifs <- DNMB:::.dnmb_rebasefinder_scan_motifs(manual_tbl)
  expect_match(manual_motifs[["SAM"]][[1]], "^present")
  expect_match(manual_motifs[["Amino-IV"]][[1]], "^present")
  expect_equal(DNMB:::.dnmb_rebasefinder_typeiii_sequence_role(manual_like), "M")
})

test_that("REBASE motif hit table lists motifs with partial and structure status", {
  tbl <- data.frame(
    locus_tag = "mod1",
    product = "DNA methyltransferase",
    translation = paste0(strrep("A", 25), "FDGTG", strrep("A", 90), "DPPY", strrep("A", 20)),
    REBASEfinder_family_id = "Type III",
    REBASEfinder_enzyme_role = "M",
    REBASEfinder_hit_label = "M.Example",
    REBASEfinder_structure_status = "structure_supported",
    REBASEfinder_structure_reference_id = "EcoP15I_ModRes_4ZCF",
    REBASEfinder_structure_role = "M",
    REBASEfinder_structure_chain_role = "M",
    stringsAsFactors = FALSE
  )
  hits <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl)

  expect_true(all(c("SAM", "Amino-IV") %in% hits$motif))
  expect_false(any(hits$motif == "SAM-like"))
  expect_true(all(hits$partial_status == "partial_or_short"))
  expect_true(any(hits$structural_verification == "structure_supported_role_consistent"))
})

test_that("REBASE motif hit table keeps non-role raw hits out of the default list", {
  tbl <- data.frame(
    locus_tag = "spec1",
    translation = paste0(strrep("A", 20), "FDGTG", strrep("A", 100)),
    REBASEfinder_family_id = "Type I",
    REBASEfinder_enzyme_role = "R",
    stringsAsFactors = FALSE
  )
  hits <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl)
  raw <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl, role_relevant_only = FALSE)

  expect_equal(nrow(hits), 0)
  expect_true(any(raw$motif == "SAM"))
  expect_true(all(!raw$role_relevant))
})

test_that("Type I specificity subunits use dual TRD architecture instead of a short SAM-like pattern", {
  tbl <- data.frame(
    locus_tag = "hsds1",
    translation = paste0(strrep("A", 80), "FTGTA", strrep("A", 360)),
    REBASEfinder_family_id = "Type I",
    REBASEfinder_enzyme_role = "S",
    `Signature accession_Pfam` = "PF01420, PF01420",
    `Signature description_Pfam` = paste(
      "Type I restriction modification DNA specificity domain",
      "Type I restriction modification DNA specificity domain",
      sep = ", "
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  motifs <- DNMB:::.dnmb_rebasefinder_scan_motifs(tbl)
  hits <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl)
  raw <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl, role_relevant_only = FALSE)

  expect_false("HsdS-FxGxA" %in% names(motifs))
  expect_true(any(hits$motif == "HsdS-2TRD"))
  expect_true(any(hits$evidence_level == "architecture_supported"))
  expect_true(any(raw$motif == "SAM" & !raw$role_relevant))
  expect_true(DNMB:::.dnmb_rebasefinder_motif_verified(tbl)[[1]])
})

test_that("C5 DNA MTases use the ordered PCQ and ENV catalytic pair without requiring a loose SAM hit", {
  seq <- paste0(strrep("A", 116), "PCQ", strrep("A", 43), "ENV", strrep("A", 322))
  tbl <- data.frame(
    locus_tag = "c5_mtase",
    product = "C5 DNA methyltransferase",
    `Signature accession_Pfam` = "PF00145",
    REBASEfinder_family_id = "Type II",
    REBASEfinder_enzyme_role = "M",
    translation = seq,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  signals <- DNMB:::.dnmb_rebasefinder_mtase_sequence_signals(seq, "PF00145 C5 DNA methyltransferase")
  hits <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl)

  expect_true(signals$verified)
  expect_identical(signals$architecture, "C5_PCQ-ENV")
  expect_setequal(hits$motif, c("N5C-PC", "N5C-ENV"))
  expect_true(all(hits$evidence_level == "supported"))
  expect_true(DNMB:::.dnmb_rebasefinder_motif_verified(tbl)[[1]])
})

test_that("MmeI-like RM fusions expose both MTase and nuclease functional evidence", {
  seq <- paste0(
    strrep("A", 40), "D", strrep("A", 10), "EAK", strrep("A", 45),
    "GAHYTS", strrep("A", 40), "FLDPACGSGNF", strrep("A", 110),
    "NPPF", strrep("A", 420)
  )
  tbl <- data.frame(
    locus_tag = "mme_rm",
    product = "Type IIG restriction enzyme/methyltransferase",
    `Signature accession_Pfam` = "PF20465, PF20466, PF20467, PF20473",
    `Signature description_Pfam` = paste(
      "MmeI helicase spacer domain", "MmeI target recognition domain",
      "MmeI C-terminal domain", "MmeI DNA-methyltransferase domain", sep = "; "
    ),
    REBASEfinder_family_id = "Type II",
    REBASEfinder_enzyme_role = "RM",
    translation = seq,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  hits <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl)
  rease_signals <- DNMB:::.dnmb_rebasefinder_rease_sequence_signals(
    seq,
    "class I SAM-dependent DNA methyltransferase"
  )
  expect_true(all(c("MmeI-X", "MmeI-I", "Amino-IV", "MmeI-PDExK", "MmeI-architecture") %in% hits$motif))
  expect_true(any(hits$expected_role == "M"))
  expect_true(any(hits$expected_role == "RM"))
  unsafe <- grepl("\n", hits$match, fixed = TRUE) |
    grepl("\r", hits$match, fixed = TRUE) |
    grepl("\t", hits$match, fixed = TRUE)
  expect_false(any(unsafe))
  expect_true(DNMB:::.dnmb_rebasefinder_motif_verified(tbl)[[1]])
  expect_true(rease_signals$mmei_fusion)
  expect_true(rease_signals$nuclease)
  expect_identical(rease_signals$signals, "MmeI-PDExK")
  expect_identical(rease_signals$architecture, "MmeI_fused_REase-MTase")
})

test_that("motif geometry detects sequence-distant residues that form a 3D active site", {
  aa3 <- c(
    A = "ALA", C = "CYS", D = "ASP", E = "GLU", F = "PHE", G = "GLY",
    H = "HIS", I = "ILE", K = "LYS", L = "LEU", M = "MET", N = "ASN",
    P = "PRO", Q = "GLN", R = "ARG", S = "SER", T = "THR", V = "VAL",
    W = "TRP", Y = "TYR"
  )
  write_model <- function(path, sequence, close = TRUE, plddt = 90) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
    xyz <- cbind(x = seq_along(chars) * 3.8, y = 0, z = 0)
    if (close) {
      xyz[80:83, ] <- cbind(x = c(39, 42, 45, 48), y = 4, z = 0)
    }
    lines <- vapply(seq_along(chars), function(i) {
      sprintf(
        "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C",
        i, aa3[[chars[[i]]]], i, xyz[i, 1], xyz[i, 2], xyz[i, 3], plddt
      )
    }, character(1))
    writeLines(c(lines, "END"), path)
  }

  sequence <- paste0(
    strrep("A", 9), "FAGAA", strrep("A", 65), "NPPY", strrep("A", 17)
  )
  motifs <- data.frame(
    locus_tag = "q1", family_id = "Type II", enzyme_role = "M",
    motif = c("SAM", "Amino-IV"), start_aa = c(10L, 80L), end_aa = c(14L, 83L),
    evidence_level = "supported", stringsAsFactors = FALSE
  )
  proteins <- data.frame(
    locus_tag = "q1", translation = sequence,
    REBASEfinder_structure_supported = TRUE,
    REBASEfinder_structure_candidate_consistent = TRUE,
    stringsAsFactors = FALSE
  )
  root <- tempfile("rebasefinder-geometry-")
  structure_dir <- file.path(root, "query_structures")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  model <- file.path(structure_dir, "q1.pdb")

  write_model(model, sequence, close = TRUE, plddt = 90)
  supported <- DNMB:::.dnmb_rebasefinder_verify_motif_geometry(
    motifs, proteins, structure_dirs = structure_dir
  )
  expect_identical(supported$pairs$combined_status, "3d_fold_supported")
  expect_lt(supported$pairs$min_ca_distance, 12)
  expect_identical(supported$summary$structural_adjacency_status, "3d_fold_supported")

  write_model(model, sequence, close = TRUE, plddt = 0.95)
  unit_scaled <- DNMB:::.dnmb_rebasefinder_verify_motif_geometry(
    motifs, proteins, structure_dirs = structure_dir
  )
  expect_identical(unit_scaled$pairs$combined_status, "3d_fold_supported")
  expect_gt(unit_scaled$pairs$motif_a_mean_plddt, 90)

  partial_sequence <- paste0(sequence, strrep("A", 60))
  proteins$translation <- partial_sequence
  write_model(model, substr(partial_sequence, 1, 83), close = TRUE, plddt = 90)
  local_support <- DNMB:::.dnmb_rebasefinder_verify_motif_geometry(
    motifs, proteins, structure_dirs = structure_dir
  )
  expect_identical(local_support$pairs$combined_status, "3d_fold_supported")
  expect_identical(local_support$summary$structure_model_scope, "partial_model")
  expect_match(local_support$summary$structural_adjacency_label, "local")

  proteins$translation <- sequence

  write_model(model, sequence, close = FALSE, plddt = 90)
  far <- DNMB:::.dnmb_rebasefinder_verify_motif_geometry(
    motifs, proteins, structure_dirs = structure_dir
  )
  expect_identical(far$pairs$geometry_status, "geometry_far")

  write_model(model, sequence, close = TRUE, plddt = 40)
  low_confidence <- DNMB:::.dnmb_rebasefinder_verify_motif_geometry(
    motifs, proteins, structure_dirs = structure_dir
  )
  expect_identical(low_confidence$pairs$geometry_status, "low_local_confidence")
})

test_that("MmeI geometry checks the MTase pocket without forcing nuclease-domain proximity", {
  motifs <- data.frame(
    locus_tag = "mme1", family_id = "Type II", enzyme_role = "RM",
    motif = c("MmeI-PDExK", "MmeI-I", "Amino-IV"),
    start_aa = c(40L, 150L, 280L), end_aa = c(54L, 160L, 283L),
    evidence_level = "supported", stringsAsFactors = FALSE
  )
  proteins <- data.frame(
    locus_tag = "mme1", translation = strrep("A", 700),
    stringsAsFactors = FALSE
  )
  result <- DNMB:::.dnmb_rebasefinder_verify_motif_geometry(
    motifs, proteins, structure_dirs = character()
  )

  expect_identical(result$pairs$pair_id, "mmei_mtase_I-IV")
  expect_false(any(grepl("PDExK", result$pairs$pair_id, fixed = TRUE)))
  expect_identical(result$pairs$combined_status, "no_structure")
})

test_that("ESM-2 motif contacts become pair outlines instead of heatmap columns", {
  root <- tempfile("rebasefinder-contacts-")
  module_dir <- file.path(root, "dnmb_module_rebasefinder")
  dir.create(module_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  contacts <- data.frame(
    locus_tag = c("strong", "weak"),
    family_id = "Type II",
    enzyme_role = "M",
    pair_id = "amino_mtase_I-IV",
    motif_a = "SAM",
    motif_b = "Amino-IV",
    max_contact_probability = c(0.63, 0.04),
    top3_mean_contact_probability = c(0.31, 0.02),
    separation_matched_percentile = c(99.9, 90),
    contact_status = c("contact_strong", "contact_weak"),
    model = "esm2_t6_8M_UR50D",
    stringsAsFactors = FALSE
  )
  write.table(
    contacts,
    file.path(module_dir, "DNMB_REBASEfinder_motif_contacts.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  result <- DNMB:::.dnmb_rebasefinder_read_motif_contacts(root)
  rownames(result$summary) <- result$summary$locus_tag

  expect_identical(result$summary["strong", "motif_contact_status"], "contact_strong")
  expect_identical(result$summary["weak", "motif_contact_status"], "contact_weak")
  expect_match(result$summary["strong", "motif_contact_label"], "CMAP\\+")
  expect_true(file.exists(result$xlsx))

  tbl <- data.frame(
    locus_tag = "strong",
    translation = paste0(strrep("A", 9), "FAGAA", strrep("A", 65), "NPPY", strrep("A", 17)),
    REBASEfinder_family_id = "Type II",
    REBASEfinder_enzyme_role = "M",
    REBASEfinder_hit_label = "M.Test",
    stringsAsFactors = FALSE
  )
  display <- DNMB:::.dnmb_rebasefinder_display_labels(tbl)
  plot <- DNMB:::.dnmb_plot_rebasefinder_motif_verification(
    tbl, display, contact_pairs = result$pairs, contact_summary = result$summary
  )
  boxes <- attr(plot, "dnmb_catalytic_pair_boxes")

  expect_equal(nrow(boxes), 1L)
  expect_identical(as.character(boxes$pair_state), "contact_supported")
  expect_identical(boxes$motif_a, "SAM")
  expect_identical(boxes$motif_b, "Amino-IV")
  expect_equal(boxes$xmin, 0.53)
  expect_equal(boxes$xmax, length(levels(plot$data$motif_label)) + 0.47)
  expect_false(any(plot$data$motif %in% c("ESM2-contact", "3D-adjacency")))
  expect_false(any(levels(plot$data$motif_label) %in% c("Contact\nmap", "3D\nmotifs")))
  expect_true(all(c("SAM", "MTase\nIV") %in% levels(plot$data$motif_label)))
})

test_that("catalytic pair outlines use aggregate 3D status and core motifs only", {
  mtase_sequence <- paste0(
    strrep("A", 9), "FAGAA", strrep("A", 65), "NPPY", strrep("A", 17)
  )
  tbl <- data.frame(
    locus_tag = c("verified", "unverified", "motor_only"),
    product = c(
      "DNA methyltransferase", "DNA methyltransferase",
      "Type I restriction enzyme R subunit"
    ),
    translation = c(mtase_sequence, mtase_sequence, strrep("A", 500)),
    REBASEfinder_family_id = c("Type II", "Type II", "Type I"),
    REBASEfinder_enzyme_role = c("M", "M", "R"),
    REBASEfinder_hit_label = c("M.Verified", "M.Unverified", "R.Motor"),
    REBASEfinder_homology_geometry_status = c(
      "homology_model_supported", "no_structure", "homology_model_motor_only"
    ),
    stringsAsFactors = FALSE
  )
  geometry_pairs <- data.frame(
    locus_tag = c("verified", "unverified", "motor_only"),
    motif_a = c("SAM", "SAM", "P-loop"),
    motif_b = c("Amino-IV", "Amino-IV", "HsdR-WB"),
    combined_status = c(
      "homology_model_supported", "no_structure", "homology_model_supported"
    ),
    stringsAsFactors = FALSE
  )
  geometry_summary <- data.frame(
    locus_tag = c("verified", "unverified", "motor_only"),
    structural_adjacency_status = c(
      "homology_model_supported", "no_structure", "homology_model_motor_only"
    ),
    stringsAsFactors = FALSE
  )

  display <- DNMB:::.dnmb_rebasefinder_display_labels(tbl)
  plot <- DNMB:::.dnmb_plot_rebasefinder_motif_verification(
    tbl,
    display,
    geometry_pairs = geometry_pairs,
    geometry_summary = geometry_summary
  )
  boxes <- attr(plot, "dnmb_catalytic_pair_boxes")
  rownames(boxes) <- boxes$locus_tag

  expect_identical(as.character(boxes["verified", "pair_state"]), "verified")
  expect_identical(as.character(boxes["unverified", "pair_state"]), "not_verified")
  expect_identical(as.character(boxes["motor_only", "pair_state"]), "not_verified")
  expect_identical(boxes["motor_only", "pair_group"], "hsdr_nuclease_motor")
  expect_identical(boxes["motor_only", "motif_a"], "HsdR-PD")
  expect_identical(boxes["motor_only", "motif_b"], "P-loop")
  expect_true(all(boxes$xmin == 0.53))
  expect_true(all(boxes$xmax == length(levels(plot$data$motif_label)) + 0.47))
  expect_false(any(boxes$pair_group %in% c("hsdr_motor_WA-WB", "hsdr_motor_WB-MIII")))
  expect_false(any(plot$data$motif %in% c("ESM2-contact", "3D-adjacency")))
})

test_that("Type IV evidence is profile-routed and never emitted as a generic Mrr sequence motif", {
  tbl <- data.frame(
    locus_tag = "gmrsd",
    product = "GmrSD modification-dependent restriction enzyme",
    `Signature accession_Pfam` = "PF03235, PF07510",
    REBASEfinder_family_id = "Type IV",
    REBASEfinder_enzyme_role = "R",
    translation = paste0(strrep("A", 100), "DABCDEFGHIJEAK", strrep("A", 500)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  hits <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl)
  raw <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl, role_relevant_only = FALSE)
  expect_false("Mrr" %in% names(DNMB:::.dnmb_rebasefinder_motif_definitions()))
  expect_identical(hits$motif, "TypeIV-profile")
  expect_identical(hits$typeiv_architecture, "GmrSD domains")
  expect_false(any(raw$motif == "Mrr"))
})

test_that("missing protein sequence yields unknown rather than motif-supported", {
  tbl <- data.frame(
    locus_tag = "missing",
    product = "restriction endonuclease",
    REBASEfinder_family_id = "Type II",
    REBASEfinder_enzyme_role = "R",
    translation = NA_character_,
    stringsAsFactors = FALSE
  )
  expect_true(is.na(DNMB:::.dnmb_rebasefinder_motif_verified(tbl)[[1]]))
})

test_that("REBASE partial status flags unusually short Type III Res proteins", {
  tbl <- data.frame(
    locus_tag = "res_short",
    translation = strrep("A", 200),
    REBASEfinder_family_id = "Type III",
    REBASEfinder_enzyme_role = "R",
    stringsAsFactors = FALSE
  )
  partial <- DNMB:::.dnmb_rebasefinder_sequence_partial_table(tbl)

  expect_equal(partial$partial_status[[1]], "partial_or_short")
  expect_match(partial$partial_reason[[1]], "below expected_min")
})

test_that("REBASE supplemental BLAST targets R/S context candidates missing BLAST", {
  hits <- data.frame(
    query = c("m1", "r1", "s1"),
    family_id = c("Type I", "Type I", "Type I"),
    enzyme_role = c("M", "R", "S"),
    evidence_mode = c("high_confidence", "operon_context", "operon_context"),
    blast_identity = c(1, NA, NA),
    stringsAsFactors = FALSE
  )
  mask <- DNMB:::.dnmb_rebasefinder_supplemental_query_mask(hits)

  expect_false(mask[[1]])
  expect_true(mask[[2]])
  expect_true(mask[[3]])
})

test_that("REBASE BLAST quality marks structural support on the BLAST point", {
  src <- c("R/plot_rebasefinder.R", "../../R/plot_rebasefinder.R")
  src <- src[file.exists(src)][1]
  skip_if_not(file.exists(src))
  text <- paste(readLines(src, warn = FALSE), collapse = "\n")
  expect_false(grepl("aes(x = 106", text, fixed = TRUE))
  expect_true(grepl("aes(x = .data$identity * 100", text, fixed = TRUE))
})

test_that("a successful integrated overview removes redundant native REBASE PDFs", {
  run_root <- tempfile("rebasefinder-integrated-pdf-")
  module_dir <- file.path(run_root, "dnmb_module_rebasefinder")
  visualization_dir <- file.path(run_root, "visualizations")
  dir.create(module_dir, recursive = TRUE)
  dir.create(visualization_dir, recursive = TRUE)
  on.exit(unlink(run_root, recursive = TRUE, force = TRUE), add = TRUE)

  root_native <- file.path(run_root, "RM_system_dotplot.pdf")
  module_native <- file.path(module_dir, "RM_system_dotplot.pdf")
  canonical <- file.path(visualization_dir, "REBASE_overview.pdf")
  writeLines("legacy", root_native)
  writeLines("legacy", module_native)
  writeLines("integrated", canonical)
  expected_removed <- normalizePath(
    c(root_native, module_native),
    winslash = "/",
    mustWork = TRUE
  )

  removed <- DNMB:::.dnmb_rebasefinder_cleanup_native_pdfs(run_root)

  expect_setequal(removed, expected_removed)
  expect_false(file.exists(root_native))
  expect_false(file.exists(module_native))
  expect_true(file.exists(canonical))
})

test_that("REBASE BLAST panel preserves context-only rows for y-axis alignment", {
  tbl <- data.frame(
    locus_tag = c("mod1", "res1"),
    contig = "ctg1",
    start = c(1, 101),
    end = c(90, 900),
    translation = c(strrep("A", 320), strrep("A", 780)),
    REBASEfinder_family_id = c("Type III", "Type III"),
    REBASEfinder_hit_label = c("M.Example", "typeIII_context:R:sequence_motif"),
    REBASEfinder_enzyme_role = c("M", "R"),
    REBASEfinder_typing_eligible = c(TRUE, FALSE),
    REBASEfinder_blast_identity = c(1, NA),
    REBASEfinder_blast_bitscore = c(500, NA),
    REBASEfinder_blast_length = c(320, NA),
    REBASEfinder_rec_seq = c("GCCAT", NA),
    REBASEfinder_operon_id = c("op1", "op1"),
    stringsAsFactors = FALSE
  )
  display_info <- DNMB:::.dnmb_rebasefinder_display_labels(tbl)
  p <- DNMB:::.dnmb_plot_rebasefinder_blast_quality(
    tbl,
    DNMB:::.dnmb_rebasefinder_role_palette(tbl),
    display_info,
    DNMB:::.dnmb_rebasefinder_palette(tbl$REBASEfinder_family_id)
  )
  y_limits <- p$scales$get_scales("y")$limits

  expect_equal(length(y_limits), 2L)
  expect_true(any(grepl("Type III R (operon)", y_limits, fixed = TRUE)))
  expect_true(any(grepl("M.Example", y_limits, fixed = TRUE)))
})

test_that("REBASE display order keeps operon-context partners adjacent", {
  tbl <- data.frame(
    locus_tag = c("modA", "resA", "specA", "modB"),
    contig = "ctg1",
    start = c(100, 1, 250, 1000),
    end = c(200, 90, 500, 1300),
    translation = c(strrep("A", 300), strrep("A", 760), strrep("A", 420), strrep("A", 310)),
    REBASEfinder_family_id = c("Type I", "Type I", "Type I", "Type II"),
    REBASEfinder_hit_label = c("M.A", "typeI_context:R:annotation", "typeI_context:S:annotation", "M.B"),
    REBASEfinder_enzyme_role = c("M", "R", "S", "M"),
    REBASEfinder_typing_eligible = c(TRUE, FALSE, FALSE, TRUE),
    REBASEfinder_blast_identity = c(1, 0.7, 0.8, 0.99),
    REBASEfinder_blast_bitscore = c(500, 300, 250, 490),
    REBASEfinder_blast_length = c(300, 760, 420, 310),
    REBASEfinder_operon_id = c("opA", NA, NA, "opB"),
    REBASEfinder_typei_context_partners = c(
      "resA:R:annotation:-1g:10bp | specA:S:annotation:+1g:10bp",
      NA, NA, NA
    ),
    stringsAsFactors = FALSE
  )
  display_info <- DNMB:::.dnmb_rebasefinder_display_labels(tbl)

  expect_equal(display_info$locus_tag[1:3], c("resA", "modA", "specA"))
  expect_equal(length(unique(display_info$operon[1:3])), 1L)
  expect_false(display_info$operon[[4]] %in% display_info$operon[1:3])
})

test_that("REBASE methylation annotations prioritize exact donor metadata and avoid ambiguous site inference", {
  cache_root <- tempfile("rebasefinder-plot-bairoch-")
  dir.create(cache_root)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)
  bairoch <- data.frame(
    enzyme_name = c("M.Donor", "R.Match", "M.Clean", "M.Amb1", "M.Amb2", "M.Unique1", "M.Unique2"),
    rec_seq = c("GATC", "CATG", "CCGG", "GGCC", "GGCC", "AAGCTT", "AAGCTT"),
    meth_type = c("m6A", "m4C", "m5C", "m5C", "m4C", "m6A", "m6A"),
    meth_pos = c("2", "3", "4", "1", "2", "3", "3"),
    meth_all = c("2(m6A)", "3(m4C)", "4(m5C)", "1(m5C)", "2(m4C)", "3(m6A)", "3(m6A)"),
    stringsAsFactors = FALSE
  )
  attr(bairoch, "dnmb_bairoch_schema") <- 2L
  saveRDS(bairoch, DNMB:::.dnmb_rebasefinder_bairoch_cache_path(cache_root))
  tbl <- data.frame(
    locus_tag = c("donor", "match", "clean", "ambiguous", "unique"),
    product = "hypothetical protein",
    translation = strrep("A", 250),
    REBASEfinder_hit_label = c("M.Clean_99", "M.Clean_99", "M.Clean_99", "typeII_context:R:test", "unknown"),
    REBASEfinder_recognition_donor = c("M.Donor", NA, NA, NA, NA),
    REBASEfinder_recognition_match = c("R.Match", "R.Match", NA, NA, NA),
    REBASEfinder_rec_seq = NA_character_,
    REBASEfinder_reference_rec_seq = c("GATC", "CATG", "CCGG", "GGCC", "AAGCTT"),
    stringsAsFactors = FALSE
  )

  out <- DNMB:::.dnmb_rebasefinder_methylation_annotations(tbl, cache_root = cache_root)
  rownames(out) <- out$locus_tag

  expect_identical(out["donor", "meth_type"], "N6A")
  expect_identical(out["donor", "meth_source"], "bairoch_recognition_donor")
  expect_identical(out["match", "meth_type"], "N4C")
  expect_identical(out["match", "meth_source"], "bairoch_recognition_match")
  expect_identical(out["clean", "meth_type"], "N5C")
  expect_identical(out["clean", "meth_source"], "bairoch_clean_hit")
  expect_true(is.na(out["ambiguous", "meth_type"]))
  expect_identical(out["unique", "meth_type"], "N6A")
  expect_identical(out["unique", "meth_pos"], "3")
  expect_identical(out["unique", "meth_source"], "bairoch_rec_seq_unique")
})

test_that("REBASE inventory excludes unassociated MTases from system counts but retains their detail scope", {
  tbl <- data.frame(
    locus_tag = c("classic_m", "orphan_m", "typeiv_r", "unresolved_r"),
    REBASEfinder_rm_association_class = c(
      "classic_RM_supported",
      "DNA_MTase_RM_association_unproven",
      "modification_dependent_restriction_supported",
      "RM_association_unproven"
    ),
    stringsAsFactors = FALSE
  )
  scope <- DNMB:::.dnmb_rebasefinder_inventory_scope(tbl)

  expect_identical(scope$system_mask, c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(scope$n_orphan_mtase, 1L)
  expect_identical(scope$n_unassociated, 2L)
  expect_identical(tbl$locus_tag[scope$system_mask], c("classic_m", "typeiv_r"))
  expect_identical(tbl$locus_tag[!scope$system_mask], c("orphan_m", "unresolved_r"))
})

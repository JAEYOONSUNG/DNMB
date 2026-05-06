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
  expect_equal(parsed["Qfok", "structure_family"], "Type IIS")
  expect_equal(parsed["Qfok", "structure_role"], "R")
  expect_equal(parsed["Qmtase", "structure_role"], "M")
  expect_equal(parsed["Qdna", "structure_status"], "structure_excluded")
  expect_false(parsed["Qdna", "structure_pass"])
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
})

test_that("Type III context rescues motif-inferred nearby missing Res partner", {
  res_seq <- paste0(
    paste(rep("A", 20), collapse = ""),
    "AGGGGGKS",
    paste(rep("A", 70), collapse = ""),
    "DEAH",
    paste(rep("A", 150), collapse = ""),
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
  expect_match(out["res1", "support"], "ResIII-WA\\+ResIII-WB\\+ResIII-PD")
})

test_that("Type III context rescues nearby long DEAD helicase Res partner", {
  helicase_seq <- paste0(
    paste(rep("A", 520), collapse = ""),
    "DEAQ",
    paste(rep("A", 120), collapse = "")
  )
  genes <- data.frame(
    locus_tag = c("res1", "mod1"),
    contig = "ctg1",
    start = c(100, 1200),
    end = c(1000, 2400),
    direction = "+",
    product = c("DEAD/DEAH box helicase", "site-specific DNA-methyltransferase"),
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

  expect_equal(out["mod1", "typeiii_context_status"], "complete_mod_res")
  expect_true("res1" %in% out$query)
  expect_equal(out["res1", "enzyme_role"], "R")
  expect_match(out["res1", "support"], "Helicase-DEAD")
})

test_that("Type I context rescues HsdR and HsdS around a Type I methylase", {
  genes <- data.frame(
    locus_tag = c("hsdR1", "hsdM1", "hsdS1"),
    contig = "ctg1",
    start = c(100, 1200, 2400),
    end = c(1000, 2200, 3400),
    direction = "+",
    gene = c(NA, NA, NA),
    product = c("DEAD/DEAH box helicase family protein",
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
})

test_that("Type IV candidates are rescued from Mrr-like restriction annotations", {
  seq <- paste0("M", "D", paste(rep("A", 8), collapse = ""), "EAK", paste(rep("A", 30), collapse = ""), "DK")
  genes <- data.frame(
    locus_tag = "typeiv1",
    contig = "ctg1",
    start = 100,
    end = 700,
    direction = "+",
    product = "restriction endonuclease-like protein",
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

test_that("REBASE overview motif scan captures Type I HsdR helicase motif III", {
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

  expect_match(motifs[["HsdR-MIII"]][[1]], "^present")
  expect_true(any(hits$motif == "HsdR-MIII"))
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

  manual_like <- paste0(strrep("A", 26), "FAGSA", strrep("A", 85), "DPPY", strrep("A", 20))
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

test_that("Type I specificity subunits use HsdS signature instead of MTase SAM in default motif list", {
  tbl <- data.frame(
    locus_tag = "hsds1",
    translation = paste0(strrep("A", 80), "FTGTA", strrep("A", 360)),
    REBASEfinder_family_id = "Type I",
    REBASEfinder_enzyme_role = "S",
    stringsAsFactors = FALSE
  )
  motifs <- DNMB:::.dnmb_rebasefinder_scan_motifs(tbl)
  hits <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl)
  raw <- DNMB:::.dnmb_rebasefinder_motif_hits_table(tbl, role_relevant_only = FALSE)

  expect_match(motifs[["HsdS-FxGxA"]][[1]], "^present")
  expect_true(any(hits$motif == "HsdS-FxGxA"))
  expect_false(any(hits$motif == "SAM-like"))
  expect_true(any(raw$motif == "SAM" & !raw$role_relevant))
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

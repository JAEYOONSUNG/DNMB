test_that("Type I-S motif parser preserves N inside half-sites", {
  parsed <- DNMB:::.dnmb_type1s_parse_motif("GNAANNNNNNNGTGG")
  expect_equal(parsed$left, "GNAA")
  expect_equal(parsed$spacer, 7L)
  expect_equal(parsed$right, "GTGG")

  parsed_right_n <- DNMB:::.dnmb_type1s_parse_motif("CAAGNNNNNNTCNNC")
  expect_equal(parsed_right_n$left, "CAAG")
  expect_equal(parsed_right_n$spacer, 6L)
  expect_equal(parsed_right_n$right, "TCNNC")

  compact <- DNMB:::.dnmb_type1s_parse_motif("CAA(6)CTC")
  expect_equal(compact$recognition, "CAANNNNNNCTC")
})

test_that("Type I-S Gold reference uses the official HTTPS endpoint", {
  expect_identical(
    DNMB:::.dnmb_type1s_gold_url(),
    "https://rebase.neb.com/rebase/Type_I_S_subunit_Gold_Standards_Protein.txt"
  )
})

test_that("IUPAC reverse complement is orientation aware", {
  expect_equal(
    unname(DNMB:::.dnmb_type1s_reverse_complement("RTAG")),
    "CTAY"
  )
  expect_equal(
    unname(DNMB:::.dnmb_type1s_reverse_complement("TCNNC")),
    "GNNGA"
  )
})

test_that("candidate selection is restricted to Type I HsdS", {
  genes <- data.frame(
    locus_tag = c(
      "hsds_1", "hsds_2", "type1_product", "refseq_subunit_s",
      "type2_product", "type4_product", "other"
    ),
    product = c(
      "HsdS specificity protein",
      "HsdS specificity protein",
      "Type I restriction enzyme specificity subunit",
      "restriction endonuclease subunit S",
      "Type II restriction enzyme specificity subunit",
      "Type IV restriction enzyme specificity subunit",
      "DNA-binding protein"
    ),
    stringsAsFactors = FALSE
  )
  rm <- data.frame(
    locus_tag = c("from_type_i", "from_type_ii", "from_type_iii"),
    rm_type = c("Type I", "Type II", "Type III"),
    subunit = c("S", "S", "S"),
    stringsAsFactors = FALSE
  )
  ids <- DNMB:::.dnmb_type1s_candidate_ids(genes, rm)
  expect_setequal(ids, c(
    "hsds_1", "hsds_2", "type1_product", "refseq_subunit_s", "from_type_i"
  ))
  expect_false(any(c("type2_product", "type4_product", "from_type_ii", "from_type_iii") %in% ids))

  expect_identical(
    DNMB:::.dnmb_type1s_candidate_rows(
      c("foo_1", "foo_2"),
      "foo_1"
    ),
    c(TRUE, FALSE)
  )
  expect_identical(
    DNMB:::.dnmb_type1s_candidate_rows(
      c("foo_1", "foo_2"),
      "gnl|foo_1"
    ),
    c(FALSE, FALSE)
  )
  expect_true(DNMB:::.dnmb_type1s_candidate_rows("foo_1", "gnl|foo_1"))
})

test_that("empty Type I-S candidate sets return a typed zero-row table", {
  empty <- DNMB:::.dnmb_type1s_empty_prediction(character())
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0L)
  expect_true(all(c(
    "locus_tag",
    "type1s_predicted_recognition",
    "type1s_prediction_eligible"
  ) %in% names(empty)))

  public_empty <- DNMB::dnmb_predict_type1s_recognition(
    data.frame(locus_tag = character(), translation = character()),
    download = FALSE,
    verbose = FALSE
  )
  expect_equal(nrow(public_empty), 0L)
  expect_identical(names(public_empty), names(empty))
  expect_identical(
    DNMB:::.dnmb_type1s_clean_id(c(" hsds_1 ", "hsds_2")),
    c("hsds_1", "hsds_2")
  )
})

test_that("cross-position TRD hits are reverse-complemented", {
  reference <- data.frame(
    reference_id = "R0001",
    source_enzyme = "S.testI",
    left_half = "CGA",
    right_half = "RTAG",
    spacer_length = 7L,
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    query_id = "Q0001_T1",
    subject_id = "R0001_T2",
    identity = 80,
    alignment_length = 100,
    query_length = 100,
    subject_length = 100,
    query_start = 1,
    query_end = 100,
    subject_start = 1,
    subject_end = 100,
    evalue = 1e-20,
    bitscore = 200,
    query_coverage = 100,
    subject_coverage = 100,
    stringsAsFactors = FALSE
  )
  selected <- DNMB:::.dnmb_type1s_select_trd_hit(hits, "Q0001_T1", "T1", reference)
  expect_equal(unname(selected$half_site), "CTAY")
  expect_equal(selected$source_position, "T2")

  hits$query_id <- "Q0001_T2"
  hits$subject_id <- "R0001_T1"
  selected_right <- DNMB:::.dnmb_type1s_select_trd_hit(
    hits,
    "Q0001_T2",
    "T2",
    reference
  )
  expect_equal(unname(selected_right$half_site), "TCG")
  expect_equal(selected_right$source_position, "T1")
})

test_that("validated same-versus-cross TRD priors are stable", {
  reference <- data.frame(
    reference_id = c("R0001", "R0002"),
    source_enzyme = c("S.sameI", "S.crossI"),
    left_half = c("AAA", "CCC"),
    right_half = c("GGG", "TTT"),
    spacer_length = c(6L, 7L),
    stringsAsFactors = FALSE
  )
  hit <- function(subject_id, bitscore) {
    data.frame(
      query_id = "Q0001_T2",
      subject_id = subject_id,
      identity = 80,
      query_coverage = 100,
      subject_coverage = 100,
      bitscore = bitscore,
      stringsAsFactors = FALSE
    )
  }
  selected_cross <- DNMB:::.dnmb_type1s_select_trd_hit(
    rbind(hit("R0001_T2", 100), hit("R0002_T1", 98)),
    "Q0001_T2",
    "T2",
    reference
  )
  expect_equal(selected_cross$source_position, "T1")

  selected_same <- DNMB:::.dnmb_type1s_select_trd_hit(
    rbind(hit("R0001_T2", 100), hit("R0002_T1", 97)),
    "Q0001_T2",
    "T2",
    reference
  )
  expect_equal(selected_same$source_position, "T2")
})

test_that("Type I references are materialized as separate TRD1 and TRD2 databases", {
  sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 20), collapse = "")
  windows <- DNMB:::.dnmb_type1s_extract_windows(sequence)
  reference <- list(table = data.frame(
    reference_id = "R0001",
    source_enzyme = "S.syntheticI",
    sequence = sequence,
    recognition = "CAANNNNNNCTC",
    left_half = "CAA",
    right_half = "CTC",
    spacer_length = 6L,
    trd1_sequence = windows$trd1,
    trd2_sequence = windows$trd2,
    stringsAsFactors = FALSE
  ))
  reference_dir <- tempfile("type1s-separate-db-")
  dir.create(reference_dir)
  on.exit(unlink(reference_dir, recursive = TRUE, force = TRUE), add = TRUE)

  files <- DNMB:::.dnmb_type1s_prepare_reference_files(reference, reference_dir)
  expect_true(all(file.exists(unlist(files))))
  expect_equal(readLines(files$trd1_fasta)[[1]], ">R0001_T1")
  expect_equal(readLines(files$trd2_fasta)[[1]], ">R0001_T2")
  expect_equal(readLines(files$scaffold_fasta)[[1]], ">R0001")

  trd1_metadata <- read.delim(files$trd1_metadata, check.names = FALSE)
  trd2_metadata <- read.delim(files$trd2_metadata, check.names = FALSE)
  expect_equal(trd1_metadata$source_position, "T1")
  expect_equal(trd2_metadata$source_position, "T2")
  expect_equal(trd1_metadata$position_half_site, "CAA")
  expect_equal(trd2_metadata$position_half_site, "CTC")
  expect_equal(trd1_metadata$canonical_half_site, "CAA")
  expect_equal(trd2_metadata$canonical_half_site, "GAG")
  expect_equal(trd2_metadata$cross_position_half_site, "GAG")
  expect_false(trd1_metadata$spacer_ambiguous)
  expect_false(trd2_metadata$spacer_ambiguous)
  expect_equal(as.character(trd1_metadata$spacer_set), "6")
  expect_equal(trd1_metadata$window_start, windows$trd1_start)
  expect_equal(trd2_metadata$window_end, windows$trd2_end)
  expect_equal(trd1_metadata$boundary_method, "gold_cv_specificity_core_window_v1")
  expect_equal(trd1_metadata$evidence, "REBASE_Type_I_S_Gold_Standard")
  scaffold_metadata <- read.delim(files$scaffold_metadata, check.names = FALSE)
  expect_equal(scaffold_metadata$scaffold_method, "conserved_scaffold_fixed_windows_v1")
  expect_equal(
    scaffold_metadata$spacer_model_version,
    "scaffold-ensemble-tael-v1.1.0"
  )
  expect_equal(scaffold_metadata$scaffold_length, nchar(windows$spacer_scaffold))

  writeLines(character(), files$trd1_fasta)
  writeLines("broken", files$trd2_metadata)
  repaired <- DNMB:::.dnmb_type1s_prepare_reference_files(reference, reference_dir)
  expect_equal(readLines(repaired$trd1_fasta)[[1]], ">R0001_T1")
  repaired_metadata <- read.delim(repaired$trd2_metadata, check.names = FALSE)
  expect_equal(nrow(repaired_metadata), 1L)
  expect_true("spacer_ambiguous" %in% names(repaired_metadata))
})

test_that("conserved HsdS scaffold captures helix and four-residue repeat features", {
  sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 20), collapse = "")
  substr(sequence, 121, 132) <- "TAELTAELTAEL"
  windows <- DNMB:::.dnmb_type1s_extract_windows(sequence)
  features <- DNMB:::.dnmb_type1s_scaffold_features(sequence)

  expect_equal(
    nchar(windows$trd1) + nchar(windows$trd2) + nchar(windows$spacer_scaffold),
    nchar(sequence)
  )
  expect_identical(windows$spacer_scaffold, features$sequence)
  expect_equal(features$central_start, 121L)
  expect_equal(features$central_end, 264L)
  expect_equal(features$tael_like_repeat_count, 3L)
  expect_gte(features$tetrapeptide_repeat_count, 3L)
  expect_true(features$helix_friendly_fraction > 0)
  expect_true(features$helix_breaker_fraction >= 0)
  expect_null(DNMB:::.dnmb_type1s_scaffold_features("ACDEFG"))
})

test_that("cached Type I reference paths are re-anchored after cache relocation", {
  cache_dir <- tempfile("type1s-relocated-cache-")
  dir.create(cache_dir, recursive = TRUE)
  on.exit(unlink(cache_dir, recursive = TRUE, force = TRUE), add = TRUE)
  gold_path <- file.path(cache_dir, DNMB:::.dnmb_type1s_gold_filename())
  writeLines("relocated gold placeholder", gold_path)
  source_md5 <- unname(tools::md5sum(gold_path))
  version <- paste(
    DNMB:::.dnmb_type1s_parser_version(),
    DNMB:::.dnmb_type1s_database_version(),
    DNMB:::.dnmb_type1s_spacer_model_version(),
    substr(source_md5, 1L, 12L),
    sep = "-"
  )
  reference_dir <- file.path(cache_dir, "type1s", version)
  dir.create(reference_dir, recursive = TRUE)
  table <- data.frame(
    reference_id = "R0001",
    source_enzyme = "S.syntheticI",
    sequence = paste(rep("A", 400), collapse = ""),
    recognition = "CAANNNNNNCTC",
    left_half = "CAA",
    right_half = "CTC",
    spacer_length = 6L,
    trd1_sequence = paste(rep("A", 75), collapse = ""),
    trd2_sequence = paste(rep("A", 64), collapse = ""),
    trd1_spacer_set = "6",
    trd1_spacer_ambiguous = FALSE,
    trd2_spacer_set = "6",
    trd2_spacer_ambiguous = FALSE,
    stringsAsFactors = FALSE
  )
  saveRDS(
    list(
      table = table,
      source_path = "/old/cache/gold.txt",
      reference_dir = "/old/cache/reference",
      version = version
    ),
    file.path(reference_dir, "reference.rds")
  )

  relocated <- testthat::with_mocked_bindings(
    DNMB:::.dnmb_type1s_reference(
      cache_dir,
      download = FALSE,
      verbose = FALSE
    ),
    .dnmb_type1s_valid_gold_file = function(...) TRUE,
    .package = "DNMB"
  )
  expect_equal(relocated$source_path, normalizePath(gold_path, winslash = "/"))
  expect_equal(relocated$reference_dir, reference_dir)
  expect_equal(relocated$database_version, DNMB:::.dnmb_type1s_database_version())
})

test_that("ambiguous TRD spacer labels fall back to the whole-HsdS neighbor", {
  hit <- list(
    spacer = 6L,
    spacer_ambiguous = FALSE,
    identity = 95,
    coverage = 90,
    subject_coverage = 90,
    source = "S.syntheticI"
  )
  full <- list(
    spacer = 7L,
    identity = 85,
    coverage = 90,
    subject_coverage = 90,
    source = "S.fullI",
    method = "whole_hsds_1nn"
  )

  consensus <- DNMB:::.dnmb_type1s_choose_spacer(hit, hit, full)
  expect_equal(consensus$spacer, 6L)
  expect_equal(consensus$method, "trd_source_consensus")

  ambiguous <- hit
  ambiguous$spacer_ambiguous <- TRUE
  fallback <- DNMB:::.dnmb_type1s_choose_spacer(ambiguous, hit, full)
  expect_identical(fallback, full)

  remote <- hit
  remote$identity <- 70
  remote_fallback <- DNMB:::.dnmb_type1s_choose_spacer(remote, remote, full)
  expect_identical(remote_fallback, full)

  last_resort <- DNMB:::.dnmb_type1s_choose_spacer(remote, remote, NULL)
  expect_equal(last_resort$spacer, 6L)
  expect_equal(last_resort$method, "trd_source_consensus")

  high_full <- full
  high_full$identity <- 95
  weighted <- full
  weighted$spacer <- 8L
  weighted$method <- "whole_hsds_weighted_knn5"
  disagreeing <- hit
  disagreeing$spacer <- 7L
  expect_identical(
    DNMB:::.dnmb_type1s_choose_spacer(hit, disagreeing, high_full, weighted),
    high_full
  )
  expect_identical(
    DNMB:::.dnmb_type1s_choose_spacer(remote, remote, full, weighted),
    weighted
  )
})

test_that("whole-HsdS spacer kNN uses reciprocal-coverage weighted voting", {
  reference <- data.frame(
    reference_id = sprintf("R%04d", 1:5),
    source_enzyme = paste0("S.vote", 1:5, "I"),
    spacer_length = c(7L, 6L, 6L, 6L, 7L),
    stringsAsFactors = FALSE
  )
  hits <- data.frame(
    subject_id = reference$reference_id,
    identity = c(100, 80, 80, 80, 60),
    query_coverage = 100,
    subject_coverage = 100,
    bitscore = 5:1,
    stringsAsFactors = FALSE
  )
  vote <- DNMB:::.dnmb_type1s_select_spacer_knn(hits, reference, k = 5L)
  expect_equal(vote$spacer, 6L)
  expect_equal(vote$neighbor_count, 5L)
  expect_equal(vote$method, "whole_hsds_weighted_knn5")
  expect_equal(vote$identity, 80)
  expect_equal(vote$source, "S.vote2I")
  expect_equal(vote$supporter_count, 3L)
  expect_equal(vote$voter_ids, paste(reference$reference_id, collapse = ","))
  expect_equal(vote$supporter_ids, paste(reference$reference_id[2:4], collapse = ","))
  expect_gt(vote$vote_support, 0.5)
  expect_gt(vote$vote_margin, 0)

  tie_reference <- reference[1:4, ]
  tie_reference$spacer_length <- c(7L, 6L, 6L, 7L)
  tied <- hits[1:4, ]
  tied$identity <- 100
  tie_vote <- DNMB:::.dnmb_type1s_select_spacer_knn(
    tied,
    tie_reference,
    k = 4L
  )
  expect_equal(tie_vote$spacer, 7L)
  expect_equal(tie_vote$method, "whole_hsds_weighted_knn4")
  expect_equal(tie_vote$vote_margin, 0)

  nearest <- DNMB:::.dnmb_type1s_select_spacer_1nn(hits, reference)
  expect_equal(nearest$spacer, 7L)
  expect_equal(nearest$source, "S.vote1I")
  expect_equal(nearest$method, "whole_hsds_1nn")
})

test_that("conserved-scaffold evidence reranks only the remote spacer vote", {
  reference <- data.frame(
    reference_id = sprintf("R%04d", 1:25),
    source_enzyme = paste0("S.vote", 1:25, "I"),
    spacer_length = c(6L, 6L, 7L, 7L, 7L, rep(7L, 20)),
    stringsAsFactors = FALSE
  )
  make_hits <- function(subject_id, identity) {
    data.frame(
      query_id = "Q0001",
      subject_id = subject_id,
      identity = identity,
      query_coverage = 100,
      subject_coverage = 100,
      bitscore = identity,
      stringsAsFactors = FALSE
    )
  }
  full_hits <- make_hits(
    reference$reference_id[1:5],
    c(100, 100, 65, 65, 65)
  )
  scaffold_hits <- make_hits(reference$reference_id[6:25], rep(90, 20))

  baseline <- DNMB:::.dnmb_type1s_select_spacer_knn(full_hits, reference, k = 5L)
  ensemble <- DNMB:::.dnmb_type1s_select_spacer_ensemble(
    full_hits,
    scaffold_hits,
    reference,
    full_k = 5L,
    scaffold_k = 20L,
    full_weight = 0.70
  )
  expect_equal(baseline$spacer, 6L)
  expect_equal(ensemble$spacer, 7L)
  expect_equal(ensemble$full_vote, 6L)
  expect_equal(ensemble$scaffold_vote, 7L)
  expect_false(ensemble$model_agreement)
  expect_false(ensemble$can_be_high)
  expect_equal(ensemble$method, "whole_hsds_scaffold_ensemble")
  expect_match(ensemble$voter_ids, "F:R0001")
  expect_match(ensemble$voter_ids, "S:R0006")

  fallback <- DNMB:::.dnmb_type1s_select_spacer_ensemble(
    full_hits,
    data.frame(),
    reference
  )
  expect_equal(fallback$spacer, baseline$spacer)
  expect_equal(fallback$method, "whole_hsds_weighted_knn5")
})

test_that("TAEL molecular ruler is restricted to calibrated canonical pairs", {
  close_reference <- list(
    spacer = 7L,
    identity = 96.8,
    coverage = 100,
    subject_coverage = 100,
    reference_id = "R0262",
    tael_like_repeat_count = 3L,
    method = "whole_hsds_1nn",
    can_be_high = TRUE
  )
  close <- DNMB:::.dnmb_type1s_apply_tael_ruler(
    close_reference,
    close_reference,
    query_repeat_count = 2L,
    mode = "close"
  )
  expect_equal(close$spacer, 6L)
  expect_equal(close$method, "whole_hsds_1nn_tael_ruler")
  expect_true(close$can_be_high)
  expect_true(close$tael_ruler_applied)
  expect_equal(close$tael_base_spacer, 7L)
  expect_equal(close$tael_reference_id, "R0262")
  expect_equal(close$tael_reference_identity, 96.8)
  expect_equal(close$tael_reference_reciprocal_coverage, 100)
  expect_equal(close$tael_base_method, "whole_hsds_1nn")

  ambiguous_left <- list(
    spacer = 7L, spacer_ambiguous = TRUE, identity = 96,
    coverage = 90, subject_coverage = 90, reference_id = "R0182",
    source = "S.leftI"
  )
  close_right <- list(
    spacer = 7L, spacer_ambiguous = FALSE, identity = 97,
    coverage = 90, subject_coverage = 90, reference_id = "R0262",
    source = "S.rightI"
  )
  expect_equal(
    DNMB:::.dnmb_type1s_choose_spacer(ambiguous_left, close_right, close)$method,
    "whole_hsds_1nn_tael_ruler"
  )
  unambiguous_left <- ambiguous_left
  unambiguous_left$spacer_ambiguous <- FALSE
  expect_equal(
    DNMB:::.dnmb_type1s_choose_spacer(unambiguous_left, close_right, close)$method,
    "trd_source_consensus"
  )

  remote_reference <- close_reference
  remote_reference$spacer <- 6L
  remote_reference$identity <- 76.8
  remote_reference$coverage <- 99
  remote_reference$subject_coverage <- 98
  remote_reference$reference_id <- "R0625"
  remote_reference$tael_like_repeat_count <- 2L
  remote <- list(
    spacer = 6L,
    identity = 80,
    coverage = 95,
    subject_coverage = 95,
    method = "whole_hsds_scaffold_ensemble",
    vote_margin = 0.39,
    can_be_high = FALSE
  )
  corrected_remote <- DNMB:::.dnmb_type1s_apply_tael_ruler(
    remote,
    remote_reference,
    query_repeat_count = 3L,
    mode = "remote"
  )
  expect_equal(corrected_remote$spacer, 7L)
  expect_equal(corrected_remote$method, "whole_hsds_scaffold_tael_ruler")
  expect_false(corrected_remote$can_be_high)
  expect_true(corrected_remote$tael_ruler_applied)
  expect_equal(corrected_remote$tael_base_method, "whole_hsds_scaffold_ensemble")

  already_correct <- remote
  already_correct$spacer <- 7L
  expect_identical(
    DNMB:::.dnmb_type1s_apply_tael_ruler(
      already_correct, remote_reference, 3L, mode = "remote"
    ),
    already_correct
  )

  high_margin <- remote
  high_margin$vote_margin <- 0.41
  expect_identical(
    DNMB:::.dnmb_type1s_apply_tael_ruler(
      high_margin, remote_reference, 3L, mode = "remote"
    ),
    high_margin
  )
  low_identity <- remote_reference
  low_identity$identity <- 69.9
  expect_identical(
    DNMB:::.dnmb_type1s_apply_tael_ruler(
      remote, low_identity, 3L, mode = "remote"
    ),
    remote
  )
  noncanonical_reference <- remote_reference
  noncanonical_reference$spacer <- 7L
  expect_identical(
    DNMB:::.dnmb_type1s_apply_tael_ruler(
      remote, noncanonical_reference, 3L, mode = "remote"
    ),
    remote
  )
  expect_identical(
    DNMB:::.dnmb_type1s_apply_tael_ruler(
      remote, remote_reference, 2L, mode = "remote"
    ),
    remote
  )
})

test_that("public Type I TRD database builder reports both position databases", {
  sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 20), collapse = "")
  windows <- DNMB:::.dnmb_type1s_extract_windows(sequence)
  reference_dir <- tempfile("type1s-public-builder-")
  dir.create(reference_dir)
  on.exit(unlink(reference_dir, recursive = TRUE, force = TRUE), add = TRUE)
  reference <- list(
    table = data.frame(
      reference_id = "R0001",
      source_enzyme = "S.syntheticI",
      sequence = sequence,
      recognition = "CAANNNNNNCTC",
      left_half = "CAA",
      right_half = "CTC",
      spacer_length = 6L,
      trd1_sequence = windows$trd1,
      trd2_sequence = windows$trd2,
      stringsAsFactors = FALSE
    ),
    version = "test-separate-trd",
    database_version = "separate-trd-v1.2.0",
    source_path = "gold.txt",
    source_md5 = "test-md5",
    reference_dir = reference_dir
  )

  built <- testthat::with_mocked_bindings(
    DNMB::dnmb_build_type1s_trd_databases(
      cache_dir = reference_dir,
      download = FALSE,
      backend = "both",
      verbose = FALSE
    ),
    .dnmb_type1s_reference = function(...) reference,
    dnmb_detect_binary = function(...) list(found = TRUE),
    .dnmb_type1s_prepare_diamond_db = function(...) list(
      trd1 = "diamond-trd1", trd2 = "diamond-trd2", full = "diamond-full"
    ),
    .dnmb_type1s_prepare_blast_db = function(...) list(
      trd1 = "blast-trd1", trd2 = "blast-trd2", full = "blast-full"
    ),
    .package = "DNMB"
  )

  expect_s3_class(built, "dnmb_type1s_databases")
  expect_equal(built$reference_records, 1L)
  expect_equal(built$statistics$trd1$conflicting_half_site_sequences, 0L)
  expect_equal(built$statistics$trd2$conflicting_half_site_sequences, 0L)
  expect_equal(built$statistics$trd1$conflicting_spacer_sequences, 0L)
  expect_equal(built$statistics$trd2$conflicting_spacer_sequences, 0L)
  expect_equal(built$databases$diamond$trd1, "diamond-trd1")
  expect_equal(built$databases$blastp$trd2, "blast-trd2")
  expect_true(file.exists(built$manifest_path))
  expect_true(all(c("artifact", "backend", "path", "size", "md5") %in% names(built$manifest)))
  expect_true(all(built$manifest$exists[built$manifest$backend == "source"]))

  auto_fallback <- testthat::with_mocked_bindings(
    DNMB::dnmb_build_type1s_trd_databases(
      cache_dir = reference_dir,
      download = FALSE,
      backend = "auto",
      verbose = FALSE
    ),
    .dnmb_type1s_reference = function(...) reference,
    dnmb_detect_binary = function(...) list(found = TRUE),
    .dnmb_type1s_prepare_diamond_db = function(...) stop("forced DIAMOND build failure"),
    .dnmb_type1s_prepare_blast_db = function(...) list(
      trd1 = "blast-trd1", trd2 = "blast-trd2", full = "blast-full"
    ),
    .package = "DNMB"
  )
  expect_null(auto_fallback$databases$diamond)
  expect_equal(auto_fallback$databases$blastp$trd1, "blast-trd1")
})

test_that("Type I-S manifests are content-stable and track artifact changes", {
  reference_dir <- tempfile("type1s-stable-manifest-")
  dir.create(reference_dir, recursive = TRUE)
  on.exit(unlink(reference_dir, recursive = TRUE, force = TRUE), add = TRUE)

  source_file <- file.path(reference_dir, "type1s_trd1.faa")
  writeLines(c(">R0001_T1", "ACDEFGHIK"), source_file)
  reference <- list(
    version = "type1s-test-reference",
    database_version = "type1s-test-database",
    source_url = DNMB:::.dnmb_type1s_gold_url(),
    source_md5 = "gold-source-md5",
    built_at = "2026-07-15T00:00:00Z",
    reference_dir = reference_dir
  )

  first <- DNMB:::.dnmb_type1s_write_manifest(
    reference,
    files = list(trd1_fasta = source_file),
    databases = list()
  )
  first_content <- readBin(first$path, what = "raw", n = file.info(first$path)$size)
  first_signature <- DNMB:::.dnmb_file_signature(first$path)
  stable_mtime <- as.POSIXct("2020-01-01T00:00:00Z", tz = "UTC")
  Sys.setFileTime(first$path, stable_mtime)
  expect_false("built_at" %in% names(first$table))
  expect_identical(
    unique(first$table$source_url),
    DNMB:::.dnmb_type1s_gold_url()
  )

  rewritten_reference <- reference
  rewritten_reference$built_at <- "2030-01-01T00:00:00Z"
  second <- DNMB:::.dnmb_type1s_write_manifest(
    rewritten_reference,
    files = list(trd1_fasta = source_file),
    databases = list()
  )
  second_content <- readBin(second$path, what = "raw", n = file.info(second$path)$size)
  second_signature <- DNMB:::.dnmb_file_signature(second$path)
  expect_identical(second_content, first_content)
  expect_identical(second_signature, first_signature)
  expect_equal(as.numeric(file.info(second$path)$mtime), as.numeric(stable_mtime))

  writeLines(c(">R0001_T1", "ACDEFGHIKL"), source_file)
  changed <- DNMB:::.dnmb_type1s_write_manifest(
    rewritten_reference,
    files = list(trd1_fasta = source_file),
    databases = list()
  )
  changed_content <- readBin(changed$path, what = "raw", n = file.info(changed$path)$size)
  changed_signature <- DNMB:::.dnmb_file_signature(changed$path)
  expect_false(identical(changed_content, first_content))
  expect_false(identical(changed_signature, first_signature))
  expect_identical(
    changed$table$md5[[1]],
    unname(tools::md5sum(source_file))
  )
})

test_that("exact Gold protein match is emitted as known", {
  sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 20), collapse = "")
  windows <- DNMB:::.dnmb_type1s_extract_windows(sequence)
  reference_table <- data.frame(
    reference_id = "R0001",
    source_enzyme = "S.syntheticI",
    sequence = sequence,
    recognition = "CAANNNNNNCTC",
    left_half = "CAA",
    right_half = "CTC",
    spacer_length = 6L,
    trd1_sequence = windows$trd1,
    trd2_sequence = windows$trd2,
    stringsAsFactors = FALSE
  )
  reference <- list(
    table = reference_table,
    version = "test-reference",
    reference_dir = tempdir()
  )
  prediction <- DNMB:::.dnmb_type1s_predict_core(
    data.frame(locus_tag = "query", translation = sequence),
    reference,
    cache_dir = tempdir(),
    verbose = FALSE
  )
  expect_equal(prediction$type1s_prediction_status, "known_gold_match")
  expect_equal(prediction$type1s_predicted_recognition, "CAANNNNNNCTC")
  expect_equal(prediction$type1s_overall_confidence, "known")
  expect_true(prediction$type1s_prediction_eligible)
  expect_equal(prediction$type1s_prediction_version, "type1s-predict-v1.4.0")
  expect_false(prediction$type1s_trd1_spacer_ambiguous)
  expect_false(prediction$type1s_trd2_spacer_ambiguous)
})

test_that("DIAMOND search failure falls back to BLASTP", {
  sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 20), collapse = "")
  query <- sequence
  substr(query, 210, 210) <- "A"
  windows <- DNMB:::.dnmb_type1s_extract_windows(sequence)
  reference <- list(
    table = data.frame(
      reference_id = "R0001",
      source_enzyme = "S.syntheticI",
      sequence = sequence,
      recognition = "CAANNNNNNCTC",
      left_half = "CAA",
      right_half = "CTC",
      spacer_length = 6L,
      trd1_sequence = windows$trd1,
      trd2_sequence = windows$trd2,
      stringsAsFactors = FALSE
    ),
    version = "test-reference",
    reference_dir = tempdir()
  )
  fake_blast <- function(query_fasta, db_prefix, output_path, cpu = 1L) {
    if (grepl("query_trd", query_fasta, fixed = TRUE)) {
      is_trd1 <- grepl("trd1", db_prefix, fixed = TRUE)
      return(data.frame(
        query_id = if (is_trd1) "Q0001_T1" else "Q0001_T2",
        subject_id = if (is_trd1) "R0001_T1" else "R0001_T2",
        identity = 95,
        query_coverage = 100,
        subject_coverage = 100,
        bitscore = 200,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      query_id = "Q0001",
      subject_id = "R0001",
      identity = 95,
      query_coverage = 100,
      subject_coverage = 100,
      bitscore = 400,
      stringsAsFactors = FALSE
    )
  }

  prediction <- testthat::with_mocked_bindings(
    DNMB:::.dnmb_type1s_predict_core(
      data.frame(locus_tag = "query", translation = query),
      reference,
      cache_dir = tempdir(),
      verbose = FALSE
    ),
    dnmb_detect_binary = function(...) list(found = TRUE),
    .dnmb_type1s_prepare_diamond_db = function(...) list(trd1 = "trd1", trd2 = "trd2", full = "full"),
    .dnmb_type1s_prepare_blast_db = function(...) list(trd1 = "trd1", trd2 = "trd2", full = "full"),
    .dnmb_type1s_run_diamond = function(...) stop("forced DIAMOND failure"),
    .dnmb_type1s_run_blast = fake_blast,
    .package = "DNMB"
  )

  expect_equal(prediction$type1s_predicted_recognition, "CAANNNNNNCTC")
  expect_equal(prediction$type1s_search_backend, "blastp")
  expect_equal(prediction$type1s_overall_confidence, "high")
})

test_that("an explicit BLASTP backend does not invoke DIAMOND", {
  sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 20), collapse = "")
  query <- sequence
  substr(query, 210, 210) <- "A"
  windows <- DNMB:::.dnmb_type1s_extract_windows(sequence)
  reference <- list(
    table = data.frame(
      reference_id = "R0001",
      source_enzyme = "S.syntheticI",
      sequence = sequence,
      recognition = "CAANNNNNNCTC",
      left_half = "CAA",
      right_half = "CTC",
      spacer_length = 6L,
      trd1_sequence = windows$trd1,
      trd2_sequence = windows$trd2,
      stringsAsFactors = FALSE
    ),
    version = "test-reference",
    reference_dir = tempdir()
  )
  fake_blast <- function(query_fasta, db_prefix, output_path, cpu = 1L) {
    if (grepl("query_trd", query_fasta, fixed = TRUE)) {
      is_trd1 <- grepl("trd1", db_prefix, fixed = TRUE)
      return(data.frame(
        query_id = if (is_trd1) "Q0001_T1" else "Q0001_T2",
        subject_id = if (is_trd1) "R0001_T1" else "R0001_T2",
        identity = 95,
        query_coverage = 100,
        subject_coverage = 100,
        bitscore = 200,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      query_id = "Q0001",
      subject_id = "R0001",
      identity = 95,
      query_coverage = 100,
      subject_coverage = 100,
      bitscore = 400,
      stringsAsFactors = FALSE
    )
  }

  prediction <- testthat::with_mocked_bindings(
    DNMB:::.dnmb_type1s_predict_core(
      data.frame(locus_tag = "query", translation = query),
      reference,
      cache_dir = tempdir(),
      backend = "blastp",
      verbose = FALSE
    ),
    dnmb_detect_binary = function(...) list(found = TRUE),
    .dnmb_type1s_prepare_diamond_db = function(...) stop("DIAMOND must not be prepared"),
    .dnmb_type1s_prepare_blast_db = function(...) list(trd1 = "trd1", trd2 = "trd2", full = "full"),
    .dnmb_type1s_run_diamond = function(...) stop("DIAMOND must not run"),
    .dnmb_type1s_run_blast = fake_blast,
    .package = "DNMB"
  )
  expect_equal(prediction$type1s_search_backend, "blastp")
  expect_equal(prediction$type1s_overall_confidence, "high")
})

test_that("search failure preserves known rows in a mixed batch", {
  sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 20), collapse = "")
  query <- sequence
  substr(query, 210, 210) <- "A"
  windows <- DNMB:::.dnmb_type1s_extract_windows(sequence)
  reference <- list(
    table = data.frame(
      reference_id = "R0001",
      source_enzyme = "S.syntheticI",
      sequence = sequence,
      recognition = "CAANNNNNNCTC",
      left_half = "CAA",
      right_half = "CTC",
      spacer_length = 6L,
      trd1_sequence = windows$trd1,
      trd2_sequence = windows$trd2,
      stringsAsFactors = FALSE
    ),
    version = "test-reference",
    reference_dir = tempdir()
  )

  prediction <- testthat::with_mocked_bindings(
    DNMB:::.dnmb_type1s_predict_core(
      data.frame(
        locus_tag = c("known", "unsearched"),
        translation = c(sequence, query)
      ),
      reference,
      cache_dir = tempdir(),
      verbose = FALSE
    ),
    dnmb_detect_binary = function(binary, ...) list(found = identical(binary, "diamond")),
    .dnmb_type1s_prepare_diamond_db = function(...) list(trd1 = "trd1", trd2 = "trd2", full = "full"),
    .dnmb_type1s_run_diamond = function(...) stop("forced search failure"),
    .package = "DNMB"
  )

  expect_identical(
    prediction$type1s_prediction_status,
    c("known_gold_match", "search_failed")
  )
  expect_equal(prediction$type1s_predicted_recognition[[1]], "CAANNNNNNCTC")
  expect_true(prediction$type1s_prediction_eligible[[1]])
  expect_false(prediction$type1s_prediction_eligible[[2]])
})

test_that("noncanonical HsdS architectures abstain", {
  short_sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", 10), collapse = "")
  expect_null(DNMB:::.dnmb_type1s_extract_windows(short_sequence))

  confidence <- DNMB:::.dnmb_type1s_confidence(
    list(identity = 95, coverage = 90, subject_coverage = 90),
    list(identity = 96, coverage = 90, subject_coverage = 90),
    list(
      identity = 95,
      coverage = 90,
      subject_coverage = 90,
      method = "whole_hsds_1nn",
      can_be_high = TRUE
    )
  )
  expect_equal(confidence$half, "high")
  expect_equal(confidence$overall, "high")
  expect_true(confidence$eligible)

  one_sided_coverage <- DNMB:::.dnmb_type1s_confidence(
    list(identity = 95, coverage = 95, subject_coverage = 70),
    list(identity = 96, coverage = 95, subject_coverage = 70),
    list(
      identity = 95,
      coverage = 95,
      subject_coverage = 70,
      method = "whole_hsds_1nn",
      can_be_high = TRUE
    )
  )
  expect_false(one_sided_coverage$overall == "high")
  expect_false(one_sided_coverage$eligible)

  weighted_vote <- DNMB:::.dnmb_type1s_confidence(
    list(identity = 95, coverage = 90, subject_coverage = 90),
    list(identity = 96, coverage = 90, subject_coverage = 90),
    list(
      identity = 95,
      coverage = 90,
      subject_coverage = 90,
      method = "whole_hsds_weighted_knn5"
    )
  )
  expect_equal(weighted_vote$spacer, "medium")
  expect_equal(weighted_vote$overall, "medium")
  expect_false(weighted_vote$eligible)

  structural_vote <- DNMB:::.dnmb_type1s_confidence(
    list(identity = 95, coverage = 90, subject_coverage = 90),
    list(identity = 96, coverage = 90, subject_coverage = 90),
    list(
      identity = 99,
      coverage = 99,
      subject_coverage = 99,
      method = "whole_hsds_scaffold_ensemble",
      can_be_high = FALSE
    )
  )
  expect_equal(structural_vote$overall, "medium")
  expect_false(structural_vote$eligible)

  close_tael_ruler <- DNMB:::.dnmb_type1s_confidence(
    list(identity = 95, coverage = 90, subject_coverage = 90),
    list(identity = 96, coverage = 90, subject_coverage = 90),
    list(
      identity = 96,
      coverage = 90,
      subject_coverage = 90,
      method = "whole_hsds_1nn_tael_ruler",
      can_be_high = TRUE
    )
  )
  expect_equal(close_tael_ruler$overall, "high")
  expect_true(close_tael_ruler$eligible)

  remote_tael_ruler <- DNMB:::.dnmb_type1s_confidence(
    list(identity = 95, coverage = 90, subject_coverage = 90),
    list(identity = 96, coverage = 90, subject_coverage = 90),
    list(
      identity = 99,
      coverage = 99,
      subject_coverage = 99,
      method = "whole_hsds_scaffold_tael_ruler",
      can_be_high = FALSE
    )
  )
  expect_equal(remote_tael_ruler$overall, "medium")
  expect_false(remote_tael_ruler$eligible)

  unrecognized_method <- DNMB:::.dnmb_type1s_confidence(
    list(identity = 95, coverage = 90, subject_coverage = 90),
    list(identity = 96, coverage = 90, subject_coverage = 90),
    list(
      identity = 99,
      coverage = 99,
      subject_coverage = 99,
      method = "future_unvalidated_model"
    )
  )
  expect_false(unrecognized_method$eligible)

  methodless_hit <- DNMB:::.dnmb_type1s_confidence(
    list(identity = 95, coverage = 90, subject_coverage = 90),
    list(identity = 96, coverage = 90, subject_coverage = 90),
    list(identity = 99, coverage = 99, subject_coverage = 99)
  )
  expect_false(methodless_hit$eligible)
})

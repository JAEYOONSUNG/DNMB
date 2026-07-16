test_that("run_DNMB exposes the REBASEfinder run-module controls", {
  option_names <- c(
    "rebasefinder_structure_tsv",
    "rebasefinder_structure_max_evalue",
    "rebasefinder_structure_min_probability",
    "rebasefinder_structure_min_tmscore",
    "rebasefinder_homology_modeling",
    "rebasefinder_homology_max_candidates",
    "rebasefinder_typeiii_context_genes"
  )
  run_defaults <- formals(DNMB::run_DNMB)[option_names]
  module_defaults <- formals(DNMB::run_module_set)[option_names]

  expect_named(run_defaults, option_names)
  expect_identical(run_defaults, module_defaults)
})

test_that("module resolver passes every REBASEfinder control unchanged", {
  captured <- NULL
  local_mocked_bindings(
    run_module_set = function(...) {
      captured <<- list(...)
      list(REBASEfinder = "captured")
    },
    .package = "DNMB"
  )
  structure_tsv <- tempfile(fileext = ".tsv")
  writeLines("query\ttarget\ne1\tt1", structure_tsv)
  on.exit(unlink(structure_tsv, force = TRUE), add = TRUE)

  result <- DNMB:::dnmb_resolve_module_results(
    module_aliases = "REBASEfinder",
    rebasefinder_structure_tsv = structure_tsv,
    rebasefinder_structure_max_evalue = 2e-4,
    rebasefinder_structure_min_probability = 0.73,
    rebasefinder_structure_min_tmscore = 0.61,
    rebasefinder_homology_modeling = FALSE,
    rebasefinder_homology_max_candidates = 9L,
    rebasefinder_typeiii_context_genes = 5L
  )

  expect_identical(result$REBASEfinder, "captured")
  expect_identical(captured$rebasefinder_structure_tsv, structure_tsv)
  expect_identical(captured$rebasefinder_structure_max_evalue, 2e-4)
  expect_identical(captured$rebasefinder_structure_min_probability, 0.73)
  expect_identical(captured$rebasefinder_structure_min_tmscore, 0.61)
  expect_false(captured$rebasefinder_homology_modeling)
  expect_identical(captured$rebasefinder_homology_max_candidates, 9L)
  expect_identical(captured$rebasefinder_typeiii_context_genes, 5L)
})

test_that("REBASEfinder stage signature tracks controls and structure TSV content", {
  local_mocked_bindings(
    .dnmb_collect_module_db_signatures = function(...) {
      list(REBASEfinder = list(installed = TRUE, version = "fixture"))
    },
    .package = "DNMB"
  )
  structure_tsv <- tempfile(fileext = ".tsv")
  on.exit(unlink(structure_tsv, force = TRUE), add = TRUE)
  writeLines("query\ttarget\nq1\tt1", structure_tsv)

  signature <- function(max_candidates = 24L) {
    DNMB:::.dnmb_module_stage_signature(
      genbank_signature = list(fixture = TRUE),
      module_aliases = "REBASEfinder",
      rebasefinder_structure_tsv = structure_tsv,
      rebasefinder_structure_max_evalue = 1e-4,
      rebasefinder_structure_min_probability = 0.72,
      rebasefinder_structure_min_tmscore = 0.57,
      rebasefinder_homology_modeling = TRUE,
      rebasefinder_homology_max_candidates = max_candidates,
      rebasefinder_typeiii_context_genes = 4L
    )
  }

  first <- signature()
  expect_identical(first$signature_version, 4L)
  expect_identical(
    first$requested$rebasefinder_homology_backend,
    "promod3-3.6.0-bundled-subclass-templates-v3"
  )
  expect_identical(first$requested$rebasefinder_homology_max_candidates, 24L)
  expect_identical(first$requested$rebasefinder_typeiii_context_genes, 4L)
  expect_equal(first$requested$rebasefinder_structure_min_tmscore, 0.57)
  expect_identical(first$related_inputs$rebasefinder_structure_tsv$file_count, 1L)

  changed_option <- signature(max_candidates = 8L)
  expect_false(identical(first$requested, changed_option$requested))

  writeLines("query\ttarget\nq1\tt2", structure_tsv)
  changed_file <- signature()
  expect_false(identical(
    first$related_inputs$rebasefinder_structure_tsv,
    changed_file$related_inputs$rebasefinder_structure_tsv
  ))
})

test_that("REBASEfinder embedded identity adds Type I-S versions without replacing structural state", {
  cache_root <- tempfile("rebasefinder-type1s-cache-")
  dir.create(cache_root, recursive = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE, force = TRUE), add = TRUE)

  active <- DNMB:::.dnmb_rebasefinder_cache_identity(cache_root = cache_root)
  expect_identical(active$type1s, list(
    signature_version = 2L,
    parser_version = DNMB:::.dnmb_type1s_parser_version(),
    database_version = DNMB:::.dnmb_type1s_database_version(),
    prediction_version = DNMB:::.dnmb_type1s_prediction_version(),
    spacer_model_version = DNMB:::.dnmb_type1s_spacer_model_version(),
    source = NULL,
    reference_key = NA_character_,
    reference_artifacts = NULL
  ))

  versions <- c(
    .dnmb_type1s_parser_version = "type1s-parser-test-v1",
    .dnmb_type1s_database_version = "type1s-db-test-v2",
    .dnmb_type1s_prediction_version = "type1s-predict-test-v3",
    .dnmb_type1s_spacer_model_version = "type1s-spacer-test-v4"
  )
  identity <- testthat::with_mocked_bindings(
    DNMB:::.dnmb_rebasefinder_cache_identity(cache_root = cache_root),
    .dnmb_optional_internal_version = function(helper_name,
                                               fallback = "unavailable") {
      value <- unname(versions[as.character(helper_name)[1]])
      if (!length(value) || is.na(value) || !nzchar(value)) fallback else value
    },
    .package = "DNMB"
  )

  expect_identical(identity$type1s, list(
    signature_version = 2L,
    parser_version = "type1s-parser-test-v1",
    database_version = "type1s-db-test-v2",
    prediction_version = "type1s-predict-test-v3",
    spacer_model_version = "type1s-spacer-test-v4",
    source = NULL,
    reference_key = NA_character_,
    reference_artifacts = NULL
  ))
  expect_identical(identity$module, "rebasefinder")
  expect_identical(identity$version, "embedded")
  expect_true(all(c(
    "files", "homology_template_manifest_md5",
    "promod3_available", "promod3_cli_mtime"
  ) %in% names(identity)))
})

test_that("Type I-S Gold and active reference changes invalidate only its cache identity", {
  wd <- tempfile("rebasefinder-type1s-artifact-stage-")
  cache_root <- tempfile("rebasefinder-type1s-artifact-db-")
  cache_dir <- file.path(cache_root, "db_modules", "rebasefinder", "cache")
  dir.create(wd, recursive = TRUE)
  dir.create(cache_dir, recursive = TRUE)
  on.exit(unlink(c(wd, cache_root), recursive = TRUE, force = TRUE), add = TRUE)

  gold_path <- file.path(
    cache_dir,
    "Type_I_S_subunit_Gold_Standards_Protein.txt"
  )
  writeLines(c(">HsdS_A AAC(6)GTT", "AAAAAAAAAA"), gold_path)
  initial_identity <- DNMB:::.dnmb_rebasefinder_cache_identity(
    cache_root = cache_root
  )
  reference_key <- initial_identity$type1s$reference_key
  expect_match(
    reference_key,
    substr(unname(tools::md5sum(gold_path)), 1L, 12L),
    fixed = TRUE
  )

  reference_dir <- file.path(cache_dir, "type1s", reference_key)
  dir.create(reference_dir, recursive = TRUE)
  reference_path <- file.path(reference_dir, "reference.rds")
  writeBin(charToRaw("reference-A"), reference_path)

  stage_signature <- function() {
    DNMB:::.dnmb_module_stage_signature(
      genbank_signature = list(fixture = TRUE),
      module_aliases = "REBASEfinder",
      module_cache_root = cache_root,
      rebasefinder_homology_modeling = FALSE
    )
  }
  baseline <- stage_signature()
  DNMB:::.dnmb_write_module_stage_cache(
    wd = wd,
    signature = baseline,
    module_results = list(REBASEfinder = data.frame(locus_tag = "fixture"))
  )

  writeBin(charToRaw("reference-B"), reference_path)
  changed_reference <- stage_signature()
  expect_identical(
    baseline$db_state$REBASEfinder$type1s$source,
    changed_reference$db_state$REBASEfinder$type1s$source
  )
  expect_false(identical(
    baseline$db_state$REBASEfinder$type1s$reference_artifacts,
    changed_reference$db_state$REBASEfinder$type1s$reference_artifacts
  ))
  expect_identical(
    baseline$db_state$REBASEfinder[
      setdiff(names(baseline$db_state$REBASEfinder), "type1s")
    ],
    changed_reference$db_state$REBASEfinder[
      setdiff(names(changed_reference$db_state$REBASEfinder), "type1s")
    ]
  )
  expect_identical(
    DNMB:::.dnmb_module_stage_cache_status(wd, changed_reference)$reason,
    "signature_changed"
  )

  writeLines(c(">HsdS_B AAC(6)GTT", "BBBBBBBBBB"), gold_path)
  changed_gold <- stage_signature()
  expect_false(identical(
    baseline$db_state$REBASEfinder$type1s$source,
    changed_gold$db_state$REBASEfinder$type1s$source
  ))
  expect_false(identical(
    baseline$db_state$REBASEfinder$type1s$reference_key,
    changed_gold$db_state$REBASEfinder$type1s$reference_key
  ))
  expect_identical(
    baseline$db_state$REBASEfinder[
      setdiff(names(baseline$db_state$REBASEfinder), "type1s")
    ],
    changed_gold$db_state$REBASEfinder[
      setdiff(names(changed_gold$db_state$REBASEfinder), "type1s")
    ]
  )
  expect_identical(
    DNMB:::.dnmb_module_stage_cache_status(wd, changed_gold)$reason,
    "signature_changed"
  )
})

test_that("Type I-S model changes invalidate cache without changing structural controls", {
  wd <- tempfile("rebasefinder-type1s-stage-")
  cache_root <- tempfile("rebasefinder-type1s-db-")
  dir.create(wd, recursive = TRUE)
  dir.create(cache_root, recursive = TRUE)
  on.exit(unlink(c(wd, cache_root), recursive = TRUE, force = TRUE), add = TRUE)

  make_signature <- function(prediction_version) {
    testthat::with_mocked_bindings(
      DNMB:::.dnmb_module_stage_signature(
        genbank_signature = list(fixture = TRUE),
        module_aliases = "REBASEfinder",
        module_cache_root = cache_root,
        rebasefinder_structure_max_evalue = 1e-4,
        rebasefinder_structure_min_probability = 0.72,
        rebasefinder_structure_min_tmscore = 0.57,
        rebasefinder_homology_modeling = TRUE,
        rebasefinder_homology_max_candidates = 24L,
        rebasefinder_typeiii_context_genes = 4L
      ),
      .dnmb_optional_internal_version = function(helper_name,
                                                 fallback = "unavailable") {
        helper_name <- as.character(helper_name)[1]
        if (identical(helper_name, ".dnmb_type1s_prediction_version")) {
          prediction_version
        } else {
          paste0(helper_name, "-stable")
        }
      },
      .package = "DNMB"
    )
  }

  original <- make_signature("type1s-predict-test-v1")
  changed <- make_signature("type1s-predict-test-v2")
  expect_identical(original$signature_version, 4L)
  expect_identical(original$requested, changed$requested)
  expect_identical(original$related_inputs, changed$related_inputs)
  expect_identical(
    original$db_state$REBASEfinder[
      setdiff(names(original$db_state$REBASEfinder), "type1s")
    ],
    changed$db_state$REBASEfinder[
      setdiff(names(changed$db_state$REBASEfinder), "type1s")
    ]
  )
  expect_identical(
    changed$db_state$REBASEfinder$type1s$prediction_version,
    "type1s-predict-test-v2"
  )

  DNMB:::.dnmb_write_module_stage_cache(
    wd = wd,
    signature = original,
    module_results = list(REBASEfinder = data.frame(locus_tag = "fixture"))
  )
  expect_identical(
    DNMB:::.dnmb_module_stage_cache_status(wd, changed)$reason,
    "signature_changed"
  )
})

test_that("missing Type I-S version helpers have a stable cache fallback", {
  expect_identical(
    DNMB:::.dnmb_optional_internal_version(
      ".dnmb_type1s_version_helper_not_present"
    ),
    "unavailable"
  )
})

test_that("REBASEfinder stage signature tracks an auto-discovered Foldseek TSV", {
  wd <- tempfile("rebasefinder-auto-structure-")
  dir.create(file.path(wd, "dnmb_module_rebasefinder"), recursive = TRUE)
  on.exit(unlink(wd, recursive = TRUE, force = TRUE), add = TRUE)
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)

  local_mocked_bindings(
    .dnmb_collect_module_db_signatures = function(...) list(),
    .package = "DNMB"
  )
  make_signature <- function() {
    DNMB:::.dnmb_module_stage_signature(
      genbank_signature = list(fixture = TRUE),
      module_aliases = "REBASEfinder",
      rebasefinder_structure_tsv = NULL
    )
  }

  before <- make_signature()
  expect_null(before$related_inputs$rebasefinder_structure_tsv)

  auto_tsv <- file.path(wd, "dnmb_module_rebasefinder", "foldseek_results.tsv")
  writeLines("query\ttarget\nq1\tt1", auto_tsv)
  after <- make_signature()
  expect_identical(after$related_inputs$rebasefinder_structure_tsv$file_count, 1L)
  expect_false(identical(
    before$related_inputs$rebasefinder_structure_tsv,
    after$related_inputs$rebasefinder_structure_tsv
  ))
})

test_that("REBASEfinder sentinel is incomplete after model failure or lost model", {
  wd <- tempfile("rebasefinder-stage-model-")
  module_dir <- file.path(wd, "dnmb_module_rebasefinder")
  model_dir <- file.path(module_dir, "promod3_query_structures")
  dir.create(model_dir, recursive = TRUE)
  on.exit(unlink(wd, recursive = TRUE, force = TRUE), add = TRUE)
  DNMB:::.dnmb_module_write_completion_sentinel("REBASEfinder", wd = wd)
  audit_path <- file.path(module_dir, "DNMB_REBASEfinder_homology_templates.tsv")
  signature <- list(requested = list(
    module_install = FALSE,
    rebasefinder_homology_modeling = TRUE
  ))

  audit <- data.frame(
    model_eligible = TRUE,
    model_status = "model_failed:failed_pm",
    model_path = NA_character_,
    model_expected_path = "/data/dnmb_module_rebasefinder/promod3_query_structures/q1.pdb"
  )
  write.table(audit, audit_path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_false(DNMB:::.dnmb_module_stage_artifacts_exist(
    "REBASEfinder", wd = wd, signature = signature
  ))

  audit$model_status <- "template_model_built"
  audit$model_path <- audit$model_expected_path
  write.table(audit, audit_path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_false(DNMB:::.dnmb_module_stage_artifacts_exist(
    "REBASEfinder", wd = wd, signature = signature
  ))

  writeLines("ATOM", file.path(model_dir, "q1.pdb"))
  expect_true(DNMB:::.dnmb_module_stage_artifacts_exist(
    "REBASEfinder", wd = wd, signature = signature
  ))

  signature$requested$rebasefinder_homology_modeling <- FALSE
  unlink(file.path(model_dir, "q1.pdb"))
  expect_true(DNMB:::.dnmb_module_stage_artifacts_exist(
    "REBASEfinder", wd = wd, signature = signature
  ))
})

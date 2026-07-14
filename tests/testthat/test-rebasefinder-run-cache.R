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

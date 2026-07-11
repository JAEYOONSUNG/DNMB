test_that("supported tool and model releases are pinned", {
  expect_identical(DNMB:::.dnmb_dbcan_required_tool_version(), "5.2.9")
  expect_identical(DNMB:::.dnmb_eggnog_required_version(), "2.1.15")
  expect_identical(DNMB:::.dnmb_defensefinder_required_tool_version(), "3.0.0")
  expect_identical(DNMB:::.dnmb_defensefinder_default_repo_ref(), "v3.0.0")
  expect_identical(DNMB:::.dnmb_defensefinder_default_models_version(), "3.1.0")
  expect_identical(
    DNMB:::.dnmb_defensefinder_default_models_repo_url(),
    "https://github.com/mdmparis/defense-finder-models.git"
  )
  expect_identical(DNMB:::.dnmb_defensefinder_default_casfinder_version(), "3.1.0")
  expect_match(DNMB:::.dnmb_dbcan_aws_s3_base(), "db_v5-2-9_5-5-2026$", perl = TRUE)
})

make_defensefinder_models_fixture <- function(path, nested = TRUE, version = "3.1.0") {
  model_dir <- if (nested) file.path(path, "defense-finder-models") else path
  dir.create(file.path(model_dir, "definitions"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(model_dir, "profiles"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("---", paste0("vers: ", version)), file.path(model_dir, "metadata.yml"))
  writeLines(
    c('<model vers="1.0">', '  <gene name="defense_test" presence="mandatory"/>', "</model>"),
    file.path(model_dir, "definitions", "DefenseTest.xml")
  )
  writeLines("HMMER3/f", file.path(model_dir, "profiles", "defense_test.hmm"))
  invisible(path)
}

make_defensefinder_cas_fixture <- function(path, version = "3.1.0") {
  dir.create(file.path(path, "definitions"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "profiles"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("---", paste0("vers: ", version)), file.path(path, "metadata.yml"))
  writeLines(
    c('<model vers="2.0">', '  <gene name="cas9_test" presence="mandatory"/>', "</model>"),
    file.path(path, "definitions", "CAS_Class2-Subtype-VI-X.xml")
  )
  writeLines("HMMER3/f", file.path(path, "profiles", "cas9_test.hmm"))
  invisible(path)
}

test_that("DefenseFinder model release is read from macsydata metadata", {
  models_dir <- tempfile("defensefinder-models-")
  dir.create(file.path(models_dir, "defense-finder-models"), recursive = TRUE)
  on.exit(unlink(models_dir, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(
    c("---", "short_desc: test fixture", "vers: 3.1.0"),
    file.path(models_dir, "defense-finder-models", "metadata.yml")
  )

  expect_identical(
    DNMB:::.dnmb_defensefinder_models_release_version(models_dir),
    "3.1.0"
  )
})

test_that("CasFinder compatibility rejects incomplete profile bundles", {
  casfinder_dir <- tempfile("casfinder-models-")
  definitions_dir <- file.path(casfinder_dir, "definitions")
  profiles_dir <- file.path(casfinder_dir, "profiles")
  dir.create(definitions_dir, recursive = TRUE)
  dir.create(profiles_dir, recursive = TRUE)
  on.exit(unlink(casfinder_dir, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(
    c('<model vers="2.0">', '  <gene name="cas9_test" presence="mandatory"/>', "</model>"),
    file.path(definitions_dir, "CAS_Class2-Subtype-VI-X.xml")
  )

  expect_false(DNMB:::.dnmb_defensefinder_casfinder_is_compatible(casfinder_dir))
  writeLines("HMMER3/f", file.path(profiles_dir, "cas9_test.hmm"))
  expect_true(DNMB:::.dnmb_defensefinder_casfinder_is_compatible(casfinder_dir))
})

test_that("DefenseFinder uses the v3 version subcommand and tagged clone", {
  cli <- tempfile("defense-finder-")
  writeLines("#!/bin/sh", cli)
  on.exit(unlink(cli, force = TRUE), add = TRUE)
  calls <- list()
  local_mocked_bindings(
    dnmb_run_external = function(command, args = character(), ...) {
      calls[[length(calls) + 1L]] <<- list(command = command, args = args)
      list(
        ok = TRUE,
        status = 0L,
        stdout = if (identical(command, cli)) "Using DefenseFinder version 3.0.0" else character(),
        stderr = character(),
        error = NULL,
        resolved_command = command
      )
    },
    .package = "DNMB"
  )

  expect_identical(
    DNMB:::.dnmb_defensefinder_tool_version(list(cli_path = cli, env_python = "")),
    "3.0.0"
  )
  expect_identical(calls[[1]]$args, "version")

  layout <- DNMB:::.dnmb_defensefinder_asset_layout(tempfile("defensefinder-cache-"))
  status <- DNMB:::.dnmb_defensefinder_prepare_repo(layout)
  expect_identical(status$status, "ok")
  clone_args <- calls[[2]]$args
  expect_true(all(c("clone", "--branch", "v3.0.0") %in% clone_args))
})

test_that("DefenseFinder model update prefers validated local model assets", {
  module_dir <- tempfile("defensefinder-module-")
  layout <- DNMB:::.dnmb_defensefinder_asset_layout(module_dir)
  models_source <- tempfile("defensefinder-model-source-")
  casfinder_source <- tempfile("casfinder-model-source-")
  make_defensefinder_models_fixture(models_source)
  make_defensefinder_cas_fixture(casfinder_source)
  on.exit(unlink(c(module_dir, models_source, casfinder_source), recursive = TRUE, force = TRUE), add = TRUE)

  local_mocked_bindings(
    dnmb_run_external = function(...) stop("local assets should avoid external installers"),
    .package = "DNMB"
  )
  status <- DNMB:::.dnmb_defensefinder_update_models(
    layout,
    asset_urls = list(models_dir = models_source, casfinder_dir = casfinder_source)
  )

  expect_identical(status$status, "ok")
  expect_identical(DNMB:::.dnmb_defensefinder_models_release_version(layout$models_dir), "3.1.0")
  expect_true(DNMB:::.dnmb_defensefinder_casfinder_is_required(file.path(layout$models_dir, "CasFinder")))
  expect_match(status$detail, "defense-finder-models <=", fixed = TRUE)
  expect_match(status$detail, "CasFinder <=", fixed = TRUE)
})

test_that("DefenseFinder falls back to the official tagged model clone", {
  module_dir <- tempfile("defensefinder-module-")
  layout <- DNMB:::.dnmb_defensefinder_asset_layout(module_dir)
  dir.create(dirname(layout$macsydata_path), recursive = TRUE, showWarnings = FALSE)
  writeLines("fixture", layout$macsydata_path)
  casfinder_source <- tempfile("casfinder-model-source-")
  make_defensefinder_cas_fixture(casfinder_source)
  on.exit(unlink(c(module_dir, casfinder_source), recursive = TRUE, force = TRUE), add = TRUE)
  calls <- list()

  local_mocked_bindings(
    dnmb_detect_binary = function(binary, required = FALSE) list(found = TRUE, path = "/test/hmmsearch"),
    dnmb_run_external = function(command, args = character(), ...) {
      calls[[length(calls) + 1L]] <<- list(command = command, args = args)
      if (identical(command, "git")) {
        make_defensefinder_models_fixture(args[[length(args)]], nested = FALSE)
        return(list(ok = TRUE, error = NULL, resolved_command = command))
      }
      list(ok = FALSE, error = "GitHub API rate limit", resolved_command = command)
    },
    .package = "DNMB"
  )
  status <- DNMB:::.dnmb_defensefinder_update_models(
    layout,
    asset_urls = list(casfinder_dir = casfinder_source)
  )

  expect_identical(status$status, "ok")
  expect_identical(calls[[1]]$command, layout$macsydata_path)
  expect_true("defense-finder-models==3.1.0" %in% calls[[1]]$args)
  clone_call <- Filter(function(x) identical(x$command, "git"), calls)[[1]]
  expect_true(all(c("clone", "--depth", "1", "--branch", "3.1.0") %in% clone_call$args))
  expect_true(DNMB:::.dnmb_defensefinder_default_models_repo_url() %in% clone_call$args)
  expect_identical(DNMB:::.dnmb_defensefinder_models_release_version(layout$models_dir), "3.1.0")
})

test_that("DefenseFinder keeps the prior model cache when staging fails", {
  module_dir <- tempfile("defensefinder-module-")
  layout <- DNMB:::.dnmb_defensefinder_asset_layout(module_dir)
  dir.create(layout$models_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(layout$models_dir, "existing-cache-marker")
  writeLines("keep", marker)
  on.exit(unlink(module_dir, recursive = TRUE, force = TRUE), add = TRUE)

  local_mocked_bindings(
    dnmb_run_external = function(command, args = character(), ...) {
      list(ok = FALSE, error = "offline", resolved_command = command)
    },
    .package = "DNMB"
  )
  status <- DNMB:::.dnmb_defensefinder_update_models(layout)

  expect_identical(status$status, "failed")
  expect_true(file.exists(marker))
  expect_identical(readLines(marker), "keep")
})

test_that("DefenseFinder installer forwards local model assets to the updater", {
  install_body <- paste(deparse(body(DNMB:::dnmb_defensefinder_install_module)), collapse = "\n")
  expect_match(
    install_body,
    "\\.dnmb_defensefinder_update_models\\(layout,\\s+asset_urls = asset_urls\\)",
    perl = TRUE
  )
})

test_that("dbCAN uses the 5.2 version subcommand", {
  calls <- list()
  local_mocked_bindings(
    dnmb_detect_binary = function(binary, required = FALSE) {
      list(found = TRUE, path = "/test/run_dbcan")
    },
    dnmb_run_external = function(command, args = character(), ...) {
      calls[[length(calls) + 1L]] <<- list(command = command, args = args)
      list(ok = TRUE, stdout = "dbCAN version: 5.2.9", stderr = character())
    },
    .package = "DNMB"
  )

  expect_identical(DNMB:::.dnmb_dbcan_tool_version(), "5.2.9")
  expect_identical(calls[[1]]$args, "version")
})

test_that("eggNOG auto-install prefers the pinned Bioconda package", {
  calls <- list()
  local_mocked_bindings(
    dnmb_detect_binary = function(binary, required = FALSE) {
      list(found = identical(binary, "mamba"), path = if (identical(binary, "mamba")) binary else "")
    },
    dnmb_run_external = function(command, args = character(), ...) {
      calls[[length(calls) + 1L]] <<- list(command = command, args = args)
      list(ok = TRUE, stdout = character(), stderr = character())
    },
    dnmb_detect_emapper = function(required = FALSE) {
      list(found = TRUE, path = "/test/emapper.py", version = "2.1.15")
    },
    .package = "DNMB"
  )

  installed <- DNMB:::.dnmb_eggnog_auto_install_emapper(verbose = FALSE)
  expect_true(installed$found)
  expect_identical(calls[[1]]$command, "mamba")
  expect_true("eggnog-mapper=2.1.15" %in% calls[[1]]$args)
})

test_that("eggNOG uses Conda release metadata over the upstream stale CLI label", {
  env_dir <- tempfile("eggnog-conda-env-")
  bin_dir <- file.path(env_dir, "bin")
  metadata_dir <- file.path(env_dir, "conda-meta")
  dir.create(bin_dir, recursive = TRUE)
  dir.create(metadata_dir, recursive = TRUE)
  on.exit(unlink(env_dir, recursive = TRUE, force = TRUE), add = TRUE)
  emapper_path <- file.path(bin_dir, "emapper.py")
  writeLines("#!/usr/bin/env python3", emapper_path)
  jsonlite::write_json(
    list(name = "eggnog-mapper", version = "2.1.15"),
    file.path(metadata_dir, "eggnog-mapper-2.1.15-test.json"),
    auto_unbox = TRUE
  )

  local_mocked_bindings(
    dnmb_detect_binary = function(binary, required = FALSE) {
      list(found = TRUE, path = emapper_path)
    },
    dnmb_run_external = function(command, args = character(), ...) {
      list(ok = TRUE, stdout = "emapper-2.1.13", stderr = character())
    },
    .package = "DNMB"
  )

  detection <- DNMB:::dnmb_detect_emapper(required = TRUE)
  expect_identical(detection$cli_version, "2.1.13")
  expect_identical(detection$package_version, "2.1.15")
  expect_identical(detection$version, "2.1.15")
})

test_that("DefenseFinder v3 duplicate system rows do not multiply gene hits", {
  systems_path <- tempfile("defensefinder-systems-", fileext = ".tsv")
  genes_path <- tempfile("defensefinder-genes-", fileext = ".tsv")
  on.exit(unlink(c(systems_path, genes_path), force = TRUE), add = TRUE)

  systems <- data.frame(
    sys_id = c("rep_RM_1", "rep_RM_1"),
    type = c("RM", "RM"),
    subtype = c("RM_Type_I", "RM_Type_I"),
    activity = c("Defense", "Defense"),
    sys_beg = c("rep_0001", "rep_0001"),
    sys_end = c("rep_0002", "rep_0002"),
    protein_in_syst = c("rep_0001,rep_0002", "rep_0001,rep_0002"),
    genes_count = c(2L, 2L),
    name_of_profiles_in_sys = c("MTase,REase", "MTase,REase")
  )
  utils::write.table(systems, systems_path, sep = "\t", quote = FALSE, row.names = FALSE)

  genes <- data.frame(
    replicon = "rep",
    hit_id = "rep_0001",
    gene_name = "MTase",
    hit_pos = 1L,
    model_fqn = "defense-finder-models/RM/RM/RM_Type_I",
    sys_id = "rep_RM_1",
    hit_status = "mandatory"
  )
  utils::write.table(genes, genes_path, sep = "\t", quote = FALSE, row.names = FALSE)
  id_map <- data.frame(
    defensefinder_id = "rep_0001",
    locus_tag = "locus_1",
    contig = "rep",
    start = 1L,
    end = 100L,
    direction = "+",
    product = "methyltransferase"
  )

  parsed_systems <- DNMB:::dnmb_defensefinder_parse_systems(systems_path)
  parsed_genes <- DNMB:::dnmb_defensefinder_parse_genes(genes_path, id_map, parsed_systems)
  expect_equal(nrow(parsed_systems), 1L)
  expect_equal(nrow(parsed_genes), 1L)
})

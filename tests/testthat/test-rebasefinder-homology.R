.homology_template_sequence <- function(template_id) {
  templates <- DNMB:::.dnmb_rebasefinder_read_homology_templates()
  row <- templates[templates$template_id == template_id, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  DNMB:::.dnmb_rebasefinder_read_single_fasta(row$fasta_path[[1]])
}

.homology_test_hit <- function(query, family, role, tier = "high") {
  data.frame(
    query = query,
    family_id = family,
    enzyme_role = role,
    final_family = family,
    final_role = role,
    curation_tier = tier,
    curation_score = if (tier == "high") 90 else 60,
    independent_evidence_axes = if (tier == "high") 3L else 2L,
    support = "",
    stringsAsFactors = FALSE
  )
}

.write_copy_template_pm <- function(path) {
  writeLines(c(
    "#!/bin/sh",
    "template=''",
    "output=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    -p) template=\"$2\"; shift 2 ;;",
    "    -o) output=\"$2\"; shift 2 ;;",
    "    *) shift ;;",
    "  esac",
    "done",
    "cp \"$template\" \"$output\""
  ), path)
  Sys.chmod(path, mode = "0755")
}

test_that("bundled R-M homology templates include HsdR and explicit noise decoys", {
  templates <- DNMB:::.dnmb_rebasefinder_read_homology_templates()

  expect_equal(nrow(templates), 16L)
  expect_equal(sum(templates$template_class == "positive"), 13L)
  expect_equal(sum(templates$template_class == "decoy"), 3L)
  expect_true("EcoR124I_HsdR_6H2J_B" %in% templates$template_id)
  expect_true(any(grepl("rRNA", templates$family_description[templates$template_class == "decoy"], ignore.case = TRUE)))
  expect_true(any(grepl("helicase", templates$family_description[templates$template_class == "decoy"], ignore.case = TRUE)))
  expect_true(any(grepl("repair", templates$description[templates$template_class == "decoy"], ignore.case = TRUE)))
  expect_true(all(file.exists(templates$pdb_path)))
  expect_true(all(file.exists(templates$fasta_path)))
})

test_that("Type I geometry rules do not match Type II families", {
  rules <- DNMB:::.dnmb_rebasefinder_structure_pair_rules()
  hsdr <- rules[grepl("^hsdr_", rules$pair_id), , drop = FALSE]
  type_ii <- rules[rules$pair_id == "mmei_mtase_I-IV", , drop = FALSE]

  expect_false(any(vapply(hsdr$family_pattern, grepl, logical(1), x = "Type II")))
  expect_true(grepl(type_ii$family_pattern, "Type II"))
})

test_that("Type II subtypes are canonicalized before homology motif-pair matching", {
  expect_identical(
    DNMB:::.dnmb_rebasefinder_canonical_structure_family(
      c("Type IIS", "Type IIG", "Type IIL", "Type ISP")
    ),
    c("Type II", "Type II", "Type II", "Type I")
  )

  cases <- data.frame(
    template_id = c("MmeI_RM_5HR4_J", "LlaBIII_RM_4XQK_A"),
    family = c("Type IIL", "Type ISP"),
    expected_pair = c("mmei_mtase_I-IV", "hsdr_motor_WB-MIII"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(cases))) {
    sequence <- .homology_template_sequence(cases$template_id[[i]])
    motif_hits <- DNMB:::.dnmb_rebasefinder_homology_motif_hits(
      cases$template_id[[i]], sequence, cases$family[[i]], "RM"
    )
    alignment <- DNMB:::.dnmb_rebasefinder_template_alignment(sequence, sequence)
    mapping <- DNMB:::.dnmb_rebasefinder_homology_motif_mapping(motif_hits, alignment)

    expect_true(mapping$complete, info = cases$template_id[[i]])
    expect_match(mapping$pairs, cases$expected_pair[[i]], fixed = TRUE)
  }
})

test_that("Type I HsdR maps to the positive template and fails softly without pm", {
  sequence <- .homology_template_sequence("EcoR124I_HsdR_6H2J_B")
  hits <- .homology_test_hit("hsdr_test", "Type I", "R")
  genes <- data.frame(
    locus_tag = "hsdr_test",
    translation = sequence,
    product = "Type I restriction enzyme R subunit",
    stringsAsFactors = FALSE
  )
  output <- tempfile("homology-no-pm-")
  cache <- tempfile("homology-cache-")
  dir.create(output)
  dir.create(cache)
  stale_dir <- file.path(output, "promod3_query_structures")
  dir.create(stale_dir)
  templates <- DNMB:::.dnmb_rebasefinder_read_homology_templates()
  stale_template <- templates$pdb_path[templates$template_id == "EcoR124I_HsdR_6H2J_B"]
  stale_model <- file.path(stale_dir, "hsdr_test.pdb")
  DNMB:::.dnmb_rebasefinder_prepare_template_pdb(stale_template, stale_model)
  expect_true(file.exists(stale_model))
  on.exit(unlink(c(output, cache), recursive = TRUE, force = TRUE), add = TRUE)

  old_path <- Sys.getenv("PATH")
  empty_bin <- tempfile("empty-bin-")
  dir.create(empty_bin)
  on.exit({
    Sys.setenv(PATH = old_path)
    unlink(empty_bin, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(PATH = empty_bin)

  result <- DNMB:::.dnmb_rebasefinder_run_homology_models(
    hits, genes, output, cache_root = cache, install = FALSE, verbose = FALSE
  )

  expect_identical(result$audit$template_id, "EcoR124I_HsdR_6H2J_B")
  expect_equal(result$audit$alignment_identity, 1)
  expect_equal(result$audit$template_coverage, 1)
  expect_identical(result$audit$selection_status, "template_mapped")
  expect_identical(result$audit$model_status, "backend_unavailable")
  expect_true(is.na(result$audit$model_path))
  expect_match(result$audit$model_expected_path, "hsdr_test[.]pdb$")
  expect_false(file.exists(stale_model))
  expect_match(result$status, "partial")
  expect_true(file.exists(result$files$tsv))
})

test_that("managed ProMod3 install falls back from a broken mamba to conda", {
  skip_on_os("windows")
  base <- tempfile("promod3-installer-")
  bin <- file.path(base, "bin")
  cache <- file.path(base, "cache")
  dir.create(bin, recursive = TRUE)
  dir.create(cache)
  on.exit(unlink(base, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines(c("#!/bin/sh", "exit 1"), file.path(bin, "mamba"))
  writeLines(c(
    "#!/bin/sh",
    "echo 'Filesystem 1024-blocks Used Available Capacity Mounted on'",
    "echo 'testfs 30000000 1000000 20000000 5% /'"
  ), file.path(bin, "df"))
  writeLines(c(
    "#!/bin/sh",
    "prefix=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  if [ \"$1\" = '-p' ]; then prefix=\"$2\"; shift 2; else shift; fi",
    "done",
    "mkdir -p \"$prefix/bin\"",
    "printf '#!/bin/sh\\necho build-model\\nexit 0\\n' > \"$prefix/bin/pm\"",
    "chmod +x \"$prefix/bin/pm\""
  ), file.path(bin, "conda"))
  Sys.chmod(file.path(bin, c("mamba", "conda", "df")), mode = "0755")

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  result <- DNMB:::.dnmb_rebasefinder_promod3_cli(
    cache_root = cache, install = TRUE, verbose = FALSE
  )

  expect_identical(result$status, "available")
  expect_true(file.exists(result$path))
  expect_match(result$detail, "conda$")
  layout <- DNMB:::.dnmb_rebasefinder_promod3_layout(cache_root = cache, create = FALSE)
  expect_true(file.exists(layout$ready))
  expect_true(DNMB:::.dnmb_rebasefinder_promod3_ready(layout))
})

test_that("managed ProMod3 requires a completion sentinel and promotes a valid legacy environment", {
  skip_on_os("windows")
  base <- paste0(tempfile("promod3-readiness-"), " with spaces")
  empty_bin <- file.path(base, "empty-bin")
  partial_cache <- file.path(base, "partial-cache")
  legacy_cache <- file.path(base, "legacy-cache")
  dir.create(empty_bin, recursive = TRUE)
  on.exit(unlink(base, recursive = TRUE, force = TRUE), add = TRUE)

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = empty_bin)
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  make_pm <- function(layout) {
    dir.create(dirname(layout$cli), recursive = TRUE)
    writeLines(c("#!/bin/sh", "echo build-model", "exit 0"), layout$cli)
    Sys.chmod(layout$cli, mode = "0755")
  }

  partial <- DNMB:::.dnmb_rebasefinder_promod3_layout(partial_cache, create = TRUE)
  make_pm(partial)
  dir.create(partial$lock)
  partial_result <- DNMB:::.dnmb_rebasefinder_promod3_cli(
    cache_root = partial_cache, install = FALSE, verbose = FALSE
  )
  expect_identical(partial_result$status, "backend_unavailable")
  expect_false(file.exists(partial$ready))
  expect_false(DNMB:::.dnmb_rebasefinder_cache_identity(partial_cache)$promod3_available)

  legacy <- DNMB:::.dnmb_rebasefinder_promod3_layout(legacy_cache, create = TRUE)
  make_pm(legacy)
  legacy_result <- DNMB:::.dnmb_rebasefinder_promod3_cli(
    cache_root = legacy_cache, install = FALSE, verbose = FALSE
  )
  expect_identical(legacy_result$status, "available")
  expect_match(legacy_result$detail, "validated legacy")
  expect_true(file.exists(legacy$ready))
  expect_true(DNMB:::.dnmb_rebasefinder_cache_identity(legacy_cache)$promod3_available)
})

test_that("an rRNA methylase decoy cannot enter ProMod3 modeling", {
  sequence <- .homology_template_sequence("KsgA_rRNA_MTase_1QYR_A")
  hits <- .homology_test_hit("ksga_noise", "Type II", "M", tier = "medium")
  genes <- data.frame(
    locus_tag = "ksga_noise",
    translation = sequence,
    product = "putative DNA methyltransferase",
    stringsAsFactors = FALSE
  )
  output <- tempfile("homology-decoy-")
  cache <- tempfile("homology-cache-")
  dir.create(output)
  dir.create(cache)
  on.exit(unlink(c(output, cache), recursive = TRUE, force = TRUE), add = TRUE)

  result <- DNMB:::.dnmb_rebasefinder_run_homology_models(
    hits, genes, output, cache_root = cache, install = FALSE, verbose = FALSE
  )

  expect_identical(result$audit$best_decoy_id, "KsgA_rRNA_MTase_1QYR_A")
  expect_lt(result$audit$decoy_margin, 0)
  expect_false(result$audit$model_eligible)
  expect_identical(result$audit$model_status, "not_eligible")
})

test_that("unclassified DNA methylases can be resolved by positive templates", {
  sequence <- .homology_template_sequence("HhaI_MTase_2HMY_B")
  hits <- .homology_test_hit("unclassified_mtase", "Unclassified DNA MTase", "M")
  genes <- data.frame(
    locus_tag = "unclassified_mtase",
    translation = sequence,
    product = "putative DNA methyltransferase",
    stringsAsFactors = FALSE
  )
  output <- tempfile("homology-unclassified-")
  cache <- tempfile("homology-cache-")
  dir.create(output)
  dir.create(cache)
  on.exit(unlink(c(output, cache), recursive = TRUE, force = TRUE), add = TRUE)

  result <- DNMB:::.dnmb_rebasefinder_run_homology_models(
    hits, genes, output, cache_root = cache, install = FALSE, verbose = FALSE
  )

  expect_identical(result$audit$template_id, "HhaI_MTase_2HMY_B")
  expect_identical(result$audit$selection_status, "template_mapped")
  expect_identical(result$audit$motif_family_used, "Type II")
  expect_true(result$audit$motif_mapping_complete)
  expect_gt(result$audit$decoy_margin, 0)
})

test_that("unclassified R candidates inherit motif family from the winning template", {
  sequence <- .homology_template_sequence("EcoR124I_HsdR_6H2J_B")
  hits <- .homology_test_hit("unknown_hsdr", "unknown", "R")
  genes <- data.frame(
    locus_tag = "unknown_hsdr", translation = sequence,
    product = "putative restriction subunit", stringsAsFactors = FALSE
  )
  output <- tempfile("homology-unknown-r-")
  cache <- tempfile("homology-cache-")
  dir.create(output)
  dir.create(cache)
  on.exit(unlink(c(output, cache), recursive = TRUE, force = TRUE), add = TRUE)

  result <- DNMB:::.dnmb_rebasefinder_run_homology_models(
    hits, genes, output, cache_root = cache, install = FALSE, verbose = FALSE
  )

  expect_identical(result$audit$template_id, "EcoR124I_HsdR_6H2J_B")
  expect_identical(result$audit$motif_family_used, "Type I")
  expect_match(result$audit$mapped_motifs, "HsdR-PD")
  expect_true(result$audit$motif_mapping_complete)
})

test_that("candidate cap reserves space for restriction roles", {
  methylases <- do.call(rbind, lapply(seq_len(24L), function(i) {
    x <- .homology_test_hit(paste0("m", i), "Type II", "M")
    x$final_role <- "M"
    x$curation_score <- 100
    x
  }))
  restriction <- .homology_test_hit("r_priority", "Type I", "R")
  restriction$final_role <- "R"
  restriction$curation_score <- 99
  hits <- rbind(methylases, restriction)

  keep <- DNMB:::.dnmb_rebasefinder_homology_candidate_keep(hits, max_candidates = 24L)

  expect_equal(sum(keep), 24L)
  expect_true(keep[hits$query == "r_priority"])
})

test_that("stage cache retries an unavailable managed backend", {
  wd <- tempfile("rebasefinder-stage-")
  module_dir <- file.path(wd, "dnmb_module_rebasefinder")
  dir.create(module_dir, recursive = TRUE)
  on.exit(unlink(wd, recursive = TRUE, force = TRUE), add = TRUE)
  DNMB:::.dnmb_module_write_completion_sentinel("REBASEfinder", wd = wd)
  write.table(
    data.frame(model_eligible = TRUE, model_status = "backend_unavailable"),
    file.path(module_dir, "DNMB_REBASEfinder_homology_templates.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  install_signature <- list(requested = list(module_install = TRUE))
  offline_signature <- list(requested = list(module_install = FALSE))
  expect_false(DNMB:::.dnmb_module_stage_artifacts_exist(
    "REBASEfinder", wd = wd, signature = install_signature
  ))
  expect_true(DNMB:::.dnmb_module_stage_artifacts_exist(
    "REBASEfinder", wd = wd, signature = offline_signature
  ))
})

test_that("user structures take precedence over managed homology models", {
  root <- tempfile("rebasefinder-structures-")
  user_dir <- file.path(root, "query_structures")
  homology_dir <- file.path(root, "dnmb_module_rebasefinder", "promod3_query_structures")
  dir.create(user_dir, recursive = TRUE)
  dir.create(homology_dir, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  dirs <- normalizePath(DNMB:::.dnmb_rebasefinder_structure_dirs(root))
  expect_lt(match(normalizePath(user_dir), dirs), match(normalizePath(homology_dir), dirs))

  managed_collision <- file.path(homology_dir, "dup_a__collision_1.pdb")
  user_model <- file.path(user_dir, "dup_a.pdb")
  writeLines("managed", managed_collision)
  writeLines("user", user_model)
  explicit <- setNames(managed_collision, "dup/a")
  expect_identical(
    DNMB:::.dnmb_rebasefinder_structure_file_for_query("dup/a", dirs, explicit),
    normalizePath(user_model)
  )
  unlink(user_model)
  expect_identical(
    DNMB:::.dnmb_rebasefinder_structure_file_for_query("dup/a", dirs, explicit),
    normalizePath(managed_collision)
  )

  prefixed <- data.frame(
    locus_tag = "dup/a",
    REBASEfinder_homology_model_path = managed_collision,
    stringsAsFactors = FALSE
  )
  expect_identical(
    unname(DNMB:::.dnmb_rebasefinder_homology_structure_path_map(prefixed)),
    managed_collision
  )
})

test_that("managed ProMod3 outputs are cleared without touching other files", {
  root <- tempfile("promod3-clean-")
  model_dir <- file.path(root, "promod3_query_structures")
  dir.create(model_dir, recursive = TRUE)
  writeLines("stale", file.path(model_dir, "old.pdb"))
  writeLines("keep", file.path(model_dir, "notes.txt"))
  writeLines("stale", file.path(root, "DNMB_REBASEfinder_homology_templates.tsv"))
  writeLines("stale", file.path(root, "DNMB_REBASEfinder_homology_templates.xlsx"))
  dir.create(file.path(root, "promod3_work"))
  writeLines("stale", file.path(root, "promod3_work", "promod3_jobs.tsv"))
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  result <- DNMB:::.dnmb_rebasefinder_clear_promod3_outputs(root)

  expect_identical(normalizePath(result), normalizePath(model_dir))
  expect_false(file.exists(file.path(model_dir, "old.pdb")))
  expect_true(file.exists(file.path(model_dir, "notes.txt")))
  expect_false(file.exists(file.path(root, "DNMB_REBASEfinder_homology_templates.tsv")))
  expect_false(file.exists(file.path(root, "DNMB_REBASEfinder_homology_templates.xlsx")))
  expect_false(dir.exists(file.path(root, "promod3_work")))
})

test_that("identical paralogs share one model job and fan out to unique paths", {
  skip_on_os("windows")
  sequence <- .homology_template_sequence("EcoR124I_HsdR_6H2J_B")
  hits <- rbind(
    .homology_test_hit("dup/a", "Type I", "R"),
    .homology_test_hit("dup_a", "Type I", "R")
  )
  hits$hit_label <- c("R.DupSlash", "R.DupUnderscore")
  genes <- data.frame(
    locus_tag = hits$query, translation = sequence,
    product = "Type I restriction enzyme R subunit", stringsAsFactors = FALSE
  )
  base <- tempfile("homology-paralogs-")
  bin <- file.path(base, "bin")
  output <- file.path(base, "output")
  cache <- file.path(base, "cache")
  dir.create(bin, recursive = TRUE)
  dir.create(output)
  dir.create(cache)
  fake_pm <- file.path(bin, "pm")
  .write_copy_template_pm(fake_pm)
  on.exit(unlink(base, recursive = TRUE, force = TRUE), add = TRUE)
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  result <- DNMB:::.dnmb_rebasefinder_run_homology_models(
    hits, genes, output, cache_root = cache, install = FALSE, cpu = 2L, verbose = FALSE
  )
  manifest <- read.delim(result$files$model_manifest, check.names = FALSE)

  expect_equal(nrow(result$audit), 2L)
  expect_true(all(result$audit$model_status == "template_model_built"))
  expect_equal(length(unique(result$audit$model_path)), 2L)
  expect_true(all(file.exists(result$audit$model_path)))
  expect_false(any(result$audit$geometry_status == "no_structure"))
  expect_true(all(!is.na(result$audit$geometry_status)))
  explicit <- setNames(result$audit$model_path, result$audit$query)
  resolved <- vapply(
    result$audit$query,
    DNMB:::.dnmb_rebasefinder_structure_file_for_query,
    character(1),
    dirs = dirname(unique(result$audit$model_path)),
    explicit_paths = explicit
  )
  expect_identical(unname(resolved), unname(result$audit$model_path))
  coverage <- DNMB:::.dnmb_rebasefinder_write_structure_coverage(
    output, genes, result$hits
  )
  expect_equal(coverage$n_queries, 2L)
  expect_equal(coverage$n_structure_files, 2L)
  expect_true(all(coverage$table$structure_file_exists))
  expect_setequal(coverage$table$structure_file, result$audit$model_path)
  expect_equal(nrow(manifest), 1L)
})

test_that("ProMod3 models use whitespace-safe staging and cross-locus cache", {
  skip_on_os("windows")
  sequence <- .homology_template_sequence("EcoR124I_HsdR_6H2J_B")
  substr(sequence, 328L, 350L) <- strrep("A", 23L)
  hits <- .homology_test_hit("hsdr_model", "Type I", "R")
  genes <- data.frame(
    locus_tag = "hsdr_model",
    translation = sequence,
    product = "Type I restriction enzyme R subunit",
    stringsAsFactors = FALSE
  )
  base <- paste0(tempfile("homology-model-"), " with spaces")
  bin <- file.path(base, "bin")
  output <- file.path(base, "output")
  cache <- file.path(base, "cache")
  dir.create(bin, recursive = TRUE)
  dir.create(output)
  dir.create(cache)
  fake_pm <- file.path(bin, "pm")
  .write_copy_template_pm(fake_pm)
  on.exit(unlink(base, recursive = TRUE, force = TRUE), add = TRUE)

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  first <- DNMB:::.dnmb_rebasefinder_run_homology_models(
    hits, genes, output, cache_root = cache, install = FALSE, verbose = FALSE
  )
  expect_identical(first$audit$model_status, "template_model_built")
  expect_false(first$audit$model_supported)
  expect_identical(first$audit$geometry_status, "homology_model_motor_only")
  expect_true(file.exists(first$audit$model_path))
  expect_identical(
    DNMB:::.dnmb_rebasefinder_structure_provenance(first$audit$model_path),
    "promod3_homology"
  )

  second_output <- file.path(base, "output_cached")
  dir.create(second_output)
  second_hits <- hits
  second_hits$query <- "hsdr_model_other"
  second_genes <- genes
  second_genes$locus_tag <- "hsdr_model_other"
  second <- DNMB:::.dnmb_rebasefinder_run_homology_models(
    second_hits, second_genes, second_output, cache_root = cache, install = FALSE, verbose = FALSE
  )
  expect_identical(second$audit$model_status, "template_model_built")
  expect_true(second$audit$model_cache_hit)
  expect_false(second$audit$model_supported)
})

test_that("homology support breaks score ties without adding an evidence axis", {
  hits <- data.frame(
    query = c("without_model", "with_model"),
    source = "rebasefinder",
    family_system = "REBASEfinder",
    family_id = "Type I",
    hit_label = c("R.A", "R.B"),
    enzyme_role = "R",
    evidence_mode = "high_confidence",
    substrate_label = NA_character_,
    support = NA_character_,
    typing_eligible = TRUE,
    blast_identity = 0.8,
    blast_evalue = 1e-30,
    blast_bitscore = 250,
    blast_length = 900,
    blast_alignment_quality = "strong",
    blast_role_compatible = TRUE,
    blast_family_compatible = TRUE,
    blast_min_coverage = 0.9,
    blast_query_coverage = 0.9,
    blast_reference_coverage = 0.9,
    partial_status = "full_length",
    mtase_motif_verified = FALSE,
    rease_operon_component_raw = "R_motor_nuclease",
    typei_operon_supported = TRUE,
    typeii_operon_supported = FALSE,
    typeiii_operon_supported = FALSE,
    reference_rec_seq = NA_character_,
    homology_model_supported = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  genes <- data.frame(
    locus_tag = hits$query,
    product = "Type I restriction enzyme R subunit",
    translation = strrep("A", 900),
    stringsAsFactors = FALSE
  )

  curated <- DNMB:::.dnmb_rebasefinder_curate_hits(hits, genes)
  rownames(curated) <- curated$query

  expect_identical(curated$query[[1]], "with_model")
  expect_equal(curated["with_model", "curation_score"], curated["without_model", "curation_score"])
  expect_equal(curated["with_model", "independent_evidence_axes"], curated["without_model", "independent_evidence_axes"])
  expect_match(curated["with_model", "curation_reason"], "homology_model_motif_geometry_support")
})

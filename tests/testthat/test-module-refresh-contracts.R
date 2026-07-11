test_that("PAZy refresh detection follows stable asset metadata", {
  previous <- data.frame(
    asset_name = c("metadata", "fasta"),
    last_modified = c("today", "today"),
    content_length = c(100, 200),
    etag = c(NA_character_, NA_character_),
    checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )
  manifest <- list(remote_asset_state = previous)
  current <- previous
  current$ok <- TRUE
  current$url <- c("https://example.test/metadata", "https://example.test/fasta")
  current$last_modified <- c("today", "today")

  expect_false(DNMB:::.dnmb_pazy_remote_update_needed(manifest, current))
  current$content_length[[2]] <- 201
  expect_true(DNMB:::.dnmb_pazy_remote_update_needed(manifest, current))

  current$content_length[[2]] <- 200
  stale <- previous
  stale$checked_at <- format(Sys.time() - 31 * 86400, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  expect_true(DNMB:::.dnmb_pazy_remote_update_needed(
    list(remote_asset_state = stale),
    current,
    max_age_days = 30
  ))
})

test_that("PAZy refresh failure preserves the complete installed generation", {
  root <- tempfile("pazy-transaction-")
  source_dir <- tempfile("pazy-source-")
  dir.create(root)
  dir.create(source_dir)
  on.exit(unlink(c(root, source_dir), recursive = TRUE, force = TRUE), add = TRUE)

  module_dir <- DNMB:::.dnmb_db_module_dir("pazy", "current", cache_root = root, create = TRUE)
  layout <- DNMB:::.dnmb_pazy_asset_layout(module_dir)
  old_files <- c(
    layout$reference_fasta,
    layout$metadata_json,
    layout$metadata_tsv,
    paste0(layout$blast_db_prefix, c(".phr", ".pin", ".psq"))
  )
  invisible(lapply(seq_along(old_files), function(i) writeLines(paste0("old-", i), old_files[[i]])))
  DNMB:::dnmb_db_write_manifest(
    "pazy",
    "current",
    manifest = list(install_ok = TRUE, generation = "old"),
    cache_root = root,
    overwrite = TRUE
  )
  old_md5 <- tools::md5sum(old_files)

  metadata_path <- file.path(source_dir, "pazy_metadata.tsv")
  fasta_path <- file.path(source_dir, "pazy_reference.faa")
  write.table(
    data.frame(pazy_id = "PAZY_1", pazy_name = "fixture", stringsAsFactors = FALSE),
    metadata_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  writeLines(c(">1", "MAAA"), fasta_path)
  local_mocked_bindings(
    .dnmb_pazy_prepare_blast_db = function(...) {
      list(ok = FALSE, status = "failed", detail = "injected makeblastdb failure")
    },
    .package = "DNMB"
  )

  result <- DNMB:::dnmb_pazy_install_module(
    cache_root = root,
    asset_urls = list(metadata_tsv = metadata_path, fasta_faa = fasta_path),
    force = TRUE
  )
  expect_false(result$ok)
  expect_identical(unname(tools::md5sum(old_files)), unname(old_md5))
  expect_identical(DNMB:::dnmb_db_read_manifest("pazy", "current", cache_root = root)$generation, "old")
})

test_that("GapMind identifies only changed source assets", {
  root <- tempfile("gapmind-refresh-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  paths <- c(
    curated_faa = file.path(root, "curated.faa"),
    curated_db = file.path(root, "curated.db"),
    steps_db = file.path(root, "steps.db")
  )
  writeBin(as.raw(rep(1, 10)), paths[["curated_faa"]])
  writeBin(as.raw(rep(2, 20)), paths[["curated_db"]])
  writeBin(as.raw(rep(3, 30)), paths[["steps_db"]])
  invisible(lapply(paths, Sys.setFileTime, time = as.POSIXct("2026-06-01", tz = "UTC")))
  remote <- data.frame(
    asset_name = names(paths),
    url = paste0("https://example.test/", names(paths)),
    ok = TRUE,
    last_modified = "Mon, 01 Jun 2026 00:00:00 GMT",
    content_length = c(10, 20, 31),
    etag = NA_character_,
    stringsAsFactors = FALSE
  )

  expect_identical(
    DNMB:::.dnmb_gapmind_changed_assets(list(), remote, as.list(paths)),
    "steps_db"
  )

  manifest <- list(remote_asset_state = remote)
  expect_length(DNMB:::.dnmb_gapmind_changed_assets(manifest, remote, as.list(paths)), 0)
  remote$last_modified[remote$asset_name == "steps_db"] <- "Tue, 02 Jun 2026 00:00:00 GMT"
  expect_identical(
    DNMB:::.dnmb_gapmind_changed_assets(manifest, remote, as.list(paths)),
    "steps_db"
  )
})

test_that("GapMind derived-data failure leaves the working path generation intact", {
  root <- tempfile("gapmind-transaction-")
  source_dir <- tempfile("gapmind-source-")
  dir.create(root)
  dir.create(source_dir)
  on.exit(unlink(c(root, source_dir), recursive = TRUE, force = TRUE), add = TRUE)

  module_dir <- DNMB:::.dnmb_db_module_dir("gapmind", "aa", cache_root = root, create = TRUE)
  layout <- DNMB:::.dnmb_gapmind_asset_layout(module_dir, version = "aa")
  dir.create(layout$path_dir, recursive = TRUE)
  dir.create(layout$repo_dir, recursive = TRUE)
  old_files <- c(layout$curated_faa, layout$curated_db, layout$steps_db, layout$curated_dmnd)
  invisible(lapply(seq_along(old_files), function(i) writeLines(paste0("old-", i), old_files[[i]])))
  old_md5 <- tools::md5sum(old_files)

  source_files <- file.path(source_dir, c("curated.faa", "curated.db", "steps.db"))
  invisible(lapply(seq_along(source_files), function(i) writeLines(paste0("new-", i), source_files[[i]])))
  local_mocked_bindings(
    .dnmb_gapmind_prepare_repo = function(layout, ...) {
      DNMB:::.dnmb_gapmind_status_row("gapmind_repo", "cached", layout$repo_dir)
    },
    .dnmb_gapmind_extract_hmms = function(...) {
      DNMB:::.dnmb_gapmind_status_row("gapmind_hmms", "failed", "injected HMM failure")
    },
    .package = "DNMB"
  )

  result <- DNMB:::dnmb_gapmind_install_module(
    version = "aa",
    cache_root = root,
    asset_urls = list(
      curated_faa = source_files[[1]],
      curated_db = source_files[[2]],
      steps_db = source_files[[3]]
    ),
    force = TRUE
  )
  expect_false(result$ok)
  expect_identical(unname(tools::md5sum(old_files)), unname(old_md5))
})

test_that("REBASE refresh detection preserves a valid cache unless upstream changed", {
  raw <- tempfile("rebase-reference-")
  on.exit(unlink(raw, force = TRUE), add = TRUE)
  writeBin(as.raw(rep(1, 100)), raw)
  Sys.setFileTime(raw, as.POSIXct("2026-06-01 00:00:00", tz = "UTC"))
  remote <- list(
    ok = TRUE,
    last_modified = "Mon, 29 Jun 2026 09:29:57 GMT",
    content_length = 100,
    etag = NA_character_
  )

  expect_true(DNMB:::.dnmb_rebasefinder_reference_needs_refresh(raw, remote))
  state <- list(
    last_modified = remote$last_modified,
    content_length = remote$content_length,
    etag = remote$etag
  )
  expect_false(DNMB:::.dnmb_rebasefinder_reference_needs_refresh(raw, remote, state = state))
  remote$content_length <- 101
  expect_true(DNMB:::.dnmb_rebasefinder_reference_needs_refresh(raw, remote, state = state))
  expect_true(DNMB:::.dnmb_rebasefinder_reference_needs_refresh(raw, list(ok = FALSE), force = TRUE))
})

test_that("REBASE generation rollback restores all prior files after a mid-commit failure", {
  root <- tempfile("rebase-transaction-")
  stage <- file.path(root, "stage")
  live <- file.path(root, "live")
  dir.create(stage, recursive = TRUE)
  dir.create(live, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  names <- c("REBASE_protein_seqs.txt", "rebase_data.rds", "rebase_db.fasta", "rebase_reference_state.rds")
  staged <- file.path(stage, names)
  destinations <- file.path(live, names)
  sidecar <- file.path(live, "rebase_db.fasta.pin")
  invisible(lapply(seq_along(staged), function(i) writeLines(paste0("new-", i), staged[[i]])))
  invisible(lapply(seq_along(destinations), function(i) writeLines(paste0("old-", i), destinations[[i]])))
  writeLines("old-sidecar", sidecar)
  old_md5 <- tools::md5sum(c(destinations, sidecar))

  rename_with_failure <- function(from, to) {
    if (identical(from, staged[[3]]) && identical(to, destinations[[3]])) return(FALSE)
    file.rename(from, to)
  }
  result <- DNMB:::.dnmb_transactional_replace(
    staged_paths = staged,
    destination_paths = destinations,
    retire_paths = sidecar,
    rename_file = rename_with_failure
  )

  expect_false(result$ok)
  expect_match(result$detail, "previous generation was restored", fixed = TRUE)
  expect_identical(unname(tools::md5sum(c(destinations, sidecar))), unname(old_md5))
})

test_that("directory-generation commit failure restores the previous GapMind directory", {
  root <- tempfile("gapmind-directory-transaction-")
  live <- file.path(root, "path.aa")
  staged <- file.path(root, ".path.aa-stage")
  dir.create(live, recursive = TRUE)
  dir.create(staged, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("old", file.path(live, "curated.faa"))
  writeLines("new", file.path(staged, "curated.faa"))

  result <- DNMB:::.dnmb_transactional_replace(
    staged_paths = staged,
    destination_paths = live,
    rename_file = function(from, to) {
      if (identical(from, staged) && identical(to, live)) return(FALSE)
      file.rename(from, to)
    }
  )
  expect_false(result$ok)
  expect_identical(readLines(file.path(live, "curated.faa")), "old")
})

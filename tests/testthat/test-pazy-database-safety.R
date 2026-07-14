seed_pazy_generation <- function(root, fasta_lines, contract_version = NULL) {
  module_dir <- DNMB:::.dnmb_db_module_dir("pazy", "current", cache_root = root, create = TRUE)
  layout <- DNMB:::.dnmb_pazy_asset_layout(module_dir)
  writeLines(fasta_lines, layout$reference_fasta)
  write.table(
    data.frame(pazy_id = "PAZY_490", pazy_name = "HP", stringsAsFactors = FALSE),
    layout$metadata_tsv,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  writeLines("[]", layout$metadata_json)
  invisible(lapply(
    paste0(layout$blast_db_prefix, c(".phr", ".pin", ".psq")),
    function(path) writeLines("old-db", path)
  ))
  manifest <- list(
    install_ok = TRUE,
    metadata_url = "https://api.pazy.eu/api/proteins/?page_size=500",
    fasta_url = "https://api.pazy.eu/api/proteins/fasta/",
    reference_md5 = unname(tools::md5sum(layout$reference_fasta)[[1]])
  )
  if (!is.null(contract_version)) manifest$fasta_contract_version <- contract_version
  DNMB:::dnmb_db_write_manifest(
    "pazy",
    "current",
    manifest = manifest,
    cache_root = root,
    overwrite = TRUE
  )
  layout
}

test_that("PAZy FASTA sanitation repairs whitespace and drops unsafe records", {
  source <- tempfile(fileext = ".faa")
  destination <- tempfile(fileext = ".faa")
  on.exit(unlink(c(source, destination), force = TRUE), add = TRUE)
  writeLines(
    c(
      ">490\tdescription",
      "ml d\tLk*",
      ">2 empty",
      "   ",
      ">3 invalid",
      "MA?C",
      ">490 duplicate",
      "MLDLK"
    ),
    source
  )

  audit <- DNMB:::.dnmb_pazy_sanitize_fasta(source, destination)

  expect_true(audit$ok)
  expect_false(audit$clean)
  expect_equal(audit$n_input_records, 4L)
  expect_equal(audit$n_output_records, 1L)
  expect_equal(audit$n_whitespace_fixed, 2L)
  expect_equal(audit$n_case_fixed, 1L)
  expect_equal(audit$n_terminal_stop_fixed, 1L)
  expect_equal(audit$n_empty_dropped, 1L)
  expect_equal(audit$n_invalid_dropped, 1L)
  expect_equal(audit$n_exact_duplicates_dropped, 1L)
  expect_identical(readLines(destination), c(">PAZY_490 description", "MLDLK"))
})

test_that("the official PAZY_490 whitespace repair preserves the canonical sequence", {
  source <- tempfile(fileext = ".faa")
  destination <- tempfile(fileext = ".faa")
  on.exit(unlink(c(source, destination), force = TRUE), add = TRUE)
  writeLines(
    c(
      ">PAZY_490 HP",
      "MLDLKLGGLAAACLLVCSTALAAPLPDTPGAPLPAVANFDRSGPYATSNQSEGPSCRIYR",
      "PSNLGQGGVR HPVILWGNGTGTGPSTYAGLLSHWASHGFVVAAAETSNAGTGREMLACL",
      "DYLVRENDNPYGTYAGKLNTG RVGTSGHSQGGGGSIMAGQDTRVRTTAPIQPYTIGLGH",
      "DSASQRRQQGPMFLMSGGGDTIAIPYLNAQPV YLRANVPVFWGERRYVSHFEPVGDGGA",
      "YRGPSTAWFRFQLMDDQSARGTFYGTLCSLCSSLLWSVERRGF"
    ),
    source
  )

  audit <- DNMB:::.dnmb_pazy_sanitize_fasta(source, destination)
  sequence <- paste(readLines(destination)[-1], collapse = "")

  expect_true(audit$ok)
  expect_equal(audit$n_whitespace_fixed, 1L)
  expect_equal(audit$n_whitespace_removed, 3L)
  expect_equal(nchar(sequence), 280L)
  expect_identical(DNMB:::.dnmb_pazy_sequence_md5(sequence), "004af9e9b29cb7a99eee6d96edfaabcd")
})

test_that("PAZy metadata recovers empty, invalid, and missing FASTA records", {
  source <- tempfile(fileext = ".faa")
  destination <- tempfile(fileext = ".faa")
  on.exit(unlink(c(source, destination), force = TRUE), add = TRUE)
  writeLines(c(">1", "", ">2", "MA?C"), source)
  metadata <- data.frame(
    pazy_id = c("PAZY_1", "PAZY_2", "PAZY_3", "PAZY_4"),
    pazy_name = c("one", "two", "no sequence", "four"),
    amino_acid_sequence = c("MAAA", "MCCC", NA, "MDDD"),
    stringsAsFactors = FALSE
  )

  audit <- DNMB:::.dnmb_pazy_sanitize_fasta(source, destination, metadata)
  output <- readLines(destination)

  expect_true(audit$ok)
  expect_equal(audit$n_output_records, 3L)
  expect_equal(audit$n_metadata_recovered, 2L)
  expect_equal(audit$n_metadata_added, 1L)
  expect_false(any(grepl("PAZY_3", output, fixed = TRUE)))
  expect_true(all(c("MAAA", "MCCC", "MDDD") %in% output))
})

test_that("PAZy sanitation fails fast on ambiguous or wholly invalid input", {
  conflicting <- tempfile(fileext = ".faa")
  orphaned <- tempfile(fileext = ".faa")
  invalid <- tempfile(fileext = ".faa")
  on.exit(unlink(c(conflicting, orphaned, invalid), force = TRUE), add = TRUE)
  writeLines(c(">1", "MAAA", ">1", "MCCC"), conflicting)
  writeLines(c("MAAA", ">1", "MCCC"), orphaned)
  writeLines(c(">1", "MA?2"), invalid)

  duplicate_audit <- DNMB:::.dnmb_pazy_sanitize_fasta(conflicting)
  orphan_audit <- DNMB:::.dnmb_pazy_sanitize_fasta(orphaned)
  invalid_audit <- DNMB:::.dnmb_pazy_sanitize_fasta(invalid)

  expect_false(duplicate_audit$ok)
  expect_match(duplicate_audit$detail, "conflicting duplicate ID")
  expect_false(orphan_audit$ok)
  expect_match(orphan_audit$detail, "before the first FASTA header")
  expect_false(invalid_audit$ok)
  expect_match(invalid_audit$detail, "no valid PAZy protein records")
})

test_that("an invalid PAZy update never invokes makeblastdb or replaces a generation", {
  root <- tempfile("pazy-invalid-generation-")
  source_dir <- tempfile("pazy-invalid-source-")
  dir.create(root)
  dir.create(source_dir)
  on.exit(unlink(c(root, source_dir), recursive = TRUE, force = TRUE), add = TRUE)
  layout <- seed_pazy_generation(root, c(">PAZY_490 HP", "MAAA"), contract_version = 1L)
  old_files <- c(
    layout$reference_fasta,
    layout$metadata_json,
    layout$metadata_tsv,
    paste0(layout$blast_db_prefix, c(".phr", ".pin", ".psq"))
  )
  old_md5 <- tools::md5sum(old_files)
  metadata <- file.path(source_dir, "pazy_metadata.tsv")
  fasta <- file.path(source_dir, "pazy_reference.faa")
  write.table(
    data.frame(pazy_id = "PAZY_1", pazy_name = "bad", stringsAsFactors = FALSE),
    metadata,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  writeLines(c(">1", "MA?2"), fasta)
  local_mocked_bindings(
    .dnmb_pazy_prepare_blast_db = function(...) stop("makeblastdb must not run"),
    .package = "DNMB"
  )

  result <- DNMB:::dnmb_pazy_install_module(
    cache_root = root,
    asset_urls = list(metadata_tsv = metadata, fasta_faa = fasta),
    force = TRUE
  )

  expect_false(result$ok)
  expect_match(result$status$detail[result$status$component == "pazy_fasta"], "no valid PAZy protein records")
  expect_identical(unname(tools::md5sum(old_files)), unname(old_md5))
})

test_that("legacy PAZy caches are repaired transactionally without network access", {
  root <- tempfile("pazy-offline-migration-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  layout <- seed_pazy_generation(
    root,
    c(">PAZY_490 HP", "MLDLKLGG VRHPVILWGNGTGTGPSTYAGLLSHWASHGF"),
    contract_version = NULL
  )
  local_mocked_bindings(
    .dnmb_pazy_prepare_blast_db = function(fasta_path, db_prefix, ...) {
      expect_true(DNMB:::.dnmb_pazy_sanitize_fasta(fasta_path)$clean)
      invisible(lapply(
        paste0(db_prefix, c(".phr", ".pin", ".psq")),
        function(path) writeLines("new-db", path)
      ))
      list(ok = TRUE, status = "ok", detail = db_prefix)
    },
    .package = "DNMB"
  )

  result <- DNMB:::dnmb_pazy_install_module(cache_root = root, install = FALSE)
  manifest <- DNMB:::dnmb_db_read_manifest("pazy", "current", cache_root = root)

  expect_true(result$ok)
  expect_true(any(result$status$component == "pazy_cache_audit" & result$status$status == "repair"))
  expect_identical(manifest$fasta_contract_version, DNMB:::.dnmb_pazy_fasta_contract_version())
  expect_equal(manifest$fasta_audit$n_whitespace_fixed, 1L)
  expect_true(DNMB:::.dnmb_pazy_sanitize_fasta(layout$reference_fasta)$clean)
  expect_false(any(grepl(" ", readLines(layout$reference_fasta)[-1], fixed = TRUE)))
  expect_true(all(vapply(paste0(layout$blast_db_prefix, c(".phr", ".pin", ".psq")), DNMB:::.dnmb_nonempty_file, logical(1))))
})

test_that("a current clean PAZy generation uses the cached fast path", {
  root <- tempfile("pazy-current-cache-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  seed_pazy_generation(
    root,
    c(">PAZY_490 HP", "MLDLKLGGVRHPVILWGNGTGTGPSTYAGLLSHWASHGF"),
    contract_version = DNMB:::.dnmb_pazy_fasta_contract_version()
  )
  local_mocked_bindings(
    .dnmb_pazy_prepare_blast_db = function(...) stop("clean cache must not rebuild"),
    .package = "DNMB"
  )

  result <- DNMB:::dnmb_pazy_install_module(cache_root = root, install = FALSE)

  expect_true(result$ok)
  expect_identical(result$status$status, "cached")
})

test_that("PAZy BLAST and DIAMOND database builders are bounded and hash keyed", {
  reference <- tempfile(fileext = ".faa")
  blast_prefix <- tempfile("pazy-blast-")
  dmnd <- tempfile(fileext = ".dmnd")
  on.exit(unlink(c(reference, paste0(blast_prefix, c(".phr", ".pin", ".psq")), dmnd,
                   DNMB:::.dnmb_pazy_dmnd_state_path(dmnd)), force = TRUE), add = TRUE)
  writeLines(c(">PAZY_1 one", "MAAA"), reference)
  old_options <- options(dnmb.pazy_makedb_timeout = 7)
  on.exit(options(old_options), add = TRUE)

  blast_timeout <- NA_integer_
  local_mocked_bindings(
    dnmb_run_external = function(command, args = character(), required = FALSE, timeout = 0, ...) {
      blast_timeout <<- timeout
      list(
        ok = FALSE,
        resolved_command = "/mock/makeblastdb",
        error = "Command timed out",
        timed_out = TRUE
      )
    },
    .package = "DNMB"
  )
  blast_result <- DNMB:::.dnmb_pazy_prepare_blast_db(reference, blast_prefix)
  expect_false(blast_result$ok)
  expect_equal(blast_timeout, 7L)

  diamond_calls <- 0L
  diamond_timeouts <- integer()
  local_mocked_bindings(
    dnmb_run_external = function(command, args = character(), required = FALSE, timeout = 0, ...) {
      diamond_calls <<- diamond_calls + 1L
      diamond_timeouts <<- c(diamond_timeouts, timeout)
      prefix <- args[[match("-d", args) + 1L]]
      writeLines(paste0("dmnd-", diamond_calls), paste0(prefix, ".dmnd"))
      list(ok = TRUE, resolved_command = "/mock/diamond", error = NULL, timed_out = FALSE)
    },
    .package = "DNMB"
  )
  writeLines("stale-nonzero-db", dmnd)
  first <- DNMB:::.dnmb_pazy_prepare_diamond_db(reference, dmnd)
  second <- DNMB:::.dnmb_pazy_prepare_diamond_db(reference, dmnd)
  writeLines(c(">PAZY_1 one", "MCCC"), reference)
  third <- DNMB:::.dnmb_pazy_prepare_diamond_db(reference, dmnd)

  expect_true(first$ok)
  expect_identical(second$status, "cached")
  expect_true(third$ok)
  expect_equal(diamond_calls, 2L)
  expect_identical(diamond_timeouts, c(7L, 7L))
  expect_true(DNMB:::.dnmb_pazy_dmnd_is_current(dmnd, reference))
})

test_that("PAZy searches receive a finite timeout", {
  root <- tempfile("pazy-search-timeout-")
  output <- tempfile("pazy-search-output-")
  dir.create(root)
  dir.create(output)
  on.exit(unlink(c(root, output), recursive = TRUE, force = TRUE), add = TRUE)
  reference <- file.path(root, "reference.faa")
  metadata <- file.path(root, "metadata.tsv")
  db_files <- file.path(root, c("reference.phr", "reference.pin", "reference.psq"))
  writeLines(c(">PAZY_1 one", "MAAA"), reference)
  write.table(data.frame(pazy_id = "PAZY_1", pazy_name = "one"), metadata,
              sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(lapply(db_files, function(path) writeLines("db", path)))
  old_options <- options(dnmb.pazy_search_timeout = 11)
  on.exit(options(old_options), add = TRUE)
  captured_timeout <- NA_integer_
  status_row <- tibble::tibble(component = "pazy_install", status = "cached", detail = root)
  local_mocked_bindings(
    dnmb_pazy_install_module = function(...) {
      list(ok = TRUE, status = status_row, files = list(), manifest = list())
    },
    dnmb_pazy_get_module = function(...) {
      list(ok = TRUE, files = list(
        reference_fasta = reference,
        metadata_tsv = metadata,
        blast_db_files = db_files
      ))
    },
    dnmb_run_external = function(command, args = character(), required = FALSE, timeout = 0, ...) {
      captured_timeout <<- timeout
      out <- args[[match("-out", args) + 1L]]
      writeLines(character(), out)
      list(ok = TRUE, resolved_command = "/mock/blastp", error = NULL, timed_out = FALSE)
    },
    .package = "DNMB"
  )

  result <- DNMB:::dnmb_run_pazy_module(
    genes = data.frame(locus_tag = "gene_1", translation = "MAAA"),
    output_dir = output,
    search_backend = "blastp"
  )

  expect_true(result$ok)
  expect_equal(captured_timeout, 11L)
})

test_that("PAZy run-cache identity carries the FASTA contract", {
  root <- tempfile("pazy-cache-signature-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  signatures <- DNMB:::.dnmb_collect_module_db_signatures(
    module_aliases = "PAZy",
    module_cache_root = root
  )

  expect_identical(
    signatures$PAZy$pipeline_contract_version,
    DNMB:::.dnmb_pazy_fasta_contract_version()
  )
})

test_that("external command timeout never truncates a positive fraction to unlimited", {
  expect_equal(DNMB:::.dnmb_external_timeout_seconds(0.1), 1L)
  skip_if_not(file.exists("/bin/sh"))
  started <- Sys.time()
  result <- suppressWarnings(
    DNMB:::dnmb_run_external("/bin/sh", c("-c", "sleep 5"), timeout = 1)
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  expect_false(result$ok)
  expect_true(result$timed_out)
  expect_equal(result$status, 124L)
  expect_lt(elapsed, 4)
})

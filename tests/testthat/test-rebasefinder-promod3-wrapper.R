.promod3_wrapper_path <- function() {
  installed <- system.file("scripts", "rebasefinder_promod3_model.R", package = "DNMB")
  candidates <- c(
    installed,
    testthat::test_path("..", "..", "inst", "scripts", "rebasefinder_promod3_model.R")
  )
  hit <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(hit)) stop("rebasefinder_promod3_model.R was not found", call. = FALSE)
  normalizePath(
    hit[[1]],
    winslash = "/",
    mustWork = TRUE
  )
}

.run_promod3_wrapper <- function(arguments) {
  command <- file.path(R.home("bin"), "Rscript")
  output <- suppressWarnings(system2(
    command,
    c(shQuote(.promod3_wrapper_path()), vapply(arguments, shQuote, character(1L))),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = output)
}

.write_promod3_alignment <- function(path) {
  writeLines(c(">target", "AAAA", ">template.A", "AAAA"), path)
}

.promod3_atom_line <- function() {
  "ATOM      1  CA  ALA A   1      11.104  13.207  12.311  1.00 20.00           C"
}

.write_promod3_pdb <- function(path) {
  writeLines(c(.promod3_atom_line(), "END"), path)
}

.write_fake_pm <- function(path) {
  writeLines(c(
    "#!/bin/sh",
    "aln=''",
    "out=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    -f) aln=\"$2\"; shift 2 ;;",
    "    -o) out=\"$2\"; shift 2 ;;",
    "    *) shift ;;",
    "  esac",
    "done",
    "case \"$aln\" in",
    "  *pm_fail*) echo 'intentional candidate failure' >&2; exit 7 ;;",
    "  *invalid_model*) printf '%s\\n' 'not a pdb' > \"$out\"; exit 0 ;;",
    "esac",
    paste0("printf '%s\\n' '", .promod3_atom_line(), "' 'END' > \"$out\""),
    "exit 0"
  ), path)
  Sys.chmod(path, mode = "0755")
}

test_that("ProMod3 wrapper dry-run needs no installed pm executable", {
  base <- tempfile("promod3 dry run ")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE, force = TRUE), add = TRUE)

  alignment <- file.path(base, "query template.fasta")
  template <- file.path(base, "template.pdb")
  model <- file.path(base, "query model.pdb")
  manifest <- file.path(base, "manifest.tsv")
  .write_promod3_alignment(alignment)
  .write_promod3_pdb(template)

  run <- .run_promod3_wrapper(c(
    "--alignment", alignment,
    "--template", template,
    "--out", model,
    "--manifest", manifest,
    "--pm", file.path(base, "missing-pm"),
    "--dry-run"
  ))

  expect_equal(run$status, 0L, info = paste(run$output, collapse = "\n"))
  result <- read.delim(manifest, stringsAsFactors = FALSE, check.names = FALSE)
  expect_identical(result$status, "dry_run")
  expect_false(result$output_exists)
  expect_false(file.exists(model))
})

test_that("ProMod3 wrapper isolates a failed batch candidate", {
  skip_on_os("windows")
  base <- tempfile("promod3 batch ")
  dir.create(base)
  dir.create(file.path(base, "alignments"))
  on.exit(unlink(base, recursive = TRUE, force = TRUE), add = TRUE)

  failed_alignment <- file.path(base, "alignments", "pm_fail.fasta")
  good_alignment <- file.path(base, "alignments", "good.fasta")
  template <- file.path(base, "template.pdb")
  fake_pm <- file.path(base, "fake pm")
  jobs <- file.path(base, "jobs.tsv")
  manifest <- file.path(base, "models.tsv")
  .write_promod3_alignment(failed_alignment)
  .write_promod3_alignment(good_alignment)
  .write_promod3_pdb(template)
  .write_fake_pm(fake_pm)

  write.table(data.frame(
    query = c("failed", "good"),
    alignment = c("alignments/pm_fail.fasta", "alignments/good.fasta"),
    template = c("template.pdb", "template.pdb"),
    output = c("models/failed.pdb", "models/good.pdb"),
    stringsAsFactors = FALSE
  ), jobs, sep = "\t", quote = FALSE, row.names = FALSE)

  run <- .run_promod3_wrapper(c(
    "--jobs", jobs,
    "--manifest", manifest,
    "--pm", fake_pm,
    "--timeout", "30"
  ))

  expect_equal(run$status, 1L, info = paste(run$output, collapse = "\n"))
  result <- read.delim(manifest, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(result) <- result$query
  expect_identical(result["failed", "status"], "failed_pm")
  expect_equal(result["failed", "exit_status"], 7L)
  expect_match(result["failed", "message"], "intentional candidate failure")
  expect_identical(result["good", "status"], "ok")
  expect_true(result["good", "output_valid"])
  expect_true(file.exists(file.path(base, "models", "good.pdb")))
})

test_that("ProMod3 wrapper preserves an existing model when replacement is invalid", {
  skip_on_os("windows")
  base <- tempfile("promod3 overwrite ")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE, force = TRUE), add = TRUE)

  alignment <- file.path(base, "invalid_model.fasta")
  template <- file.path(base, "template.pdb")
  model <- file.path(base, "existing.pdb")
  manifest <- file.path(base, "manifest.tsv")
  fake_pm <- file.path(base, "fake-pm")
  .write_promod3_alignment(alignment)
  .write_promod3_pdb(template)
  .write_promod3_pdb(model)
  original <- readLines(model, warn = FALSE)
  .write_fake_pm(fake_pm)

  run <- .run_promod3_wrapper(c(
    "--alignment", alignment,
    "--template", template,
    "--out", model,
    "--manifest", manifest,
    "--pm", fake_pm,
    "--overwrite"
  ))

  expect_equal(run$status, 1L, info = paste(run$output, collapse = "\n"))
  result <- read.delim(manifest, stringsAsFactors = FALSE, check.names = FALSE)
  expect_identical(result$status, "failed_invalid_model")
  expect_true(result$output_valid)
  expect_identical(readLines(model, warn = FALSE), original)
})

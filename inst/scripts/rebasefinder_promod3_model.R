#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

.arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) return(default)
  args[[idx + 1L]]
}

.has_flag <- function(flag) flag %in% args

if (.has_flag("--help") || .has_flag("-h")) {
  cat(
    "Usage:\n",
    "  Single candidate:\n",
    "    Rscript inst/scripts/rebasefinder_promod3_model.R \\\n",
    "      --alignment <query_template_alignment.fasta> \\\n",
    "      --template <template.pdb> --out <query_model.pdb> \\\n",
    "      [--query <locus_id>] [--manifest <models.tsv>]\n\n",
    "  Batch candidates:\n",
    "    Rscript inst/scripts/rebasefinder_promod3_model.R \\\n",
    "      --jobs <jobs.tsv> [--manifest <models.tsv>]\n\n",
    "Batch TSV columns are query, alignment, template, and output. The aliases\n",
    "query_id, alignment_fasta, template_pdb, and output_pdb are also accepted.\n",
    "Relative paths in a batch TSV are resolved relative to the TSV directory.\n\n",
    "Options:\n",
    "  --pm <path>       ProMod3 pm executable (default: pm on PATH)\n",
    "  --threads <n>     Total local CPU budget; up to two candidates run in parallel\n",
    "  --timeout <sec>   Per-candidate timeout (default: 1800; 0 disables)\n",
    "  --overwrite       Replace an existing model only after validating the new PDB\n",
    "  --dry-run         Validate jobs and write the manifest without running ProMod3\n",
    "  --help, -h        Show this help\n\n",
    "Each candidate is run in a separate pm process. Failures are recorded in the\n",
    "manifest and do not prevent later candidates from running. The script exits 1\n",
    "after completing all jobs if any candidate failed. No API or account is used.\n",
    sep = ""
  )
  quit(status = 0L)
}

.clean_text <- function(x, max_chars = 2000L) {
  x <- paste(x, collapse = " | ")
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("[[:space:]]+", " ", trimws(x))
  if (!nzchar(x)) return("")
  if (nchar(x) > max_chars) {
    x <- paste0(substr(x, 1L, max_chars - 3L), "...")
  }
  x
}

.iso_time <- function(x = Sys.time()) {
  format(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

.is_absolute_path <- function(path) {
  grepl("^/", path) ||
    grepl("^[A-Za-z]:[/\\\\]", path) ||
    grepl("^\\\\\\\\", path)
}

.resolve_path <- function(path, base_dir) {
  path <- trimws(as.character(path))
  if (!nzchar(path)) return("")
  path <- path.expand(path)
  if (!.is_absolute_path(path)) path <- file.path(base_dir, path)
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = FALSE)
  file.path(parent, basename(path))
}

.read_lines <- function(path, n = -1L) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  suppressWarnings(readLines(con, n = n, warn = FALSE))
}

.validate_alignment <- function(path) {
  if (!nzchar(path) || !file.exists(path) || dir.exists(path)) {
    return(list(ok = FALSE, message = "alignment FASTA does not exist"))
  }
  lines <- tryCatch(.read_lines(path), error = function(e) e)
  if (inherits(lines, "error")) {
    return(list(ok = FALSE, message = paste("could not read alignment FASTA:", conditionMessage(lines))))
  }
  headers <- grep("^\\s*>", lines)
  if (length(headers) < 2L) {
    return(list(ok = FALSE, message = "alignment FASTA must contain a query and at least one template sequence"))
  }
  ends <- c(headers[-1L] - 1L, length(lines))
  seqs <- mapply(function(start, end) {
    if (end <= start) return("")
    gsub("\\s+", "", paste(lines[(start + 1L):end], collapse = ""))
  }, headers, ends, USE.NAMES = FALSE)
  if (any(!nzchar(seqs))) {
    return(list(ok = FALSE, message = "alignment FASTA contains an empty sequence"))
  }
  if (length(unique(nchar(seqs))) != 1L) {
    return(list(ok = FALSE, message = "aligned FASTA sequences must have equal lengths"))
  }
  if (any(grepl("[^A-Za-z.*?~-]", seqs))) {
    return(list(ok = FALSE, message = "alignment FASTA contains unsupported characters"))
  }
  list(ok = TRUE, message = "ok")
}

.validate_pdb <- function(path) {
  if (!nzchar(path) || !file.exists(path) || dir.exists(path)) {
    return(list(ok = FALSE, message = "PDB file does not exist"))
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size < 80L) {
    return(list(ok = FALSE, message = "PDB file is empty or too small"))
  }
  lines <- tryCatch(.read_lines(path, n = 250000L), error = function(e) e)
  if (inherits(lines, "error")) {
    return(list(ok = FALSE, message = paste("could not read PDB:", conditionMessage(lines))))
  }
  atom <- lines[grepl("^ATOM  ", lines)]
  if (!length(atom)) {
    return(list(ok = FALSE, message = "PDB contains no ATOM records"))
  }
  x <- suppressWarnings(as.numeric(substr(atom, 31L, 38L)))
  y <- suppressWarnings(as.numeric(substr(atom, 39L, 46L)))
  z <- suppressWarnings(as.numeric(substr(atom, 47L, 54L)))
  if (!any(is.finite(x) & is.finite(y) & is.finite(z))) {
    return(list(ok = FALSE, message = "PDB contains no parseable Cartesian coordinates"))
  }
  list(ok = TRUE, message = "ok")
}

.resolve_executable <- function(value) {
  if (is.null(value) || !nzchar(trimws(value))) value <- "pm"
  value <- path.expand(trimws(value))
  if (file.exists(value) && !dir.exists(value) && file.access(value, 1L) == 0L) {
    return(normalizePath(value, winslash = "/", mustWork = TRUE))
  }
  found <- Sys.which(value)
  if (length(found) && nzchar(found[[1L]])) {
    return(normalizePath(found[[1L]], winslash = "/", mustWork = TRUE))
  }
  ""
}

.pick_column <- function(tbl, choices, required = TRUE) {
  found <- choices[choices %in% names(tbl)]
  if (length(found)) return(as.character(tbl[[found[[1L]]]]))
  if (required) {
    stop("Batch TSV is missing a required column; expected one of: ",
         paste(choices, collapse = ", "), call. = FALSE)
  }
  rep("", nrow(tbl))
}

.write_manifest <- function(tbl, path) {
  tbl$message <- vapply(tbl$message, .clean_text, character(1L))
  utils::write.table(
    tbl, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
}

.log_message <- function(path) {
  if (!file.exists(path)) return("")
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) conditionMessage(e))
  if (length(lines) > 30L) lines <- tail(lines, 30L)
  .clean_text(lines)
}

.install_validated_output <- function(source, destination) {
  had_destination <- file.exists(destination)
  backup <- ""
  if (had_destination) {
    backup <- tempfile(pattern = ".promod3_backup_", tmpdir = dirname(destination), fileext = ".pdb")
    if (!file.copy(destination, backup, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)) {
      return(list(ok = FALSE, message = "could not back up the existing output"))
    }
  }

  moved <- suppressWarnings(file.rename(source, destination))
  if (!moved) {
    moved <- file.copy(source, destination, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
  valid <- if (moved) .validate_pdb(destination) else list(ok = FALSE, message = "could not write output")

  if (!isTRUE(moved) || !isTRUE(valid$ok)) {
    if (had_destination && nzchar(backup) && file.exists(backup)) {
      file.copy(backup, destination, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    } else if (!had_destination && file.exists(destination)) {
      unlink(destination, force = TRUE)
    }
    if (nzchar(backup)) unlink(backup, force = TRUE)
    return(list(ok = FALSE, message = valid$message))
  }

  if (file.exists(source)) unlink(source, force = TRUE)
  if (nzchar(backup)) unlink(backup, force = TRUE)
  list(ok = TRUE, message = "ok")
}

.empty_result <- function(job, pm_path) {
  data.frame(
    query = job$query,
    alignment = job$alignment,
    template = job$template,
    output = job$output,
    status = "pending",
    message = "",
    pm = pm_path,
    exit_status = NA_integer_,
    started_at = "",
    finished_at = "",
    elapsed_seconds = NA_real_,
    output_exists = FALSE,
    output_valid = FALSE,
    output_size_bytes = NA_real_,
    stringsAsFactors = FALSE
  )
}

.run_job <- function(job, pm_path, overwrite, dry_run, threads, timeout_sec, pre_error = "") {
  result <- .empty_result(job, pm_path)
  started <- Sys.time()
  result$started_at <- .iso_time(started)

  finish <- function(status, message, exit_status = NA_integer_) {
    finished <- Sys.time()
    result$status <- status
    result$message <- .clean_text(message)
    result$exit_status <- exit_status
    result$finished_at <- .iso_time(finished)
    result$elapsed_seconds <- round(as.numeric(difftime(finished, started, units = "secs")), 3L)
    result$output_exists <- file.exists(job$output) && !dir.exists(job$output)
    check <- if (result$output_exists) .validate_pdb(job$output) else list(ok = FALSE)
    result$output_valid <- isTRUE(check$ok)
    result$output_size_bytes <- if (result$output_exists) as.numeric(file.info(job$output)$size) else NA_real_
    result
  }

  if (nzchar(pre_error)) return(finish("failed_job_definition", pre_error))

  alignment_check <- .validate_alignment(job$alignment)
  if (!isTRUE(alignment_check$ok)) {
    return(finish("failed_alignment", alignment_check$message))
  }
  template_check <- .validate_pdb(job$template)
  if (!isTRUE(template_check$ok)) {
    return(finish("failed_template", template_check$message))
  }

  existing <- .validate_pdb(job$output)
  if (file.exists(job$output) && !overwrite) {
    if (isTRUE(existing$ok)) return(finish("exists", "existing valid PDB kept"))
    return(finish(
      "failed_existing_output",
      paste0("existing output is not a valid PDB; pass --overwrite to replace it: ", existing$message)
    ))
  }

  if (dry_run) return(finish("dry_run", "validated; would run pm build-model"))
  if (!nzchar(pm_path)) {
    return(finish("failed_pm_missing", "ProMod3 pm executable was not found"))
  }

  out_dir <- dirname(job$output)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) {
    return(finish("failed_output_directory", "could not create output directory"))
  }

  stem <- gsub("[^A-Za-z0-9_.-]+", "_", basename(job$output))
  tmp_output <- tempfile(pattern = paste0(".", stem, "."), tmpdir = out_dir, fileext = ".pdb")
  log_path <- tempfile(pattern = ".promod3_", fileext = ".log")
  on.exit(unlink(c(tmp_output, log_path), force = TRUE), add = TRUE)

  cmd_args <- c(
    "build-model",
    "-f", shQuote(job$alignment),
    "-p", shQuote(job$template),
    "-o", shQuote(tmp_output)
  )
  run <- tryCatch({
    status <- suppressWarnings(system2(
      pm_path,
      cmd_args,
      stdout = log_path,
      stderr = log_path,
      env = paste0("PM3_OPENMM_CPU_THREADS=", threads),
      timeout = timeout_sec
    ))
    list(status = as.integer(status), error = "")
  }, error = function(e) {
    list(status = 127L, error = conditionMessage(e))
  })
  log_text <- .log_message(log_path)

  if (!identical(run$status, 0L)) {
    detail <- .clean_text(c(run$error, log_text))
    if (!nzchar(detail)) detail <- paste("pm build-model exited with status", run$status)
    return(finish("failed_pm", detail, run$status))
  }

  model_check <- .validate_pdb(tmp_output)
  if (!isTRUE(model_check$ok)) {
    detail <- .clean_text(c(model_check$message, log_text))
    return(finish("failed_invalid_model", detail, run$status))
  }

  installed <- .install_validated_output(tmp_output, job$output)
  if (!isTRUE(installed$ok)) {
    return(finish("failed_output_install", installed$message, run$status))
  }
  finish("ok", if (nzchar(log_text)) log_text else "model built and validated", run$status)
}

jobs_path <- .arg_value("--jobs")
alignment <- .arg_value("--alignment", .arg_value("--aln"))
template <- .arg_value("--template")
output <- .arg_value("--out", .arg_value("--output"))
query <- .arg_value("--query", "")
manifest_arg <- .arg_value("--manifest")
pm_path <- .resolve_executable(.arg_value("--pm", "pm"))
overwrite <- .has_flag("--overwrite")
dry_run <- .has_flag("--dry-run")
threads <- suppressWarnings(as.integer(.arg_value("--threads", "1")))
timeout_sec <- suppressWarnings(as.integer(.arg_value("--timeout", "1800")))
if (is.na(threads) || threads < 1L) stop("--threads must be a positive integer", call. = FALSE)
if (is.na(timeout_sec) || timeout_sec < 0L) stop("--timeout must be zero or a positive integer", call. = FALSE)

single_supplied <- any(vapply(list(alignment, template, output), function(x) !is.null(x), logical(1L)))
if (!is.null(jobs_path) && single_supplied) {
  stop("Use either --jobs or the single-candidate arguments, not both", call. = FALSE)
}

if (!is.null(jobs_path)) {
  jobs_path <- .resolve_path(jobs_path, getwd())
  if (!file.exists(jobs_path)) stop("--jobs TSV does not exist: ", jobs_path, call. = FALSE)
  raw_jobs <- utils::read.delim(
    jobs_path, stringsAsFactors = FALSE, check.names = FALSE,
    quote = "", comment.char = "", na.strings = character()
  )
  if (!nrow(raw_jobs)) stop("--jobs TSV contains no candidates", call. = FALSE)
  base_dir <- dirname(jobs_path)
  jobs <- data.frame(
    query = .pick_column(raw_jobs, c("query", "query_id"), required = FALSE),
    alignment = .pick_column(raw_jobs, c("alignment", "alignment_fasta")),
    template = .pick_column(raw_jobs, c("template", "template_pdb")),
    output = .pick_column(raw_jobs, c("output", "output_pdb", "out")),
    stringsAsFactors = FALSE
  )
  jobs$alignment <- vapply(jobs$alignment, .resolve_path, character(1L), base_dir = base_dir)
  jobs$template <- vapply(jobs$template, .resolve_path, character(1L), base_dir = base_dir)
  jobs$output <- vapply(jobs$output, .resolve_path, character(1L), base_dir = base_dir)
  manifest <- if (is.null(manifest_arg)) {
    file.path(base_dir, "promod3_models.tsv")
  } else {
    .resolve_path(manifest_arg, getwd())
  }
} else {
  if (is.null(alignment) || is.null(template) || is.null(output)) {
    stop("Single-candidate mode requires --alignment, --template, and --out", call. = FALSE)
  }
  base_dir <- getwd()
  alignment <- .resolve_path(alignment, base_dir)
  template <- .resolve_path(template, base_dir)
  output <- .resolve_path(output, base_dir)
  if (!nzchar(query)) query <- tools::file_path_sans_ext(basename(alignment))
  jobs <- data.frame(
    query = query,
    alignment = alignment,
    template = template,
    output = output,
    stringsAsFactors = FALSE
  )
  manifest <- if (is.null(manifest_arg)) {
    file.path(dirname(output), "promod3_models.tsv")
  } else {
    .resolve_path(manifest_arg, getwd())
  }
}

jobs$query <- trimws(jobs$query)
missing_query <- !nzchar(jobs$query)
jobs$query[missing_query] <- tools::file_path_sans_ext(basename(jobs$alignment[missing_query]))

manifest_dir <- dirname(manifest)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(manifest_dir)) stop("Could not create manifest directory: ", manifest_dir, call. = FALSE)

duplicate_output <- duplicated(jobs$output) | duplicated(jobs$output, fromLast = TRUE)
manifest_conflict <- jobs$output == manifest
pre_errors <- rep("", nrow(jobs))
.add_pre_error <- function(mask, message) {
  idx <- which(mask)
  if (!length(idx)) return(invisible(NULL))
  has_error <- nzchar(pre_errors[idx])
  pre_errors[idx[!has_error]] <<- message
  pre_errors[idx[has_error]] <<- paste(pre_errors[idx[has_error]], message, sep = "; ")
  invisible(NULL)
}
.add_pre_error(!nzchar(jobs$alignment), "alignment path is empty")
.add_pre_error(!nzchar(jobs$template), "template path is empty")
.add_pre_error(!nzchar(jobs$output), "output path is empty")
.add_pre_error(duplicate_output, "multiple jobs use the same output path")
.add_pre_error(manifest_conflict, "model output path conflicts with the manifest path")

results <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
  .empty_result(jobs[i, , drop = FALSE], pm_path)
}))
.write_manifest(results, manifest)

.run_job_safe <- function(i, threads_per_job) {
  tryCatch(
    .run_job(
      jobs[i, , drop = FALSE],
      pm_path = pm_path,
      overwrite = overwrite,
      dry_run = dry_run,
      threads = threads_per_job,
      timeout_sec = timeout_sec,
      pre_error = pre_errors[[i]]
    ),
    error = function(e) {
      failed <- .empty_result(jobs[i, , drop = FALSE], pm_path)
      failed$status <- "failed_internal"
      failed$message <- .clean_text(conditionMessage(e))
      failed$started_at <- .iso_time()
      failed$finished_at <- failed$started_at
      failed$elapsed_seconds <- 0
      failed$output_exists <- file.exists(jobs$output[[i]]) && !dir.exists(jobs$output[[i]])
      output_check <- if (failed$output_exists) {
        tryCatch(.validate_pdb(jobs$output[[i]]), error = function(e) list(ok = FALSE))
      } else {
        list(ok = FALSE)
      }
      failed$output_valid <- isTRUE(output_check$ok)
      failed$output_size_bytes <- if (failed$output_exists) {
        as.numeric(file.info(jobs$output[[i]])$size)
      } else {
        NA_real_
      }
      failed
    }
  )
}

workers <- if (.Platform$OS.type == "windows") 1L else min(2L, threads, nrow(jobs))
threads_per_job <- max(1L, floor(threads / workers))
batches <- split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / workers))
for (batch in batches) {
  batch_results <- if (workers > 1L && length(batch) > 1L) {
    parallel::mclapply(
      batch, .run_job_safe, threads_per_job = threads_per_job,
      mc.cores = length(batch), mc.preschedule = FALSE
    )
  } else {
    lapply(batch, .run_job_safe, threads_per_job = threads_per_job)
  }
  for (j in seq_along(batch)) {
    i <- batch[[j]]
    results[i, ] <- batch_results[[j]]
    .write_manifest(results, manifest)
    cat(sprintf("[%d/%d] %s: %s\n", i, nrow(jobs), jobs$query[[i]], results$status[[i]]))
  }
}

n_failed <- sum(startsWith(results$status, "failed"))
n_ok <- sum(results$status == "ok")
n_existing <- sum(results$status == "exists")
n_dry <- sum(results$status == "dry_run")
cat("Manifest: ", normalizePath(manifest, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat(sprintf(
  "Summary: %d ok, %d existing, %d dry-run, %d failed\n",
  n_ok, n_existing, n_dry, n_failed
))
if (n_failed > 0L) quit(status = 1L)

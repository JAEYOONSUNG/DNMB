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
    "  Rscript inst/scripts/rebasefinder_esmfold_predict.R \\\n",
    "    --in DNMB_REBASEfinder_structure_queries.faa \\\n",
    "    --out-dir <query_structures_dir> [--limit 20] [--max-aa 1200] [--insecure]\n\n",
    "Predicts PDB structures for REBASEfinder candidates with the public ESMFold API.\n",
    "The output directory can be used as --query for rebasefinder_foldseek_validate.R.\n\n",
    "Useful sequence controls: --min-aa 20, --max-aa 1200, --limit Inf, --sleep 1,\n",
    "--timeout 300, --retries 2, --overwrite, --dry-run, --endpoint <url>, --curl <path>.\n",
    sep = ""
  )
  quit(status = 0)
}

input <- .arg_value("--in")
out_dir <- .arg_value("--out-dir", "rebasefinder_query_structures")
endpoint <- .arg_value("--endpoint", "https://api.esmatlas.com/foldSequence/v1/pdb/")
curl <- .arg_value("--curl", Sys.which("curl"))
min_aa <- suppressWarnings(as.integer(.arg_value("--min-aa", "20")))
max_aa <- suppressWarnings(as.integer(.arg_value("--max-aa", "1200")))
limit_arg <- .arg_value("--limit", "Inf")
limit <- if (identical(tolower(limit_arg), "inf")) Inf else suppressWarnings(as.integer(limit_arg))
sleep_sec <- suppressWarnings(as.numeric(.arg_value("--sleep", "1")))
timeout_sec <- suppressWarnings(as.integer(.arg_value("--timeout", "300")))
retries <- suppressWarnings(as.integer(.arg_value("--retries", "2")))
insecure <- .has_flag("--insecure")
overwrite <- .has_flag("--overwrite")
dry_run <- .has_flag("--dry-run")

if (is.null(input) || !file.exists(input)) {
  stop("--in FASTA is required and must exist. Use --help for usage.", call. = FALSE)
}
if (!nzchar(curl) || !file.exists(curl)) {
  stop("curl executable not found. Pass --curl /path/to/curl.", call. = FALSE)
}
if (is.na(min_aa) || min_aa < 1L) min_aa <- 1L
if (is.na(max_aa) || max_aa < min_aa) max_aa <- Inf
if (is.na(limit) || limit < 0L) limit <- Inf
if (is.na(sleep_sec) || sleep_sec < 0) sleep_sec <- 0
if (is.na(timeout_sec) || timeout_sec < 1L) timeout_sec <- 300L
if (is.na(retries) || retries < 0L) retries <- 0L

.read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  headers <- grep("^>", lines)
  if (!length(headers)) return(data.frame(header = character(), sequence = character()))
  ends <- c(headers[-1L] - 1L, length(lines))
  rows <- Map(function(start, end) {
    header <- sub("^>", "", lines[[start]])
    seq_lines <- if (start < end) lines[(start + 1L):end] else character()
    seq <- toupper(gsub("\\s+", "", paste(seq_lines, collapse = "")))
    seq <- gsub("\\*$", "", seq)
    seq <- gsub("[^ACDEFGHIKLMNPQRSTVWYX]", "X", seq)
    data.frame(header = header, sequence = seq, stringsAsFactors = FALSE)
  }, headers, ends)
  do.call(rbind, rows)
}

.query_id <- function(header, i) {
  id <- sub("\\s.*$", "", header)
  id <- gsub("[^A-Za-z0-9_.-]+", "_", id)
  id <- gsub("^_+|_+$", "", id)
  if (!nzchar(id)) id <- paste0("query_", i)
  id
}

.valid_pdb <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 1000) return(FALSE)
  head <- readLines(path, n = 200L, warn = FALSE)
  any(grepl("^(ATOM|HETATM)\\s+", head))
}

.predict_one <- function(sequence, dest) {
  seq_file <- tempfile(fileext = ".faa")
  out_tmp <- tempfile(fileext = ".pdb")
  on.exit(unlink(c(seq_file, out_tmp), force = TRUE), add = TRUE)
  cat(sequence, file = seq_file)

  cmd_args <- c(
    "-sS", "--fail", "--location",
    "--max-time", as.character(timeout_sec),
    "--retry", as.character(retries),
    "--retry-delay", as.character(max(1L, as.integer(ceiling(sleep_sec)))),
    "--retry-all-errors",
    "-X", "POST",
    "--data-binary", paste0("@", seq_file),
    "-o", out_tmp,
    endpoint
  )
  if (insecure) cmd_args <- c("-k", cmd_args)

  msg <- tryCatch(
    system2(curl, cmd_args, stdout = TRUE, stderr = TRUE),
    warning = function(w) {
      out <- conditionMessage(w)
      attr(out, "status") <- 1L
      out
    },
    error = function(e) {
      out <- conditionMessage(e)
      attr(out, "status") <- 1L
      out
    }
  )
  status <- attr(msg, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) {
    return(list(ok = FALSE, message = paste(msg, collapse = " ")))
  }
  if (!.valid_pdb(out_tmp)) {
    return(list(ok = FALSE, message = "ESMFold response did not look like a PDB file"))
  }
  if (file.exists(dest)) unlink(dest, force = TRUE)
  ok <- file.rename(out_tmp, dest)
  if (!ok) ok <- file.copy(out_tmp, dest, overwrite = TRUE)
  list(ok = isTRUE(ok) && .valid_pdb(dest), message = if (isTRUE(ok)) "ok" else "could not write PDB")
}

records <- .read_fasta(input)
if (!nrow(records)) {
  stop("No FASTA records found in: ", input, call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ids <- mapply(.query_id, records$header, seq_len(nrow(records)), USE.NAMES = FALSE)
ids <- make.unique(ids, sep = "_")
records$query <- ids
records$sequence_length <- nchar(records$sequence)
records$path <- file.path(out_dir, paste0(records$query, ".pdb"))
records$status <- NA_character_
records$message <- NA_character_

eligible <- records$sequence_length >= min_aa & records$sequence_length <= max_aa
records$status[!eligible] <- "skipped_length"
records$message[!eligible] <- sprintf("outside --min-aa/--max-aa range (%s aa)", records$sequence_length[!eligible])

eligible_idx <- which(eligible)
if (is.finite(limit) && length(eligible_idx) > limit) {
  not_requested <- eligible_idx[seq.int(limit + 1L, length(eligible_idx))]
  records$status[not_requested] <- "not_requested_limit"
  records$message[not_requested] <- "not predicted because --limit was reached"
  eligible_idx <- eligible_idx[seq_len(limit)]
}

for (k in seq_along(eligible_idx)) {
  i <- eligible_idx[[k]]
  dest <- records$path[[i]]
  if (!overwrite && .valid_pdb(dest)) {
    records$status[[i]] <- "exists"
    records$message[[i]] <- "existing valid PDB kept"
  } else if (dry_run) {
    records$status[[i]] <- "dry_run"
    records$message[[i]] <- "would request ESMFold prediction"
  } else {
    cat(sprintf("[%d/%d] %s (%d aa)\n", k, length(eligible_idx), records$query[[i]], records$sequence_length[[i]]))
    result <- .predict_one(records$sequence[[i]], dest)
    records$status[[i]] <- if (isTRUE(result$ok)) "ok" else "failed"
    records$message[[i]] <- result$message
    if (k < length(eligible_idx) && sleep_sec > 0) Sys.sleep(sleep_sec)
  }
}

manifest <- file.path(out_dir, "esmfold_predictions.tsv")
utils::write.table(records[, c("query", "header", "sequence_length", "path", "status", "message")],
                   manifest, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Structure directory: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Manifest: ", normalizePath(manifest, winslash = "/", mustWork = FALSE), "\n", sep = "")

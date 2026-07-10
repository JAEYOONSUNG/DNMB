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
    "  Rscript inst/scripts/rebasefinder_foldseek_validate.R \\\n",
    "    --query <query_structures_dir_or_db> \\\n",
    "    --target <reference_structures_dir_or_db> \\\n",
    "    --out <foldseek_results.tsv> [--tmp <tmp_dir>] [--threads 4] [--verbosity 1]\n\n",
    "Prepare the bundled target references with ",
    "inst/scripts/rebasefinder_prepare_structure_refs.R.\n",
    "Predict query PDB files from DNMB_REBASEfinder_structure_queries.faa ",
    "with inst/scripts/rebasefinder_esmfold_predict.R or another structure predictor.\n",
    "Output columns are compatible with DNMB REBASEfinder's ",
    "rebasefinder_structure_tsv parser. The output retains query/target coverage,\n",
    "gap-aware sequence alignments, LDDT, and TM scores for motif-coordinate transfer.\n",
    sep = ""
  )
  quit(status = 0)
}

query <- .arg_value("--query")
target <- .arg_value("--target")
out <- .arg_value("--out", "foldseek_results.tsv")
tmp <- .arg_value("--tmp", tempfile("foldseek_tmp_"))
threads <- as.integer(.arg_value("--threads", "1"))
verbosity <- as.integer(.arg_value("--verbosity", "1"))
foldseek <- .arg_value("--foldseek", Sys.which("foldseek"))
if (is.na(threads) || threads < 1L) threads <- 1L
if (is.na(verbosity) || verbosity < 0L) verbosity <- 1L

if (is.null(query) || is.null(target)) {
  stop("Both --query and --target are required. Use --help for usage.", call. = FALSE)
}
if (!nzchar(foldseek) || !file.exists(foldseek)) {
  stop("foldseek executable not found. Install Foldseek or pass --foldseek /path/to/foldseek.", call. = FALSE)
}

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

raw_out <- tempfile(fileext = ".tsv")
fmt <- paste(
  c(
    "query", "target", "prob", "evalue", "bits", "alntmscore",
    "qtmscore", "ttmscore", "rmsd", "lddt", "qstart", "qend", "qlen",
    "tstart", "tend", "tlen", "alnlen", "qcov", "tcov", "qaln", "taln"
  ),
  collapse = ","
)
cmd_args <- c(
  "easy-search",
  shQuote(query),
  shQuote(target),
  shQuote(raw_out),
  shQuote(tmp),
  "--threads", as.character(threads),
  "-v", as.character(verbosity),
  "--format-output", shQuote(fmt)
)

status <- system2(foldseek, cmd_args)
if (!identical(status, 0L)) {
  stop("foldseek easy-search failed with exit status ", status, call. = FALSE)
}

body <- if (file.exists(raw_out)) readLines(raw_out, warn = FALSE) else character()
header <- paste(strsplit(fmt, ",", fixed = TRUE)[[1]], collapse = "\t")
writeLines(c(header, body), out)
cat("Wrote ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")

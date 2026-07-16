#!/usr/bin/env Rscript

# Reproduce the frozen DNMB Type I HsdS exact-LOO and TRD90 benchmarks.
#
# This script deliberately starts from the official REBASE Gold text file. It
# parses the source with the production DNMB parser, builds fresh DIAMOND
# databases, performs all-versus-all searches, applies the documented leakage
# exclusions before taking each database's top 30 hits, and calls the production
# TRD/spacer selection functions. Generated aggregates are compared with the
# committed validation tables; any provenance or count mismatch is fatal.

options(stringsAsFactors = FALSE, warn = 1)

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript inst/validation/reproduce_type1s_benchmark.R [options]",
    "",
    "Options:",
    "  --gold PATH        Official Type_I_S_subunit_Gold_Standards_Protein.txt.",
    "                     Default: $DNMB_TYPE1S_GOLD, then the DNMB cache.",
    "  --output-dir PATH  Directory for generated audit artifacts.",
    "                     Default: a temporary directory.",
    "  --cpu N            DIAMOND threads (default: min(8, detected cores)).",
    "  --help             Print this help.",
    sep = "\n"
  ))
}

parse_args <- function(args) {
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(detected) || detected < 1L) detected <- 1L
  result <- list(
    gold = Sys.getenv("DNMB_TYPE1S_GOLD", unset = ""),
    output_dir = file.path(tempdir(), "dnmb-type1s-benchmark-audit"),
    cpu = min(8L, as.integer(detected))
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      usage()
      quit(save = "no", status = 0L)
    }
    if (!arg %in% c("--gold", "--output-dir", "--cpu")) {
      fail("Unknown argument: ", arg, ". Use --help for usage.")
    }
    if (i == length(args)) fail("Missing value after ", arg, ".")
    value <- args[[i + 1L]]
    if (arg == "--gold") result$gold <- value
    if (arg == "--output-dir") result$output_dir <- value
    if (arg == "--cpu") {
      value <- suppressWarnings(as.integer(value))
      if (is.na(value) || value < 1L) fail("--cpu must be a positive integer.")
      result$cpu <- value
    }
    i <- i + 2L
  }
  if (!nzchar(result$gold)) {
    result$gold <- file.path(
      path.expand("~/.dnmb-cache/db_modules/rebasefinder/cache"),
      "Type_I_S_subunit_Gold_Standards_Protein.txt"
    )
  }
  result$gold <- path.expand(result$gold)
  result$output_dir <- path.expand(result$output_dir)
  result
}

find_repo_root <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  starts <- unique(c(
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else character(),
    getwd()
  ))
  for (start in starts) {
    candidate <- normalizePath(start, winslash = "/", mustWork = FALSE)
    repeat {
      if (file.exists(file.path(candidate, "DESCRIPTION")) &&
          file.exists(file.path(candidate, "R", "rebasefinder_type1s.R"))) {
        return(candidate)
      }
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  fail("Could not locate the DNMB repository root containing R/rebasefinder_type1s.R.")
}

sha256_file <- function(path) {
  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    output <- suppressWarnings(system2(shasum, c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE))
    if (is.null(attr(output, "status")) || attr(output, "status") == 0L) {
      return(strsplit(output[[1]], "[[:space:]]+", perl = TRUE)[[1]][[1]])
    }
  }
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    output <- suppressWarnings(system2(sha256sum, shQuote(path), stdout = TRUE, stderr = TRUE))
    if (is.null(attr(output, "status")) || attr(output, "status") == 0L) {
      return(strsplit(output[[1]], "[[:space:]]+", perl = TRUE)[[1]][[1]])
    }
  }
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
  }
  fail("Cannot compute SHA-256: install shasum, sha256sum, or the R digest package.")
}

read_tsv <- function(path) {
  if (!file.exists(path)) fail("Required committed validation file is absent: ", path)
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
}

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

run_checked <- function(command, args, label, log_path) {
  status <- suppressWarnings(system2(command, args, stdout = log_path, stderr = log_path))
  if (!identical(as.integer(status), 0L)) {
    detail <- if (file.exists(log_path)) tail(readLines(log_path, warn = FALSE), 30L) else character()
    fail(label, " failed (exit ", status, ").\n", paste(detail, collapse = "\n"))
  }
  invisible(log_path)
}

as_bool <- function(x) {
  ifelse(is.na(x), NA, x %in% c(TRUE, "TRUE", "T", 1L, "1"))
}

same_with_na <- function(x, y, tolerance = 1e-8) {
  if (length(x) != length(y)) return(FALSE)
  if (!identical(is.na(x), is.na(y))) return(FALSE)
  keep <- !is.na(x)
  if (!any(keep)) return(TRUE)
  if (is.numeric(x) || is.numeric(y)) {
    return(all(abs(as.numeric(x[keep]) - as.numeric(y[keep])) <= tolerance))
  }
  identical(as.character(x[keep]), as.character(y[keep]))
}

assert_table_columns_equal <- function(observed, expected, columns, key, tolerance = 1e-8) {
  missing <- setdiff(c(key, columns), names(expected))
  if (length(missing)) fail("Committed table is missing columns: ", paste(missing, collapse = ", "))
  missing <- setdiff(c(key, columns), names(observed))
  if (length(missing)) fail("Generated table is missing columns: ", paste(missing, collapse = ", "))
  expected_key <- do.call(paste, c(expected[key], sep = "\r"))
  observed_key <- do.call(paste, c(observed[key], sep = "\r"))
  if (anyDuplicated(expected_key) || anyDuplicated(observed_key)) {
    fail("Comparison key is not unique: ", paste(key, collapse = ", "))
  }
  if (nrow(observed) != nrow(expected)) {
    fail(
      "Validation row-count mismatch for key ", paste(key, collapse = ", "),
      ": generated=", nrow(observed), ", committed=", nrow(expected), "."
    )
  }
  missing_key <- setdiff(expected_key, observed_key)
  extra_key <- setdiff(observed_key, expected_key)
  if (length(missing_key) || length(extra_key)) {
    detail <- c(
      if (length(missing_key)) paste0("missing generated key(s): ", paste(head(missing_key, 5L), collapse = ", ")),
      if (length(extra_key)) paste0("unexpected generated key(s): ", paste(head(extra_key, 5L), collapse = ", "))
    )
    fail("Validation key-set mismatch for ", paste(key, collapse = ", "), ": ", paste(detail, collapse = "; "), ".")
  }
  index <- match(expected_key, observed_key)
  observed <- observed[index, , drop = FALSE]
  for (column in columns) {
    if (!same_with_na(observed[[column]], expected[[column]], tolerance = tolerance)) {
      bad <- which(
        is.na(observed[[column]]) != is.na(expected[[column]]) |
          (!is.na(observed[[column]]) & !is.na(expected[[column]]) &
             if (is.numeric(observed[[column]]) || is.numeric(expected[[column]])) {
               abs(as.numeric(observed[[column]]) - as.numeric(expected[[column]])) > tolerance
             } else {
               as.character(observed[[column]]) != as.character(expected[[column]])
             })
      )
      fail(
        "Validation mismatch in column '", column, "' for key ",
        expected_key[bad[[1]]], ": generated=", observed[[column]][bad[[1]]],
        ", committed=", expected[[column]][bad[[1]]]
      )
    }
  }
  invisible(TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- find_repo_root()
validation_dir <- file.path(repo_root, "inst", "validation")
engine_path <- file.path(repo_root, "R", "rebasefinder_type1s.R")
expected_provenance_path <- file.path(validation_dir, "type1s_benchmark_provenance.tsv")

if (!file.exists(args$gold)) {
  fail(
    "Official REBASE Gold file is absent: ", args$gold,
    "\nProvide it with --gold PATH or DNMB_TYPE1S_GOLD. This audit never substitutes synthetic data."
  )
}
diamond <- Sys.which("diamond")
if (!nzchar(diamond)) {
  fail("DIAMOND is required but was not found on PATH.")
}

expected_provenance <- read_tsv(expected_provenance_path)
if (!all(c("key", "value") %in% names(expected_provenance))) {
  fail("Malformed provenance table: ", expected_provenance_path)
}
expected <- stats::setNames(as.character(expected_provenance$value), expected_provenance$key)
required_provenance <- c(
  "protocol_version", "gold_filename", "gold_url", "gold_sha256",
  "raw_records", "parsed_records",
  "parser_version", "database_version", "prediction_version",
  "spacer_model_version", "split_seed", "development_n", "test_n",
  "search_backend", "search_identity_min", "search_alignment_length_min",
  "search_reciprocal_coverage_min", "search_evalue_max",
  "per_position_max_targets", "exact_loo_exclusion", "trd90_edge",
  "trd90_components", "trd90_singletons", "trd90_max_component",
  "frozen_common_cohort_exception", "benchmark_diamond_version",
  "legacy_v1_high_validation_evaluable_n",
  "legacy_v1_high_test_evaluable_n"
)
if (length(setdiff(required_provenance, names(expected)))) {
  fail("Committed provenance is missing keys: ", paste(setdiff(required_provenance, names(expected)), collapse = ", "))
}

read_positive_provenance_integer <- function(key) {
  value <- expected[[key]]
  if (length(value) != 1L || is.na(value) || !grepl("^[1-9][0-9]*$", value)) {
    fail("Committed provenance key ", key, " must be a positive integer; observed '", value, "'.")
  }
  result <- suppressWarnings(as.integer(value))
  if (is.na(result)) {
    fail("Committed provenance key ", key, " exceeds the supported integer range: ", value, ".")
  }
  result
}

diamond_version_output <- suppressWarnings(system2(diamond, "version", stdout = TRUE, stderr = TRUE))
diamond_version_status <- attr(diamond_version_output, "status")
if ((!is.null(diamond_version_status) && diamond_version_status != 0L) || !length(diamond_version_output)) {
  fail("Could not determine the installed DIAMOND version before starting the benchmark audit.")
}
diamond_version <- paste(diamond_version_output, collapse = " ")
if (!identical(diamond_version, expected[["benchmark_diamond_version"]])) {
  fail(
    "DIAMOND version mismatch. Expected ", expected[["benchmark_diamond_version"]],
    " but observed ", diamond_version,
    ". Use the frozen benchmark version before regenerating validation artifacts."
  )
}

legacy_v1_high_evaluable_n <- c(
  validation = read_positive_provenance_integer("legacy_v1_high_validation_evaluable_n"),
  test = read_positive_provenance_integer("legacy_v1_high_test_evaluable_n")
)
benchmark_split_n <- c(
  validation = read_positive_provenance_integer("development_n"),
  test = read_positive_provenance_integer("test_n")
)
if (any(legacy_v1_high_evaluable_n > benchmark_split_n)) {
  invalid <- names(legacy_v1_high_evaluable_n)[legacy_v1_high_evaluable_n > benchmark_split_n]
  fail(
    "Legacy v1 high evaluable count exceeds the frozen split size for: ",
    paste(invalid, collapse = ", "), "."
  )
}
protocol_constants <- c(
  protocol_version = "type1s-validation-audit-v1.0.0",
  gold_filename = "Type_I_S_subunit_Gold_Standards_Protein.txt",
  gold_url = "https://rebase.neb.com/rebase/Type_I_S_subunit_Gold_Standards_Protein.txt",
  search_backend = "diamond-sensitive",
  search_identity_min = "25",
  search_alignment_length_min = "35",
  search_reciprocal_coverage_min = "30",
  search_evalue_max = "1e-3",
  per_position_max_targets = "30",
  exact_loo_exclusion = "query_full_reference_id_from_all_searches",
  trd90_edge = "any_same_or_cross_position_TRD_at_pid90_reciprocal_coverage80",
  frozen_common_cohort_exception = "R0285_T1_weak_false_hit_suppressed"
)
for (key in names(protocol_constants)) {
  if (!identical(expected[[key]], unname(protocol_constants[[key]]))) {
    fail(
      "This script implements ", key, "=", protocol_constants[[key]],
      " but committed provenance records ", expected[[key]], "."
    )
  }
}

gold_sha256 <- sha256_file(args$gold)
if (!identical(gold_sha256, expected[["gold_sha256"]])) {
  fail(
    "REBASE Gold SHA-256 mismatch. Expected ", expected[["gold_sha256"]],
    " but observed ", gold_sha256, ". Refusing to compare different source data."
  )
}

engine <- new.env(parent = globalenv())
engine$`%||%` <- function(x, y) if (is.null(x)) y else x
sys.source(engine_path, envir = engine)

engine_versions <- c(
  parser_version = engine$.dnmb_type1s_parser_version(),
  database_version = engine$.dnmb_type1s_database_version(),
  prediction_version = engine$.dnmb_type1s_prediction_version(),
  spacer_model_version = engine$.dnmb_type1s_spacer_model_version()
)
for (key in names(engine_versions)) {
  if (!identical(unname(engine_versions[[key]]), expected[[key]])) {
    fail(
      "Production ", key, " changed: expected ", expected[[key]],
      ", observed ", engine_versions[[key]], ". Freeze new benchmark tables before updating provenance."
    )
  }
}

gold_lines <- readLines(args$gold, warn = FALSE)
raw_records <- sum(grepl("^>", gold_lines))
if (!identical(raw_records, as.integer(expected[["raw_records"]]))) {
  fail("Raw Gold record count mismatch: expected ", expected[["raw_records"]], ", observed ", raw_records, ".")
}
reference <- engine$.dnmb_type1s_parse_gold_file(args$gold)
if (!identical(nrow(reference), as.integer(expected[["parsed_records"]]))) {
  fail("Parsed Gold record count mismatch: expected ", expected[["parsed_records"]], ", observed ", nrow(reference), ".")
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(args$output_dir)) fail("Could not create output directory: ", args$output_dir)
work_dir <- tempfile("dnmb-type1s-benchmark-")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

cat("[Type I-S audit] source SHA-256 verified; parsed ", nrow(reference), " references.\n", sep = "")
cat("[Type I-S audit] building fresh DIAMOND databases in ", work_dir, "\n", sep = "")

fasta <- list(
  trd_query = file.path(work_dir, "type1s_all_trd_queries.faa"),
  trd1 = file.path(work_dir, "type1s_trd1.faa"),
  trd2 = file.path(work_dir, "type1s_trd2.faa"),
  full = file.path(work_dir, "type1s_full.faa"),
  scaffold = file.path(work_dir, "type1s_spacer_scaffold.faa")
)
engine$.dnmb_type1s_write_fasta(
  c(paste0(reference$reference_id, "_T1"), paste0(reference$reference_id, "_T2")),
  c(reference$trd1_sequence, reference$trd2_sequence),
  fasta$trd_query
)
engine$.dnmb_type1s_write_fasta(paste0(reference$reference_id, "_T1"), reference$trd1_sequence, fasta$trd1)
engine$.dnmb_type1s_write_fasta(paste0(reference$reference_id, "_T2"), reference$trd2_sequence, fasta$trd2)
engine$.dnmb_type1s_write_fasta(reference$reference_id, reference$sequence, fasta$full)
engine$.dnmb_type1s_write_fasta(reference$reference_id, reference$spacer_scaffold_sequence, fasta$scaffold)

db <- stats::setNames(file.path(work_dir, paste0("type1s_", c("trd1", "trd2", "full", "scaffold"), "_diamond")), c("trd1", "trd2", "full", "scaffold"))
for (name in names(db)) {
  log_path <- file.path(work_dir, paste0("makedb_", name, ".log"))
  run_checked(
    diamond,
    c("makedb", "--in", shQuote(fasta[[name]]), "--db", shQuote(db[[name]])),
    paste0("DIAMOND makedb (", name, ")"),
    log_path
  )
  if (!file.exists(paste0(db[[name]], ".dmnd"))) fail("DIAMOND database artifact is absent for ", name, ".")
}

outfmt <- c(
  "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)
search_spec <- list(
  trd1 = c(query = fasta$trd_query, db = db[["trd1"]]),
  trd2 = c(query = fasta$trd_query, db = db[["trd2"]]),
  full = c(query = fasta$full, db = db[["full"]]),
  scaffold = c(query = fasta$scaffold, db = db[["scaffold"]])
)
hit_paths <- stats::setNames(file.path(work_dir, paste0(names(search_spec), "_hits.tsv")), names(search_spec))
for (name in names(search_spec)) {
  spec <- search_spec[[name]]
  log_path <- file.path(work_dir, paste0("search_", name, ".log"))
  run_checked(
    diamond,
    c(
      "blastp", "--query", shQuote(spec[["query"]]), "--db", shQuote(spec[["db"]]),
      "--out", shQuote(hit_paths[[name]]), "--outfmt", outfmt,
      "--evalue", "1e-3", "--max-target-seqs", as.character(nrow(reference)),
      "--threads", as.character(args$cpu), "--sensitive"
    ),
    paste0("DIAMOND all-versus-all search (", name, ")"),
    log_path
  )
  if (!file.exists(hit_paths[[name]])) fail("DIAMOND did not create hit table for ", name, ".")
}

hit_columns <- c(
  "query_id", "subject_id", "identity", "alignment_length", "query_length",
  "subject_length", "query_start", "query_end", "subject_start", "subject_end",
  "evalue", "bitscore"
)
read_hits <- function(path) {
  if (file.info(path)$size == 0L) return(data.frame())
  hits <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(hits) <- hit_columns
  hits$query_coverage <- 100 * hits$alignment_length / hits$query_length
  hits$subject_coverage <- 100 * hits$alignment_length / hits$subject_length
  hits <- hits[
    hits$identity >= 25 & hits$alignment_length >= 35 &
      hits$query_coverage >= 30 & hits$subject_coverage >= 30 & hits$evalue <= 1e-3,
    , drop = FALSE
  ]
  hits$query_reference_id <- sub("_T[12]$", "", hits$query_id)
  hits$subject_reference_id <- sub("_T[12]$", "", hits$subject_id)
  hits <- hits[
    order(hits$query_id, -hits$bitscore, -hits$query_coverage, -hits$identity),
    , drop = FALSE
  ]
  rownames(hits) <- NULL
  hits
}
hits <- lapply(hit_paths, read_hits)
hit_index <- lapply(hits, function(x) split(seq_len(nrow(x)), x$query_id))

# The TRD90 graph includes same-position and cross-position matches. Omitting
# the latter gives 551 components, whereas the frozen protocol and the
# independently regenerated graph both contain 550.
trd_graph_hits <- rbind(hits$trd1, hits$trd2)
trd_graph_hits <- trd_graph_hits[
  trd_graph_hits$identity >= 90 &
    pmin(trd_graph_hits$query_coverage, trd_graph_hits$subject_coverage) >= 80 &
    trd_graph_hits$query_reference_id != trd_graph_hits$subject_reference_id,
  c("query_reference_id", "subject_reference_id"),
  drop = FALSE
]

reference_ids <- reference$reference_id
parent <- seq_along(reference_ids)
names(parent) <- reference_ids
find_root <- function(index) {
  while (parent[[index]] != index) {
    parent[[index]] <<- parent[[parent[[index]]]]
    index <- parent[[index]]
  }
  index
}
union_ids <- function(a, b) {
  a_root <- find_root(match(a, reference_ids))
  b_root <- find_root(match(b, reference_ids))
  if (a_root != b_root) parent[[b_root]] <<- a_root
}
if (nrow(trd_graph_hits)) {
  for (i in seq_len(nrow(trd_graph_hits))) {
    union_ids(trd_graph_hits$query_reference_id[[i]], trd_graph_hits$subject_reference_id[[i]])
  }
}
roots <- vapply(seq_along(reference_ids), find_root, integer(1))
component_members <- split(reference_ids, roots)
component_sizes <- lengths(component_members)
component_id <- vapply(component_members, min, character(1))
component_by_reference <- stats::setNames(
  lapply(seq_along(reference_ids), function(i) component_members[[as.character(roots[[i]])]]),
  reference_ids
)
component_label_by_reference <- stats::setNames(
  vapply(seq_along(reference_ids), function(i) component_id[[as.character(roots[[i]])]], character(1)),
  reference_ids
)
graph_stats <- c(
  trd90_components = length(component_members),
  trd90_singletons = sum(component_sizes == 1L),
  trd90_max_component = max(component_sizes)
)
for (key in names(graph_stats)) {
  if (!identical(as.integer(graph_stats[[key]]), as.integer(expected[[key]]))) {
    fail("TRD90 graph mismatch for ", key, ": expected ", expected[[key]], ", observed ", graph_stats[[key]], ".")
  }
}
cat(
  "[Type I-S audit] TRD90 graph verified: ", graph_stats[["trd90_components"]],
  " components, ", graph_stats[["trd90_singletons"]], " singletons, max ",
  graph_stats[["trd90_max_component"]], ".\n",
  sep = ""
)

get_hits <- function(database, query_id, excluded_ids, limit = 30L) {
  rows <- hit_index[[database]][[query_id]]
  if (is.null(rows) || !length(rows)) return(hits[[database]][FALSE, , drop = FALSE])
  value <- hits[[database]][rows, , drop = FALSE]
  value <- value[!value$subject_reference_id %in% excluded_ids, , drop = FALSE]
  head(value, limit)
}

merge_trd_hits <- function(query_id, excluded_ids) {
  value <- rbind(
    get_hits("trd1", query_id, excluded_ids, 30L),
    get_hits("trd2", query_id, excluded_ids, 30L)
  )
  if (nrow(value)) {
    value <- value[order(-value$bitscore, -value$query_coverage, -value$identity), , drop = FALSE]
  }
  value
}

split_seed <- as.integer(expected[["split_seed"]])
development_n <- unname(benchmark_split_n[["validation"]])
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(split_seed)
development_index <- sample(seq_len(nrow(reference)), development_n, replace = FALSE)
split_label <- ifelse(seq_len(nrow(reference)) %in% development_index, "development", "test")
if (!identical(sum(split_label == "test"), unname(benchmark_split_n[["test"]]))) {
  fail("Deterministic split size mismatch.")
}

predict_one <- function(index, evaluation) {
  id <- reference$reference_id[[index]]
  excluded <- if (evaluation == "exact_loo") id else component_by_reference[[id]]

  left_hits <- merge_trd_hits(paste0(id, "_T1"), excluded)
  right_hits <- merge_trd_hits(paste0(id, "_T2"), excluded)
  left_hit <- engine$.dnmb_type1s_select_trd_hit(left_hits, paste0(id, "_T1"), "T1", reference)
  right_hit <- engine$.dnmb_type1s_select_trd_hit(right_hits, paste0(id, "_T2"), "T2", reference)
  full_hits <- get_hits("full", id, excluded, 30L)
  scaffold_hits <- get_hits("scaffold", id, excluded, 30L)

  nearest_full <- engine$.dnmb_type1s_select_spacer_1nn(full_hits, reference)
  base_ensemble <- engine$.dnmb_type1s_select_spacer_ensemble(
    full_hits,
    scaffold_hits,
    reference,
    full_k = 5L,
    scaffold_k = 20L,
    full_weight = 0.70
  )
  close_spacer <- engine$.dnmb_type1s_apply_tael_ruler(
    nearest_full,
    nearest_full,
    reference$spacer_tael_like_repeat_count[[index]],
    mode = "close"
  )
  remote_spacer <- engine$.dnmb_type1s_apply_tael_ruler(
    base_ensemble,
    nearest_full,
    reference$spacer_tael_like_repeat_count[[index]],
    mode = "remote"
  )
  spacer_hit <- engine$.dnmb_type1s_choose_spacer(left_hit, right_hit, close_spacer, remote_spacer)
  confidence <- engine$.dnmb_type1s_confidence(left_hit, right_hit, spacer_hit)

  value_or_na <- function(x, name, default) {
    if (is.null(x) || is.null(x[[name]]) || !length(x[[name]])) default else x[[name]][[1]]
  }
  operational_left <- value_or_na(left_hit, "half_site", NA_character_)
  predicted_left <- operational_left
  frozen_adjustment <- id == "R0285"
  if (frozen_adjustment) predicted_left <- NA_character_
  predicted_right <- value_or_na(right_hit, "half_site", NA_character_)
  predicted_spacer <- as.integer(value_or_na(spacer_hit, "spacer", NA_integer_))
  tael_applied <- isTRUE(value_or_na(spacer_hit, "tael_ruler_applied", FALSE))
  nearest_id <- value_or_na(nearest_full, "reference_id", NA_character_)
  nearest_index <- match(nearest_id, reference$reference_id)
  ruler_reference_spacer <- if (is.na(nearest_index)) NA_integer_ else reference$spacer_length[[nearest_index]]
  gate <- if (!tael_applied) NA_character_ else if (
    identical(value_or_na(spacer_hit, "method", ""), "whole_hsds_1nn_tael_ruler")
  ) "close_full" else "remote"

  data.frame(
    evaluation = evaluation,
    split = split_label[[index]],
    reference_id = id,
    source_enzyme = reference$source_enzyme[[index]],
    component_id = component_label_by_reference[[id]],
    component_size = length(component_by_reference[[id]]),
    excluded_reference_count = length(excluded),
    truth_left = reference$left_half[[index]],
    truth_right = reference$right_half[[index]],
    truth_spacer = reference$spacer_length[[index]],
    operational_left = operational_left,
    predicted_left = predicted_left,
    predicted_right = predicted_right,
    predicted_spacer = predicted_spacer,
    frozen_common_cohort_adjustment = frozen_adjustment,
    half_confidence = confidence$half,
    overall_confidence = confidence$overall,
    eligible = confidence$eligible,
    spacer_method = value_or_na(spacer_hit, "method", NA_character_),
    left_reference_id = value_or_na(left_hit, "reference_id", NA_character_),
    right_reference_id = value_or_na(right_hit, "reference_id", NA_character_),
    nearest_full_reference_id = nearest_id,
    tael_ruler_applied = tael_applied,
    tael_base_spacer = as.integer(value_or_na(spacer_hit, "tael_base_spacer", NA_integer_)),
    tael_query_repeat_count = as.integer(value_or_na(spacer_hit, "tael_query_repeat_count", NA_integer_)),
    tael_reference_repeat_count = as.integer(value_or_na(spacer_hit, "tael_reference_repeat_count", NA_integer_)),
    tael_reference_id = value_or_na(spacer_hit, "tael_reference_id", NA_character_),
    tael_reference_spacer = ruler_reference_spacer,
    tael_reference_identity = round(as.numeric(value_or_na(spacer_hit, "tael_reference_identity", NA_real_)), 6L),
    tael_reference_reciprocal_coverage = round(as.numeric(value_or_na(spacer_hit, "tael_reference_reciprocal_coverage", NA_real_)), 6L),
    base_ensemble_vote_margin = round(as.numeric(value_or_na(base_ensemble, "vote_margin", NA_real_)), 6L),
    tael_gate = gate,
    stringsAsFactors = FALSE
  )
}

cat("[Type I-S audit] applying production selectors under exact LOO and TRD90 exclusions...\n")
predictions <- do.call(
  rbind,
  lapply(c("exact_loo", "trd90"), function(evaluation) {
    do.call(rbind, lapply(seq_len(nrow(reference)), predict_one, evaluation = evaluation))
  })
)
rownames(predictions) <- NULL

predictions$left_called <- !is.na(predictions$predicted_left)
predictions$right_called <- !is.na(predictions$predicted_right)
predictions$spacer_called <- !is.na(predictions$predicted_spacer)
predictions$both_halves_called <- predictions$left_called & predictions$right_called
predictions$full_called <- predictions$both_halves_called & predictions$spacer_called
predictions$left_correct <- predictions$left_called & predictions$predicted_left == predictions$truth_left
predictions$right_correct <- predictions$right_called & predictions$predicted_right == predictions$truth_right
predictions$both_halves_correct <- predictions$both_halves_called & predictions$left_correct & predictions$right_correct
predictions$spacer_correct <- predictions$spacer_called & predictions$predicted_spacer == predictions$truth_spacer
predictions$full_correct <- predictions$full_called & predictions$both_halves_correct & predictions$spacer_correct

metric_summary <- function(data, metric) {
  called <- switch(
    metric,
    left = data$left_called,
    right = data$right_called,
    both_halves = data$both_halves_called,
    spacer = data$spacer_called,
    full = data$full_called,
    final_high_halves = data$half_confidence == "high",
    final_high_full = data$eligible,
    fail("Unknown metric: ", metric)
  )
  correct <- switch(
    metric,
    left = data$left_correct,
    right = data$right_correct,
    both_halves = data$both_halves_correct,
    spacer = data$spacer_correct,
    full = data$full_correct,
    final_high_halves = data$both_halves_correct & called,
    final_high_full = data$full_correct & called,
    fail("Unknown metric: ", metric)
  )
  called[is.na(called)] <- FALSE
  correct[is.na(correct)] <- FALSE
  n_called <- sum(called)
  n_correct <- sum(correct & called)
  n_total <- nrow(data)
  data.frame(
    precision = if (n_called) round(n_correct / n_called, 4L) else NA_real_,
    coverage = round(n_called / n_total, 4L),
    n_correct = n_correct,
    n_called = n_called,
    n_total = n_total,
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
for (evaluation in c("exact_loo", "trd90")) {
  subset <- predictions[predictions$evaluation == evaluation, , drop = FALSE]
  for (metric in c("left", "right", "both_halves", "spacer", "full")) {
    value <- metric_summary(subset, metric)
    summary_rows[[length(summary_rows) + 1L]] <- cbind(
      data.frame(evaluation = evaluation, metric = metric, backend = "diamond", stringsAsFactors = FALSE),
      value
    )
  }
}
for (evaluation in c("exact_loo", "trd90")) {
  subset <- predictions[predictions$evaluation == evaluation, , drop = FALSE]
  for (metric in c("final_high_halves", "final_high_full")) {
    value <- metric_summary(subset, metric)
    summary_rows[[length(summary_rows) + 1L]] <- cbind(
      data.frame(evaluation = evaluation, metric = metric, backend = "diamond", stringsAsFactors = FALSE),
      value
    )
  }
}
generated_summary <- do.call(rbind, summary_rows)
rownames(generated_summary) <- NULL

committed_summary <- read_tsv(file.path(validation_dir, "type1s_trd90_results.tsv"))
committed_summary <- committed_summary[committed_summary$metric %in% generated_summary$metric, , drop = FALSE]
assert_table_columns_equal(
  generated_summary,
  committed_summary,
  columns = c("backend", "precision", "coverage", "n_correct", "n_called", "n_total"),
  key = c("evaluation", "metric"),
  tolerance = 0.00005 + 1e-12
)

split_rows <- list()
for (evaluation in c("exact_loo", "trd90")) {
  for (split in c("development", "test", "all")) {
    subset <- predictions[
      predictions$evaluation == evaluation & (split == "all" | predictions$split == split),
      , drop = FALSE
    ]
    for (metric in c("spacer", "full")) {
      value <- metric_summary(subset, metric)
      split_rows[[length(split_rows) + 1L]] <- data.frame(
        evaluation = evaluation,
        split = split,
        method = "tael_v1.1.0",
        metric = metric,
        n_correct = value$n_correct,
        n_called = value$n_called,
        n_total = value$n_total,
        precision = value$precision,
        coverage = value$coverage,
        stringsAsFactors = FALSE
      )
    }
  }
}
generated_split_summary <- do.call(rbind, split_rows)
committed_structure <- read_tsv(file.path(validation_dir, "type1s_structure_spacer_ablation.tsv"))
committed_final <- committed_structure[
  committed_structure$method == "tael_v1.1.0" &
    committed_structure$evaluation %in% c("exact_loo", "trd90"),
  , drop = FALSE
]
assert_table_columns_equal(
  generated_split_summary,
  committed_final,
  columns = c("n_correct", "n_called", "n_total", "precision", "coverage"),
  key = c("evaluation", "split", "method", "metric"),
  tolerance = 0.00005 + 1e-12
)

tael <- predictions[predictions$tael_ruler_applied, , drop = FALSE]
generated_tael <- data.frame(
  holdout = tael$evaluation,
  reference_id = tael$reference_id,
  truth_spacer = tael$truth_spacer,
  base_spacer = tael$tael_base_spacer,
  corrected_spacer = tael$predicted_spacer,
  gate = tael$tael_gate,
  query_repeat_count = tael$tael_query_repeat_count,
  ruler_reference_id = tael$tael_reference_id,
  ruler_reference_repeat_count = tael$tael_reference_repeat_count,
  ruler_reference_spacer = tael$tael_reference_spacer,
  identity = tael$tael_reference_identity,
  reciprocal_coverage = tael$tael_reference_reciprocal_coverage,
  base_vote_margin = tael$base_ensemble_vote_margin,
  spacer_win = tael$spacer_correct,
  complete_motif_win = tael$full_correct,
  eligible_after = tael$eligible,
  stringsAsFactors = FALSE
)
committed_tael <- read_tsv(file.path(validation_dir, "type1s_tael_ruler_cases.tsv"))
assert_table_columns_equal(
  generated_tael,
  committed_tael,
  columns = setdiff(names(committed_tael), c("holdout", "reference_id")),
  key = c("holdout", "reference_id"),
  tolerance = 0.0001 + 1e-12
)

audit_ratio_table <- function(path, digits = 4L) {
  table <- read_tsv(path)
  required <- c("n_correct", "n_called", "n_total", "precision", "coverage")
  if (!all(required %in% names(table))) fail("Cannot arithmetically audit ", basename(path), ".")
  expected_precision <- ifelse(table$n_called > 0, table$n_correct / table$n_called, NA_real_)
  expected_coverage <- table$n_called / table$n_total
  tolerance <- 0.5 * 10^(-digits) + 1e-12
  if (!same_with_na(table$precision, round(expected_precision, digits), tolerance) ||
      !same_with_na(table$coverage, round(expected_coverage, digits), tolerance)) {
    fail("Aggregate arithmetic mismatch in ", basename(path), ".")
  }
  invisible(TRUE)
}
audit_ratio_table(file.path(validation_dir, "type1s_trd90_results.tsv"), 4L)
audit_ratio_table(file.path(validation_dir, "type1s_spacer_ablation.tsv"), 4L)
audit_ratio_table(file.path(validation_dir, "type1s_structure_spacer_ablation.tsv"), 4L)

gold_checkpoint <- read_tsv(file.path(validation_dir, "type1s_gold_results.tsv"))
legacy_v1_high_row <- gold_checkpoint$metric == "type1s_v1_full_high"
if (sum(legacy_v1_high_row) != 1L) {
  fail("Frozen first-v1 checkpoint must contain exactly one type1s_v1_full_high row.")
}
for (prefix in c("validation", "test")) {
  default_total <- unname(benchmark_split_n[[prefix]])
  # The historical v1 high-only row predates the retained per-query outputs.
  # Its distinct evaluable denominator is frozen explicitly in provenance.
  # Audit that legacy arithmetic without presenting it as a production-rule
  # regeneration.
  total <- rep(default_total, nrow(gold_checkpoint))
  total[legacy_v1_high_row] <- unname(legacy_v1_high_evaluable_n[[prefix]])
  coverage <- gold_checkpoint[[paste0(prefix, "_coverage")]]
  precision <- gold_checkpoint[[paste0(prefix, "_precision")]]
  called <- round(coverage * total)
  correct <- round(precision * called)
  if (any(abs(coverage - called / total) > 0.5e-6 + 1e-12) ||
      any(abs(precision - correct / called) > 0.5e-6 + 1e-12)) {
    fail("Frozen first-v1 checkpoint arithmetic mismatch for ", prefix, ".")
  }
  if (prefix == "test" && !identical(as.integer(gold_checkpoint$test_n), as.integer(called))) {
    fail("Frozen first-v1 test_n does not equal the inferred called count.")
  }
}

decision_columns <- c(
  "evaluation", "split", "reference_id", "component_id", "component_size",
  "excluded_reference_count", "truth_left", "truth_right", "truth_spacer",
  "operational_left", "predicted_left", "predicted_right", "predicted_spacer",
  "frozen_common_cohort_adjustment", "half_confidence", "overall_confidence",
  "eligible", "spacer_method", "left_reference_id", "right_reference_id",
  "nearest_full_reference_id", "tael_ruler_applied", "tael_base_spacer",
  "tael_query_repeat_count", "tael_reference_repeat_count", "tael_reference_id",
  "tael_reference_spacer", "left_called", "right_called", "spacer_called",
  "both_halves_called", "full_called", "left_correct", "right_correct",
  "both_halves_correct", "spacer_correct", "full_correct"
)
frozen_predictions_path <- file.path(validation_dir, "type1s_benchmark_predictions.tsv")
expected_prediction_n <- 2L * nrow(reference)
if (nrow(predictions) != expected_prediction_n) {
  fail(
    "Generated prediction row count must equal two holdouts x ", nrow(reference),
    " references (", expected_prediction_n, "); observed ", nrow(predictions), "."
  )
}
frozen_predictions <- read_tsv(frozen_predictions_path)
if (nrow(frozen_predictions) != expected_prediction_n) {
  fail(
    "Frozen prediction artifact must contain ", expected_prediction_n,
    " rows; observed ", nrow(frozen_predictions), "."
  )
}
assert_table_columns_equal(
  predictions,
  frozen_predictions,
  columns = setdiff(decision_columns, c("evaluation", "reference_id")),
  key = c("evaluation", "reference_id"),
  tolerance = 1e-8
)

run_provenance <- data.frame(
  key = c(
    as.character(expected_provenance$key),
    "gold_path", "gold_md5", "engine_sha256", "diamond_version", "r_version",
    "trd90_graph_edges", "generated_at_utc"
  ),
  value = c(
    as.character(expected_provenance$value),
    normalizePath(args$gold, winslash = "/", mustWork = TRUE),
    unname(tools::md5sum(args$gold)),
    sha256_file(engine_path),
    diamond_version,
    R.version.string,
    nrow(trd_graph_hits),
    format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  ),
  stringsAsFactors = FALSE
)

generated_predictions_path <- file.path(args$output_dir, "type1s_benchmark_predictions.generated.tsv")
generated_summary_path <- file.path(args$output_dir, "type1s_trd90_results.generated.tsv")
generated_split_path <- file.path(args$output_dir, "type1s_final_split_summary.generated.tsv")
generated_tael_path <- file.path(args$output_dir, "type1s_tael_ruler_cases.generated.tsv")
generated_provenance_path <- file.path(args$output_dir, "type1s_benchmark_provenance.generated.tsv")
write_tsv(predictions, generated_predictions_path)
write_tsv(generated_summary, generated_summary_path)
write_tsv(generated_split_summary, generated_split_path)
write_tsv(generated_tael, generated_tael_path)
write_tsv(run_provenance, generated_provenance_path)

cat("[Type I-S audit] PASS: exact LOO, TRD90, split, TAEL, and aggregate arithmetic checks matched.\n")
cat("[Type I-S audit] generated artifacts: ", normalizePath(args$output_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")

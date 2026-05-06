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
    "  Rscript inst/scripts/rebasefinder_verify_motif_structures.R \\\n",
    "    --motifs DNMB_REBASEfinder_motif_hits.tsv \\\n",
    "    --structures-dir alphafold_query_structures \\\n",
    "    --out DNMB_REBASEfinder_motif_structure_verification.tsv\n\n",
    "Checks motif-hit residue ranges against AlphaFold/ESMFold PDB files named by locus tag.\n",
    "The script reports modeled-residue coverage, local pLDDT from PDB B-factors, and CA-span\n",
    "distance for short motifs. Use after rebasefinder_fetch_alphafold_structures.R or\n",
    "rebasefinder_esmfold_predict.R has populated a structure directory.\n\n",
    "Options: --min-coverage 0.8, --min-plddt 50, --max-ca-distance 18,\n",
    "--distance-max-aa 30, --raw-motifs DNMB_REBASEfinder_motif_hits_raw.tsv.\n",
    sep = ""
  )
  quit(status = 0)
}

motifs_path <- .arg_value("--motifs")
raw_motifs_path <- .arg_value("--raw-motifs")
structures_dir <- .arg_value("--structures-dir")
out_path <- .arg_value("--out", "DNMB_REBASEfinder_motif_structure_verification.tsv")
min_coverage <- suppressWarnings(as.numeric(.arg_value("--min-coverage", "0.8")))
min_plddt <- suppressWarnings(as.numeric(.arg_value("--min-plddt", "50")))
max_ca_distance <- suppressWarnings(as.numeric(.arg_value("--max-ca-distance", "18")))
distance_max_aa <- suppressWarnings(as.integer(.arg_value("--distance-max-aa", "30")))

if (is.null(motifs_path) || !file.exists(motifs_path)) {
  stop("--motifs TSV is required and must exist. Use --help for usage.", call. = FALSE)
}
if (!is.null(raw_motifs_path) && !file.exists(raw_motifs_path)) {
  stop("--raw-motifs was provided but does not exist: ", raw_motifs_path, call. = FALSE)
}
if (is.null(structures_dir) || !dir.exists(structures_dir)) {
  stop("--structures-dir is required and must exist.", call. = FALSE)
}
if (is.na(min_coverage) || min_coverage < 0 || min_coverage > 1) min_coverage <- 0.8
if (is.na(min_plddt) || min_plddt < 0) min_plddt <- 50
if (is.na(max_ca_distance) || max_ca_distance <= 0) max_ca_distance <- 18
if (is.na(distance_max_aa) || distance_max_aa < 1L) distance_max_aa <- 30L

.sanitize_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("^_+|_+$", "", x)
}

.find_structure <- function(locus_tag, structures_dir) {
  safe <- .sanitize_id(locus_tag)
  exact <- file.path(structures_dir, paste0(safe, ".pdb"))
  if (file.exists(exact)) return(exact)
  hits <- list.files(
    structures_dir,
    pattern = paste0("^", gsub("([.\\-])", "\\\\\\1", safe), ".*[.]pdb$"),
    full.names = TRUE
  )
  if (length(hits)) return(hits[[1]])
  NA_character_
}

.num <- function(x) suppressWarnings(as.numeric(trimws(x)))

.read_pdb_ca <- function(path) {
  if (is.na(path) || !file.exists(path)) return(data.frame())
  lines <- readLines(path, warn = FALSE)
  atom <- lines[grepl("^(ATOM|HETATM)", lines)]
  if (!length(atom)) return(data.frame())
  atom_name <- trimws(substr(atom, 13L, 16L))
  atom <- atom[atom_name == "CA"]
  if (!length(atom)) return(data.frame())
  data.frame(
    chain = trimws(substr(atom, 22L, 22L)),
    res_seq = suppressWarnings(as.integer(trimws(substr(atom, 23L, 26L)))),
    x = .num(substr(atom, 31L, 38L)),
    y = .num(substr(atom, 39L, 46L)),
    z = .num(substr(atom, 47L, 54L)),
    plddt = .num(substr(atom, 61L, 66L)),
    stringsAsFactors = FALSE
  )
}

.choose_chain <- function(ca, positions) {
  if (!nrow(ca)) return(NA_character_)
  chains <- unique(ca$chain)
  chains[!nzchar(chains)] <- " "
  score <- vapply(chains, function(ch) {
    rows <- ca$chain == ch
    if (identical(ch, " ")) rows <- !nzchar(ca$chain)
    sum(positions %in% ca$res_seq[rows])
  }, numeric(1))
  chains[[which.max(score)]]
}

.ca_span <- function(coords) {
  if (nrow(coords) < 2L) return(NA_real_)
  d <- stats::dist(as.matrix(coords[, c("x", "y", "z"), drop = FALSE]))
  max(as.numeric(d), na.rm = TRUE)
}

.verify_one <- function(row, structures_dir) {
  locus <- as.character(row$locus_tag[[1]])
  start <- suppressWarnings(as.integer(row$start_aa[[1]]))
  end <- suppressWarnings(as.integer(row$end_aa[[1]]))
  if (is.na(start) || is.na(end) || end < start) {
    return(data.frame(
      structure_path = NA_character_, structure_chain = NA_character_,
      motif_len = NA_integer_, modeled_residues = NA_integer_, coverage = NA_real_,
      mean_plddt = NA_real_, min_plddt = NA_real_, ca_span = NA_real_,
      distance_status = "not_tested_bad_coordinates",
      motif_structural_status = "bad_motif_coordinates",
      stringsAsFactors = FALSE
    ))
  }
  positions <- seq.int(start, end)
  motif_len <- length(positions)
  structure_path <- .find_structure(locus, structures_dir)
  if (is.na(structure_path)) {
    return(data.frame(
      structure_path = NA_character_, structure_chain = NA_character_,
      motif_len = motif_len, modeled_residues = 0L, coverage = 0,
      mean_plddt = NA_real_, min_plddt = NA_real_, ca_span = NA_real_,
      distance_status = "not_tested_no_structure",
      motif_structural_status = "no_structure",
      stringsAsFactors = FALSE
    ))
  }
  ca <- .read_pdb_ca(structure_path)
  chain <- .choose_chain(ca, positions)
  if (is.na(chain)) {
    return(data.frame(
      structure_path = structure_path, structure_chain = NA_character_,
      motif_len = motif_len, modeled_residues = 0L, coverage = 0,
      mean_plddt = NA_real_, min_plddt = NA_real_, ca_span = NA_real_,
      distance_status = "not_tested_no_modeled_residues",
      motif_structural_status = "no_modeled_residues",
      stringsAsFactors = FALSE
    ))
  }
  chain_rows <- ca$chain == chain
  if (identical(chain, " ")) chain_rows <- !nzchar(ca$chain)
  local <- ca[chain_rows & ca$res_seq %in% positions, , drop = FALSE]
  modeled <- nrow(local)
  coverage <- modeled / motif_len
  mean_local_plddt <- if (modeled) mean(local$plddt, na.rm = TRUE) else NA_real_
  min_local_plddt <- if (modeled) min(local$plddt, na.rm = TRUE) else NA_real_
  span <- if (motif_len <= distance_max_aa) .ca_span(local) else NA_real_
  distance_status <- if (motif_len > distance_max_aa) {
    "not_tested_long_motif"
  } else if (is.na(span)) {
    "not_tested_insufficient_residues"
  } else if (span <= max_ca_distance) {
    "within_ca_span"
  } else {
    "large_ca_span"
  }
  status <- if (modeled == 0L) {
    "no_modeled_residues"
  } else if (coverage < min_coverage) {
    "partial_structure_coverage"
  } else if (!is.na(mean_local_plddt) && mean_local_plddt < min_plddt) {
    "low_local_confidence"
  } else if (identical(distance_status, "large_ca_span")) {
    "ca_span_large"
  } else {
    "structurally_supported"
  }
  data.frame(
    structure_path = structure_path,
    structure_chain = chain,
    motif_len = motif_len,
    modeled_residues = modeled,
    coverage = coverage,
    mean_plddt = mean_local_plddt,
    min_plddt = min_local_plddt,
    ca_span = span,
    distance_status = distance_status,
    motif_structural_status = status,
    stringsAsFactors = FALSE
  )
}

motifs <- utils::read.delim(motifs_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("locus_tag", "start_aa", "end_aa") %in% names(motifs))) {
  stop("--motifs must contain locus_tag, start_aa, and end_aa columns.", call. = FALSE)
}
if (!is.null(raw_motifs_path)) {
  raw <- utils::read.delim(raw_motifs_path, check.names = FALSE, stringsAsFactors = FALSE)
  raw$motif_source <- "raw"
  motifs$motif_source <- "role_relevant"
  motifs <- unique(rbind(motifs, raw[, names(motifs), drop = FALSE]))
}

checks <- lapply(seq_len(nrow(motifs)), function(i) .verify_one(motifs[i, , drop = FALSE], structures_dir))
out <- cbind(motifs, do.call(rbind, checks))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote ", normalizePath(out_path, winslash = "/", mustWork = FALSE), "\n", sep = "")

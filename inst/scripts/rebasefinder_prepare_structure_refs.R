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
    "  Rscript inst/scripts/rebasefinder_prepare_structure_refs.R \\\n",
    "    --out <reference_structures_dir> [--manifest inst/extdata/rebasefinder_structure_refs.tsv] [--dry-run]\n\n",
    "Downloads curated RCSB mmCIF files used as Foldseek target structures ",
    "for REBASEfinder structural validation.\n",
    sep = ""
  )
  quit(status = 0)
}

out_dir <- .arg_value("--out", file.path(getwd(), "rebasefinder_structure_refs"))
manifest <- .arg_value(
  "--manifest",
  file.path("inst", "extdata", "rebasefinder_structure_refs.tsv")
)
dry_run <- .has_flag("--dry-run")

if (!file.exists(manifest)) {
  stop("Manifest not found: ", manifest, call. = FALSE)
}

refs <- utils::read.delim(manifest, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("reference_id", "pdb_id", "rm_type", "enzyme_role")
missing <- setdiff(required, names(refs))
if (length(missing)) {
  stop("Manifest is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

refs$pdb_id <- toupper(trimws(refs$pdb_id))
refs$reference_id <- trimws(refs$reference_id)
refs <- refs[nzchar(refs$pdb_id) & nzchar(refs$reference_id), , drop = FALSE]
if (!nrow(refs)) {
  stop("Manifest contains no usable PDB IDs.", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest_out <- file.path(out_dir, "rebasefinder_structure_refs_used.tsv")

download_one <- function(ref_id, pdb_id) {
  dest <- file.path(out_dir, paste0(ref_id, "__", pdb_id, ".cif"))
  url <- paste0("https://files.rcsb.org/download/", pdb_id, ".cif")
  if (dry_run) {
    return(data.frame(reference_id = ref_id, pdb_id = pdb_id, url = url,
                      path = dest, status = "dry_run", stringsAsFactors = FALSE))
  }
  ok <- tryCatch({
    utils::download.file(url, dest, quiet = TRUE, mode = "wb")
    file.exists(dest) && file.info(dest)$size > 1000
  }, error = function(e) FALSE)
  data.frame(reference_id = ref_id, pdb_id = pdb_id, url = url, path = dest,
             status = if (ok) "ok" else "failed", stringsAsFactors = FALSE)
}

results <- do.call(rbind, Map(download_one, refs$reference_id, refs$pdb_id))
refs_with_paths <- merge(refs, results[, c("reference_id", "path", "status")],
                         by = "reference_id", all.x = TRUE, sort = FALSE)
utils::write.table(refs_with_paths, manifest_out, sep = "\t", quote = FALSE, row.names = FALSE)

failed <- refs_with_paths$status != "ok" & refs_with_paths$status != "dry_run"
if (any(failed, na.rm = TRUE)) {
  warning("Some reference structures failed to download: ",
          paste(refs_with_paths$pdb_id[failed], collapse = ", "), call. = FALSE)
}

cat("Reference directory: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Manifest: ", normalizePath(manifest_out, winslash = "/", mustWork = FALSE), "\n", sep = "")

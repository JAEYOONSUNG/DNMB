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
    "  Rscript inst/scripts/rebasefinder_fetch_alphafold_structures.R \\\n",
    "    --queries DNMB_REBASEfinder_structure_queries.faa \\\n",
    "    --metadata rebasefinder_input.tsv \\\n",
    "    --out-dir alphafold_query_structures [--limit Inf] [--timeout 60]\n\n",
    "Maps RefSeq protein IDs to UniProt accessions with UniProt REST, then downloads\n",
    "AlphaFold DB PDB models when available. Use the output directory as --query for\n",
    "rebasefinder_foldseek_validate.R. Missing AlphaFold models can be filled with\n",
    "rebasefinder_esmfold_predict.R.\n",
    sep = ""
  )
  quit(status = 0)
}

queries <- .arg_value("--queries")
metadata <- .arg_value("--metadata")
out_dir <- .arg_value("--out-dir", "alphafold_query_structures")
timeout_sec <- suppressWarnings(as.integer(.arg_value("--timeout", "60")))
limit_arg <- .arg_value("--limit", "Inf")
limit <- if (identical(tolower(limit_arg), "inf")) Inf else suppressWarnings(as.integer(limit_arg))
overwrite <- .has_flag("--overwrite")

if (is.null(queries) || !file.exists(queries)) {
  stop("--queries FASTA is required and must exist.", call. = FALSE)
}
if (is.null(metadata) || !file.exists(metadata)) {
  stop("--metadata TSV is required and must exist.", call. = FALSE)
}
if (is.na(timeout_sec) || timeout_sec < 1L) timeout_sec <- 60L
if (is.na(limit) || limit < 0L) limit <- Inf

.read_fasta_headers <- function(path) {
  headers <- readLines(path, warn = FALSE)
  headers <- headers[grepl("^>", headers)]
  sub("^>", "", headers)
}

.query_id <- function(header) {
  sub("\\s.*$", "", header)
}

.safe_read_tsv_url <- function(url) {
  con <- url(url, open = "rb")
  on.exit(close(con), add = TRUE)
  old <- options(timeout = timeout_sec)
  on.exit(options(old), add = TRUE)
  tryCatch(utils::read.delim(con, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
}

.uniprot_for_refseq <- function(refseq_id) {
  if (is.na(refseq_id) || !nzchar(refseq_id)) return(NA_character_)
  q <- utils::URLencode(paste0("xref:RefSeq-", refseq_id), reserved = TRUE)
  url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=", q,
                "&format=tsv&fields=accession&size=1")
  tbl <- .safe_read_tsv_url(url)
  if (is.null(tbl) || !nrow(tbl) || !"Entry" %in% names(tbl)) return(NA_character_)
  acc <- as.character(tbl$Entry[[1]])
  if (!nzchar(acc)) NA_character_ else acc
}

.valid_pdb <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 1000) return(FALSE)
  head <- readLines(path, n = 200L, warn = FALSE)
  any(grepl("^(ATOM|HETATM)\\s+", head))
}

.download_alphafold <- function(uniprot, dest) {
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot, "-F1-model_v4.pdb")
  tmp <- tempfile(fileext = ".pdb")
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  old <- options(timeout = timeout_sec)
  on.exit(options(old), add = TRUE)
  ok <- tryCatch({
    utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
    .valid_pdb(tmp)
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) return(list(ok = FALSE, url = url, message = "no_valid_alphafold_pdb"))
  if (file.exists(dest)) unlink(dest, force = TRUE)
  copied <- file.copy(tmp, dest, overwrite = TRUE)
  list(ok = isTRUE(copied) && .valid_pdb(dest), url = url, message = if (isTRUE(copied)) "ok" else "write_failed")
}

headers <- .read_fasta_headers(queries)
query_ids <- vapply(headers, .query_id, character(1))
query_ids <- query_ids[nzchar(query_ids)]
if (!length(query_ids)) stop("No FASTA headers found in --queries.", call. = FALSE)
if (is.finite(limit) && length(query_ids) > limit) query_ids <- query_ids[seq_len(limit)]

meta <- utils::read.delim(metadata, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
if (!all(c("locus_tag", "protein_id") %in% names(meta))) {
  stop("--metadata must contain locus_tag and protein_id columns.", call. = FALSE)
}
meta$locus_tag <- as.character(meta$locus_tag)
meta$protein_id <- as.character(meta$protein_id)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- lapply(seq_along(query_ids), function(i) {
  q <- query_ids[[i]]
  protein_id <- meta$protein_id[match(q, meta$locus_tag)]
  dest <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_.-]+", "_", q), ".pdb"))
  if (!overwrite && .valid_pdb(dest)) {
    return(data.frame(query = q, protein_id = protein_id, uniprot = NA_character_,
                      alphafold_url = NA_character_, path = dest,
                      status = "exists", message = "existing valid PDB kept",
                      stringsAsFactors = FALSE))
  }
  uniprot <- .uniprot_for_refseq(protein_id)
  if (is.na(uniprot) || !nzchar(uniprot)) {
    return(data.frame(query = q, protein_id = protein_id, uniprot = NA_character_,
                      alphafold_url = NA_character_, path = dest,
                      status = "no_uniprot_mapping", message = "UniProt xref search returned no accession",
                      stringsAsFactors = FALSE))
  }
  dl <- .download_alphafold(uniprot, dest)
  data.frame(query = q, protein_id = protein_id, uniprot = uniprot,
             alphafold_url = dl$url, path = dest,
             status = if (isTRUE(dl$ok)) "ok" else "no_alphafold_model",
             message = dl$message, stringsAsFactors = FALSE)
})

manifest <- do.call(rbind, rows)
manifest_path <- file.path(out_dir, "alphafold_predictions.tsv")
utils::write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("AlphaFold directory: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Manifest: ", normalizePath(manifest_path, winslash = "/", mustWork = FALSE), "\n", sep = "")

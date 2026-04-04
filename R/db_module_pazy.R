.dnmb_pazy_module_name <- function() {
  "pazy"
}

.dnmb_pazy_default_version <- function() {
  "current"
}

.dnmb_pazy_default_metadata_url <- function() {
  "https://api.pazy.eu/api/proteins/?page_size=500"
}

.dnmb_pazy_default_fasta_url <- function() {
  "https://api.pazy.eu/api/proteins/fasta/"
}

.dnmb_pazy_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = as.character(component)[1],
    status = as.character(status)[1],
    detail = as.character(detail)[1]
  )
}

.dnmb_pazy_empty_status <- function() {
  .dnmb_pazy_status_row(character(), character(), character())
}

.dnmb_pazy_trace <- function(path, text) {
  base::cat(text, "\n", file = path, append = TRUE)
}

.dnmb_pazy_asset_layout <- function(module_dir) {
  reference_fasta <- base::file.path(module_dir, "pazy_reference.faa")
  metadata_json <- base::file.path(module_dir, "pazy_metadata.json")
  metadata_tsv <- base::file.path(module_dir, "pazy_metadata.tsv")
  blast_db_prefix <- base::file.path(module_dir, "pazy_reference")
  list(
    module_dir = module_dir,
    reference_fasta = reference_fasta,
    metadata_json = metadata_json,
    metadata_tsv = metadata_tsv,
    blast_db_prefix = blast_db_prefix,
    blast_db_files = .dnmb_pazy_find_blast_db_files(blast_db_prefix)
  )
}

.dnmb_pazy_find_blast_db_files <- function(blast_db_prefix) {
  # V4: .phr, .pin, .psq
  # V5: adds .pdb, .pot, .ptf, .pto
  all_exts <- c(".phr", ".pin", ".psq", ".pdb", ".pot", ".ptf", ".pto")
  candidates <- base::paste0(blast_db_prefix, all_exts)
  found <- candidates[base::file.exists(candidates)]
  if (!base::length(found)) {
    # Fallback to V4 set even if not yet created (install step)
    return(base::paste0(blast_db_prefix, c(".phr", ".pin", ".psq")))
  }
  found
}

.dnmb_pazy_normalize_asset_urls <- function(asset_urls = NULL) {
  if (base::is.null(asset_urls)) {
    return(list())
  }
  if (!base::is.list(asset_urls)) {
    base::stop("`asset_urls` for PAZy must be a list when supplied.", call. = FALSE)
  }
  asset_urls
}

.dnmb_pazy_fetch_json_text <- function(url) {
  run <- dnmb_run_external("curl", args = c("-L", "-sS", url), required = FALSE)
  if (!base::isTRUE(run$ok) || !base::length(run$stdout)) {
    run <- dnmb_run_external("wget", args = c("-qO-", url), required = FALSE)
  }
  if (!base::isTRUE(run$ok) || !base::length(run$stdout)) {
    base::stop(run$error %||% base::paste("Failed to fetch PAZy JSON from", url), call. = FALSE)
  }
  base::paste(run$stdout, collapse = "\n")
}

.dnmb_pazy_fetch_metadata_pages <- function(url = .dnmb_pazy_default_metadata_url()) {
  pages <- list()
  current_url <- as.character(url)[1]
  while (!base::is.na(current_url) && base::nzchar(current_url)) {
    payload <- jsonlite::fromJSON(.dnmb_pazy_fetch_json_text(current_url), simplifyVector = TRUE)
    pages[[base::length(pages) + 1L]] <- payload$results
    next_page_url <- payload[["next"]]
    current_url <- if (base::is.null(next_page_url) || !base::nzchar(next_page_url)) "" else base::as.character(next_page_url)[1]
  }
  dplyr::bind_rows(pages)
}

.dnmb_pazy_flatten_metadata <- function(metadata_tbl) {
  if (base::is.null(metadata_tbl) || !base::is.data.frame(metadata_tbl) || !base::nrow(metadata_tbl)) {
    return(data.frame())
  }

  rows <- lapply(base::seq_len(base::nrow(metadata_tbl)), function(i) {
    row <- metadata_tbl[i, , drop = FALSE]
    substrates <- row$substrates[[1]]
    accessions <- row$accessions[[1]]
    organism <- row$organism[[1]]
    literature <- row$literature[[1]]
    tibble::tibble(
      pazy_id = base::paste0("PAZY_", base::as.character(row$id[[1]])),
      pazy_name = base::as.character(row$name[[1]]),
      verified = base::isTRUE(row$verified[[1]]),
      substrate_names = if (base::is.data.frame(substrates) && base::nrow(substrates)) base::paste(substrates$name, collapse = "; ") else NA_character_,
      substrate_abbreviations = if (base::is.data.frame(substrates) && base::nrow(substrates)) base::paste(substrates$abbreviation, collapse = "; ") else NA_character_,
      organism_name = if (base::is.list(organism)) base::as.character(organism$scientific_name %||% NA_character_) else NA_character_,
      domain = if (base::is.list(organism)) base::as.character(organism$domain %||% NA_character_) else NA_character_,
      phylum = if (base::is.list(organism)) base::as.character(organism$phylum %||% NA_character_) else NA_character_,
      accession_list = if (base::is.data.frame(accessions) && base::nrow(accessions)) base::paste(accessions$accession, collapse = "; ") else NA_character_,
      accession_db_list = if (base::is.data.frame(accessions) && base::nrow(accessions)) base::paste(accessions$database, collapse = "; ") else NA_character_,
      literature_doi = if (base::is.data.frame(literature) && base::nrow(literature)) base::paste(literature$doi, collapse = "; ") else NA_character_
    )
  })

  dplyr::bind_rows(rows)
}

.dnmb_pazy_rewrite_fasta_with_ids <- function(source_fasta, dest_fasta, metadata_flat) {
  lines <- base::readLines(source_fasta, warn = FALSE)
  out <- character()
  current_id <- NA_character_
  for (line in lines) {
    if (base::startsWith(line, ">")) {
      header <- base::sub("^>", "", line)
      source_id <- base::strsplit(header, " ", fixed = TRUE)[[1]][1]
      source_id <- base::trimws(source_id)
      current_id <- if (base::startsWith(source_id, "PAZY_")) source_id else base::paste0("PAZY_", source_id)
      label <- metadata_flat$pazy_name[base::match(current_id, metadata_flat$pazy_id)]
      label <- ifelse(base::is.na(label) | !base::nzchar(label), current_id, label)
      out <- base::c(out, base::paste0(">", current_id, " ", label))
    } else {
      out <- base::c(out, line)
    }
  }
  base::writeLines(out, dest_fasta)
  invisible(dest_fasta)
}

.dnmb_pazy_prepare_blast_db <- function(fasta_path, db_prefix, trace_log = NULL) {
  output_files <- base::paste0(db_prefix, c(".phr", ".pin", ".psq"))
  if (base::all(base::file.exists(output_files))) {
    return(list(ok = TRUE, status = "cached", detail = db_prefix))
  }
  stage_dir <- tempfile("dnmb-pazy-makeblastdb-")
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  staged_faa <- base::file.path(stage_dir, "pazy_reference.faa")
  staged_db <- base::file.path(stage_dir, "pazy_reference")
  base::file.copy(fasta_path, staged_faa, overwrite = TRUE)
  args <- c("-in", staged_faa, "-dbtype", "prot", "-out", staged_db)
  if (!base::is.null(trace_log) && base::nzchar(trace_log)) {
    .dnmb_pazy_trace(trace_log, base::sprintf("[%s] makeblastdb %s", base::Sys.time(), .dnmb_format_command("makeblastdb", args)))
  }
  run <- dnmb_run_external("makeblastdb", args = args, required = FALSE)
  staged_files <- base::paste0(staged_db, c(".phr", ".pin", ".psq"))
  ok <- base::isTRUE(run$ok) && base::all(base::file.exists(staged_files))
  if (ok) {
    # Copy all BLAST DB files including v5 LMDB files (.pdb, .pot, .ptf, .pto)
    all_staged <- base::list.files(stage_dir, pattern = base::paste0("^", base::basename(staged_db), "\\."),
                                   full.names = TRUE)
    dest_dir <- base::dirname(db_prefix)
    for (sf in all_staged) {
      dest_name <- base::file.path(dest_dir, sub(base::basename(staged_db), base::basename(db_prefix), base::basename(sf), fixed = TRUE))
      base::file.copy(sf, dest_name, overwrite = TRUE)
    }
  }
  list(
    ok = ok,
    status = if (ok) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    detail = if (ok) db_prefix else (run$error %||% db_prefix),
    command = run
  )
}

dnmb_pazy_install_module <- function(version = .dnmb_pazy_default_version(),
                                     cache_root = NULL,
                                     install = TRUE,
                                     metadata_url = .dnmb_pazy_default_metadata_url(),
                                     fasta_url = .dnmb_pazy_default_fasta_url(),
                                     asset_urls = NULL) {
  module <- .dnmb_pazy_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_pazy_asset_layout(module_dir)
  trace_log <- base::file.path(module_dir, "pazy_install_trace.log")
  asset_urls <- .dnmb_pazy_normalize_asset_urls(asset_urls)
  metadata_source <- asset_urls$metadata_tsv %||% asset_urls$metadata_json %||% metadata_url
  fasta_source <- asset_urls$fasta_faa %||% asset_urls$reference_fasta %||% fasta_url
  status <- .dnmb_pazy_empty_status()

  if (base::all(base::file.exists(c(layout$reference_fasta, layout$metadata_tsv, layout$metadata_json, layout$blast_db_files)))) {
    manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
    if (!base::is.null(manifest) && base::isTRUE(manifest$install_ok)) {
      return(list(
        ok = TRUE,
        status = .dnmb_pazy_status_row("pazy_install", "cached", module_dir),
        files = c(list(reference_fasta = layout$reference_fasta, metadata_tsv = layout$metadata_tsv), stats::setNames(as.list(layout$blast_db_files), base::basename(layout$blast_db_files))),
        manifest = manifest
      ))
    }
  }

  if (!base::isTRUE(install) && !base::all(base::file.exists(c(layout$reference_fasta, layout$metadata_tsv)))) {
    return(list(
      ok = FALSE,
      status = .dnmb_pazy_status_row("pazy_install", "missing", "PAZy reference is missing and module_install is FALSE."),
      files = list(),
      manifest = NULL
    ))
  }

  if (base::interactive() && base::isTRUE(install) && !base::all(base::file.exists(c(layout$reference_fasta, layout$metadata_tsv)))) {
    answer <- utils::askYesNo("PAZy reference data are missing. Download the current PAZy protein set from the official API?")
    if (!base::isTRUE(answer)) {
      return(list(
        ok = FALSE,
        status = .dnmb_pazy_status_row("pazy_install", "declined", "User declined PAZy download."),
        files = list(),
        manifest = NULL
      ))
    }
  }

  metadata_flat <- if (base::dir.exists(metadata_source)) {
    utils::read.delim(base::file.path(metadata_source, "pazy_metadata.tsv"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  } else if (base::file.exists(metadata_source) && grepl("\\.tsv$", metadata_source, ignore.case = TRUE)) {
    utils::read.delim(metadata_source, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  } else if (base::file.exists(metadata_source) && grepl("\\.json$", metadata_source, ignore.case = TRUE)) {
    jsonlite::fromJSON(metadata_source, simplifyVector = TRUE)
  } else {
    .dnmb_pazy_flatten_metadata(.dnmb_pazy_fetch_metadata_pages(metadata_source))
  }
  utils::write.table(metadata_flat, file = layout$metadata_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  jsonlite::write_json(metadata_flat, path = layout$metadata_json, auto_unbox = TRUE, pretty = TRUE)
  status <- dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_metadata", "ok", layout$metadata_tsv))

  raw_fasta <- base::file.path(module_dir, "pazy_raw.faa")
  if (base::dir.exists(fasta_source)) {
    base::file.copy(base::file.path(fasta_source, "pazy_reference.faa"), raw_fasta, overwrite = TRUE)
  } else if (base::file.exists(fasta_source)) {
    base::file.copy(fasta_source, raw_fasta, overwrite = TRUE)
  } else {
    download <- .dnmb_download_asset(fasta_source, raw_fasta, insecure = FALSE)
    if (!base::isTRUE(download$ok) || !base::file.exists(raw_fasta)) {
      return(list(
        ok = FALSE,
        status = dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_fasta", "failed", download$error %||% fasta_source)),
        files = list(),
        manifest = NULL
      ))
    }
  }
  .dnmb_pazy_rewrite_fasta_with_ids(raw_fasta, layout$reference_fasta, metadata_flat)
  status <- dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_fasta", "ok", layout$reference_fasta))

  blast_prepare <- .dnmb_pazy_prepare_blast_db(layout$reference_fasta, layout$blast_db_prefix, trace_log = trace_log)
  status <- dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_prepare", blast_prepare$status, blast_prepare$detail))
  if (!base::isTRUE(blast_prepare$ok)) {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  manifest <- list(
    install_ok = TRUE,
    module = module,
    version = version,
    module_dir = module_dir,
    metadata_url = metadata_source,
    fasta_url = fasta_source,
    reference_fasta = layout$reference_fasta,
    metadata_tsv = layout$metadata_tsv,
    blast_db_prefix = layout$blast_db_prefix
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)

  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_pazy_status_row("pazy_install", "ok", module_dir)),
    files = c(
      list(reference_fasta = layout$reference_fasta, metadata_tsv = layout$metadata_tsv, trace_log = trace_log),
      stats::setNames(as.list(layout$blast_db_files), base::basename(layout$blast_db_files))
    ),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_pazy_get_module <- function(version = .dnmb_pazy_default_version(),
                                 cache_root = NULL,
                                 required = FALSE) {
  manifest <- dnmb_db_read_manifest(.dnmb_pazy_module_name(), version, cache_root = cache_root, required = FALSE)
  layout <- .dnmb_pazy_asset_layout(.dnmb_db_module_dir(.dnmb_pazy_module_name(), version, cache_root = cache_root, create = FALSE))
  ok <- !base::is.null(manifest) && base::isTRUE(manifest$install_ok) && base::all(base::file.exists(c(layout$reference_fasta, layout$metadata_tsv, layout$blast_db_files)))
  if (required && !ok) {
    base::stop("PAZy module is not installed for version `", version, "`.", call. = FALSE)
  }
  list(
    ok = ok,
    manifest = manifest,
    files = list(
      reference_fasta = layout$reference_fasta,
      metadata_tsv = layout$metadata_tsv,
      blast_db_prefix = layout$blast_db_prefix,
      blast_db_files = layout$blast_db_files
    )
  )
}

dnmb_pazy_parse_metadata <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

dnmb_pazy_parse_blast_tabular <- function(path, metadata) {
  if (!base::file.exists(path) || !base::isTRUE(base::file.info(path)$size > 0)) {
    return(data.frame())
  }
  tbl <- utils::read.delim(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (base::ncol(tbl) < 12L) {
    return(data.frame())
  }
  tbl <- tbl[, base::seq_len(base::min(13L, base::ncol(tbl))), drop = FALSE]
  base::colnames(tbl) <- c(
    "query", "subject_id", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    if (base::ncol(tbl) >= 13L) "subject_title" else NULL
  )
  tbl$pident <- base::as.numeric(tbl$pident)
  tbl$length <- base::as.integer(tbl$length)
  tbl$qlen <- base::as.integer(tbl$qlen)
  tbl$slen <- base::as.integer(tbl$slen)
  tbl$evalue <- base::as.numeric(tbl$evalue)
  tbl$bitscore <- base::as.numeric(tbl$bitscore)
  tbl$qcov <- tbl$length / tbl$qlen
  tbl$scov <- tbl$length / tbl$slen
  tbl <- dplyr::left_join(tbl, metadata, by = c("subject_id" = "pazy_id"))
  tbl
}

dnmb_pazy_normalize_hits <- function(hits) {
  if (base::is.null(hits) || !base::is.data.frame(hits) || !base::nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  support <- base::paste0(
    "pazy_id=", hits$subject_id,
    "; evalue=", base::signif(hits$evalue, 4),
    "; bitscore=", base::signif(hits$bitscore, 4),
    "; pident=", base::sprintf("%.1f", hits$pident),
    "; qcov=", base::sprintf("%.3f", hits$qcov),
    "; scov=", base::sprintf("%.3f", hits$scov),
    "; verified=", ifelse(hits$verified %in% TRUE, "TRUE", "FALSE"),
    "; substrate=", hits$substrate_abbreviations
  )

  out <- base::data.frame(
    query = .dnmb_module_clean_annotation_key(hits$query),
    source = "pazy",
    family_system = "PAZy",
    family_id = base::as.character(hits$pazy_name),
    hit_label = base::as.character(hits$pazy_name),
    enzyme_role = "plastic-active enzyme",
    evidence_mode = "direct",
    substrate_label = base::as.character(hits$substrate_abbreviations),
    support = support,
    pazy_id = base::as.character(hits$subject_id),
    verified = base::as.logical(hits$verified),
    accession_list = base::as.character(hits$accession_list),
    organism_name = base::as.character(hits$organism_name),
    pident = base::as.numeric(hits$pident),
    alignment_length = base::as.integer(hits$length),
    query_length = base::as.integer(hits$qlen),
    subject_length = base::as.integer(hits$slen),
    evalue = base::as.numeric(hits$evalue),
    bitscore = base::as.numeric(hits$bitscore),
    qcov = base::as.numeric(hits$qcov),
    scov = base::as.numeric(hits$scov),
    typing_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  out <- out[, c(.dnmb_module_optional_long_columns(), base::setdiff(base::names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.dnmb_pazy_append_evalue_threshold <- function() {
  1e-10
}

.dnmb_pazy_append_qcov_threshold <- function() {
  0.5
}

.dnmb_pazy_append_scov_threshold <- function() {
  0.5
}

.dnmb_pazy_hits_for_output <- function(hits,
                                       evalue_threshold = .dnmb_pazy_append_evalue_threshold(),
                                       qcov_threshold = .dnmb_pazy_append_qcov_threshold(),
                                       scov_threshold = .dnmb_pazy_append_scov_threshold()) {
  if (base::is.null(hits) || !base::is.data.frame(hits) || !base::nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  out <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::all(c("evalue", "qcov", "scov") %in% base::names(out))) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out$evalue <- suppressWarnings(base::as.numeric(out$evalue))
  out$qcov <- suppressWarnings(base::as.numeric(out$qcov))
  out$scov <- suppressWarnings(base::as.numeric(out$scov))
  out <- out[
    !base::is.na(out$evalue) &
      !base::is.na(out$qcov) &
      !base::is.na(out$scov) &
      out$evalue <= base::as.numeric(evalue_threshold)[1] &
      out$qcov >= base::as.numeric(qcov_threshold)[1] &
      out$scov >= base::as.numeric(scov_threshold)[1],
    ,
    drop = FALSE
  ]
  base::rownames(out) <- NULL
  if (!base::nrow(out)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out
}

.dnmb_pazy_output_table <- function(genes,
                                    hits,
                                    evalue_threshold = .dnmb_pazy_append_evalue_threshold(),
                                    qcov_threshold = .dnmb_pazy_append_qcov_threshold(),
                                    scov_threshold = .dnmb_pazy_append_scov_threshold()) {
  filtered_hits <- .dnmb_pazy_hits_for_output(
    hits = hits,
    evalue_threshold = evalue_threshold,
    qcov_threshold = qcov_threshold,
    scov_threshold = scov_threshold
  )
  out <- .dnmb_module_output_table(genes = genes, hits = filtered_hits)
  drop_cols <- base::intersect(
    c("enzyme_role", "evidence_mode", "support", "typing_eligible", "verified"),
    base::names(out)
  )
  if (base::length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  base_cols <- base::intersect(dnmb_backbone_columns(), base::names(out))
  pazy_cols <- c("family_id", "hit_label", "substrate_label", base::setdiff(base::names(out), c(base_cols, "family_id", "hit_label", "substrate_label")))
  out[, c(base_cols, pazy_cols), drop = FALSE]
}

dnmb_run_pazy_module <- function(genes,
                                 output_dir,
                                 version = .dnmb_pazy_default_version(),
                                 cache_root = NULL,
                                 install = TRUE,
                                 metadata_url = .dnmb_pazy_default_metadata_url(),
                                 fasta_url = .dnmb_pazy_default_fasta_url(),
                                 asset_urls = NULL,
                                 cpu = 1L,
                                 genbank = NULL,
                                 search_backend = c("diamond", "blastp")) {
  search_backend <- base::match.arg(search_backend)
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- base::file.path(output_dir, "pazy_module_trace.log")
  fasta_path <- base::file.path(output_dir, "pazy_query_proteins.faa")
  existing_faa <- dnmb_resolve_query_faa(genbank = genbank, output_dir = output_dir, fallback_filename = base::basename(fasta_path))
  if (!base::is.null(existing_faa) && .dnmb_can_reuse_query_fasta(existing_faa, genes)) {
    proteins <- .dnmb_prepare_query_proteins(genes)
    fasta <- list(path = existing_faa, n = base::nrow(proteins), proteins = proteins)
    fasta_path <- existing_faa
  } else {
    fasta <- .dnmb_write_query_fasta(genes, fasta_path)
  }

  status <- .dnmb_pazy_status_row("pazy_query_fasta", if (fasta$n) "ok" else "empty", base::paste0("proteins=", fasta$n))
  if (!fasta$n) {
    return(list(ok = TRUE, status = status, files = list(query_fasta = fasta_path), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  install_result <- dnmb_pazy_install_module(
    version = version,
    cache_root = cache_root,
    install = install,
    metadata_url = metadata_url,
    fasta_url = fasta_url,
    asset_urls = asset_urls
  )
  status <- dplyr::bind_rows(status, install_result$status)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = status, files = install_result$files, hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  module <- dnmb_pazy_get_module(version = version, cache_root = cache_root, required = TRUE)
  blast_out <- base::file.path(output_dir, "pazy_blastp.tsv")
  command <- NULL

  # --- diamond backend ---
  if (base::identical(search_backend, "diamond")) {
    diamond_check <- dnmb_run_external("diamond", args = "version", required = FALSE)
    if (base::nzchar(diamond_check$resolved_command)) {
      dmnd_db <- base::file.path(output_dir, "pazy_reference.dmnd")
      if (!base::file.exists(dmnd_db)) {
        makedb_run <- dnmb_run_external("diamond", args = c("makedb", "--in", module$files$reference_fasta, "-d", dmnd_db), required = FALSE)
        if (!base::isTRUE(makedb_run$ok)) {
          .dnmb_pazy_trace(trace_log, base::sprintf("[%s] diamond makedb failed, falling back to blastp", base::Sys.time()))
          search_backend <- "blastp"
        }
      }
      if (base::identical(search_backend, "diamond")) {
        stage_out <- base::file.path(output_dir, "pazy_blastp.tsv")
        d_args <- c(
          "blastp", "-q", fasta_path, "-d", dmnd_db,
          "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle",
          "-k", "5", "--threads", base::as.character(base::as.integer(cpu)[1]),
          "--sensitive", "-o", stage_out
        )
        .dnmb_pazy_trace(trace_log, base::sprintf("[%s] diamond blastp --sensitive threads=%s", base::Sys.time(), cpu))
        command <- dnmb_run_external("diamond", args = d_args, required = FALSE)
      }
    } else {
      .dnmb_pazy_trace(trace_log, base::sprintf("[%s] diamond not found, falling back to blastp", base::Sys.time()))
      search_backend <- "blastp"
    }
  }

  # --- blastp backend (fallback) ---
  if (base::identical(search_backend, "blastp") || base::is.null(command)) {
    stage_dir <- tempfile("dnmb-pazy-blastp-")
    base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
    stage_query <- base::file.path(stage_dir, "query.faa")
    stage_db <- base::file.path(stage_dir, "pazy_reference")
    stage_out <- base::file.path(stage_dir, "pazy_blastp.tsv")
    base::file.copy(fasta_path, stage_query, overwrite = TRUE)
    base::file.copy(module$files$reference_fasta, base::file.path(stage_dir, "pazy_reference.faa"), overwrite = TRUE)
    db_src <- module$files$blast_db_files
    if (base::is.character(db_src)) {
      db_src <- db_src[base::file.exists(db_src)]
      db_dest <- base::paste0(stage_db, ".", tools::file_ext(db_src))
      base::file.copy(db_src, db_dest, overwrite = TRUE)
    }
    args <- c(
      "-query", stage_query, "-db", stage_db,
      "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle",
      "-max_target_seqs", "5",
      "-num_threads", base::as.character(base::as.integer(cpu)[1]),
      "-out", stage_out
    )
    .dnmb_pazy_trace(trace_log, base::sprintf("[%s] blastp %s", base::Sys.time(), .dnmb_format_command("blastp", args)))
    command <- dnmb_run_external("blastp", args = args, required = FALSE)
    if (base::isTRUE(command$ok) && base::file.exists(stage_out)) {
      base::file.copy(stage_out, blast_out, overwrite = TRUE)
    }
  }

  if (base::isTRUE(command$ok) && !base::file.exists(blast_out)) {
    base::writeLines(character(0), blast_out)
  }
  metadata <- dnmb_pazy_parse_metadata(module$files$metadata_tsv)
  base::names(metadata)[base::match("pazy_id", base::names(metadata))] <- "pazy_id"
  parsed <- if (base::isTRUE(command$ok) && base::file.exists(blast_out)) dnmb_pazy_parse_blast_tabular(blast_out, metadata = metadata) else data.frame()
  hits <- dnmb_pazy_normalize_hits(parsed)
  status <- dplyr::bind_rows(
    status,
    .dnmb_pazy_status_row(
      "pazy_blastp",
      if (base::isTRUE(command$ok)) if (base::nrow(hits)) "ok" else "empty" else if (!base::nzchar(command$resolved_command)) "missing" else "failed",
      if (base::isTRUE(command$ok)) blast_out else (command$error %||% "PAZy blastp failed")
    )
  )

  out_table <- .dnmb_pazy_output_table(genes = genes, hits = hits)

  # Write pazy_merged.tsv for downstream visualization (Panel A/B/C pub figure)
  merged_path <- base::file.path(output_dir, "pazy_merged.tsv")
  base::tryCatch({
    if (base::is.data.frame(parsed) && base::nrow(parsed) > 0L) {
      backbone <- c("locus_tag", "gene", "product", "protein_id", "contig", "start", "end", "direction")
      base_genes <- genes[, base::intersect(backbone, base::names(genes)), drop = FALSE]
      # Rename parsed columns to match _pub expectations
      p <- parsed
      p$locus_tag <- .dnmb_module_clean_annotation_key(p$query)
      col_map <- c(pazy_name = "family_id", subject_id = "pazy_id",
                   substrate_abbreviations = "substrate_label",
                   length = "alignment_length", qlen = "query_length", slen = "subject_length")
      for (old_nm in base::names(col_map)) {
        if (old_nm %in% base::names(p)) base::names(p)[base::names(p) == old_nm] <- col_map[[old_nm]]
      }
      # Best hit per locus_tag
      p <- p[base::order(p$evalue, -p$bitscore), , drop = FALSE]
      p <- p[!base::duplicated(p$locus_tag), , drop = FALSE]
      merge_cols <- base::setdiff(base::names(p), c("query", backbone))
      merged <- base::merge(base_genes, p[, c("locus_tag", merge_cols), drop = FALSE],
                            by = "locus_tag", all.x = TRUE)
      rename_idx <- !base::names(merged) %in% backbone
      base::names(merged)[rename_idx] <- base::paste0("PAZy_", base::names(merged)[rename_idx])
      utils::write.table(merged, merged_path, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }, error = function(e) NULL)

  list(
    ok = base::isTRUE(command$ok),
    status = status,
    files = c(install_result$files, list(query_fasta = fasta_path, blast_tsv = blast_out, trace_log = trace_log)),
    hits = hits,
    output_table = out_table,
    module_result = parsed
  )
}
#' Internal PAZy module helpers
#'
#' Installation, search, parsing, and output-building helpers for the DNMB
#' PAZy module workflow.
#'
#' @name dnmb_internal_pazy_module
#' @keywords internal
#' @noRd
NULL

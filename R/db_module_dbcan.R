.dnmb_dbcan_module_name <- function() {
  "dbcan"
}

.dnmb_dbcan_default_version <- function() {
  "current"
}

.dnmb_dbcan_required_tool_version <- function() {
  "5.2.9"
}

.dnmb_dbcan_pipeline_contract_version <- function() {
  2L
}

.dnmb_dbcan_default_base_url <- function() {
  "https://pro.unl.edu/dbCAN2/download/Databases"
}

.dnmb_dbcan_release_page_url <- function() {
  "https://pro.unl.edu/dbCAN2/"
}

.dnmb_dbcan_status_row <- function(component, status, detail) {
  tibble::tibble(
    component = as.character(component),
    status = as.character(status),
    detail = as.character(detail)
  )
}

.dnmb_dbcan_empty_status <- function() {
  .dnmb_dbcan_status_row(character(), character(), character())
}

.dnmb_dbcan_trace <- function(path, text) {
  cat(text, "\n", file = path, append = TRUE)
}

.dnmb_dbcan_asset_layout <- function(module_dir) {
  hmm_path <- file.path(module_dir, "dbCAN.txt")
  mapping_path <- file.path(module_dir, "fam-substrate-mapping.tsv")
  list(
    module_dir = module_dir,
    hmm_path = hmm_path,
    mapping_path = mapping_path,
    hmmpress_files = paste0(hmm_path, c(".h3f", ".h3i", ".h3m", ".h3p"))
  )
}

.dnmb_dbcan_fallback_release_info <- function() {
  list(
    version = "V14",
    release_note = "Fallback release info",
    source = "fallback",
    dbcan_release_date = NA_character_
  )
}

.dnmb_dbcan_release_info <- function(page_url = .dnmb_dbcan_release_page_url()) {
  html <- character()
  curl_detection <- dnmb_detect_binary("curl", required = FALSE)
  if (isTRUE(curl_detection$found)) {
    run <- dnmb_run_external("curl", c("-k", "-L", "-sS", page_url), required = FALSE)
    if (isTRUE(run$ok) && length(run$stdout)) {
      html <- run$stdout
    }
  }
  if (!length(html)) {
    wget_detection <- dnmb_detect_binary("wget", required = FALSE)
    if (isTRUE(wget_detection$found)) {
      run <- dnmb_run_external("wget", c("--no-check-certificate", "-qO-", page_url), required = FALSE)
      if (isTRUE(run$ok) && length(run$stdout)) {
        html <- run$stdout
      }
    }
  }
  if (!length(html)) {
    return(.dnmb_dbcan_fallback_release_info())
  }

  pattern <- "([0-9]{1,2}/[0-9]{1,2}/[0-9]{4}).*dbCAN HMMdb v([0-9]+) is released \\(based on CAZyDB ([0-9/]+)\\)"
  match_line <- html[grepl(pattern, html, perl = TRUE)]
  if (!length(match_line)) {
    return(.dnmb_dbcan_fallback_release_info())
  }

  groups <- regmatches(match_line[[1]], regexec(pattern, match_line[[1]], perl = TRUE))[[1]]
  if (length(groups) < 4L) {
    return(.dnmb_dbcan_fallback_release_info())
  }

  list(
    version = paste0("V", groups[[3]]),
    release_note = trimws(gsub("<[^>]+>", "", match_line[[1]])),
    source = page_url,
    dbcan_release_date = groups[[2]],
    cazydb_release_date = groups[[4]]
  )
}

.dnmb_dbcan_aws_s3_base <- function() {
  # dbCAN's canonical AWS S3 distribution (see
  # dbcan/constants/databases_constants.py::AWS_S3_URL in the
  # `dbcan` python package). The legacy HTTP mirror under
  # `https://bcb.unl.edu/dbCAN2/download/...` is dead — every path
  # 302-redirects to a landing page and downloads save HTML instead
  # of the actual HMM/DIAMOND/TSV files.
  "https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026"
}

.dnmb_dbcan_asset_urls <- function(base_url = .dnmb_dbcan_default_base_url(),
                                   version = .dnmb_dbcan_default_version(),
                                   asset_urls = NULL,
                                   release_info = NULL) {
  if (!is.null(asset_urls)) {
    urls <- unlist(asset_urls, use.names = TRUE)
    urls <- stats::setNames(as.character(urls), names(urls))
    if (!length(urls) || is.null(names(urls)) || any(!nzchar(names(urls)))) {
      stop("`asset_urls` must be a named character vector or list.", call. = FALSE)
    }
    return(urls)
  }

  # Ignore the legacy versioned HTTP layout entirely — the S3 mirror
  # ships a single `dbCAN.hmm` (no per-version HMMdb-Vxx.txt split) and
  # the same `fam-substrate-mapping.tsv`. DNMB saves the HMM locally
  # as `dbCAN.txt` which hmmsearch still happily accepts.
  s3 <- .dnmb_dbcan_aws_s3_base()
  stats::setNames(
    c(
      paste0(s3, "/dbCAN.hmm"),
      paste0(s3, "/fam-substrate-mapping.tsv")
    ),
    c("dbCAN.txt", "fam-substrate-mapping.tsv")
  )
}

.dnmb_dbcan_is_hmm <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(FALSE)
  }
  first_lines <- tryCatch(readLines(path, n = 5L, warn = FALSE), error = function(...) character())
  any(grepl("^HMMER", first_lines))
}

.dnmb_dbcan_is_mapping <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(FALSE)
  }
  first_lines <- tryCatch(readLines(path, n = 5L, warn = FALSE), error = function(...) character())
  !any(grepl("^<!DOCTYPE html|^<html", first_lines, ignore.case = TRUE)) && any(grepl("\t", first_lines))
}

.dnmb_dbcan_hmmpress_summary <- function(hmmpress, args = character()) {
  paste(
    c(
      paste0("path=", if (isTRUE(hmmpress$found) && nzchar(hmmpress$path)) hmmpress$path else "missing"),
      if (!is.null(hmmpress$version)) paste0("version=", hmmpress$version),
      if (length(args)) paste0("args=", .dnmb_format_command("hmmpress", args))
    ),
    collapse = "; "
  )
}

.dnmb_dbcan_detect_hmmer <- function(binary) {
  detection <- dnmb_detect_binary(binary, required = FALSE)
  version <- NA_character_
  if (isTRUE(detection$found)) {
    run <- dnmb_run_external(binary, args = "-h", required = FALSE)
    version <- .dnmb_parse_tool_version(c(run$stdout, run$stderr))
  }
  list(
    binary = binary,
    path = detection$path,
    found = isTRUE(detection$found),
    version = version
  )
}

.dnmb_dbcan_tool_version <- function() {
  detection <- dnmb_detect_binary("run_dbcan", required = FALSE)
  if (!isTRUE(detection$found)) return(NA_character_)
  run <- dnmb_run_external(detection$path, args = "version", required = FALSE)
  .dnmb_parse_tool_version(c(run$stdout, run$stderr))
}

.dnmb_dbcan_prepare_hmm_db <- function(hmm_path, force = FALSE, trace_log = NULL) {
  pressed_files <- paste0(hmm_path, c(".h3f", ".h3i", ".h3m", ".h3p"))
  hmmpress <- .dnmb_dbcan_detect_hmmer("hmmpress")
  args <- c("-f", hmm_path)
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_dbcan_trace(trace_log, sprintf("[%s] hmmpress %s", Sys.time(), .dnmb_dbcan_hmmpress_summary(hmmpress, args)))
  }

  if (!isTRUE(force) && all(file.exists(pressed_files))) {
    return(list(
      ok = TRUE,
      status = "cached",
      detail = paste("dbCAN HMMER database files already exist", .dnmb_dbcan_hmmpress_summary(hmmpress, args), sep = "; "),
      hmm_path = hmm_path,
      hmmpress_files = pressed_files,
      command = NULL,
      hmmpress = hmmpress,
      args = args
    ))
  }

  run <- dnmb_run_external("hmmpress", args = args, required = FALSE)
  ok <- isTRUE(run$ok) && all(file.exists(pressed_files))
  detail <- if (ok) {
    paste(hmm_path, .dnmb_dbcan_hmmpress_summary(hmmpress, args), sep = "; ")
  } else {
    stderr_text <- .dnmb_compact_output(c(run$stderr, run$stdout), max_lines = 6L)
    paste(
      c(
        paste0("path=", if (isTRUE(hmmpress$found) && nzchar(hmmpress$path)) hmmpress$path else "missing"),
        paste0("version=", hmmpress$version %||% "unknown"),
        paste0("args=", .dnmb_format_command("hmmpress", args)),
        if (nzchar(stderr_text)) paste0("stderr=", stderr_text)
      ),
      collapse = "; "
    )
  }
  status <- if (ok) "ok" else if (!nzchar(run$resolved_command)) "missing" else "failed"
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_dbcan_trace(trace_log, sprintf("[%s] hmmpress_result status=%s detail=%s", Sys.time(), status, detail))
  }

  list(
    ok = ok,
    status = status,
    detail = detail,
    hmm_path = hmm_path,
    hmmpress_files = pressed_files[file.exists(pressed_files)],
    command = run,
    hmmpress = hmmpress,
    args = args
  )
}

.dnmb_dbcan_prepare_pul_blast_db <- function(pul_faa, trace_log = NULL) {
  output_files <- paste0(pul_faa, c(".phr", ".pin", ".psq"))
  if (all(file.exists(output_files))) {
    return(list(
      ok = TRUE,
      status = "cached",
      detail = pul_faa,
      files = output_files
    ))
  }

  stage_dir <- tempfile("dnmb-dbcan-pul-blast-")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  staged_faa <- file.path(stage_dir, "PUL.faa")
  file.copy(pul_faa, staged_faa, overwrite = TRUE)
  args <- c("-in", staged_faa, "-dbtype", "prot", "-out", staged_faa)
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_dbcan_trace(trace_log, sprintf("[%s] pul_makeblastdb args=%s", Sys.time(), .dnmb_format_command("makeblastdb", args)))
  }
  run <- dnmb_run_external("makeblastdb", args = args, required = FALSE)
  staged_files <- paste0(staged_faa, c(".phr", ".pin", ".psq"))
  ok <- isTRUE(run$ok) && all(file.exists(staged_files))
  if (ok) {
    file.copy(staged_files, output_files, overwrite = TRUE)
  }
  list(
    ok = ok,
    status = if (ok) "ok" else if (!nzchar(run$resolved_command)) "missing" else "failed",
    detail = if (ok) pul_faa else (run$error %||% pul_faa),
    files = output_files[file.exists(output_files)],
    command = run,
    args = args
  )
}

.dnmb_dbcan_manifest_diagnostics <- function(layout, prepare_result = NULL, release_info = NULL) {
  hmmpress <- prepare_result$hmmpress %||% .dnmb_dbcan_detect_hmmer("hmmpress")
  args <- prepare_result$args %||% c("-f", layout$hmm_path)
  info <- release_info %||% .dnmb_dbcan_fallback_release_info()

  list(
    expected_tool_version = .dnmb_dbcan_required_tool_version(),
    tool_version = .dnmb_dbcan_tool_version(),
    hmmpress_path = if (isTRUE(hmmpress$found) && nzchar(hmmpress$path)) hmmpress$path else "",
    hmmpress_version = hmmpress$version %||% NA_character_,
    hmmpress_args = as.character(args),
    resolved_release_version = if (identical(info$source %||% "", "fallback")) NA_character_ else info$version %||% NA_character_,
    release_info_source = info$source %||% NA_character_,
    dbcan_release_date = info$dbcan_release_date %||% NA_character_,
    cazydb_release_date = info$cazydb_release_date %||% NA_character_,
    release_note = info$release_note %||% NA_character_
  )
}

.dnmb_dbcan_remote_asset_state <- function(urls, enabled = TRUE) {
  if (!isTRUE(enabled) || !length(urls)) {
    return(NULL)
  }

  rows <- lapply(names(urls), function(asset_name) {
    metadata <- .dnmb_remote_asset_metadata(urls[[asset_name]], insecure = TRUE)
    tibble::tibble(
      asset_name = asset_name,
      url = metadata$url %||% NA_character_,
      ok = isTRUE(metadata$ok),
      last_modified = metadata$last_modified %||% NA_character_,
      content_length = metadata$content_length %||% NA_real_,
      etag = metadata$etag %||% NA_character_
    )
  })

  dplyr::bind_rows(rows)
}

.dnmb_dbcan_remote_update_needed <- function(manifest, remote_state) {
  if (is.null(manifest) || is.null(remote_state) || !is.data.frame(remote_state) || !nrow(remote_state)) {
    return(FALSE)
  }
  if (!all(remote_state$ok %in% TRUE)) {
    return(FALSE)
  }

  previous <- manifest$remote_asset_state
  if (is.null(previous) || !is.data.frame(previous) || !nrow(previous)) {
    return(TRUE)
  }

  previous <- previous[, intersect(c("asset_name", "last_modified", "content_length", "etag"), names(previous)), drop = FALSE]
  current <- remote_state[, c("asset_name", "last_modified", "content_length", "etag"), drop = FALSE]
  previous <- previous[order(previous$asset_name), , drop = FALSE]
  current <- current[order(current$asset_name), , drop = FALSE]
  rownames(previous) <- NULL
  rownames(current) <- NULL

  !identical(previous, current)
}

.dnmb_dbcan_profile_family <- function(profile_id) {
  clean <- trimws(as.character(profile_id)[1])
  clean <- sub("\\.hmm$", "", clean, ignore.case = TRUE)
  clean <- sub("\\.txt$", "", clean, ignore.case = TRUE)
  if (grepl("^GT2_", clean)) {
    return("GT2")
  }
  clean
}

.dnmb_dbcan_clean_query_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^>", "", x)
  x <- sub("[[:space:]].*$", "", x)
  x <- gsub("^gnl\\|", "", x)
  x <- gsub("\\|", ":", x)
  x <- gsub(
    "^lcl:(AC|NC|BG|NT|NW|NZ)_([A-Za-z]+)?[0-9]+\\.[0-9]+_prot_",
    "",
    x
  )
  x <- gsub("^extdb:", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

.dnmb_dbcan_missing_value <- function(x) {
  x <- trimws(as.character(x))
  is.na(x) | !nzchar(x) | tolower(x) %in% c("-", "na", "n/a", "none", "null")
}

.dnmb_dbcan_contig_keys <- function(genes) {
  n <- nrow(genes)
  display <- if ("contig" %in% names(genes)) as.character(genes$contig) else rep(NA_character_, n)
  key <- rep(NA_character_, n)

  if ("contig_number" %in% names(genes)) {
    number <- suppressWarnings(as.integer(genes$contig_number))
    use <- !is.na(number)
    key[use] <- sprintf("replicon_%04d", number[use])
  }

  accession_columns <- c("record_accession", "contig_accession", "accession", "sequence_id")
  for (column_name in accession_columns[accession_columns %in% names(genes)]) {
    value <- .dnmb_dbcan_clean_query_id(genes[[column_name]])
    use <- is.na(key) & !is.na(value) & nzchar(value)
    key[use] <- value[use]
  }

  fallback <- gsub("[^A-Za-z0-9_.-]+", "_", display)
  fallback[is.na(fallback) | !nzchar(fallback)] <- "unknown_replicon"
  key[is.na(key)] <- fallback[is.na(key)]
  list(key = key, display = display)
}

.dnmb_dbcan_build_gene_map <- function(genes) {
  if (!is.data.frame(genes)) {
    stop("`genes` must be a data frame.", call. = FALSE)
  }
  n <- nrow(genes)
  value_or_na <- function(column_name) {
    if (column_name %in% names(genes)) as.character(genes[[column_name]]) else rep(NA_character_, n)
  }

  locus_tag <- .dnmb_dbcan_clean_query_id(value_or_na("locus_tag"))
  protein_id <- .dnmb_dbcan_clean_query_id(value_or_na("protein_id"))
  translation <- if ("translation" %in% names(genes)) {
    .dnmb_normalize_translation(genes$translation)
  } else {
    rep(NA_character_, n)
  }
  contigs <- .dnmb_dbcan_contig_keys(genes)
  start <- suppressWarnings(as.numeric(value_or_na("start")))
  end <- suppressWarnings(as.numeric(value_or_na("end")))
  direction <- value_or_na("direction")

  mapping_source <- ifelse(!is.na(locus_tag), "locus_tag", ifelse(!is.na(protein_id), "protein_id", "coordinate"))
  coordinate_id <- paste0(
    contigs$key,
    ":",
    ifelse(is.finite(start), format(start, scientific = FALSE, trim = TRUE), "NA"),
    "-",
    ifelse(is.finite(end), format(end, scientific = FALSE, trim = TRUE), "NA"),
    ":",
    ifelse(is.na(direction) | !nzchar(direction), ".", direction)
  )
  mapping_id <- ifelse(!is.na(locus_tag), locus_tag, ifelse(!is.na(protein_id), protein_id, coordinate_id))

  out <- tibble::tibble(
    dbcan_gene_row = seq_len(n),
    dbcan_query_id = sprintf("DNMBCAZY_%07d", seq_len(n)),
    mapping_id = mapping_id,
    mapping_source = mapping_source,
    locus_tag = locus_tag,
    protein_id = protein_id,
    dbcan_contig_key = contigs$key,
    contig = contigs$display,
    start = start,
    end = end,
    direction = direction,
    translation = translation
  )
  out
}

.dnmb_dbcan_write_gene_map <- function(id_map, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  columns <- setdiff(names(id_map), "translation")
  utils::write.table(
    as.data.frame(id_map[, columns, drop = FALSE], stringsAsFactors = FALSE),
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
  invisible(path)
}

.dnmb_dbcan_write_mapped_query_fasta <- function(id_map, path) {
  id_map <- id_map[!is.na(id_map$translation) & nzchar(id_map$translation), , drop = FALSE]
  fasta <- tibble::tibble(
    protein_label = id_map$dbcan_query_id,
    protein_seq = id_map$translation
  )
  .dnmb_write_protein_fasta(fasta, path)
  list(path = path, proteins = id_map, n = nrow(id_map))
}

.dnmb_dbcan_restore_query_map <- function(tbl, id_map) {
  if (is.null(tbl) || !is.data.frame(tbl) || !nrow(tbl) ||
      is.null(id_map) || !is.data.frame(id_map) || !nrow(id_map) ||
      !"query" %in% names(tbl)) {
    return(tbl)
  }
  out <- as.data.frame(tbl, stringsAsFactors = FALSE)
  query_id <- .dnmb_dbcan_clean_query_id(out$query)
  idx <- match(query_id, id_map$dbcan_query_id)
  out$dbcan_query_id <- query_id
  out$dbcan_gene_row <- id_map$dbcan_gene_row[idx]
  out$dbcan_mapping_source <- id_map$mapping_source[idx]
  out$dbcan_contig_key <- id_map$dbcan_contig_key[idx]
  mapped <- !is.na(idx)
  out$query[mapped] <- id_map$mapping_id[idx[mapped]]
  out$query[!mapped] <- query_id[!mapped]
  tibble::as_tibble(out)
}

.dnmb_dbcan_raw_domain_tokens <- function(call) {
  x <- trimws(as.character(call)[1])
  if (.dnmb_dbcan_missing_value(x)) return(character())
  tokens <- trimws(unlist(strsplit(x, "[|+;]", perl = TRUE), use.names = FALSE))
  tokens <- sub("\\([^)]*\\)$", "", tokens)
  tokens <- tokens[! .dnmb_dbcan_missing_value(tokens)]
  tokens
}

.dnmb_dbcan_family_tokens <- function(call) {
  tokens <- .dnmb_dbcan_raw_domain_tokens(call)
  if (!length(tokens)) return(character())
  tokens <- sub("\\.hmm$", "", tokens, ignore.case = TRUE)
  tokens <- sub("_e[0-9]+$", "", tokens, ignore.case = TRUE)
  family <- sub(
    "^((?:AA|CBM|CE|GH|GT|PL)[0-9]+(?:_[0-9]+)?).*$",
    "\\1",
    tokens,
    perl = TRUE,
    ignore.case = TRUE
  )
  valid <- grepl("^(AA|CBM|CE|GH|GT|PL)[0-9]+(?:_[0-9]+)?$", family, perl = TRUE, ignore.case = TRUE)
  unique(toupper(family[valid]))
}

.dnmb_dbcan_primary_family <- function(families) {
  families <- as.character(families)
  families <- families[!is.na(families) & nzchar(families)]
  if (!length(families)) return(NA_character_)
  catalytic <- families[grepl("^(AA|CE|GH|GT|PL)", families)]
  if (length(catalytic)) catalytic[[1]] else families[[1]]
}

.dnmb_dbcan_normalize_substrate <- function(x) {
  x <- trimws(as.character(x)[1])
  if (.dnmb_dbcan_missing_value(x)) return(NA_character_)
  values <- trimws(unlist(strsplit(x, "[;,]", perl = TRUE), use.names = FALSE))
  values <- values[! .dnmb_dbcan_missing_value(values)]
  values <- gsub("^hostglycan$", "host glycan", values, ignore.case = TRUE)
  values <- unique(values)
  if (!length(values)) NA_character_ else paste(values, collapse = "; ")
}

.dnmb_dbcan_empty_overview <- function() {
  tibble::tibble(
    query = character(),
    dbcan_ec = character(),
    dbcan_hmm_call = character(),
    dbcan_subfamily_call = character(),
    dbcan_diamond_call = character(),
    dbcan_tool_count = integer(),
    dbcan_recommended = character(),
    dbcan_overview_substrate = character(),
    dbcan_evidence_sources = character(),
    dbcan_evidence_tier = character(),
    dbcan_consensus = logical(),
    dbcan_domain_architecture = character(),
    dbcan_all_families = character(),
    dbcan_primary_family = character(),
    dbcan_family_count = integer()
  )
}

dnmb_dbcan_parse_overview <- function(path, id_map = NULL) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(.dnmb_dbcan_empty_overview())
  }
  tbl <- tryCatch(
    utils::read.delim(
      path,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      check.names = FALSE,
      quote = "",
      comment.char = ""
    ),
    error = function(e) NULL
  )
  if (is.null(tbl) || !nrow(tbl)) return(.dnmb_dbcan_empty_overview())

  column <- function(candidates, default = NA_character_) {
    found <- candidates[candidates %in% names(tbl)]
    if (length(found)) as.character(tbl[[found[[1]]]]) else rep(default, nrow(tbl))
  }
  query <- .dnmb_dbcan_clean_query_id(column(c("Gene ID", "Gene.ID", "query", "Protein ID")))
  ec <- column(c("EC#", "EC", "EC number"))
  hmm <- column(c("dbCAN_hmm", "HMMER"))
  subfamily <- column(c("dbCAN_sub", "dbCAN-sub"))
  diamond <- column(c("DIAMOND", "diamond"))
  recommended <- column(c("Recommend Results", "Recommended Results", "Recommend.Results"))
  substrate <- column(c("Substrate", "substrate"))
  tool_count <- suppressWarnings(as.integer(column(c("#ofTools", "ofTools", "Tools"))))

  has_hmm <- !.dnmb_dbcan_missing_value(hmm)
  has_sub <- !.dnmb_dbcan_missing_value(subfamily)
  has_diamond <- !.dnmb_dbcan_missing_value(diamond)
  inferred_count <- as.integer(has_hmm) + as.integer(has_sub) + as.integer(has_diamond)
  tool_count[is.na(tool_count)] <- inferred_count[is.na(tool_count)]

  chosen <- recommended
  for (candidate in list(subfamily, hmm, diamond)) {
    take <- .dnmb_dbcan_missing_value(chosen) & !.dnmb_dbcan_missing_value(candidate)
    chosen[take] <- candidate[take]
  }
  family_list <- lapply(chosen, .dnmb_dbcan_family_tokens)
  domain_list <- lapply(chosen, .dnmb_dbcan_raw_domain_tokens)
  all_families <- vapply(
    family_list,
    function(x) if (length(x)) paste(x, collapse = "; ") else NA_character_,
    character(1)
  )
  architecture <- vapply(
    domain_list,
    function(x) if (length(x)) paste(x, collapse = " | ") else NA_character_,
    character(1)
  )
  primary <- vapply(family_list, .dnmb_dbcan_primary_family, character(1))
  family_count <- vapply(family_list, length, integer(1))
  sources <- vapply(seq_len(nrow(tbl)), function(i) {
    source <- c("HMM", "dbCAN-sub", "DIAMOND")[c(has_hmm[[i]], has_sub[[i]], has_diamond[[i]])]
    if (length(source)) paste(source, collapse = "+") else NA_character_
  }, character(1))
  consensus <- tool_count >= 2L & !.dnmb_dbcan_missing_value(recommended)
  tier <- ifelse(
    tool_count >= 3L & consensus,
    "very_high",
    ifelse(
      consensus,
      "high",
      ifelse(has_hmm, "medium", ifelse(has_sub, "low", "audit"))
    )
  )

  out <- tibble::tibble(
    query = query,
    dbcan_ec = vapply(ec, function(x) {
      value <- trimws(as.character(x))
      if (.dnmb_dbcan_missing_value(value)) NA_character_ else value
    }, character(1)),
    dbcan_hmm_call = ifelse(has_hmm, hmm, NA_character_),
    dbcan_subfamily_call = ifelse(has_sub, subfamily, NA_character_),
    dbcan_diamond_call = ifelse(has_diamond, diamond, NA_character_),
    dbcan_tool_count = tool_count,
    dbcan_recommended = ifelse(.dnmb_dbcan_missing_value(recommended), NA_character_, recommended),
    dbcan_overview_substrate = vapply(substrate, .dnmb_dbcan_normalize_substrate, character(1)),
    dbcan_evidence_sources = sources,
    dbcan_evidence_tier = tier,
    dbcan_consensus = consensus,
    dbcan_domain_architecture = architecture,
    dbcan_all_families = all_families,
    dbcan_primary_family = primary,
    dbcan_family_count = family_count
  )
  out <- out[!is.na(out$query) & nzchar(out$query) & !is.na(out$dbcan_all_families), , drop = FALSE]
  .dnmb_dbcan_restore_query_map(out, id_map)
}

.dnmb_dbcan_family_substrate_map <- function(path) {
  empty <- tibble::tibble(
    family = character(),
    dbcan_family_substrate_prior = character()
  )
  if (is.null(path) || !file.exists(path) || !isTRUE(file.info(path)$size > 0)) return(empty)
  tbl <- tryCatch(
    utils::read.delim(
      path,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      check.names = FALSE,
      comment.char = "#",
      quote = "\""
    ),
    error = function(e) NULL
  )
  if (is.null(tbl) || !nrow(tbl)) return(empty)
  names(tbl) <- trimws(names(tbl))
  required <- c("Family", "Substrate_high_level")
  if (!all(required %in% names(tbl))) return(empty)
  removed <- if ("type" %in% names(tbl)) toupper(trimws(as.character(tbl$type))) == "REMOVED" else rep(FALSE, nrow(tbl))
  preferred <- if ("new_Substrate_high_level" %in% names(tbl)) as.character(tbl$new_Substrate_high_level) else rep(NA_character_, nrow(tbl))
  original <- as.character(tbl$Substrate_high_level)
  substrate <- ifelse(!.dnmb_dbcan_missing_value(preferred), preferred, original)
  family <- toupper(trimws(as.character(tbl$Family)))
  keep <- !removed & !.dnmb_dbcan_missing_value(family) & !.dnmb_dbcan_missing_value(substrate)
  if (!any(keep)) return(empty)
  data.frame(family = family[keep], substrate = substrate[keep], stringsAsFactors = FALSE) |>
    dplyr::group_by(.data$family) |>
    dplyr::summarise(
      dbcan_family_substrate_prior = paste(sort(unique(.data$substrate)), collapse = "; "),
      .groups = "drop"
    )
}

.dnmb_dbcan_add_family_substrate_prior <- function(overview, mapping_path) {
  if (is.null(overview) || !is.data.frame(overview) || !nrow(overview)) return(overview)
  mapping <- .dnmb_dbcan_family_substrate_map(mapping_path)
  out <- as.data.frame(overview, stringsAsFactors = FALSE)
  out$dbcan_family_substrate_prior <- NA_character_
  if (!nrow(mapping)) return(tibble::as_tibble(out))
  prior_map <- stats::setNames(mapping$dbcan_family_substrate_prior, mapping$family)
  out$dbcan_family_substrate_prior <- vapply(out$dbcan_all_families, function(value) {
    families <- .dnmb_dbcan_family_tokens(value)
    if (!length(families)) return(NA_character_)
    substrates <- unname(prior_map[families])
    missing <- is.na(substrates) | !nzchar(substrates)
    if (any(missing)) {
      parent <- sub("^([A-Z]+[0-9]+).*$", "\\1", families[missing])
      substrates[missing] <- unname(prior_map[parent])
    }
    substrates <- unlist(strsplit(substrates[!is.na(substrates) & nzchar(substrates)], ";\\s*", perl = TRUE), use.names = FALSE)
    substrates <- sort(unique(trimws(substrates[nzchar(trimws(substrates))])))
    if (length(substrates)) paste(substrates, collapse = "; ") else NA_character_
  }, character(1))
  tibble::as_tibble(out)
}

.dnmb_dbcan_merge_overview_hits <- function(hits, overview) {
  if (is.null(overview) || !is.data.frame(overview) || !nrow(overview)) return(hits)
  hits <- as.data.frame(hits, stringsAsFactors = FALSE)
  overview <- as.data.frame(overview, stringsAsFactors = FALSE)

  hit_idx <- if (nrow(hits) && "dbcan_gene_row" %in% names(hits) && "dbcan_gene_row" %in% names(overview)) {
    match(hits$dbcan_gene_row, overview$dbcan_gene_row)
  } else if (nrow(hits)) {
    match(.dnmb_dbcan_clean_query_id(hits$query), .dnmb_dbcan_clean_query_id(overview$query))
  } else integer()
  if (nrow(hits)) {
    metadata_columns <- c(
      "dbcan_tool_count", "dbcan_evidence_sources", "dbcan_evidence_tier",
      "dbcan_consensus", "dbcan_all_families", "dbcan_primary_family",
      "dbcan_domain_architecture", "dbcan_ec", "dbcan_overview_substrate",
      "dbcan_family_substrate_prior"
    )
    for (column_name in metadata_columns[metadata_columns %in% names(overview)]) {
      hits[[column_name]] <- overview[[column_name]][hit_idx]
    }
    direct_substrate <- if ("dbcan_overview_substrate" %in% names(overview)) overview$dbcan_overview_substrate[hit_idx] else NA_character_
    hits$substrate_label <- ifelse(!is.na(direct_substrate) & nzchar(direct_substrate), direct_substrate, hits$substrate_label)
    tier <- if ("dbcan_evidence_tier" %in% names(hits)) hits$dbcan_evidence_tier else NA_character_
    hits$typing_eligible <- tier %in% c("very_high", "high", "medium")
  }

  represented <- if (nrow(hits) && "dbcan_gene_row" %in% names(hits)) {
    unique(hits$dbcan_gene_row[!is.na(hits$dbcan_gene_row)])
  } else if (nrow(hits)) {
    unique(.dnmb_dbcan_clean_query_id(hits$query))
  } else {
    integer()
  }
  missing <- if ("dbcan_gene_row" %in% names(overview) && length(represented) && is.numeric(represented)) {
    !overview$dbcan_gene_row %in% represented
  } else if ("dbcan_gene_row" %in% names(overview) && !length(represented)) {
    rep(TRUE, nrow(overview))
  } else {
    !.dnmb_dbcan_clean_query_id(overview$query) %in% represented
  }

  synthetic <- lapply(which(missing), function(i) {
    families <- .dnmb_dbcan_family_tokens(overview$dbcan_all_families[[i]])
    if (!length(families)) return(NULL)
    substrate <- overview$dbcan_overview_substrate[[i]]
    if (is.na(substrate) || !nzchar(substrate)) substrate <- overview$dbcan_family_substrate_prior[[i]] %||% NA_character_
    data.frame(
      query = rep(overview$query[[i]], length(families)),
      source = rep("dbcan_overview", length(families)),
      family_system = rep("dbCAN", length(families)),
      family_id = families,
      hit_label = families,
      enzyme_role = rep("CAZyme", length(families)),
      evidence_mode = rep(if (isTRUE(overview$dbcan_consensus[[i]])) "consensus" else "single_tool", length(families)),
      substrate_label = rep(substrate, length(families)),
      support = rep(paste0(
        "tools=", overview$dbcan_tool_count[[i]],
        "; tier=", overview$dbcan_evidence_tier[[i]],
        "; sources=", overview$dbcan_evidence_sources[[i]]
      ), length(families)),
      typing_eligible = rep(overview$dbcan_evidence_tier[[i]] %in% c("very_high", "high", "medium"), length(families)),
      dbcan_query_id = rep(overview$dbcan_query_id[[i]] %||% NA_character_, length(families)),
      dbcan_gene_row = rep(overview$dbcan_gene_row[[i]] %||% NA_integer_, length(families)),
      dbcan_mapping_source = rep(overview$dbcan_mapping_source[[i]] %||% NA_character_, length(families)),
      dbcan_contig_key = rep(overview$dbcan_contig_key[[i]] %||% NA_character_, length(families)),
      dbcan_tool_count = rep(overview$dbcan_tool_count[[i]], length(families)),
      dbcan_evidence_sources = rep(overview$dbcan_evidence_sources[[i]], length(families)),
      dbcan_evidence_tier = rep(overview$dbcan_evidence_tier[[i]], length(families)),
      dbcan_consensus = rep(overview$dbcan_consensus[[i]], length(families)),
      dbcan_all_families = rep(overview$dbcan_all_families[[i]], length(families)),
      dbcan_primary_family = rep(overview$dbcan_primary_family[[i]], length(families)),
      dbcan_domain_architecture = rep(overview$dbcan_domain_architecture[[i]], length(families)),
      dbcan_ec = rep(overview$dbcan_ec[[i]], length(families)),
      dbcan_overview_substrate = rep(overview$dbcan_overview_substrate[[i]], length(families)),
      stringsAsFactors = FALSE
    )
  })
  tibble::as_tibble(dplyr::bind_rows(hits, dplyr::bind_rows(synthetic)))
}

.dnmb_dbcan_empty_hits <- function() {
  tibble::tibble(
    query = character(),
    profile_id = character(),
    profile_length = integer(),
    gene_length = integer(),
    evalue = numeric(),
    profile_start = integer(),
    profile_end = integer(),
    gene_start = integer(),
    gene_end = integer(),
    coverage = numeric(),
    family_id = character(),
    hit_label = character(),
    substrate_label = character()
  )
}

.dnmb_dbcan_filter_overlaps <- function(hits) {
  if (is.null(hits) || !nrow(hits)) {
    return(.dnmb_dbcan_empty_hits())
  }

  split_hits <- split(as.data.frame(hits, stringsAsFactors = FALSE), hits$query)
  filtered <- lapply(split_hits, function(tbl) {
    tbl <- tbl[order(tbl$gene_start, tbl$gene_end, tbl$evalue), , drop = FALSE]
    if (nrow(tbl) <= 1L) {
      return(tbl)
    }

    i <- 1L
    while (i < nrow(tbl)) {
      len1 <- tbl$gene_end[[i]] - tbl$gene_start[[i]] + 1L
      len2 <- tbl$gene_end[[i + 1L]] - tbl$gene_start[[i + 1L]] + 1L
      overlap <- min(tbl$gene_end[[i]], tbl$gene_end[[i + 1L]]) -
        max(tbl$gene_start[[i]], tbl$gene_start[[i + 1L]]) + 1L
      if (!is.na(overlap) && overlap > 0 && (
        (!is.na(len1) && len1 > 0 && overlap / len1 > 0.5) ||
          (!is.na(len2) && len2 > 0 && overlap / len2 > 0.5)
      )) {
        score_i <- c(
          ifelse(is.na(tbl$evalue[[i]]), Inf, tbl$evalue[[i]]),
          -ifelse(is.na(tbl$coverage[[i]]), -Inf, tbl$coverage[[i]])
        )
        score_j <- c(
          ifelse(is.na(tbl$evalue[[i + 1L]]), Inf, tbl$evalue[[i + 1L]]),
          -ifelse(is.na(tbl$coverage[[i + 1L]]), -Inf, tbl$coverage[[i + 1L]])
        )
        drop_index <- if (score_i[[1]] < score_j[[1]] ||
          (identical(score_i[[1]], score_j[[1]]) && score_i[[2]] <= score_j[[2]])) i + 1L else i
        tbl <- tbl[-drop_index, , drop = FALSE]
        if (i > 1L) {
          i <- i - 1L
        }
      } else {
        i <- i + 1L
      }
    }
    tbl
  })

  tibble::as_tibble(dplyr::bind_rows(filtered))
}

dnmb_dbcan_parse_domtblout <- function(path,
                                       evalue_threshold = 1e-15,
                                       coverage_threshold = 0.35) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(.dnmb_dbcan_empty_hits())
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  if (!length(lines)) {
    return(.dnmb_dbcan_empty_hits())
  }

  parsed <- lapply(lines, function(line) {
    fields <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(fields) < 19L) {
      return(NULL)
    }

    profile_id <- fields[[4]]
    profile_length <- suppressWarnings(as.integer(fields[[6]]))
    gene_id <- fields[[1]]
    gene_length <- suppressWarnings(as.integer(fields[[3]]))
    evalue <- suppressWarnings(as.numeric(fields[[13]]))
    profile_start <- suppressWarnings(as.integer(fields[[16]]))
    profile_end <- suppressWarnings(as.integer(fields[[17]]))
    gene_start <- suppressWarnings(as.integer(fields[[18]]))
    gene_end <- suppressWarnings(as.integer(fields[[19]]))
    coverage <- if (is.na(profile_length) || profile_length <= 0L || is.na(profile_start) || is.na(profile_end)) {
      NA_real_
    } else {
      (profile_end - profile_start + 1L) / profile_length
    }
    family_id <- .dnmb_dbcan_profile_family(profile_id)

    tibble::tibble(
      query = .dnmb_dbcan_clean_query_id(gene_id),
      profile_id = profile_id,
      profile_length = profile_length,
      gene_length = gene_length,
      evalue = evalue,
      profile_start = profile_start,
      profile_end = profile_end,
      gene_start = gene_start,
      gene_end = gene_end,
      coverage = coverage,
      family_id = family_id,
      hit_label = family_id,
      substrate_label = NA_character_
    )
  })

  hits <- dplyr::bind_rows(parsed)
  if (!nrow(hits)) {
    return(.dnmb_dbcan_empty_hits())
  }

  hits <- hits[!is.na(hits$evalue) & hits$evalue <= evalue_threshold & !is.na(hits$coverage) & hits$coverage >= coverage_threshold, , drop = FALSE]
  hits <- .dnmb_dbcan_filter_overlaps(hits)
  rownames(hits) <- NULL
  if (!nrow(hits)) {
    return(.dnmb_dbcan_empty_hits())
  }

  hits <- hits[, names(.dnmb_dbcan_empty_hits()), drop = FALSE]
  tibble::as_tibble(hits)
}

dnmb_dbcan_normalize_hits <- function(hits) {
  if (is.null(hits) || !nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  required <- c("query", "family_id", "hit_label", "evalue", "coverage", "profile_length", "gene_length", "profile_start", "profile_end", "gene_start", "gene_end")
  missing <- setdiff(required, names(hits))
  if (length(missing)) {
    stop("`hits` is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  support <- paste0(
    "evalue=", signif(hits$evalue, 4),
    "; coverage=", sprintf("%.3f", hits$coverage),
    "; profile_length=", hits$profile_length,
    "; gene_length=", hits$gene_length,
    "; profile_start=", hits$profile_start,
    "; profile_end=", hits$profile_end,
    "; gene_start=", hits$gene_start,
    "; gene_end=", hits$gene_end
  )

  out <- data.frame(
    query = .dnmb_dbcan_clean_query_id(hits$query),
    source = rep("dbcan", nrow(hits)),
    family_system = rep("dbCAN", nrow(hits)),
    family_id = as.character(hits$family_id),
    hit_label = as.character(hits$hit_label),
    enzyme_role = rep("CAZyme", nrow(hits)),
    evidence_mode = rep("direct", nrow(hits)),
    substrate_label = as.character(hits$substrate_label),
    support = support,
    profile_id = as.character(hits$profile_id),
    evalue = as.numeric(hits$evalue),
    coverage = as.numeric(hits$coverage),
    profile_length = as.integer(hits$profile_length),
    gene_length = as.integer(hits$gene_length),
    profile_start = as.integer(hits$profile_start),
    profile_end = as.integer(hits$profile_end),
    gene_start = as.integer(hits$gene_start),
    gene_end = as.integer(hits$gene_end),
    typing_eligible = rep(TRUE, nrow(hits)),
    stringsAsFactors = FALSE
  )
  mapping_columns <- intersect(
    c("dbcan_query_id", "dbcan_gene_row", "dbcan_mapping_source", "dbcan_contig_key"),
    names(hits)
  )
  for (column_name in mapping_columns) {
    out[[column_name]] <- hits[[column_name]]
  }
  out <- out[, c(.dnmb_module_optional_long_columns(), setdiff(names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
  hit_order <- order(out$query, hits$evalue, -hits$coverage)
  out <- out[hit_order, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.dnmb_dbcan_domain_summary <- function(hits) {
  if (is.null(hits) || !is.data.frame(hits) || !nrow(hits)) return(tibble::tibble())
  hits <- as.data.frame(hits, stringsAsFactors = FALSE)
  if ("source" %in% names(hits)) {
    hits <- hits[hits$source == "dbcan", , drop = FALSE]
  }
  required <- c("family_id", "gene_start", "gene_end", "evalue", "coverage")
  if (!nrow(hits) || !all(required %in% names(hits))) return(tibble::tibble())
  key <- if ("dbcan_gene_row" %in% names(hits) && any(!is.na(hits$dbcan_gene_row))) {
    paste0("row:", hits$dbcan_gene_row)
  } else {
    paste0("query:", .dnmb_dbcan_clean_query_id(hits$query))
  }
  groups <- split(seq_len(nrow(hits)), key)
  rows <- lapply(groups, function(idx) {
    tbl <- hits[idx, , drop = FALSE]
    tbl <- tbl[order(tbl$gene_start, tbl$gene_end, tbl$evalue, na.last = TRUE), , drop = FALSE]
    families <- as.character(tbl$family_id)
    families <- families[!is.na(families) & nzchar(families)]
    architecture <- if (length(families)) {
      paste0(
        as.character(tbl$family_id),
        "@",
        as.character(tbl$gene_start),
        "-",
        as.character(tbl$gene_end)
      )
    } else {
      character()
    }
    best_order <- order(tbl$evalue, -tbl$coverage, na.last = TRUE)
    best <- tbl[best_order[[1]], , drop = FALSE]
    out <- data.frame(
      query = as.character(best$query[[1]]),
      dbcan_gene_row = if ("dbcan_gene_row" %in% names(best)) suppressWarnings(as.integer(best$dbcan_gene_row[[1]])) else NA_integer_,
      dbcan_hmm_families = if (length(families)) paste(unique(families), collapse = "; ") else NA_character_,
      dbcan_hmm_domain_architecture = if (length(architecture)) paste(architecture, collapse = " | ") else NA_character_,
      dbcan_hmm_domain_count = length(families),
      stringsAsFactors = FALSE
    )
    best_columns <- intersect(
      c("profile_id", "evalue", "coverage", "profile_length", "gene_length", "profile_start", "profile_end", "gene_start", "gene_end", "support"),
      names(best)
    )
    for (column_name in best_columns) out[[column_name]] <- best[[column_name]][[1]]
    out
  })
  tibble::as_tibble(dplyr::bind_rows(rows))
}

.dnmb_dbcan_output_table <- function(genes, hits, cgc_genes = NULL,
                                      substrate = NULL, overview = NULL,
                                      id_map = NULL) {
  out <- .dnmb_module_output_table(
    genes = genes,
    hits = hits,
    key_normalizer = .dnmb_dbcan_clean_query_id
  )
  out$dbcan_gene_row <- seq_len(nrow(out))

  attach_rows <- function(target, source, exclude = character()) {
    if (is.null(source) || !is.data.frame(source) || !nrow(source)) return(target)
    source <- as.data.frame(source, stringsAsFactors = FALSE)
    idx <- rep(NA_integer_, nrow(target))
    if ("dbcan_gene_row" %in% names(source) && any(!is.na(source$dbcan_gene_row))) {
      idx <- match(target$dbcan_gene_row, suppressWarnings(as.integer(source$dbcan_gene_row)))
    }
    unresolved <- is.na(idx)
    if (any(unresolved) && "query" %in% names(source)) {
      query <- .dnmb_dbcan_clean_query_id(source$query)
      locus <- .dnmb_dbcan_clean_query_id(target$locus_tag)
      idx[unresolved] <- match(locus[unresolved], query)
      unresolved <- is.na(idx)
      if (any(unresolved) && "protein_id" %in% names(target)) {
        protein <- .dnmb_dbcan_clean_query_id(target$protein_id)
        idx[unresolved] <- match(protein[unresolved], query)
      }
    }
    include <- !is.na(idx)
    join_cols <- setdiff(names(source), c("query", "dbcan_gene_row", exclude))
    for (column_name in join_cols) {
      if (!column_name %in% names(target)) {
        target[[column_name]] <- .dnmb_na_vector_like(source[[column_name]], nrow(target))
      }
      target[[column_name]][include] <- source[[column_name]][idx[include]]
    }
    target
  }

  if (!is.null(id_map) && is.data.frame(id_map) && nrow(id_map)) {
    map_columns <- intersect(
      c("dbcan_gene_row", "dbcan_query_id", "mapping_source", "dbcan_contig_key"),
      names(id_map)
    )
    map <- id_map[, map_columns, drop = FALSE]
    if ("mapping_source" %in% names(map)) names(map)[names(map) == "mapping_source"] <- "dbcan_mapping_source"
    out <- attach_rows(out, map)
  }

  domain_summary <- .dnmb_dbcan_domain_summary(hits)
  out <- attach_rows(out, domain_summary)
  out <- attach_rows(out, overview)
  out <- attach_rows(out, cgc_genes)

  if (!is.null(substrate) && is.data.frame(substrate) && nrow(substrate) && "dbcan_cgc_id" %in% names(out)) {
    substrate <- as.data.frame(substrate, stringsAsFactors = FALSE)
    substrate <- substrate[!duplicated(substrate$dbcan_cgc_id), , drop = FALSE]
    match_rows <- match(out$dbcan_cgc_id, substrate$dbcan_cgc_id)
    include <- !is.na(match_rows)
    join_cols <- setdiff(names(substrate), "dbcan_cgc_id")
    for (column_name in join_cols) {
      if (!column_name %in% names(out)) {
        out[[column_name]] <- .dnmb_na_vector_like(substrate[[column_name]], nrow(out))
      }
      out[[column_name]][include] <- substrate[[column_name]][match_rows[include]]
    }
  }

  prefer_text <- function(...) {
    values <- list(...)
    result <- rep(NA_character_, nrow(out))
    for (value in values) {
      if (is.null(value)) next
      value <- as.character(value)
      take <- (is.na(result) | !nzchar(result)) & !is.na(value) & nzchar(value) & value != "-"
      result[take] <- value[take]
    }
    result
  }

  overview_primary <- if ("dbcan_primary_family" %in% names(out)) out$dbcan_primary_family else NULL
  hmm_primary <- if ("dbcan_hmm_families" %in% names(out)) {
    vapply(
      out$dbcan_hmm_families,
      function(value) .dnmb_dbcan_primary_family(.dnmb_dbcan_family_tokens(value)),
      character(1)
    )
  } else {
    NULL
  }
  out$family_id <- prefer_text(overview_primary, hmm_primary, out$family_id)
  all_families <- if ("dbcan_all_families" %in% names(out)) out$dbcan_all_families else NULL
  hmm_families <- if ("dbcan_hmm_families" %in% names(out)) out$dbcan_hmm_families else NULL
  out$dbcan_hit <- prefer_text(all_families, hmm_families, out$family_id)

  pul_substrate <- if ("dbcan_pul_substrate" %in% names(out)) out$dbcan_pul_substrate else rep(NA_character_, nrow(out))
  sub_substrate <- if ("dbcan_sub_substrate" %in% names(out)) out$dbcan_sub_substrate else rep(NA_character_, nrow(out))
  overview_substrate <- if ("dbcan_overview_substrate" %in% names(out)) out$dbcan_overview_substrate else rep(NA_character_, nrow(out))
  family_prior <- if ("dbcan_family_substrate_prior" %in% names(out)) out$dbcan_family_substrate_prior else rep(NA_character_, nrow(out))
  out$substrate_label <- prefer_text(pul_substrate, sub_substrate, overview_substrate, family_prior)
  out$dbcan_substrate_source <- ifelse(
    !is.na(pul_substrate) & nzchar(pul_substrate),
    "dbCAN-PUL",
    ifelse(
      !is.na(sub_substrate) & nzchar(sub_substrate),
      "dbCAN-sub",
      ifelse(
        !is.na(overview_substrate) & nzchar(overview_substrate),
        "overview",
        ifelse(!is.na(family_prior) & nzchar(family_prior), "family_prior", NA_character_)
      )
    )
  )

  drop_cols <- intersect(
    c("hit_label", "enzyme_role", "evidence_mode", "typing_eligible"),
    names(out)
  )
  if (length(drop_cols)) out[drop_cols] <- NULL

  base_cols <- intersect(dnmb_backbone_columns(), names(out))
  dbcan_cols <- c("dbcan_hit", "family_id", setdiff(names(out), c(base_cols, "dbcan_hit", "family_id")))
  out[, c(base_cols, dbcan_cols), drop = FALSE]
}

.dnmb_dbcan_standalone_layout <- function(module_dir) {
  file.path(module_dir, "standalone_bundle")
}

# Legacy URL list removed — S3 download is handled by `run_dbcan database`

.dnmb_dbcan_supporting_paths <- function(bundle_dir) {
  list(
    bundle_dir = bundle_dir,
    hmm_txt = file.path(bundle_dir, "dbCAN.hmm"),
    fam_substrate_mapping = file.path(bundle_dir, "fam-substrate-mapping.tsv"),
    pul_dmnd = file.path(bundle_dir, "PUL.dmnd"),
    pul_excel = file.path(bundle_dir, "dbCAN-PUL.xlsx"),
    pul_tar = file.path(bundle_dir, "dbCAN-PUL.tar.gz"),
    pul_dir = file.path(bundle_dir, "dbCAN-PUL"),
    tcdb_dmnd = file.path(bundle_dir, "TCDB.dmnd"),
    tf_hmm = file.path(bundle_dir, "TF.hmm"),
    tf_dmnd = file.path(bundle_dir, "TF.dmnd"),
    stp_hmm = file.path(bundle_dir, "STP.hmm"),
    dbcan_sub_hmm = file.path(bundle_dir, "dbCAN-sub.hmm"),
    dbcan_sub_hmm_real = file.path(bundle_dir, "dbCAN_sub.hmm"),
    cazy_dmnd = file.path(bundle_dir, "CAZy.dmnd"),
    peptidase_dmnd = file.path(bundle_dir, "peptidase_db.dmnd"),
    sulfatlas_dmnd = file.path(bundle_dir, "sulfatlas_db.dmnd")
  )
}

.dnmb_dbcan_standalone_ready <- function(paths) {
  required <- c(
    paths$hmm_txt,
    paste0(paths$hmm_txt, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paths$fam_substrate_mapping,
    paths$pul_dmnd,
    paths$pul_excel,
    paths$tcdb_dmnd,
    paths$tf_hmm,
    paths$tf_dmnd,
    paths$stp_hmm,
    paste0(paths$stp_hmm, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paste0(paths$tf_hmm, c(".h3f", ".h3i", ".h3m", ".h3p")),
    paths$cazy_dmnd,
    paths$peptidase_dmnd,
    paths$sulfatlas_dmnd,
    paths$dbcan_sub_hmm,
    paste0(paths$dbcan_sub_hmm, c(".h3f", ".h3i", ".h3m", ".h3p"))
  )
  all(file.exists(required)) && dir.exists(paths$pul_dir)
}

.dnmb_dbcan_stage_hmm_for_standalone <- function(module_dir, bundle_dir) {
  # Stage the module's dbCAN.hmm (stored as dbCAN.txt) into the bundle
  # as dbCAN.hmm so run_dbcan can find it.
  layout <- .dnmb_dbcan_asset_layout(module_dir)
  dest <- file.path(bundle_dir, "dbCAN.hmm")
  if (file.exists(dest)) {
    return(list(ok = TRUE, files = dest))
  }
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(layout$hmm_path, dest, overwrite = TRUE)
  # Also copy hmmpress indices if they exist
  for (ext in c(".h3f", ".h3i", ".h3m", ".h3p")) {
    src <- paste0(layout$hmm_path, ext)
    if (file.exists(src)) file.copy(src, paste0(dest, ext), overwrite = TRUE)
  }
  list(ok = ok, files = dest)
}

.dnmb_dbcan_prepare_reference_indexes <- function(paths) {

  # Create symlink dbCAN-sub.hmm -> dbCAN_sub.hmm if the S3-layout

  # uses underscore but run_dbcan expects hyphen
  if (!file.exists(paths$dbcan_sub_hmm) && file.exists(paths$dbcan_sub_hmm_real)) {
    suppressWarnings(file.symlink(
      normalizePath(paths$dbcan_sub_hmm_real, winslash = "/", mustWork = TRUE),
      paths$dbcan_sub_hmm
    ))
  }

  # hmmpress every .hmm that lacks its .h3f index
  hmm_files <- c(paths$hmm_txt, paths$tf_hmm, paths$stp_hmm, paths$dbcan_sub_hmm)
  for (hmm in hmm_files) {
    if (file.exists(hmm) && !file.exists(paste0(hmm, ".h3f"))) {
      dnmb_run_external("hmmpress", c("-f", hmm), required = FALSE)
    }
  }
}

.dnmb_dbcan_install_supporting_bundle <- function(module_dir,
                                                  trace_log,
                                                  base_url = .dnmb_dbcan_default_base_url(),
                                                  force = FALSE) {
  bundle_dir <- .dnmb_dbcan_standalone_layout(module_dir)
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- .dnmb_dbcan_supporting_paths(bundle_dir)
  status <- .dnmb_dbcan_empty_status()

  if (!isTRUE(force) && .dnmb_dbcan_standalone_ready(paths)) {
    return(list(ok = TRUE, status = .dnmb_dbcan_status_row("dbcan_support_bundle", "cached", bundle_dir), paths = paths))
  }

  # Delegate to `run_dbcan database --aws_s3 --cgc` for downloading
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] run_dbcan database --aws_s3 --cgc --db_dir %s", Sys.time(), bundle_dir))
  db_run <- dnmb_run_external("run_dbcan", c("database", "--aws_s3", "--cgc", "--db_dir", bundle_dir), required = FALSE)
  db_ok <- isTRUE(db_run$ok)
  status <- dplyr::bind_rows(
    status,
    .dnmb_dbcan_status_row(
      "dbcan_support_download",
      if (db_ok) "ok" else "failed",
      if (db_ok) bundle_dir else .dnmb_compact_output(c(db_run$stderr, db_run$stdout), max_lines = 6L)
    )
  )
  if (!db_ok) {
    return(list(ok = FALSE, status = status, paths = paths))
  }

  # Stage the module's dbCAN.hmm into the bundle if run_dbcan did not
  # provide it (the S3 download already includes dbCAN.hmm, but the
  # module may have a more recent pressed copy).
  if (!file.exists(paths$hmm_txt)) {
    hmm_stage <- .dnmb_dbcan_stage_hmm_for_standalone(module_dir = module_dir, bundle_dir = bundle_dir)
    status <- dplyr::bind_rows(
      status,
      .dnmb_dbcan_status_row(
        "dbcan_support_bundle:hmm_stage",
        if (isTRUE(hmm_stage$ok)) "ok" else "failed",
        paths$hmm_txt
      )
    )
    if (!isTRUE(hmm_stage$ok)) {
      return(list(ok = FALSE, status = status, paths = paths))
    }
  }

  # Symlink + hmmpress for HMM files
  .dnmb_dbcan_prepare_reference_indexes(paths)
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] dbcan_support_bundle=%s", Sys.time(), bundle_dir))

  ready <- .dnmb_dbcan_standalone_ready(paths)
  list(
    ok = ready,
    status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_support_bundle", if (ready) "ok" else "incomplete", bundle_dir)),
    paths = paths
  )
}

.dnmb_dbcan_can_run_standalone <- function(genes) {
  required <- c("start", "end", "direction")
  if (!is.data.frame(genes) || !all(required %in% names(genes))) return(FALSE)
  has_replicon <- any(c("contig_number", "record_accession", "contig_accession", "accession", "sequence_id", "contig") %in% names(genes))
  start <- suppressWarnings(as.numeric(genes$start))
  end <- suppressWarnings(as.numeric(genes$end))
  has_replicon && any(is.finite(start) & is.finite(end) & start <= end)
}

.dnmb_dbcan_safe_contig_map <- function(genes) {
  keys <- .dnmb_dbcan_contig_keys(as.data.frame(genes, stringsAsFactors = FALSE))
  contigs <- unique(keys$key)
  contigs <- contigs[!is.na(contigs) & nzchar(contigs)]
  tibble::tibble(
    contig = contigs,
    safe_contig = paste0("ctg", seq_along(contigs))
  )
}

.dnmb_dbcan_write_cluster_tsv <- function(genes, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- as.data.frame(genes, stringsAsFactors = FALSE)
  out$locus_tag <- .dnmb_dbcan_clean_query_id(out$locus_tag)
  out <- out[!is.na(out$locus_tag) & nzchar(out$locus_tag), c("contig", "locus_tag", "start", "end", "direction"), drop = FALSE]
  mapping <- .dnmb_dbcan_safe_contig_map(out)
  contig_key <- .dnmb_dbcan_contig_keys(out)$key
  out$contig <- mapping$safe_contig[match(contig_key, mapping$contig)]
  utils::write.table(out, file = path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  list(path = path, contig_map = mapping)
}

.dnmb_dbcan_write_query_gff <- function(genes, path, id_map = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(id_map)) id_map <- .dnmb_dbcan_build_gene_map(genes)
  out <- as.data.frame(id_map, stringsAsFactors = FALSE)
  out <- out[!is.na(out$dbcan_query_id) & nzchar(out$dbcan_query_id) &
             !is.na(out$dbcan_contig_key) & nzchar(out$dbcan_contig_key) &
             is.finite(out$start) & is.finite(out$end) & out$start <= out$end, , drop = FALSE]

  mapping <- tibble::tibble(
    contig = unique(out$dbcan_contig_key),
    safe_contig = paste0("ctg", seq_along(unique(out$dbcan_contig_key)))
  )
  out$safe_contig <- mapping$safe_contig[match(out$dbcan_contig_key, mapping$contig)]

  out$start <- as.numeric(out$start)
  out$end <- as.numeric(out$end)
  out <- out[order(out$safe_contig, out$start), ]

  strand <- ifelse(
    grepl("complement|minus|^-|^-1$", as.character(out$direction), ignore.case = TRUE),
    "-",
    "+"
  )
  attrs <- paste0("ID=", out$dbcan_query_id, ";Name=", out$dbcan_query_id)
  body <- paste(
    out$safe_contig,
    "DNMB",
    "CDS",
    as.integer(out$start),
    as.integer(out$end),
    ".",
    strand,
    "0",
    attrs,
    sep = "\t"
  )
  lines <- c("##gff-version 3", body)
  writeLines(lines, path)
  list(path = path, contig_map = mapping, id_map = id_map, n = nrow(out))
}

.dnmb_dbcan_restore_cgc_ids <- function(tbl, contig_map) {
  if (is.null(tbl) || !is.data.frame(tbl) || !nrow(tbl) || is.null(contig_map) || !nrow(contig_map) || !"dbcan_cgc_id" %in% names(tbl)) {
    return(tbl)
  }
  out <- as.data.frame(tbl, stringsAsFactors = FALSE)
  out$dbcan_cgc_id <- as.character(out$dbcan_cgc_id)
  safe_contig <- sub("\\|.*$", "", out$dbcan_cgc_id)
  suffix <- sub("^[^|]+", "", out$dbcan_cgc_id)
  idx <- match(safe_contig, contig_map$safe_contig)
  mapped <- !is.na(idx)
  out$dbcan_cgc_id[mapped] <- paste0(contig_map$contig[idx[mapped]], suffix[mapped])
  if ("dbcan_contig_key" %in% names(out)) {
    contig_idx <- match(out$dbcan_contig_key, contig_map$safe_contig)
    keep <- !is.na(contig_idx)
    out$dbcan_contig_key[keep] <- contig_map$contig[contig_idx[keep]]
  }
  out
}

.dnmb_dbcan_stage_dir_for_synteny <- function(source_dir) {
  stage_parent <- tempfile("dnmb-dbcan-synteny-db-")
  dir.create(stage_parent, recursive = TRUE, showWarnings = FALSE)
  stage_dir <- file.path(stage_parent, "db")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  source_pul_dir <- normalizePath(file.path(source_dir, "dbCAN-PUL"), winslash = "/", mustWork = TRUE)
  target_pul_dir <- file.path(stage_dir, "dbCAN-PUL")
  linked <- suppressWarnings(file.symlink(source_pul_dir, target_pul_dir))
  if (!isTRUE(linked)) {
    dir.create(target_pul_dir, recursive = TRUE, showWarnings = FALSE)
    copy_ok <- file.copy(list.files(source_pul_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE), target_pul_dir, recursive = TRUE)
    if (!all(copy_ok)) {
      stop("Failed to stage dbCAN-PUL database directory for synteny plotting.", call. = FALSE)
    }
  }
  list(
    db_dir = stage_dir,
    cleanup = function() unlink(stage_parent, recursive = TRUE, force = TRUE)
  )
}

.dnmb_dbcan_copy_dir_contents <- function(src_dir, dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(src_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(entries)) {
    return(invisible(dest_dir))
  }
  ok <- file.copy(entries, dest_dir, recursive = TRUE, overwrite = TRUE)
  if (!all(ok)) {
    stop("Failed to copy staged run_dbcan outputs.", call. = FALSE)
  }
  invisible(dest_dir)
}

.dnmb_dbcan_stage_for_run <- function(query_fasta, cluster_path, bundle_dir) {
  stage_parent <- tempfile("dnmb-dbcan-run-")
  dir.create(stage_parent, recursive = TRUE, showWarnings = FALSE)
  stage_query <- file.path(stage_parent, "query.faa")
  stage_cluster <- file.path(stage_parent, "cluster.tsv")
  stage_db <- file.path(stage_parent, "db")
  stage_out <- file.path(stage_parent, "run_dbcan")
  dir.create(stage_out, recursive = TRUE, showWarnings = FALSE)
  file.copy(query_fasta, stage_query, overwrite = TRUE)
  file.copy(cluster_path, stage_cluster, overwrite = TRUE)
  linked <- suppressWarnings(file.symlink(normalizePath(bundle_dir, winslash = "/", mustWork = TRUE), stage_db))
  if (!isTRUE(linked)) {
    .dnmb_dbcan_copy_dir_contents(bundle_dir, stage_db)
  }
  list(
    query_fasta = stage_query,
    cluster_tsv = stage_cluster,
    db_dir = stage_db,
    out_dir = stage_out,
    cleanup = function() unlink(stage_parent, recursive = TRUE, force = TRUE)
  )
}

.dnmb_dbcan_synteny_output_dir <- function(run_dbcan_dir) {
  candidates <- file.path(run_dbcan_dir, c("synteny_pdf", "synteny.pdf"))
  existing <- candidates[dir.exists(candidates)]
  if (length(existing)) existing[[1]] else candidates[[1]]
}

.dnmb_dbcan_run_syntenic_plot <- function(run_dbcan_dir, bundle_dir, trace_log = NULL) {
  required_files <- c(
    file.path(run_dbcan_dir, "PUL_blast.out"),
    file.path(run_dbcan_dir, "cgc_standard_out.tsv"),
    file.path(run_dbcan_dir, "substrate_prediction.tsv")
  )
  if (!all(file.exists(required_files))) {
    return(list(
      ok = FALSE,
      status = .dnmb_dbcan_status_row("dbcan_synteny", "missing_inputs", run_dbcan_dir),
      command = NULL,
      files = list()
    ))
  }

  stage <- .dnmb_dbcan_stage_dir_for_synteny(bundle_dir)
  on.exit(stage$cleanup(), add = TRUE)
  args <- c(
    "syntenic_plot",
    "-b", "PUL_blast.out",
    "--cgc", "cgc_standard_out.tsv",
    "-i", "substrate_prediction.tsv",
    "--db", stage$db_dir
  )
  if (!is.null(trace_log) && nzchar(trace_log)) {
    .dnmb_dbcan_trace(trace_log, sprintf("[%s] syntenic_plot args=%s", Sys.time(), .dnmb_format_command("syntenic_plot", args)))
  }
  command <- dnmb_run_external("syntenic_plot", args = args, wd = run_dbcan_dir, required = FALSE)
  synteny_dir <- .dnmb_dbcan_synteny_output_dir(run_dbcan_dir)
  ok <- isTRUE(command$ok) && dir.exists(synteny_dir)
  status <- .dnmb_dbcan_status_row(
    "dbcan_synteny",
    if (ok) "ok" else if (!nzchar(command$resolved_command)) "missing" else "failed",
    if (ok) synteny_dir else (command$error %||% synteny_dir)
  )

  list(
    ok = ok,
    status = status,
    command = command,
    files = if (dir.exists(synteny_dir)) list(synteny_dir = synteny_dir) else list()
  )
}

dnmb_dbcan_parse_hmmer_table <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(.dnmb_dbcan_empty_hits())
  }
  lines <- readLines(path, warn = FALSE)
  if (length(lines) <= 1L) {
    return(.dnmb_dbcan_empty_hits())
  }
  lines <- lines[-1L]
  parsed <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 10L) {
      return(NULL)
    }
    tibble::tibble(
      query = .dnmb_dbcan_clean_query_id(fields[[3]]),
      profile_id = fields[[1]],
      profile_length = suppressWarnings(as.integer(fields[[2]])),
      gene_length = suppressWarnings(as.integer(fields[[4]])),
      evalue = suppressWarnings(as.numeric(fields[[5]])),
      profile_start = suppressWarnings(as.integer(fields[[6]])),
      profile_end = suppressWarnings(as.integer(fields[[7]])),
      gene_start = suppressWarnings(as.integer(fields[[8]])),
      gene_end = suppressWarnings(as.integer(fields[[9]])),
      coverage = suppressWarnings(as.numeric(fields[[10]])),
      family_id = .dnmb_dbcan_profile_family(fields[[1]]),
      hit_label = .dnmb_dbcan_profile_family(fields[[1]]),
      substrate_label = NA_character_
    )
  })
  hits <- dplyr::bind_rows(parsed)
  if (!nrow(hits)) {
    return(.dnmb_dbcan_empty_hits())
  }
  hits <- hits[, names(.dnmb_dbcan_empty_hits()), drop = FALSE]
  tibble::as_tibble(hits)
}

dnmb_dbcan_parse_cgc_standard <- function(path) {
  empty <- tibble::tibble(
    query = character(), dbcan_cgc_id = character(),
    dbcan_contig_key = character(), dbcan_cgc_gene_type = character(),
    dbcan_cgc_protein_family = character()
  )
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) return(empty)

  tbl <- tryCatch(
    utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                      fill = TRUE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tbl) || !nrow(tbl)) return(empty)

  # Map column names to DNMB's internal names
  col_map <- c(
    "Protein ID" = "query",
    "CGC#" = "cgc_num",
    "Contig ID" = "contig_id",
    "Gene Type" = "dbcan_cgc_gene_type",
    "Gene Annotation" = "dbcan_cgc_protein_family"
  )
  for (old_name in names(col_map)) {
    if (old_name %in% names(tbl)) {
      names(tbl)[names(tbl) == old_name] <- col_map[old_name]
    }
  }

  # Build dbcan_cgc_id from contig + CGC#
  if ("contig_id" %in% names(tbl) && "cgc_num" %in% names(tbl)) {
    tbl$dbcan_cgc_id <- paste(tbl$contig_id, tbl$cgc_num, sep = "|")
    tbl$dbcan_contig_key <- as.character(tbl$contig_id)
  } else {
    tbl$dbcan_cgc_id <- NA_character_
    tbl$dbcan_contig_key <- NA_character_
  }

  # Clean query
  if ("query" %in% names(tbl)) {
    tbl$query <- .dnmb_dbcan_clean_query_id(tbl$query)
  }

  result <- tbl[, intersect(c("query", "dbcan_cgc_id", "dbcan_contig_key", "dbcan_cgc_gene_type", "dbcan_cgc_protein_family"), names(tbl)), drop = FALSE]
  tibble::as_tibble(result)
}

dnmb_dbcan_parse_substrate_table <- function(path) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    return(tibble::tibble(
      dbcan_cgc_id = character(),
      dbcan_pul_id = character(),
      dbcan_pul_substrate = character(),
      dbcan_pul_bitscore = numeric(),
      dbcan_pul_signature_pairs = character(),
      dbcan_sub_substrate = character(),
      dbcan_sub_score = character()
    ))
  }

  tbl <- utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    fill = TRUE,
    check.names = FALSE
  )
  if (!nrow(tbl)) {
    return(tibble::tibble(
      dbcan_cgc_id = character(),
      dbcan_pul_id = character(),
      dbcan_pul_substrate = character(),
      dbcan_pul_bitscore = numeric(),
      dbcan_pul_signature_pairs = character(),
      dbcan_sub_substrate = character(),
      dbcan_sub_score = character()
    ))
  }
  names(tbl)[seq_len(min(7L, ncol(tbl)))] <- c(
    "dbcan_cgc_id",
    "dbcan_pul_id",
    "dbcan_pul_substrate",
    "dbcan_pul_bitscore",
    "dbcan_pul_signature_pairs",
    "dbcan_sub_substrate",
    "dbcan_sub_score"
  )[seq_len(min(7L, ncol(tbl)))]
  out <- tibble::as_tibble(tbl)
  out$dbcan_pul_bitscore <- suppressWarnings(as.numeric(out$dbcan_pul_bitscore))
  out[, c("dbcan_cgc_id", "dbcan_pul_id", "dbcan_pul_substrate", "dbcan_pul_bitscore", "dbcan_pul_signature_pairs", "dbcan_sub_substrate", "dbcan_sub_score")]
}

.dnmb_dbcan_augment_hits_with_standalone <- function(hits, cgc_genes = NULL, substrate = NULL) {
  hits <- as.data.frame(hits, stringsAsFactors = FALSE)
  if (!is.null(cgc_genes) && is.data.frame(cgc_genes) && nrow(cgc_genes)) {
    hits <- dplyr::left_join(hits, cgc_genes, by = c("query" = "query"))
  }
  if (!is.null(substrate) && is.data.frame(substrate) && nrow(substrate) && "dbcan_cgc_id" %in% names(hits)) {
    hits <- dplyr::left_join(hits, substrate, by = "dbcan_cgc_id")
  }
  hits
}

dnmb_dbcan_run_standalone <- function(query_fasta,
                                      genes,
                                      output_dir,
                                      version = .dnmb_dbcan_default_version(),
                                      cache_root = NULL,
                                      cpu = 1L,
                                      install = FALSE,
                                      base_url = .dnmb_dbcan_default_base_url(),
                                      id_map = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_dbcan_empty_status()
  trace_log <- file.path(output_dir, "dbcan_trace.log")
  module_dir <- .dnmb_db_module_dir(.dnmb_dbcan_module_name(), version, cache_root = cache_root, create = TRUE)
  bundle_dir <- .dnmb_dbcan_standalone_layout(module_dir)

  module <- dnmb_dbcan_get_module(version = version, cache_root = cache_root, required = FALSE)
  if (!isTRUE(module$ok) && isTRUE(install)) {
    install_result <- dnmb_dbcan_install_module(version = version, cache_root = cache_root, base_url = base_url, force = FALSE, prepare = TRUE)
    status <- dplyr::bind_rows(status, install_result$status)
    module <- dnmb_dbcan_get_module(version = version, cache_root = cache_root, required = FALSE)
  }
  if (!isTRUE(module$ok)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_module", "missing", paste0("dbCAN module not installed for version ", version))),
      files = list(query_fasta = query_fasta),
      command = NULL,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      cgc_genes = tibble::tibble(),
      substrate = tibble::tibble(),
      manifest = module$manifest
    ))
  }

  support_bundle <- .dnmb_dbcan_install_supporting_bundle(
    module_dir = module_dir,
    trace_log = trace_log,
    base_url = base_url,
    force = FALSE
  )
  status <- dplyr::bind_rows(status, support_bundle$status)
  if (!isTRUE(support_bundle$ok)) {
    return(list(
      ok = FALSE,
      status = status,
      files = list(query_fasta = query_fasta, trace_log = trace_log),
      command = NULL,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      cgc_genes = tibble::tibble(),
      substrate = tibble::tibble(),
      manifest = module$manifest
    ))
  }

  if (is.null(id_map)) id_map <- .dnmb_dbcan_build_gene_map(genes)
  gff_info <- .dnmb_dbcan_write_query_gff(
    genes,
    file.path(output_dir, "dbcan_query.gff"),
    id_map = id_map
  )

  run_dbcan_output <- file.path(output_dir, "run_dbcan")
  if (dir.exists(run_dbcan_output)) {
    unlink(run_dbcan_output, recursive = TRUE, force = TRUE)
  }
  dir.create(run_dbcan_output, recursive = TRUE, showWarnings = FALSE)

  # Stage: symlink db_dir for run_dbcan
  stage_parent <- tempfile("dnmb-dbcan-run-")
  dir.create(stage_parent, recursive = TRUE, showWarnings = FALSE)
  stage_db <- file.path(stage_parent, "db")
  linked <- suppressWarnings(file.symlink(normalizePath(bundle_dir, winslash = "/", mustWork = TRUE), stage_db))
  if (!isTRUE(linked)) {
    .dnmb_dbcan_copy_dir_contents(bundle_dir, stage_db)
  }
  on.exit(unlink(stage_parent, recursive = TRUE, force = TRUE), add = TRUE)

  # Call run_dbcan easy_substrate (new 5.x subcommand CLI)
  args <- c(
    "easy_substrate",
    "--db_dir", stage_db,
    "--mode", "protein",
    "--output_dir", run_dbcan_output,
    "--input_raw_data", query_fasta,
    "--input_gff", gff_info$path,
    "--additional_genes", "TC",
    "--additional_genes", "TF",
    "--additional_genes", "STP",
    "--additional_logic", "any",
    "--num_null_gene", "2",
    "--use_null_genes",
    "--threads", as.character(max(1L, as.integer(cpu)[1]))
  )
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] run_dbcan args=%s", Sys.time(), .dnmb_format_command("run_dbcan", args)))
  command <- dnmb_run_external("run_dbcan", args = args, required = FALSE)

  hmmer_path <- file.path(run_dbcan_output, "dbCAN_hmm_results.tsv")
  overview_path <- file.path(run_dbcan_output, "overview.tsv")
  cgc_standard_path <- file.path(run_dbcan_output, "cgc_standard_out.tsv")
  substrate_path <- file.path(run_dbcan_output, "substrate_prediction.tsv")

  raw_hits <- .dnmb_dbcan_restore_query_map(
    dnmb_dbcan_parse_hmmer_table(hmmer_path),
    id_map
  )
  overview <- dnmb_dbcan_parse_overview(overview_path, id_map = id_map)
  cgc_genes <- .dnmb_dbcan_restore_cgc_ids(
    dnmb_dbcan_parse_cgc_standard(cgc_standard_path),
    contig_map = gff_info$contig_map
  )
  cgc_genes <- .dnmb_dbcan_restore_query_map(cgc_genes, id_map)
  substrate <- .dnmb_dbcan_restore_cgc_ids(
    dnmb_dbcan_parse_substrate_table(substrate_path),
    contig_map = gff_info$contig_map
  )
  hits <- dnmb_dbcan_normalize_hits(raw_hits)
  mapped_domains_path <- file.path(run_dbcan_output, "dnmb_mapped_domains.tsv")
  gene_summary_path <- file.path(run_dbcan_output, "dnmb_cazy_gene_summary.tsv")
  utils::write.table(
    as.data.frame(hits, stringsAsFactors = FALSE),
    file = mapped_domains_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  utils::write.table(
    as.data.frame(overview, stringsAsFactors = FALSE),
    file = gene_summary_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  synteny <- .dnmb_dbcan_run_syntenic_plot(
    run_dbcan_dir = run_dbcan_output,
    bundle_dir = bundle_dir,
    trace_log = trace_log
  )

  ok <- isTRUE(command$ok) && file.exists(overview_path)
  status <- dplyr::bind_rows(
    status,
    .dnmb_dbcan_status_row("run_dbcan", if (ok) "ok" else if (isTRUE(command$ok)) "partial" else "failed", run_dbcan_output),
    .dnmb_dbcan_status_row("dbcan_overview", if (file.exists(overview_path)) "ok" else "missing", overview_path),
    .dnmb_dbcan_status_row("dbcan_cgc", if (file.exists(cgc_standard_path)) "ok" else "missing", cgc_standard_path),
    .dnmb_dbcan_status_row("dbcan_substrate", if (file.exists(substrate_path)) "ok" else "missing", substrate_path),
    synteny$status
  )

  list(
    ok = ok,
    status = status,
    files = c(
      list(
        query_fasta = query_fasta,
        query_gff = gff_info$path,
        run_dbcan_dir = run_dbcan_output,
        hmmer_out = hmmer_path,
        overview = overview_path,
        cgc_standard = cgc_standard_path,
        substrate = substrate_path,
        mapped_domains = mapped_domains_path,
        gene_summary = gene_summary_path,
        trace_log = trace_log
      ),
      synteny$files,
      module$files
    ),
    command = command,
    raw_hits = raw_hits,
    hits = hits,
    overview = overview,
    cgc_genes = cgc_genes,
    substrate = substrate,
    manifest = module$manifest
  )
}

dnmb_dbcan_install_module <- function(version = .dnmb_dbcan_default_version(),
                                      cache_root = NULL,
                                      base_url = .dnmb_dbcan_default_base_url(),
                                      asset_urls = NULL,
                                      force = FALSE,
                                      prepare = TRUE) {
  module <- .dnmb_dbcan_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_dbcan_asset_layout(module_dir)
  install_trace_log <- file.path(module_dir, "dbcan_install_trace.log")
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  release_info <- .dnmb_dbcan_release_info()
  urls <- .dnmb_dbcan_asset_urls(base_url = base_url, version = version, asset_urls = asset_urls, release_info = release_info)
  status <- .dnmb_dbcan_empty_status()
  remote_asset_state <- .dnmb_dbcan_remote_asset_state(
    urls = urls,
    enabled = identical(trimws(as.character(version)[1]), "current") && is.null(asset_urls)
  )
  remote_refresh_needed <- .dnmb_dbcan_remote_update_needed(manifest, remote_asset_state)
  if (!is.null(remote_asset_state) && nrow(remote_asset_state)) {
    remote_detail <- paste(
      sprintf(
        "%s(last_modified=%s; content_length=%s)",
        remote_asset_state$asset_name,
        ifelse(is.na(remote_asset_state$last_modified), "NA", remote_asset_state$last_modified),
        ifelse(is.na(remote_asset_state$content_length), "NA", as.character(remote_asset_state$content_length))
      ),
      collapse = "; "
    )
    status <- dplyr::bind_rows(
      status,
      .dnmb_dbcan_status_row(
        "dbcan_remote_state",
        if (all(remote_asset_state$ok %in% TRUE)) {
          if (remote_refresh_needed) "update_detected" else "ok"
        } else {
          "unavailable"
        },
        remote_detail
      )
    )
    .dnmb_dbcan_trace(install_trace_log, sprintf("[%s] remote_asset_state %s", Sys.time(), remote_detail))
  }
  if (isTRUE(remote_refresh_needed) && !isTRUE(force)) {
    force <- TRUE
    .dnmb_dbcan_trace(install_trace_log, sprintf("[%s] forcing refresh because current dbCAN metadata changed", Sys.time()))
  }

  required_assets <- c("dbCAN.txt")
  required_ready <- file.exists(layout$hmm_path) && .dnmb_dbcan_is_hmm(layout$hmm_path)
  manifest_ready <- !is.null(manifest)
  if (manifest_ready && required_ready && !isTRUE(force)) {
    prepare_result <- if (isTRUE(prepare)) {
      .dnmb_dbcan_prepare_hmm_db(layout$hmm_path, force = FALSE, trace_log = install_trace_log)
    } else {
      NULL
    }
    if (!is.null(prepare_result)) {
      status <- dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_prepare", prepare_result$status, prepare_result$detail))
    }
    cached_ok <- if (isTRUE(prepare) && !is.null(prepare_result)) {
      isTRUE(prepare_result$ok)
    } else {
      isTRUE(!is.null(manifest$install_ok) && manifest$install_ok) || .dnmb_dbcan_is_hmm(layout$hmm_path)
    }
    asset_identity <- .dnmb_db_cached_file_md5(
      c(layout$hmm_path, layout$mapping_path),
      manifest = manifest
    )
    manifest_fields <- unclass(manifest)
    manifest_fields[c("module", "version", "module_dir", "manifest_path", "written_at")] <- NULL
    manifest_fields <- utils::modifyList(
      manifest_fields,
      c(
        .dnmb_dbcan_manifest_diagnostics(layout, prepare_result = prepare_result, release_info = release_info),
        list(
          hmmpress_files = layout$hmmpress_files[file.exists(layout$hmmpress_files)],
          prepared_with_hmmpress = isTRUE(!is.null(prepare_result) && prepare_result$status %in% c("ok", "cached") && length(prepare_result$hmmpress_files) == 4L),
          asset_state = asset_identity$asset_state,
          asset_md5 = asset_identity$asset_md5,
          remote_asset_state = remote_asset_state,
          install_ok = cached_ok
        )
      )
    )
    dnmb_db_write_manifest(
      module = module,
      version = version,
      cache_root = cache_root,
      manifest = manifest_fields,
      overwrite = TRUE
    )
    .dnmb_db_autoprune_default_versions(
      module = module,
      version = version,
      default_version = .dnmb_dbcan_default_version(),
      cache_root = cache_root
    )
    manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
    return(list(
      ok = cached_ok,
      status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_install", if (cached_ok) "cached" else "failed", module_dir)),
      manifest = manifest,
      files = c(
        list(hmm_path = layout$hmm_path, trace_log = install_trace_log),
        if (file.exists(layout$mapping_path) && .dnmb_dbcan_is_mapping(layout$mapping_path)) list(mapping_path = layout$mapping_path),
        stats::setNames(as.list(layout$hmmpress_files[file.exists(layout$hmmpress_files)]), basename(layout$hmmpress_files[file.exists(layout$hmmpress_files)]))
      )
    ))
  }

  download_results <- list()
  for (asset_name in names(urls)) {
    asset_path <- switch(
      asset_name,
      "dbCAN.txt" = layout$hmm_path,
      "fam-substrate-mapping.tsv" = layout$mapping_path,
      file.path(module_dir, asset_name)
    )
    result <- .dnmb_download_asset(urls[[asset_name]], asset_path, insecure = TRUE)
    download_results[[asset_name]] <- result
    required_asset <- asset_name %in% required_assets
    asset_ok <- if (identical(asset_name, "dbCAN.txt")) {
      .dnmb_dbcan_is_hmm(asset_path)
    } else if (identical(asset_name, "fam-substrate-mapping.tsv")) {
      .dnmb_dbcan_is_mapping(asset_path)
    } else {
      isTRUE(result$ok) || file.exists(asset_path)
    }
    status <- dplyr::bind_rows(
      status,
      .dnmb_dbcan_status_row(
        paste0("dbcan_download:", asset_name),
        if (asset_ok) "ok" else if (required_asset) "failed" else "optional_missing",
        if (asset_ok) asset_path else (result$error %||% asset_name)
      )
    )
  }

  ok <- .dnmb_dbcan_is_hmm(layout$hmm_path)
  prepare_result <- NULL
  if (ok && isTRUE(prepare)) {
    prepare_result <- .dnmb_dbcan_prepare_hmm_db(layout$hmm_path, force = isTRUE(force), trace_log = install_trace_log)
    status <- dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_prepare", prepare_result$status, prepare_result$detail))
    ok <- isTRUE(prepare_result$ok)
  }

  mapping_available <- file.exists(layout$mapping_path) && .dnmb_dbcan_is_mapping(layout$mapping_path)
  asset_identity <- .dnmb_db_cached_file_md5(
    c(layout$hmm_path, layout$mapping_path),
    manifest = manifest
  )
  manifest_payload <- c(
    list(
      source = "dbcan",
      base_url = base_url,
      asset_urls = urls,
      asset_files = names(urls),
      required_assets = required_assets,
      hmm_path = layout$hmm_path,
      mapping_path = if (mapping_available) layout$mapping_path else NA_character_,
      hmmpress_files = layout$hmmpress_files[file.exists(layout$hmmpress_files)],
      prepared_with_hmmpress = isTRUE(!is.null(prepare_result) && prepare_result$status %in% c("ok", "cached") && length(prepare_result$hmmpress_files) == 4L),
      mapping_available = mapping_available,
      asset_state = asset_identity$asset_state,
      asset_md5 = asset_identity$asset_md5,
      remote_asset_state = remote_asset_state
    ),
    .dnmb_dbcan_manifest_diagnostics(layout, prepare_result = prepare_result, release_info = release_info),
    list(
      download_methods = vapply(download_results, function(x) x$method %||% NA_character_, character(1)),
      install_ok = ok
    )
  )
  manifest_path <- dnmb_db_write_manifest(
    module = module,
    version = version,
    cache_root = cache_root,
    manifest = manifest_payload,
    overwrite = TRUE
  )
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_dbcan_default_version(),
    cache_root = cache_root
  )
  manifest <- readRDS(manifest_path)
  class(manifest) <- c("dnmb_db_manifest", class(manifest))

  list(
    ok = ok,
    status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_install", if (ok) "ok" else "failed", module_dir)),
    manifest = manifest,
    files = c(
      list(hmm_path = layout$hmm_path, trace_log = install_trace_log),
      if (mapping_available) list(mapping_path = layout$mapping_path),
      stats::setNames(as.list(layout$hmmpress_files[file.exists(layout$hmmpress_files)]), basename(layout$hmmpress_files[file.exists(layout$hmmpress_files)]))
    )
  )
}

dnmb_dbcan_get_module <- function(version = .dnmb_dbcan_default_version(),
                                  cache_root = NULL,
                                  required = FALSE) {
  module <- .dnmb_dbcan_module_name()
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = FALSE)
  layout <- .dnmb_dbcan_asset_layout(module_dir)
  ok <- !is.null(manifest) && .dnmb_dbcan_is_hmm(layout$hmm_path)

  if (!ok && isTRUE(required)) {
    stop("dbCAN module is not installed for version `", version, "`.", call. = FALSE)
  }

  list(
    ok = ok,
    module = module,
    version = version,
    module_dir = module_dir,
    manifest = manifest,
    files = c(
      list(hmm_path = layout$hmm_path),
      if (file.exists(layout$mapping_path) && .dnmb_dbcan_is_mapping(layout$mapping_path)) list(mapping_path = layout$mapping_path),
      stats::setNames(as.list(layout$hmmpress_files[file.exists(layout$hmmpress_files)]), basename(layout$hmmpress_files[file.exists(layout$hmmpress_files)]))
    )
  )
}

dnmb_dbcan_run_hmmsearch <- function(query_fasta,
                                     output_dir,
                                     version = .dnmb_dbcan_default_version(),
                                     cache_root = NULL,
                                     cpu = 1L,
                                     install = FALSE,
                                     base_url = .dnmb_dbcan_default_base_url(),
                                     asset_urls = NULL,
                                     evalue_threshold = 1e-15,
                                     coverage_threshold = 0.35,
                                     id_map = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_dbcan_empty_status()
  trace_log <- file.path(output_dir, "dbcan_trace.log")
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] enter dnmb_dbcan_run_hmmsearch query=%s", Sys.time(), query_fasta))

  module <- dnmb_dbcan_get_module(version = version, cache_root = cache_root, required = FALSE)
  if (!isTRUE(module$ok) && isTRUE(install)) {
    install_result <- dnmb_dbcan_install_module(
      version = version,
      cache_root = cache_root,
      base_url = base_url,
      asset_urls = asset_urls,
      force = FALSE,
      prepare = TRUE
    )
    status <- dplyr::bind_rows(status, install_result$status)
    module <- dnmb_dbcan_get_module(version = version, cache_root = cache_root, required = FALSE)
  }

  if (!isTRUE(module$ok)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_module", "missing", paste0("dbCAN module not installed for version ", version))),
      files = list(query_fasta = query_fasta),
      command = NULL,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      manifest = module$manifest
    ))
  }

  prepare_result <- .dnmb_dbcan_prepare_hmm_db(module$files$hmm_path, force = FALSE, trace_log = trace_log)
  status <- dplyr::bind_rows(status, .dnmb_dbcan_status_row("dbcan_prepare", prepare_result$status, prepare_result$detail))
  if (!isTRUE(prepare_result$ok)) {
    return(list(
      ok = FALSE,
      status = status,
      files = c(list(query_fasta = query_fasta), module$files),
      command = prepare_result$command,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      manifest = module$manifest
    ))
  }

  domtblout <- file.path(output_dir, "dbcan_hmmsearch.domtblout")
  stdout_log <- file.path(output_dir, "dbcan_hmmsearch.stdout.log")
  stderr_log <- file.path(output_dir, "dbcan_hmmsearch.stderr.log")
  hmm_stdout <- file.path(output_dir, "dbcan_hmmsearch.out")
  args <- c(
    "--domtblout", domtblout,
    "--cpu", as.character(as.integer(cpu)),
    "-o", hmm_stdout,
    module$files$hmm_path,
    query_fasta
  )
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] hmmsearch args=%s", Sys.time(), .dnmb_format_command("hmmsearch", args)))
  command <- dnmb_run_external("hmmsearch", args = args, required = FALSE)
  writeLines(command$stdout, con = stdout_log)
  writeLines(command$stderr, con = stderr_log)

  raw_hits <- dnmb_dbcan_parse_domtblout(
    domtblout,
    evalue_threshold = evalue_threshold,
    coverage_threshold = coverage_threshold
  )
  raw_hits <- .dnmb_dbcan_restore_query_map(raw_hits, id_map)
  hits <- dnmb_dbcan_normalize_hits(raw_hits)
  status <- dplyr::bind_rows(
    status,
    .dnmb_dbcan_status_row(
      "hmmsearch",
      if (isTRUE(command$ok)) if (nrow(raw_hits)) "ok" else "empty" else if (!nzchar(command$resolved_command)) "missing" else "failed",
      if (isTRUE(command$ok)) domtblout else (command$error %||% "hmmsearch failed")
    )
  )

  list(
    ok = isTRUE(command$ok),
    status = status,
    files = c(
      list(
        query_fasta = query_fasta,
        domtblout = domtblout,
        hmm_stdout = hmm_stdout,
        stdout_log = stdout_log,
        stderr_log = stderr_log,
        trace_log = trace_log
      ),
      module$files
    ),
    command = command,
    raw_hits = raw_hits,
    hits = hits,
    manifest = module$manifest
  )
}

dnmb_run_dbcan_module <- function(genes,
                                  output_dir,
                                  version = .dnmb_dbcan_default_version(),
                                  cache_root = NULL,
                                  install = TRUE,
                                  base_url = .dnmb_dbcan_default_base_url(),
                                  asset_urls = NULL,
                                  cpu = 1L,
                                  genbank = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  standalone_output <- file.path(output_dir, "run_dbcan")
  if (dir.exists(standalone_output)) {
    unlink(standalone_output, recursive = TRUE, force = TRUE)
  }
  trace_log <- file.path(output_dir, "dbcan_module_trace.log")
  fasta_path <- file.path(output_dir, "dbcan_query_proteins.faa")
  id_map <- .dnmb_dbcan_build_gene_map(genes)
  id_map_path <- file.path(output_dir, "dbcan_gene_id_map.tsv")
  .dnmb_dbcan_write_gene_map(id_map, id_map_path)
  fasta <- .dnmb_dbcan_write_mapped_query_fasta(id_map, fasta_path)
  status <- .dnmb_dbcan_status_row(
    "dbcan_query_fasta",
    if (fasta$n) "ok" else "empty",
    paste0("proteins=", fasta$n)
  )
  .dnmb_dbcan_trace(trace_log, sprintf("[%s] query_faa=%s synthetic_ids=true proteins=%s", Sys.time(), fasta_path, fasta$n))

  if (!fasta$n) {
    return(list(
      ok = TRUE,
      status = status,
      files = list(query_fasta = fasta_path, gene_id_map = id_map_path, trace_log = trace_log),
      manifest = NULL,
      raw_hits = .dnmb_dbcan_empty_hits(),
      hits = .dnmb_module_empty_optional_long_table(),
      overview = .dnmb_dbcan_empty_overview(),
      output_table = .dnmb_dbcan_output_table(
        genes = genes,
        hits = .dnmb_module_empty_optional_long_table(),
        id_map = id_map
      ),
      query_proteins = fasta$proteins
    ))
  }

  standalone_inputs_ok <- .dnmb_dbcan_can_run_standalone(genes)
  run_dbcan_detection <- dnmb_detect_binary("run_dbcan", required = FALSE)
  standalone_possible <- standalone_inputs_ok && isTRUE(run_dbcan_detection$found)
  if (standalone_possible) {
    search <- dnmb_dbcan_run_standalone(
      query_fasta = fasta_path,
      genes = genes,
      output_dir = output_dir,
      version = version,
      cache_root = cache_root,
      cpu = cpu,
      install = install,
      base_url = base_url,
      id_map = id_map
    )
    if (!isTRUE(search$ok)) {
      standalone_status <- search$status
      if (dir.exists(standalone_output)) {
        unlink(standalone_output, recursive = TRUE, force = TRUE)
      }
      fallback <- dnmb_dbcan_run_hmmsearch(
        query_fasta = fasta_path,
        output_dir = output_dir,
        version = version,
        cache_root = cache_root,
        cpu = cpu,
        install = install,
        base_url = base_url,
        asset_urls = asset_urls,
        id_map = id_map
      )
      fallback$status <- dplyr::bind_rows(
        standalone_status,
        .dnmb_dbcan_status_row(
          "run_dbcan_fallback",
          if (isTRUE(fallback$ok)) "ok" else "failed",
          "standalone dbCAN failed; retried with HMMER"
        ),
        fallback$status
      )
      search <- fallback
    }
  } else {
    missing_cols <- setdiff(c("start", "end", "direction"), names(genes))
    skip_detail <- if (!standalone_inputs_ok) {
      paste0("missing required GenBank columns for CGC mode: ", paste(missing_cols, collapse = ", "))
    } else {
      run_dbcan_detection$message
    }
    search <- dnmb_dbcan_run_hmmsearch(
      query_fasta = fasta_path,
      output_dir = output_dir,
      version = version,
      cache_root = cache_root,
      cpu = cpu,
      install = install,
      base_url = base_url,
      asset_urls = asset_urls,
      id_map = id_map
    )
    search$status <- dplyr::bind_rows(
      .dnmb_dbcan_status_row("run_dbcan", "skipped", skip_detail),
      search$status
    )
  }

  mapping_path <- search$files$mapping_path %||% NULL
  if (!is.null(mapping_path) && length(mapping_path)) {
    search$overview <- .dnmb_dbcan_add_family_substrate_prior(
      search$overview %||% .dnmb_dbcan_empty_overview(),
      as.character(mapping_path[[1]])
    )
    gene_summary_path <- file.path(output_dir, "run_dbcan", "dnmb_cazy_gene_summary.tsv")
    if (dir.exists(dirname(gene_summary_path))) {
      utils::write.table(
        as.data.frame(search$overview, stringsAsFactors = FALSE),
        file = gene_summary_path,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE,
        na = ""
      )
    }
  }
  search$hits <- .dnmb_dbcan_merge_overview_hits(
    search$hits,
    search$overview %||% .dnmb_dbcan_empty_overview()
  )

  list(
    ok = isTRUE(search$ok),
    status = dplyr::bind_rows(status, search$status),
    files = c(search$files, list(gene_id_map = id_map_path, module_trace_log = trace_log)),
    manifest = search$manifest,
    raw_hits = search$raw_hits,
    hits = search$hits,
    overview = search$overview %||% .dnmb_dbcan_empty_overview(),
    cgc_genes = search$cgc_genes %||% tibble::tibble(),
    substrate = search$substrate %||% tibble::tibble(),
    output_table = .dnmb_dbcan_output_table(
      genes = genes,
      hits = search$hits,
      cgc_genes = search$cgc_genes %||% NULL,
      substrate = search$substrate %||% NULL,
      overview = search$overview %||% NULL,
      id_map = id_map
    ),
    query_proteins = fasta$proteins
  )
}

.dnmb_dbapis_module_name <- function() {
  "dbapis"
}

.dnmb_dbapis_default_version <- function() {
  "current"
}

.dnmb_dbapis_default_repo_url <- function() {
  "https://github.com/azureycy/dbAPIS.git"
}

.dnmb_dbapis_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_dbapis_empty_status <- function() {
  .dnmb_dbapis_status_row(character(), character(), character())
}

.dnmb_dbapis_asset_layout <- function(module_dir) {
  data_dir <- base::file.path(module_dir, "data_download")
  list(
    module_dir = module_dir,
    repo_dir = base::file.path(module_dir, "dbAPIS"),
    data_dir = data_dir,
    hmm_zip_path = base::file.path(data_dir, "dbAPIS.hmm.zip"),
    hmm_path = base::file.path(data_dir, "dbAPIS.hmm"),
    mapping_path = base::file.path(data_dir, "seed_family_mapping.tsv"),
    info_path = base::file.path(data_dir, "seed_and_familyrep_all_infor.tsv"),
    fasta_path = base::file.path(data_dir, "anti_defense.pep"),
    readme_path = base::file.path(data_dir, "readme.txt")
  )
}

.dnmb_dbapis_normalize_asset_urls <- function(asset_urls = NULL) {
  if (base::is.null(asset_urls)) {
    return(list())
  }
  if (!base::is.list(asset_urls)) {
    base::stop("`asset_urls` for dbAPIS must be a list when supplied.", call. = FALSE)
  }
  asset_urls
}

.dnmb_dbapis_copy_download_assets <- function(src_data_dir, layout) {
  required_files <- c("dbAPIS.hmm.zip", "seed_family_mapping.tsv", "seed_and_familyrep_all_infor.tsv", "anti_defense.pep", "readme.txt")
  if (!base::dir.exists(src_data_dir)) {
    base::stop("dbAPIS data_download directory not found: ", src_data_dir, call. = FALSE)
  }
  base::dir.create(layout$data_dir, recursive = TRUE, showWarnings = FALSE)
  for (file_name in required_files) {
    src <- base::file.path(src_data_dir, file_name)
    dest <- base::file.path(layout$data_dir, file_name)
    if (!base::file.exists(src)) {
      base::stop("dbAPIS asset is missing: ", src, call. = FALSE)
    }
    ok <- base::file.copy(src, dest, overwrite = TRUE)
    if (!base::isTRUE(ok)) {
      base::stop("Failed to copy dbAPIS asset: ", src, call. = FALSE)
    }
  }
  utils::unzip(layout$hmm_zip_path, exdir = layout$data_dir, overwrite = TRUE)
  macosx_dir <- base::file.path(layout$data_dir, "__MACOSX")
  if (base::dir.exists(macosx_dir)) {
    unlink(macosx_dir, recursive = TRUE, force = TRUE)
  }
  invisible(layout$data_dir)
}

.dnmb_dbapis_prepare_repo <- function(layout,
                                      repo_source = .dnmb_dbapis_default_repo_url(),
                                      force = FALSE) {
  repo_source <- base::as.character(repo_source)[1]
  if (base::dir.exists(layout$repo_dir) && !base::isTRUE(force)) {
    return(.dnmb_dbapis_status_row("dbapis_repo", "cached", layout$repo_dir))
  }

  if (base::dir.exists(layout$repo_dir) && base::isTRUE(force)) {
    unlink(layout$repo_dir, recursive = TRUE, force = TRUE)
  }

  if (base::dir.exists(repo_source)) {
    .dnmb_defensefinder_copy_dir_contents(repo_source, layout$repo_dir)
    return(.dnmb_dbapis_status_row("dbapis_repo", "ok", layout$repo_dir))
  }

  run <- dnmb_run_external(
    "git",
    args = c("clone", "--depth", "1", repo_source, layout$repo_dir),
    required = FALSE
  )
  .dnmb_dbapis_status_row(
    "dbapis_repo",
    if (base::isTRUE(run$ok)) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    if (base::isTRUE(run$ok)) layout$repo_dir else (run$error %||% repo_source)
  )
}

.dnmb_dbapis_install_ready <- function(layout, manifest = NULL) {
  !base::is.null(manifest) &&
    base::isTRUE(manifest$install_ok) &&
    base::file.exists(layout$hmm_path) &&
    base::file.exists(layout$mapping_path) &&
    base::file.exists(layout$info_path) &&
    base::file.exists(layout$fasta_path)
}

dnmb_dbapis_install_module <- function(version = .dnmb_dbapis_default_version(),
                                       cache_root = NULL,
                                       install = TRUE,
                                       repo_url = .dnmb_dbapis_default_repo_url(),
                                       asset_urls = NULL,
                                       force = FALSE) {
  module <- .dnmb_dbapis_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_dbapis_asset_layout(module_dir)
  asset_urls <- .dnmb_dbapis_normalize_asset_urls(asset_urls)
  status <- .dnmb_dbapis_empty_status()
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)

  if (.dnmb_dbapis_install_ready(layout, manifest = manifest) && !base::isTRUE(force)) {
    return(list(
      ok = TRUE,
      status = dplyr::bind_rows(status, .dnmb_dbapis_status_row("dbapis_install", "cached", module_dir)),
      files = list(hmm = layout$hmm_path, mapping = layout$mapping_path, info = layout$info_path, fasta = layout$fasta_path),
      manifest = manifest
    ))
  }

  if (!base::isTRUE(install)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_dbapis_status_row("dbapis_install", "missing", "dbAPIS is missing and module_install is FALSE.")),
      files = list(),
      manifest = manifest
    ))
  }

  repo_status <- .dnmb_dbapis_prepare_repo(
    layout,
    repo_source = asset_urls$repo_dir %||% asset_urls$repo_url %||% repo_url,
    force = force
  )
  status <- dplyr::bind_rows(status, repo_status)
  if (!repo_status$status %in% c("ok", "cached")) {
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }

  asset_status <- tryCatch({
    .dnmb_dbapis_copy_download_assets(base::file.path(layout$repo_dir, "data_download"), layout)
    .dnmb_dbapis_status_row("dbapis_assets", "ok", layout$data_dir)
  }, error = function(e) {
    .dnmb_dbapis_status_row("dbapis_assets", "failed", conditionMessage(e))
  })
  status <- dplyr::bind_rows(status, asset_status)
  if (asset_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }

  manifest <- list(
    install_ok = TRUE,
    repo_dir = layout$repo_dir,
    data_dir = layout$data_dir,
    hmm_path = layout$hmm_path,
    mapping_path = layout$mapping_path,
    info_path = layout$info_path,
    fasta_path = layout$fasta_path,
    readme_path = layout$readme_path
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_dbapis_default_version(),
    cache_root = cache_root
  )

  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_dbapis_status_row("dbapis_install", "ok", module_dir)),
    files = list(hmm = layout$hmm_path, mapping = layout$mapping_path, info = layout$info_path, fasta = layout$fasta_path),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_dbapis_get_module <- function(version = .dnmb_dbapis_default_version(),
                                   cache_root = NULL,
                                   required = TRUE) {
  manifest <- dnmb_db_read_manifest(.dnmb_dbapis_module_name(), version, cache_root = cache_root, required = FALSE)
  if (base::is.null(manifest)) {
    if (base::isTRUE(required)) {
      base::stop("dbAPIS module is not installed.", call. = FALSE)
    }
    return(list(ok = FALSE, manifest = NULL))
  }
  list(
    ok = TRUE,
    manifest = manifest,
    files = list(
      data_dir = manifest$data_dir,
      hmm = manifest$hmm_path,
      mapping = manifest$mapping_path,
      info = manifest$info_path,
      fasta = manifest$fasta_path
    )
  )
}

.dnmb_dbapis_prepare_hmm_db <- function(hmm_path, force = FALSE) {
  hmm_path <- base::path.expand(base::as.character(hmm_path)[1])
  required_files <- paste0(hmm_path, c(".h3m", ".h3i", ".h3f", ".h3p"))
  if (!base::file.exists(hmm_path)) {
    return(list(ok = FALSE, status = "missing", detail = hmm_path))
  }
  if (!base::isTRUE(force) && base::all(base::file.exists(required_files))) {
    return(list(ok = TRUE, status = "ok", detail = hmm_path, files = required_files))
  }
  run <- dnmb_run_external("hmmpress", args = c(hmm_path), required = FALSE)
  list(
    ok = base::isTRUE(run$ok),
    status = if (base::isTRUE(run$ok)) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    detail = if (base::isTRUE(run$ok)) hmm_path else (run$error %||% hmm_path),
    files = required_files
  )
}

.dnmb_dbapis_family_metadata <- function(info_path, mapping_path, include_acr = FALSE) {
  info_tbl <- utils::read.delim(info_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  mapping_tbl <- utils::read.delim(mapping_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  info_tbl$family_id <- base::as.character(info_tbl[["APIS families"]])
  info_tbl$hit_label <- base::as.character(info_tbl[["APIS genes"]])
  info_tbl$defense_type <- base::as.character(info_tbl[["Defense systems"]])
  info_tbl$representative_protein <- base::as.character(info_tbl[["Representative protein"]])
  info_tbl$description <- base::as.character(info_tbl[["Description"]])
  info_tbl$clan_id <- base::as.character(info_tbl[["clan ID"]])
  info_tbl$clan_defense_type <- base::as.character(info_tbl[["clan_system"]])
  info_tbl$profile_name <- base::sub("\\.hmm$", "", base::as.character(info_tbl[["Family HMM"]]))

  mapping_summary <- mapping_tbl |>
    dplyr::group_by(.data$family_ID) |>
    dplyr::summarise(
      mapping_defense_type = dplyr::first(.data[["inhibited_defense_system(inferred from seed proteins in the family)"]]),
      mapping_clan_id = dplyr::first(.data$clan_ID),
      mapping_clan_defense_type = dplyr::first(.data[["clan_inhibited_defense_system(inferred from families with seed proteins)"]]),
      .groups = "drop"
    )

  meta <- info_tbl |>
    dplyr::left_join(mapping_summary, by = c("family_id" = "family_ID")) |>
    dplyr::mutate(
      defense_type = dplyr::coalesce(.data$defense_type, .data$mapping_defense_type),
      clan_id = dplyr::coalesce(.data$clan_id, .data$mapping_clan_id),
      clan_defense_type = dplyr::coalesce(.data$clan_defense_type, .data$mapping_clan_defense_type)
    ) |>
    dplyr::select(.data$family_id, .data$profile_name, .data$hit_label, .data$defense_type, .data$representative_protein, .data$description, .data$clan_id, .data$clan_defense_type)

  meta$family_id <- ifelse(base::is.na(meta$family_id) | !base::nzchar(meta$family_id), meta$profile_name, meta$family_id)
  meta <- meta[!base::is.na(meta$family_id) & base::nzchar(meta$family_id), , drop = FALSE]
  if (!base::isTRUE(include_acr)) {
    meta <- meta[grepl("^APIS", meta$family_id), , drop = FALSE]
  }
  meta <- meta[!base::duplicated(meta$family_id), , drop = FALSE]
  base::rownames(meta) <- NULL
  meta
}

.dnmb_dbapis_empty_hits <- function() {
  tibble::tibble(
    query = character(),
    family_id = character(),
    hit_label = character(),
    defense_type = character(),
    clan_id = character(),
    clan_defense_type = character(),
    representative_protein = character(),
    description = character(),
    i_evalue = numeric(),
    score = numeric(),
    profile_length = integer(),
    gene_length = integer(),
    hmm_from = integer(),
    hmm_to = integer(),
    ali_from = integer(),
    ali_to = integer(),
    hmm_coverage = numeric(),
    query_coverage = numeric()
  )
}

dnmb_dbapis_parse_domtblout <- function(path,
                                        metadata,
                                        evalue_threshold = 1e-5,
                                        hmm_coverage_threshold = 0.35,
                                        include_acr = FALSE) {
  if (!base::file.exists(path) || !isTRUE(base::file.info(path)$size > 0)) {
    return(.dnmb_dbapis_empty_hits())
  }
  lines <- base::readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  if (!base::length(lines)) {
    return(.dnmb_dbapis_empty_hits())
  }

  meta_index <- metadata
  base::rownames(meta_index) <- meta_index$family_id

  parsed <- lapply(lines, function(line) {
    fields <- strsplit(base::trimws(line), "\\s+")[[1]]
    if (base::length(fields) < 22L) {
      return(NULL)
    }
    family_id <- base::as.character(fields[[1]])
    if (!base::isTRUE(include_acr) && grepl("^Acr", family_id)) {
      return(NULL)
    }
    meta_row <- meta_index[family_id, , drop = FALSE]
    if (!base::nrow(meta_row)) {
      meta_row <- data.frame(
        family_id = family_id,
        hit_label = family_id,
        defense_type = NA_character_,
        clan_id = NA_character_,
        clan_defense_type = NA_character_,
        representative_protein = NA_character_,
        description = NA_character_,
        stringsAsFactors = FALSE
      )
    }
    profile_length <- suppressWarnings(base::as.integer(fields[[3]]))
    gene_id <- base::as.character(fields[[4]])
    gene_length <- suppressWarnings(base::as.integer(fields[[6]]))
    i_evalue <- suppressWarnings(base::as.numeric(fields[[13]]))
    score <- suppressWarnings(base::as.numeric(fields[[14]]))
    hmm_from <- suppressWarnings(base::as.integer(fields[[16]]))
    hmm_to <- suppressWarnings(base::as.integer(fields[[17]]))
    ali_from <- suppressWarnings(base::as.integer(fields[[18]]))
    ali_to <- suppressWarnings(base::as.integer(fields[[19]]))
    hmm_coverage <- if (base::is.na(profile_length) || profile_length <= 0L || base::is.na(hmm_from) || base::is.na(hmm_to)) {
      NA_real_
    } else {
      (hmm_to - hmm_from + 1) / profile_length
    }
    query_coverage <- if (base::is.na(gene_length) || gene_length <= 0L || base::is.na(ali_from) || base::is.na(ali_to)) {
      NA_real_
    } else {
      (ali_to - ali_from + 1) / gene_length
    }

    tibble::tibble(
      query = .dnmb_module_clean_annotation_key(gene_id),
      family_id = family_id,
      hit_label = base::as.character(meta_row$hit_label[[1]] %||% family_id),
      defense_type = base::as.character(meta_row$defense_type[[1]] %||% NA_character_),
      clan_id = base::as.character(meta_row$clan_id[[1]] %||% NA_character_),
      clan_defense_type = base::as.character(meta_row$clan_defense_type[[1]] %||% NA_character_),
      representative_protein = base::as.character(meta_row$representative_protein[[1]] %||% NA_character_),
      description = base::as.character(meta_row$description[[1]] %||% NA_character_),
      i_evalue = i_evalue,
      score = score,
      profile_length = profile_length,
      gene_length = gene_length,
      hmm_from = hmm_from,
      hmm_to = hmm_to,
      ali_from = ali_from,
      ali_to = ali_to,
      hmm_coverage = hmm_coverage,
      query_coverage = query_coverage
    )
  })

  hits <- dplyr::bind_rows(parsed)
  if (!base::nrow(hits)) {
    return(.dnmb_dbapis_empty_hits())
  }
  hits <- hits[
    !base::is.na(hits$i_evalue) &
      hits$i_evalue <= evalue_threshold &
      !base::is.na(hits$hmm_coverage) &
      hits$hmm_coverage >= hmm_coverage_threshold,
    ,
    drop = FALSE
  ]
  if (!base::nrow(hits)) {
    return(.dnmb_dbapis_empty_hits())
  }
  hits <- hits[base::order(hits$query, hits$i_evalue, -hits$score, -hits$hmm_coverage), , drop = FALSE]
  hits <- hits[!base::duplicated(hits[, c("query", "family_id")]), , drop = FALSE]
  base::rownames(hits) <- NULL
  tibble::as_tibble(hits)
}

dnmb_dbapis_normalize_hits <- function(hits) {
  if (base::is.null(hits) || !base::nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  support <- base::paste0(
    "family=", hits$family_id,
    "; defense_type=", hits$defense_type,
    "; i_evalue=", signif(hits$i_evalue, 4),
    "; score=", base::round(hits$score, 2),
    "; hmm_cov=", sprintf("%.3f", hits$hmm_coverage),
    "; query_cov=", sprintf("%.3f", hits$query_coverage)
  )
  out <- base::data.frame(
    query = .dnmb_module_clean_annotation_key(hits$query),
    source = base::rep("dbapis", base::nrow(hits)),
    family_system = base::rep("dbAPIS", base::nrow(hits)),
    family_id = base::as.character(hits$family_id),
    hit_label = base::as.character(hits$hit_label),
    enzyme_role = base::as.character(hits$defense_type),
    evidence_mode = base::rep("direct", base::nrow(hits)),
    substrate_label = base::as.character(hits$clan_defense_type),
    support = support,
    defense_type = base::as.character(hits$defense_type),
    clan_id = base::as.character(hits$clan_id),
    clan_defense_type = base::as.character(hits$clan_defense_type),
    representative_protein = base::as.character(hits$representative_protein),
    description = base::as.character(hits$description),
    i_evalue = base::as.numeric(hits$i_evalue),
    hit_score = base::as.numeric(hits$score),
    hmm_coverage = base::as.numeric(hits$hmm_coverage),
    query_coverage = base::as.numeric(hits$query_coverage),
    profile_length = base::as.integer(hits$profile_length),
    gene_length = base::as.integer(hits$gene_length),
    hmm_from = base::as.integer(hits$hmm_from),
    hmm_to = base::as.integer(hits$hmm_to),
    ali_from = base::as.integer(hits$ali_from),
    ali_to = base::as.integer(hits$ali_to),
    typing_eligible = base::rep(TRUE, base::nrow(hits)),
    stringsAsFactors = FALSE
  )
  out <- out[base::order(out$query, out$i_evalue, -out$hmm_coverage), , drop = FALSE]
  base::rownames(out) <- NULL
  out[, c(.dnmb_module_optional_long_columns(), base::setdiff(base::names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
}

.dnmb_dbapis_output_table <- function(genes, hits) {
  out <- .dnmb_module_output_table(genes = genes, hits = hits)
  ordered <- c(
    base::intersect(dnmb_backbone_columns(), base::names(out)),
    base::intersect(
      c("family_id", "hit_label", "defense_type", "clan_id", "clan_defense_type", "representative_protein", "description", "i_evalue", "hit_score", "hmm_coverage", "query_coverage", "support"),
      base::names(out)
    ),
    base::setdiff(base::names(out), c(base::intersect(dnmb_backbone_columns(), base::names(out)), "family_id", "hit_label", "defense_type", "clan_id", "clan_defense_type", "representative_protein", "description", "i_evalue", "hit_score", "hmm_coverage", "query_coverage", "support"))
  )
  out[, ordered, drop = FALSE]
}

dnmb_run_dbapis_module <- function(genes,
                                   output_dir,
                                   version = .dnmb_dbapis_default_version(),
                                   cache_root = NULL,
                                   install = TRUE,
                                   repo_url = .dnmb_dbapis_default_repo_url(),
                                   asset_urls = NULL,
                                   cpu = 1L,
                                   evalue_threshold = 1e-5,
                                   hmm_coverage_threshold = 0.35,
                                   include_acr = FALSE) {
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_dbapis_empty_status()
  trace_log <- base::file.path(output_dir, "dbapis_module_trace.log")

  install_result <- dnmb_dbapis_install_module(
    version = version,
    cache_root = cache_root,
    install = install,
    repo_url = repo_url,
    asset_urls = asset_urls
  )
  status <- dplyr::bind_rows(status, install_result$status)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), raw_hits = .dnmb_dbapis_empty_hits(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  module <- dnmb_dbapis_get_module(version = version, cache_root = cache_root, required = TRUE)
  prepare_result <- .dnmb_dbapis_prepare_hmm_db(module$files$hmm, force = FALSE)
  status <- dplyr::bind_rows(status, .dnmb_dbapis_status_row("dbapis_prepare", prepare_result$status, prepare_result$detail))
  if (!base::isTRUE(prepare_result$ok)) {
    return(list(ok = FALSE, status = status, files = c(list(trace_log = trace_log), module$files), raw_hits = .dnmb_dbapis_empty_hits(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  query_info <- .dnmb_write_query_fasta(genes, base::file.path(output_dir, "dbapis_query_proteins.faa"))
  domtblout <- base::file.path(output_dir, "dbapis_hmmscan.domtblout")
  stdout_log <- base::file.path(output_dir, "dbapis_hmmscan.out")
  stderr_log <- base::file.path(output_dir, "dbapis_hmmscan.err")
  raw_hits_path <- base::file.path(output_dir, "dbapis_hits.tsv")

  cmd <- c(
    "--domtblout", domtblout,
    "--cpu", base::as.character(base::as.integer(cpu)[1]),
    "-o", stdout_log,
    module$files$hmm,
    query_info$path
  )
  base::cat(base::paste0("[", base::Sys.time(), "] ", .dnmb_format_command("hmmscan", cmd), "\n"), file = trace_log, append = TRUE)
  run <- dnmb_run_external("hmmscan", args = cmd, required = FALSE)
  base::writeLines(run$stderr, con = stderr_log)

  metadata <- .dnmb_dbapis_family_metadata(module$files$info, module$files$mapping, include_acr = include_acr)
  raw_hits <- dnmb_dbapis_parse_domtblout(
    domtblout,
    metadata = metadata,
    evalue_threshold = evalue_threshold,
    hmm_coverage_threshold = hmm_coverage_threshold,
    include_acr = include_acr
  )
  if (base::nrow(raw_hits)) {
    utils::write.table(raw_hits, file = raw_hits_path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  hits <- dnmb_dbapis_normalize_hits(raw_hits)
  output_table <- .dnmb_dbapis_output_table(genes = genes, hits = hits)

  status <- dplyr::bind_rows(
    status,
    .dnmb_dbapis_status_row(
      "dbapis_run",
      if (base::isTRUE(run$ok)) if (base::nrow(raw_hits)) "ok" else "empty" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
      if (base::isTRUE(run$ok)) domtblout else (run$error %||% domtblout)
    )
  )

  list(
    ok = base::isTRUE(run$ok),
    status = status,
    files = list(
      trace_log = trace_log,
      query_fasta = query_info$path,
      domtblout = domtblout,
      stdout_log = stdout_log,
      stderr_log = stderr_log,
      hits = raw_hits_path
    ),
    raw_hits = raw_hits,
    hits = hits,
    output_table = output_table
  )
}

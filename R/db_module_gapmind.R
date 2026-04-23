.dnmb_gapmind_module_name <- function() {
  "gapmind"
}

.dnmb_gapmind_default_version <- function() {
  "aa"
}

.dnmb_gapmind_supported_versions <- function() {
  c("aa", "carbon")
}

.dnmb_gapmind_validate_version <- function(version) {
  value <- trimws(as.character(version)[1])
  if (is.na(value) || !nzchar(value) || !value %in% .dnmb_gapmind_supported_versions()) {
    base::stop(
      "`version` for GapMind must be one of: ",
      paste(.dnmb_gapmind_supported_versions(), collapse = ", "),
      call. = FALSE
    )
  }
  value
}

.dnmb_gapmind_database_name <- function(version) {
  version <- .dnmb_gapmind_validate_version(version)
  if (identical(version, "aa")) {
    return("GapMindAA")
  }
  "GapMindCarbon"
}

.dnmb_gapmind_main_score_threshold <- function() {
  2
}

.dnmb_gapmind_main_require_best_path <- function() {
  TRUE
}

dnmb_gapmind_carbohydrate_pathways <- function() {
  c(
    "arabinose",
    "cellobiose",
    "deoxyribonate",
    "deoxyribose",
    "fructose",
    "fucose",
    "galactose",
    "galacturonate",
    "gluconate",
    "glucose",
    "glucose-6-P",
    "glucosamine",
    "glucuronate",
    "lactose",
    "maltose",
    "mannitol",
    "mannose",
    "myoinositol",
    "NAG",
    "rhamnose",
    "ribose",
    "sorbitol",
    "sucrose",
    "trehalose",
    "xylitol",
    "xylose"
  )
}

.dnmb_gapmind_default_repo_url <- function() {
  "https://github.com/morgannprice/PaperBLAST.git"
}

.dnmb_gapmind_default_resource_base_url <- function(version = .dnmb_gapmind_default_version()) {
  version <- .dnmb_gapmind_validate_version(version)
  base::paste0("https://papers.genomics.lbl.gov/tmp/path.", version)
}

.dnmb_gapmind_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_gapmind_empty_status <- function() {
  .dnmb_gapmind_status_row(character(), character(), character())
}

.dnmb_gapmind_trace <- function(path, text) {
  base::cat(base::paste0(text, "\n"), file = path, append = TRUE)
}

.dnmb_gapmind_asset_layout <- function(module_dir, version = .dnmb_gapmind_default_version()) {
  version <- .dnmb_gapmind_validate_version(version)
  path_dir <- base::file.path(module_dir, base::paste0("path.", version))
  list(
    module_dir = module_dir,
    version = version,
    repo_dir = base::file.path(module_dir, "PaperBLAST"),
    path_dir = path_dir,
    curated_faa = base::file.path(path_dir, "curated.faa"),
    curated_db = base::file.path(path_dir, "curated.db"),
    steps_db = base::file.path(path_dir, "steps.db"),
    curated_dmnd = base::file.path(path_dir, "curated.faa.dmnd"),
    install_trace_log = base::file.path(module_dir, "gapmind_install_trace.log")
  )
}

.dnmb_gapmind_normalize_asset_urls <- function(asset_urls = NULL) {
  if (base::is.null(asset_urls)) {
    return(list())
  }
  if (!base::is.list(asset_urls)) {
    base::stop("`asset_urls` for GapMind must be a list when supplied.", call. = FALSE)
  }
  asset_urls
}

.dnmb_gapmind_copy_dir_contents <- function(src_dir, dest_dir) {
  base::dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  entries <- base::list.files(src_dir, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  if (!base::length(entries)) {
    return(invisible(dest_dir))
  }
  ok <- base::file.copy(entries, dest_dir, recursive = TRUE, copy.mode = TRUE)
  if (!base::all(ok)) {
    base::stop("Failed to copy GapMind assets from ", src_dir, " to ", dest_dir, call. = FALSE)
  }
  invisible(dest_dir)
}

.dnmb_gapmind_perl_candidates <- function() {
  candidates <- c(
    "/opt/biotools/bin/perl",
    dnmb_detect_binary("perl", required = FALSE)$path,
    base::Sys.which("perl"),
    "/usr/bin/perl"
  )
  candidates <- unique(base::path.expand(base::trimws(base::as.character(candidates))))
  candidates[base::nzchar(candidates) & base::file.exists(candidates)]
}

.dnmb_gapmind_find_perl <- function(required_modules = c("DBI", "DBD::SQLite")) {
  required_modules <- base::as.character(required_modules)
  required_modules <- required_modules[!base::is.na(required_modules) & base::nzchar(required_modules)]
  perl_candidates <- .dnmb_gapmind_perl_candidates()
  if (!base::length(perl_candidates)) {
    return("")
  }

  args <- c(
    if (base::length(required_modules)) base::paste0("-M", required_modules) else character(),
    "-e",
    "print qq{ok\\n}"
  )
  for (perl_path in perl_candidates) {
    run <- dnmb_run_external(perl_path, args = args, required = FALSE)
    if (base::isTRUE(run$ok)) {
      return(perl_path)
    }
  }
  ""
}

.dnmb_gapmind_runtime_env <- function(perl_path = "", stage_root = "") {
  path_entries <- c(
    if (base::nzchar(stage_root)) base::file.path(stage_root, "bin") else character(),
    if (base::nzchar(perl_path)) base::dirname(perl_path) else character(),
    strsplit(base::Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  )
  path_entries <- unique(path_entries[base::nzchar(path_entries) & base::file.exists(path_entries)])
  c(PATH = base::paste(path_entries, collapse = .Platform$path.sep))
}

.dnmb_gapmind_prepare_repo <- function(layout,
                                       repo_source = .dnmb_gapmind_default_repo_url(),
                                       force = FALSE) {
  repo_source <- base::as.character(repo_source)[1]
  if (base::dir.exists(layout$repo_dir) && !base::isTRUE(force)) {
    return(.dnmb_gapmind_status_row("gapmind_repo", "cached", layout$repo_dir))
  }

  if (base::dir.exists(layout$repo_dir) && base::isTRUE(force)) {
    unlink(layout$repo_dir, recursive = TRUE, force = TRUE)
  }

  if (base::dir.exists(repo_source)) {
    .dnmb_gapmind_copy_dir_contents(repo_source, layout$repo_dir)
    return(.dnmb_gapmind_status_row("gapmind_repo", "ok", layout$repo_dir))
  }

  run <- dnmb_run_external(
    "git",
    args = c("clone", "--depth", "1", repo_source, layout$repo_dir),
    required = FALSE
  )
  .dnmb_gapmind_status_row(
    "gapmind_repo",
    if (base::isTRUE(run$ok)) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    if (base::isTRUE(run$ok)) layout$repo_dir else (run$error %||% repo_source)
  )
}

.dnmb_gapmind_required_paths <- function(layout) {
  c(layout$curated_faa, layout$curated_db, layout$steps_db, layout$curated_dmnd)
}

.dnmb_gapmind_bundled_table <- function(layout) {
  base::file.path(layout$repo_dir, "gaps", layout$version, base::paste0(layout$version, ".table"))
}

.dnmb_gapmind_bundled_requires <- function(layout) {
  base::file.path(layout$repo_dir, "gaps", layout$version, "requires.tsv")
}

.dnmb_gapmind_extract_hmms <- function(layout) {
  existing <- base::list.files(layout$path_dir, pattern = "[.]hmm$", full.names = TRUE)
  if (base::length(existing)) {
    return(.dnmb_gapmind_status_row("gapmind_hmms", "cached", base::paste0("count=", base::length(existing))))
  }

  perl_path <- .dnmb_gapmind_find_perl(required_modules = c("DBI", "DBD::SQLite"))
  if (!base::nzchar(perl_path)) {
    return(.dnmb_gapmind_status_row(
      "gapmind_hmms",
      "missing",
      "No perl executable with DBI and DBD::SQLite was found for GapMind."
    ))
  }
  script_path <- base::file.path(layout$repo_dir, "bin", "extractHmms.pl")
  run <- dnmb_run_external(
    perl_path,
    args = c(script_path, layout$steps_db, layout$path_dir),
    env = .dnmb_gapmind_runtime_env(perl_path = perl_path),
    required = FALSE
  )
  existing <- base::list.files(layout$path_dir, pattern = "[.]hmm$", full.names = TRUE)
  .dnmb_gapmind_status_row(
    "gapmind_hmms",
    if (base::isTRUE(run$ok) && base::length(existing)) "ok" else "failed",
    if (base::length(existing)) base::paste0("count=", base::length(existing)) else (run$error %||% layout$path_dir)
  )
}

.dnmb_gapmind_prepare_diamond_db <- function(layout) {
  if (base::file.exists(layout$curated_dmnd)) {
    return(.dnmb_gapmind_status_row("gapmind_diamond", "cached", layout$curated_dmnd))
  }

  run <- dnmb_run_external(
    "diamond",
    args = c("makedb", "--in", layout$curated_faa, "-d", sub("[.]dmnd$", "", layout$curated_dmnd)),
    required = FALSE
  )
  .dnmb_gapmind_status_row(
    "gapmind_diamond",
    if (base::isTRUE(run$ok) && base::file.exists(layout$curated_dmnd)) "ok" else "failed",
    if (base::file.exists(layout$curated_dmnd)) layout$curated_dmnd else (run$error %||% layout$curated_faa)
  )
}

dnmb_gapmind_install_module <- function(version = .dnmb_gapmind_default_version(),
                                        cache_root = NULL,
                                        install = TRUE,
                                        repo_url = .dnmb_gapmind_default_repo_url(),
                                        resource_base_url = .dnmb_gapmind_default_resource_base_url(version),
                                        asset_urls = NULL,
                                        force = FALSE) {
  version <- .dnmb_gapmind_validate_version(version)
  module <- .dnmb_gapmind_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_gapmind_asset_layout(module_dir, version = version)
  base::dir.create(layout$path_dir, recursive = TRUE, showWarnings = FALSE)
  asset_urls <- .dnmb_gapmind_normalize_asset_urls(asset_urls)
  status <- .dnmb_gapmind_empty_status()

  if (base::all(base::file.exists(.dnmb_gapmind_required_paths(layout))) &&
      base::file.exists(.dnmb_gapmind_bundled_table(layout)) &&
      base::dir.exists(layout$repo_dir) &&
      !base::isTRUE(force)) {
    manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
    if (!base::is.null(manifest) && base::isTRUE(manifest$install_ok)) {
      return(list(
        ok = TRUE,
        status = .dnmb_gapmind_status_row("gapmind_install", "cached", module_dir),
        files = list(
          repo_dir = layout$repo_dir,
          path_dir = layout$path_dir,
          curated_faa = layout$curated_faa,
          curated_db = layout$curated_db,
          steps_db = layout$steps_db,
          curated_dmnd = layout$curated_dmnd,
          pathway_table = .dnmb_gapmind_bundled_table(layout),
          aa_table = .dnmb_gapmind_bundled_table(layout),
          requires_tsv = .dnmb_gapmind_bundled_requires(layout)
        ),
        manifest = manifest
      ))
    }
  }

  if (!base::isTRUE(install)) {
    return(list(
      ok = FALSE,
      status = .dnmb_gapmind_status_row("gapmind_install", "missing", "GapMind resources are missing and module_install is FALSE."),
      files = list(),
      manifest = NULL
    ))
  }

  if (base::interactive()) {
    answer <- utils::askYesNo(base::paste0("GapMind reference data are missing. Download the official GapMind ", version, " resources?"))
    if (!base::isTRUE(answer)) {
      return(list(
        ok = FALSE,
        status = .dnmb_gapmind_status_row("gapmind_install", "declined", "User declined GapMind download."),
        files = list(),
        manifest = NULL
      ))
    }
  }

  repo_status <- .dnmb_gapmind_prepare_repo(
    layout = layout,
    repo_source = asset_urls$repo_dir %||% asset_urls$repo_url %||% repo_url,
    force = force
  )
  status <- dplyr::bind_rows(status, repo_status)
  if (!repo_status$status %in% c("ok", "cached")) {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  url_map <- list(
    curated_faa = asset_urls$curated_faa %||% base::paste0(resource_base_url, "/curated.faa"),
    curated_db = asset_urls$curated_db %||% base::paste0(resource_base_url, "/curated.db"),
    steps_db = asset_urls$steps_db %||% base::paste0(resource_base_url, "/steps.db")
  )
  dest_map <- list(
    curated_faa = layout$curated_faa,
    curated_db = layout$curated_db,
    steps_db = layout$steps_db
  )

  for (asset_name in base::names(url_map)) {
    dest <- dest_map[[asset_name]]
    source <- url_map[[asset_name]]
    ok <- FALSE
    detail <- dest
    if (base::file.exists(dest) && !base::isTRUE(force)) {
      ok <- TRUE
      status <- dplyr::bind_rows(status, .dnmb_gapmind_status_row(asset_name, "cached", dest))
      next
    }
    if (base::file.exists(source)) {
      ok <- base::file.copy(source, dest, overwrite = TRUE)
      detail <- if (ok) dest else source
    } else {
      download <- .dnmb_download_asset(source, dest, insecure = FALSE)
      ok <- base::isTRUE(download$ok) && base::file.exists(dest)
      detail <- if (ok) dest else (download$error %||% source)
    }
    status <- dplyr::bind_rows(status, .dnmb_gapmind_status_row(asset_name, if (ok) "ok" else "failed", detail))
    if (!ok) {
      return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
    }
  }

  hmm_status <- .dnmb_gapmind_extract_hmms(layout)
  status <- dplyr::bind_rows(status, hmm_status)
  if (!hmm_status$status %in% c("ok", "cached")) {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  diamond_status <- .dnmb_gapmind_prepare_diamond_db(layout)
  status <- dplyr::bind_rows(status, diamond_status)
  if (!diamond_status$status %in% c("ok", "cached")) {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  manifest <- list(
    install_ok = TRUE,
    module = module,
    version = version,
    module_dir = module_dir,
    repo_dir = layout$repo_dir,
    path_dir = layout$path_dir,
    curated_faa = layout$curated_faa,
    curated_db = layout$curated_db,
    steps_db = layout$steps_db,
    curated_dmnd = layout$curated_dmnd,
    pathway_table = .dnmb_gapmind_bundled_table(layout),
    requires_tsv = .dnmb_gapmind_bundled_requires(layout),
    resource_base_url = resource_base_url,
    repo_url = repo_url
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)

  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_gapmind_status_row("gapmind_install", "ok", module_dir)),
    files = list(
      repo_dir = layout$repo_dir,
      path_dir = layout$path_dir,
      curated_faa = layout$curated_faa,
      curated_db = layout$curated_db,
      steps_db = layout$steps_db,
      curated_dmnd = layout$curated_dmnd,
      pathway_table = .dnmb_gapmind_bundled_table(layout),
      requires_tsv = .dnmb_gapmind_bundled_requires(layout)
    ),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_gapmind_get_module <- function(version = .dnmb_gapmind_default_version(),
                                    cache_root = NULL,
                                    required = TRUE) {
  version <- .dnmb_gapmind_validate_version(version)
  manifest <- dnmb_db_read_manifest(.dnmb_gapmind_module_name(), version, cache_root = cache_root, required = FALSE)
  if (base::is.null(manifest)) {
    if (base::isTRUE(required)) {
      base::stop("GapMind module is not installed for version `", version, "`.", call. = FALSE)
    }
    return(list(ok = FALSE, manifest = NULL))
  }
  list(
    ok = TRUE,
    manifest = manifest,
    files = list(
      repo_dir = manifest$repo_dir,
      path_dir = manifest$path_dir,
      curated_faa = manifest$curated_faa,
      curated_db = manifest$curated_db,
      steps_db = manifest$steps_db,
      curated_dmnd = manifest$curated_dmnd,
      pathway_table = manifest$pathway_table %||% manifest$aa_table,
      requires_tsv = manifest$requires_tsv
    )
  )
}

dnmb_gapmind_parse_pathway_table <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame(pathway_id = character(), pathway_desc = character(), stringsAsFactors = FALSE))
  }
  tbl <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!base::nrow(tbl)) {
    return(data.frame(pathway_id = character(), pathway_desc = character(), stringsAsFactors = FALSE))
  }
  base::names(tbl)[base::match("pathwayId", base::names(tbl))] <- "pathway_id"
  base::names(tbl)[base::match("desc", base::names(tbl))] <- "pathway_desc"
  tbl[, c("pathway_id", "pathway_desc"), drop = FALSE]
}

.dnmb_gapmind_parse_org_locus <- function(x) {
  x <- base::as.character(x)
  x <- trimws(x)
  x <- sub("^[^:]+:", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

dnmb_gapmind_parse_steps <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame(pathway = character(), step = character(), locus_tag = character(), on_best_path = integer(), stringsAsFactors = FALSE))
  }
  tbl <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!base::nrow(tbl)) {
    return(data.frame(pathway = character(), step = character(), locus_tag = character(), on_best_path = integer(), stringsAsFactors = FALSE))
  }
  tbl$locus_tag <- .dnmb_gapmind_parse_org_locus(tbl$locusId)
  tbl$on_best_path <- suppressWarnings(base::as.integer(tbl$onBestPath))
  tbl[, c("pathway", "step", "locus_tag", "on_best_path"), drop = FALSE]
}

dnmb_gapmind_parse_candidates <- function(path,
                                          steps_path = NULL,
                                          pathway_table = NULL) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  tbl <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!base::nrow(tbl)) {
    return(data.frame())
  }
  tbl$locus_tag <- .dnmb_gapmind_parse_org_locus(tbl$locusId)
  if ("locusId2" %in% base::names(tbl)) {
    tbl$locus_tag2 <- .dnmb_gapmind_parse_org_locus(tbl$locusId2)
  } else {
    tbl$locus_tag2 <- NA_character_
  }
  tbl$score <- suppressWarnings(base::as.numeric(tbl$score))
  tbl$blastBits <- suppressWarnings(base::as.numeric(tbl$blastBits))
  tbl$identity <- suppressWarnings(base::as.numeric(tbl$identity))
  tbl$blastCoverage <- suppressWarnings(base::as.numeric(tbl$blastCoverage))
  tbl$blastScore <- suppressWarnings(base::as.numeric(tbl$blastScore))
  tbl$hmmBits <- suppressWarnings(base::as.numeric(tbl$hmmBits))
  tbl$hmmCoverage <- suppressWarnings(base::as.numeric(tbl$hmmCoverage))
  tbl$hmmScore <- suppressWarnings(base::as.numeric(tbl$hmmScore))
  tbl$otherBits <- suppressWarnings(base::as.numeric(tbl$otherBits))
  tbl$otherIdentity <- suppressWarnings(base::as.numeric(tbl$otherIdentity))
  tbl$otherCoverage <- suppressWarnings(base::as.numeric(tbl$otherCoverage))

  if (!base::is.null(steps_path) && base::file.exists(steps_path)) {
    step_tbl <- dnmb_gapmind_parse_steps(steps_path)
    if (base::nrow(step_tbl)) {
      tbl <- dplyr::left_join(tbl, step_tbl, by = c("pathway" = "pathway", "step" = "step", "locus_tag" = "locus_tag"))
    }
  }
  if (base::is.null(tbl$on_best_path)) {
    tbl$on_best_path <- NA_integer_
  }
  if (!base::is.null(pathway_table) && base::is.data.frame(pathway_table) && base::nrow(pathway_table)) {
    tbl <- dplyr::left_join(tbl, pathway_table, by = c("pathway" = "pathway_id"))
  } else {
    tbl$pathway_desc <- NA_character_
  }
  tbl
}

.dnmb_gapmind_score_label <- function(score) {
  score <- suppressWarnings(base::as.numeric(score))
  out <- ifelse(
    base::is.na(score),
    NA_character_,
    ifelse(score >= 2, "high", ifelse(score >= 1, "medium", "low"))
  )
  base::as.character(out)
}

.dnmb_gapmind_support_string <- function(tbl) {
  n <- base::nrow(tbl)
  if (!n) {
    return(character())
  }
  parts <- base::vector("list", n)
  add_part <- function(label, values, digits = 4) {
    if (base::is.null(values)) {
      return()
    }
    vals <- values
    if (base::is.numeric(vals)) {
      vals <- ifelse(base::is.na(vals), NA_character_, base::format(signif(vals, digits), trim = TRUE))
    } else {
      vals <- base::as.character(vals)
    }
    for (i in seq_len(n)) {
      if (!base::is.na(vals[[i]]) && base::nzchar(trimws(vals[[i]]))) {
        parts[[i]] <<- c(parts[[i]], base::paste0(label, "=", vals[[i]]))
      }
    }
  }
  add_part("score", tbl$score, digits = 2)
  add_part("blastBits", tbl$blastBits)
  add_part("identity", tbl$identity)
  add_part("blastCoverage", tbl$blastCoverage)
  add_part("hmmBits", tbl$hmmBits)
  add_part("hmmCoverage", tbl$hmmCoverage)
  add_part("otherBits", tbl$otherBits)
  add_part("otherIdentity", tbl$otherIdentity)
  add_part("otherCoverage", tbl$otherCoverage)
  add_part("curated", tbl$curatedIds)
  add_part("hmm", tbl$hmmId)
  add_part("other", tbl$otherIds)
  add_part("split_locus", tbl$locus_tag2)
  add_part("on_best_path", tbl$on_best_path)
  vapply(parts, function(x) if (base::length(x)) base::paste(x, collapse = "; ") else NA_character_, character(1))
}

dnmb_gapmind_normalize_hits <- function(candidates, version = .dnmb_gapmind_default_version()) {
  if (base::is.null(candidates) || !base::is.data.frame(candidates) || !base::nrow(candidates)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  version <- .dnmb_gapmind_validate_version(version)

  tbl <- base::as.data.frame(candidates, stringsAsFactors = FALSE)
  for (column_name in c("pathway_desc", "curatedIds", "curatedDesc", "hmmId", "hmmDesc", "otherIds", "locus_tag2")) {
    if (base::is.null(tbl[[column_name]])) {
      tbl[[column_name]] <- base::rep(NA_character_, base::nrow(tbl))
    }
  }
  tbl$query <- .dnmb_module_clean_annotation_key(tbl$locus_tag)
  tbl$confidence <- .dnmb_gapmind_score_label(tbl$score)
  tbl$reference_id <- ifelse(
    !base::is.na(tbl$curatedIds) & base::nzchar(tbl$curatedIds),
    tbl$curatedIds,
    tbl$hmmId
  )
  tbl$support <- .dnmb_gapmind_support_string(tbl)
  tbl$candidate_type <- ifelse(
    !base::is.na(tbl$blastBits) & !base::is.na(tbl$hmmBits),
    "blast+hmm",
    ifelse(!base::is.na(tbl$blastBits), "blast", ifelse(!base::is.na(tbl$hmmBits), "hmm", NA_character_))
  )

  out <- data.frame(
    query = tbl$query,
    source = "gapmind",
    family_system = if (identical(version, "carbon")) "GapMind-Carbon" else "GapMind-AA",
    family_id = base::as.character(tbl$pathway),
    hit_label = base::as.character(tbl$step),
    enzyme_role = base::as.character(tbl$pathway_desc),
    evidence_mode = base::as.character(tbl$confidence),
    substrate_label = NA_character_,
    support = base::as.character(tbl$support),
    typing_eligible = tbl$on_best_path %in% 1L,
    pathway_id = base::as.character(tbl$pathway),
    pathway_desc = base::as.character(tbl$pathway_desc),
    step_id = base::as.character(tbl$step),
    step_score = suppressWarnings(base::as.numeric(tbl$score)),
    confidence = base::as.character(tbl$confidence),
    on_best_path = tbl$on_best_path %in% 1L,
    reference_id = base::as.character(tbl$reference_id),
    candidate_type = base::as.character(tbl$candidate_type),
    locus_tag2 = base::as.character(tbl$locus_tag2),
    curated_ids = base::as.character(tbl$curatedIds),
    curated_desc = base::as.character(tbl$curatedDesc),
    blast_bits = suppressWarnings(base::as.numeric(tbl$blastBits)),
    identity = suppressWarnings(base::as.numeric(tbl$identity)),
    blast_coverage = suppressWarnings(base::as.numeric(tbl$blastCoverage)),
    blast_score = suppressWarnings(base::as.numeric(tbl$blastScore)),
    hmm_bits = suppressWarnings(base::as.numeric(tbl$hmmBits)),
    hmm_id = base::as.character(tbl$hmmId),
    hmm_coverage = suppressWarnings(base::as.numeric(tbl$hmmCoverage)),
    hmm_score = suppressWarnings(base::as.numeric(tbl$hmmScore)),
    hmm_desc = base::as.character(tbl$hmmDesc),
    other_ids = base::as.character(tbl$otherIds),
    other_bits = suppressWarnings(base::as.numeric(tbl$otherBits)),
    other_identity = suppressWarnings(base::as.numeric(tbl$otherIdentity)),
    other_coverage = suppressWarnings(base::as.numeric(tbl$otherCoverage)),
    stringsAsFactors = FALSE
  )
  out <- out[!base::is.na(out$query) & base::nzchar(out$query), , drop = FALSE]
  base::rownames(out) <- NULL
  out
}

.dnmb_gapmind_hits_for_output <- function(hits,
                                          version = .dnmb_gapmind_default_version(),
                                          score_threshold = .dnmb_gapmind_main_score_threshold(),
                                          require_best_path = .dnmb_gapmind_main_require_best_path()) {
  if (base::is.null(hits) || !base::is.data.frame(hits) || !base::nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  version <- .dnmb_gapmind_validate_version(version)
  out <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  out$step_score <- suppressWarnings(base::as.numeric(out$step_score))
  out$priority_bits <- pmax(
    suppressWarnings(base::as.numeric(out$blast_bits)),
    suppressWarnings(base::as.numeric(out$hmm_bits)),
    na.rm = TRUE
  )
  out$priority_bits[!base::is.finite(out$priority_bits)] <- NA_real_
  out <- out[
    !base::is.na(out$step_score) & out$step_score >= score_threshold,
    ,
    drop = FALSE
  ]
  if (base::isTRUE(require_best_path) && "on_best_path" %in% base::names(out)) {
    keep_best <- out$on_best_path %in% TRUE
    out <- out[keep_best, , drop = FALSE]
  }
  if (identical(version, "carbon") && "pathway_id" %in% base::names(out)) {
    out <- out[out$pathway_id %in% dnmb_gapmind_carbohydrate_pathways(), , drop = FALSE]
  }
  if (!base::nrow(out)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out <- out[order(
    out$query,
    -(out$on_best_path %in% TRUE),
    -out$step_score,
    -ifelse(base::is.na(out$priority_bits), -Inf, out$priority_bits),
    out$pathway_id,
    out$step_id
  ), , drop = FALSE]
  out <- out[!duplicated(out$query), , drop = FALSE]
  if ("priority_bits" %in% base::names(out)) {
    out$priority_bits <- NULL
  }
  base::rownames(out) <- NULL
  out
}

.dnmb_gapmind_output_table <- function(genes,
                                       hits,
                                       version = .dnmb_gapmind_default_version(),
                                       score_threshold = .dnmb_gapmind_main_score_threshold(),
                                       require_best_path = .dnmb_gapmind_main_require_best_path()) {
  selected_hits <- .dnmb_gapmind_hits_for_output(
    hits,
    version = version,
    score_threshold = score_threshold,
    require_best_path = require_best_path
  )
  out <- .dnmb_module_output_table(genes = genes, hits = selected_hits)
  drop_cols <- base::intersect(
    c(
      "enzyme_role",
      "evidence_mode",
      "substrate_label",
      "typing_eligible",
      "blast_bits",
      "identity",
      "blast_coverage",
      "blast_score",
      "hmm_bits",
      "hmm_id",
      "hmm_coverage",
      "hmm_score",
      "hmm_desc",
      "other_ids",
      "other_bits",
      "other_identity",
      "other_coverage"
    ),
    base::names(out)
  )
  if (base::length(drop_cols)) {
    out[drop_cols] <- NULL
  }

  if (!"pathway_id" %in% base::names(out)) {
    out$pathway_id <- if ("family_id" %in% base::names(out)) base::as.character(out$family_id) else NA_character_
  } else if ("family_id" %in% base::names(out)) {
    missing_path <- base::is.na(out$pathway_id) | !base::nzchar(base::as.character(out$pathway_id))
    out$pathway_id[missing_path] <- base::as.character(out$family_id[missing_path])
  }

  if (!"step_id" %in% base::names(out)) {
    out$step_id <- if ("hit_label" %in% base::names(out)) base::as.character(out$hit_label) else NA_character_
  } else if ("hit_label" %in% base::names(out)) {
    missing_step <- base::is.na(out$step_id) | !base::nzchar(base::as.character(out$step_id))
    out$step_id[missing_step] <- base::as.character(out$hit_label[missing_step])
  }

  generic_drop <- base::intersect(c("family_id", "hit_label"), base::names(out))
  if (base::length(generic_drop)) {
    out[generic_drop] <- NULL
  }

  if ("pathway_desc" %in% base::names(out)) {
    out$pathway_desc <- base::as.character(out$pathway_desc)
  } else {
    out$pathway_desc <- NA_character_
  }
  ordered <- c(
    base::intersect(dnmb_backbone_columns(), base::names(out)),
    base::intersect(
      c("pathway_id", "pathway_desc", "step_id", "step_score", "confidence", "on_best_path", "reference_id", "candidate_type", "support", "locus_tag2"),
      base::names(out)
    ),
    base::setdiff(base::names(out), c(base::intersect(dnmb_backbone_columns(), base::names(out)), "pathway_id", "pathway_desc", "step_id", "step_score", "confidence", "on_best_path", "reference_id", "candidate_type", "support", "locus_tag2"))
  )
  out[, ordered, drop = FALSE]
}

.dnmb_gapmind_sanitize_id <- function(x, default = "dnmb") {
  value <- base::as.character(x)[1]
  value <- gsub("[^A-Za-z0-9_.-]", "_", value)
  value <- gsub("_+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  if (base::is.na(value) || !base::nzchar(value)) {
    value <- default
  }
  value
}

.dnmb_gapmind_org_prefix <- function(output_dir, genes, genbank = NULL) {
  prefix <- if (!base::is.null(genbank) && base::nzchar(base::as.character(genbank)[1])) {
    tools::file_path_sans_ext(base::basename(base::as.character(genbank)[1]))
  } else {
    "dnmb_gapmind"
  }
  prefix <- .dnmb_gapmind_sanitize_id(prefix, default = "dnmb_gapmind")
  base::file.path(output_dir, "gapmind_orgs", prefix)
}

.dnmb_gapmind_write_org_inputs <- function(genes, prefix, genome_name = "DNMB query genome") {
  proteins <- .dnmb_prepare_query_proteins(genes)
  org_id <- "local__DNMB"
  base::dir.create(base::dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  faa_path <- base::paste0(prefix, ".faa")
  org_path <- base::paste0(prefix, ".org")

  con <- base::file(faa_path, open = "w")
  on.exit(base::close(con), add = TRUE)
  for (i in seq_len(base::nrow(proteins))) {
    locus_id <- proteins$locus_tag[[i]]
    desc <- ""
    if ("product" %in% base::names(proteins)) {
      desc <- base::as.character(proteins$product[[i]])
    }
    desc <- gsub("[\r\n\t]+", " ", desc, perl = TRUE)
    desc <- trimws(gsub("[[:space:]]{2,}", " ", desc, perl = TRUE))
    base::writeLines(base::paste0(">", org_id, ":", locus_id, " ", locus_id, " ", desc), con = con)
    base::writeLines(proteins$translation[[i]], con = con)
  }

  org_tbl <- data.frame(
    orgId = org_id,
    gdb = "local",
    gid = "DNMB",
    genomeName = base::as.character(genome_name)[1],
    nProteins = base::nrow(proteins),
    stringsAsFactors = FALSE
  )
  utils::write.table(org_tbl, file = org_path, sep = "\t", row.names = FALSE, quote = FALSE)
  list(prefix = prefix, faa = faa_path, org = org_path, org_id = org_id, proteins = proteins)
}

.dnmb_gapmind_link_or_copy <- function(source, dest) {
  if (base::file.exists(dest) || base::dir.exists(dest)) {
    unlink(dest, recursive = TRUE, force = TRUE)
  }
  ok <- base::file.symlink(source, dest)
  if (!base::isTRUE(ok)) {
    ok <- base::file.copy(source, dest, recursive = TRUE, copy.mode = TRUE)
  }
  if (!base::isTRUE(ok)) {
    base::stop("Failed to stage GapMind asset from ", source, " to ", dest, call. = FALSE)
  }
  invisible(dest)
}

.dnmb_gapmind_stage_runtime <- function(module, output_dir, perl_path = "") {
  stage_root <- base::tempfile("dnmb-gapmind-run-")
  base::dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)
  base::dir.create(base::file.path(stage_root, "bin"), recursive = TRUE, showWarnings = FALSE)
  base::dir.create(base::file.path(stage_root, "lib"), recursive = TRUE, showWarnings = FALSE)
  base::dir.create(base::file.path(stage_root, "tmp", "path.aa"), recursive = TRUE, showWarnings = FALSE)

  .dnmb_gapmind_copy_dir_contents(base::file.path(module$files$repo_dir, "lib"), base::file.path(stage_root, "lib"))
  for (script_name in c("gapsearch.pl", "gaprevsearch.pl", "gapsummary.pl", "checkGapRequirements.pl")) {
    script_src <- base::file.path(module$files$repo_dir, "bin", script_name)
    script_dest <- base::file.path(stage_root, "bin", script_name)
    ok <- base::file.copy(script_src, script_dest, overwrite = TRUE, copy.mode = TRUE)
    if (!base::isTRUE(ok)) {
      base::stop("Failed to copy GapMind script from ", script_src, " to ", script_dest, call. = FALSE)
    }
  }

  for (binary_name in c("diamond", "hmmsearch", "hmmfetch")) {
    binary_path <- dnmb_detect_binary(binary_name, required = TRUE)$path
    .dnmb_gapmind_link_or_copy(binary_path, base::file.path(stage_root, "bin", binary_name))
  }
  if (base::nzchar(perl_path) && base::file.exists(perl_path)) {
    .dnmb_gapmind_link_or_copy(perl_path, base::file.path(stage_root, "bin", "perl"))
  }

  for (file_name in base::list.files(module$files$path_dir, full.names = FALSE)) {
    .dnmb_gapmind_link_or_copy(
      base::file.path(module$files$path_dir, file_name),
      base::file.path(stage_root, "tmp", "path.aa", file_name)
    )
  }

  list(
    stage_root = stage_root,
    path_dir = base::file.path(stage_root, "tmp", "path.aa"),
    cleanup = function() unlink(stage_root, recursive = TRUE, force = TRUE)
  )
}

dnmb_run_gapmind_module <- function(genes,
                                    output_dir,
                                    version = .dnmb_gapmind_default_version(),
                                    cache_root = NULL,
                                    install = TRUE,
                                    repo_url = .dnmb_gapmind_default_repo_url(),
                                    resource_base_url = .dnmb_gapmind_default_resource_base_url(version),
                                    asset_urls = NULL,
                                    cpu = 1L,
                                    genbank = NULL) {
  version <- .dnmb_gapmind_validate_version(version)
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- base::file.path(output_dir, "gapmind_module_trace.log")
  status <- .dnmb_gapmind_empty_status()

  install_result <- dnmb_gapmind_install_module(
    version = version,
    cache_root = cache_root,
    install = install,
    repo_url = repo_url,
    resource_base_url = resource_base_url,
    asset_urls = asset_urls
  )
  status <- dplyr::bind_rows(status, install_result$status)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), hits = .dnmb_module_empty_optional_long_table(), candidates = data.frame(), steps = data.frame(), rules = data.frame(), output_table = data.frame()))
  }

  module <- dnmb_gapmind_get_module(version = version, cache_root = cache_root, required = TRUE)
  perl_path <- .dnmb_gapmind_find_perl(required_modules = c("DBI", "DBD::SQLite"))
  status <- dplyr::bind_rows(
    status,
    .dnmb_gapmind_status_row(
      "gapmind_perl",
      if (base::nzchar(perl_path)) "ok" else "missing",
      if (base::nzchar(perl_path)) perl_path else "No perl executable with DBI and DBD::SQLite was found for GapMind."
    )
  )
  if (!base::nzchar(perl_path)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), hits = .dnmb_module_empty_optional_long_table(), candidates = data.frame(), steps = data.frame(), rules = data.frame(), output_table = data.frame()))
  }

  stage <- .dnmb_gapmind_stage_runtime(module, output_dir = output_dir, perl_path = perl_path)
  on.exit(stage$cleanup(), add = TRUE)
  gapmind_env <- .dnmb_gapmind_runtime_env(perl_path = perl_path, stage_root = stage$stage_root)

  org_prefix <- .dnmb_gapmind_org_prefix(stage$stage_root, genes = genes, genbank = genbank)
  genome_name <- if (!base::is.null(genbank) && base::nzchar(base::as.character(genbank)[1])) base::basename(base::as.character(genbank)[1]) else "DNMB query genome"
  org_inputs <- .dnmb_gapmind_write_org_inputs(genes = genes, prefix = org_prefix, genome_name = genome_name)

  org_dmnd <- base::paste0(org_prefix, ".faa.dmnd")
  diamond_make <- dnmb_run_external(
    "diamond",
    args = c("makedb", "--in", org_inputs$faa, "-d", sub("[.]dmnd$", "", org_dmnd)),
    required = FALSE
  )
  status <- dplyr::bind_rows(status, .dnmb_gapmind_status_row("gapmind_org_diamond", if (base::isTRUE(diamond_make$ok) && base::file.exists(org_dmnd)) "ok" else "failed", if (base::file.exists(org_dmnd)) org_dmnd else (diamond_make$error %||% org_inputs$faa)))
  if (!base::file.exists(org_dmnd)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), hits = .dnmb_module_empty_optional_long_table(), candidates = data.frame(), steps = data.frame(), rules = data.frame(), output_table = data.frame()))
  }

  hits_prefix <- base::file.path(stage$stage_root, version)
  gapsearch_args <- c(
    base::file.path(stage$stage_root, "bin", "gapsearch.pl"),
    "-diamond",
    "-orgs", org_prefix,
    "-set", version,
    "-dir", stage$path_dir,
    "-out", base::paste0(hits_prefix, ".hits"),
    "-nCPU", base::as.character(base::as.integer(cpu)[1])
  )
  .dnmb_gapmind_trace(trace_log, base::sprintf("[%s] %s", base::Sys.time(), .dnmb_format_command(perl_path, gapsearch_args)))
  gapsearch_run <- dnmb_run_external(perl_path, args = gapsearch_args, wd = stage$stage_root, env = gapmind_env, required = FALSE)
  status <- dplyr::bind_rows(status, .dnmb_gapmind_status_row("gapmind_search", if (base::isTRUE(gapsearch_run$ok)) "ok" else "failed", gapsearch_run$error %||% base::paste0(hits_prefix, ".hits")))
  if (!base::isTRUE(gapsearch_run$ok)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), hits = .dnmb_module_empty_optional_long_table(), candidates = data.frame(), steps = data.frame(), rules = data.frame(), output_table = data.frame()))
  }

  revhits_path <- base::paste0(hits_prefix, ".revhits")
  gaprev_args <- c(
    base::file.path(stage$stage_root, "bin", "gaprevsearch.pl"),
    "-diamond",
    "-hitsFile", base::paste0(hits_prefix, ".hits"),
    "-orgs", org_prefix,
    "-curated", base::file.path(stage$path_dir, "curated.faa.dmnd"),
    "-out", revhits_path,
    "-nCPU", base::as.character(base::as.integer(cpu)[1])
  )
  .dnmb_gapmind_trace(trace_log, base::sprintf("[%s] %s", base::Sys.time(), .dnmb_format_command(perl_path, gaprev_args)))
  gaprev_run <- dnmb_run_external(perl_path, args = gaprev_args, wd = stage$stage_root, env = gapmind_env, required = FALSE)
  status <- dplyr::bind_rows(status, .dnmb_gapmind_status_row("gapmind_revsearch", if (base::isTRUE(gaprev_run$ok)) "ok" else "failed", gaprev_run$error %||% revhits_path))
  if (!base::isTRUE(gaprev_run$ok)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), hits = .dnmb_module_empty_optional_long_table(), candidates = data.frame(), steps = data.frame(), rules = data.frame(), output_table = data.frame()))
  }

  summary_prefix <- base::file.path(stage$stage_root, version, "aa.sum")
  base::dir.create(base::dirname(summary_prefix), recursive = TRUE, showWarnings = FALSE)
  gapsummary_args <- c(
    base::file.path(stage$stage_root, "bin", "gapsummary.pl"),
    "-set", version,
    "-orgs", org_prefix,
    "-dbDir", stage$path_dir,
    "-hits", base::paste0(hits_prefix, ".hits"),
    "-revhits", revhits_path,
    "-out", summary_prefix
  )
  .dnmb_gapmind_trace(trace_log, base::sprintf("[%s] %s", base::Sys.time(), .dnmb_format_command(perl_path, gapsummary_args)))
  gapsummary_run <- dnmb_run_external(perl_path, args = gapsummary_args, wd = stage$stage_root, env = gapmind_env, required = FALSE)
  status <- dplyr::bind_rows(status, .dnmb_gapmind_status_row("gapmind_summary", if (base::isTRUE(gapsummary_run$ok)) "ok" else "failed", gapsummary_run$error %||% summary_prefix))

  cand_path <- base::paste0(summary_prefix, ".cand")
  steps_path <- base::paste0(summary_prefix, ".steps")
  rules_path <- base::paste0(summary_prefix, ".rules")
  out_files <- c(hits = base::paste0(hits_prefix, ".hits"), revhits = revhits_path, cand = cand_path, steps = steps_path, rules = rules_path)
  for (name in base::names(out_files)) {
    if (base::file.exists(out_files[[name]])) {
      base::file.copy(out_files[[name]], base::file.path(output_dir, base::basename(out_files[[name]])), overwrite = TRUE)
    }
  }

  pathway_tbl <- dnmb_gapmind_parse_pathway_table(module$files$pathway_table)
  candidates <- dnmb_gapmind_parse_candidates(cand_path, steps_path = steps_path, pathway_table = pathway_tbl)
  steps_tbl <- if (base::file.exists(steps_path)) utils::read.delim(steps_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
  rules_tbl <- if (base::file.exists(rules_path)) utils::read.delim(rules_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
  hits <- dnmb_gapmind_normalize_hits(candidates, version = version)
  output_table <- .dnmb_gapmind_output_table(genes = genes, hits = hits, version = version)

  list(
    ok = base::isTRUE(gapsummary_run$ok) && base::file.exists(cand_path) && base::file.exists(steps_path),
    status = status,
    files = list(
      trace_log = trace_log,
      hits = base::file.path(output_dir, base::basename(base::paste0(hits_prefix, ".hits"))),
      revhits = base::file.path(output_dir, base::basename(revhits_path)),
      cand = base::file.path(output_dir, base::basename(cand_path)),
      steps = base::file.path(output_dir, base::basename(steps_path)),
      rules = base::file.path(output_dir, base::basename(rules_path))
    ),
    candidates = candidates,
    steps = steps_tbl,
    rules = rules_tbl,
    hits = hits,
    output_table = output_table
  )
}

# ---------------------------------------------------------------------------
# Annotation-based supplementation for GapMind AA pathway steps
# ---------------------------------------------------------------------------

#' Mapping of amino acid biosynthesis pathway steps to EC numbers, gene names,
#' and product keywords. Used by annotation supplementation to rescue steps
#' that GapMind missed or scored low.
#' @return A data frame with columns: pathway, step, ec, gene_pattern, product_pattern.
#' @noRd
.dnmb_gapmind_aa_step_annotation_map <- function() {
 rows <- list(
  # --- Histidine ---
  list("his","hisA","5.3.1.16","^hisA$","phosphoribosylformimino.*AICAR.*isomerase|HisA"),
  list("his","hisB","3.6.1.31|3.1.3.15","^hisB$","imidazole glycerol.phosphate dehydratase|histidinol.phosphatase"),
  list("his","hisC","2.6.1.9","^hisC$","histidinol.phosphate aminotransferase"),
  list("his","hisD","1.1.1.23","^hisD$","histidinol dehydrogenase"),
  list("his","hisE","2.4.2.17","^hisE$","phosphoribosyl.ATP pyrophosphohydrolase"),
  list("his","hisF","4.1.3.-","^hisF$","imidazole glycerol.phosphate synthase|cyclase.*HisF|HisF"),
  list("his","hisG","2.4.2.17","^hisG$","ATP phosphoribosyltransferase|hisG"),
  list("his","hisH","2.4.2.17","^hisH$","glutamine amidotransferase.*HisH|amidotransferase.*imidazole"),
  list("his","hisI","3.5.4.19|3.6.1.31","^hisI$","phosphoribosyl.AMP cyclohydrolase|bifunctional.*hisIE"),
  list("his","hisN","3.1.3.15","^hisN$","histidinol.phosphatase|PHP domain"),
  list("his","prs","2.7.6.1","^prs$|^prsA$","ribose.phosphate pyrophosphokinase|PRPP synthetase"),

  # --- Serine / Glycine / Cysteine ---
  list("ser","serA","1.1.1.95","^serA$","phosphoglycerate dehydrogenase"),
  list("ser","serB","3.1.3.3","^serB$","phosphoserine phosphatase"),
  list("ser","serC","2.6.1.52","^serC$","phosphoserine aminotransferase"),
  list("gly","glyA","2.1.2.1","^glyA$","serine hydroxymethyltransferase"),
  list("cys","cysE","2.3.1.30","^cysE$","serine O-acetyltransferase"),
  list("cys","cysK","2.5.1.47","^cysK$","cysteine synthase|O-acetylserine.*thiol lyase"),
  list("cys","cysM","2.5.1.47","^cysM$","cysteine synthase B|O-acetylserine.*lyase"),

  # --- Chorismate / Aromatic ---
  list("chorismate","aroA","2.5.1.19","^aroA$","EPSPS|5-enolpyruvylshikimate-3-phosphate synthase"),
  list("chorismate","aroB","4.2.3.4","^aroB$","3-dehydroquinate synthase"),
  list("chorismate","aroC","4.2.1.10","^aroC$","chorismate synthase"),
  list("chorismate","aroD","4.2.1.10","^aroD$","3-dehydroquinate dehydratase"),
  list("chorismate","aroE","1.1.1.25","^aroE$","shikimate dehydrogenase"),
  list("chorismate","aroK","2.7.1.71","^aroK$|^aroL$","shikimate kinase"),
  list("phe","pheA","4.2.1.51|5.4.99.5","^pheA$","prephenate dehydratase|chorismate mutase"),
  list("tyr","tyrA","1.3.1.12|5.4.99.5","^tyrA$","prephenate dehydrogenase|arogenate dehydrogenase"),
  list("trp","trpA","4.2.1.20","^trpA$","tryptophan synthase.*alpha"),
  list("trp","trpB","4.2.1.20","^trpB$","tryptophan synthase.*beta"),
  list("trp","trpC","5.3.1.24","^trpC$","indole-3-glycerol.phosphate synthase"),
  list("trp","trpD","2.4.2.18","^trpD$","anthranilate phosphoribosyltransferase"),
  list("trp","trpE","4.1.3.27","^trpE$","anthranilate synthase"),

  # --- Branched-chain (shared: ilvB/ilvC/ilvD/ilvE across val/ile/leu) ---
  list("val","ilvB","2.2.1.6","^ilvB$|^ilvI$|^ilvG$|^ilvH$","acetolactate synthase|acetohydroxy acid synthase"),
  list("val","ilvC","1.1.1.86","^ilvC$","ketol-acid reductoisomerase"),
  list("val","ilvD","4.2.1.9","^ilvD$","dihydroxy-acid dehydratase"),
  list("val","ilvE","2.6.1.42","^ilvE$","branched.chain.*amino.*transferase"),
  list("leu","leuA","2.3.3.13","^leuA$","2-isopropylmalate synthase"),
  list("leu","leuB","1.1.1.85","^leuB$","3-isopropylmalate dehydrogenase"),
  list("leu","leuC","4.2.1.33","^leuC$","isopropylmalate isomerase.*large"),
  list("leu","leuD","4.2.1.33","^leuD$","isopropylmalate isomerase.*small"),
  list("leu","ilvE","2.6.1.42","^ilvE$","branched.chain.*amino.*transferase"),
  list("ile","ilvA","4.3.1.19","^ilvA$|^tdcB$","threonine dehydratase|threonine deaminase"),
  list("ile","ilvB","2.2.1.6","^ilvB$|^ilvI$|^ilvG$|^ilvH$","acetolactate synthase|acetohydroxy acid synthase"),
  list("ile","ilvC","1.1.1.86","^ilvC$","ketol-acid reductoisomerase"),
  list("ile","ilvD","4.2.1.9","^ilvD$","dihydroxy-acid dehydratase"),
  list("ile","ilvE","2.6.1.42","^ilvE$","branched.chain.*amino.*transferase"),

  # --- Aspartate family ---
  list("thr","thrA","1.1.1.3","^thrA$","aspartate kinase|homoserine dehydrogenase"),
  list("thr","thrB","2.7.1.39","^thrB$","homoserine kinase"),
  list("thr","thrC","4.2.3.1","^thrC$","threonine synthase"),
  list("met","metE","2.1.1.14","^metE$","methionine synthase.*B12-independent"),
  list("met","metH","2.1.1.13","^metH$","methionine synthase.*B12-dependent|cobalamin.*methionine"),
  list("lys","dapA","4.3.3.7","^dapA$","dihydrodipicolinate synthase"),
  list("lys","dapB","1.17.1.8","^dapB$","dihydrodipicolinate reductase"),
  list("lys","lysA","4.1.1.20","^lysA$","diaminopimelate decarboxylase"),
  list("asn","asnA","6.3.1.1","^asnA$","asparagine synthetase.*ammonia"),
  list("asn","asnB","6.3.5.4","^asnB$","asparagine synthetase.*glutamine"),

  # --- Glutamate family ---
  list("gln","glnA","6.3.1.2","^glnA$","glutamine synthetase"),
  list("pro","proA","2.7.2.11","^proA$","glutamate 5-kinase|gamma-glutamyl kinase"),
  list("pro","proB","2.7.2.11","^proB$","glutamate 5-kinase"),
  list("pro","proC","1.5.1.2","^proC$","pyrroline-5-carboxylate reductase"),
  list("arg","argA","2.3.1.1","^argA$","N-acetylglutamate synthase"),
  list("arg","argB","2.7.2.8","^argB$","acetylglutamate kinase"),
  list("arg","argC","1.2.1.38","^argC$","N-acetyl-gamma-glutamyl-phosphate reductase"),
  list("arg","argD","2.6.1.11","^argD$","acetylornithine aminotransferase"),
  list("arg","argF","2.1.3.3","^argF$|^argI$","ornithine carbamoyltransferase"),
  list("arg","argG","6.3.4.5","^argG$","argininosuccinate synthase"),
  list("arg","argH","4.3.2.1","^argH$","argininosuccinate lyase")
 )
 df <- data.frame(
   pathway = vapply(rows, `[[`, "", 1L),
   step    = vapply(rows, `[[`, "", 2L),
   ec      = vapply(rows, `[[`, "", 3L),
   gene_pattern   = vapply(rows, `[[`, "", 4L),
   product_pattern = vapply(rows, `[[`, "", 5L),
   stringsAsFactors = FALSE
 )
 df
}

#' Supplement GapMind step status with GenBank annotation evidence
#'
#' For each pathway step that is missing or low-confidence in GapMind results,
#' search the genbank_table for matching EC numbers, gene names, and product
#' descriptions. Returns additional step assignments with confidence = "annotation".
#'
#' @param genbank_table Data frame with columns: locus_tag, gene, product, EC_number,
#'   start, end, direction, contig.
#' @param step_status Existing GapMind step status data frame (pathway_id, step_id,
#'   confidence, locus_tag). Can be NULL.
#' @return Data frame with same columns as step_status, containing supplementary
#'   annotation-based assignments. Only includes steps not already high/medium.
#' @noRd
.dnmb_gapmind_supplement_by_annotation <- function(genbank_table, step_status = NULL) {
  annot_map <- .dnmb_gapmind_aa_step_annotation_map()
  gt <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)

  # Determine which steps already have good evidence
  good_keys <- character(0)
  if (!is.null(step_status) && nrow(step_status) > 0L) {
    good <- step_status[step_status$confidence %in% c("high", "medium"), , drop = FALSE]
    good_keys <- paste0(good$pathway_id, "::", good$step_id)
  }

  results <- list()
  for (i in seq_len(nrow(annot_map))) {
    pw   <- annot_map$pathway[i]
    step <- annot_map$step[i]
    key  <- paste0(pw, "::", step)
    if (key %in% good_keys) next

    ec_pat   <- annot_map$ec[i]
    gene_pat <- annot_map$gene_pattern[i]
    prod_pat <- annot_map$product_pattern[i]

    # Search by EC number — exact matching (split on | for alternatives,
    # escape dots, add word boundaries to prevent partial matches like
    # EC 1.1.1.3 matching 2.1.1.193)
    ec_hits <- character(0)
    if (nzchar(ec_pat) && "EC_number" %in% names(gt)) {
      ec_col <- gt$EC_number
      ec_col[is.na(ec_col)] <- ""
      ec_alternatives <- strsplit(ec_pat, "\\|")[[1]]
      ec_regex <- paste0("(^|[, ])", gsub("\\.", "\\\\.", ec_alternatives), "($|[, ])")
      ec_regex_combined <- paste(ec_regex, collapse = "|")
      ec_match <- grepl(ec_regex_combined, ec_col)
      ec_hits <- gt$locus_tag[ec_match]
    }

    # Search by gene name
    gene_hits <- character(0)
    if (nzchar(gene_pat) && "gene" %in% names(gt)) {
      gene_col <- gt$gene
      gene_col[is.na(gene_col)] <- ""
      gene_match <- grepl(gene_pat, gene_col, ignore.case = TRUE)
      gene_hits <- gt$locus_tag[gene_match]
    }

    # Search by product description
    prod_hits <- character(0)
    if (nzchar(prod_pat) && "product" %in% names(gt)) {
      prod_col <- gt$product
      prod_col[is.na(prod_col)] <- ""
      prod_match <- grepl(prod_pat, prod_col, ignore.case = TRUE)
      prod_hits <- gt$locus_tag[prod_match]
    }

    all_hits <- unique(c(ec_hits, gene_hits, prod_hits))
    # Remove hits that are not protein-coding (no protein_id)
    if ("protein_id" %in% names(gt) && length(all_hits) > 0L) {
      pid <- gt$protein_id[match(all_hits, gt$locus_tag)]
      all_hits <- all_hits[!is.na(pid) & nzchar(pid)]
    }

    if (length(all_hits) > 0L) {
      # Pick the best hit: prefer gene name match > EC match > product match
      best <- all_hits[1]
      if (length(gene_hits) > 0L) {
        best <- gene_hits[1]
      } else if (length(ec_hits) > 0L) {
        best <- ec_hits[1]
      }

      # Determine evidence source
      sources <- character(0)
      if (best %in% gene_hits) sources <- c(sources, "gene_name")
      if (best %in% ec_hits)   sources <- c(sources, "EC_number")
      if (best %in% prod_hits) sources <- c(sources, "product_desc")

      results[[length(results) + 1L]] <- data.frame(
        pathway_id = pw,
        step_id    = step,
        confidence = "annotation",
        locus_tag  = best,
        evidence_source = paste(sources, collapse = "+"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) == 0L) {
    return(data.frame(
      pathway_id = character(0), step_id = character(0),
      confidence = character(0), locus_tag = character(0),
      evidence_source = character(0), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, results)
}

#' Detect putative operons from genomic coordinates
#'
#' Groups genes into putative operons based on proximity (max intergenic distance)
#' and co-directionality on the same contig.
#'
#' @param genbank_table Data frame with locus_tag, contig, start, end, direction.
#' @param max_gap Maximum intergenic distance in bp to consider genes in the same
#'   operon (default 300 bp).
#' @return A data frame with columns: locus_tag, operon_id (integer).
#' @noRd
.dnmb_detect_operons <- function(genbank_table, max_gap = 300L) {
  gt <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  required <- c("locus_tag", "contig", "start", "end", "direction")
  if (!all(required %in% names(gt))) return(NULL)

  gt$start <- as.numeric(gt$start)
  gt$end   <- as.numeric(gt$end)
  gt <- gt[!is.na(gt$start) & !is.na(gt$end), , drop = FALSE]
  gt <- gt[order(gt$contig, gt$start), , drop = FALSE]

  operon_id <- integer(nrow(gt))
  current_operon <- 1L
  operon_id[1] <- current_operon

  for (i in seq_len(nrow(gt))[-1]) {
    same_contig <- identical(gt$contig[i], gt$contig[i - 1L])
    same_strand <- identical(gt$direction[i], gt$direction[i - 1L])
    gap <- gt$start[i] - gt$end[i - 1L]
    if (same_contig && same_strand && gap >= 0L && gap <= max_gap) {
      operon_id[i] <- current_operon
    } else {
      current_operon <- current_operon + 1L
      operon_id[i] <- current_operon
    }
  }

  data.frame(locus_tag = gt$locus_tag, operon_id = operon_id, stringsAsFactors = FALSE)
}

#' Boost annotation evidence using operon context
#'
#' For annotation-supplemented steps, check if the candidate gene shares an
#' operon with other genes in the same pathway (from GapMind or annotation
#' evidence). Genes in a pathway operon get confidence upgraded to "operon".
#'
#' @param step_status Combined step status (GapMind + annotation supplement).
#' @param genbank_table Data frame with genomic coordinates.
#' @param max_gap Maximum intergenic distance for operon detection.
#' @return Updated step_status with operon-boosted confidence.
#' @noRd
.dnmb_gapmind_operon_boost <- function(step_status, genbank_table, max_gap = 300L) {
  if (is.null(step_status) || nrow(step_status) == 0L) return(step_status)

  operons <- .dnmb_detect_operons(genbank_table, max_gap = max_gap)
  if (is.null(operons)) return(step_status)

  # Merge operon IDs onto step status
  ss <- merge(step_status, operons, by = "locus_tag", all.x = TRUE)

  # For each pathway, check if multiple steps share an operon
  for (pw in unique(ss$pathway_id)) {
    pw_rows <- which(ss$pathway_id == pw & !is.na(ss$operon_id))
    if (length(pw_rows) < 2L) next

    pw_operons <- ss$operon_id[pw_rows]
    # Find operons with 2+ pathway genes
    op_counts <- table(pw_operons)
    clustered_ops <- names(op_counts)[op_counts >= 2L]

    if (length(clustered_ops) == 0L) next

    # Boost annotation-only or low confidence genes in these operons
    for (r in pw_rows) {
      if (ss$operon_id[r] %in% as.integer(clustered_ops) &&
          ss$confidence[r] %in% c("annotation", "low")) {
        ss$confidence[r] <- "operon"
      }
    }
  }

  ss$operon_id <- NULL
  ss
}

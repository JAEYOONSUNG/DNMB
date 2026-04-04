.dnmb_prophage_module_name <- function() {
  "prophage"
}

.dnmb_prophage_default_backend <- function() {
  "phispy"
}

.dnmb_prophage_supported_backends <- function() {
  c("phispy", "virsorter2", "pide")
}

.dnmb_prophage_normalize_backend <- function(backend = .dnmb_prophage_default_backend()) {
  backend <- tolower(trimws(as.character(backend)[1]))
  if (!nzchar(backend)) {
    backend <- .dnmb_prophage_default_backend()
  }
  if (backend %in% c("phispy", "phi_spy", "phi-spy")) {
    return("phispy")
  }
  if (backend %in% c("virsorter2", "virsorter", "vs2")) {
    return("virsorter2")
  }
  if (backend %in% c("pide")) {
    return("pide")
  }
  stop(
    "Unsupported prophage backend: ",
    backend,
    ". Supported backends: ",
    paste(.dnmb_prophage_supported_backends(), collapse = ", "),
    call. = FALSE
  )
}

.dnmb_prophage_default_version <- function(backend = .dnmb_prophage_default_backend()) {
  .dnmb_prophage_normalize_backend(backend)
}

.dnmb_prophage_default_randomforest_trees <- function() {
  100L
}

.dnmb_prophage_default_repo_url <- function() {
  "https://github.com/linsalrob/PhiSpy.git"
}

.dnmb_prophage_default_virsorter2_repo_url <- function() {
  "https://github.com/jiarong/VirSorter2.git"
}

.dnmb_prophage_default_pide_repo_url <- function() {
  "https://github.com/chyghy/PIDE.git"
}

.dnmb_prophage_default_pide_model_url <- function() {
  "https://zenodo.org/records/12759619/files/PIDE.model.tar.gz"
}

.dnmb_prophage_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_prophage_empty_status <- function() {
  .dnmb_prophage_status_row(character(), character(), character())
}

.dnmb_prophage_trace <- function(path, text) {
  base::cat(base::paste0(text, "\n"), file = path, append = TRUE)
}

.dnmb_prophage_asset_layout <- function(module_dir) {
  list(
    module_dir = module_dir,
    repo_dir = base::file.path(module_dir, "PhiSpy"),
    env_dir = base::file.path(module_dir, "venv"),
    env_python = base::file.path(module_dir, "venv", "bin", "python"),
    env_pip = base::file.path(module_dir, "venv", "bin", "pip"),
    cli_path = base::file.path(module_dir, "venv", "bin", "PhiSpy.py"),
    install_trace_log = base::file.path(module_dir, "prophage_install_trace.log")
  )
}

.dnmb_prophage_resolve_cli <- function(layout) {
  candidates <- c(layout$cli_path, dnmb_detect_binary("PhiSpy.py", required = FALSE)$path)
  candidates <- unique(candidates[base::nzchar(candidates)])
  for (candidate in candidates) {
    if (!base::file.exists(candidate)) {
      next
    }
    run <- dnmb_run_external(candidate, args = "--help", required = FALSE)
    if (base::isTRUE(run$ok)) {
      return(candidate)
    }
  }
  ""
}

.dnmb_prophage_candidate_python <- function() {
  # Prefer Homebrew python on Apple Silicon to avoid conda x86 compiler issues
  homebrew_candidates <- c(
    "/opt/homebrew/bin/python3",
    "/opt/homebrew/bin/python3.13",
    "/opt/homebrew/bin/python3.12",
    "/opt/homebrew/bin/python3.11"
  )
  for (hb in homebrew_candidates) {
    if (base::file.exists(hb)) return(hb)
  }
  candidates <- c(
    dnmb_detect_binary("python3.12", required = FALSE)$path,
    dnmb_detect_binary("python3.11", required = FALSE)$path,
    dnmb_detect_binary("python3.13", required = FALSE)$path,
    dnmb_detect_binary("python3.10", required = FALSE)$path,
    dnmb_detect_binary("python3", required = FALSE)$path
  )
  candidates <- unique(candidates[nzchar(candidates)])
  if (!base::length(candidates)) return("")
  candidates[[1]]
}

.dnmb_prophage_candidate_conda <- function() {
  conda <- dnmb_detect_binary("conda", required = FALSE)
  if (isTRUE(conda$found)) {
    return(conda$path)
  }
  ""
}

.dnmb_prophage_no_space_root <- function() {
  root <- base::file.path(path.expand("~"), ".dnmb-external")
  base::dir.create(root, recursive = TRUE, showWarnings = FALSE)
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

.dnmb_prophage_resolve_virsorter2_db_dir <- function(db_dir) {
  db_dir <- normalizePath(db_dir, winslash = "/", mustWork = FALSE)
  nested_dir <- base::file.path(db_dir, "db")
  candidates <- c(db_dir, nested_dir)
  for (candidate in candidates) {
    if (base::dir.exists(base::file.path(candidate, "group")) &&
        base::dir.exists(base::file.path(candidate, "hmm"))) {
      return(candidate)
    }
  }
  db_dir
}

.dnmb_prophage_hmm_file_complete <- function(path) {
  if (!base::file.exists(path)) {
    return(FALSE)
  }
  size <- suppressWarnings(base::file.info(path)$size)
  if (!base::is.finite(size) || size < 1024) {
    return(FALSE)
  }
  con <- base::file(path, open = "rb")
  on.exit(base::close(con), add = TRUE)
  head_raw <- base::readBin(con, what = "raw", n = 256L)
  if (!length(head_raw)) {
    return(FALSE)
  }
  head_txt <- rawToChar(head_raw)
  if (!grepl("^HMMER3/f", head_txt)) {
    return(FALSE)
  }
  seek_from_end <- max(0, size - 4096)
  base::seek(con, where = seek_from_end, origin = "start")
  tail_raw <- base::readBin(con, what = "raw", n = 4096L)
  if (!length(tail_raw)) {
    return(FALSE)
  }
  tail_txt <- rawToChar(tail_raw)
  grepl("//\\s*$", tail_txt)
}

.dnmb_prophage_virsorter2_db_ready <- function(db_dir) {
  effective_db_dir <- .dnmb_prophage_resolve_virsorter2_db_dir(db_dir)
  hmm_path <- base::file.path(effective_db_dir, "hmm", "viral", "combined.hmm")
  ready <- base::dir.exists(base::file.path(effective_db_dir, "group")) &&
    base::dir.exists(base::file.path(effective_db_dir, "hmm")) &&
    .dnmb_prophage_hmm_file_complete(hmm_path)
  list(ok = ready, db_dir = effective_db_dir, hmm_path = hmm_path)
}

.dnmb_prophage_reset_virsorter2_db_extract <- function(db_dir) {
  db_dir <- normalizePath(db_dir, winslash = "/", mustWork = FALSE)
  base::unlink(base::file.path(db_dir, "db"), recursive = TRUE, force = TRUE)
  for (subdir in c("group", "hmm", "rbs", "conda_envs")) {
    base::unlink(base::file.path(db_dir, subdir), recursive = TRUE, force = TRUE)
  }
  invisible(db_dir)
}

.dnmb_prophage_virsorter2_setup_assets_ready <- function(db_dir) {
  required <- c(
    "db.tgz",
    "Pfam-A-acc2desc.tsv",
    "Pfam-A-Archaea.hmm",
    "Pfam-A-Bacteria.hmm",
    "Pfam-A-Eukaryota.hmm",
    "Pfam-A-Mixed.hmm",
    "Pfam-A-Viruses.hmm",
    "combined.hmm.gz.split_00",
    "combined.hmm.gz.split_01",
    "combined.hmm.gz.split_02"
  )
  base::all(base::file.exists(base::file.path(db_dir, required)))
}

.dnmb_prophage_virsorter2_setup_asset_specs <- function() {
  zenodo <- "https://zenodo.org/record/4269607/files"
  specs <- list(
    list(dest = "db.tgz", url = paste0(zenodo, "/db.tgz?download=1"), compressed = FALSE),
    list(dest = "combined.hmm.gz.split_00", url = paste0(zenodo, "/combined.hmm.gz.split_00?download=1"), compressed = FALSE),
    list(dest = "combined.hmm.gz.split_01", url = paste0(zenodo, "/combined.hmm.gz.split_01?download=1"), compressed = FALSE),
    list(dest = "combined.hmm.gz.split_02", url = paste0(zenodo, "/combined.hmm.gz.split_02?download=1"), compressed = FALSE),
    list(dest = "Pfam-A-acc2desc.tsv", url = paste0(zenodo, "/Pfam-A-acc2desc.tsv.gz?download=1"), compressed = TRUE),
    list(dest = "Pfam-A-Archaea.hmm", url = paste0(zenodo, "/Pfam-A-Archaea.hmm.gz?download=1"), compressed = TRUE),
    list(dest = "Pfam-A-Bacteria.hmm", url = paste0(zenodo, "/Pfam-A-Bacteria.hmm.gz?download=1"), compressed = TRUE),
    list(dest = "Pfam-A-Eukaryota.hmm", url = paste0(zenodo, "/Pfam-A-Eukaryota.hmm.gz?download=1"), compressed = TRUE),
    list(dest = "Pfam-A-Mixed.hmm", url = paste0(zenodo, "/Pfam-A-Mixed.hmm.gz?download=1"), compressed = TRUE),
    list(dest = "Pfam-A-Viruses.hmm", url = paste0(zenodo, "/Pfam-A-Viruses.hmm.gz?download=1"), compressed = TRUE)
  )
  stats::setNames(specs, vapply(specs, `[[`, "", "dest"))
}

.dnmb_prophage_download_virsorter2_setup_assets <- function(db_dir) {
  specs <- .dnmb_prophage_virsorter2_setup_asset_specs()
  base::dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)
  for (spec in specs) {
    dest <- base::file.path(db_dir, spec$dest)
    if (base::file.exists(dest) && base::file.info(dest)$size > 0) {
      next
    }
    if (isTRUE(spec$compressed)) {
      gz_dest <- paste0(dest, ".gz")
      dl <- .dnmb_download_asset(spec$url, gz_dest, insecure = FALSE)
      if (!base::isTRUE(dl$ok)) {
        return(list(ok = FALSE, detail = dl$error %||% spec$url))
      }
      unzip_run <- dnmb_run_external("gunzip", args = c("-f", gz_dest), required = FALSE)
      if (!base::isTRUE(unzip_run$ok) || !base::file.exists(dest)) {
        return(list(ok = FALSE, detail = unzip_run$error %||% paste("Failed to gunzip", gz_dest)))
      }
    } else {
      dl <- .dnmb_download_asset(spec$url, dest, insecure = FALSE)
      if (!base::isTRUE(dl$ok)) {
        return(list(ok = FALSE, detail = dl$error %||% spec$url))
      }
    }
  }
  list(ok = TRUE, detail = db_dir)
}

.dnmb_prophage_finalize_virsorter2_db <- function(layout) {
  if (!.dnmb_prophage_virsorter2_setup_assets_ready(layout$db_dir)) {
    return(list(ok = FALSE, detail = "VirSorter2 setup assets are incomplete."))
  }
  .dnmb_prophage_reset_virsorter2_db_extract(layout$db_dir)
  script <- paste(
    "set -euo pipefail",
    paste("cd", shQuote(layout$db_dir)),
    "rm -rf group hmm rbs db",
    "tar -xzf db.tgz",
    "mv db/* .",
    "rmdir db",
    "mkdir -p hmm/pfam hmm/viral",
    "mv Pfam-A-*.hmm hmm/pfam",
    "mv Pfam-A-acc2desc.tsv hmm/pfam",
    "cat combined.hmm.gz.split_* | gunzip -c > hmm/viral/combined.hmm",
    "chmod -R 755 group hmm rbs",
    sep = " && "
  )
  run <- dnmb_run_external("bash", args = c("-lc", script), required = FALSE)
  if (!base::isTRUE(run$ok)) {
    return(list(ok = FALSE, detail = run$error %||% "VirSorter2 DB finalize failed."))
  }
  db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
  if (!isTRUE(db_check$ok)) {
    return(list(ok = FALSE, detail = "VirSorter2 DB finalize completed but DB validation failed."))
  }
  config_run <- dnmb_run_external(
    layout$cli_path,
    args = c("config", "--init-source", "--db-dir", db_check$db_dir),
    env = c(
      PATH = paste(unique(c(layout$env_bin, Sys.getenv("PATH"))), collapse = .Platform$path.sep),
      CONDA_PKGS_DIRS = layout$pkgs_dir
    ),
    required = FALSE
  )
  if (!base::isTRUE(config_run$ok)) {
    return(list(ok = FALSE, detail = config_run$error %||% "VirSorter2 config update failed."))
  }
  list(ok = TRUE, detail = db_check$db_dir)
}

.dnmb_prophage_parse_genbank_records <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  lines <- base::readLines(path, warn = FALSE)
  if (!base::length(lines)) {
    return(data.frame())
  }
  loci <- base::grep("^LOCUS", lines)
  if (!base::length(loci)) {
    return(data.frame())
  }
  ends <- c(loci[-1L] - 1L, base::length(lines))
  records <- lapply(seq_along(loci), function(i) {
    rec <- lines[loci[[i]]:ends[[i]]]
    locus_id <- sub("^LOCUS\\s+([^ ]+).*$", "\\1", rec[[1]])
    acc_idx <- base::grep("^ACCESSION", rec)
    accession <- if (base::length(acc_idx)) {
      trimws(sub("^ACCESSION\\s+", "", rec[[acc_idx[[1]]]]))
    } else {
      locus_id
    }
    accession <- strsplit(accession, "\\s+")[[1]][[1]]
    def_idx <- base::grep("^DEFINITION", rec)
    definition <- NA_character_
    if (base::length(def_idx)) {
      definition <- sub("^DEFINITION\\s+", "", rec[[def_idx[[1]]]])
      j <- def_idx[[1]] + 1L
      while (j <= base::length(rec) && grepl("^\\s{12}\\S", rec[[j]])) {
        definition <- paste(definition, trimws(rec[[j]]))
        j <- j + 1L
      }
      definition <- gsub("\\.$", "", trimws(gsub("\\s+", " ", definition)))
    }
    origin_idx <- base::grep("^ORIGIN", rec)
    end_idx <- base::grep("^//", rec)
    sequence <- NA_character_
    if (base::length(origin_idx) && base::length(end_idx) && end_idx[[1]] > origin_idx[[1]]) {
      sequence <- paste(rec[(origin_idx[[1]] + 1L):(end_idx[[1]] - 1L)], collapse = "")
      sequence <- gsub("[^ACGTacgt]", "", sequence)
      sequence <- toupper(sequence)
    }
    data.frame(
      locus = locus_id,
      accession = accession,
      definition = definition,
      sequence = sequence,
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(records)
  rownames(out) <- NULL
  out
}

.dnmb_prophage_write_contig_fasta <- function(records, out_fasta) {
  if (!base::nrow(records)) {
    return(out_fasta)
  }
  con <- base::file(out_fasta, open = "w")
  on.exit(base::close(con), add = TRUE)
  for (i in seq_len(nrow(records))) {
    if (is.na(records$sequence[[i]]) || !nzchar(records$sequence[[i]])) {
      next
    }
    base::writeLines(paste0(">", records$accession[[i]]), con = con)
    base::writeLines(records$sequence[[i]], con = con)
  }
  out_fasta
}

.dnmb_prophage_match_contig_label <- function(contig_id, records, genes_tbl) {
  contig_id <- as.character(contig_id)[1]
  genes_contigs <- unique(as.character(genes_tbl$contig))
  if (contig_id %in% genes_contigs) {
    return(contig_id)
  }
  if (base::nrow(records)) {
    idx <- match(contig_id, records$accession)
    if (!is.na(idx) && !is.na(records$definition[[idx]]) && records$definition[[idx]] %in% genes_contigs) {
      return(records$definition[[idx]])
    }
    idx <- match(contig_id, records$locus)
    if (!is.na(idx) && !is.na(records$definition[[idx]]) && records$definition[[idx]] %in% genes_contigs) {
      return(records$definition[[idx]])
    }
  }
  contig_id
}

.dnmb_prophage_hits_from_regions <- function(region_tbl,
                                             genes,
                                             genbank_records = data.frame(),
                                             backend = .dnmb_prophage_default_backend()) {
  if (base::is.null(region_tbl) || !base::is.data.frame(region_tbl) || !base::nrow(region_tbl)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  backend <- .dnmb_prophage_normalize_backend(backend)
  genes_tbl <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  genes_tbl$locus_tag <- .dnmb_module_clean_annotation_key(genes_tbl$locus_tag)
  genes_tbl$start <- suppressWarnings(base::as.numeric(genes_tbl$start))
  genes_tbl$end <- suppressWarnings(base::as.numeric(genes_tbl$end))
  rows <- lapply(seq_len(nrow(region_tbl)), function(i) {
    contig_label <- .dnmb_prophage_match_contig_label(region_tbl$contig[[i]], genbank_records, genes_tbl)
    keep <- genes_tbl$contig == contig_label &
      !is.na(genes_tbl$start) &
      !is.na(genes_tbl$end) &
      pmin(genes_tbl$start, genes_tbl$end) >= region_tbl$prophage_start[[i]] &
      pmax(genes_tbl$start, genes_tbl$end) <= region_tbl$prophage_end[[i]]
    if (!any(keep)) {
      return(NULL)
    }
    support_bits <- c(
      paste0("backend=", backend),
      if ("prophage_score" %in% names(region_tbl) && !is.na(region_tbl$prophage_score[[i]])) paste0("score=", region_tbl$prophage_score[[i]]) else NULL,
      if ("prophage_group" %in% names(region_tbl) && !is.na(region_tbl$prophage_group[[i]])) paste0("group=", region_tbl$prophage_group[[i]]) else NULL,
      if ("partial" %in% names(region_tbl) && !is.na(region_tbl$partial[[i]])) paste0("partial=", region_tbl$partial[[i]]) else NULL,
      if ("hallmark_cnt" %in% names(region_tbl) && !is.na(region_tbl$hallmark_cnt[[i]])) paste0("hallmark_cnt=", region_tbl$hallmark_cnt[[i]]) else NULL,
      paste0("region=", region_tbl$prophage_start[[i]], "-", region_tbl$prophage_end[[i]])
    )
    data.frame(
      query = genes_tbl$locus_tag[keep],
      source = "prophage",
      family_system = backend,
      family_id = as.character(region_tbl$prophage_id[[i]]),
      hit_label = as.character(region_tbl$prophage_id[[i]]),
      enzyme_role = "prophage_region",
      evidence_mode = backend,
      substrate_label = if ("prophage_group" %in% names(region_tbl)) as.character(region_tbl$prophage_group[[i]]) else NA_character_,
      support = paste(support_bits, collapse = "; "),
      typing_eligible = TRUE,
      prophage_id = as.character(region_tbl$prophage_id[[i]]),
      prophage_rank = if ("prophage_score" %in% names(region_tbl)) suppressWarnings(as.numeric(region_tbl$prophage_score[[i]])) else NA_real_,
      prophage_my_status = NA_real_,
      prophage_pp = if ("hallmark_cnt" %in% names(region_tbl)) suppressWarnings(as.numeric(region_tbl$hallmark_cnt[[i]])) else NA_real_,
      prophage_start = suppressWarnings(as.numeric(region_tbl$prophage_start[[i]])),
      prophage_end = suppressWarnings(as.numeric(region_tbl$prophage_end[[i]])),
      prophage_backend = backend,
      prophage_score = if ("prophage_score" %in% names(region_tbl)) suppressWarnings(as.numeric(region_tbl$prophage_score[[i]])) else NA_real_,
      prophage_group = if ("prophage_group" %in% names(region_tbl)) as.character(region_tbl$prophage_group[[i]]) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(base::Filter(base::Negate(base::is.null), rows))
  if (!base::nrow(out)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  out
}

.dnmb_prophage_copy_dir_contents <- function(src_dir, dest_dir) {
  base::dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  entries <- base::list.files(src_dir, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  if (!base::length(entries)) {
    return(invisible(dest_dir))
  }
  ok <- base::file.copy(entries, dest_dir, recursive = TRUE, copy.mode = TRUE)
  if (!base::all(ok)) {
    base::stop("Failed to copy PhiSpy assets from ", src_dir, " to ", dest_dir, call. = FALSE)
  }
  invisible(dest_dir)
}

.dnmb_prophage_prepare_repo <- function(layout,
                                        repo_source = .dnmb_prophage_default_repo_url(),
                                        force = FALSE) {
  repo_source <- base::as.character(repo_source)[1]
  if (base::dir.exists(layout$repo_dir) && !base::isTRUE(force)) {
    run <- dnmb_run_external("git", args = c("-C", layout$repo_dir, "pull", "--ff-only"), required = FALSE)
    return(.dnmb_prophage_status_row(
      "prophage_repo",
      if (base::isTRUE(run$ok)) "updated" else "cached",
      if (base::isTRUE(run$ok)) layout$repo_dir else (run$error %||% layout$repo_dir)
    ))
  }

  if (base::dir.exists(layout$repo_dir) && base::isTRUE(force)) {
    unlink(layout$repo_dir, recursive = TRUE, force = TRUE)
  }

  if (base::dir.exists(repo_source)) {
    .dnmb_prophage_copy_dir_contents(repo_source, layout$repo_dir)
    return(.dnmb_prophage_status_row("prophage_repo", "ok", layout$repo_dir))
  }

  run <- dnmb_run_external(
    "git",
    args = c("clone", "--depth", "1", repo_source, layout$repo_dir),
    env = c(GIT_LFS_SKIP_SMUDGE = "1"),
    required = FALSE
  )
  .dnmb_prophage_status_row(
    "prophage_repo",
    if (base::isTRUE(run$ok)) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    if (base::isTRUE(run$ok)) layout$repo_dir else (run$error %||% repo_source)
  )
}

.dnmb_prophage_prepare_env <- function(layout, force = FALSE) {
  py <- .dnmb_prophage_candidate_python()
  if (!base::nzchar(py)) {
    return(.dnmb_prophage_status_row("prophage_python", "missing", "No python3 executable found in PATH."))
  }
  if (base::dir.exists(layout$env_dir) && base::isTRUE(force)) {
    unlink(layout$env_dir, recursive = TRUE, force = TRUE)
  }
  # Check if existing venv works (may be from different platform/python version)
  if (base::dir.exists(layout$env_dir)) {
    if (!base::file.exists(layout$env_python)) {
      unlink(layout$env_dir, recursive = TRUE, force = TRUE)
    } else {
      test <- dnmb_run_external(layout$env_python, args = c("-c", "print('ok')"), required = FALSE)
      if (!base::isTRUE(test$ok)) {
        unlink(layout$env_dir, recursive = TRUE, force = TRUE)
      }
    }
  }
  if (!base::file.exists(layout$env_python)) {
    run <- dnmb_run_external(py, args = c("-m", "venv", layout$env_dir), required = FALSE)
    if (!base::isTRUE(run$ok)) {
      return(.dnmb_prophage_status_row("prophage_python", "failed", run$error %||% py))
    }
  }
  .dnmb_prophage_status_row("prophage_python", "ok", layout$env_python)
}

.dnmb_prophage_install_cli <- function(layout) {
  # Clean compiler env to prevent conda x86 flags from poisoning arm64 builds
  clean_env <- c(
    CC = "/usr/bin/cc", CXX = "/usr/bin/c++",
    CFLAGS = "", CXXFLAGS = "", LDFLAGS = "",
    ARCHFLAGS = paste0("-arch ", Sys.info()[["machine"]])
  )
  # pip upgrade is optional — don't fail if it doesn't work
  dnmb_run_external(layout$env_python, args = c("-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"), env = clean_env, required = FALSE)
  # Install PhiSpy
  run <- dnmb_run_external(layout$env_pip, args = c("install", "--no-cache-dir", "-U", layout$repo_dir), env = clean_env, required = FALSE)
  if (!base::isTRUE(run$ok)) {
    return(.dnmb_prophage_status_row("prophage_cli", "failed", run$error %||% layout$env_pip))
  }
  .dnmb_prophage_status_row(
    "prophage_cli",
    if (base::file.exists(layout$cli_path)) "ok" else "failed",
    if (base::file.exists(layout$cli_path)) layout$cli_path else layout$env_dir
  )
}

dnmb_prophage_install_module <- function(version = .dnmb_prophage_default_version(),
                                         cache_root = NULL,
                                         install = TRUE,
                                         repo_url = .dnmb_prophage_default_repo_url(),
                                         asset_urls = NULL,
                                         force = FALSE) {
  module <- .dnmb_prophage_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_prophage_asset_layout(module_dir)
  status <- .dnmb_prophage_empty_status()

  if (!base::isTRUE(install)) {
    manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
    cli_path <- .dnmb_prophage_resolve_cli(layout)
    ready <- !base::is.null(manifest) && base::isTRUE(manifest$install_ok) && base::nzchar(cli_path)
    if (ready) {
      return(list(ok = TRUE, status = .dnmb_prophage_status_row("prophage_install", "cached", module_dir), files = list(cli = cli_path, env_python = layout$env_python), manifest = manifest))
    }
    return(list(ok = FALSE, status = .dnmb_prophage_status_row("prophage_install", "missing", "PhiSpy is missing and module_install is FALSE."), files = list(), manifest = NULL))
  }

  repo_status <- .dnmb_prophage_prepare_repo(layout, repo_source = repo_url, force = force)
  status <- dplyr::bind_rows(status, repo_status)
  if (!repo_status$status %in% c("ok", "cached", "updated")) {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }
  env_status <- .dnmb_prophage_prepare_env(layout, force = force)
  status <- dplyr::bind_rows(status, env_status)
  if (env_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }
  cli_status <- .dnmb_prophage_install_cli(layout)
  status <- dplyr::bind_rows(status, cli_status)
  if (cli_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = NULL))
  }

  manifest <- list(
    install_ok = TRUE,
    repo_dir = layout$repo_dir,
    env_python = layout$env_python,
    cli_path = layout$cli_path
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_prophage_status_row("prophage_install", "ok", module_dir)),
    files = list(cli = layout$cli_path, env_python = layout$env_python),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_prophage_get_module <- function(version = .dnmb_prophage_default_version(),
                                     cache_root = NULL,
                                     required = TRUE) {
  manifest <- dnmb_db_read_manifest(.dnmb_prophage_module_name(), version, cache_root = cache_root, required = FALSE)
  if (base::is.null(manifest)) {
    if (base::isTRUE(required)) {
      base::stop("PhiSpy module is not installed.", call. = FALSE)
    }
    return(list(ok = FALSE, manifest = NULL))
  }
  layout <- .dnmb_prophage_asset_layout(.dnmb_db_module_dir(.dnmb_prophage_module_name(), version, cache_root = cache_root, create = FALSE))
  cli_path <- .dnmb_prophage_resolve_cli(layout)
  list(
    ok = TRUE,
    manifest = manifest,
    files = list(
      cli = cli_path,
      env_python = layout$env_python,
      repo_dir = layout$repo_dir
    )
  )
}

dnmb_prophage_parse_coordinates <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  out <- tryCatch(
    utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  )
  if (!base::is.data.frame(out) || !base::nrow(out)) {
    out <- tryCatch(
      utils::read.delim(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) data.frame()
    )
  }
  out
}

dnmb_prophage_parse_information <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

.dnmb_prophage_standardize_coordinates <- function(tbl) {
  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !base::nrow(tbl)) {
    return(data.frame())
  }
  out <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  names_lc <- tolower(names(out))
  names(out) <- names_lc
  data.frame(
    prophage_id = base::as.character(out[[1]]),
    contig = base::as.character(out[[2]]),
    prophage_start = suppressWarnings(base::as.numeric(out[[3]])),
    prophage_end = suppressWarnings(base::as.numeric(out[[4]])),
    stringsAsFactors = FALSE
  )
}

.dnmb_prophage_standardize_information <- function(tbl) {
  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !base::nrow(tbl)) {
    return(data.frame())
  }
  out <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  names_lc <- gsub("[^a-z0-9]+", "_", tolower(names(out)))
  names(out) <- names_lc
  col_or_na <- function(name) if (name %in% names(out)) out[[name]] else base::rep(NA, nrow(out))
  data.frame(
    gene_id = base::as.character(col_or_na("identifier")),
    product = base::as.character(col_or_na("function")),
    contig = base::as.character(col_or_na("contig")),
    gene_start = suppressWarnings(base::as.numeric(col_or_na("start"))),
    gene_end = suppressWarnings(base::as.numeric(col_or_na("stop"))),
    gene_position = suppressWarnings(base::as.numeric(col_or_na("position"))),
    rank = suppressWarnings(base::as.numeric(col_or_na("rank"))),
    my_status = suppressWarnings(base::as.numeric(col_or_na("my_status"))),
    pp = suppressWarnings(base::as.numeric(col_or_na("pp"))),
    final_status = base::as.character(col_or_na("final_status")),
    attl_start = suppressWarnings(base::as.numeric(col_or_na("start_of_attl"))),
    attl_end = suppressWarnings(base::as.numeric(col_or_na("end_of_attl"))),
    attr_start = suppressWarnings(base::as.numeric(col_or_na("start_of_attr"))),
    attr_end = suppressWarnings(base::as.numeric(col_or_na("end_of_attr"))),
    stringsAsFactors = FALSE
  )
}

dnmb_prophage_gene_hits <- function(info_tbl, coord_tbl, genes) {
  if ((!base::is.data.frame(info_tbl) || !base::nrow(info_tbl)) && (!base::is.data.frame(coord_tbl) || !base::nrow(coord_tbl))) {
    return(data.frame())
  }
  genes_tbl <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  genes_tbl$locus_tag <- .dnmb_module_clean_annotation_key(genes_tbl$locus_tag)
  genes_tbl$start <- suppressWarnings(base::as.numeric(genes_tbl$start))
  genes_tbl$end <- suppressWarnings(base::as.numeric(genes_tbl$end))
  if (base::is.data.frame(info_tbl) && base::nrow(info_tbl)) {
    info_tbl$gene_id <- .dnmb_module_clean_annotation_key(info_tbl$gene_id)
    by_locus <- dplyr::left_join(
      genes_tbl,
      info_tbl,
      by = c("locus_tag" = "gene_id")
    )
    if ("contig.x" %in% names(by_locus) && !"contig" %in% names(by_locus)) {
      by_locus$contig <- by_locus$contig.x
    }
    if ("start.x" %in% names(by_locus) && !"start" %in% names(by_locus)) {
      by_locus$start <- by_locus$start.x
    }
    if ("end.x" %in% names(by_locus) && !"end" %in% names(by_locus)) {
      by_locus$end <- by_locus$end.x
    }
    if ("final_status" %in% names(by_locus) && any(!is.na(by_locus$final_status) & by_locus$final_status != "0")) {
      if (base::is.data.frame(coord_tbl) && base::nrow(coord_tbl)) {
        by_locus <- dplyr::left_join(by_locus, coord_tbl, by = c("final_status" = "prophage_id", "contig" = "contig"))
      }
      return(by_locus)
    }
    joined <- dplyr::left_join(
      genes_tbl,
      info_tbl,
      by = c("contig" = "contig", "start" = "gene_start", "end" = "gene_end")
    )
    if (base::is.data.frame(coord_tbl) && base::nrow(coord_tbl)) {
      joined <- dplyr::left_join(joined, coord_tbl, by = c("final_status" = "prophage_id", "contig" = "contig"))
    }
    return(joined)
  }
  out <- genes_tbl
  out$prophage_id <- NA_character_
  out$prophage_start <- NA_real_
  out$prophage_end <- NA_real_
  for (i in seq_len(nrow(coord_tbl))) {
    keep <- out$contig == coord_tbl$contig[[i]] &
      suppressWarnings(base::as.numeric(out$start)) >= coord_tbl$prophage_start[[i]] &
      suppressWarnings(base::as.numeric(out$end)) <= coord_tbl$prophage_end[[i]]
    out$prophage_id[keep] <- coord_tbl$prophage_id[[i]]
    out$prophage_start[keep] <- coord_tbl$prophage_start[[i]]
    out$prophage_end[keep] <- coord_tbl$prophage_end[[i]]
  }
  out
}

dnmb_prophage_normalize_hits <- function(gene_hits, backend = "phispy") {
  backend <- .dnmb_prophage_normalize_backend(backend)
  if (base::is.null(gene_hits) || !base::is.data.frame(gene_hits) || !base::nrow(gene_hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  tbl <- base::as.data.frame(gene_hits, stringsAsFactors = FALSE)
  for (column_name in c("final_status", "rank", "my_status", "pp", "prophage_id", "prophage_start", "prophage_end")) {
    if (base::is.null(tbl[[column_name]])) {
      tbl[[column_name]] <- base::rep(NA, base::nrow(tbl))
    }
  }
  tbl$final_status <- trimws(base::as.character(tbl$final_status))
  tbl$prophage_id <- trimws(base::as.character(tbl$prophage_id))
  missing_region <- base::is.na(tbl$prophage_id) | !base::nzchar(base::as.character(tbl$prophage_id))
  tbl$prophage_id[missing_region] <- base::as.character(tbl$final_status[missing_region])
  tbl$query <- .dnmb_module_clean_annotation_key(tbl$locus_tag)
  tbl$has_prophage <- !base::is.na(tbl$final_status) & base::nzchar(tbl$final_status) & tbl$final_status != "0"
  if ("prophage_id" %in% names(tbl)) {
    tbl$has_prophage <- tbl$has_prophage | (!base::is.na(tbl$prophage_id) & base::nzchar(tbl$prophage_id) & tbl$prophage_id != "0")
  }
  tbl <- tbl[tbl$has_prophage, , drop = FALSE]
  if (!base::nrow(tbl)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  prophage_label <- if ("prophage_id" %in% names(tbl)) base::as.character(tbl$prophage_id) else base::as.character(tbl$final_status)
  data.frame(
    query = tbl$query,
    source = "prophage",
    family_system = backend,
    family_id = prophage_label,
    hit_label = prophage_label,
    enzyme_role = "prophage_region",
    evidence_mode = backend,
    substrate_label = NA_character_,
    support = base::paste0(
      "rank=", base::as.character(tbl$rank),
      "; my_status=", base::as.character(tbl$my_status),
      "; pp=", base::as.character(tbl$pp),
      "; region=", ifelse(base::is.na(tbl$prophage_start), "NA", base::as.character(tbl$prophage_start)),
      "-",
      ifelse(base::is.na(tbl$prophage_end), "NA", base::as.character(tbl$prophage_end))
    ),
    typing_eligible = TRUE,
    prophage_id = prophage_label,
    prophage_rank = suppressWarnings(base::as.numeric(tbl$rank)),
    prophage_my_status = suppressWarnings(base::as.numeric(tbl$my_status)),
    prophage_pp = suppressWarnings(base::as.numeric(tbl$pp)),
    prophage_start = suppressWarnings(base::as.numeric(tbl$prophage_start)),
    prophage_end = suppressWarnings(base::as.numeric(tbl$prophage_end)),
    prophage_backend = backend,
    prophage_score = NA_real_,
    prophage_group = NA_character_,
    stringsAsFactors = FALSE
  )
}

.dnmb_prophage_output_table <- function(genes, hits) {
  out <- .dnmb_module_output_table(genes = genes, hits = hits)
  drop_cols <- base::intersect(c("hit_label", "enzyme_role", "evidence_mode", "substrate_label", "typing_eligible"), base::names(out))
  if (base::length(drop_cols)) {
    out[drop_cols] <- NULL
  }
  ordered <- c(
    base::intersect(dnmb_backbone_columns(), base::names(out)),
    base::intersect(c("prophage_id", "prophage_rank", "prophage_my_status", "prophage_pp", "prophage_start", "prophage_end", "prophage_backend", "prophage_score", "prophage_group", "support"), base::names(out)),
    base::setdiff(base::names(out), c(base::intersect(dnmb_backbone_columns(), base::names(out)), "prophage_id", "prophage_rank", "prophage_my_status", "prophage_pp", "prophage_start", "prophage_end", "prophage_backend", "prophage_score", "prophage_group", "support"))
  )
  out[, ordered, drop = FALSE]
}

.dnmb_prophage_virsorter2_asset_layout <- function(module_dir) {
  external_root <- base::file.path(.dnmb_prophage_no_space_root(), "prophage", "virsorter2")
  list(
    module_dir = module_dir,
    external_root = external_root,
    env_dir = base::file.path(external_root, "conda_env"),
    env_python = base::file.path(external_root, "conda_env", "bin", "python"),
    env_bin = base::file.path(external_root, "conda_env", "bin"),
    cli_path = base::file.path(external_root, "conda_env", "bin", "virsorter"),
    db_dir = base::file.path(external_root, "db"),
    pkgs_dir = base::file.path(external_root, "pkgs")
  )
}

.dnmb_prophage_pide_asset_layout <- function(module_dir) {
  external_root <- base::file.path(.dnmb_prophage_no_space_root(), "prophage", "pide")
  list(
    module_dir = module_dir,
    external_root = external_root,
    repo_dir = base::file.path(external_root, "PIDE"),
    env_dir = base::file.path(external_root, "conda_env"),
    env_python = base::file.path(external_root, "conda_env", "bin", "python"),
    model_dir = base::file.path(external_root, "model"),
    model_archive = base::file.path(external_root, "PIDE.model.tar.gz"),
    pkgs_dir = base::file.path(external_root, "pkgs")
  )
}

.dnmb_prophage_conda_create <- function(env_dir, packages, channels) {
  conda <- .dnmb_prophage_candidate_conda()
  if (!base::nzchar(conda)) {
    return(list(ok = FALSE, detail = "conda not found"))
  }
  external_root <- dirname(env_dir)
  pkgs_dir <- base::file.path(external_root, "pkgs")
  base::dir.create(external_root, recursive = TRUE, showWarnings = FALSE)
  base::dir.create(pkgs_dir, recursive = TRUE, showWarnings = FALSE)
  if (!base::dir.exists(env_dir)) {
    args <- c("create", "-y", "-p", env_dir)
    if (length(channels)) {
      for (ch in channels) {
        args <- c(args, "-c", ch)
      }
    }
    args <- c(args, packages)
    run <- dnmb_run_external(
      conda,
      args = args,
      env = c(CONDA_PKGS_DIRS = pkgs_dir),
      required = FALSE
    )
    if (!base::isTRUE(run$ok)) {
      return(list(ok = FALSE, detail = run$error %||% "conda create failed"))
    }
  }
  snakemake_path <- base::file.path(env_dir, "bin", "snakemake")
  if (!base::file.exists(snakemake_path)) {
    args <- c("install", "-y", "-p", env_dir)
    if (length(channels)) {
      for (ch in channels) {
        args <- c(args, "-c", ch)
      }
    }
    args <- c(args, "snakemake-minimal")
    run <- dnmb_run_external(
      conda,
      args = args,
      env = c(CONDA_PKGS_DIRS = pkgs_dir),
      required = FALSE
    )
    if (!base::isTRUE(run$ok)) {
      return(list(ok = FALSE, detail = run$error %||% "conda install snakemake-minimal failed"))
    }
  }
  list(ok = TRUE, detail = env_dir)
}

.dnmb_prophage_vs2_step_conda_env <- function(layout) {
  create_run <- .dnmb_prophage_conda_create(
    env_dir = layout$env_dir,
    packages = c("python=3.10", "virsorter=2", "snakemake-minimal"),
    channels = c("conda-forge", "bioconda")
  )
  if (!isTRUE(create_run$ok)) {
    return(list(ok = FALSE, detail = create_run$detail))
  }
  if (!base::file.exists(layout$cli_path)) {
    # Partial install: env exists but virsorter binary missing — reinstall
    conda <- .dnmb_prophage_candidate_conda()
    if (base::nzchar(conda)) {
      repair_run <- dnmb_run_external(
        conda,
        args = c("install", "-y", "-p", layout$env_dir, "-c", "conda-forge", "-c", "bioconda", "virsorter=2"),
        env = c(CONDA_PKGS_DIRS = layout$pkgs_dir),
        required = FALSE
      )
      if (base::isTRUE(repair_run$ok) && base::file.exists(layout$cli_path)) {
        return(list(ok = TRUE, detail = layout$env_dir))
      }
    }
    return(list(ok = FALSE, detail = "virsorter executable not found after install"))
  }
  list(ok = TRUE, detail = layout$env_dir)
}

.dnmb_prophage_vs2_step_python_deps <- function(layout) {
  deps_ok <- dnmb_run_external(
    layout$env_python,
    args = c("-c", "import screed, pandas, numpy, sklearn, Bio, joblib"),
    required = FALSE
  )
  if (base::isTRUE(deps_ok$ok)) {
    return(list(ok = TRUE, detail = "dependencies already present"))
  }
  conda <- .dnmb_prophage_candidate_conda()
  dep_run <- dnmb_run_external(
    conda,
    args = c(
      "install", "-y", "-p", layout$env_dir,
      "-c", "conda-forge",
      "-c", "bioconda",
      "screed", "pandas", "numpy<1.24", "scikit-learn",
      "imbalanced-learn", "seaborn", "biopython",
      "prodigal", "hmmer", "last",
      "ncbi-genome-download", "ruamel.yaml"
    ),
    env = c(CONDA_PKGS_DIRS = layout$pkgs_dir),
    required = FALSE
  )
  if (!base::isTRUE(dep_run$ok)) {
    return(list(ok = FALSE, detail = dep_run$error %||% "virsorter dependency install failed"))
  }
  list(ok = TRUE, detail = layout$env_dir)
}

.dnmb_prophage_vs2_run_config <- function(layout, db_dir) {
  dnmb_run_external(
    layout$cli_path,
    args = c("config", "--init-source", "--db-dir", db_dir),
    env = c(
      PATH = paste(unique(c(layout$env_bin, Sys.getenv("PATH"))), collapse = .Platform$path.sep),
      CONDA_PKGS_DIRS = layout$pkgs_dir
    ),
    required = FALSE
  )
}

.dnmb_prophage_vs2_try_archive <- function(layout) {
  archive_path <- base::file.path(layout$db_dir, "db.tgz")
  if (!base::file.exists(archive_path)) {
    return(list(ok = FALSE, detail = "archive not found"))
  }
  .dnmb_prophage_reset_virsorter2_db_extract(layout$db_dir)
  untar_run <- dnmb_run_external("tar", args = c("-xzf", archive_path, "-C", dirname(layout$db_dir)), required = FALSE)
  if (!base::isTRUE(untar_run$ok)) {
    return(list(ok = FALSE, detail = untar_run$error %||% "archive extraction failed"))
  }
  db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
  .dnmb_prophage_vs2_run_config(layout, db_check$db_dir)
  db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
  list(ok = isTRUE(db_check$ok), db_dir = db_check$db_dir, detail = if (isTRUE(db_check$ok)) db_check$db_dir else "archive DB validation failed")
}

.dnmb_prophage_vs2_step_setup_db <- function(layout, cpu = 1L) {
  db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
  if (isTRUE(db_check$ok)) {
    return(list(ok = TRUE, db_dir = db_check$db_dir))
  }
  # Strategy 1: download assets and finalize
  if (!.dnmb_prophage_virsorter2_setup_assets_ready(layout$db_dir)) {
    download_run <- .dnmb_prophage_download_virsorter2_setup_assets(layout$db_dir)
    if (isTRUE(download_run$ok)) {
      finalize_run <- .dnmb_prophage_finalize_virsorter2_db(layout)
      if (isTRUE(finalize_run$ok)) {
        db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
        if (isTRUE(db_check$ok)) return(list(ok = TRUE, db_dir = db_check$db_dir))
      }
    }
  }
  # Strategy 2: finalize from existing assets
  if (.dnmb_prophage_virsorter2_setup_assets_ready(layout$db_dir)) {
    finalize_run <- .dnmb_prophage_finalize_virsorter2_db(layout)
    if (isTRUE(finalize_run$ok)) {
      db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
      if (isTRUE(db_check$ok)) return(list(ok = TRUE, db_dir = db_check$db_dir))
    }
  }
  # Strategy 3: extract from archive
  archive_result <- .dnmb_prophage_vs2_try_archive(layout)
  if (isTRUE(archive_result$ok)) {
    return(list(ok = TRUE, db_dir = archive_result$db_dir))
  }
  # Strategy 4: run virsorter setup command
  setup_run <- dnmb_run_external(
    layout$cli_path,
    args = c("setup", "-d", layout$db_dir, "-j", as.character(base::as.integer(cpu)[1]), "--skip-deps-install"),
    env = c(
      PATH = paste(unique(c(layout$env_bin, Sys.getenv("PATH"))), collapse = .Platform$path.sep),
      CONDA_PKGS_DIRS = layout$pkgs_dir
    ),
    required = FALSE
  )
  db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
  if (isTRUE(db_check$ok)) {
    return(list(ok = TRUE, db_dir = db_check$db_dir))
  }
  # Strategy 5: retry finalize/archive after setup
  if (.dnmb_prophage_virsorter2_setup_assets_ready(layout$db_dir)) {
    finalize_run <- .dnmb_prophage_finalize_virsorter2_db(layout)
    if (isTRUE(finalize_run$ok)) {
      db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
      if (isTRUE(db_check$ok)) return(list(ok = TRUE, db_dir = db_check$db_dir))
    }
  }
  archive_result <- .dnmb_prophage_vs2_try_archive(layout)
  if (isTRUE(archive_result$ok)) {
    return(list(ok = TRUE, db_dir = archive_result$db_dir))
  }
  list(ok = FALSE, db_dir = layout$db_dir, detail = setup_run$error %||% "VirSorter2 DB setup failed after all strategies")
}

.dnmb_prophage_install_virsorter2 <- function(version = "virsorter2",
                                              cache_root = NULL,
                                              install = TRUE,
                                              cpu = 1L,
                                              asset_urls = NULL) {
  module_dir <- .dnmb_db_module_dir(.dnmb_prophage_module_name(), version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_prophage_virsorter2_asset_layout(module_dir)
  if (!base::isTRUE(install)) {
    manifest <- dnmb_db_read_manifest(.dnmb_prophage_module_name(), version, cache_root = cache_root, required = FALSE)
    db_check <- .dnmb_prophage_virsorter2_db_ready(layout$db_dir)
    ready <- !base::is.null(manifest) && base::isTRUE(manifest$install_ok) && base::file.exists(layout$cli_path) && base::isTRUE(db_check$ok)
    if (ready) {
      return(list(ok = TRUE, files = list(cli = layout$cli_path, db_dir = db_check$db_dir, env_python = layout$env_python), manifest = manifest, status = .dnmb_prophage_status_row("prophage_install", "cached", module_dir)))
    }
    return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "missing", "VirSorter2 is missing and module_install is FALSE.")))
  }
  # Step 1: conda environment
  env_result <- .dnmb_prophage_vs2_step_conda_env(layout)
  if (!isTRUE(env_result$ok)) {
    return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "failed", env_result$detail)))
  }
  # Step 2: Python dependencies
  dep_result <- .dnmb_prophage_vs2_step_python_deps(layout)
  if (!isTRUE(dep_result$ok)) {
    return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "failed", dep_result$detail)))
  }
  # Step 3: database setup
  db_result <- .dnmb_prophage_vs2_step_setup_db(layout, cpu = cpu)
  if (!isTRUE(db_result$ok)) {
    return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "failed", db_result$detail %||% "VirSorter2 DB is incomplete after setup.")))
  }
  effective_db_dir <- db_result$db_dir
  manifest <- list(install_ok = TRUE, cli_path = layout$cli_path, db_dir = effective_db_dir, env_python = layout$env_python)
  dnmb_db_write_manifest(.dnmb_prophage_module_name(), version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  list(ok = TRUE, files = list(cli = layout$cli_path, db_dir = effective_db_dir, env_python = layout$env_python), manifest = manifest, status = .dnmb_prophage_status_row("prophage_install", "ok", module_dir))
}

.dnmb_prophage_install_pide <- function(version = "pide",
                                        cache_root = NULL,
                                        install = TRUE,
                                        repo_url = .dnmb_prophage_default_pide_repo_url(),
                                        model_url = .dnmb_prophage_default_pide_model_url()) {
  module_dir <- .dnmb_db_module_dir(.dnmb_prophage_module_name(), version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_prophage_pide_asset_layout(module_dir)
  if (!base::isTRUE(install)) {
    manifest <- dnmb_db_read_manifest(.dnmb_prophage_module_name(), version, cache_root = cache_root, required = FALSE)
    ready <- !base::is.null(manifest) && base::isTRUE(manifest$install_ok) && base::file.exists(layout$env_python) && base::dir.exists(layout$repo_dir) && base::dir.exists(layout$model_dir)
    if (ready) {
      return(list(ok = TRUE, files = list(env_python = layout$env_python, repo_dir = layout$repo_dir, model_dir = layout$model_dir), manifest = manifest, status = .dnmb_prophage_status_row("prophage_install", "cached", module_dir)))
    }
    return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "missing", "PIDE is missing and module_install is FALSE.")))
  }
  create_run <- .dnmb_prophage_conda_create(
    env_dir = layout$env_dir,
    packages = c("python=3.10", "pandas", "biopython", "prodigal", "pip"),
    channels = c("conda-forge", "bioconda")
  )
  if (!isTRUE(create_run$ok)) {
    return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "failed", create_run$detail)))
  }
  for (pkg in c("torch", "fair-esm")) {
    run <- dnmb_run_external(layout$env_python, args = c("-m", "pip", "install", pkg), required = FALSE)
    if (!base::isTRUE(run$ok)) {
      return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "failed", run$error %||% paste("Failed to install", pkg))))
    }
  }
  if (!base::dir.exists(layout$repo_dir)) {
    repo_run <- dnmb_run_external("git", args = c("clone", "--depth", "1", repo_url, layout$repo_dir), required = FALSE)
    if (!base::isTRUE(repo_run$ok)) {
      return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "failed", repo_run$error %||% repo_url)))
    }
  }
  if (!base::dir.exists(layout$model_dir) || !base::length(base::list.files(layout$model_dir, all.files = TRUE, no.. = TRUE))) {
    base::dir.create(layout$model_dir, recursive = TRUE, showWarnings = FALSE)
    dl <- .dnmb_download_asset(model_url, layout$model_archive, insecure = FALSE)
    if (!base::isTRUE(dl$ok)) {
      return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "failed", dl$error %||% model_url)))
    }
    untar_run <- dnmb_run_external("tar", args = c("-xzf", layout$model_archive, "-C", layout$model_dir), required = FALSE)
    if (!base::isTRUE(untar_run$ok)) {
      return(list(ok = FALSE, files = list(), manifest = NULL, status = .dnmb_prophage_status_row("prophage_install", "failed", untar_run$error %||% "Failed to extract PIDE model")))
    }
  }
  manifest <- list(install_ok = TRUE, env_python = layout$env_python, repo_dir = layout$repo_dir, model_dir = layout$model_dir)
  dnmb_db_write_manifest(.dnmb_prophage_module_name(), version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  list(ok = TRUE, files = list(env_python = layout$env_python, repo_dir = layout$repo_dir, model_dir = layout$model_dir), manifest = manifest, status = .dnmb_prophage_status_row("prophage_install", "ok", module_dir))
}

.dnmb_prophage_guess_pide_model <- function(model_dir) {
  candidates <- base::list.files(model_dir, recursive = TRUE, full.names = TRUE)
  candidates <- candidates[grepl("\\.(pt|pth|model|ckpt)$", candidates, ignore.case = TRUE)]
  if (!base::length(candidates)) {
    return("")
  }
  candidates[[1]]
}

.dnmb_prophage_parse_virsorter2_boundary <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

.dnmb_prophage_parse_virsorter2_scores <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

.dnmb_prophage_standardize_virsorter2_regions <- function(boundary_tbl, score_tbl = data.frame()) {
  if (!base::is.data.frame(boundary_tbl) || !base::nrow(boundary_tbl)) {
    return(data.frame())
  }
  b <- base::as.data.frame(boundary_tbl, stringsAsFactors = FALSE)
  names(b) <- gsub("[^a-z0-9]+", "_", tolower(names(b)))
  if (base::is.data.frame(score_tbl) && base::nrow(score_tbl)) {
    s <- base::as.data.frame(score_tbl, stringsAsFactors = FALSE)
    names(s) <- gsub("[^a-z0-9]+", "_", tolower(names(s)))
    key_b <- if ("seqname" %in% names(b)) "seqname" else names(b)[[1]]
    key_s <- if ("seqname" %in% names(s)) "seqname" else names(s)[[1]]
    b <- dplyr::left_join(b, s, by = stats::setNames(key_s, key_b))
  }
  seqname_col <- if ("seqname" %in% names(b)) "seqname" else names(b)[[1]]
  start_col <- if ("trim_bp_start" %in% names(b)) "trim_bp_start" else if ("start" %in% names(b)) "start" else NULL
  end_col <- if ("trim_bp_end" %in% names(b)) "trim_bp_end" else if ("end" %in% names(b)) "end" else NULL
  if (is.null(start_col) || is.null(end_col)) {
    return(data.frame())
  }
  score_col <- if ("trim_pr" %in% names(b)) "trim_pr" else if ("max_score" %in% names(b)) "max_score" else NULL
  group_col <- if ("max_score_group" %in% names(b)) "max_score_group" else if ("group" %in% names(b)) "group" else NULL
  hallmark_col <- if ("hallmark_cnt" %in% names(b)) "hallmark_cnt" else NULL
  partial_col <- if ("partial" %in% names(b)) "partial" else NULL
  seqname <- as.character(b[[seqname_col]])
  data.frame(
    prophage_id = seqname,
    contig = sub("\\|\\|.*$", "", seqname),
    prophage_start = suppressWarnings(as.numeric(b[[start_col]])),
    prophage_end = suppressWarnings(as.numeric(b[[end_col]])),
    prophage_score = if (!is.null(score_col)) suppressWarnings(as.numeric(b[[score_col]])) else NA_real_,
    prophage_group = if (!is.null(group_col)) as.character(b[[group_col]]) else NA_character_,
    hallmark_cnt = if (!is.null(hallmark_col)) suppressWarnings(as.numeric(b[[hallmark_col]])) else NA_real_,
    partial = if (!is.null(partial_col)) as.character(b[[partial_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

.dnmb_prophage_parse_pide_clusters <- function(path) {
  if (!base::file.exists(path)) {
    return(data.frame())
  }
  header <- base::strsplit(base::readLines(path, n = 1L, warn = FALSE), ",", fixed = TRUE)[[1]]
  expected <- c(header, "Total_ORFs", "B_count", "B_ratio")
  utils::read.table(
    path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    fill = TRUE,
    quote = "",
    comment.char = "",
    skip = 1L,
    col.names = expected
  )
}

.dnmb_prophage_standardize_pide_regions <- function(cluster_tbl) {
  if (!base::is.data.frame(cluster_tbl) || !base::nrow(cluster_tbl)) {
    return(data.frame())
  }
  tbl <- base::as.data.frame(cluster_tbl, stringsAsFactors = FALSE)
  names(tbl) <- gsub("[^a-z0-9]+", "_", tolower(names(tbl)))
  contig_col <- if ("contig" %in% names(tbl)) "contig" else names(tbl)[[1]]
  start_col <- if ("start" %in% names(tbl)) "start" else NULL
  end_col <- if ("end" %in% names(tbl)) "end" else NULL
  score_col <- if ("score" %in% names(tbl)) "score" else NULL
  hallmark_col <- if ("total_orfs" %in% names(tbl)) "total_orfs" else NULL
  if (is.null(start_col) || is.null(end_col)) {
    return(data.frame())
  }
  region_id <- paste0(as.character(tbl[[contig_col]]), "_PI", seq_len(nrow(tbl)))
  data.frame(
    prophage_id = region_id,
    contig = as.character(tbl[[contig_col]]),
    prophage_start = suppressWarnings(as.numeric(tbl[[start_col]])),
    prophage_end = suppressWarnings(as.numeric(tbl[[end_col]])),
    prophage_score = if (!is.null(score_col)) suppressWarnings(as.numeric(tbl[[score_col]])) else NA_real_,
    prophage_group = "PIDE_cluster",
    hallmark_cnt = if (!is.null(hallmark_col)) suppressWarnings(as.numeric(tbl[[hallmark_col]])) else NA_real_,
    partial = NA_character_,
    stringsAsFactors = FALSE
  )
}

.dnmb_run_prophage_backend_virsorter2 <- function(genes,
                                                  output_dir,
                                                  cache_root = NULL,
                                                  install = TRUE,
                                                  cpu = 1L,
                                                  genbank = NULL) {
  genbank_records <- .dnmb_prophage_parse_genbank_records(genbank)
  if (!base::nrow(genbank_records)) {
    return(list(ok = FALSE, status = .dnmb_prophage_status_row("prophage_genbank", "failed", "Failed to parse GenBank sequences for VirSorter2."), files = list(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }
  install_result <- .dnmb_prophage_install_virsorter2(cache_root = cache_root, install = install, cpu = cpu)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = install_result$status, files = list(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  input_fna <- base::file.path(output_dir, "virsorter2_input.fna")
  .dnmb_prophage_write_contig_fasta(genbank_records, input_fna)
  stage_dir <- base::tempfile("dnmb-virsorter2-run-")
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  stage_fna <- base::file.path(stage_dir, "input.fna")
  base::file.copy(input_fna, stage_fna, overwrite = TRUE)
  stage_out <- base::file.path(stage_dir, "out")
  run <- dnmb_run_external(
    install_result$files$cli,
    args = c("run", "-w", stage_out, "-d", install_result$files$db_dir, "-i", stage_fna, "--min-length", "1500", "-j", as.character(base::as.integer(cpu)[1]), "--include-groups", "dsDNAphage,ssDNA", "--use-conda-off", "all"),
    env = c(
      PATH = paste(unique(c(dirname(install_result$files$cli), dirname(install_result$files$env_python), Sys.getenv("PATH"))), collapse = .Platform$path.sep),
      PYTHONNOUSERSITE = "1"
    ),
    required = FALSE
  )
  boundary_path <- base::file.path(stage_out, "final-viral-boundary.tsv")
  score_path <- base::file.path(stage_out, "final-viral-score.tsv")
  if (base::file.exists(boundary_path)) base::file.copy(boundary_path, base::file.path(output_dir, "virsorter2_boundary.tsv"), overwrite = TRUE)
  if (base::file.exists(score_path)) base::file.copy(score_path, base::file.path(output_dir, "virsorter2_score.tsv"), overwrite = TRUE)
  region_tbl <- .dnmb_prophage_standardize_virsorter2_regions(
    .dnmb_prophage_parse_virsorter2_boundary(boundary_path),
    .dnmb_prophage_parse_virsorter2_scores(score_path)
  )
  hits <- .dnmb_prophage_hits_from_regions(region_tbl, genes = genes, genbank_records = genbank_records, backend = "virsorter2")
  output_table <- .dnmb_prophage_output_table(genes = genes, hits = hits)
  vs2_ok <- base::isTRUE(run$ok)
  vs2_status <- if (!vs2_ok) "failed" else if (base::nrow(region_tbl) > 0L) "ok" else "empty"
  vs2_detail <- if (!vs2_ok) (run$error %||% "VirSorter2 execution failed") else base::paste0(base::nrow(region_tbl), " prophage regions detected")
  list(
    ok = vs2_ok,
    status = .dnmb_prophage_status_row("prophage_run", vs2_status, vs2_detail),
    files = list(boundary = base::file.path(output_dir, "virsorter2_boundary.tsv"), score = base::file.path(output_dir, "virsorter2_score.tsv")),
    regions = region_tbl,
    hits = hits,
    output_table = output_table
  )
}

.dnmb_run_prophage_backend_pide <- function(genes,
                                            output_dir,
                                            cache_root = NULL,
                                            install = TRUE,
                                            cpu = 1L,
                                            genbank = NULL) {
  genbank_records <- .dnmb_prophage_parse_genbank_records(genbank)
  if (!base::nrow(genbank_records)) {
    return(list(ok = FALSE, status = .dnmb_prophage_status_row("prophage_genbank", "failed", "Failed to parse GenBank sequences for PIDE."), files = list(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }
  install_result <- .dnmb_prophage_install_pide(cache_root = cache_root, install = install)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = install_result$status, files = list(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }
  model_path <- .dnmb_prophage_guess_pide_model(install_result$files$model_dir)
  if (!base::nzchar(model_path) || !base::file.exists(model_path)) {
    return(list(ok = FALSE, status = .dnmb_prophage_status_row("prophage_model", "missing", "PIDE model file not found after installation."), files = list(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  input_fna <- base::file.path(output_dir, "pide_input.fna")
  .dnmb_prophage_write_contig_fasta(genbank_records, input_fna)
  stage_dir <- base::tempfile("dnmb-pide-run-")
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  stage_fna <- base::file.path(stage_dir, "input.fna")
  base::file.copy(input_fna, stage_fna, overwrite = TRUE)
  stage_out <- base::file.path(stage_dir, "out")
  base::dir.create(stage_out, recursive = TRUE, showWarnings = FALSE)
  classify_path <- base::file.path(install_result$files$repo_dir, "classification.py")
  run <- dnmb_run_external(
    install_result$files$env_python,
    args = c(classify_path, stage_fna, model_path, "-o", stage_out),
    env = c(
      PATH = paste(unique(c(dirname(install_result$files$env_python), Sys.getenv("PATH"))), collapse = .Platform$path.sep),
      KMP_DUPLICATE_LIB_OK = "TRUE"
    ),
    required = FALSE
  )
  cluster_path <- base::file.path(stage_out, "cluster.csv")
  virus_contig_path <- base::file.path(stage_out, "virus_contig.txt")
  if (base::file.exists(cluster_path)) {
    base::file.copy(cluster_path, base::file.path(output_dir, "pide_cluster.csv"), overwrite = TRUE)
  }
  if (base::file.exists(virus_contig_path)) {
    base::file.copy(virus_contig_path, base::file.path(output_dir, "pide_virus_contig.txt"), overwrite = TRUE)
  }
  region_tbl <- .dnmb_prophage_standardize_pide_regions(.dnmb_prophage_parse_pide_clusters(cluster_path))
  hits <- .dnmb_prophage_hits_from_regions(region_tbl, genes = genes, genbank_records = genbank_records, backend = "pide")
  output_table <- .dnmb_prophage_output_table(genes = genes, hits = hits)
  pide_ok <- base::isTRUE(run$ok)
  pide_status <- if (!pide_ok) "failed" else if (base::nrow(region_tbl) > 0L) "ok" else "empty"
  pide_detail <- if (!pide_ok) (run$error %||% "PIDE execution failed") else base::paste0(base::nrow(region_tbl), " prophage regions detected")
  list(
    ok = pide_ok,
    status = .dnmb_prophage_status_row("prophage_run", pide_status, pide_detail),
    files = list(cluster = base::file.path(output_dir, "pide_cluster.csv"), virus_contig = base::file.path(output_dir, "pide_virus_contig.txt")),
    regions = region_tbl,
    hits = hits,
    output_table = output_table
  )
}

dnmb_run_prophage_module <- function(genes,
                                     output_dir,
                                     version = .dnmb_prophage_default_version(),
                                     backend = .dnmb_prophage_default_backend(),
                                     cache_root = NULL,
                                     install = TRUE,
                                     repo_url = .dnmb_prophage_default_repo_url(),
                                     cpu = 1L,
                                     genbank = NULL) {
  backend <- .dnmb_prophage_normalize_backend(backend)
  version <- .dnmb_prophage_default_version(backend)
  genbank <- genbank %||% .dnmb_module_detect_genbank(base::getwd())
  if (base::is.null(genbank) || !base::file.exists(genbank)) {
    return(list(ok = FALSE, status = .dnmb_prophage_status_row("prophage_genbank", "missing", "GenBank file is required for prophage analysis."), files = list(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }
  if (identical(backend, "virsorter2")) {
    return(.dnmb_run_prophage_backend_virsorter2(genes = genes, output_dir = output_dir, cache_root = cache_root, install = install, cpu = cpu, genbank = genbank))
  }
  if (identical(backend, "pide")) {
    return(.dnmb_run_prophage_backend_pide(genes = genes, output_dir = output_dir, cache_root = cache_root, install = install, cpu = cpu, genbank = genbank))
  }

  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trace_log <- base::file.path(output_dir, "prophage_module_trace.log")
  status <- .dnmb_prophage_empty_status()
  install_result <- dnmb_prophage_install_module(version = version, cache_root = cache_root, install = install, repo_url = repo_url)
  status <- dplyr::bind_rows(status, install_result$status)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  module <- dnmb_prophage_get_module(version = version, cache_root = cache_root, required = TRUE)
  stage_dir <- base::tempfile("dnmb-phispy-run-")
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  stage_genbank <- base::file.path(stage_dir, base::basename(genbank))
  base::file.copy(genbank, stage_genbank, overwrite = TRUE)
  stage_out <- base::file.path(stage_dir, "out")
  base::dir.create(stage_out, recursive = TRUE, showWarnings = FALSE)

  args <- c(
    stage_genbank,
    "-o", stage_out,
    "--output_choice", "9",
    "--threads", base::as.character(base::as.integer(cpu)[1]),
    "--randomforest_trees", base::as.character(.dnmb_prophage_default_randomforest_trees())
  )
  .dnmb_prophage_trace(trace_log, base::sprintf("[%s] %s", base::Sys.time(), .dnmb_format_command(module$files$cli, args)))
  run <- dnmb_run_external(module$files$cli, args = args, required = FALSE)
  coord_path <- base::file.path(stage_out, "prophage_coordinates.tsv")
  info_path <- base::file.path(stage_out, "prophage_information.tsv")
  if (base::file.exists(coord_path)) base::file.copy(coord_path, base::file.path(output_dir, "prophage_coordinates.tsv"), overwrite = TRUE)
  if (base::file.exists(info_path)) base::file.copy(info_path, base::file.path(output_dir, "prophage_information.tsv"), overwrite = TRUE)

  coord_tbl <- .dnmb_prophage_standardize_coordinates(dnmb_prophage_parse_coordinates(coord_path))
  info_tbl <- .dnmb_prophage_standardize_information(dnmb_prophage_parse_information(info_path))
  gene_hits <- dnmb_prophage_gene_hits(info_tbl, coord_tbl, genes)
  hits <- dnmb_prophage_normalize_hits(gene_hits, backend = "phispy")
  output_table <- .dnmb_prophage_output_table(genes = genes, hits = hits)
  phispy_ok <- base::isTRUE(run$ok)
  phispy_status <- if (!phispy_ok) "failed" else if (base::nrow(coord_tbl) > 0L || base::nrow(info_tbl) > 0L) "ok" else "empty"
  phispy_detail <- if (!phispy_ok) (run$error %||% "PhiSpy execution failed") else base::paste0(base::nrow(coord_tbl), " prophage regions detected")
  status <- dplyr::bind_rows(status, .dnmb_prophage_status_row("prophage_run", phispy_status, phispy_detail))

  list(
    ok = phispy_ok,
    status = status,
    files = list(
      trace_log = trace_log,
      coordinates = base::file.path(output_dir, "prophage_coordinates.tsv"),
      information = base::file.path(output_dir, "prophage_information.tsv")
    ),
    coordinates = coord_tbl,
    information = info_tbl,
    gene_hits = gene_hits,
    hits = hits,
    output_table = output_table
  )
}

.dnmb_prophage_build_summary_table <- function(genbank_table, gbff_path = NULL, output_dir = NULL) {
  # --- Edge case: empty or missing input ---
  if (is.null(genbank_table) || !is.data.frame(genbank_table) || nrow(genbank_table) == 0L) {
    return(data.frame())
  }
  required_cols <- c("Prophage_prophage_id", "locus_tag", "contig", "start", "end", "product")
  missing_cols <- setdiff(required_cols, colnames(genbank_table))
  if (length(missing_cols)) {
    warning("Missing required columns for prophage summary: ",
            paste(missing_cols, collapse = ", "), call. = FALSE)
    return(data.frame())
  }

  # --- 1. Filter to prophage-positive rows ---
  pid <- as.character(genbank_table$Prophage_prophage_id)
  keep <- !is.na(pid) & nzchar(pid) & pid != "0"
  tbl <- genbank_table[keep, , drop = FALSE]
  if (nrow(tbl) == 0L) {
    return(data.frame())
  }

  # --- Helper: safely extract a column or return NA ---
  .safe_col <- function(df, col, default = NA_character_) {
    if (col %in% colnames(df)) df[[col]] else rep(default, nrow(df))
  }

  # --- 2. Per-region summary ---
  ids <- unique(as.character(tbl$Prophage_prophage_id))
  summaries <- lapply(ids, function(rid) {
    sub <- tbl[as.character(tbl$Prophage_prophage_id) == rid, , drop = FALSE]
    gene_count <- nrow(sub)

    # Region boundaries
    prophage_start_vals <- suppressWarnings(as.numeric(.safe_col(sub, "Prophage_prophage_start")))
    prophage_end_vals   <- suppressWarnings(as.numeric(.safe_col(sub, "Prophage_prophage_end")))
    gene_start_vals     <- suppressWarnings(as.numeric(sub$start))
    gene_end_vals       <- suppressWarnings(as.numeric(sub$end))

    region_start <- if (any(!is.na(prophage_start_vals))) min(prophage_start_vals, na.rm = TRUE) else min(gene_start_vals, na.rm = TRUE)
    region_end   <- if (any(!is.na(prophage_end_vals)))   max(prophage_end_vals,   na.rm = TRUE) else max(gene_end_vals,   na.rm = TRUE)
    length_bp    <- region_end - region_start + 1L

    # CDS count (protein_id present = CDS feature)
    cds_count <- if ("protein_id" %in% names(sub)) {
      sum(!is.na(sub$protein_id) & nzchar(as.character(sub$protein_id)), na.rm = TRUE)
    } else {
      gene_count
    }

    # Gene categorization
    category <- .dnmb_gene_arrow_category_prophage(sub$product)
    integrase_present   <- any(category == "Integration", na.rm = TRUE)
    structural_count    <- sum(category %in% c("Head/packaging", "Tail"), na.rm = TRUE)
    lysis_count         <- sum(category == "Lysis", na.rm = TRUE)
    replication_count   <- sum(category == "DNA replication", na.rm = TRUE)
    recombination_count <- sum(category == "Recombination", na.rm = TRUE)
    anti_defense_count  <- sum(category == "Anti-defense", na.rm = TRUE)
    regulation_count    <- sum(category == "Regulation", na.rm = TRUE)
    hyp_denom <- if (cds_count > 0L) cds_count else gene_count
    hypothetical_frac   <- round(sum(category == "Hypothetical", na.rm = TRUE) / hyp_denom, 3)

    # Rank and phage probability
    rank_vals <- suppressWarnings(as.numeric(.safe_col(sub, "Prophage_prophage_rank")))
    mean_rank <- if (all(is.na(rank_vals))) NA_real_ else mean(rank_vals, na.rm = TRUE)

    pp_vals <- suppressWarnings(as.numeric(.safe_col(sub, "Prophage_prophage_pp")))
    phage_gene_frac <- if (all(is.na(pp_vals))) {
      NA_real_
    } else {
      round(sum(pp_vals > 0, na.rm = TRUE) / gene_count, 3)
    }

    # Support field parsing
    support_vals <- .safe_col(sub, "Prophage_support")
    support_all  <- paste(support_vals[!is.na(support_vals)], collapse = " ")
    att_site_detected <- grepl("attL|attR", support_all, ignore.case = TRUE)
    backend_match <- regmatches(support_all, regexpr("backend=([^ ;,]+)", support_all))
    backend <- if (length(backend_match) && nzchar(backend_match)) {
      sub("^backend=", "", backend_match)
    } else {
      NA_character_
    }

    # Completeness proxy
    modules_present <- sum(c(integrase_present, structural_count > 0L, lysis_count > 0L))
    completeness_proxy <- if (modules_present >= 3L) {
      "intact"
    } else if (modules_present == 2L) {
      "questionable"
    } else {
      "incomplete"
    }

    # Confidence tier
    pgf <- if (is.na(phage_gene_frac)) 0 else phage_gene_frac
    confidence_tier <- if (pgf >= 0.5 && structural_count >= 3L) {
      "high"
    } else if (pgf >= 0.3 || structural_count >= 2L) {
      "medium"
    } else {
      "low"
    }

    # Build a human-readable functional summary
    parts <- character()
    if (integrase_present) parts <- c(parts, "integrase")
    if (structural_count > 0L) parts <- c(parts, paste0(structural_count, " structural"))
    if (lysis_count > 0L) parts <- c(parts, paste0(lysis_count, " lysis"))
    if (replication_count > 0L) parts <- c(parts, paste0(replication_count, " replication"))
    if (recombination_count > 0L) parts <- c(parts, paste0(recombination_count, " recombination"))
    if (anti_defense_count > 0L) parts <- c(parts, paste0(anti_defense_count, " anti-defense"))
    if (regulation_count > 0L) parts <- c(parts, paste0(regulation_count, " regulation"))
    functional_summary <- if (length(parts)) paste(parts, collapse = "; ") else "no functional modules detected"

    data.frame(
      region_id           = rid,
      contig              = as.character(sub$contig[1L]),
      start               = region_start,
      end                 = region_end,
      length_bp           = length_bp,
      gene_count          = gene_count,
      cds_count           = cds_count,
      integrase_present   = integrase_present,
      structural_count    = structural_count,
      lysis_count         = lysis_count,
      replication_count   = replication_count,
      recombination_count = recombination_count,
      anti_defense_count  = anti_defense_count,
      regulation_count    = regulation_count,
      hypothetical_frac   = hypothetical_frac,
      mean_rank           = mean_rank,
      phage_gene_frac     = phage_gene_frac,
      att_site_detected   = att_site_detected,
      completeness_proxy  = completeness_proxy,
      backend             = backend,
      confidence_tier     = confidence_tier,
      functional_summary  = functional_summary,
      stringsAsFactors    = FALSE
    )
  })
  summary_df <- dplyr::bind_rows(summaries)

  # --- 3. GC content computation ---
  if (!is.null(gbff_path) && file.exists(gbff_path)) {
    records <- .dnmb_prophage_parse_genbank_records(gbff_path)
    if (is.data.frame(records) && nrow(records) > 0L && "sequence" %in% colnames(records)) {
      # Helper to compute GC fraction from a sequence string
      .gc_frac <- function(seq_str) {
        if (is.na(seq_str) || !nzchar(seq_str)) return(NA_real_)
        chars <- strsplit(seq_str, "")[[1]]
        gc <- sum(chars %in% c("G", "C"))
        round(gc / length(chars), 4)
      }
      # Host GC: weighted average across all contigs
      all_seqs <- records$sequence[!is.na(records$sequence)]
      host_gc <- if (length(all_seqs)) {
        full <- paste(all_seqs, collapse = "")
        .gc_frac(full)
      } else {
        NA_real_
      }

      # Match contigs using accession or locus
      region_gc <- vapply(seq_len(nrow(summary_df)), function(i) {
        ctg <- summary_df$contig[i]
        ctg_norm <- .dnmb_normalize_contig_label(ctg)
        rec_idx <- which(
          .dnmb_normalize_contig_label(records$accession) == ctg_norm |
          .dnmb_normalize_contig_label(records$locus) == ctg_norm |
          .dnmb_normalize_contig_label(records$definition) == ctg_norm
        )
        if (!length(rec_idx)) return(NA_real_)
        seq_full <- records$sequence[rec_idx[1L]]
        if (is.na(seq_full) || !nzchar(seq_full)) return(NA_real_)
        s <- max(1L, summary_df$start[i])
        e <- min(nchar(seq_full), summary_df$end[i])
        if (s > e) return(NA_real_)
        subseq <- substr(seq_full, s, e)
        .gc_frac(subseq)
      }, numeric(1))

      summary_df$region_gc <- region_gc
      summary_df$host_gc   <- host_gc
      summary_df$delta_gc  <- round(region_gc - host_gc, 4)
    }
  }

  # Add GC columns as NA if not computed
  if (!"region_gc" %in% colnames(summary_df)) {
    summary_df$region_gc <- NA_real_
    summary_df$host_gc   <- NA_real_
    summary_df$delta_gc  <- NA_real_
  }

  # --- 4. Integration site analysis ---
  full_tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  full_tbl$gene_start_num <- suppressWarnings(base::as.numeric(full_tbl$start))
  full_tbl$gene_end_num   <- suppressWarnings(base::as.numeric(full_tbl$end))
  has_anticodon <- "anticodon" %in% names(full_tbl)

  integration_cols <- vapply(seq_len(nrow(summary_df)), function(i) {
    ctg <- summary_df$contig[i]
    pp_start <- summary_df$start[i]
    pp_end   <- summary_df$end[i]
    contig_genes <- full_tbl[full_tbl$contig == ctg, , drop = FALSE]
    contig_genes <- contig_genes[order(contig_genes$gene_start_num), , drop = FALSE]
    flank_bp <- 5000L

    # Integrase info
    pp_genes <- contig_genes[
      contig_genes$gene_start_num >= pp_start & contig_genes$gene_end_num <= pp_end, , drop = FALSE]
    int_genes <- pp_genes[grepl("integrase|excisionase", tolower(pp_genes$product)), , drop = FALSE]
    integrase_family <- if (nrow(int_genes)) {
      paste(unique(trimws(int_genes$product)), collapse = "; ")
    } else {
      NA_character_
    }
    integrase_position <- if (nrow(int_genes)) {
      int_mid <- (as.numeric(int_genes$gene_start_num[1]) + as.numeric(int_genes$gene_end_num[1])) / 2
      region_mid <- (pp_start + pp_end) / 2
      if (int_mid < pp_start + (pp_end - pp_start) * 0.25) "left_end"
      else if (int_mid > pp_start + (pp_end - pp_start) * 0.75) "right_end"
      else "internal"
    } else {
      NA_character_
    }

    # tRNA near boundaries
    trna_genes <- if (has_anticodon) {
      contig_genes[!is.na(contig_genes$anticodon) & nzchar(contig_genes$anticodon), , drop = FALSE]
    } else {
      contig_genes[grepl("tRNA", contig_genes$product, ignore.case = TRUE), , drop = FALSE]
    }
    trna_near_left <- trna_genes[
      abs(trna_genes$gene_start_num - pp_start) < flank_bp | abs(trna_genes$gene_end_num - pp_start) < flank_bp,
      , drop = FALSE]
    trna_near_right <- trna_genes[
      abs(trna_genes$gene_start_num - pp_end) < flank_bp | abs(trna_genes$gene_end_num - pp_end) < flank_bp,
      , drop = FALSE]
    trna_near_all <- unique(rbind(trna_near_left, trna_near_right))

    trna_target <- if (nrow(trna_near_all)) {
      labels <- vapply(seq_len(nrow(trna_near_all)), function(j) {
        ac <- if (has_anticodon && !is.na(trna_near_all$anticodon[j])) trna_near_all$anticodon[j] else ""
        prod <- trna_near_all$product[j]
        dist_left <- min(abs(trna_near_all$gene_start_num[j] - pp_start), abs(trna_near_all$gene_end_num[j] - pp_start))
        dist_right <- min(abs(trna_near_all$gene_start_num[j] - pp_end), abs(trna_near_all$gene_end_num[j] - pp_end))
        boundary <- if (dist_left < dist_right) "L" else "R"
        dist <- min(dist_left, dist_right)
        inside <- trna_near_all$gene_start_num[j] >= pp_start & trna_near_all$gene_end_num[j] <= pp_end
        paste0(prod, if (nzchar(ac)) paste0("(", ac, ")") else "",
               " ", boundary, ":", dist, "bp",
               if (inside) " DISRUPTED" else "")
      }, character(1))
      paste(labels, collapse = "; ")
    } else {
      NA_character_
    }
    trna_disrupted <- if (nrow(trna_near_all)) {
      any(trna_near_all$gene_start_num >= pp_start & trna_near_all$gene_end_num <= pp_end)
    } else {
      FALSE
    }

    # Flanking host genes (closest 3 each side)
    pp_idx <- which(contig_genes$gene_start_num >= pp_start & contig_genes$gene_end_num <= pp_end)
    left_flank <- if (length(pp_idx) && min(pp_idx) > 1) {
      idx <- max(1, min(pp_idx) - 3):(min(pp_idx) - 1)
      paste(contig_genes$product[idx], collapse = "; ")
    } else {
      NA_character_
    }
    right_flank <- if (length(pp_idx) && max(pp_idx) < nrow(contig_genes)) {
      idx <- (max(pp_idx) + 1):min(nrow(contig_genes), max(pp_idx) + 3)
      paste(contig_genes$product[idx], collapse = "; ")
    } else {
      NA_character_
    }

    c(integrase_family = integrase_family,
      integrase_position = integrase_position,
      trna_target = trna_target,
      trna_disrupted = as.character(trna_disrupted),
      left_flank_genes = left_flank,
      right_flank_genes = right_flank)
  }, character(6))

  summary_df$integrase_family    <- integration_cols["integrase_family", ]
  summary_df$integrase_position  <- integration_cols["integrase_position", ]
  summary_df$trna_target         <- integration_cols["trna_target", ]
  summary_df$trna_disrupted      <- as.logical(integration_cols["trna_disrupted", ])
  summary_df$left_flank_genes    <- integration_cols["left_flank_genes", ]
  summary_df$right_flank_genes   <- integration_cols["right_flank_genes", ]

  # --- 5. att site detection (multi-layer) ---
  att_results <- lapply(seq_len(nrow(summary_df)), function(i) {
    tryCatch(
      .dnmb_prophage_detect_att_sites(
        genome_seq = if (exists("records") && is.data.frame(records) && nrow(records)) {
          ctg_norm <- .dnmb_normalize_contig_label(summary_df$contig[i])
          rec_idx <- which(
            .dnmb_normalize_contig_label(records$accession) == ctg_norm |
            .dnmb_normalize_contig_label(records$locus) == ctg_norm |
            .dnmb_normalize_contig_label(records$definition) == ctg_norm
          )
          if (length(rec_idx)) records$sequence[[rec_idx[1]]] else ""
        } else "",
        pp_start = summary_df$start[i],
        pp_end = summary_df$end[i],
        genbank_table = genbank_table,
        contig = summary_df$contig[i],
        cache_dir = tryCatch({
          plot_dir <- file.path(output_dir %||% getwd(), "visualizations", "_cache",
                                paste0("prophage_ref_", summary_df$region_id[i]))
          if (dir.exists(plot_dir)) plot_dir else NULL
        }, error = function(e) NULL),
        comparative = FALSE
      ),
      error = function(e) .dnmb_prophage_empty_att_result()
    )
  })
  summary_df$att_core_seq      <- vapply(att_results, function(r) r$att_core_sequence %||% NA_character_, character(1))
  summary_df$att_core_length   <- vapply(att_results, function(r) as.integer(r$att_core_length %||% NA_integer_), integer(1))
  summary_df$att_evidence      <- vapply(att_results, function(r) paste(r$evidence_layers, collapse=";"), character(1))
  summary_df$att_confidence    <- vapply(att_results, function(r) r$confidence %||% NA_character_, character(1))
  summary_df$attL_seq          <- vapply(att_results, function(r) r$attL_sequence %||% NA_character_, character(1))
  summary_df$attR_seq          <- vapply(att_results, function(r) r$attR_sequence %||% NA_character_, character(1))
  summary_df$attB_seq          <- vapply(att_results, function(r) r$attB_sequence %||% NA_character_, character(1))
  summary_df$att_related_genome <- vapply(att_results, function(r) r$related_genome_used %||% NA_character_, character(1))

  rownames(summary_df) <- NULL
  summary_df
}

.dnmb_prophage_empty_att_result <- function() {
  list(
    att_core_sequence = NA_character_,
    att_core_length = NA_integer_,
    attL_sequence = NA_character_, attL_start = NA_integer_, attL_end = NA_integer_,
    attR_sequence = NA_character_, attR_start = NA_integer_, attR_end = NA_integer_,
    attB_sequence = NA_character_,
    evidence_layers = character(),
    confidence = "none",
    related_genome_used = NA_character_
  )
}

.dnmb_prophage_clean_dna <- function(seq) {
  seq <- base::toupper(base::gsub("[^ACGT]", "", base::as.character(seq)[1]))
  if (is.na(seq)) {
    return("")
  }
  seq
}

.dnmb_prophage_revcomp <- function(seq) {
  seq <- .dnmb_prophage_clean_dna(seq)
  if (!nzchar(seq)) {
    return("")
  }
  base::paste(rev(strsplit(chartr("ACGT", "TGCA", seq), "", fixed = TRUE)[[1]]), collapse = "")
}

.dnmb_prophage_find_inverted_repeats <- function(left_seq, right_seq, min_len = 8L, max_len = 50L) {
  # Inverted repeats: left boundary matches reverse complement of right boundary
  left_seq <- .dnmb_prophage_clean_dna(left_seq)
  right_rc <- .dnmb_prophage_revcomp(right_seq)
  if (!nzchar(left_seq) || !nzchar(right_rc)) {
    return(list(seq = "", len = 0L, left_pos = NA_integer_, right_pos = NA_integer_))
  }
  best <- list(seq = "", len = 0L, left_pos = NA_integer_, right_pos = NA_integer_)
  for (k in seq.int(min(max_len, nchar(left_seq), nchar(right_rc)), min_len, by = -1L)) {
    left_kmers <- vapply(seq_len(nchar(left_seq) - k + 1L), function(i) substr(left_seq, i, i + k - 1L), character(1))
    right_kmers <- vapply(seq_len(nchar(right_rc) - k + 1L), function(i) substr(right_rc, i, i + k - 1L), character(1))
    shared <- which(left_kmers %in% right_kmers)
    if (length(shared)) {
      best <- list(
        seq = left_kmers[shared[1]],
        len = k,
        left_pos = shared[1],
        right_pos = match(left_kmers[shared[1]], right_kmers)
      )
      break
    }
  }
  best
}

.dnmb_prophage_find_tandem_repeats <- function(seq, min_unit = 2L, max_unit = 8L, min_copies = 3L) {
  # Find short tandem repeats (STRs) in a sequence
  seq <- .dnmb_prophage_clean_dna(seq)
  if (!nzchar(seq)) {
    return(data.frame(unit = character(), unit_len = integer(), copies = integer(),
                      total_len = integer(), position = integer(), stringsAsFactors = FALSE))
  }
  results <- list()
  for (unit_len in min_unit:max_unit) {
    for (pos in seq_len(nchar(seq) - unit_len * min_copies + 1L)) {
      unit <- substr(seq, pos, pos + unit_len - 1L)
      copies <- 1L
      check_pos <- pos + unit_len
      while (check_pos + unit_len - 1L <= nchar(seq) &&
             substr(seq, check_pos, check_pos + unit_len - 1L) == unit) {
        copies <- copies + 1L
        check_pos <- check_pos + unit_len
      }
      if (copies >= min_copies) {
        results[[length(results) + 1L]] <- data.frame(
          unit = unit, unit_len = unit_len, copies = copies,
          total_len = unit_len * copies, position = pos, stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(results)) {
    return(data.frame(unit = character(), unit_len = integer(), copies = integer(),
                      total_len = integer(), position = integer(), stringsAsFactors = FALSE))
  }
  out <- dplyr::bind_rows(results)
  out <- out[!duplicated(paste(out$unit, out$position)), , drop = FALSE]
  out[order(-out$total_len, out$position), , drop = FALSE]
}

.dnmb_prophage_auto_find_reference <- function(genome_seq, pp_start, pp_end,
                                               flank_size = 1000L,
                                               candidate_fastas = NULL) {
  # Concatenate flanks to create "simulated empty site"
  # BLAST this against candidate genomes or nt — a hit spanning the junction = genome without prophage
  left_flank <- substr(genome_seq, max(1, pp_start - flank_size), pp_start - 1)
  right_flank <- substr(genome_seq, pp_end + 1, min(nchar(genome_seq), pp_end + flank_size))
  simulated_attB <- paste0(left_flank, right_flank)
  junction_pos <- nchar(left_flank)

  tmp_dir <- base::tempdir()
  query_path <- file.path(tmp_dir, "simulated_attB.fna")
  writeLines(c(">simulated_attB", simulated_attB), query_path)

  # --- Helper: parse BLAST outfmt 6 result and find best reference ---
  .parse_blast_hits <- function(result_path, junction_pos, accession_label = NULL) {
    if (!file.exists(result_path) || file.info(result_path)$size == 0) {
      return(NULL)
    }
    hits <- tryCatch(
      utils::read.table(result_path, sep = "\t", header = FALSE, quote = "",
                        comment.char = "", stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(hits) || !nrow(hits)) return(NULL)

    # outfmt "6" default columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    if (ncol(hits) >= 12L) {
      names(hits)[1:12] <- c("qseqid", "saccver", "pident", "length",
                             "mismatch", "gapopen", "qstart", "qend",
                             "sstart", "send", "evalue", "bitscore")
    } else {
      return(NULL)
    }

    # Find hits that span the junction (= both flanks align to the same genome)
    spanning <- hits[hits$qstart < junction_pos & hits$qend > junction_pos, , drop = FALSE]
    if (nrow(spanning)) {
      best <- spanning[which.max(spanning$length), ]
      return(list(
        accession = if (!is.null(accession_label)) accession_label else best$saccver,
        title = best$saccver,
        identity = best$pident,
        length = best$length,
        ref_start = min(best$sstart, best$send),
        ref_end = max(best$sstart, best$send),
        spanning = TRUE
      ))
    }

    # Try finding pairs: left hit + right hit from same subject
    left_hits <- hits[hits$qend <= junction_pos + 100, , drop = FALSE]
    right_hits <- hits[hits$qstart >= junction_pos - 100, , drop = FALSE]
    common_acc <- intersect(left_hits$saccver, right_hits$saccver)
    if (length(common_acc)) {
      best_acc <- common_acc[1]
      lh <- left_hits[left_hits$saccver == best_acc, , drop = FALSE][1, ]
      rh <- right_hits[right_hits$saccver == best_acc, , drop = FALSE][1, ]
      return(list(
        accession = if (!is.null(accession_label)) accession_label else best_acc,
        title = best_acc,
        left_identity = lh$pident,
        right_identity = rh$pident,
        left_ref_pos = c(lh$sstart, lh$send),
        right_ref_pos = c(rh$sstart, rh$send),
        spanning = FALSE
      ))
    }
    NULL
  }

  # --- 1. Try local BLAST against each candidate FASTA ---
  if (!is.null(candidate_fastas) && length(candidate_fastas)) {
    candidate_fastas <- candidate_fastas[file.exists(candidate_fastas)]
    best_local <- NULL
    best_local_identity <- 0
    best_local_gap <- NA_real_

    for (cand_path in candidate_fastas) {
      result_path <- file.path(tmp_dir, paste0("attB_blast_", basename(cand_path), ".tsv"))
      # Copy to tempdir to avoid path-with-spaces issues in BLAST
      safe_cand <- file.path(tmp_dir, paste0("cand_", basename(cand_path)))
      if (!file.exists(safe_cand)) base::file.copy(cand_path, safe_cand, overwrite = TRUE)
      db_prefix <- file.path(tmp_dir, paste0("cand_db_", basename(cand_path)))
      mk_run <- dnmb_run_external("makeblastdb", args = c(
        "-in", safe_cand, "-dbtype", "nucl", "-out", db_prefix
      ), required = FALSE)
      if (!isTRUE(mk_run$ok)) next

      run <- dnmb_run_external("blastn", args = c(
        "-task", "megablast",
        "-db", db_prefix,
        "-query", query_path,
        "-max_target_seqs", "10",
        "-evalue", "1e-10",
        "-outfmt", "6",
        "-out", result_path
      ), required = FALSE)

      if (!isTRUE(run$ok)) next

      # Extract real accession from FASTA header or filename
      cand_first_line <- readLines(cand_path, n = 1L, warn = FALSE)
      cand_label <- if (grepl("^>", cand_first_line)) {
        trimws(sub("^>", "", sub("\\s.*$", "", cand_first_line)))
      } else {
        sub("^ref_", "", sub("\\.[^.]*$", "", basename(cand_path)))
      }
      hit <- .parse_blast_hits(result_path, junction_pos, accession_label = cand_label)
      if (!is.null(hit)) {
        # Check gap size: small gap = empty site (preferred), large gap = another prophage
        hit_gap <- NA_real_
        if (isTRUE(hit$spanning)) {
          hit_gap <- 0  # spanning = empty site
        } else if (!is.null(hit$left_ref_pos) && !is.null(hit$right_ref_pos)) {
          hit_gap <- abs(min(hit$right_ref_pos) - max(hit$left_ref_pos))
        }
        # Prefer empty sites (gap < 5kb) over occupied sites (gap > 5kb)
        hit_id <- if (isTRUE(hit$spanning)) hit$identity else mean(c(hit$left_identity, hit$right_identity))
        is_empty <- !is.na(hit_gap) && hit_gap < 5000
        prev_empty <- !is.na(best_local_gap) && best_local_gap < 5000
        if (is.null(best_local) || (is_empty && !prev_empty) || (is_empty == prev_empty && hit_id > best_local_identity)) {
          best_local <- hit
          best_local_identity <- hit_id
          best_local_gap <- hit_gap
        }
      }
      # Clean up temporary database files
      unlink(list.files(tmp_dir, pattern = paste0("^cand_db_", basename(cand_path)), full.names = TRUE))
    }

    # Return if empty site found; otherwise try right-only for a better reference
    if (!is.null(best_local) && !is.na(best_local_gap) && best_local_gap < 5000) return(best_local)

    # Try right-flank-only match (tRNA conserved across species, left may diverge)
    right_query <- file.path(tmp_dir, "right_flank_only.fna")
    writeLines(c(">right_flank", right_flank), right_query)
    for (cand_path in candidate_fastas) {
      safe_cand <- file.path(tmp_dir, paste0("cand_", basename(cand_path)))
      if (!file.exists(safe_cand)) base::file.copy(cand_path, safe_cand, overwrite = TRUE)
      db_prefix <- file.path(tmp_dir, paste0("cand_rdb_", basename(cand_path)))
      mk_run <- dnmb_run_external("makeblastdb", args = c("-in", safe_cand, "-dbtype", "nucl", "-out", db_prefix), required = FALSE)
      if (!isTRUE(mk_run$ok)) next
      result_path <- file.path(tmp_dir, paste0("attB_right_", basename(cand_path), ".tsv"))
      run <- dnmb_run_external("blastn", args = c("-task", "megablast", "-db", db_prefix, "-query", right_query,
                               "-evalue", "1e-10", "-outfmt", "6", "-out", result_path), required = FALSE)
      if (!isTRUE(run$ok) || !file.exists(result_path) || file.info(result_path)$size == 0) next
      hits <- tryCatch(utils::read.table(result_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(hits) || !nrow(hits) || ncol(hits) < 10) next
      best_hit <- hits[which.max(as.integer(hits$V4)), ]
      if (as.numeric(best_hit$V3) >= 85) {
        cand_first_line <- readLines(cand_path, n = 1L, warn = FALSE)
        cand_label <- if (grepl("^>", cand_first_line)) trimws(sub("^>", "", sub("\\s.*$", "", cand_first_line)))
                      else sub("^ref_", "", sub("\\.[^.]*$", "", basename(cand_path)))
        return(list(
          accession = cand_label,
          title = cand_label,
          right_identity = as.numeric(best_hit$V3),
          right_ref_pos = c(as.integer(best_hit$V9), as.integer(best_hit$V10)),
          spanning = FALSE,
          right_only = TRUE
        ))
      }
    }
  }

  NULL
}

.dnmb_prophage_boundary_windows <- function(genome_seq, pp_start, pp_end, flank_bp = 500L) {
  genome_seq <- .dnmb_prophage_clean_dna(genome_seq)
  seq_len <- nchar(genome_seq)
  left_start <- max(1L, as.integer(pp_start) - as.integer(flank_bp))
  left_end <- min(seq_len, as.integer(pp_start) + as.integer(flank_bp))
  right_start <- max(1L, as.integer(pp_end) - as.integer(flank_bp))
  right_end <- min(seq_len, as.integer(pp_end) + as.integer(flank_bp))
  list(
    left = list(
      seq = substr(genome_seq, left_start, left_end),
      window_start = left_start,
      window_end = left_end,
      boundary = as.integer(pp_start)
    ),
    right = list(
      seq = substr(genome_seq, right_start, right_end),
      window_start = right_start,
      window_end = right_end,
      boundary = as.integer(pp_end)
    )
  )
}

.dnmb_prophage_locate_pattern <- function(pattern, window) {
  pattern <- .dnmb_prophage_clean_dna(pattern)
  if (!nzchar(pattern) || !length(window$seq) || !nzchar(window$seq)) {
    return(NULL)
  }
  hits <- gregexpr(pattern, window$seq, fixed = TRUE)[[1]]
  if (!length(hits) || identical(hits[[1]], -1L)) {
    return(NULL)
  }
  centers <- window$window_start + hits - 1L + floor((nchar(pattern) - 1L) / 2L)
  best_idx <- which.min(abs(centers - window$boundary))
  start <- window$window_start + hits[[best_idx]] - 1L
  list(start = start, end = start + nchar(pattern) - 1L)
}

.dnmb_prophage_find_direct_repeats <- function(left_seq,
                                               right_seq,
                                               min_len = 10L,
                                               max_len = 50L,
                                               max_boundary_distance = 100L) {
  left_seq <- .dnmb_prophage_clean_dna(left_seq)
  right_seq <- .dnmb_prophage_clean_dna(right_seq)
  best <- list(
    seq = "",
    len = 0L,
    identity = 0,
    score = 0,
    left_pos = NA_integer_,
    right_pos = NA_integer_
  )
  if (!nzchar(left_seq) || !nzchar(right_seq)) {
    return(best)
  }
  max_len <- min(as.integer(max_len), nchar(left_seq), nchar(right_seq))
  min_len <- min(as.integer(min_len), max_len)
  if (max_len < min_len) {
    return(best)
  }
  for (k in seq.int(max_len, min_len, by = -1L)) {
    left_idx <- seq_len(nchar(left_seq) - k + 1L)
    left_kmers <- vapply(left_idx, function(i) substr(left_seq, i, i + k - 1L), character(1))
    right_idx <- seq_len(nchar(right_seq) - k + 1L)
    right_kmers <- vapply(right_idx, function(i) substr(right_seq, i, i + k - 1L), character(1))
    shared <- intersect(unique(left_kmers), unique(right_kmers))
    if (!length(shared)) {
      next
    }
    left_center <- (nchar(left_seq) + 1) / 2
    right_center <- (nchar(right_seq) + 1) / 2
    ranked <- lapply(shared, function(kmer) {
      left_matches <- which(left_kmers == kmer)
      right_matches <- which(right_kmers == kmer)
      lp <- left_matches[which.min(abs((left_matches + (k - 1) / 2) - left_center))]
      rp <- right_matches[which.min(abs((right_matches + (k - 1) / 2) - right_center))]
      left_dist <- abs((lp + (k - 1) / 2) - left_center)
      right_dist <- abs((rp + (k - 1) / 2) - right_center)
      data.frame(
        kmer = kmer,
        left_pos = lp,
        right_pos = rp,
        left_dist = left_dist,
        right_dist = right_dist,
        dist = left_dist + right_dist,
        stringsAsFactors = FALSE
      )
    })
    ranked <- dplyr::bind_rows(ranked)
    ranked <- ranked[
      ranked$left_dist <= as.integer(max_boundary_distance) &
        ranked$right_dist <= as.integer(max_boundary_distance),
      ,
      drop = FALSE
    ]
    if (!nrow(ranked)) {
      next
    }
    ranked <- ranked[order(ranked$dist, ranked$left_pos, ranked$right_pos), , drop = FALSE]
    best <- list(
      seq = ranked$kmer[[1]],
      len = as.integer(k),
      identity = 100,
      score = as.integer(k) * 100,
      left_pos = as.integer(ranked$left_pos[[1]]),
      right_pos = as.integer(ranked$right_pos[[1]])
    )
    break
  }
  best
}

.dnmb_prophage_best_core_from_source <- function(source_seq,
                                                 windows,
                                                 min_len = 10L,
                                                 max_len = 50L,
                                                 max_boundary_distance = 100L) {
  source_seq <- .dnmb_prophage_clean_dna(source_seq)
  best <- list(seq = "", len = 0L, score = 0L)
  if (!nzchar(source_seq)) {
    return(best)
  }
  max_len <- min(as.integer(max_len), nchar(source_seq))
  min_len <- min(as.integer(min_len), max_len)
  if (max_len < min_len) {
    return(best)
  }
  source_center <- (nchar(source_seq) + 1) / 2
  for (k in seq.int(max_len, min_len, by = -1L)) {
    idx <- seq_len(nchar(source_seq) - k + 1L)
    idx <- idx[order(abs((idx + (k - 1) / 2) - source_center))]
    seen <- character()
    for (i in idx) {
      candidate <- substr(source_seq, i, i + k - 1L)
      if (candidate %in% seen) {
        next
      }
      seen <- c(seen, candidate)
      left_hit <- .dnmb_prophage_locate_pattern(candidate, windows$left)
      right_hit <- .dnmb_prophage_locate_pattern(candidate, windows$right)
      if (is.null(left_hit) || is.null(right_hit)) {
        next
      }
      left_mid <- (left_hit$start + left_hit$end) / 2
      right_mid <- (right_hit$start + right_hit$end) / 2
      if (abs(left_mid - windows$left$boundary) > as.integer(max_boundary_distance) ||
          abs(right_mid - windows$right$boundary) > as.integer(max_boundary_distance)) {
        next
      }
      return(list(seq = candidate, len = as.integer(k), score = as.integer(k) * 100L))
    }
  }
  best
}

.dnmb_prophage_make_candidate <- function(layer,
                                          seq,
                                          windows,
                                          score = NULL,
                                          attB_sequence = NA_character_,
                                          related_genome = NA_character_,
                                          source_note = NA_character_) {
  seq <- .dnmb_prophage_clean_dna(seq)
  if (!nzchar(seq)) {
    return(NULL)
  }
  left_hit <- .dnmb_prophage_locate_pattern(seq, windows$left)
  right_hit <- .dnmb_prophage_locate_pattern(seq, windows$right)
  if (is.null(left_hit) || is.null(right_hit)) {
    return(NULL)
  }
  if (is.null(score) || !length(score) || is.na(score[[1]])) {
    score <- nchar(seq) * 100
  }
  list(
    layer = base::as.character(layer)[1],
    seq = seq,
    len = as.integer(nchar(seq)),
    score = as.numeric(score)[1],
    attL_start = as.integer(left_hit$start),
    attL_end = as.integer(left_hit$end),
    attR_start = as.integer(right_hit$start),
    attR_end = as.integer(right_hit$end),
    attB_sequence = if (!is.na(attB_sequence[[1]]) && nzchar(attB_sequence[[1]])) {
      .dnmb_prophage_clean_dna(attB_sequence[[1]])
    } else {
      NA_character_
    },
    related_genome = base::as.character(related_genome)[1],
    source_note = base::as.character(source_note)[1]
  )
}

.dnmb_prophage_extract_feature_sequence <- function(feature_row, genome_seq) {
  genome_seq <- .dnmb_prophage_clean_dna(genome_seq)
  if (!nzchar(genome_seq)) {
    return("")
  }
  start <- suppressWarnings(as.integer(min(feature_row$gs[[1]], feature_row$ge[[1]], na.rm = TRUE)))
  end <- suppressWarnings(as.integer(max(feature_row$gs[[1]], feature_row$ge[[1]], na.rm = TRUE)))
  if (!is.finite(start) || !is.finite(end) || start < 1L || end > nchar(genome_seq) || start > end) {
    return("")
  }
  seq <- substr(genome_seq, start, end)
  strand <- if ("strand" %in% names(feature_row)) feature_row$strand[[1]] else if ("direction" %in% names(feature_row)) feature_row$direction[[1]] else "+"
  if (!is.na(strand) && identical(as.character(strand), "-")) {
    seq <- .dnmb_prophage_revcomp(seq)
  }
  seq
}

.dnmb_prophage_is_trna_feature <- function(tbl) {
  feature_types <- if ("feature_types" %in% names(tbl)) tbl$feature_types else ""
  products <- if ("product" %in% names(tbl)) tbl$product else ""
  anticodon <- if ("anticodon" %in% names(tbl)) tbl$anticodon else ""
  (!is.na(anticodon) & nzchar(as.character(anticodon))) |
    grepl("tRNA", feature_types, ignore.case = TRUE) |
    grepl("tRNA", products, ignore.case = TRUE)
}

.dnmb_prophage_parse_location_simple <- function(raw_location) {
  nums <- suppressWarnings(as.integer(unlist(regmatches(raw_location, gregexpr("[0-9]+", raw_location)))))
  nums <- nums[!is.na(nums)]
  if (!length(nums)) {
    return(list(start = NA_integer_, end = NA_integer_, strand = "+"))
  }
  list(
    start = min(nums),
    end = max(nums),
    strand = if (grepl("complement", raw_location, fixed = TRUE)) "-" else "+"
  )
}

.dnmb_prophage_reference_att_annotations <- function(gbk_file,
                                                     ref_seq,
                                                     terminal_bp = 500L) {
  if (!file.exists(gbk_file) || !nzchar(ref_seq)) {
    return(character())
  }
  lines <- readLines(gbk_file, warn = FALSE)
  feat_idx <- grep("^\\s{5}\\S", lines)
  if (!length(feat_idx)) {
    return(character())
  }
  candidates <- character()
  for (i in seq_along(feat_idx)) {
    start_idx <- feat_idx[[i]]
    end_idx <- if (i < length(feat_idx)) feat_idx[[i + 1L]] - 1L else length(lines)
    block <- lines[start_idx:end_idx]
    feature_type <- trimws(substr(block[[1]], 6, 20))
    qualifier_blob <- paste(block, collapse = " ")
    is_att_feature <- grepl("attP|attachment", qualifier_blob, ignore.case = TRUE)
    is_repeat_feature <- identical(feature_type, "repeat_region")
    if (!is_att_feature && !is_repeat_feature) {
      next
    }
    loc_lines <- block[!grepl("^\\s{21}/", block)]
    loc_text <- paste(trimws(loc_lines), collapse = "")
    loc <- .dnmb_prophage_parse_location_simple(loc_text)
    if (!is.finite(loc$start) || !is.finite(loc$end)) {
      next
    }
    if (!(loc$start <= terminal_bp || loc$end >= (nchar(ref_seq) - terminal_bp))) {
      next
    }
    seq <- substr(ref_seq, loc$start, loc$end)
    if (identical(loc$strand, "-")) {
      seq <- .dnmb_prophage_revcomp(seq)
    }
    seq <- .dnmb_prophage_clean_dna(seq)
    if (nchar(seq) >= 10L) {
      candidates <- c(candidates, seq)
    }
  }
  unique(candidates)
}

.dnmb_prophage_write_single_fasta_local <- function(header, sequence, path) {
  con <- base::file(path, "w")
  on.exit(base::close(con), add = TRUE)
  base::writeLines(paste0(">", header), con = con)
  base::writeLines(sequence, con = con)
  invisible(path)
}

.dnmb_prophage_blast_subject <- function(query_seq,
                                         subject_fasta,
                                         out_tsv,
                                         min_pident = 85,
                                         min_length = 200L) {
  query_seq <- .dnmb_prophage_clean_dna(query_seq)
  if (!nzchar(query_seq) || !file.exists(subject_fasta)) {
    return(data.frame())
  }
  blastn <- dnmb_detect_binary("blastn", required = FALSE)
  if (!isTRUE(blastn$found)) {
    return(data.frame())
  }
  query_fasta <- tempfile("dnmb-att-query-", fileext = ".fna")
  .dnmb_prophage_write_single_fasta_local("query", query_seq, query_fasta)
  run <- dnmb_run_external(
    blastn$path,
    args = c(
      "-task", "megablast",
      "-query", query_fasta,
      "-subject", subject_fasta,
      "-dust", "no",
      "-max_hsps", "50",
      "-outfmt", "6 sseqid sstart send pident length qstart qend sstrand bitscore",
      "-out", out_tsv
    ),
    required = FALSE
  )
  unlink(query_fasta, force = TRUE)
  if (!isTRUE(run$ok) || !file.exists(out_tsv) || !isTRUE(file.info(out_tsv)$size > 0)) {
    return(data.frame())
  }
  tbl <- utils::read.table(
    out_tsv,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  if (!nrow(tbl)) {
    return(data.frame())
  }
  names(tbl) <- c("sseqid", "sstart", "send", "pident", "length", "qstart", "qend", "sstrand", "bitscore")
  tbl$sstart <- suppressWarnings(as.integer(tbl$sstart))
  tbl$send <- suppressWarnings(as.integer(tbl$send))
  tbl$pident <- suppressWarnings(as.numeric(tbl$pident))
  tbl$length <- suppressWarnings(as.integer(tbl$length))
  tbl$qstart <- suppressWarnings(as.integer(tbl$qstart))
  tbl$qend <- suppressWarnings(as.integer(tbl$qend))
  tbl$bitscore <- suppressWarnings(as.numeric(tbl$bitscore))
  tbl$hit_start <- pmin(tbl$sstart, tbl$send)
  tbl$hit_end <- pmax(tbl$sstart, tbl$send)
  tbl$strand <- ifelse(tbl$send >= tbl$sstart, "+", "-")
  tbl <- tbl[!is.na(tbl$pident) & tbl$pident >= min_pident &
               !is.na(tbl$length) & tbl$length >= min(as.integer(min_length), nchar(query_seq)),
             , drop = FALSE]
  if (!nrow(tbl)) {
    return(data.frame())
  }
  tbl[order(tbl$sseqid, -tbl$bitscore, -tbl$length, -tbl$pident), , drop = FALSE]
}

.dnmb_prophage_select_best_related_pair <- function(left_hits,
                                                    right_hits,
                                                    max_empty_bp = 5000L,
                                                    max_overlap_bp = 200L) {
  if (!is.data.frame(left_hits) || !nrow(left_hits) || !is.data.frame(right_hits) || !nrow(right_hits)) {
    return(NULL)
  }
  best <- NULL
  best_score <- -Inf
  for (i in seq_len(nrow(left_hits))) {
    for (j in seq_len(nrow(right_hits))) {
      if (!identical(left_hits$sseqid[[i]], right_hits$sseqid[[j]])) {
        next
      }
      if (!identical(left_hits$strand[[i]], right_hits$strand[[j]])) {
        next
      }
      strand <- left_hits$strand[[i]]
      if (identical(strand, "+")) {
        if (left_hits$hit_start[[i]] >= right_hits$hit_start[[j]]) {
          next
        }
        empty_start <- left_hits$hit_end[[i]] + 1L
        empty_end <- right_hits$hit_start[[j]] - 1L
      } else {
        if (left_hits$hit_start[[i]] <= right_hits$hit_start[[j]]) {
          next
        }
        empty_start <- right_hits$hit_end[[j]] + 1L
        empty_end <- left_hits$hit_start[[i]] - 1L
      }
      gap_bp <- empty_end - empty_start + 1L
      if (gap_bp > as.integer(max_empty_bp) || gap_bp < (-as.integer(max_overlap_bp))) {
        next
      }
      score <- sum(
        left_hits$bitscore[[i]],
        right_hits$bitscore[[j]],
        na.rm = TRUE
      ) - max(gap_bp, 0L) / 10
      if (is.na(score) || score <= best_score) {
        next
      }
      best_score <- score
      best <- list(
        sseqid = left_hits$sseqid[[i]],
        strand = strand,
        left_start = left_hits$hit_start[[i]],
        left_end = left_hits$hit_end[[i]],
        right_start = right_hits$hit_start[[j]],
        right_end = right_hits$hit_end[[j]],
        empty_start = empty_start,
        empty_end = empty_end,
        gap_bp = gap_bp,
        score = score
      )
    }
  }
  best
}

.dnmb_prophage_match_record_by_id <- function(record_id, records) {
  record_id <- as.character(record_id)[1]
  if (!nzchar(record_id) || !nrow(records)) {
    return(NA_integer_)
  }
  idx <- match(record_id, records$accession)
  if (!is.na(idx)) {
    return(idx)
  }
  idx <- match(record_id, records$locus)
  if (!is.na(idx)) {
    return(idx)
  }
  idx <- match(record_id, records$definition)
  if (!is.na(idx)) {
    return(idx)
  }
  NA_integer_
}

.dnmb_prophage_extract_related_junction <- function(records, pair, pad_bp = 80L) {
  idx <- .dnmb_prophage_match_record_by_id(pair$sseqid, records)
  if (is.na(idx)) {
    return(NULL)
  }
  seq_full <- .dnmb_prophage_clean_dna(records$sequence[[idx]])
  if (!nzchar(seq_full)) {
    return(NULL)
  }
  seq_len <- nchar(seq_full)
  overlap_start <- max(pair$left_start, pair$right_start)
  overlap_end <- min(pair$left_end, pair$right_end)
  if (identical(pair$strand, "+")) {
    junction_start <- max(1L, pair$left_end - as.integer(pad_bp))
    junction_end <- min(seq_len, pair$right_start + as.integer(pad_bp))
    junction_seq <- substr(seq_full, junction_start, junction_end)
    attb_raw <- if (pair$empty_start <= pair$empty_end) {
      substr(seq_full, max(1L, pair$empty_start), min(seq_len, pair$empty_end))
    } else if (overlap_start <= overlap_end) {
      substr(seq_full, overlap_start, overlap_end)
    } else {
      ""
    }
  } else {
    junction_start <- max(1L, pair$right_end - as.integer(pad_bp))
    junction_end <- min(seq_len, pair$left_start + as.integer(pad_bp))
    junction_seq <- .dnmb_prophage_revcomp(substr(seq_full, junction_start, junction_end))
    attb_raw <- if (pair$empty_start <= pair$empty_end) {
      .dnmb_prophage_revcomp(substr(seq_full, max(1L, pair$empty_start), min(seq_len, pair$empty_end)))
    } else if (overlap_start <= overlap_end) {
      .dnmb_prophage_revcomp(substr(seq_full, overlap_start, overlap_end))
    } else {
      ""
    }
  }
  list(
    accession = records$accession[[idx]] %||% records$locus[[idx]],
    definition = records$definition[[idx]],
    junction_seq = junction_seq,
    attB_raw = attb_raw
  )
}

.dnmb_prophage_guess_organism_name <- function(genbank_table, contig = NULL) {
  if (is.null(genbank_table) || !is.data.frame(genbank_table) || !nrow(genbank_table)) {
    return(NA_character_)
  }
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  if (!is.null(contig) && "contig" %in% names(tbl)) {
    tbl <- tbl[tbl$contig == contig, , drop = FALSE]
  }
  if (!nrow(tbl)) {
    return(NA_character_)
  }
  for (col in c("organism", "source_organism", "species", "definition")) {
    if (!col %in% names(tbl)) {
      next
    }
    vals <- unique(trimws(as.character(tbl[[col]])))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (!length(vals)) {
      next
    }
    picked <- vals[[1]]
    tokens <- strsplit(picked, "\\s+")[[1]]
    if (length(tokens) >= 2L) {
      return(paste(tokens[[1]], tokens[[2]]))
    }
    return(picked)
  }
  NA_character_
}

.dnmb_prophage_recursive_pick <- function(x, keys) {
  if (is.null(x)) {
    return(NA_character_)
  }
  if (!is.list(x)) {
    return(NA_character_)
  }
  nms <- names(x)
  if (!is.null(nms) && length(nms)) {
    hit <- which(tolower(nms) %in% tolower(keys))
    if (length(hit)) {
      val <- x[[hit[[1]]]]
      if (length(val)) {
        return(as.character(val[[1]]))
      }
    }
  }
  for (item in x) {
    out <- .dnmb_prophage_recursive_pick(item, keys)
    if (!is.na(out) && nzchar(out)) {
      return(out)
    }
  }
  NA_character_
}

.dnmb_prophage_ncbi_fetch_json <- function(url,
                                           email = Sys.getenv("NCBI_EMAIL"),
                                           api_key = Sys.getenv("NCBI_API_KEY")) {
  params <- c(tool = "DNMB")
  if (!is.na(email) && nzchar(email)) {
    params[["email"]] <- email
  }
  if (!is.na(api_key) && nzchar(api_key)) {
    params[["api_key"]] <- api_key
  }
  if (length(params)) {
    url <- paste0(
      url,
      if (grepl("\\?", url, fixed = TRUE)) "&" else "?",
      paste(
        paste0(names(params), "=", utils::URLencode(unname(params), reserved = TRUE)),
        collapse = "&"
      )
    )
  }
  dest <- tempfile("dnmb-ncbi-", fileext = ".json")
  dl <- .dnmb_download_asset(url, dest, insecure = FALSE)
  if (!isTRUE(dl$ok) || !file.exists(dest) || !isTRUE(file.info(dest)$size > 0)) {
    return(list(ok = FALSE, data = NULL))
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(dest, simplifyVector = FALSE),
    error = function(e) NULL
  )
  list(ok = !is.null(parsed), data = parsed)
}

.dnmb_prophage_ncbi_related_assemblies <- function(organism_name,
                                                   max_records = 10L) {
  organism_name <- trimws(as.character(organism_name)[1])
  if (is.na(organism_name) || !nzchar(organism_name)) {
    return(tibble::tibble())
  }
  term <- paste0(
    "(", organism_name, "[Organism]) AND latest[filter] AND ",
    "(\"complete genome\"[Assembly Level] OR chromosome[Assembly Level])"
  )
  search_url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "?db=assembly&retmode=json&retmax=",
    as.integer(max_records),
    "&term=",
    utils::URLencode(term, reserved = TRUE)
  )
  search <- .dnmb_prophage_ncbi_fetch_json(search_url)
  ids <- search$data$esearchresult$idlist %||% list()
  ids <- unlist(ids, use.names = FALSE)
  ids <- ids[nzchar(ids)]
  if (!length(ids)) {
    return(tibble::tibble())
  }
  summary_url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
    "?db=assembly&retmode=json&id=",
    paste(ids, collapse = ",")
  )
  summary <- .dnmb_prophage_ncbi_fetch_json(summary_url)
  if (!isTRUE(summary$ok) || is.null(summary$data$result)) {
    return(tibble::tibble())
  }
  result <- summary$data$result
  rows <- lapply(ids, function(id) {
    rec <- result[[id]]
    if (is.null(rec)) {
      return(NULL)
    }
    tibble::tibble(
      assembly_uid = id,
      assembly_accession = .dnmb_prophage_recursive_pick(rec, c("assemblyaccession", "assembly_accession")),
      organism_name = .dnmb_prophage_recursive_pick(rec, c("organism", "organism_name", "speciesname")),
      assembly_level = .dnmb_prophage_recursive_pick(rec, c("assemblylevel", "assembly_level")),
      ftp_path_refseq = .dnmb_prophage_recursive_pick(rec, c("ftppath_refseq", "ftp_path_refseq")),
      ftp_path_genbank = .dnmb_prophage_recursive_pick(rec, c("ftppath_genbank", "ftp_path_genbank"))
    )
  })
  out <- dplyr::bind_rows(rows)
  out <- out[!is.na(out$assembly_accession) & nzchar(out$assembly_accession), , drop = FALSE]
  out <- dplyr::distinct(out, .data$assembly_accession, .keep_all = TRUE)
  out
}

.dnmb_prophage_download_related_genbank <- function(assembly_row, out_dir) {
  ftp_path <- assembly_row$ftp_path_refseq[[1]]
  if (is.na(ftp_path) || !nzchar(ftp_path)) {
    ftp_path <- assembly_row$ftp_path_genbank[[1]]
  }
  if (is.na(ftp_path) || !nzchar(ftp_path)) {
    return("")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ftp_http <- sub("^ftp://", "https://", ftp_path)
  ftp_base <- basename(ftp_path)
  gz_url <- paste0(ftp_http, "/", ftp_base, "_genomic.gbff.gz")
  dest_gz <- file.path(out_dir, paste0(ftp_base, "_genomic.gbff.gz"))
  dest_gbff <- sub("\\.gz$", "", dest_gz)
  if (file.exists(dest_gbff) && isTRUE(file.info(dest_gbff)$size > 0)) {
    return(dest_gbff)
  }
  dl <- .dnmb_download_asset(gz_url, dest_gz, insecure = FALSE)
  if (!isTRUE(dl$ok) || !file.exists(dest_gz)) {
    return("")
  }
  gunzip_run <- dnmb_run_external("gunzip", args = c("-f", dest_gz), required = FALSE)
  if (!isTRUE(gunzip_run$ok) || !file.exists(dest_gbff)) {
    lines <- tryCatch(readLines(gzfile(dest_gz), warn = FALSE), error = function(e) character())
    if (length(lines)) {
      writeLines(lines, dest_gbff)
    }
  }
  if (file.exists(dest_gbff) && isTRUE(file.info(dest_gbff)$size > 0)) {
    dest_gbff
  } else {
    ""
  }
}

.dnmb_prophage_fetch_related_genomes <- function(organism_name,
                                                 max_genomes = 3L,
                                                 cache_dir = NULL) {
  organism_name <- trimws(as.character(organism_name)[1])
  if (is.na(organism_name) || !nzchar(organism_name)) {
    return(character())
  }
  dest_dir <- if (!is.null(cache_dir) && nzchar(cache_dir)) {
    file.path(cache_dir, "related_genomes")
  } else {
    file.path(tempdir(), "dnmb-related-genomes")
  }
  tbl <- .dnmb_prophage_ncbi_related_assemblies(organism_name, max_records = max(6L, as.integer(max_genomes) * 3L))
  if (!nrow(tbl)) {
    return(character())
  }
  tbl <- tbl[order(
    !is.na(tbl$ftp_path_refseq) & nzchar(tbl$ftp_path_refseq),
    !is.na(tbl$assembly_level) & tbl$assembly_level == "Complete Genome",
    decreasing = TRUE
  ), , drop = FALSE]
  paths <- character()
  for (i in seq_len(nrow(tbl))) {
    path <- .dnmb_prophage_download_related_genbank(tbl[i, , drop = FALSE], dest_dir)
    if (nzchar(path)) {
      paths <- c(paths, path)
    }
    if (length(paths) >= as.integer(max_genomes)) {
      break
    }
  }
  unique(paths[nzchar(paths)])
}

.dnmb_prophage_comparative_attB_from_related <- function(left_flank,
                                                         right_flank,
                                                         windows,
                                                         related_genbank,
                                                         min_len = 10L,
                                                         max_len = 50L) {
  records <- .dnmb_prophage_parse_genbank_records(related_genbank)
  if (!nrow(records)) {
    return(NULL)
  }
  subject_fasta <- tempfile("dnmb-att-related-", fileext = ".fna")
  .dnmb_prophage_write_contig_fasta(records, subject_fasta)
  left_hits <- .dnmb_prophage_blast_subject(
    left_flank,
    subject_fasta,
    out_tsv = tempfile("dnmb-att-left-", fileext = ".tsv"),
    min_pident = 85,
    min_length = max(80L, floor(nchar(left_flank) * 0.25))
  )
  right_hits <- .dnmb_prophage_blast_subject(
    right_flank,
    subject_fasta,
    out_tsv = tempfile("dnmb-att-right-", fileext = ".tsv"),
    min_pident = 85,
    min_length = max(80L, floor(nchar(right_flank) * 0.25))
  )
  unlink(subject_fasta, force = TRUE)
  pair <- .dnmb_prophage_select_best_related_pair(left_hits, right_hits, max_empty_bp = 5000L, max_overlap_bp = 250L)
  if (is.null(pair)) {
    return(NULL)
  }
  junction <- .dnmb_prophage_extract_related_junction(records, pair, pad_bp = 80L)
  if (is.null(junction) || (!nzchar(junction$junction_seq) && !nzchar(junction$attB_raw))) {
    return(NULL)
  }
  core <- .dnmb_prophage_best_core_from_source(
    paste(junction$attB_raw, junction$junction_seq, sep = ""),
    windows = windows,
    min_len = min_len,
    max_len = max_len
  )
  if (!nzchar(core$seq)) {
    return(NULL)
  }
  list(
    attB_sequence = if (nzchar(junction$attB_raw)) junction$attB_raw else core$seq,
    att_core_sequence = core$seq,
    att_core_length = core$len,
    related_genome = junction$accession,
    related_genome_title = junction$definition,
    gap_bp = pair$gap_bp,
    score = pair$score
  )
}

.dnmb_prophage_comparative_attB_remote_nt <- function(left_flank,
                                                      right_flank,
                                                      windows,
                                                      organism_name = NULL) {
  blastn <- dnmb_detect_binary("blastn", required = FALSE)
  if (!isTRUE(blastn$found)) {
    return(NULL)
  }
  query_locus <- paste0(.dnmb_prophage_clean_dna(left_flank), .dnmb_prophage_clean_dna(right_flank))
  if (nchar(query_locus) < 200L) {
    return(NULL)
  }
  tmp_dir <- tempdir()
  query_fasta <- file.path(tmp_dir, "attB_query.fna")
  result_tsv <- file.path(tmp_dir, "attB_blast_result.tsv")
  .dnmb_prophage_write_single_fasta_local("query_locus", query_locus, query_fasta)
  run <- dnmb_run_external(
    blastn$path,
    args = c(
      "-task", "megablast",
      "-db", "nt",
      "-remote",
      "-query", query_fasta,
      "-max_target_seqs", "10",
      "-evalue", "1e-20",
      "-outfmt", "6 saccver sstart send pident length qstart qend qcovs stitle",
      "-out", result_tsv
    ),
    required = FALSE
  )
  if (!isTRUE(run$ok) || !file.exists(result_tsv) || !isTRUE(file.info(result_tsv)$size > 0)) {
    return(NULL)
  }
  hits <- utils::read.table(
    result_tsv,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  if (!nrow(hits)) {
    return(NULL)
  }
  names(hits) <- c("saccver", "sstart", "send", "pident", "length", "qstart", "qend", "qcovs", "stitle")
  hits$pident <- suppressWarnings(as.numeric(hits$pident))
  hits$length <- suppressWarnings(as.integer(hits$length))
  hits$qstart <- suppressWarnings(as.integer(hits$qstart))
  hits$qend <- suppressWarnings(as.integer(hits$qend))
  if (!is.null(organism_name) && nzchar(organism_name)) {
    keep <- grepl(organism_name, hits$stitle, ignore.case = TRUE, fixed = FALSE)
    if (any(keep)) {
      hits <- hits[keep, , drop = FALSE]
    }
  }
  if (!nrow(hits)) {
    return(NULL)
  }
  junction <- nchar(left_flank)
  hits <- hits[hits$qstart <= junction & hits$qend >= junction, , drop = FALSE]
  if (!nrow(hits)) {
    return(NULL)
  }
  hits <- hits[order(-hits$length, -hits$pident), , drop = FALSE]
  best <- hits[1, , drop = FALSE]
  ref_start <- min(best$sstart[[1]], best$send[[1]])
  query_junction_offset <- junction - best$qstart[[1]] + 1L
  ref_junction <- ref_start + query_junction_offset - 1L
  ref_url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    "?db=nuccore&id=",
    best$saccver[[1]],
    "&rettype=fasta&seq_start=",
    max(1L, ref_junction - 80L),
    "&seq_stop=",
    ref_junction + 80L
  )
  ref_fasta <- file.path(tmp_dir, "attB_ref.fna")
  dl <- .dnmb_download_asset(ref_url, ref_fasta, insecure = FALSE)
  if (!isTRUE(dl$ok) || !file.exists(ref_fasta) || !isTRUE(file.info(ref_fasta)$size > 0)) {
    return(NULL)
  }
  ref_lines <- readLines(ref_fasta, warn = FALSE)
  ref_seq <- .dnmb_prophage_clean_dna(paste(ref_lines[!grepl("^>", ref_lines)], collapse = ""))
  core <- .dnmb_prophage_best_core_from_source(ref_seq, windows = windows, min_len = 10L, max_len = 50L)
  if (!nzchar(core$seq)) {
    return(NULL)
  }
  list(
    attB_sequence = core$seq,
    att_core_sequence = core$seq,
    att_core_length = core$len,
    related_genome = as.character(best$saccver[[1]]),
    related_genome_title = as.character(best$stitle[[1]])
  )
}

.dnmb_prophage_comparative_attB <- function(genome_seq,
                                            pp_start,
                                            pp_end,
                                            contig,
                                            genbank_table = NULL,
                                            related_genbanks = NULL,
                                            organism_name = NULL,
                                            cache_dir = NULL,
                                            max_related_genomes = 3L,
                                            allow_remote_blast = TRUE) {
  genome_seq <- .dnmb_prophage_clean_dna(genome_seq)
  if (!nzchar(genome_seq)) {
    return(NULL)
  }
  left_flank <- substr(genome_seq, max(1L, as.integer(pp_start) - 2000L), max(1L, as.integer(pp_start) - 1L))
  right_flank <- substr(genome_seq, min(nchar(genome_seq), as.integer(pp_end) + 1L), min(nchar(genome_seq), as.integer(pp_end) + 2000L))
  if (nchar(left_flank) < 50L || nchar(right_flank) < 50L) {
    return(NULL)
  }
  windows <- .dnmb_prophage_boundary_windows(genome_seq, pp_start = pp_start, pp_end = pp_end, flank_bp = 500L)
  if (is.null(organism_name) || !nzchar(organism_name)) {
    organism_name <- .dnmb_prophage_guess_organism_name(genbank_table, contig = contig)
  }
  related_genbanks <- as.character(related_genbanks %||% character())
  related_genbanks <- related_genbanks[file.exists(related_genbanks)]
  if (!length(related_genbanks) && !is.null(organism_name) && nzchar(organism_name)) {
    related_genbanks <- .dnmb_prophage_fetch_related_genomes(
      organism_name = organism_name,
      max_genomes = max_related_genomes,
      cache_dir = cache_dir
    )
  }
  if (length(related_genbanks)) {
    for (path in related_genbanks[seq_len(min(length(related_genbanks), as.integer(max_related_genomes)))]) {
      local_hit <- tryCatch(
        .dnmb_prophage_comparative_attB_from_related(
          left_flank = left_flank,
          right_flank = right_flank,
          windows = windows,
          related_genbank = path
        ),
        error = function(e) NULL
      )
      if (!is.null(local_hit) && nzchar(local_hit$att_core_sequence %||% "")) {
        local_hit$method <- "local_related_genbank"
        return(local_hit)
      }
    }
  }
  if (!isTRUE(allow_remote_blast)) {
    return(NULL)
  }
  remote_hit <- tryCatch(
    .dnmb_prophage_comparative_attB_remote_nt(
      left_flank = left_flank,
      right_flank = right_flank,
      windows = windows,
      organism_name = organism_name
    ),
    error = function(e) NULL
  )
  if (!is.null(remote_hit)) {
    remote_hit$method <- "remote_nt_bootstrap"
  }
  remote_hit
}

.dnmb_prophage_detect_att_sites <- function(genome_seq,
                                            pp_start,
                                            pp_end,
                                            genbank_table = NULL,
                                            contig = NULL,
                                            cache_dir = NULL,
                                            related_genbanks = NULL,
                                            organism_name = NULL,
                                            direct_repeat_flank_bp = 500L,
                                            comparative = TRUE,
                                            comparative_max_genomes = 3L,
                                            allow_remote_blast = TRUE) {
  result <- .dnmb_prophage_empty_att_result()
  genome_seq <- .dnmb_prophage_clean_dna(genome_seq)
  pp_start <- suppressWarnings(as.integer(pp_start)[1])
  pp_end <- suppressWarnings(as.integer(pp_end)[1])
  if (!nzchar(genome_seq) || !is.finite(pp_start) || !is.finite(pp_end) || pp_start >= pp_end) {
    return(result)
  }
  windows <- .dnmb_prophage_boundary_windows(
    genome_seq,
    pp_start = pp_start,
    pp_end = pp_end,
    flank_bp = as.integer(direct_repeat_flank_bp)
  )
  candidates <- list()
  evidence <- character()
  add_candidate <- function(candidate) {
    if (is.null(candidate) || !nzchar(candidate$seq %||% "")) {
      return(invisible(NULL))
    }
    candidates[[length(candidates) + 1L]] <<- candidate
    evidence <<- unique(c(evidence, candidate$layer))
    invisible(candidate)
  }

  # Layer 1: direct repeats within +/-500 bp of each boundary.
  direct_repeat <- .dnmb_prophage_find_direct_repeats(
    windows$left$seq,
    windows$right$seq,
    min_len = 7L,
    max_len = 50L
  )
  if (direct_repeat$len >= 7L) {
    add_candidate(.dnmb_prophage_make_candidate(
      layer = "direct_repeat",
      seq = direct_repeat$seq,
      windows = windows,
      score = direct_repeat$score
    ))
  }

  # Layer 1b: Inverted repeats at boundaries
  ir <- .dnmb_prophage_find_inverted_repeats(windows$left$seq, windows$right$seq, min_len = 8L, max_len = 50L)
  if (ir$len >= 8L) {
    add_candidate(list(
      layer = "inverted_repeat",
      seq = ir$seq,
      len = ir$len,
      left_pos = if (!is.na(ir$left_pos)) windows$left$genome_start + ir$left_pos - 1L else NA_integer_,
      right_pos = if (!is.na(ir$right_pos)) windows$right$genome_start + ir$right_pos - 1L else NA_integer_,
      score = ir$len * 8L  # IR bonus
    ))
  }

  # Layer 1c: Tandem repeats at boundaries
  for (boundary_seq in list(windows$left$seq, windows$right$seq)) {
    tr <- .dnmb_prophage_find_tandem_repeats(boundary_seq, min_unit = 2L, max_unit = 8L, min_copies = 3L)
    if (nrow(tr)) {
      best_tr <- tr[1, ]  # longest tandem repeat
      add_candidate(list(
        layer = "tandem_repeat",
        seq = paste0("(", best_tr$unit, ")x", best_tr$copies),
        len = best_tr$total_len,
        score = best_tr$total_len * 5L
      ))
    }
  }

  # Layer 2: tRNA overlap/adjacency, using the tRNA 3' end as the candidate core.
  if (!is.null(genbank_table) && is.data.frame(genbank_table) && nrow(genbank_table) && !is.null(contig)) {
    gt <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
    gt$gs <- suppressWarnings(as.integer(gt$start))
    gt$ge <- suppressWarnings(as.integer(gt$end))
    if (!"strand" %in% names(gt) && "direction" %in% names(gt)) {
      gt$strand <- gt$direction
    }
    trna_tbl <- gt[gt$contig == contig & .dnmb_prophage_is_trna_feature(gt), , drop = FALSE]
    if (nrow(trna_tbl)) {
      near_boundary <- vapply(seq_len(nrow(trna_tbl)), function(i) {
        overlap_left <- trna_tbl$gs[[i]] <= pp_start && trna_tbl$ge[[i]] >= pp_start
        overlap_right <- trna_tbl$gs[[i]] <= pp_end && trna_tbl$ge[[i]] >= pp_end
        adjacent_left <- min(abs(trna_tbl$gs[[i]] - pp_start), abs(trna_tbl$ge[[i]] - pp_start)) <= 500L
        adjacent_right <- min(abs(trna_tbl$gs[[i]] - pp_end), abs(trna_tbl$ge[[i]] - pp_end)) <= 500L
        overlap_left || overlap_right || adjacent_left || adjacent_right
      }, logical(1))
      trna_tbl <- trna_tbl[near_boundary, , drop = FALSE]
      if (nrow(trna_tbl)) {
        for (i in seq_len(nrow(trna_tbl))) {
          trna_seq <- .dnmb_prophage_extract_feature_sequence(trna_tbl[i, , drop = FALSE], genome_seq)
          if (!nzchar(trna_seq)) {
            next
          }
          trna_3prime <- substr(trna_seq, max(1L, nchar(trna_seq) - 24L), nchar(trna_seq))
          core <- .dnmb_prophage_best_core_from_source(trna_3prime, windows = windows, min_len = 10L, max_len = 25L)
          if (!nzchar(core$seq)) {
            next
          }
          boundary_bonus <- if (
            (trna_tbl$gs[[i]] <= pp_start && trna_tbl$ge[[i]] >= pp_start) ||
            (trna_tbl$gs[[i]] <= pp_end && trna_tbl$ge[[i]] >= pp_end)
          ) {
            200
          } else {
            0
          }
          add_candidate(.dnmb_prophage_make_candidate(
            layer = "tRNA_overlap",
            seq = core$seq,
            windows = windows,
            score = core$score + boundary_bonus,
            attB_sequence = core$seq,
            source_note = trna_tbl$product[[i]] %||% NA_character_
          ))
        }
      }
    }
  }

  # Layer 3: cached reference phages, preferring annotated attP/repeat regions.
  if (!is.null(cache_dir) && dir.exists(cache_dir)) {
    ref_files <- list.files(
      cache_dir,
      pattern = "\\.(gbk|gbff|gb)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    for (ref_file in ref_files) {
      ref_records <- .dnmb_prophage_parse_genbank_records(ref_file)
      if (!nrow(ref_records)) {
        next
      }
      ref_seq <- .dnmb_prophage_clean_dna(ref_records$sequence[[1]])
      if (!nzchar(ref_seq)) {
        next
      }
      ref_candidates <- .dnmb_prophage_reference_att_annotations(ref_file, ref_seq, terminal_bp = 500L)
      if (!length(ref_candidates)) {
        termini_repeat <- .dnmb_prophage_find_direct_repeats(
          substr(ref_seq, 1L, min(500L, nchar(ref_seq))),
          substr(ref_seq, max(1L, nchar(ref_seq) - 499L), nchar(ref_seq)),
          min_len = 10L,
          max_len = 50L
        )
        if (termini_repeat$len >= 10L) {
          ref_candidates <- termini_repeat$seq
        }
      }
      if (!length(ref_candidates)) {
        next
      }
      for (candidate_seq in unique(ref_candidates)) {
        core <- .dnmb_prophage_best_core_from_source(candidate_seq, windows = windows, min_len = 10L, max_len = 50L)
        if (!nzchar(core$seq)) {
          next
        }
        add_candidate(.dnmb_prophage_make_candidate(
          layer = "reference_attP",
          seq = core$seq,
          windows = windows,
          score = core$score + 100,
          source_note = basename(ref_file)
        ))
      }
    }
  }

  # Layer 4: comparative empty-site recovery; local related genomes are preferred.
  if (isTRUE(comparative)) {
    comparative_hit <- tryCatch(
      .dnmb_prophage_comparative_attB(
        genome_seq = genome_seq,
        pp_start = pp_start,
        pp_end = pp_end,
        contig = contig,
        genbank_table = genbank_table,
        related_genbanks = related_genbanks,
        organism_name = organism_name,
        cache_dir = cache_dir,
        max_related_genomes = comparative_max_genomes,
        allow_remote_blast = allow_remote_blast
      ),
      error = function(e) NULL
    )
    if (!is.null(comparative_hit) && nzchar(comparative_hit$att_core_sequence %||% "")) {
      comparative_layer <- if (identical(comparative_hit$method %||% "", "local_related_genbank")) {
        "comparative_genomics"
      } else {
        "comparative_remote_bootstrap"
      }
      add_candidate(.dnmb_prophage_make_candidate(
        layer = comparative_layer,
        seq = comparative_hit$att_core_sequence,
        windows = windows,
        score = (comparative_hit$att_core_length %||% nchar(comparative_hit$att_core_sequence)) * 140,
        attB_sequence = comparative_hit$attB_sequence %||% comparative_hit$att_core_sequence,
        related_genome = comparative_hit$related_genome %||% NA_character_,
        source_note = comparative_hit$related_genome_title %||% NA_character_
      ))
    }
  }

  result$evidence_layers <- unique(evidence)
  if (!length(candidates)) {
    return(result)
  }

  layer_priority <- c(
    comparative_genomics = 5,
    reference_attP = 4,
    tRNA_overlap = 3,
    direct_repeat = 2,
    comparative_remote_bootstrap = 1
  )
  candidate_tbl <- lapply(seq_along(candidates), function(i) {
    cand <- candidates[[i]]
    support_layers <- unique(vapply(candidates, function(other) {
      same_core <- identical(cand$seq, other$seq) ||
        grepl(cand$seq, other$seq, fixed = TRUE) ||
        grepl(other$seq, cand$seq, fixed = TRUE)
      if (isTRUE(same_core)) other$layer else NA_character_
    }, character(1)))
    support_layers <- support_layers[!is.na(support_layers) & nzchar(support_layers)]
    data.frame(
      idx = i,
      layer = cand$layer,
      seq = cand$seq,
      len = cand$len,
      score = cand$score,
      support_count = length(support_layers),
      priority = unname(layer_priority[cand$layer] %||% 0),
      stringsAsFactors = FALSE
    )
  })
  candidate_tbl <- dplyr::bind_rows(candidate_tbl)
  candidate_tbl <- candidate_tbl[order(
    -candidate_tbl$support_count,
    -candidate_tbl$priority,
    -candidate_tbl$score,
    -candidate_tbl$len,
    candidate_tbl$idx
  ), , drop = FALSE]
  best <- candidates[[candidate_tbl$idx[[1]]]]

  result$att_core_sequence <- best$seq
  result$att_core_length <- as.integer(best$len)
  result$attL_sequence <- best$seq
  result$attL_start <- as.integer(best$attL_start)
  result$attL_end <- as.integer(best$attL_end)
  result$attR_sequence <- best$seq
  result$attR_start <- as.integer(best$attR_start)
  result$attR_end <- as.integer(best$attR_end)

  attb_candidates <- vapply(candidates, function(x) x$attB_sequence %||% NA_character_, character(1))
  attb_candidates <- unique(attb_candidates[!is.na(attb_candidates) & nzchar(attb_candidates)])
  if (length(attb_candidates)) {
    result$attB_sequence <- attb_candidates[[1]]
  }
  related_candidates <- vapply(candidates, function(x) x$related_genome %||% NA_character_, character(1))
  related_candidates <- unique(related_candidates[!is.na(related_candidates) & nzchar(related_candidates)])
  if (length(related_candidates)) {
    result$related_genome_used <- related_candidates[[1]]
  }

  has_local_comparative <- "comparative_genomics" %in% result$evidence_layers
  has_remote_bootstrap <- "comparative_remote_bootstrap" %in% result$evidence_layers
  if (has_local_comparative && candidate_tbl$support_count[[1]] >= 2L) {
    result$confidence <- "high"
  } else if (has_local_comparative || length(result$evidence_layers) >= 2L) {
    result$confidence <- "medium"
  } else if (has_remote_bootstrap || length(result$evidence_layers) == 1L) {
    result$confidence <- "low"
  } else {
    result$confidence <- "none"
  }
  result
}
#' Internal prophage module helpers
#'
#' Backend orchestration, parsing, summarization, and cache helpers for the
#' DNMB prophage workflow.
#'
#' @name dnmb_internal_prophage_module
#' @keywords internal
#' @noRd
NULL

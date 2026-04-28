.dnmb_acrfinder_module_name <- function() {
  "acrfinder"
}

.dnmb_acrfinder_default_version <- function() {
  "current"
}

.dnmb_acrfinder_default_repo_url <- function() {
  "https://github.com/HaidYi/acrfinder.git"
}

.dnmb_acrfinder_status_row <- function(component, status, detail = NA_character_) {
  tibble::tibble(
    component = base::as.character(component)[1],
    status = base::as.character(status)[1],
    detail = base::as.character(detail)[1]
  )
}

.dnmb_acrfinder_empty_status <- function() {
  .dnmb_acrfinder_status_row(character(), character(), character())
}

.dnmb_acrfinder_trace <- function(path, text) {
  base::cat(base::paste0(text, "\n"), file = path, append = TRUE)
}

.dnmb_acrfinder_trace_result <- function(path, result, label = "acrfinder_run") {
  if (base::is.null(result) || !base::is.list(result)) {
    return(invisible(FALSE))
  }
  lines <- c(
    base::sprintf("[%s] %s status=%s ok=%s", base::Sys.time(), label, result$status %||% NA_integer_, base::isTRUE(result$ok)),
    if (!base::is.null(result$error) && base::nzchar(base::as.character(result$error))) {
      base::paste0("error: ", base::as.character(result$error))
    } else {
      character()
    },
    if (base::length(result$stderr)) {
      c("stderr:", result$stderr)
    } else {
      character()
    },
    if (base::length(result$stdout)) {
      c("stdout:", result$stdout)
    } else {
      character()
    }
  )
  .dnmb_acrfinder_trace(path, base::paste(lines, collapse = "\n"))
  invisible(TRUE)
}

.dnmb_acrfinder_patch_biopython_compat <- function(repo_dir) {
  mask_path <- base::file.path(repo_dir, "mask_fna_with_spacers.py")
  changed <- FALSE
  if (base::file.exists(mask_path)) {
    lines <- base::readLines(mask_path, warn = FALSE)
    import_idx <- grep("^from Bio\\.Alphabet import SingleLetterAlphabet$", lines)
    if (base::length(import_idx)) {
      lines <- lines[-import_idx]
      changed <- TRUE
    }
    seq_idx <- grep("SingleLetterAlphabet\\(\\)", lines, fixed = FALSE)
    if (base::length(seq_idx)) {
      lines[seq_idx] <- gsub("Seq\\(([^,]+), SingleLetterAlphabet\\(\\)\\)", "Seq(\\1)", lines[seq_idx])
      changed <- TRUE
    }
    if (changed) {
      base::writeLines(lines, con = mask_path)
    }
  }

  runner_path <- base::file.path(repo_dir, "acr_aca_cri_runner.py")
  if (base::file.exists(runner_path)) {
    lines <- base::readLines(runner_path, warn = FALSE)
    if (!base::any(lines == "from shutil import which as shutil_which")) {
      import_anchor <- grep("^from os import path as os_path$", lines)
      if (base::length(import_anchor)) {
        lines <- append(lines, "from shutil import which as shutil_which", after = import_anchor[[1]])
        changed <- TRUE
      }
    }
    if (!base::any(lines == "def _dnmb_rpsblast_binary():")) {
      import_idx <- grep("^from shutil import which as shutil_which$", lines)
      if (base::length(import_idx)) {
        lines <- append(
          lines,
          c(
            "",
            "def _dnmb_rpsblast_binary():",
            "    return shutil_which('rpsblast+') or shutil_which('rpsblast')"
          ),
          after = import_idx[[1]]
        )
        changed <- TRUE
      }
    }
    selector_idx <- grep("CRISPR_CAS_FINDER_EXECUTABLE, CRISPR_CAS_FINDER_SO, NUM_CPUS = ", lines, fixed = TRUE)
    if (base::length(selector_idx) && base::grepl("dependencies/CRISPRCasFinder/sel392v2\\.so", lines[[selector_idx[[1]]]])) {
      lines[[selector_idx[[1]]]] <- "CRISPR_CAS_FINDER_EXECUTABLE, CRISPR_CAS_FINDER_SO, NUM_CPUS = 'dependencies/CRISPRCasFinder/CRISPRCasFinder.pl', '/usr/lib/vmatch/SELECT/sel392.so', '4'  # CRISPRCasFinder files"
      changed <- TRUE
    }
    rpsblast_exec_idx <- grep("execute\\(\\['rpsblast\\+',", lines)
    if (base::length(rpsblast_exec_idx)) {
      lines[rpsblast_exec_idx] <- sub("execute\\(\\['rpsblast\\+',", "execute([_dnmb_rpsblast_binary(),", lines[rpsblast_exec_idx])
      changed <- TRUE
    }
    cdd_block_start <- grep("^\\tskip_cdd_filter = ", lines)
    if (!base::length(cdd_block_start)) {
      legacy_start <- grep("^\\t\\t*fasta_file, rpsblast_file = generate_filter_file\\(ACR_ACA_FILE, CDD_DBFILE, INTERMEDIATES\\)$", lines)
      if (base::length(legacy_start)) {
        cdd_block_start <- legacy_start[[1]]
      }
    } else {
      cdd_block_start <- cdd_block_start[[1]]
    }
    cdd_block_end <- grep("^\\twith open\\(ACR_ACA_FILE\\) as acrFile:", lines)
    if (base::length(cdd_block_start) == 1L && base::length(cdd_block_end) == 1L && cdd_block_end[[1]] > cdd_block_start) {
      replacement <- c(
        "\tskip_cdd_filter = (GENOME_TYPE == 'V') or (_dnmb_rpsblast_binary() is None) or (not os_path.exists(ESCAPE_DBFILE)) or (not (os_path.exists(CDD_DBFILE) or os_path.exists(CDD_DBFILE + '.pn')))",
        "\tif skip_cdd_filter:",
        "\t\tfasta_file, rpsblast_file = None, None",
        "\t\tescape_set = set()",
        "\telse:",
        "\t\tfasta_file, rpsblast_file = generate_filter_file(ACR_ACA_FILE, CDD_DBFILE, INTERMEDIATES)",
        "\t\t# obtain the escape list",
        "\t\tescape_set = set()",
        "\t\tfor line in open(ESCAPE_DBFILE).readlines():",
        "\t\t\tline = line.rstrip()",
        "\t\t\tif line:",
        "\t\t\t\tescape_set.add(line)"
      )
      current_block <- lines[cdd_block_start:(cdd_block_end[[1]] - 1L)]
      if (!identical(current_block, replacement)) {
        lines <- c(
          if (cdd_block_start > 1L) lines[seq_len(cdd_block_start - 1L)] else character(),
          replacement,
          lines[cdd_block_end[[1]]:base::length(lines)]
        )
        changed <- TRUE
      }
    }

    assign_idx <- grep("^\\t\\t\\tis_cddfilter =", lines)
    if (base::length(assign_idx)) {
      corrected_assignment <- "\t\t\tis_cddfilter = (GENOME_TYPE != 'V') and (not skip_cdd_filter)"
      if (!identical(lines[[assign_idx[[1]]]], corrected_assignment)) {
        lines[[assign_idx[[1]]]] <- corrected_assignment
        changed <- TRUE
      }
      drop_idx <- assign_idx[-1L]
      if (base::length(drop_idx)) {
        lines <- lines[-drop_idx]
        changed <- TRUE
      }
    }

    if (changed) {
      base::writeLines(lines, con = runner_path)
    }
  }

  crisprcas_path <- base::file.path(repo_dir, "dependencies", "CRISPRCasFinder", "CRISPRCasFinder.pl")
  if (base::file.exists(crisprcas_path)) {
    lines <- base::readLines(crisprcas_path, warn = FALSE)
    so_idx <- grep('^my \\$so = "\\./sel392v2\\.so";', lines)
    if (base::length(so_idx)) {
      lines[[so_idx[[1]]]] <- 'my $so = $ENV{"DNMB_VMATCH_SEL392"} || "/usr/lib/vmatch/SELECT/sel392.so"; # path to shared object (.so) file (former name: $pathSoFile)'
      changed <- TRUE
    }
    trunc_idx <- grep('elsif\\(\\$ARGV\\[\\$i\\]=~/-truncDR/ or \\$ARGV\\[\\$i\\]=~/-t/\\)\\{', lines)
    if (base::length(trunc_idx)) {
      lines[[trunc_idx[[1]]]] <- '    elsif($ARGV[$i]=~/^-truncDR$/ or $ARGV[$i]=~/^-t$/){'
      changed <- TRUE
    }
    if (base::isTRUE(changed)) {
      base::writeLines(lines, con = crisprcas_path)
    }
  }
  changed
}

.dnmb_acrfinder_asset_layout <- function(module_dir) {
  repo_dir <- base::file.path(module_dir, "acrfinder")
  data_dir <- base::file.path(repo_dir, "dependencies", "diamond_query")
  list(
    module_dir = module_dir,
    repo_dir = repo_dir,
    env_dir = base::file.path(module_dir, "venv"),
    env_python = base::file.path(module_dir, "venv", "bin", "python"),
    env_pip = base::file.path(module_dir, "venv", "bin", "pip"),
    acr_fasta = base::file.path(data_dir, "known-acr.faa"),
    aca_fasta = base::file.path(data_dir, "401-aca.faa")
  )
}

.dnmb_acrfinder_normalize_asset_urls <- function(asset_urls = NULL) {
  if (base::is.null(asset_urls)) {
    return(list())
  }
  if (!base::is.list(asset_urls)) {
    base::stop("`asset_urls` for AcrFinder must be a list when supplied.", call. = FALSE)
  }
  asset_urls
}

.dnmb_acrfinder_normalize_genome_type <- function(genome_type) {
  value <- base::toupper(base::trimws(base::as.character(genome_type)[1]))
  if (base::is.na(value) || !base::nzchar(value) || identical(value, "AUTO")) {
    return(NA_character_)
  }
  if (!value %in% c("A", "B", "V")) {
    base::stop("`genome_type` for AcrFinder must be one of: A, B, V, or AUTO.", call. = FALSE)
  }
  value
}

.dnmb_acrfinder_guess_genome_type <- function(genbank = NULL) {
  genbank_path <- base::trimws(base::as.character(genbank)[1])
  if (base::is.na(genbank_path) || !base::nzchar(genbank_path) || !base::file.exists(genbank_path)) {
    return("B")
  }

  header_lines <- tryCatch(
    base::readLines(genbank_path, warn = FALSE, n = 160L),
    error = function(e) character()
  )
  header_text <- base::paste(header_lines, collapse = "\n")
  if (!base::nzchar(header_text)) {
    return("B")
  }

  if (grepl("\\bArchaea\\b", header_text, ignore.case = TRUE, perl = TRUE)) {
    return("A")
  }
  if (grepl("\\bViruses\\b|\\bvirus\\b|\\bphage\\b", header_text, ignore.case = TRUE, perl = TRUE)) {
    return("V")
  }
  "B"
}

.dnmb_acrfinder_resolve_genome_type <- function(genome_type = NULL, genbank = NULL) {
  explicit_type <- .dnmb_acrfinder_normalize_genome_type(genome_type)
  if (!base::is.na(explicit_type) && base::nzchar(explicit_type)) {
    return(explicit_type)
  }
  .dnmb_acrfinder_guess_genome_type(genbank = genbank)
}

.dnmb_acrfinder_cached_install_state <- function(layout, manifest = NULL) {
  if (base::is.null(manifest) || !base::isTRUE(manifest$install_ok)) {
    return(list(ok = FALSE, detail = "AcrFinder manifest is missing or marked incomplete."))
  }

  if (!base::dir.exists(layout$repo_dir)) {
    return(list(ok = FALSE, detail = base::paste0("AcrFinder repo is missing: ", layout$repo_dir)))
  }

  required_paths <- c(layout$env_python, layout$acr_fasta, layout$aca_fasta)
  missing <- required_paths[!base::file.exists(required_paths)]
  if (base::length(missing)) {
    return(list(ok = FALSE, detail = base::paste0("AcrFinder runtime is missing: ", missing[[1]])))
  }

  py_run <- dnmb_run_external(
    layout$env_python,
    args = c("-c", "import sys, Bio; print(sys.executable)"),
    required = FALSE
  )
  if (!base::isTRUE(py_run$ok)) {
    return(list(ok = FALSE, detail = py_run$error %||% layout$env_python))
  }

  list(ok = TRUE, detail = layout$env_python)
}

.dnmb_acrfinder_install_ready <- function(layout, manifest = NULL) {
  base::isTRUE(.dnmb_acrfinder_cached_install_state(layout, manifest = manifest)$ok)
}

.dnmb_acrfinder_prepare_env <- function(layout, force = FALSE) {
  py <- .dnmb_defensefinder_candidate_python()
  if (!base::nzchar(py)) {
    return(.dnmb_acrfinder_status_row("acrfinder_python", "missing", "No python3 executable found in PATH."))
  }
  if (base::dir.exists(layout$env_dir) && base::isTRUE(force)) {
    unlink(layout$env_dir, recursive = TRUE, force = TRUE)
  }
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
    create_run <- dnmb_run_external(py, args = c("-m", "venv", layout$env_dir), required = FALSE)
    if (!base::isTRUE(create_run$ok)) {
      return(.dnmb_acrfinder_status_row("acrfinder_python", "failed", create_run$error %||% py))
    }
  }
  pip_run <- dnmb_run_external(layout$env_python, args = c("-m", "pip", "install", "--upgrade", "pip"), required = FALSE)
  bio_run <- dnmb_run_external(layout$env_python, args = c("-m", "pip", "install", "biopython"), required = FALSE)
  if (!base::isTRUE(bio_run$ok)) {
    return(.dnmb_acrfinder_status_row("acrfinder_python", "failed", bio_run$error %||% layout$env_python))
  }
  .dnmb_acrfinder_status_row(
    "acrfinder_python",
    if (base::isTRUE(pip_run$ok) || base::isTRUE(bio_run$ok)) "ok" else "failed",
    layout$env_python
  )
}

.dnmb_acrfinder_prepare_repo <- function(layout,
                                         repo_source = .dnmb_acrfinder_default_repo_url(),
                                         force = FALSE) {
  repo_source <- base::as.character(repo_source)[1]
  if (base::dir.exists(layout$repo_dir) && !base::isTRUE(force)) {
    return(.dnmb_acrfinder_status_row("acrfinder_repo", "cached", layout$repo_dir))
  }

  if (base::dir.exists(layout$repo_dir) && base::isTRUE(force)) {
    unlink(layout$repo_dir, recursive = TRUE, force = TRUE)
  }

  if (base::dir.exists(repo_source)) {
    .dnmb_defensefinder_copy_dir_contents(repo_source, layout$repo_dir)
    return(.dnmb_acrfinder_status_row("acrfinder_repo", "ok", layout$repo_dir))
  }

  run <- dnmb_run_external(
    "git",
    args = c("clone", "--depth", "1", repo_source, layout$repo_dir),
    required = FALSE
  )
  .dnmb_acrfinder_status_row(
    "acrfinder_repo",
    if (base::isTRUE(run$ok)) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    if (base::isTRUE(run$ok)) layout$repo_dir else (run$error %||% repo_source)
  )
}

dnmb_acrfinder_install_module <- function(version = .dnmb_acrfinder_default_version(),
                                          cache_root = NULL,
                                          install = TRUE,
                                          repo_url = .dnmb_acrfinder_default_repo_url(),
                                          asset_urls = NULL,
                                          force = FALSE) {
  module <- .dnmb_acrfinder_module_name()
  module_dir <- .dnmb_db_module_dir(module, version, cache_root = cache_root, create = TRUE)
  layout <- .dnmb_acrfinder_asset_layout(module_dir)
  asset_urls <- .dnmb_acrfinder_normalize_asset_urls(asset_urls)
  status <- .dnmb_acrfinder_empty_status()
  manifest <- dnmb_db_read_manifest(module, version, cache_root = cache_root, required = FALSE)
  cache_state <- .dnmb_acrfinder_cached_install_state(layout, manifest = manifest)

  if (base::dir.exists(layout$repo_dir)) {
    patch_changed <- .dnmb_acrfinder_patch_biopython_compat(layout$repo_dir)
    if (base::isTRUE(patch_changed)) {
      status <- dplyr::bind_rows(
        status,
        .dnmb_acrfinder_status_row("acrfinder_repo_patch", "ok", "Applied Biopython/CDD compatibility patch.")
      )
    }
  }

  if (!base::isTRUE(cache_state$ok) && !base::isTRUE(force) && !base::is.null(manifest) && base::isTRUE(manifest$install_ok)) {
    status <- dplyr::bind_rows(
      status,
      .dnmb_acrfinder_status_row(
        "acrfinder_install",
        "stale",
        cache_state$detail %||% "Cached AcrFinder runtime is not executable in the current environment."
      )
    )
  }

  if (base::isTRUE(cache_state$ok) && !base::isTRUE(force)) {
    return(list(
      ok = TRUE,
      status = dplyr::bind_rows(status, .dnmb_acrfinder_status_row("acrfinder_install", "cached", module_dir)),
      files = list(acr_fasta = layout$acr_fasta, aca_fasta = layout$aca_fasta),
      manifest = manifest
    ))
  }

  if (!base::isTRUE(install)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_acrfinder_status_row("acrfinder_install", "missing", "AcrFinder is missing and module_install is FALSE.")),
      files = list(),
      manifest = manifest
    ))
  }

  repo_status <- .dnmb_acrfinder_prepare_repo(
    layout,
    repo_source = asset_urls$repo_dir %||% asset_urls$repo_url %||% repo_url,
    force = force
  )
  status <- dplyr::bind_rows(status, repo_status)
  if (!repo_status$status %in% c("ok", "cached")) {
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }

  if (!base::file.exists(layout$acr_fasta) || !base::file.exists(layout$aca_fasta)) {
    return(list(
      ok = FALSE,
      status = dplyr::bind_rows(status, .dnmb_acrfinder_status_row("acrfinder_assets", "failed", "AcrFinder diamond query FASTA files are missing.")),
      files = list(),
      manifest = manifest
    ))
  }

  .dnmb_acrfinder_patch_biopython_compat(layout$repo_dir)

  env_status <- .dnmb_acrfinder_prepare_env(layout, force = force)
  status <- dplyr::bind_rows(status, env_status)
  if (env_status$status != "ok") {
    return(list(ok = FALSE, status = status, files = list(), manifest = manifest))
  }

  manifest <- list(
    install_ok = TRUE,
    repo_dir = layout$repo_dir,
    env_python = layout$env_python,
    acr_fasta = layout$acr_fasta,
    aca_fasta = layout$aca_fasta
  )
  dnmb_db_write_manifest(module, version, manifest = manifest, cache_root = cache_root, overwrite = TRUE)
  .dnmb_db_autoprune_default_versions(
    module = module,
    version = version,
    default_version = .dnmb_acrfinder_default_version(),
    cache_root = cache_root
  )

  list(
    ok = TRUE,
    status = dplyr::bind_rows(status, .dnmb_acrfinder_status_row("acrfinder_install", "ok", module_dir)),
    files = list(acr_fasta = layout$acr_fasta, aca_fasta = layout$aca_fasta),
    manifest = dnmb_db_read_manifest(module, version, cache_root = cache_root, required = TRUE)
  )
}

dnmb_acrfinder_get_module <- function(version = .dnmb_acrfinder_default_version(),
                                      cache_root = NULL,
                                      required = TRUE) {
  manifest <- dnmb_db_read_manifest(.dnmb_acrfinder_module_name(), version, cache_root = cache_root, required = FALSE)
  if (base::is.null(manifest)) {
    if (base::isTRUE(required)) {
      base::stop("AcrFinder module is not installed.", call. = FALSE)
    }
    return(list(ok = FALSE, manifest = NULL))
  }
  list(
    ok = TRUE,
    manifest = manifest,
    files = list(
      env_python = manifest$env_python,
      acr_fasta = manifest$acr_fasta,
      aca_fasta = manifest$aca_fasta
    )
  )
}

.dnmb_acrfinder_header_label <- function(x) {
  x <- base::sub("^>", "", base::as.character(x))
  x <- base::trimws(x)
  x <- x[1]
  if (base::is.na(x) || !base::nzchar(x)) {
    return(NA_character_)
  }
  strsplit(x, "[[:space:]]+")[[1]][1]
}

.dnmb_acrfinder_db_headers <- function(path) {
  lines <- base::readLines(path, warn = FALSE)
  headers <- lines[grepl("^>", lines)]
  headers <- vapply(headers, .dnmb_acrfinder_header_label, character(1))
  headers[!base::is.na(headers) & base::nzchar(headers)]
}

.dnmb_acrfinder_make_db <- function(fasta_path, db_prefix) {
  run <- dnmb_run_external(
    "diamond",
    args = c("makedb", "--in", fasta_path, "--db", db_prefix),
    required = FALSE
  )
  list(
    ok = base::isTRUE(run$ok),
    status = if (base::isTRUE(run$ok)) "ok" else if (!base::nzchar(run$resolved_command)) "missing" else "failed",
    detail = if (base::isTRUE(run$ok)) db_prefix else (run$error %||% db_prefix)
  )
}

.dnmb_acrfinder_empty_hits <- function() {
  tibble::tibble(
    query = character(),
    acr_id = character(),
    pident = numeric(),
    qcov = numeric(),
    evalue = numeric(),
    bitscore = numeric(),
    align_len = integer()
  )
}

dnmb_acrfinder_parse_diamond <- function(path,
                                         evalue_threshold = 1e-5,
                                         identity_threshold = 30,
                                         coverage_threshold = 0.8) {
  if (!base::file.exists(path) || !isTRUE(base::file.info(path)$size > 0)) {
    return(.dnmb_acrfinder_empty_hits())
  }
  tbl <- tryCatch(
    utils::read.delim(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (!base::nrow(tbl)) {
    return(.dnmb_acrfinder_empty_hits())
  }
  base::colnames(tbl) <- c("acr_id", "query", "pident", "align_len", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcov")
  tbl$query <- .dnmb_module_clean_annotation_key(tbl$query)
  tbl$pident <- suppressWarnings(base::as.numeric(tbl$pident))
  tbl$qcov <- suppressWarnings(base::as.numeric(tbl$qcov))
  tbl$evalue <- suppressWarnings(base::as.numeric(tbl$evalue))
  tbl$bitscore <- suppressWarnings(base::as.numeric(tbl$bitscore))
  tbl$align_len <- suppressWarnings(base::as.integer(tbl$align_len))
  tbl <- tbl[
    !base::is.na(tbl$evalue) &
      tbl$evalue <= evalue_threshold &
      !base::is.na(tbl$pident) &
      tbl$pident >= identity_threshold &
      !base::is.na(tbl$qcov) &
      tbl$qcov >= (coverage_threshold * 100),
    ,
    drop = FALSE
  ]
  if (!base::nrow(tbl)) {
    return(.dnmb_acrfinder_empty_hits())
  }
  tbl <- tbl[base::order(tbl$query, tbl$evalue, -tbl$bitscore, -tbl$qcov), , drop = FALSE]
  tbl <- tbl[!base::duplicated(tbl$query), c("query", "acr_id", "pident", "qcov", "evalue", "bitscore", "align_len"), drop = FALSE]
  base::rownames(tbl) <- NULL
  tibble::as_tibble(tbl)
}

dnmb_acrfinder_normalize_hits <- function(hits) {
  if (base::is.null(hits) || !base::nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  support <- base::paste0(
    "acr=", hits$acr_id,
    "; pident=", base::round(hits$pident, 2),
    "; qcov=", base::round(hits$qcov, 2),
    "; evalue=", signif(hits$evalue, 4),
    "; bitscore=", base::round(hits$bitscore, 2)
  )
  out <- base::data.frame(
    query = .dnmb_module_clean_annotation_key(hits$query),
    source = base::rep("acrfinder", base::nrow(hits)),
    family_system = base::rep("AcrFinder", base::nrow(hits)),
    family_id = base::as.character(hits$acr_id),
    hit_label = base::as.character(hits$acr_id),
    enzyme_role = base::rep("CRISPR-Cas", base::nrow(hits)),
    evidence_mode = base::rep("direct", base::nrow(hits)),
    substrate_label = base::rep("anti-CRISPR", base::nrow(hits)),
    support = support,
    acr_id = base::as.character(hits$acr_id),
    pident = base::as.numeric(hits$pident),
    qcov = base::as.numeric(hits$qcov),
    evalue = base::as.numeric(hits$evalue),
    bitscore = base::as.numeric(hits$bitscore),
    align_len = base::as.integer(hits$align_len),
    typing_eligible = base::rep(TRUE, base::nrow(hits)),
    stringsAsFactors = FALSE
  )
  out[, c(.dnmb_module_optional_long_columns(), base::setdiff(base::names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
}

.dnmb_acrfinder_parse_hit_string <- function(x) {
  x <- base::as.character(x)[1]
  if (base::is.na(x) || !base::nzchar(x) || identical(x, "-") || identical(x, "---")) {
    return(list(hit_label = NA_character_, pident = NA_real_))
  }
  parts <- strsplit(x, "\\|", fixed = FALSE)[[1]]
  list(
    hit_label = parts[[1]],
    pident = if (base::length(parts) >= 2L) suppressWarnings(base::as.numeric(parts[[2]])) else NA_real_
  )
}

.dnmb_acrfinder_read_table <- function(path) {
  if (!base::file.exists(path) || !isTRUE(base::file.info(path)$size > 0)) {
    return(data.frame())
  }
  lines <- base::readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*$", lines)]
  if (!base::length(lines)) {
    return(data.frame())
  }
  utils::read.delim(text = base::paste(lines, collapse = "\n"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "")
}

.dnmb_acrfinder_parse_outputs <- function(homology_path,
                                          gba_path,
                                          identity_threshold = 30) {
  homology_tbl <- .dnmb_acrfinder_read_table(homology_path)
  gba_tbl <- .dnmb_acrfinder_read_table(gba_path)
  out <- list(
    homology = tibble::tibble(),
    gba = tibble::tibble()
  )

  if (base::nrow(homology_tbl)) {
    names(homology_tbl) <- trimws(names(homology_tbl))
    parsed_hit <- lapply(homology_tbl[["Acr_Hit|pident"]], .dnmb_acrfinder_parse_hit_string)
    homology_tbl$acr_id <- vapply(parsed_hit, `[[`, character(1), "hit_label")
    homology_tbl$pident <- vapply(parsed_hit, `[[`, numeric(1), "pident")
    homology_tbl$query <- .dnmb_module_clean_annotation_key(homology_tbl[["Protein ID"]])
    homology_tbl$evidence_mode <- "homology"
    homology_tbl <- homology_tbl[!base::is.na(homology_tbl$acr_id) & !base::is.na(homology_tbl$query), , drop = FALSE]
    if (base::nrow(homology_tbl)) {
      out$homology <- tibble::as_tibble(homology_tbl)
    }
  }

  if (base::nrow(gba_tbl)) {
    names(gba_tbl) <- trimws(names(gba_tbl))
    parsed_hit <- lapply(gba_tbl[["Acr_Hit|pident"]], .dnmb_acrfinder_parse_hit_string)
    gba_tbl$acr_id <- vapply(parsed_hit, `[[`, character(1), "hit_label")
    gba_tbl$pident <- vapply(parsed_hit, `[[`, numeric(1), "pident")
    gba_tbl$query <- .dnmb_module_clean_annotation_key(gba_tbl[["Protein ID"]])
    gba_tbl$evidence_mode <- base::tolower(gsub("[[:space:]]+", "_", gba_tbl$Classification))
    keep <- (!base::is.na(gba_tbl$acr_id) & base::nzchar(gba_tbl$acr_id)) |
      base::as.character(gba_tbl[["Acr/Aca"]]) == "Acr"
    gba_tbl <- gba_tbl[keep & !base::is.na(gba_tbl$query), , drop = FALSE]
    if (base::nrow(gba_tbl)) {
      missing_acr <- base::is.na(gba_tbl$acr_id) | !base::nzchar(gba_tbl$acr_id)
      gba_tbl$acr_id[missing_acr] <- "AcrFinder_candidate"
      out$gba <- tibble::as_tibble(gba_tbl)
    }
  }

  out
}

dnmb_acrfinder_full_normalize_hits <- function(parsed_outputs) {
  tables <- list(parsed_outputs$homology, parsed_outputs$gba)
  tables <- tables[vapply(tables, function(x) base::is.data.frame(x) && base::nrow(x), logical(1))]
  if (!base::length(tables)) {
    return(.dnmb_module_empty_optional_long_table())
  }
  tbl <- dplyr::bind_rows(tables)
  support <- vapply(seq_len(nrow(tbl)), function(i) {
    parts <- c(
      if (!base::is.na(tbl$acr_id[[i]]) && base::nzchar(tbl$acr_id[[i]])) base::paste0("acr=", tbl$acr_id[[i]]) else NA_character_,
      if (!base::is.na(tbl$pident[[i]])) base::paste0("pident=", base::round(tbl$pident[[i]], 2)) else NA_character_,
      if ("Classification" %in% names(tbl) && !base::is.na(tbl$Classification[[i]]) && base::nzchar(tbl$Classification[[i]])) base::paste0("classification=", tbl$Classification[[i]]) else NA_character_,
      if ("Self Target w/in 5000 BP" %in% names(tbl) && !base::is.na(tbl[["Self Target w/in 5000 BP"]][[i]]) && base::nzchar(tbl[["Self Target w/in 5000 BP"]][[i]])) base::paste0("self_target_near=", tbl[["Self Target w/in 5000 BP"]][[i]]) else NA_character_,
      if ("Self Target Outside 5000 BP" %in% names(tbl) && !base::is.na(tbl[["Self Target Outside 5000 BP"]][[i]]) && base::nzchar(tbl[["Self Target Outside 5000 BP"]][[i]])) base::paste0("self_target_far=", tbl[["Self Target Outside 5000 BP"]][[i]]) else NA_character_,
      if ("MGE/Prophage MetaData" %in% names(tbl) && !base::is.na(tbl[["MGE/Prophage MetaData"]][[i]]) && base::nzchar(tbl[["MGE/Prophage MetaData"]][[i]]) && tbl[["MGE/Prophage MetaData"]][[i]] != "-") base::paste0("mge=", tbl[["MGE/Prophage MetaData"]][[i]]) else NA_character_
    )
    parts <- parts[!base::is.na(parts)]
    if (!base::length(parts)) NA_character_ else base::paste(parts, collapse = "; ")
  }, character(1))
  out <- base::data.frame(
    query = .dnmb_module_clean_annotation_key(tbl$query),
    source = base::rep("acrfinder", base::nrow(tbl)),
    family_system = base::rep("AcrFinder", base::nrow(tbl)),
    family_id = base::as.character(tbl$acr_id),
    hit_label = base::as.character(tbl$acr_id),
    enzyme_role = base::rep("CRISPR-Cas", base::nrow(tbl)),
    evidence_mode = base::as.character(tbl$evidence_mode),
    substrate_label = base::rep("anti-CRISPR", base::nrow(tbl)),
    support = support,
    typing_eligible = base::rep(TRUE, base::nrow(tbl)),
    classification = if ("Classification" %in% names(tbl)) base::as.character(tbl$Classification) else NA_character_,
    pident = if ("pident" %in% names(tbl)) base::as.numeric(tbl$pident) else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- out[base::order(out$query, base::match(out$evidence_mode, c("high_confidence", "medium_confidence", "low_confidence", "homology"))), , drop = FALSE]
  out <- out[!base::duplicated(out$query), , drop = FALSE]
  base::rownames(out) <- NULL
  out[, c(.dnmb_module_optional_long_columns(), base::setdiff(base::names(out), .dnmb_module_optional_long_columns())), drop = FALSE]
}

.dnmb_acrfinder_prepare_full_input <- function(genes, output_dir, genbank) {
  if (base::is.null(genbank) || !base::nzchar(base::trimws(base::as.character(genbank)[1])) || !base::file.exists(base::as.character(genbank)[1])) {
    base::stop("AcrFinder full workflow requires a GenBank file path.", call. = FALSE)
  }
  input_dir <- base::file.path(output_dir, "subjects")
  base::dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- tools::file_path_sans_ext(base::basename(base::as.character(genbank)[1]))
  padloc_input <- .dnmb_write_padloc_input(genes, output_dir = input_dir, prefix = prefix)
  genbank_records <- .dnmb_prophage_parse_genbank_records(genbank)
  fna_path <- base::file.path(input_dir, paste0(prefix, ".fna"))
  .dnmb_prophage_write_contig_fasta(genbank_records, fna_path)
  list(
    prefix = prefix,
    fna = fna_path,
    gff = padloc_input$gff,
    faa = padloc_input$faa
  )
}

.dnmb_acrfinder_find_on_path <- function(binary, path_entries = character()) {
  binary <- base::as.character(binary)[1]
  path_entries <- base::as.character(path_entries)
  path_entries <- path_entries[!base::is.na(path_entries) & base::nzchar(path_entries)]
  candidates <- base::file.path(path_entries, binary)
  existing <- candidates[base::file.exists(candidates)]
  if (base::length(existing)) {
    return(existing[[1]])
  }
  detected <- base::Sys.which(binary)
  if (base::length(detected) && base::nzchar(detected[[1]])) {
    return(base::unname(detected[[1]]))
  }
  ""
}

.dnmb_acrfinder_find_rpsblast <- function(path_entries = character()) {
  for (binary in c("rpsblast+", "rpsblast")) {
    found <- .dnmb_acrfinder_find_on_path(binary, path_entries = path_entries)
    if (base::nzchar(found)) {
      return(found)
    }
  }
  ""
}

.dnmb_acrfinder_runtime_status <- function(path_entries = character()) {
  required_bins <- c(
    "diamond", "blastn", "perl", "clustalw", "muscle", "fuzznuc", "needle",
    "vmatch2", "mkvtree2", "vsubseqselect2"
  )
  rows <- lapply(required_bins, function(binary) {
    found <- .dnmb_acrfinder_find_on_path(binary, path_entries = path_entries)
    .dnmb_acrfinder_status_row(
      paste0("acrfinder_dependency_", binary),
      if (base::nzchar(found)) "ok" else "missing",
      if (base::nzchar(found)) found else paste0(binary, " not found in PATH.")
    )
  })
  status <- dplyr::bind_rows(rows)

  perl_path <- .dnmb_acrfinder_find_on_path("perl", path_entries = path_entries)
  if (base::nzchar(perl_path)) {
    perl_modules <- c(
      bioperl = "Bio::SeqIO",
      bioperl_run_clustalw = "Bio::Tools::Run::Alignment::Clustalw",
      date_calc = "Date::Calc",
      json_parse = "JSON::Parse"
    )
    perl_rows <- lapply(names(perl_modules), function(label) {
      module <- perl_modules[[label]]
      check <- dnmb_run_external(perl_path, args = c(paste0("-M", module), "-e", "1"), required = FALSE)
      .dnmb_acrfinder_status_row(
        paste0("acrfinder_dependency_", label),
        if (base::isTRUE(check$ok)) "ok" else "missing",
        if (base::isTRUE(check$ok)) module else paste0("Perl module ", module, " is missing.")
      )
    })
    status <- dplyr::bind_rows(status, dplyr::bind_rows(perl_rows))
  }

  rpsblast <- .dnmb_acrfinder_find_rpsblast(path_entries = path_entries)
  status <- dplyr::bind_rows(
    status,
    .dnmb_acrfinder_status_row(
      "acrfinder_dependency_rpsblast",
      if (base::nzchar(rpsblast)) "ok" else "skipped",
      if (base::nzchar(rpsblast)) rpsblast else "rpsblast/rpsblast+ not found; AcrFinder CDD filtering will be skipped."
    )
  )
  sel392 <- base::Sys.getenv("DNMB_VMATCH_SEL392", unset = "/usr/lib/vmatch/SELECT/sel392.so")
  status <- dplyr::bind_rows(
    status,
    .dnmb_acrfinder_status_row(
      "acrfinder_dependency_vmatch_sel392",
      if (base::file.exists(sel392)) "ok" else "missing",
      if (base::file.exists(sel392)) sel392 else "Vmatch selector sel392.so is missing."
    )
  )
  status
}

.dnmb_acrfinder_output_prefixes <- function(prefix) {
  prefix <- base::trimws(base::as.character(prefix)[1])
  if (base::is.na(prefix) || !base::nzchar(prefix)) {
    return(character())
  }
  unique(c(prefix, sub("\\.[0-9]+$", "", prefix)))
}

.dnmb_acrfinder_first_existing_output <- function(output_dir, prefixes, suffix) {
  candidates <- base::file.path(output_dir, base::paste0(prefixes, suffix))
  existing <- candidates[base::file.exists(candidates)]
  if (base::length(existing)) {
    return(existing[[1]])
  }
  glob <- base::list.files(output_dir, pattern = base::paste0(suffix, "$"), full.names = TRUE)
  if (base::length(glob)) {
    return(glob[[1]])
  }
  candidates[[1]]
}

.dnmb_acrfinder_output_table <- function(genes, hits) {
  genes_df <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  hits_df <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (base::nrow(hits_df) && all(c("locus_tag", "protein_id") %in% base::names(genes_df)) && "query" %in% base::names(hits_df)) {
    protein_keys <- .dnmb_module_clean_annotation_key(genes_df$protein_id)
    locus_keys <- .dnmb_module_clean_annotation_key(genes_df$locus_tag)
    hit_keys <- .dnmb_module_clean_annotation_key(hits_df$query)
    matched <- base::match(hit_keys, protein_keys)
    replace <- !base::is.na(matched) & !base::is.na(locus_keys[matched]) & base::nzchar(locus_keys[matched])
    if (base::any(replace)) {
      hits_df$query_protein_id <- hits_df$query
      hits_df$query[replace] <- locus_keys[matched[replace]]
    }
  }
  out <- .dnmb_module_output_table(genes = genes_df, hits = hits_df)
  ordered <- c(
    base::intersect(dnmb_backbone_columns(), base::names(out)),
    base::intersect(c("family_id", "hit_label", "acr_id", "pident", "qcov", "evalue", "bitscore", "align_len", "support"), base::names(out)),
    base::setdiff(base::names(out), c(base::intersect(dnmb_backbone_columns(), base::names(out)), "family_id", "hit_label", "acr_id", "pident", "qcov", "evalue", "bitscore", "align_len", "support"))
  )
  out[, ordered, drop = FALSE]
}

dnmb_run_acrfinder_module <- function(genes,
                                      output_dir,
                                      version = .dnmb_acrfinder_default_version(),
                                      cache_root = NULL,
                                      install = TRUE,
                                      repo_url = .dnmb_acrfinder_default_repo_url(),
                                      asset_urls = NULL,
                                      cpu = 1L,
                                      genbank = NULL,
                                      evalue_threshold = 1e-5,
                                      identity_threshold = 30,
                                      coverage_threshold = 0.8,
                                      genome_type = "AUTO",
                                      min_cdd_proteins = 2L,
                                      use_prophage_db = TRUE,
                                      use_gi_db = TRUE) {
  base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  status <- .dnmb_acrfinder_empty_status()
  trace_log <- base::file.path(output_dir, "acrfinder_module_trace.log")
  resolved_genome_type <- .dnmb_acrfinder_resolve_genome_type(genome_type = genome_type, genbank = genbank)
  status <- dplyr::bind_rows(
    status,
    .dnmb_acrfinder_status_row("acrfinder_genome_type", "ok", resolved_genome_type)
  )

  install_result <- dnmb_acrfinder_install_module(
    version = version,
    cache_root = cache_root,
    install = install,
    repo_url = repo_url,
    asset_urls = asset_urls
  )
  status <- dplyr::bind_rows(status, install_result$status)
  if (!base::isTRUE(install_result$ok)) {
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), raw_hits = .dnmb_acrfinder_empty_hits(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  module <- dnmb_acrfinder_get_module(version = version, cache_root = cache_root, required = TRUE)
  full_input <- tryCatch(
    .dnmb_acrfinder_prepare_full_input(genes = genes, output_dir = output_dir, genbank = genbank),
    error = function(e) e
  )
  if (inherits(full_input, "error")) {
    status <- dplyr::bind_rows(status, .dnmb_acrfinder_status_row("acrfinder_input", "failed", conditionMessage(full_input)))
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), raw_hits = .dnmb_acrfinder_empty_hits(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }

  final_dir <- base::file.path(output_dir, "acrfinder_output")
  if (base::dir.exists(final_dir)) {
    base::unlink(final_dir, recursive = TRUE, force = TRUE)
  }
  base::dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  cmd <- c(
    base::file.path(module$manifest$repo_dir, "acr_aca_cri_runner.py"),
    "-n", full_input$fna,
    "-f", full_input$gff,
    "-a", full_input$faa,
    "-o", final_dir,
    "-z", resolved_genome_type,
    "-c", base::as.character(base::as.integer(min_cdd_proteins)[1]),
    "-p", if (base::isTRUE(use_prophage_db)) "true" else "false",
    "-g", if (base::isTRUE(use_gi_db)) "true" else "false",
    "--identity", base::as.character(identity_threshold),
    "--coverage", base::as.character(coverage_threshold),
    "--e_value", base::as.character(evalue_threshold),
    "--num_threads", base::as.character(base::as.integer(cpu)[1])
  )
  base::cat(base::paste0("[", base::Sys.time(), "] ", .dnmb_format_command(module$files$env_python, cmd), "\n"), file = trace_log, append = TRUE)
  # Order priority for AcrFinder's tool resolution:
  #   1. AcrFinder's own bin
  #   2. AcrFinder python venv bin
  #   3. Distro system bins (/usr/bin, /usr/local/bin, /usr/sbin)
  #      explicitly hoisted ahead of the inherited PATH so an
  #      `/opt/biotools/bin` conda env never wins over the apt-installed
  #      perl that has Bio::SeqIO / Date::Calc / JSON::Parse on its @INC,
  #      and the apt clustalw / muscle / vmatch links resolve cleanly.
  #   4. The rest of inherited PATH (covers conda envs, /opt/biotools,
  #      etc. for tools that genuinely live there).
  #   5. CRISPRCasFinder bundled bin as last-resort fallback for bare-
  #      metal users without system clustalw2/vmatch2/muscle binaries.
  system_first <- c("/usr/bin", "/usr/local/bin", "/usr/sbin")
  inherited <- strsplit(base::Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  acr_path_entries <- c(
    base::file.path(module$manifest$repo_dir, "bin"),
    base::dirname(module$files$env_python),
    system_first,
    inherited,
    base::file.path(module$manifest$repo_dir, "dependencies", "CRISPRCasFinder", "bin")
  )
  acr_path_entries <- unique(acr_path_entries[base::nzchar(acr_path_entries) & base::file.exists(acr_path_entries)])
  dependency_status <- .dnmb_acrfinder_runtime_status(acr_path_entries)
  status <- dplyr::bind_rows(status, dependency_status)
  blocking_dependency <- dependency_status$status %in% c("missing", "failed")
  if (base::any(blocking_dependency)) {
    detail <- base::paste(dependency_status$detail[blocking_dependency], collapse = "; ")
    status <- dplyr::bind_rows(status, .dnmb_acrfinder_status_row("acrfinder_run", "failed", detail))
    return(list(ok = FALSE, status = status, files = list(trace_log = trace_log), raw_hits = .dnmb_acrfinder_empty_hits(), hits = .dnmb_module_empty_optional_long_table(), output_table = data.frame()))
  }
  run <- dnmb_run_external(
    module$files$env_python,
    args = cmd,
    wd = module$manifest$repo_dir,
    env = c(
      PATH = base::paste(acr_path_entries, collapse = .Platform$path.sep),
      MACSY_HOME = base::file.path(module$manifest$repo_dir, "dependencies", "CRISPRCasFinder", "macsyfinder-1.0.5")
    ),
    required = FALSE
  )
  .dnmb_acrfinder_trace_result(trace_log, run, label = "acrfinder_run")

  output_prefixes <- .dnmb_acrfinder_output_prefixes(full_input$prefix)
  homology_path <- .dnmb_acrfinder_first_existing_output(final_dir, output_prefixes, "_homology_based.out")
  gba_candidates <- c(
    vapply(
      output_prefixes,
      function(prefix) .dnmb_acrfinder_first_existing_output(final_dir, prefix, "_guilt_by_association.out"),
      character(1),
      USE.NAMES = FALSE
    ),
    vapply(
      output_prefixes,
      function(prefix) .dnmb_acrfinder_first_existing_output(final_dir, prefix, "_guilt-by-association.out"),
      character(1),
      USE.NAMES = FALSE
    )
  )
  gba_existing <- gba_candidates[base::file.exists(gba_candidates)]
  gba_path <- if (base::length(gba_existing)) gba_existing[[1]] else gba_candidates[[1]]
  parsed_outputs <- .dnmb_acrfinder_parse_outputs(
    homology_path = homology_path,
    gba_path = gba_path,
    identity_threshold = identity_threshold
  )
  hits <- dnmb_acrfinder_full_normalize_hits(parsed_outputs)
  output_table <- .dnmb_acrfinder_output_table(genes = genes, hits = hits)

  raw_hits_n <- 0L
  if (base::is.data.frame(parsed_outputs$homology)) raw_hits_n <- raw_hits_n + base::nrow(parsed_outputs$homology)
  if (base::is.data.frame(parsed_outputs$gba)) raw_hits_n <- raw_hits_n + base::nrow(parsed_outputs$gba)
  effective_ok <- base::isTRUE(run$ok) || raw_hits_n > 0L
  status <- dplyr::bind_rows(
    status,
    .dnmb_acrfinder_status_row(
      "acrfinder_run",
      if (base::isTRUE(run$ok)) {
        if (raw_hits_n > 0L) "ok" else "empty"
      } else if (raw_hits_n > 0L) {
        "partial"
      } else if (!base::nzchar(run$resolved_command)) {
        "missing"
      } else {
        "failed"
      },
      if (base::isTRUE(run$ok)) final_dir else (run$error %||% final_dir)
    )
  )

  list(
    ok = effective_ok,
    status = status,
    files = list(
      trace_log = trace_log,
      input_fna = full_input$fna,
      input_gff = full_input$gff,
      input_faa = full_input$faa,
      homology = homology_path,
      gba = gba_path
    ),
    raw_hits = dplyr::bind_rows(parsed_outputs$homology, parsed_outputs$gba),
    hits = hits,
    output_table = output_table
  )
}

#' Run InterProScan on protein sequences
#'
#' Executes InterProScan via \code{interproscan.sh} and returns the path to the
#' output TSV.
#' This function is designed for Linux environments (including Docker) where
#' InterProScan is installed.
#'
#' @param fasta_path Path to a FASTA file of protein sequences.
#'   When \code{NULL}, proteins are extracted from the active
#'   \code{genbank_table}.
#' @param output_dir Directory for InterProScan output files.
#' @param applications Character vector of InterProScan analyses to run (e.g.,
#'   \code{c("Pfam", "TIGRFAM", "CDD")}). \code{NULL} runs all available.
#' @param cpu Number of threads. Defaults to 80\% of available cores.
#' @param extra_args Additional arguments passed to \code{interproscan.sh}.
#' @param interproscan_path Path to the \code{interproscan.sh} script. When
#'   \code{NULL}, searches \code{INTERPROSCAN_HOME} env var, then PATH.
#' @param genbank_table Optional gene table used when \code{fasta_path} is
#'   \code{NULL}.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list with components:
#'   \describe{
#'     \item{tsv}{Path to the main \code{.tsv} output file.}
#'     \item{tsv_sites}{Path to the \code{.tsv.sites} file, or \code{NULL}.}
#'     \item{exit_code}{Integer exit code of the InterProScan process.}
#'   }
#' @export
run_interproscan <- function(fasta_path = NULL,
                             output_dir = NULL,
                             applications = NULL,
                             cpu = .dnmb_default_cpu(),
                             extra_args = character(),
                             interproscan_path = NULL,
                             genbank_table = NULL,
                             verbose = TRUE) {

  # ── Locate interproscan.sh ──────────────────────────────────
  ipr_bin <- .dnmb_find_interproscan(interproscan_path)

  # ── Prepare output directory ────────────────────────────────
  if (is.null(output_dir) || !nzchar(trimws(output_dir))) {
    output_dir <- file.path(getwd(), "dnmb_interproscan")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  # ── Prepare input FASTA ─────────────────────────────────────
  if (is.null(fasta_path)) {
    genes <- .dnmb_resolve_run_genes(genbank_table = genbank_table)
    fasta_path <- file.path(output_dir, "query_proteins.faa")
    prep <- .dnmb_write_query_fasta(genes, fasta_path)
    if (prep$n == 0L) {
      stop("No protein sequences available for InterProScan.", call. = FALSE)
    }
    if (isTRUE(verbose)) {
      message("[InterProScan] Wrote ", prep$n, " protein sequences to ", fasta_path)
    }
  }
  fasta_path <- normalizePath(fasta_path, winslash = "/", mustWork = TRUE)

  # ── Build command arguments ─────────────────────────────────
  output_base <- file.path(output_dir, "interproscan_results")
  args <- c(
    "-i", fasta_path,
    "-o", paste0(output_base, ".tsv"),
    "-f", "tsv",
    "-cpu", as.character(as.integer(cpu)[1]),
    "--goterms",
    "--pathways",
    "--iprlookup",
    "-dra"
  )

  if (!is.null(applications) && length(applications) > 0L) {
    args <- c(args, "-appl", paste(applications, collapse = ","))
  }

  if (length(extra_args) > 0L) {
    args <- c(args, as.character(extra_args))
  }

  # ── Run InterProScan ────────────────────────────────────────
  if (isTRUE(verbose)) {
    message("[InterProScan] Running: ", ipr_bin, " ", paste(args, collapse = " "))
  }

  exit_code <- system2(ipr_bin, args = args, stdout = "", stderr = "")

  tsv_path <- paste0(output_base, ".tsv")
  tsv_sites_path <- paste0(output_base, ".tsv.sites")

  if (exit_code != 0L) {
    warning("InterProScan exited with code ", exit_code, call. = FALSE)
  }

  if (!file.exists(tsv_path)) {
    stop("InterProScan output file not found: ", tsv_path, call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message("[InterProScan] Done. Output: ", tsv_path)
  }

  list(
    tsv = tsv_path,
    tsv_sites = if (file.exists(tsv_sites_path)) tsv_sites_path else NULL,
    exit_code = exit_code
  )
}


#' Run InterProScan and integrate results into genbank_table
#'
#' Convenience wrapper that runs InterProScan, parses the output via
#' \code{InterProScan_annotations()}, and merges the result onto
#' \code{genbank_table} in the global environment.
#'
#' @inheritParams run_interproscan
#' @return Invisibly returns the updated \code{genbank_table}.
#' @export
run_interproscan_pipeline <- function(fasta_path = NULL,
                                      output_dir = NULL,
                                      applications = NULL,
                                      cpu = .dnmb_default_cpu(),
                                      extra_args = character(),
                                      interproscan_path = NULL,
                                      genbank_table = NULL,
                                      verbose = TRUE) {

  result <- run_interproscan(
    fasta_path = fasta_path,
    output_dir = output_dir,
    applications = applications,
    cpu = cpu,
    extra_args = extra_args,
    interproscan_path = interproscan_path,
    genbank_table = genbank_table,
    verbose = verbose
  )

  # Parse with existing organizer
  InterProScan_annotations(InterProScan_dir = dirname(result$tsv))

  # Merge onto genbank_table
  ipr_table <- get0("InterProScan_table", envir = .GlobalEnv, inherits = FALSE)
  if (!is.null(ipr_table)) {
    gt <- get0("genbank_table", envir = .GlobalEnv, inherits = FALSE)
    if (!is.null(gt) && "locus_tag" %in% colnames(gt) && "query" %in% colnames(ipr_table)) {
      gt <- merge(gt, ipr_table, by.x = "locus_tag", by.y = "query", all.x = TRUE)
      assign("genbank_table", gt, envir = .GlobalEnv)
      if (isTRUE(verbose)) {
        message("[InterProScan] Merged InterProScan results into genbank_table.")
      }
    }
  }

  invisible(get0("genbank_table", envir = .GlobalEnv, inherits = FALSE))
}


# ── Internal helpers ──────────────────────────────────────────

.dnmb_interproscan_fallback_version <- function() "5.72-103.0"

.dnmb_interproscan_normalize_version <- function(version) {
  value <- trimws(as.character(version)[1])
  if (is.na(value) || !nzchar(value)) {
    return("")
  }
  sub("^v", "", value)
}

.dnmb_interproscan_latest_version <- function(default = .dnmb_interproscan_fallback_version()) {
  override <- .dnmb_interproscan_normalize_version(Sys.getenv("DNMB_INTERPROSCAN_VERSION", unset = ""))
  if (nzchar(override)) {
    return(override)
  }

  remote_info <- tryCatch(
    .dnmb_db_remote_check_interproscan(list(version = .dnmb_interproscan_normalize_version(default))),
    error = function(e) NULL
  )
  remote_version <- .dnmb_interproscan_normalize_version(
    if (!is.null(remote_info)) remote_info$remote_version else ""
  )

  if (nzchar(remote_version)) {
    remote_version
  } else {
    .dnmb_interproscan_normalize_version(default)
  }
}

.dnmb_interproscan_installed_script <- function(cache_root = NULL) {
  ipr_root <- file.path(.dnmb_db_cache_root(cache_root = cache_root, create = FALSE), "interproscan")
  if (!dir.exists(ipr_root)) {
    return("")
  }

  version_dirs <- list.dirs(ipr_root, full.names = TRUE, recursive = FALSE)
  if (!length(version_dirs)) {
    return("")
  }

  scripts <- file.path(version_dirs, "interproscan.sh")
  scripts <- scripts[file.exists(scripts)]
  if (!length(scripts)) {
    return("")
  }

  info <- file.info(scripts)
  scripts <- scripts[order(info$mtime, decreasing = TRUE, na.last = TRUE)]
  normalizePath(scripts[[1]], winslash = "/", mustWork = TRUE)
}

.dnmb_interproscan_installed_version <- function(cache_root = NULL) {
  script_path <- .dnmb_interproscan_installed_script(cache_root = cache_root)
  if (!nzchar(script_path)) {
    return("")
  }
  basename(dirname(script_path))
}

.dnmb_interproscan_default_version <- function(cache_root = NULL) {
  .dnmb_interproscan_latest_version(default = .dnmb_interproscan_fallback_version())
}

.dnmb_ensure_interproscan <- function(cache_root = NULL, version = NULL) {
  version <- .dnmb_interproscan_normalize_version(
    if (is.null(version) || !nzchar(trimws(as.character(version)[1]))) {
      .dnmb_interproscan_default_version(cache_root = cache_root)
    } else {
      version
    }
  )
  ipr_dir <- file.path(.dnmb_db_cache_root(cache_root = cache_root, create = TRUE), "interproscan", version)
  ipr_sh <- file.path(ipr_dir, "interproscan.sh")
  if (file.exists(ipr_sh)) {
    if (is.null(dnmb_db_read_manifest("interproscan", version, cache_root = cache_root, required = FALSE))) {
      try(
        dnmb_db_write_manifest(
          "interproscan",
          version,
          manifest = list(source_url = NA_character_, tarball = NA_character_),
          cache_root = cache_root,
          overwrite = FALSE
        ),
        silent = TRUE
      )
    }
    return(ipr_dir)
  }

  message("[InterProScan] Downloading InterProScan ", version, " to cache...")
  message("[InterProScan] This is a one-time download (~15GB). Please wait.")
  dir.create(ipr_dir, recursive = TRUE, showWarnings = FALSE)

  tarball <- paste0("interproscan-", version, "-64-bit.tar.gz")
  url <- paste0("https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/", version, "/", tarball)
  dest <- file.path(tempdir(), tarball)

  dl <- dnmb_run_external("wget", args = c("-q", "--no-check-certificate", url, "-O", dest), required = FALSE)
  if (!isTRUE(dl$ok) || !file.exists(dest)) {
    stop("[InterProScan] Download failed. Check internet connection.", call. = FALSE)
  }

  message("[InterProScan] Extracting...")
  ex <- dnmb_run_external("tar", args = c("-xzf", dest, "--strip-components=1", "-C", ipr_dir), required = FALSE)
  unlink(dest, force = TRUE)
  if (!isTRUE(ex$ok)) {
    stop("[InterProScan] Extraction failed.", call. = FALSE)
  }

  # Setup
  setup_py <- file.path(ipr_dir, "setup.py")
  if (file.exists(setup_py)) {
    message("[InterProScan] Running setup...")
    dnmb_run_external("python3", args = c(setup_py, "-f", file.path(ipr_dir, "interproscan.properties")), wd = ipr_dir, required = FALSE)
  }

  if (file.exists(ipr_sh)) {
    try(
      dnmb_db_write_manifest(
        "interproscan",
        version,
        manifest = list(source_url = url, tarball = tarball),
        cache_root = cache_root,
        overwrite = TRUE
      ),
      silent = TRUE
    )
    message("[InterProScan] Installation complete: ", ipr_dir)
    return(ipr_dir)
  }
  stop("[InterProScan] Installation failed.", call. = FALSE)
}

.dnmb_find_interproscan <- function(path = NULL, cache_root = NULL) {
  # 1. User-supplied path
  if (!is.null(path) && nzchar(trimws(path))) {
    path <- path.expand(trimws(path))
    if (file.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
    stop("interproscan.sh not found at: ", path, call. = FALSE)
  }

  # 2. INTERPROSCAN_HOME environment variable
  home <- Sys.getenv("INTERPROSCAN_HOME", "")
  if (nzchar(home)) {
    candidate <- file.path(home, "interproscan.sh")
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  # 3. PATH lookup
  on_path <- Sys.which("interproscan.sh")
  if (nzchar(on_path)) {
    return(normalizePath(on_path, winslash = "/", mustWork = TRUE))
  }

  # 4. Auto-download to cache
  cached_candidate <- .dnmb_interproscan_installed_script(cache_root = cache_root)
  ipr_dir <- tryCatch(
    .dnmb_ensure_interproscan(cache_root = cache_root),
    error = function(e) {
      if (nzchar(cached_candidate)) {
        warning(
          "[InterProScan] Could not prepare the latest cached release (",
          conditionMessage(e),
          "). Falling back to existing cached installation at ",
          cached_candidate,
          call. = FALSE
        )
      }
      NULL
    }
  )
  if (!is.null(ipr_dir)) {
    candidate <- file.path(ipr_dir, "interproscan.sh")
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  if (nzchar(cached_candidate)) {
    return(cached_candidate)
  }

  stop(
    "interproscan.sh not found. ",
    "Set the INTERPROSCAN_HOME environment variable or pass the path explicitly. ",
    "InterProScan is only available on Linux; use the DNMB Docker image for macOS/Windows.",
    call. = FALSE
  )
}

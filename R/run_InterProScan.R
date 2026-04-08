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

  # Copy output to working directory so InterProScan_annotations() finds it
  tsv_dest <- file.path(getwd(), basename(result$tsv))
  file.copy(result$tsv, tsv_dest, overwrite = TRUE)
  if (!is.null(result$tsv_sites)) {
    file.copy(result$tsv_sites, file.path(getwd(), basename(result$tsv_sites)), overwrite = TRUE)
  }

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

.dnmb_interproscan_pinned_version <- function() "5.77-108.0"

.dnmb_interproscan_latest_version <- function() {
  latest <- tryCatch({
    con <- url("https://api.github.com/repos/ebi-pf-team/interproscan/releases/latest")
    on.exit(close(con))
    lines <- readLines(con, warn = FALSE)
    tag <- sub('.*"tag_name"\\s*:\\s*"([^"]+)".*', "\\1", grep('"tag_name"', lines, value = TRUE)[1])
    if (is.na(tag) || !nzchar(tag)) NA_character_ else trimws(tag)
  }, error = function(e) NA_character_)

  if (!is.na(latest) && nzchar(latest)) {
    latest
  } else {
    .dnmb_interproscan_pinned_version()
  }
}

.dnmb_interproscan_default_version <- function() .dnmb_interproscan_latest_version()

.dnmb_interproscan_cache_root <- function(cache_root = NULL, create = FALSE) {
  root <- file.path(.dnmb_db_cache_root(cache_root = cache_root, create = create), "interproscan")
  if (isTRUE(create)) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
  }
  if (dir.exists(root)) {
    normalizePath(root, winslash = "/", mustWork = FALSE)
  } else {
    root
  }
}

.dnmb_interproscan_cached_versions <- function(cache_root = NULL) {
  root <- .dnmb_interproscan_cache_root(cache_root = cache_root, create = FALSE)
  if (!dir.exists(root)) {
    return(character())
  }
  dirs <- list.dirs(root, full.names = FALSE, recursive = FALSE)
  dirs <- dirs[nzchar(dirs)]
  if (!length(dirs)) {
    return(character())
  }
  ok <- vapply(dirs, function(x) {
    !inherits(try(base::package_version(gsub("-", ".", x, fixed = TRUE)), silent = TRUE), "try-error")
  }, logical(1))
  dirs <- dirs[ok]
  if (!length(dirs)) {
    return(character())
  }
  order_key <- vapply(dirs, function(x) as.character(base::package_version(gsub("-", ".", x, fixed = TRUE))), character(1))
  dirs[order(base::package_version(order_key), decreasing = TRUE)]
}

.dnmb_interproscan_existing_path <- function(cache_root = NULL, version = NULL) {
  if (!is.null(version) && nzchar(trimws(version))) {
    candidate <- file.path(.dnmb_interproscan_cache_root(cache_root = cache_root, create = FALSE), trimws(version), "interproscan.sh")
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  versions <- .dnmb_interproscan_cached_versions(cache_root = cache_root)
  for (ver in versions) {
    candidate <- file.path(.dnmb_interproscan_cache_root(cache_root = cache_root, create = FALSE), ver, "interproscan.sh")
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  NULL
}

.dnmb_ensure_interproscan <- function(cache_root = NULL, version = .dnmb_interproscan_default_version()) {
  ipr_dir <- file.path(.dnmb_interproscan_cache_root(cache_root = cache_root, create = TRUE), version)
  ipr_sh <- file.path(ipr_dir, "interproscan.sh")
  if (file.exists(ipr_sh)) {
    .dnmb_db_autoprune_default_versions(
      module = "interproscan",
      version = version,
      default_version = .dnmb_interproscan_default_version(),
      cache_root = cache_root,
      preserve = character()
    )
    return(ipr_dir)
  }

  # Incomplete previous download/extraction can leave an empty or partial
  # version directory behind. Remove it before retrying the install.
  if (dir.exists(ipr_dir) && !file.exists(ipr_sh)) {
    unlink(ipr_dir, recursive = TRUE, force = TRUE)
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
    .dnmb_db_autoprune_default_versions(
      module = "interproscan",
      version = version,
      default_version = .dnmb_interproscan_default_version(),
      cache_root = cache_root,
      preserve = character()
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

  # 4. Existing cached install
  cached <- .dnmb_interproscan_existing_path(
    cache_root = cache_root,
    version = .dnmb_interproscan_default_version()
  )
  if (!is.null(cached) && file.exists(cached)) {
    return(cached)
  }

  # 5. Auto-download to cache
  ensure_error <- NULL
  ipr_dir <- tryCatch(
    .dnmb_ensure_interproscan(cache_root = cache_root),
    error = function(e) {
      ensure_error <<- conditionMessage(e)
      NULL
    }
  )
  if (!is.null(ipr_dir)) {
    candidate <- file.path(ipr_dir, "interproscan.sh")
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  if (!is.null(ensure_error) && nzchar(ensure_error)) {
    stop(
      "InterProScan setup failed: ",
      ensure_error,
      call. = FALSE
    )
  } else {
    stop(
      "interproscan.sh not found. ",
      "Set the INTERPROSCAN_HOME environment variable or pass the path explicitly. ",
      "InterProScan is only available on Linux; use the DNMB Docker image for macOS/Windows.",
      call. = FALSE
    )
  }
}

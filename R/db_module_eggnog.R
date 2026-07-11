#' EggNOG-mapper module for DNMB
#'
#' Runs eggnog-mapper (emapper.py) on protein sequences and returns
#' standardized module output with KO, COG, EC, and other functional
#' annotations.
#'
#' @name dnmb_eggnog_module
NULL

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

.dnmb_eggnog_module_name <- function() "EggNOG"

.dnmb_eggnog_default_db <- function() "auto"

.dnmb_eggnog_detect_tax_scope <- function(genbank = NULL) {
  if (is.null(genbank) || !file.exists(genbank)) return(NULL)
  lines <- tryCatch(readLines(genbank, n = 30L, warn = FALSE), error = function(e) character())
  org_line <- grep("^\\s+(Bacteria|Archaea|Eukaryota)", lines, value = TRUE)
  if (!length(org_line)) return(NULL)
  domain <- trimws(sub(";.*", "", trimws(org_line[[1]])))
  scope <- switch(domain,
    Bacteria = "Bacteria",
    Archaea = "Archaea",
    Eukaryota = "Eukaryota",
    NULL
  )
  scope
}

.dnmb_eggnog_tsv_columns <- function() {

  c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
    "max_annot_lvl", "COG_category", "Description", "Preferred_name",
    "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
    "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
    "CAZy", "BiGG_Reaction", "PFAMs")
}

.dnmb_eggnog_cog_colors <- function() {
  c(
    J = "#ff0000", A = "#c2af58", K = "#ff9900", L = "#ffff00",
    B = "#ffc600", D = "#99ff00", Y = "#493126", V = "#ff008a",
    T = "#0000ff", M = "#9ec928", N = "#006633", Z = "#660099",
    W = "#336699", U = "#33cc99", O = "#00ffff", C = "#9900ff",
    G = "#805642", E = "#ff00ff", F = "#99334d", H = "#727dcc",
    I = "#5c5a1b", P = "#0099ff", Q = "#ffcc99", R = "#ff9999",
    S = "#d6aadf"
  )
}

.dnmb_eggnog_cog_legends <- function() {
  c(
    J = "[J] Translation, ribosomal structure and biogenesis",
    A = "[A] RNA processing and modification",
    K = "[K] Transcription",
    L = "[L] Replication, recombination and repair",
    B = "[B] Chromatin structure and dynamics",
    D = "[D] Cell cycle control, cell division, chromosome partitioning",
    Y = "[Y] Nuclear structure",
    V = "[V] Defense mechanisms",
    T = "[T] Signal transduction mechanisms",
    M = "[M] Cell wall/membrane/envelope biogenesis",
    N = "[N] Cell motility",
    Z = "[Z] Cytoskeleton",
    W = "[W] Extracellular structures",
    U = "[U] Intracellular trafficking, secretion, and vesicular transport",
    O = "[O] Posttranslational modification, protein turnover, chaperones",
    C = "[C] Energy production and conversion",
    G = "[G] Carbohydrate transport and metabolism",
    E = "[E] Amino acid transport and metabolism",
    F = "[F] Nucleotide transport and metabolism",
    H = "[H] Coenzyme transport and metabolism",
    I = "[I] Lipid transport and metabolism",
    P = "[P] Inorganic ion transport and metabolism",
    Q = "[Q] Secondary metabolites biosynthesis, transport and catabolism",
    R = "[R] General function prediction only",
    S = "[S] Function unknown"
  )
}

# ---------------------------------------------------------------------------
# Detection & installation
# ---------------------------------------------------------------------------

.dnmb_eggnog_required_version <- function() {
  "2.1.15"
}

.dnmb_eggnog_database_release <- function() {
  "5.0.2"
}

.dnmb_eggnog_version_is_supported <- function(version) {
  base::identical(base::as.character(version)[1], .dnmb_eggnog_required_version())
}

.dnmb_eggnog_conda_package_version <- function(emapper_path) {
  emapper_path <- base::normalizePath(emapper_path, winslash = "/", mustWork = FALSE)
  metadata_dir <- base::file.path(base::dirname(base::dirname(emapper_path)), "conda-meta")
  metadata_files <- base::list.files(
    metadata_dir,
    pattern = "^eggnog-mapper-.*\\.json$",
    full.names = TRUE
  )
  if (!base::length(metadata_files)) {
    return(NA_character_)
  }
  for (metadata_path in base::rev(base::sort(metadata_files))) {
    metadata <- tryCatch(
      jsonlite::fromJSON(metadata_path, simplifyVector = TRUE),
      error = function(e) NULL
    )
    if (!base::is.null(metadata) && base::identical(metadata$name, "eggnog-mapper") &&
        !base::is.null(metadata$version) && base::nzchar(base::as.character(metadata$version)[1])) {
      return(base::as.character(metadata$version)[1])
    }
  }
  NA_character_
}

#' Detect eggnog-mapper (emapper.py)
#'
#' Searches PATH and common conda environments for emapper.py. Returns a list
#' with `found`, `path`, and `version`.
#'
#' @param required If TRUE, stops on failure.
#' @return A list with components `found`, `path`, `version`.
#' @export
dnmb_detect_emapper <- function(required = FALSE) {
  result <- list(
    found = FALSE,
    path = NA_character_,
    version = NA_character_,
    cli_version = NA_character_,
    package_version = NA_character_
  )

  # Try emapper.py on current PATH
  detection <- dnmb_detect_binary("emapper.py", required = FALSE)

  # If not found, search common conda env paths
  if (!isTRUE(detection$found)) {
    conda_roots <- c(
      file.path(Sys.getenv("HOME"), "miniforge3"),
      file.path(Sys.getenv("HOME"), "miniconda3"),
      file.path(Sys.getenv("HOME"), "anaconda3"),
      "/opt/conda",
      "/opt/biotools"
    )
    for (root in conda_roots) {
      env_dirs <- list.dirs(file.path(root, "envs"), recursive = FALSE, full.names = TRUE)
      candidates <- c(
        file.path(root, "bin", "emapper.py"),
        file.path(env_dirs, "bin", "emapper.py")
      )
      found_path <- candidates[file.exists(candidates)]
      if (length(found_path)) {
        detection <- list(found = TRUE, path = found_path[1])
        break
      }
    }
  }

  if (!isTRUE(detection$found)) {
    if (isTRUE(required)) {
      stop(
        "emapper.py not found. Install eggnog-mapper: ",
        "conda install -c bioconda eggnog-mapper",
        call. = FALSE
      )
    }
    return(result)
  }

  result$found <- TRUE
  result$path <- detection$path

  # Get version (use full path to avoid PATH dependency)
  ver_run <- dnmb_run_external(result$path, args = "--version", required = FALSE)
  ver_lines <- .dnmb_nonempty_lines(c(ver_run$stdout, ver_run$stderr))
  ver_match <- grep("emapper-([0-9.]+)", ver_lines, value = TRUE)
  if (length(ver_match)) {
    m <- regmatches(ver_match[1], regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", ver_match[1]))
    if (length(m)) result$version <- m
  }
  result$cli_version <- result$version
  result$package_version <- .dnmb_eggnog_conda_package_version(result$path)
  if (!base::is.na(result$package_version)) {
    # Conda metadata is the authoritative release identity for packaged builds.
    result$version <- result$package_version
  }

  result
}

.dnmb_eggnog_database_release_from_db <- function(data_dir) {
  annotation_db <- base::file.path(data_dir, "eggnog.db")
  if (!.dnmb_nonempty_file(annotation_db)) return(NA_character_)

  query <- "SELECT version FROM version LIMIT 1;"
  sqlite <- dnmb_detect_binary("sqlite3", required = FALSE)
  if (base::isTRUE(sqlite$found)) {
    run <- dnmb_run_external(sqlite$path, args = c(annotation_db, query), required = FALSE)
    value <- .dnmb_nonempty_lines(run$stdout)
    if (base::isTRUE(run$ok) && base::length(value)) return(base::trimws(value[[1]]))
  }

  python_query <- paste0(
    "import sqlite3,sys; ",
    "row=sqlite3.connect(sys.argv[1]).execute('SELECT version FROM version LIMIT 1').fetchone(); ",
    "print(row[0] if row else '')"
  )
  python_paths <- unique(c(
    dnmb_detect_binary("python3", required = FALSE)$path,
    dnmb_detect_binary("python", required = FALSE)$path
  ))
  python_paths <- python_paths[!base::is.na(python_paths) & base::nzchar(python_paths)]
  for (python_path in python_paths) {
    run <- dnmb_run_external(
      python_path,
      args = c("-c", python_query, annotation_db),
      required = FALSE
    )
    value <- .dnmb_nonempty_lines(run$stdout)
    if (base::isTRUE(run$ok) && base::length(value)) return(base::trimws(value[[1]]))
  }
  NA_character_
}

#' Check or download eggnog-mapper databases
#'
#' Verifies the eggnog-mapper database directory exists; optionally downloads
#' databases via `download_eggnog_data.py`.
#'
#' @param data_dir Database directory. If NULL, uses the emapper default
#'   (`~/.local/share/eggnog-mapper/`).
#' @param db_name Database name (default: "bact").
#' @param install Whether to download missing databases.
#' @param cache_root Optional DNMB cache root.
#' @return A list with `ok`, `data_dir`, `db_name`.
#' @export
dnmb_eggnog_ensure_database <- function(data_dir = NULL,
                                        db_name = .dnmb_eggnog_default_db(),
                                        install = TRUE,
                                        cache_root = NULL) {
  if (is.null(data_dir)) {
    # Check EGGNOG_DATA_DIR env var first
    env_dir <- Sys.getenv("EGGNOG_DATA_DIR", unset = "")
    if (nzchar(env_dir) && dir.exists(env_dir)) {
      data_dir <- env_dir
    } else {
      # Use DNMB cache root → db_modules/eggnog/ (same level as other modules)
      db_root <- .dnmb_db_cache_root(cache_root = cache_root, create = TRUE)
      data_dir <- file.path(db_root, "eggnog", "data")
      dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  # Check if database files exist
  diamond_db <- file.path(data_dir, "eggnog_proteins.dmnd")
  taxa_db <- file.path(data_dir, "eggnog.taxa.db")
  annotation_db <- file.path(data_dir, "eggnog.db")
  required_db_files <- c(diamond_db, annotation_db, taxa_db)
  db_ok <- base::all(base::vapply(required_db_files, .dnmb_nonempty_file, logical(1)))

  if (!db_ok && isTRUE(install)) {
    message("[DNMB EggNOG] Downloading eggnog-mapper databases to: ", data_dir)
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

    # Try download_eggnog_data.py first
    dl_cmd <- "download_eggnog_data.py"
    emapper <- dnmb_detect_emapper(required = FALSE)
    if (isTRUE(emapper$found) && !is.na(emapper$path)) {
      dl_candidate <- file.path(dirname(emapper$path), "download_eggnog_data.py")
      if (file.exists(dl_candidate)) dl_cmd <- dl_candidate
    }
    dl_run <- dnmb_run_external(dl_cmd, args = c("-y", "--data_dir", data_dir), required = FALSE)

    # Re-check after download_eggnog_data.py
    db_ok <- base::all(base::vapply(required_db_files, .dnmb_nonempty_file, logical(1)))

    # Fallback: direct download (eggnogdb.embl.de sometimes down, use eggnog5 mirror)
    if (!db_ok) {
      message("[DNMB EggNOG] Official downloader failed, trying mirror (eggnog5.embl.de)...")
      db_ok <- .dnmb_eggnog_fallback_download(data_dir)
    }

    if (!db_ok) {
      return(list(ok = FALSE, data_dir = data_dir, db_name = db_name,
                  error = "Database download failed. Run manually: download_eggnog_data.py -y --data_dir <path>"))
    }
  }

  manifest <- NULL
  if (base::isTRUE(db_ok)) {
    emapper <- dnmb_detect_emapper(required = FALSE)
    db_files <- c(
      diamond_db,
      taxa_db,
      annotation_db
    )
    db_files <- db_files[base::file.exists(db_files)]
    manifest <- list(
      install_ok = TRUE,
      data_dir = base::normalizePath(data_dir, winslash = "/", mustWork = FALSE),
      tool_version = emapper$version %||% NA_character_,
      emapper_version = emapper$version %||% NA_character_,
      expected_tool_version = .dnmb_eggnog_required_version(),
      expected_db_release = .dnmb_eggnog_database_release(),
      db_release = .dnmb_eggnog_database_release_from_db(data_dir),
      asset_state = .dnmb_db_local_asset_state(db_files)
    )
    dnmb_db_write_manifest(
      "eggnog", "data",
      manifest = manifest,
      cache_root = cache_root,
      overwrite = TRUE
    )
    manifest <- dnmb_db_read_manifest("eggnog", "data", cache_root = cache_root, required = TRUE)
  }

  list(ok = db_ok, data_dir = data_dir, db_name = db_name, manifest = manifest)
}

# ---------------------------------------------------------------------------
# Fallback DB download (mirror)
# ---------------------------------------------------------------------------

.dnmb_eggnog_fallback_download <- function(data_dir) {
  base_urls <- c(
    "http://eggnog5.embl.de/download/emapperdb-5.0.2",
    "http://eggnogdb.embl.de/download/emapperdb-5.0.2"
  )

  files_needed <- list(
    list(name = "eggnog_proteins.dmnd.gz", final = "eggnog_proteins.dmnd", decompress = "gunzip"),
    list(name = "eggnog.db.gz",            final = "eggnog.db",            decompress = "gunzip"),
    list(name = "eggnog.taxa.tar.gz",      final = "eggnog.taxa.db",       decompress = "tar")
  )

  for (finfo in files_needed) {
    final_path <- file.path(data_dir, finfo$final)
    if (file.exists(final_path) && file.info(final_path)$size > 0) next

    gz_path <- file.path(data_dir, finfo$name)
    downloaded <- FALSE

    for (base_url in base_urls) {
      url <- paste0(base_url, "/", finfo$name)
      message("[DNMB EggNOG] Downloading ", finfo$name, " from ", base_url)
      # Prefer wget/curl for large files (better resume + timeout handling)
      dl <- 1L
      wget_check <- dnmb_detect_binary("wget", required = FALSE)
      curl_check <- dnmb_detect_binary("curl", required = FALSE)
      if (isTRUE(wget_check$found)) {
        run <- dnmb_run_external("wget", args = c("-c", "-q", "--show-progress", "-O", gz_path, url), required = FALSE)
        if (isTRUE(run$ok)) dl <- 0L
      } else if (isTRUE(curl_check$found)) {
        run <- dnmb_run_external("curl", args = c("-L", "-C", "-", "-o", gz_path, url), required = FALSE)
        if (isTRUE(run$ok)) dl <- 0L
      } else {
        old_timeout <- getOption("timeout")
        options(timeout = 7200)
        dl <- tryCatch(
          utils::download.file(url, gz_path, mode = "wb", quiet = FALSE, method = "libcurl"),
          error = function(e) 1L
        )
        options(timeout = old_timeout)
      }
      if (dl == 0L && file.exists(gz_path) && file.info(gz_path)$size > 1000) {
        downloaded <- TRUE
        break
      }
    }

    if (!downloaded) {
      message("[DNMB EggNOG] Failed to download ", finfo$name)
      return(FALSE)
    }

    # Decompress
    message("[DNMB EggNOG] Decompressing ", finfo$name)
    if (identical(finfo$decompress, "gunzip")) {
      run <- dnmb_run_external("gunzip", args = c("-f", gz_path), required = FALSE)
    } else if (identical(finfo$decompress, "tar")) {
      run <- dnmb_run_external("tar", args = c("-xzf", gz_path, "-C", data_dir), required = FALSE)
      unlink(gz_path, force = TRUE)
    }
  }

  # Verify
  required <- c("eggnog_proteins.dmnd", "eggnog.db", "eggnog.taxa.db")
  base::all(base::vapply(
    base::file.path(data_dir, required),
    .dnmb_nonempty_file,
    logical(1)
  ))
}

# ---------------------------------------------------------------------------
# Auto-install emapper.py
# ---------------------------------------------------------------------------

.dnmb_eggnog_auto_install_emapper <- function(verbose = TRUE) {
  required_version <- .dnmb_eggnog_required_version()
  conda_spec <- base::paste0("eggnog-mapper=", required_version)
  installers <- list(
    list(name = "mamba", cmd = "mamba", args = c("install", "-y", "-c", "bioconda", "-c", "conda-forge", conda_spec)),
    list(name = "conda", cmd = "conda", args = c("install", "-y", "-c", "bioconda", "-c", "conda-forge", conda_spec))
  )

  for (inst in installers) {
    check <- dnmb_detect_binary(inst$cmd, required = FALSE)
    if (!isTRUE(check$found)) next

    if (isTRUE(verbose)) message("[DNMB EggNOG] Trying ", inst$name, " for eggNOG-mapper ", required_version)
    run <- dnmb_run_external(inst$cmd, args = inst$args, required = FALSE)
    if (isTRUE(run$ok)) {
      # Re-detect after install
      emapper <- dnmb_detect_emapper(required = FALSE)
      if (isTRUE(emapper$found) && .dnmb_eggnog_version_is_supported(emapper$version)) {
        if (isTRUE(verbose)) message("[DNMB EggNOG] Installed via ", inst$name, " (v", emapper$version %||% "?", ")")
        return(emapper)
      }
    }
  }

  if (isTRUE(verbose)) {
    message(
      "[DNMB EggNOG] Auto-install failed. Please install manually: conda install -c bioconda -c conda-forge ",
      conda_spec
    )
  }
  list(found = FALSE, path = NA_character_, version = NA_character_)
}

# ---------------------------------------------------------------------------
# Core module function
# ---------------------------------------------------------------------------

#' Run eggnog-mapper module
#'
#' Executes emapper.py on the query proteome and parses the output into a
#' standardized DNMB module result.
#'
#' @param genes A data frame (genbank_table) with at least `locus_tag`.
#' @param output_dir Output directory for module files.
#' @param data_dir eggnog-mapper database directory. If NULL, auto-detected.
#' @param db_name Taxonomic scope database (default: "bact").
#' @param cpu Number of threads (default: 1).
#' @param genbank Path to the GenBank file (for protein FASTA extraction).
#' @param install Whether to install missing databases.
#' @param sensmode Sensitivity mode for diamond ("fast", "mid-sensitive",
#'   "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive").
#'   Default: "mid-sensitive".
#' @param evalue E-value threshold (default: 1e-3).
#' @param verbose Print progress messages.
#' @param cache_root Optional DNMB cache root.
#'
#' @return A list with `ok`, `status`, `output_table`, `hits`, `files`.
#' @export
dnmb_run_eggnog_module <- function(genes,
                                   output_dir = NULL,
                                   data_dir = NULL,
                                   db_name = .dnmb_eggnog_default_db(),
                                   cpu = 1L,
                                   genbank = NULL,
                                   install = TRUE,
                                   sensmode = "mid-sensitive",
                                   evalue = 1e-3,
                                   verbose = TRUE,
                                   cache_root = NULL) {

  module_name <- .dnmb_eggnog_module_name()

  # Resolve output directory
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), paste0("dnmb_module_", tolower(module_name)))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Status tracking
  status <- list()
  .add_status <- function(component, ok, detail = NULL) {
    status[[length(status) + 1L]] <<- data.frame(
      component = component,
      ok = ok,
      detail = detail %||% "",
      stringsAsFactors = FALSE
    )
  }

  # 0. Auto-detect tax_scope from GenBank if db_name is "auto"
  if (identical(db_name, "auto")) {
    detected <- .dnmb_eggnog_detect_tax_scope(genbank)
    if (!is.null(detected)) {
      db_name <- detected
      if (isTRUE(verbose)) message("[DNMB EggNOG] Auto-detected tax_scope: ", db_name)
    } else {
      db_name <- "bact"
      if (isTRUE(verbose)) message("[DNMB EggNOG] Could not detect taxonomy, using default: ", db_name)
    }
  }

  # 1. Detect emapper.py — auto-install if missing and install=TRUE
  emapper <- dnmb_detect_emapper(required = FALSE)
  emapper_supported <- isTRUE(emapper$found) && .dnmb_eggnog_version_is_supported(emapper$version)
  if (!isTRUE(emapper_supported) && isTRUE(install)) {
    if (isTRUE(verbose)) {
      message(
        "[DNMB EggNOG] eggNOG-mapper ", .dnmb_eggnog_required_version(),
        " is required; attempting a pinned install (found ", emapper$version %||% "none", ")."
      )
    }
    emapper <- .dnmb_eggnog_auto_install_emapper(verbose = verbose)
  }
  emapper_supported <- isTRUE(emapper$found) && .dnmb_eggnog_version_is_supported(emapper$version)
  if (!isTRUE(emapper_supported)) {
    .add_status("emapper_detection", FALSE,
                paste0(
                  "eggNOG-mapper ", .dnmb_eggnog_required_version(), " is required; found ",
                  emapper$version %||% "none",
                  if (isTRUE(install)) "; pinned auto-install failed" else ""
                ))
    return(list(
      ok = FALSE,
      status = do.call(rbind, status),
      output_table = NULL, hits = NULL, files = list()
    ))
  }
  .add_status("emapper_detection", TRUE,
              paste0("emapper.py v", emapper$version %||% "unknown"))

  # 2. Check databases
  db_check <- dnmb_eggnog_ensure_database(
    data_dir = data_dir, db_name = db_name, install = install, cache_root = cache_root
  )
  if (!isTRUE(db_check$ok)) {
    .add_status("database", FALSE, db_check$error %||% "databases not found")
    return(list(
      ok = FALSE,
      status = do.call(rbind, status),
      output_table = NULL, hits = NULL, files = list()
    ))
  }
  .add_status("database", TRUE, paste0("data_dir=", db_check$data_dir))

  # 3. Prepare query protein FASTA from genes (locus_tag as header)
  fasta_path <- file.path(output_dir, "eggnog_query_proteins.faa")
  existing_faa <- dnmb_resolve_query_faa(
    genbank = genbank, output_dir = output_dir,
    fallback_filename = basename(fasta_path)
  )
  if (!is.null(existing_faa) && .dnmb_can_reuse_query_fasta(existing_faa, genes)) {
    proteins <- .dnmb_prepare_query_proteins(genes)
    fasta <- list(path = existing_faa, n = nrow(proteins), proteins = proteins)
    query_faa <- existing_faa
  } else {
    fasta <- .dnmb_write_query_fasta(genes, fasta_path)
    query_faa <- fasta_path
  }
  if (fasta$n == 0L) {
    .add_status("query_fasta", FALSE, "no protein sequences in genes")
    return(list(
      ok = TRUE,
      status = do.call(rbind, status),
      output_table = data.frame(), hits = .dnmb_eggnog_empty_hits(),
      files = list(query_fasta = query_faa)
    ))
  }
  .add_status("query_fasta", TRUE, paste0("proteins=", fasta$n))

  # 4. Run emapper.py
  output_prefix <- file.path(output_dir, "eggnog_out")
  annotations_file <- paste0(output_prefix, ".emapper.annotations")

  # Check if cached result exists and is non-empty
  if (file.exists(annotations_file) && file.info(annotations_file)$size > 0) {
    if (isTRUE(verbose)) {
      message("[DNMB EggNOG] Using cached result: ", basename(annotations_file))
    }
    .add_status("emapper_run", TRUE, "cached")
  } else {
    if (isTRUE(verbose)) {
      message("[DNMB EggNOG] Running emapper.py (cpu=", cpu, ", db=", db_name, ")")
    }

    emapper_args <- c(
      "-i", query_faa,
      "--output", output_prefix,
      "-m", "diamond",
      "--cpu", as.character(as.integer(cpu)),
      "--mp_start_method", "forkserver",
      "--data_dir", db_check$data_dir,
      "--dmnd_ignore_warnings",
      "--sensmode", "fast",
      "--evalue", format(evalue, scientific = TRUE),
      "--score", "60",
      "--pident", "40",
      "--query_cover", "20",
      "--subject_cover", "20",
      "--itype", "proteins",
      "--target_orthologs", "all",
      "--go_evidence", "non-electronic",
      "--pfam_realign", "none",
      "--override",
      "--no_file_comments",
      "--excel"
    )

    # Add taxonomic scope (auto-detected or user-specified)
    if (!is.null(db_name) && nzchar(db_name)) {
      emapper_args <- c(emapper_args, "--tax_scope", db_name)
    }

    emapper_cmd <- if (!is.na(emapper$path) && nzchar(emapper$path)) emapper$path else "emapper.py"
    emapper_path <- base::normalizePath(emapper_cmd, winslash = "/", mustWork = FALSE)
    emapper_env_path <- base::paste(
      unique(c(
        base::dirname(emapper_path),
        base::strsplit(base::Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
      )),
      collapse = .Platform$path.sep
    )
    run <- dnmb_run_external(
      emapper_cmd,
      args = emapper_args,
      env = c(PATH = emapper_env_path),
      required = FALSE,
      stream_stderr = TRUE
    )

    if (!isTRUE(run$ok) || !file.exists(annotations_file)) {
      stderr_msg <- paste(utils::tail(run$stderr, 5), collapse = "; ")
      .add_status("emapper_run", FALSE, paste0("emapper.py failed: ", stderr_msg))
      return(list(
        ok = FALSE,
        status = do.call(rbind, status),
        output_table = NULL, hits = NULL,
        files = list(query_fasta = query_faa)
      ))
    }
    .add_status("emapper_run", TRUE, "completed")
  }

  # 5. Parse output
  hits <- .dnmb_eggnog_parse_annotations(annotations_file)
  .add_status("parse", TRUE, paste0(nrow(hits), " annotations"))

  # 6. Build output table (join with genes by locus_tag)
  output_table <- .dnmb_eggnog_build_output_table(genes, hits)
  .add_status("output_table", TRUE, paste0(nrow(output_table), " rows"))

  list(
    ok = TRUE,
    status = do.call(rbind, status),
    output_table = output_table,
    hits = hits,
    files = list(
      query_fasta = query_faa,
      annotations = annotations_file,
      orthologs = paste0(output_prefix, ".emapper.seed_orthologs"),
      hits_file = paste0(output_prefix, ".emapper.hits")
    )
  )
}

# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

.dnmb_eggnog_parse_annotations <- function(annotations_file) {
  # emapper.py with --no_file_comments outputs a clean TSV

  # With comments: lines starting with # are header/footer
  lines <- readLines(annotations_file, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(trimws(lines))]

  if (length(lines) == 0L) {
    return(.dnmb_eggnog_empty_hits())
  }

  # Parse TSV
  con <- textConnection(paste(lines, collapse = "\n"))
  on.exit(close(con), add = TRUE)

  expected_cols <- .dnmb_eggnog_tsv_columns()
  hits <- utils::read.delim(con, header = FALSE, stringsAsFactors = FALSE,
                            quote = "", comment.char = "",
                            col.names = if (length(lines) > 0) {
                              # Read first line to determine column count
                              ncols <- length(strsplit(lines[1], "\t")[[1]])
                              if (ncols == length(expected_cols)) {
                                expected_cols
                              } else if (ncols < length(expected_cols)) {
                                expected_cols[seq_len(ncols)]
                              } else {
                                c(expected_cols, paste0("X", seq_len(ncols - length(expected_cols))))
                              }
                            } else expected_cols)

  # Clean query column: extract protein_id / locus_tag
  if ("query" %in% names(hits)) {
    hits$query <- .dnmb_eggnog_clean_query(hits$query)
  }

  # Coerce numeric columns
  if ("evalue" %in% names(hits)) {
    hits$evalue <- suppressWarnings(as.numeric(hits$evalue))
  }
  if ("score" %in% names(hits)) {
    hits$score <- suppressWarnings(as.numeric(hits$score))
  }

  # Replace "-" with NA for annotation columns
  dash_cols <- intersect(
    c("COG_category", "Description", "Preferred_name", "GOs", "EC",
      "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "PFAMs"),
    names(hits)
  )
  for (col in dash_cols) {
    hits[[col]][hits[[col]] == "-"] <- NA_character_
  }

  hits
}

.dnmb_eggnog_empty_hits <- function() {
  cols <- .dnmb_eggnog_tsv_columns()
  df <- data.frame(matrix(ncol = length(cols), nrow = 0), stringsAsFactors = FALSE)
  names(df) <- cols
  df
}

.dnmb_eggnog_clean_query <- function(query) {
  # Handle NCBI-style protein IDs from GenBank FASTA extraction
  # e.g., "gnl|extdb|WP_014789012.1" -> "WP_014789012.1"
  # e.g., "lcl|NC_123456.1_prot_WP_014789012.1_123" -> "WP_014789012.1"
  q <- query
  q <- gsub("^gnl\\|", "", q)
  q <- gsub("^extdb:", "", q)
  q <- gsub("\\|", ":", q)

  # Handle lcl:NC_xxx_prot_PROTID_NNN pattern
  lcl_pat <- "^lcl:(AC|NC|BG|NT|NW|NZ)_([a-zA-Z]+)?[0-9]+\\.[0-9]+_prot_"
  has_lcl <- grepl(lcl_pat, q)
  if (any(has_lcl)) {
    q[has_lcl] <- gsub(lcl_pat, "", q[has_lcl])
    q[has_lcl] <- gsub("_[0-9]{1,4}$", "", q[has_lcl])
  }

  q
}

# ---------------------------------------------------------------------------
# Output table builder
# ---------------------------------------------------------------------------

.dnmb_eggnog_build_output_table <- function(genes, hits) {
  if (nrow(hits) == 0L) {
    out <- genes[, "locus_tag", drop = FALSE]
    for (col in c("seed_ortholog", "evalue", "score", "eggNOG_OGs",
                  "COG_category", "Description", "Preferred_name",
                  "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
                  "PFAMs", "COG_color", "COG_legend")) {
      out[[col]] <- NA_character_
    }
    return(out)
  }

  # Match hits to genes by locus_tag (DNMB FASTA uses locus_tag as header)
  hit_key <- hits$query
  match_idx <- match(genes$locus_tag, hit_key)

  # Fallback: if query looks like protein_id (e.g. WP_xxx), try protein_id
  if (all(is.na(match_idx)) && "protein_id" %in% names(genes)) {
    match_idx <- match(genes$protein_id, hit_key)
  }

  # Build output columns
  out_cols <- c("seed_ortholog", "evalue", "score", "eggNOG_OGs",
                "COG_category", "Description", "Preferred_name",
                "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
                "PFAMs")

  out <- genes[, "locus_tag", drop = FALSE]

  for (col in out_cols) {
    if (col %in% names(hits)) {
      out[[col]] <- hits[[col]][match_idx]
    } else {
      out[[col]] <- NA_character_
    }
  }

  # Add COG color and legend
  cog_first <- substr(out$COG_category, 1, 1)
  cog_first[is.na(cog_first) | cog_first == ""] <- NA_character_

  cog_colors <- .dnmb_eggnog_cog_colors()
  cog_legends <- .dnmb_eggnog_cog_legends()

  out$COG_color <- cog_colors[cog_first]
  out$COG_legend <- cog_legends[cog_first]

  out
}

# ---------------------------------------------------------------------------
# Output table for module_runner (fallback)
# ---------------------------------------------------------------------------

.dnmb_eggnog_output_table <- function(genes, hits) {
  .dnmb_eggnog_build_output_table(genes, hits)
}

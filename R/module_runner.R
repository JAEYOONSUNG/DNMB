.dnmb_module_clean_annotation_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^gnl\\|", "", x)
  x <- gsub("\\|", ":", x)
  x <- gsub("lcl:(AC|NC|BG|NT|NW|NZ)_([A-Za-z]+)?[0-9]+\\.[0-9]+_prot_", "", x)
  x <- gsub("_[0-9]{1,4}$", "", x)
  x <- gsub("^extdb:", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

.dnmb_module_default_prophage_backend <- function() {
  backend_fn <- get0(".dnmb_prophage_default_backend", mode = "function", inherits = TRUE)
  if (is.function(backend_fn)) {
    return(backend_fn())
  }
  "phispy"
}

.dnmb_module_gapmind_name <- function(version) {
  gapmind_name_fn <- get0(".dnmb_gapmind_database_name", mode = "function", inherits = TRUE)
  if (is.function(gapmind_name_fn)) {
    return(gapmind_name_fn(version))
  }
  if (identical(as.character(version)[1], "aa")) {
    return("GapMindAA")
  }
  "GapMindCarbon"
}

.dnmb_module_optional_long_columns <- function() {
  c(
    "query",
    "source",
    "family_system",
    "family_id",
    "hit_label",
    "enzyme_role",
    "evidence_mode",
    "substrate_label",
    "support",
    "typing_eligible"
  )
}

.dnmb_module_empty_optional_long_table <- function() {
  data.frame(
    query = character(),
    source = character(),
    family_system = character(),
    family_id = character(),
    hit_label = character(),
    enzyme_role = character(),
    evidence_mode = character(),
    substrate_label = character(),
    support = character(),
    typing_eligible = logical(),
    stringsAsFactors = FALSE
  )
}

.dnmb_module_detect_genbank <- function(path = getwd()) {
  candidates <- list.files(
    path,
    pattern = "\\.gbk$|\\.gb$|\\.gbff$",
    full.names = TRUE
  )
  candidates <- candidates[!grepl("(^|/|\\\\)~\\$", candidates)]
  if (!length(candidates)) {
    return(NULL)
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

.dnmb_write_protein_fasta <- function(proteins, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!nrow(proteins)) {
    writeLines(character(), con = path)
    return(invisible(path))
  }

  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(nrow(proteins))) {
    seq_clean <- gsub("[^A-Za-z*]", "", proteins$protein_seq[[i]])
    seq_clean <- toupper(seq_clean)
    seq_clean <- gsub("\\*", "", seq_clean)
    if (!nzchar(seq_clean)) {
      next
    }
    writeLines(paste0(">", proteins$protein_label[[i]]), con = con)
    writeLines(seq_clean, con = con)
  }
  invisible(path)
}

.dnmb_normalize_translation <- function(x) {
  x <- as.character(x)
  x <- gsub("[[:space:]]+", "", x)
  x <- toupper(x)
  x <- gsub("\\*", "", x)
  x <- gsub("[^A-Z]", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

.dnmb_prepare_query_proteins <- function(genes) {
  if (!is.data.frame(genes)) {
    stop("`genes` must be a data frame.", call. = FALSE)
  }
  required_columns <- c("locus_tag", "translation")
  missing_columns <- setdiff(required_columns, names(genes))
  if (length(missing_columns)) {
    stop(
      "`genes` must contain columns: ",
      paste(required_columns, collapse = ", "),
      ". Missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  proteins <- tibble::as_tibble(genes)
  proteins$locus_tag <- .dnmb_module_clean_annotation_key(proteins$locus_tag)
  proteins$translation <- .dnmb_normalize_translation(proteins$translation)
  proteins <- proteins[!is.na(proteins$locus_tag) & nzchar(proteins$locus_tag) & !is.na(proteins$translation), , drop = FALSE]
  proteins <- proteins[!duplicated(proteins$locus_tag), , drop = FALSE]
  rownames(proteins) <- NULL
  proteins
}

.dnmb_write_query_fasta <- function(genes, path) {
  proteins <- .dnmb_prepare_query_proteins(genes)
  if (!nrow(proteins)) {
    writeLines(character(), con = path)
    return(list(path = path, proteins = proteins, n = 0L))
  }

  fasta_tbl <- tibble::tibble(
    protein_label = proteins$locus_tag,
    protein_seq = proteins$translation
  )
  .dnmb_write_protein_fasta(fasta_tbl, path)
  list(path = path, proteins = proteins, n = nrow(proteins))
}

.dnmb_fasta_headers <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(character())
  }

  headers <- readLines(path, warn = FALSE)
  headers <- headers[grepl("^>", headers)]
  headers <- sub("^>", "", headers)
  headers <- trimws(headers)
  headers <- headers[nzchar(headers)]
  headers
}

.dnmb_can_reuse_query_fasta <- function(path, genes) {
  proteins <- .dnmb_prepare_query_proteins(genes)
  headers <- .dnmb_fasta_headers(path)
  if (!length(headers) || length(headers) != nrow(proteins)) {
    return(FALSE)
  }

  identical(headers, proteins$locus_tag)
}

.dnmb_best_module_hits <- function(hits) {
  if (is.null(hits) || !is.data.frame(hits) || !nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  required <- .dnmb_module_optional_long_columns()
  extra_columns <- setdiff(names(hits), required)
  missing <- setdiff(required, names(hits))
  if (length(missing)) {
    stop("`hits` is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  hits <- as.data.frame(hits, stringsAsFactors = FALSE)
  hits$query <- .dnmb_module_clean_annotation_key(hits$query)
  ordered_columns <- c(required, extra_columns)
  hits <- hits[!is.na(hits$query), ordered_columns, drop = FALSE]
  if (!nrow(hits)) {
    return(.dnmb_module_empty_optional_long_table())
  }

  hit_rank <-
    ifelse(!is.na(hits$family_id) & nzchar(hits$family_id), 100L, 0L) +
    ifelse(as.character(hits$evidence_mode) == "direct", 10L, 0L) +
    ifelse(isTRUE(hits$typing_eligible) | hits$typing_eligible %in% TRUE, 1L, 0L)
  hit_order <- seq_len(nrow(hits))
  hits <- hits[order(hits$query, -hit_rank, hit_order), , drop = FALSE]
  hits <- hits[!duplicated(hits$query), ordered_columns, drop = FALSE]
  rownames(hits) <- NULL
  hits
}

.dnmb_na_vector_like <- function(x, n) {
  if (is.logical(x)) {
    return(rep(NA, n))
  }
  if (is.integer(x)) {
    return(rep(NA_integer_, n))
  }
  if (is.numeric(x)) {
    return(rep(NA_real_, n))
  }
  rep(NA_character_, n)
}

.dnmb_module_output_table <- function(genes, hits) {
  genes <- as.data.frame(genes, stringsAsFactors = FALSE)
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  base_cols <- intersect(
    c("locus_tag", "gene", "product", "protein_id", "contig", "start", "end", "direction"),
    names(genes)
  )
  out <- genes[, base_cols, drop = FALSE]
  out$family_id <- NA_character_
  out$hit_label <- NA_character_
  out$enzyme_role <- NA_character_
  out$evidence_mode <- NA_character_
  out$substrate_label <- NA_character_
  out$support <- NA_character_
  out$typing_eligible <- NA

  best_hits <- .dnmb_best_module_hits(hits)
  if (!nrow(best_hits)) {
    return(out)
  }

  extra_columns <- setdiff(names(best_hits), .dnmb_module_optional_long_columns())
  if (length(extra_columns)) {
    for (column_name in extra_columns) {
      out[[column_name]] <- .dnmb_na_vector_like(best_hits[[column_name]], nrow(out))
    }
  }

  matches <- match(out$locus_tag, best_hits$query)
  include <- !is.na(matches)
  out$family_id[include] <- as.character(best_hits$family_id[matches[include]])
  out$hit_label[include] <- as.character(best_hits$hit_label[matches[include]])
  out$enzyme_role[include] <- as.character(best_hits$enzyme_role[matches[include]])
  out$evidence_mode[include] <- as.character(best_hits$evidence_mode[matches[include]])
  out$substrate_label[include] <- as.character(best_hits$substrate_label[matches[include]])
  out$support[include] <- as.character(best_hits$support[matches[include]])
  out$typing_eligible[include] <- best_hits$typing_eligible[matches[include]]
  if (length(extra_columns)) {
    for (column_name in extra_columns) {
      out[[column_name]][include] <- best_hits[[column_name]][matches[include]]
    }
  }
  out
}

.dnmb_resolve_run_genes <- function(genbank_table = NULL) {
  genes <- genbank_table
  if (is.null(genes)) {
    genes <- get0("genbank_table", envir = .GlobalEnv, inherits = FALSE)
  }
  if (is.null(genes) || !is.data.frame(genes)) {
    stop("No active genbank_table is available for module execution.", call. = FALSE)
  }
  genes
}

.dnmb_module_output_dir <- function(db, output_dir = NULL) {
  if (!is.null(output_dir) && nzchar(trimws(as.character(output_dir)[1]))) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(output_dir, winslash = "/", mustWork = FALSE))
  }

  out <- file.path(getwd(), paste0("dnmb_module_", tolower(db)))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  normalizePath(out, winslash = "/", mustWork = FALSE)
}

.dnmb_module_status_detail <- function(status_table) {
  if (is.null(status_table) || !is.data.frame(status_table) || !nrow(status_table)) {
    return(NA_character_)
  }

  parts <- paste0(
    as.character(status_table$component),
    "=",
    as.character(status_table$status),
    ifelse(
      is.na(status_table$detail) | !nzchar(trimws(as.character(status_table$detail))),
      "",
      paste0("(", trimws(as.character(status_table$detail)), ")")
    )
  )
  paste(parts, collapse = "; ")
}

.dnmb_module_try_run <- function(module_name, run_fn) {
  message("[DNMB module] ", module_name, " started")
  start_time <- proc.time()
  result <- tryCatch(
    run_fn(),
    error = function(e) {
      list(
        ok = FALSE,
        status = tibble::tibble(
          component = as.character(module_name)[1],
          status = "failed",
          detail = conditionMessage(e)
        ),
        output_table = data.frame(),
        hits = .dnmb_module_empty_optional_long_table()
      )
    }
  )
  elapsed <- round((proc.time() - start_time)[["elapsed"]], 1)
  if (isTRUE(result$ok)) {
    n_hits <- if (!is.null(result$hits) && is.data.frame(result$hits)) nrow(result$hits) else 0L
    message("[DNMB module] ", module_name, " completed (", n_hits, " hits, ", elapsed, "s)")
  } else {
    fail_detail <- if (!is.null(result$status) && is.data.frame(result$status) && "detail" %in% names(result$status)) {
      paste(result$status$detail[nzchar(result$status$detail)], collapse = "; ")
    } else ""
    message("[DNMB module] ", module_name, " failed (", elapsed, "s)", if (nzchar(fail_detail)) paste0(": ", fail_detail) else "")
  }

  if (!is.list(result) || is.null(result$ok)) {
    return(list(
      ok = FALSE,
      status = tibble::tibble(
        component = as.character(module_name)[1],
        status = "failed",
        detail = "Module returned an invalid result object."
      ),
      output_table = data.frame(),
      hits = .dnmb_module_empty_optional_long_table()
    ))
  }

  result
}

#' Run DNMB database modules
#'
#' Executes live DNMB module workflows against the active `genbank_table` or a
#' supplied gene table and returns per-module results or a merged locus-level
#' table.
#'
#' Live execution is currently implemented for `dbCAN`, `MEROPS`, `CLEAN`, `PAZy`, `GapMind`, `DefenseFinder`, `PADLOC`, `DefensePredictor`, `ISelement`, and `Prophage`.
#'
#' @param db Optional character vector of module names to run.
#' @param module_dbCAN,module_MEROPS,module_CLEAN,module_PAZy,module_GapMind,module_DefenseFinder,module_PADLOC,module_DefensePredictor,module_ISelement,module_Prophage Logical toggles used when
#'   `db` is not supplied.
#' @param module_Prophage_backend Character backend for `Prophage`; one of
#'   `phispy`, `virsorter2`, or `pide`.
#' @param genbank_table Optional gene table to use instead of the global
#'   `genbank_table`.
#' @param genbank Optional GenBank path used for query FASTA reuse.
#' @param output_dir Optional output directory for module artifacts.
#' @param module_version Optional module version string.
#' @param module_cache_root Optional cache root override.
#' @param module_install Logical; install module assets when missing.
#' @param module_base_url Optional module asset base URL override.
#' @param module_asset_urls Optional named local paths or URLs for module
#'   assets.
#' @param module_cpu Integer thread count for external tools.
#' @param merge Logical; when `TRUE`, return a merged locus-level table instead
#'   of a list of per-module runs.
#' @param verbose Logical; emit progress messages.
#'
#' @return A named list of `dnmb_module_run` objects, or a merged data frame
#'   when `merge = TRUE`.
#' @export
run_module_set <- function(db = NULL,
                           module_dbCAN = TRUE,
                           module_MEROPS = TRUE,
                           module_CLEAN = .dnmb_cuda_default_module(),
                           module_PAZy = TRUE,
                           module_GapMind = TRUE,
                           module_DefenseFinder = TRUE,
                           module_PADLOC = TRUE,
                           module_DefensePredictor = TRUE,
                           module_REBASEfinder = TRUE,
                           module_ISelement = TRUE,
                           module_Prophage = FALSE,
                           module_PhiSpy = TRUE,
                           module_VirSorter2 = FALSE,
                           module_PIDE = .dnmb_cuda_default_module(),
                           module_EggNOG = TRUE,
                           module_Prophage_backend = .dnmb_module_default_prophage_backend(),
                           genbank_table = NULL,
                           genbank = NULL,
                           output_dir = NULL,
                           module_version = NULL,
                           module_cache_root = NULL,
                           module_install = TRUE,
                           module_base_url = NULL,
                           module_asset_urls = NULL,
                           module_cpu = .dnmb_default_cpu(),
                           merge = FALSE,
                           verbose = TRUE,
                           iselement_analysis_depth = "standard",
                           iselement_related_genbanks = NULL,
                           iselement_related_metadata = NULL,
                           iselement_auto_discover_related = TRUE,
                           iselement_max_related = 5L) {
  # Legacy `module_Prophage` shim: route onto a specific backend flag based on
  # module_Prophage_backend, so older callers keep working without sprouting a
  # duplicate Prophage branch.
  if (isTRUE(module_Prophage)) {
    legacy_backend <- tryCatch(
      .dnmb_prophage_normalize_backend(module_Prophage_backend),
      error = function(e) "phispy"
    )
    warning(
      "`module_Prophage = TRUE` is deprecated; use `module_PhiSpy`, ",
      "`module_VirSorter2`, or `module_PIDE` instead.",
      call. = FALSE
    )
    if (identical(legacy_backend, "phispy")) {
      module_PhiSpy <- TRUE
    } else if (identical(legacy_backend, "virsorter2")) {
      module_VirSorter2 <- TRUE
    } else if (identical(legacy_backend, "pide")) {
      module_PIDE <- TRUE
    }
    module_Prophage <- FALSE
  }

  if (is.null(db)) {
    db <- names(which(c(
      EggNOG = isTRUE(module_EggNOG),
      CLEAN = isTRUE(module_CLEAN),
      DefenseFinder = isTRUE(module_DefenseFinder),
      PADLOC = isTRUE(module_PADLOC),
      DefensePredictor = isTRUE(module_DefensePredictor),
      REBASEfinder = isTRUE(module_REBASEfinder),
      GapMind = isTRUE(module_GapMind),
      MEROPS = isTRUE(module_MEROPS),
      dbCAN = isTRUE(module_dbCAN),
      PAZy = isTRUE(module_PAZy),
      ISelement = isTRUE(module_ISelement),
      PhiSpy = isTRUE(module_PhiSpy),
      VirSorter2 = isTRUE(module_VirSorter2),
      PIDE = isTRUE(module_PIDE)
    )))
  } else {
    db <- as.character(db)
    if ("Prophage" %in% db) {
      legacy_backend <- tryCatch(
        .dnmb_prophage_normalize_backend(module_Prophage_backend),
        error = function(e) "phispy"
      )
      replacement <- switch(
        legacy_backend,
        phispy = "PhiSpy",
        virsorter2 = "VirSorter2",
        pide = "PIDE",
        "PhiSpy"
      )
      warning(
        "Database alias 'Prophage' is deprecated; use 'PhiSpy', ",
        "'VirSorter2', or 'PIDE' instead (routed to '", replacement, "').",
        call. = FALSE
      )
      db <- unique(c(setdiff(db, "Prophage"), replacement))
    }
  }

  db <- unique(as.character(db))
  db <- db[!is.na(db) & nzchar(trimws(db))]
  if ("GapMind" %in% db) {
    gapmind_index <- match("GapMind", db)
    before <- if (gapmind_index > 1L) db[seq_len(gapmind_index - 1L)] else character()
    after <- if (gapmind_index < length(db)) db[seq.int(gapmind_index + 1L, length(db))] else character()
    db <- c(
      before,
      "GapMindAA",
      "GapMindCarbon",
      after
    )
    db <- db[nzchar(db)]
  }
  if (!length(db)) {
    if (isTRUE(merge)) {
      return(data.frame())
    }
    return(list())
  }

  unsupported <- setdiff(db, c("dbCAN", "MEROPS", "CLEAN", "PAZy", "GapMindAA", "GapMindCarbon", "DefenseFinder", "PADLOC", "DefensePredictor", "REBASEfinder", "ISelement", "PhiSpy", "VirSorter2", "PIDE", "EggNOG"))
  if (length(unsupported)) {
    stop(
      "Live module execution is currently implemented only for dbCAN, MEROPS, CLEAN, PAZy, GapMind, DefenseFinder, PADLOC, DefensePredictor, ISelement, PhiSpy, VirSorter2, PIDE, and EggNOG. Unsupported module(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }

  genes <- .dnmb_resolve_run_genes(genbank_table = genbank_table)
  genbank <- genbank %||% .dnmb_module_detect_genbank(getwd())
  runs <- list()

  # DB freshness check for enabled modules
  if (isTRUE(verbose)) {
    db_version_map <- list(
      EggNOG = list(module = "eggnog", version = "data"),
      CLEAN = list(module = "clean", version = module_version %||% "split100"),
      DefenseFinder = list(module = "defensefinder", version = module_version %||% "current"),
      PADLOC = list(module = "padloc", version = module_version %||% "current"),
      DefensePredictor = list(module = "defensepredictor", version = module_version %||% "current"),
      GapMindAA = list(module = "gapmind", version = "aa"),
      GapMindCarbon = list(module = "gapmind", version = "carbon"),
      MEROPS = list(module = "merops", version = module_version %||% "current"),
      dbCAN = list(module = "dbcan", version = module_version %||% "current"),
      PAZy = list(module = "pazy", version = module_version %||% "current"),
      ISelement = list(module = "iselement", version = module_version %||% "current"),
      PhiSpy = list(module = "prophage", version = "phispy"),
      VirSorter2 = list(module = "prophage", version = "virsorter2"),
      PIDE = list(module = "prophage", version = "pide")
    )
    for (mod_name in intersect(db, names(db_version_map))) {
      info <- db_version_map[[mod_name]]
      dnmb_db_check_freshness(info$module, info$version, cache_root = module_cache_root, verbose = TRUE)
    }
  }

  ## --- 1. EggNOG ---
  if ("EggNOG" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── EggNOG ──")
    }
    eggnog_result <- .dnmb_module_try_run("EggNOG", function() {
      dnmb_run_eggnog_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("EggNOG", output_dir = output_dir),
        cpu = as.integer(module_cpu)[1],
        genbank = genbank,
        install = isTRUE(module_install),
        verbose = verbose
      )
    })
    if (!isTRUE(eggnog_result$ok)) {
      warning(
        "EggNOG module execution failed: ",
        .dnmb_module_status_detail(eggnog_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    runs$EggNOG <- structure(
      list(
        database = "EggNOG",
        output_table = eggnog_result$output_table %||% .dnmb_eggnog_output_table(genes = genes, hits = eggnog_result$hits),
        hits = eggnog_result$hits,
        module_result = eggnog_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 2. CLEAN ---
  if ("CLEAN" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── CLEAN ──")
    }
    clean_result <- .dnmb_module_try_run("CLEAN", function() {
      dnmb_run_clean_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("CLEAN", output_dir = output_dir),
        version = module_version %||% .dnmb_clean_default_version(),
        cache_root = module_cache_root,
        install = isTRUE(module_install),
        base_url = module_base_url %||% .dnmb_clean_default_base_url(),
        asset_urls = if (is.list(module_asset_urls) && "CLEAN" %in% names(module_asset_urls)) module_asset_urls[["CLEAN"]] else module_asset_urls,
        cpu = as.integer(module_cpu)[1],
        genbank = genbank
      )
    })
    if (!isTRUE(clean_result$ok)) {
      warning(
        "CLEAN module execution failed: ",
        .dnmb_module_status_detail(clean_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    clean_output_table <- clean_result$output_table %||% .dnmb_clean_output_table(genes = genes, hits = clean_result$hits)
    runs$CLEAN <- structure(
      list(
        database = "CLEAN",
        output_table = clean_output_table,
        hits = clean_result$hits,
        module_result = clean_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 3. DefenseFinder ---
  if ("DefenseFinder" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── DefenseFinder ──")
    }
    defensefinder_result <- .dnmb_module_try_run("DefenseFinder", function() {
      dnmb_run_defensefinder_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("DefenseFinder", output_dir = output_dir),
        version = module_version %||% .dnmb_defensefinder_default_version(),
        cache_root = module_cache_root,
        install = isTRUE(module_install),
        repo_url = .dnmb_defensefinder_default_repo_url(),
        asset_urls = if (is.list(module_asset_urls) && "DefenseFinder" %in% names(module_asset_urls)) module_asset_urls[["DefenseFinder"]] else module_asset_urls,
        cpu = as.integer(module_cpu)[1],
        genbank = genbank
      )
    })
    if (!isTRUE(defensefinder_result$ok)) {
      warning(
        "DefenseFinder module execution failed: ",
        .dnmb_module_status_detail(defensefinder_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    runs$DefenseFinder <- structure(
      list(
        database = "DefenseFinder",
        output_table = defensefinder_result$output_table %||% .dnmb_defensefinder_output_table(genes = genes, hits = defensefinder_result$hits),
        hits = defensefinder_result$hits,
        module_result = defensefinder_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 3a. PADLOC ---
  if ("PADLOC" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── PADLOC ──")
    }
    padloc_result <- .dnmb_module_try_run("PADLOC", function() {
      dnmb_run_padloc_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("PADLOC", output_dir = output_dir),
        version = module_version %||% .dnmb_padloc_default_version(),
        cache_root = module_cache_root,
        install = isTRUE(module_install),
        cpu = as.integer(module_cpu)[1]
      )
    })
    if (!isTRUE(padloc_result$ok)) {
      warning(
        "PADLOC module execution failed: ",
        .dnmb_module_status_detail(padloc_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
      runs$PADLOC <- structure(
        list(
          database = "PADLOC",
          output_table = padloc_result$output_table %||% .dnmb_padloc_output_table(genes = genes, hits = padloc_result$hits),
          hits = padloc_result$hits,
          module_result = padloc_result
        ),
        class = "dnmb_module_run"
      )
    }
  }

  ## --- 3aa. DefensePredictor ---
  if ("DefensePredictor" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── DefensePredictor ──")
    }
    defensepredictor_result <- .dnmb_module_try_run("DefensePredictor", function() {
      dnmb_run_defensepredictor_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("DefensePredictor", output_dir = output_dir),
        version = module_version %||% .dnmb_defensepredictor_default_version(),
        cache_root = module_cache_root,
        install = isTRUE(module_install)
      )
    })
    if (!isTRUE(defensepredictor_result$ok)) {
      warning(
        "DefensePredictor module execution failed: ",
        .dnmb_module_status_detail(defensepredictor_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
      runs$DefensePredictor <- structure(
        list(
          database = "DefensePredictor",
          output_table = defensepredictor_result$output_table %||% .dnmb_defensepredictor_output_table(genes = genes, hits = defensepredictor_result$hits),
          hits = defensepredictor_result$hits,
          module_result = defensepredictor_result
        ),
        class = "dnmb_module_run"
      )
    }
  }

  ## --- 3b. REBASEfinder ---
  if ("REBASEfinder" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── REBASEfinder ──")
    }
    rebasefinder_result <- .dnmb_module_try_run("REBASEfinder", function() {
      dnmb_run_rebasefinder_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("REBASEfinder", output_dir = output_dir),
        genbank = genbank,
        compare_rebase = TRUE,
        search_motifs = TRUE,
        install = isTRUE(module_install),
        cpu = as.integer(module_cpu)[1],
        verbose = verbose
      )
    })
    if (!isTRUE(rebasefinder_result$ok)) {
      warning(
        "REBASEfinder module execution failed: ",
        .dnmb_module_status_detail(rebasefinder_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    runs$REBASEfinder <- structure(
      list(
        database = "REBASEfinder",
        output_table = rebasefinder_result$output_table %||% .dnmb_rebasefinder_output_table(genes = genes, hits = rebasefinder_result$hits),
        hits = rebasefinder_result$hits,
        module_result = rebasefinder_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 4. GapMind (AA + Carbon) ---
  for (gapmind_version in c("aa", "carbon")) {
    gapmind_name <- .dnmb_module_gapmind_name(gapmind_version)
    if (!gapmind_name %in% db) {
      next
    }
    if (isTRUE(verbose)) {
      message("[DNMB] ── ", gapmind_name, " ──")
    }
    gapmind_asset_urls <- module_asset_urls
    if (is.list(module_asset_urls) && "GapMind" %in% names(module_asset_urls)) {
      gapmind_asset_urls <- module_asset_urls[["GapMind"]]
    }
    if (is.list(module_asset_urls) && gapmind_name %in% names(module_asset_urls)) {
      gapmind_asset_urls <- module_asset_urls[[gapmind_name]]
    }
    gapmind_result <- .dnmb_module_try_run(gapmind_name, function() {
      dnmb_run_gapmind_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir(gapmind_name, output_dir = output_dir),
        version = gapmind_version,
        cache_root = module_cache_root,
        install = isTRUE(module_install),
        repo_url = .dnmb_gapmind_default_repo_url(),
        resource_base_url = .dnmb_gapmind_default_resource_base_url(gapmind_version),
        asset_urls = gapmind_asset_urls,
        cpu = as.integer(module_cpu)[1],
        genbank = genbank
      )
    })
    if (!isTRUE(gapmind_result$ok)) {
      warning(
        gapmind_name,
        " module execution failed: ",
        .dnmb_module_status_detail(gapmind_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    runs[[gapmind_name]] <- structure(
      list(
        database = gapmind_name,
        output_table = gapmind_result$output_table %||% .dnmb_gapmind_output_table(genes = genes, hits = gapmind_result$hits, version = gapmind_version),
        hits = gapmind_result$hits,
        module_result = gapmind_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 5. MEROPS ---
  if ("MEROPS" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── MEROPS ──")
    }
    merops_result <- .dnmb_module_try_run("MEROPS", function() {
      dnmb_run_merops_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("MEROPS", output_dir = output_dir),
        version = module_version %||% .dnmb_merops_default_version(),
        cache_root = module_cache_root,
        install = isTRUE(module_install),
        base_url = module_base_url %||% .dnmb_merops_default_base_url(),
        asset_urls = if (is.list(module_asset_urls) && "MEROPS" %in% names(module_asset_urls)) module_asset_urls[["MEROPS"]] else module_asset_urls,
        cpu = as.integer(module_cpu)[1],
        genbank = genbank
      )
    })
    if (!isTRUE(merops_result$ok)) {
      warning(
        "MEROPS module execution failed: ",
        .dnmb_module_status_detail(merops_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    runs$MEROPS <- structure(
      list(
        database = "MEROPS",
        output_table = .dnmb_merops_output_table(genes = genes, hits = merops_result$hits),
        hits = merops_result$hits,
        module_result = merops_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 6. dbCAN ---
  if ("dbCAN" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── dbCAN ──")
    }
    dbcan_result <- .dnmb_module_try_run("dbCAN", function() {
      dnmb_run_dbcan_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("dbCAN", output_dir = output_dir),
        version = module_version %||% .dnmb_dbcan_default_version(),
        cache_root = module_cache_root,
        install = isTRUE(module_install),
        base_url = module_base_url %||% .dnmb_dbcan_default_base_url(),
        asset_urls = if (is.list(module_asset_urls) && "dbCAN" %in% names(module_asset_urls)) module_asset_urls[["dbCAN"]] else module_asset_urls,
        cpu = as.integer(module_cpu)[1],
        genbank = genbank
      )
    })
    if (!isTRUE(dbcan_result$ok)) {
      warning(
        "dbCAN module execution failed: ",
        .dnmb_module_status_detail(dbcan_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    dbcan_output_table <- dbcan_result$output_table %||% .dnmb_dbcan_output_table(genes = genes, hits = dbcan_result$hits)
    runs$dbCAN <- structure(
      list(
        database = "dbCAN",
        output_table = dbcan_output_table,
        hits = dbcan_result$hits,
        module_result = dbcan_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 7. PAZy ---
  if ("PAZy" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── PAZy ──")
    }
    pazy_result <- .dnmb_module_try_run("PAZy", function() {
      dnmb_run_pazy_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("PAZy", output_dir = output_dir),
        version = module_version %||% .dnmb_pazy_default_version(),
        cache_root = module_cache_root,
        install = isTRUE(module_install),
        metadata_url = .dnmb_pazy_default_metadata_url(),
        fasta_url = .dnmb_pazy_default_fasta_url(),
        asset_urls = if (is.list(module_asset_urls) && "PAZy" %in% names(module_asset_urls)) module_asset_urls[["PAZy"]] else module_asset_urls,
        cpu = as.integer(module_cpu)[1],
        genbank = genbank
      )
    })
    if (!isTRUE(pazy_result$ok)) {
      warning(
        "PAZy module execution failed: ",
        .dnmb_module_status_detail(pazy_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    runs$PAZy <- structure(
      list(
        database = "PAZy",
        output_table = pazy_result$output_table %||% .dnmb_pazy_output_table(genes = genes, hits = pazy_result$hits),
        hits = pazy_result$hits,
        module_result = pazy_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 8. ISelement ---
  if ("ISelement" %in% db) {
    if (isTRUE(verbose)) {
      message("[DNMB] ── ISelement ──")
    }
    iselement_result <- .dnmb_module_try_run("ISelement", function() {
      dnmb_run_iselement_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir("ISelement", output_dir = output_dir),
        genbank = genbank,
        detection_mode = "hybrid",
        analysis_depth = iselement_analysis_depth %||% "standard",
        related_genbanks = iselement_related_genbanks,
        related_metadata = iselement_related_metadata,
        auto_discover_related = iselement_auto_discover_related %||% TRUE,
        max_related = iselement_max_related %||% 5L,
        verbose = verbose
      )
    })
    if (!isTRUE(iselement_result$ok)) {
      warning(
        "ISelement module execution failed: ",
        .dnmb_module_status_detail(iselement_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
    runs$ISelement <- structure(
      list(
        database = "ISelement",
        output_table = iselement_result$output_table %||% .dnmb_iselement_output_table(genes = genes, hits = iselement_result$hits),
        hits = iselement_result$hits,
        module_result = iselement_result
      ),
      class = "dnmb_module_run"
    )
    }
  }

  ## --- 9. Prophage backends (PhiSpy / VirSorter2 / PIDE) ---
  prophage_backend_specs <- list(
    PhiSpy     = list(backend = "phispy",     label = "PhiSpy"),
    VirSorter2 = list(backend = "virsorter2", label = "VirSorter2"),
    PIDE       = list(backend = "pide",       label = "PIDE")
  )
  for (prophage_db in intersect(names(prophage_backend_specs), db)) {
    spec <- prophage_backend_specs[[prophage_db]]
    if (isTRUE(verbose)) {
      message("[DNMB] ── ", spec$label, " ──")
    }
    prophage_result <- .dnmb_module_try_run(spec$label, function() {
      dnmb_run_prophage_module(
        genes = genes,
        output_dir = .dnmb_module_output_dir(prophage_db, output_dir = output_dir),
        genbank = genbank,
        backend = spec$backend,
        cache_root = module_cache_root,
        install = isTRUE(module_install),
        cpu = as.integer(module_cpu)[1]
      )
    })
    if (!isTRUE(prophage_result$ok)) {
      warning(
        spec$label, " module execution failed: ",
        .dnmb_module_status_detail(prophage_result$status) %||% "unknown error",
        call. = FALSE
      )
    } else {
      runs[[prophage_db]] <- structure(
        list(
          database = spec$label,
          output_table = prophage_result$output_table %||% .dnmb_prophage_output_table(genes = genes, hits = prophage_result$hits),
          hits = prophage_result$hits,
          module_result = prophage_result
        ),
        class = "dnmb_module_run"
      )
    }
  }

  if (isTRUE(merge)) {
    if (!length(runs)) {
      warning(
        "No requested modules completed successfully; returning the input gene table unchanged.",
        call. = FALSE
      )
      return(as.data.frame(genes, stringsAsFactors = FALSE))
    }
    merge_module_results(runs)
  } else {
    runs
  }
}

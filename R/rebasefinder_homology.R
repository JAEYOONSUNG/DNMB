.dnmb_rebasefinder_homology_template_root <- function() {
  installed <- base::system.file(
    "extdata", "rebasefinder_structure_templates",
    package = "DNMB"
  )
  source_file <- tryCatch(
    utils::getSrcFilename(.dnmb_rebasefinder_homology_template_root, full.names = TRUE),
    error = function(e) ""
  )
  source_root <- if (base::length(source_file) && base::nzchar(source_file[[1]])) {
    base::dirname(base::dirname(source_file[[1]]))
  } else {
    ""
  }
  candidates <- base::unique(c(
    installed,
    if (base::nzchar(source_root)) base::file.path(source_root, "inst", "extdata", "rebasefinder_structure_templates") else "",
    base::file.path(base::getwd(), "inst", "extdata", "rebasefinder_structure_templates"),
    base::file.path(base::dirname(base::getwd()), "inst", "extdata", "rebasefinder_structure_templates")
  ))
  candidates <- candidates[base::nzchar(candidates) & base::dir.exists(candidates)]
  has_manifest <- base::vapply(candidates, function(path) {
    base::any(base::file.exists(base::file.path(
      path, c("manifest.tsv", "templates.tsv", "template_manifest.tsv")
    ))) || base::length(base::list.files(path, pattern = "manifest.*[.]tsv$", ignore.case = TRUE)) > 0L
  }, logical(1))
  candidates <- candidates[has_manifest]
  if (!base::length(candidates)) return(NA_character_)
  base::normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

.dnmb_rebasefinder_homology_manifest_path <- function(root = .dnmb_rebasefinder_homology_template_root()) {
  if (base::is.na(root) || !base::dir.exists(root)) return(NA_character_)
  preferred <- base::file.path(root, c("manifest.tsv", "templates.tsv", "template_manifest.tsv"))
  hit <- preferred[base::file.exists(preferred)]
  if (!base::length(hit)) {
    hit <- base::list.files(root, pattern = "manifest.*[.]tsv$", full.names = TRUE, ignore.case = TRUE)
  }
  if (!base::length(hit)) return(NA_character_)
  base::normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

.dnmb_rebasefinder_read_homology_templates <- function(root = .dnmb_rebasefinder_homology_template_root()) {
  manifest_path <- .dnmb_rebasefinder_homology_manifest_path(root)
  if (base::is.na(manifest_path)) return(base::data.frame())
  manifest <- tryCatch(
    utils::read.delim(
      manifest_path, stringsAsFactors = FALSE, check.names = FALSE,
      quote = "", comment.char = ""
    ),
    error = function(e) base::data.frame()
  )
  required <- c(
    "template_id", "template_class", "rm_type", "enzyme_role",
    "pdb_file", "fasta_file"
  )
  if (!base::nrow(manifest) || !base::all(required %in% base::names(manifest))) {
    return(base::data.frame())
  }
  resolve <- function(path) {
    path <- base::as.character(path)
    absolute <- base::grepl("^/|^[A-Za-z]:[/\\\\]", path)
    path[!absolute] <- base::file.path(root, path[!absolute])
    base::normalizePath(path, winslash = "/", mustWork = FALSE)
  }
  manifest$pdb_path <- resolve(manifest$pdb_file)
  manifest$fasta_path <- resolve(manifest$fasta_file)
  manifest$template_class <- base::tolower(base::trimws(manifest$template_class))
  manifest$family_description <- if ("family" %in% base::names(manifest)) {
    base::as.character(manifest$family)
  } else {
    NA_character_
  }
  # rm_type is the canonical compatibility key; the manifest's family field is
  # a descriptive fold/subfamily label (for example, C5 DNA MTase).
  manifest$family <- base::as.character(manifest$rm_type)
  manifest$manifest_path <- manifest_path
  manifest <- manifest[
    manifest$template_class %in% c("positive", "decoy") &
      base::file.exists(manifest$pdb_path) & base::file.exists(manifest$fasta_path),
    , drop = FALSE
  ]
  base::rownames(manifest) <- NULL
  manifest
}

.dnmb_rebasefinder_read_single_fasta <- function(path) {
  if (base::is.na(path) || !base::file.exists(path)) return(NA_character_)
  lines <- tryCatch(base::readLines(path, warn = FALSE), error = function(e) character())
  lines <- lines[!base::grepl("^>", lines) & base::nzchar(base::trimws(lines))]
  .dnmb_rebasefinder_normalize_protein(base::paste(lines, collapse = ""))
}

.dnmb_rebasefinder_template_alignment <- function(query_sequence, template_sequence) {
  query_sequence <- .dnmb_rebasefinder_normalize_protein(query_sequence)
  template_sequence <- .dnmb_rebasefinder_normalize_protein(template_sequence)
  if (base::is.na(query_sequence) || base::is.na(template_sequence)) return(NULL)
  alignment <- tryCatch(
    base::suppressWarnings(Biostrings::pairwiseAlignment(
      Biostrings::AAString(query_sequence),
      Biostrings::AAString(template_sequence),
      type = "global"
    )),
    error = function(e) NULL
  )
  if (base::is.null(alignment)) return(NULL)
  aligned_query_text <- base::as.character(Biostrings::alignedPattern(alignment))
  aligned_template_text <- base::as.character(Biostrings::alignedSubject(alignment))
  aligned_query <- base::strsplit(aligned_query_text, "", fixed = TRUE)[[1]]
  aligned_template <- base::strsplit(aligned_template_text, "", fixed = TRUE)[[1]]
  paired <- aligned_query != "-" & aligned_template != "-"
  if (!base::any(paired)) return(NULL)
  query_start <- Biostrings::start(Biostrings::pattern(alignment))[[1]]
  query_end <- Biostrings::end(Biostrings::pattern(alignment))[[1]]
  template_start <- Biostrings::start(Biostrings::subject(alignment))[[1]]
  template_end <- Biostrings::end(Biostrings::subject(alignment))[[1]]
  query_position <- query_start - 1L
  template_position <- template_start - 1L
  map_query <- integer()
  map_template <- integer()
  for (i in base::seq_along(aligned_query)) {
    if (aligned_query[[i]] != "-") query_position <- query_position + 1L
    if (aligned_template[[i]] != "-") template_position <- template_position + 1L
    if (aligned_query[[i]] != "-" && aligned_template[[i]] != "-") {
      map_query <- c(map_query, query_position)
      map_template <- c(map_template, template_position)
    }
  }
  identical_residue <- aligned_query[paired] == aligned_template[paired]
  aligned_length <- base::sum(paired)
  list(
    aligned_query = aligned_query_text,
    aligned_template = aligned_template_text,
    query_start = base::as.integer(query_start),
    query_end = base::as.integer(query_end),
    template_start = base::as.integer(template_start),
    template_end = base::as.integer(template_end),
    aligned_length = base::as.integer(aligned_length),
    identity = base::mean(identical_residue),
    query_coverage = aligned_length / base::nchar(query_sequence),
    query_span_coverage = aligned_length / (query_end - query_start + 1L),
    template_coverage = aligned_length / base::nchar(template_sequence),
    score = base::as.numeric(Biostrings::score(alignment)),
    mapping = base::data.frame(
      query_pos = map_query,
      template_pos = map_template,
      stringsAsFactors = FALSE
    )
  )
}

.dnmb_rebasefinder_homology_motif_hits <- function(query, sequence, family, role) {
  sequence <- .dnmb_rebasefinder_normalize_protein(sequence)
  if (base::is.na(sequence)) return(base::data.frame())
  catalog <- .dnmb_rebasefinder_motif_catalog()
  rows <- list()
  n <- 0L
  for (motif in base::names(catalog)) {
    definition <- catalog[[motif]]
    pattern <- definition$pattern
    if (base::is.null(pattern) || base::is.na(pattern) || !base::nzchar(pattern)) next
    expected_role <- definition$expected_role %||% NA_character_
    role_ok <- base::is.na(role) || !base::nzchar(role) ||
      role %in% c(expected_role, "RM") || expected_role == "RM"
    expected_family <- definition$expected_family %||% NA_character_
    family_ok <- if (base::is.na(expected_family) || !base::nzchar(expected_family)) {
      TRUE
    } else if (base::is.na(family) || !base::nzchar(family)) {
      FALSE
    } else {
      base::identical(
        .dnmb_rebasefinder_canonical_structure_family(expected_family),
        .dnmb_rebasefinder_canonical_structure_family(family)
      )
    }
    if (!role_ok || !family_ok) next
    matches <- base::gregexpr(pattern, sequence, perl = TRUE)[[1]]
    if (!base::length(matches) || matches[[1]] < 1L) next
    lengths <- base::attr(matches, "match.length")
    for (i in base::seq_along(matches)) {
      start <- base::as.integer(matches[[i]])
      end <- start + base::as.integer(lengths[[i]]) - 1L
      range <- definition$pos_range
      if (!base::is.null(range) && base::length(range) == 2L) {
        fraction <- start / base::nchar(sequence)
        if (fraction < range[[1]] || fraction > range[[2]]) next
      }
      n <- n + 1L
      rows[[n]] <- base::data.frame(
        locus_tag = query,
        family_id = family,
        enzyme_role = role,
        motif = motif,
        start_aa = start,
        end_aa = end,
        matched_sequence = base::substr(sequence, start, end),
        evidence_level = "sequence_hint",
        stringsAsFactors = FALSE
      )
    }
  }
  if (!base::length(rows)) return(base::data.frame())
  base::do.call(base::rbind, rows)
}

.dnmb_rebasefinder_homology_motif_mapping <- function(motif_hits, alignment, flank = 3L) {
  if (base::is.null(alignment) || !base::nrow(motif_hits)) {
    return(list(complete = FALSE, pairs = NA_character_, mapped_motifs = NA_character_))
  }
  mapped <- base::vapply(base::seq_len(base::nrow(motif_hits)), function(i) {
    start <- base::as.integer(motif_hits$start_aa[[i]])
    end <- base::as.integer(motif_hits$end_aa[[i]])
    required <- base::seq.int(base::max(1L, start - flank), end + flank)
    base::all(required %in% alignment$mapping$query_pos)
  }, logical(1))
  mapped_motifs <- base::unique(base::as.character(motif_hits$motif[mapped]))
  rules <- .dnmb_rebasefinder_structure_pair_rules()
  family <- base::as.character(motif_hits$family_id[[1]])
  rule_family <- .dnmb_rebasefinder_canonical_structure_family(family)
  role <- base::as.character(motif_hits$enzyme_role[[1]])
  applicable <- rules[
    base::vapply(base::seq_len(base::nrow(rules)), function(i) {
      family_ok <- base::is.na(rule_family) || !base::nzchar(rule_family) ||
        base::grepl(rules$family_pattern[[i]], rule_family, ignore.case = TRUE)
      role_ok <- base::is.na(role) || !base::nzchar(role) ||
        role %in% base::trimws(base::strsplit(rules$roles[[i]], ",", fixed = TRUE)[[1]])
      family_ok && role_ok
    }, logical(1)),
    , drop = FALSE
  ]
  complete_pairs <- if (base::nrow(applicable)) {
    applicable$pair_id[
      applicable$motif_a %in% mapped_motifs & applicable$motif_b %in% mapped_motifs
    ]
  } else {
    character()
  }
  list(
    complete = base::length(complete_pairs) > 0L,
    pairs = if (base::length(complete_pairs)) base::paste(base::unique(complete_pairs), collapse = ",") else NA_character_,
    mapped_motifs = if (base::length(mapped_motifs)) base::paste(mapped_motifs, collapse = ",") else NA_character_
  )
}

.dnmb_rebasefinder_homology_candidate_keep <- function(hits, max_candidates = 24L) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  if (!base::nrow(hits)) return(logical())
  max_candidates <- base::suppressWarnings(base::as.integer(max_candidates)[[1]])
  if (base::is.na(max_candidates) || max_candidates < 1L) max_candidates <- 24L
  tier <- if ("curation_tier" %in% base::names(hits)) {
    base::as.character(hits$curation_tier)
  } else {
    base::rep("review", base::nrow(hits))
  }
  axes <- if ("independent_evidence_axes" %in% base::names(hits)) {
    base::suppressWarnings(base::as.integer(hits$independent_evidence_axes))
  } else {
    base::rep(0L, base::nrow(hits))
  }
  axes[base::is.na(axes)] <- 0L
  role <- if ("final_role" %in% base::names(hits)) {
    base::as.character(hits$final_role)
  } else if ("enzyme_role" %in% base::names(hits)) {
    base::as.character(hits$enzyme_role)
  } else {
    base::rep(NA_character_, base::nrow(hits))
  }
  keep <- !tier %in% c("excluded_noise", "other_defense") &
    role %in% c("M", "R", "RM") &
    (tier %in% c("high", "medium") | axes >= 2L)
  keep[base::is.na(keep)] <- FALSE
  selected <- base::which(keep)
  if (base::length(selected) > max_candidates) {
    tier_rank <- base::match(tier[selected], c("high", "medium", "review"))
    score <- if ("curation_score" %in% base::names(hits)) {
      base::suppressWarnings(base::as.numeric(hits$curation_score[selected]))
    } else {
      base::rep(0, base::length(selected))
    }
    score[base::is.na(score)] <- 0
    selected <- selected[base::order(tier_rank, -score, -axes[selected], selected)]
    restriction <- selected[role[selected] %in% c("R", "RM")]
    restriction_quota <- base::min(
      base::length(restriction),
      base::max(1L, base::ceiling(max_candidates * 0.25))
    )
    priority <- if (restriction_quota > 0L) {
      restriction[base::seq_len(restriction_quota)]
    } else {
      integer()
    }
    selected <- c(priority, selected[!selected %in% priority])
    keep[] <- FALSE
    keep[selected[base::seq_len(max_candidates)]] <- TRUE
  }
  keep
}

.dnmb_rebasefinder_hash_text <- function(text) {
  path <- base::tempfile(fileext = ".txt")
  on.exit(base::unlink(path, force = TRUE), add = TRUE)
  base::writeLines(base::enc2utf8(base::as.character(text)), path, useBytes = TRUE)
  base::unname(tools::md5sum(path)[[1]])
}

.dnmb_rebasefinder_promod3_layout <- function(cache_root = NULL, create = TRUE) {
  root <- base::file.path(
    .dnmb_db_cache_root(cache_root = cache_root, create = create),
    "rebasefinder", "promod3", "3.6.0"
  )
  if (base::isTRUE(create)) base::dir.create(root, recursive = TRUE, showWarnings = FALSE)
  list(
    root = root,
    env = base::file.path(root, "env"),
    cli = base::file.path(root, "env", "bin", "pm"),
    ready = base::file.path(root, ".install-complete"),
    lock = base::file.path(root, ".install-lock"),
    models = base::file.path(root, "models"),
    packages = base::file.path(root, "pkgs")
  )
}

.dnmb_rebasefinder_promod3_cli_usable <- function(path, timeout = 30L) {
  path <- base::as.character(path)[[1]]
  if (base::is.na(path) || !base::nzchar(path) || !base::file.exists(path) ||
      base::file.access(path, mode = 1L) != 0L) {
    return(FALSE)
  }
  output <- tryCatch(
    base::suppressWarnings(base::system2(
      path,
      args = "help",
      stdout = TRUE,
      stderr = TRUE,
      timeout = base::as.integer(timeout)
    )),
    error = function(e) base::structure(character(), status = 127L)
  )
  status <- base::attr(output, "status")
  if (base::is.null(status)) status <- 0L
  base::identical(base::as.integer(status), 0L) &&
    base::any(base::grepl("build-model", output, fixed = TRUE))
}

.dnmb_rebasefinder_promod3_ready <- function(layout) {
  !base::is.null(layout) && base::file.exists(layout$ready) &&
    base::file.exists(layout$cli) && base::file.access(layout$cli, mode = 1L) == 0L
}

.dnmb_rebasefinder_mark_promod3_ready <- function(layout) {
  base::dir.create(layout$root, recursive = TRUE, showWarnings = FALSE)
  marker <- base::tempfile(".install-complete-", tmpdir = layout$root)
  on.exit(base::unlink(marker, force = TRUE), add = TRUE)
  base::writeLines(
    c(
      "promod3=3.6.0",
      base::paste0("completed_at=", base::format(base::Sys.time(), tz = "UTC", usetz = TRUE))
    ),
    marker,
    useBytes = TRUE
  )
  base::unlink(layout$ready, force = TRUE)
  base::isTRUE(base::file.rename(marker, layout$ready)) && base::file.exists(layout$ready)
}

.dnmb_rebasefinder_candidate_conda <- function() {
  candidates <- c(
    base::Sys.which("micromamba"), base::Sys.which("mamba"), base::Sys.which("conda"),
    base::file.path(base::path.expand("~"), "miniforge3", "bin", "conda"),
    base::file.path(base::path.expand("~"), "mambaforge", "bin", "mamba"),
    "/opt/homebrew/Caskroom/miniforge/base/bin/conda",
    "/opt/conda/bin/conda"
  )
  candidates <- candidates[base::nzchar(candidates) & base::file.exists(candidates)]
  if (!base::length(candidates)) return(character())
  base::unique(base::normalizePath(candidates, winslash = "/", mustWork = TRUE))
}

.dnmb_rebasefinder_free_disk_bytes <- function(path) {
  if (.Platform$OS.type == "windows" || !base::nzchar(base::Sys.which("df"))) return(NA_real_)
  base::dir.create(path, recursive = TRUE, showWarnings = FALSE)
  output <- tryCatch(
    base::system2("df", args = c("-Pk", base::shQuote(path)), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (base::length(output) < 2L) return(NA_real_)
  fields <- base::strsplit(base::trimws(output[[base::length(output)]]), "[[:space:]]+")[[1]]
  if (base::length(fields) < 4L) return(NA_real_)
  available_kb <- base::suppressWarnings(base::as.numeric(fields[[4]]))
  if (!base::is.finite(available_kb)) NA_real_ else available_kb * 1024
}

.dnmb_rebasefinder_promod3_cli <- function(cache_root = NULL,
                                           install = TRUE,
                                           verbose = TRUE) {
  path_cli <- base::Sys.which("pm")
  if (base::nzchar(path_cli)) {
    return(list(path = path_cli, status = "available", detail = "pm found on PATH"))
  }
  layout <- .dnmb_rebasefinder_promod3_layout(cache_root = cache_root, create = TRUE)
  if (.dnmb_rebasefinder_promod3_ready(layout)) {
    return(list(path = layout$cli, status = "available", detail = "managed ProMod3 environment"))
  }
  if (.Platform$OS.type == "windows") {
    return(list(path = "", status = "backend_unavailable", detail = "ProMod3 has no native Windows build; use WSL2 or Docker"))
  }
  lock_dir <- layout$lock
  lock_deadline <- base::Sys.time() + 1800
  waiting_announced <- FALSE
  repeat {
    if (.dnmb_rebasefinder_promod3_ready(layout)) {
      return(list(path = layout$cli, status = "available", detail = "managed ProMod3 environment installed by another process"))
    }
    lock_acquired <- base::dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)
    if (base::isTRUE(lock_acquired)) break
    if (!base::isTRUE(install)) {
      return(list(
        path = "", status = "backend_unavailable",
        detail = "managed ProMod3 installation is incomplete or in progress and module_install=FALSE"
      ))
    }
    lock_info <- base::file.info(lock_dir)
    lock_age <- base::as.numeric(base::difftime(base::Sys.time(), lock_info$mtime[[1]], units = "secs"))
    if (base::is.finite(lock_age) && lock_age > 7200) {
      base::unlink(lock_dir, recursive = TRUE, force = TRUE)
      next
    }
    if (base::Sys.time() >= lock_deadline) {
      return(list(
        path = "", status = "backend_unavailable",
        detail = "timed out waiting for another ProMod3 installation"
      ))
    }
    if (base::isTRUE(verbose) && !waiting_announced) {
      base::message("[REBASEfinder] Waiting for another ProMod3 installation...")
      waiting_announced <- TRUE
    }
    base::Sys.sleep(2)
  }
  base::on.exit(base::unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  if (.dnmb_rebasefinder_promod3_ready(layout)) {
    return(list(path = layout$cli, status = "available", detail = "managed ProMod3 environment"))
  }
  if (base::file.exists(layout$cli)) {
    if (.dnmb_rebasefinder_promod3_cli_usable(layout$cli) &&
        .dnmb_rebasefinder_mark_promod3_ready(layout)) {
      return(list(
        path = layout$cli, status = "available",
        detail = "validated legacy managed ProMod3 environment"
      ))
    }
    base::unlink(c(layout$env, layout$packages, layout$ready), recursive = TRUE, force = TRUE)
  }
  if (!base::isTRUE(install)) {
    return(list(path = "", status = "backend_unavailable", detail = "pm not installed and module_install=FALSE"))
  }
  free_bytes <- .dnmb_rebasefinder_free_disk_bytes(layout$root)
  minimum_bytes <- 5 * 1024^3
  if (base::is.finite(free_bytes) && free_bytes < minimum_bytes) {
    return(list(
      path = "", status = "backend_unavailable",
      detail = base::sprintf(
        "insufficient disk space for managed ProMod3 install (%.2f GiB free; 5.0 GiB required)",
        free_bytes / 1024^3
      )
    ))
  }
  conda_candidates <- .dnmb_rebasefinder_candidate_conda()
  if (!base::length(conda_candidates)) {
    return(list(path = "", status = "backend_unavailable", detail = "pm and conda/mamba were not found"))
  }
  base::dir.create(layout$packages, recursive = TRUE, showWarnings = FALSE)
  if (base::isTRUE(verbose)) {
    base::message("[REBASEfinder] Installing the local ProMod3 backend (first run only)...")
  }
  failures <- character()
  for (conda in conda_candidates) {
    base::unlink(c(layout$env, layout$ready), recursive = TRUE, force = TRUE)
    run <- dnmb_run_external(
      conda,
      args = c(
        "create", "-y", "-p", layout$env,
        "-c", "conda-forge", "-c", "bioconda", "promod3=3.6.0"
      ),
      env = c(CONDA_PKGS_DIRS = layout$packages),
      required = FALSE
    )
    if (base::isTRUE(run$ok) && .dnmb_rebasefinder_promod3_cli_usable(layout$cli) &&
        .dnmb_rebasefinder_mark_promod3_ready(layout)) {
      base::unlink(layout$packages, recursive = TRUE, force = TRUE)
      return(list(
        path = layout$cli, status = "available",
        detail = base::sprintf("managed ProMod3 environment installed with %s", base::basename(conda))
      ))
    }
    failure <- run$error %||% "managed ProMod3 installation failed"
    failure <- base::substr(base::gsub("[\r\n]+", " ", failure), 1L, 1200L)
    failures <- c(failures, base::sprintf("%s: %s", conda, failure))
  }
  base::unlink(c(layout$env, layout$packages, layout$ready), recursive = TRUE, force = TRUE)
  list(
    path = "", status = "backend_unavailable",
    detail = base::paste(failures, collapse = " | ")
  )
}

.dnmb_rebasefinder_promod3_wrapper_path <- function() {
  installed <- base::system.file("scripts", "rebasefinder_promod3_model.R", package = "DNMB")
  source_file <- tryCatch(
    utils::getSrcFilename(.dnmb_rebasefinder_promod3_wrapper_path, full.names = TRUE),
    error = function(e) ""
  )
  source_root <- if (base::length(source_file) && base::nzchar(source_file[[1]])) {
    base::dirname(base::dirname(source_file[[1]]))
  } else {
    ""
  }
  candidates <- base::unique(c(
    installed,
    if (base::nzchar(source_root)) base::file.path(source_root, "inst", "scripts", "rebasefinder_promod3_model.R") else "",
    base::file.path(base::getwd(), "inst", "scripts", "rebasefinder_promod3_model.R"),
    base::file.path(base::dirname(base::getwd()), "inst", "scripts", "rebasefinder_promod3_model.R")
  ))
  hit <- candidates[base::nzchar(candidates) & base::file.exists(candidates)]
  if (!base::length(hit)) return(NA_character_)
  base::normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

.dnmb_rebasefinder_prepare_template_pdb <- function(path, destination) {
  if (base::is.na(path) || !base::file.exists(path)) return(NA_character_)
  base::dir.create(base::dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (base::grepl("[.]gz$", path, ignore.case = TRUE)) {
    input <- base::gzfile(path, open = "rt")
    output <- base::file(destination, open = "wt")
    ok <- tryCatch({
      repeat {
        lines <- base::readLines(input, n = 10000L, warn = FALSE)
        if (!base::length(lines)) break
        base::writeLines(lines, output)
      }
      TRUE
    }, error = function(e) FALSE)
    try(base::close(input), silent = TRUE)
    try(base::close(output), silent = TRUE)
    if (!ok) return(NA_character_)
  } else if (!base::file.copy(path, destination, overwrite = TRUE)) {
    return(NA_character_)
  }
  if (!base::file.exists(destination)) NA_character_ else destination
}

.dnmb_rebasefinder_write_promod3_alignment <- function(path, alignment) {
  base::dir.create(base::dirname(path), recursive = TRUE, showWarnings = FALSE)
  base::writeLines(c(
    ">target", alignment$aligned_query,
    ">template", alignment$aligned_template
  ), path)
  path
}

.dnmb_rebasefinder_write_homology_audit <- function(output_dir, table) {
  tsv <- base::file.path(output_dir, "DNMB_REBASEfinder_homology_templates.tsv")
  utils::write.table(table, tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  xlsx <- base::file.path(output_dir, "DNMB_REBASEfinder_homology_templates.xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    tryCatch(openxlsx::write.xlsx(table, xlsx, overwrite = TRUE), error = function(e) NULL)
  }
  list(
    tsv = base::normalizePath(tsv, winslash = "/", mustWork = FALSE),
    xlsx = if (base::file.exists(xlsx)) base::normalizePath(xlsx, winslash = "/", mustWork = FALSE) else NA_character_
  )
}

.dnmb_rebasefinder_merge_homology <- function(hits, audit) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  audit <- base::as.data.frame(audit, stringsAsFactors = FALSE)
  if (!base::nrow(hits) || !base::nrow(audit) || !"query" %in% base::names(audit)) return(hits)
  fields <- base::setdiff(base::names(audit), c(
    "family_id", "enzyme_role", "curation_tier", "curation_score",
    "aligned_query", "aligned_template", "template_pdb_path", "template_fasta_path"
  ))
  fields <- base::setdiff(fields, "query")
  idx <- base::match(.dnmb_module_clean_annotation_key(hits$query), .dnmb_module_clean_annotation_key(audit$query))
  ok <- !base::is.na(idx)
  for (field in fields) {
    target <- base::paste0("homology_", field)
    prototype <- audit[[field]]
    if (base::is.logical(prototype)) {
      value <- base::rep(NA, base::nrow(hits))
    } else if (base::is.numeric(prototype) || base::is.integer(prototype)) {
      value <- base::rep(NA_real_, base::nrow(hits))
    } else {
      value <- base::rep(NA_character_, base::nrow(hits))
    }
    value[ok] <- prototype[idx[ok]]
    hits[[target]] <- value
  }
  hits$homology_model_supported <- hits$homology_model_supported %in% TRUE
  support_idx <- base::which(hits$homology_model_supported)
  if (base::length(support_idx) && "support" %in% base::names(hits)) {
    for (i in support_idx) {
      note <- base::sprintf(
        "ProMod3 homology model: %s; identity=%.3f; template_cov=%.3f; geometry=%s",
        hits$homology_template_id[[i]], hits$homology_alignment_identity[[i]],
        hits$homology_template_coverage[[i]], hits$homology_geometry_status[[i]]
      )
      old <- base::as.character(hits$support[[i]])
      hits$support[[i]] <- if (base::is.na(old) || !base::nzchar(old)) note else base::paste(old, note, sep = "; ")
    }
  }
  hits
}

.dnmb_rebasefinder_clear_promod3_outputs <- function(output_dir) {
  model_dir <- base::file.path(output_dir, "promod3_query_structures")
  base::dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  stale <- base::list.files(model_dir, pattern = "[.]pdb$", full.names = TRUE, ignore.case = TRUE)
  if (base::length(stale)) base::unlink(stale, force = TRUE)
  base::unlink(
    c(
      base::file.path(output_dir, "DNMB_REBASEfinder_homology_templates.tsv"),
      base::file.path(output_dir, "DNMB_REBASEfinder_homology_templates.xlsx"),
      base::file.path(output_dir, "promod3_work")
    ),
    recursive = TRUE,
    force = TRUE
  )
  model_dir
}

.dnmb_rebasefinder_run_homology_models <- function(hits,
                                                   genes,
                                                   output_dir,
                                                   cache_root = NULL,
                                                   install = TRUE,
                                                   cpu = 1L,
                                                   max_candidates = 24L,
                                                   min_identity = 0.30,
                                                   min_template_coverage = 0.60,
                                                   min_aligned_length = 80L,
                                                   min_decoy_margin = 0.02,
                                                   min_decoy_relative_margin = 0.15,
                                                   verbose = TRUE) {
  hits <- base::as.data.frame(hits, stringsAsFactors = FALSE)
  genes <- base::as.data.frame(genes, stringsAsFactors = FALSE)
  empty <- list(
    hits = hits, audit = base::data.frame(), status = "not_run",
    detail = "No candidates", files = list()
  )
  output_model_dir <- .dnmb_rebasefinder_clear_promod3_outputs(output_dir)
  templates <- .dnmb_rebasefinder_read_homology_templates()
  if (!base::nrow(hits) || !base::nrow(genes) || !base::nrow(templates)) {
    empty$detail <- if (!base::nrow(templates)) "Bundled R-M template manifest is unavailable" else "No REBASEfinder candidates"
    return(empty)
  }
  genes$locus_tag <- .dnmb_module_clean_annotation_key(genes$locus_tag)
  selected <- .dnmb_rebasefinder_homology_candidate_keep(hits, max_candidates = max_candidates)
  candidate_hits <- hits[selected, , drop = FALSE]
  if (!base::nrow(candidate_hits)) return(empty)

  template_sequences <- base::vapply(templates$fasta_path, .dnmb_rebasefinder_read_single_fasta, character(1))
  template_ok <- !base::is.na(template_sequences) & base::nzchar(template_sequences)
  templates <- templates[template_ok, , drop = FALSE]
  template_sequences <- template_sequences[template_ok]
  if (!base::nrow(templates)) {
    empty$detail <- "Bundled template FASTA files are unreadable"
    return(empty)
  }

  audit_rows <- list()
  alignment_objects <- list()
  motif_tables <- list()
  for (candidate_index in base::seq_len(base::nrow(candidate_hits))) {
    candidate <- candidate_hits[candidate_index, , drop = FALSE]
    query <- .dnmb_module_clean_annotation_key(candidate$query[[1]])
    gene_index <- base::match(query, genes$locus_tag)
    if (base::is.na(gene_index) || !"translation" %in% base::names(genes)) next
    sequence <- .dnmb_rebasefinder_normalize_protein(genes$translation[[gene_index]])
    if (base::is.na(sequence)) next
    family <- if ("final_family" %in% base::names(candidate)) candidate$final_family[[1]] else candidate$family_id[[1]]
    role <- if ("final_role" %in% base::names(candidate)) candidate$final_role[[1]] else candidate$enzyme_role[[1]]
    family_unspecified <- base::is.na(family) || !base::nzchar(base::trimws(family)) ||
      base::grepl("unclassified|unknown|unspecified|undetermined", family, ignore.case = TRUE)
    one_rows <- list()
    one_alignments <- list()
    one_motif_tables <- list()
    for (template_index in base::seq_len(base::nrow(templates))) {
      template <- templates[template_index, , drop = FALSE]
      family_ok <- if (template$template_class[[1]] == "decoy") {
        TRUE
      } else {
        family_unspecified ||
          .dnmb_rebasefinder_structure_family_compatible(family, template$family[[1]])[[1]]
      }
      role_ok <- if (template$template_class[[1]] == "decoy") {
        TRUE
      } else {
        .dnmb_rebasefinder_structure_role_compatible(role, template$enzyme_role[[1]])[[1]]
      }
      # Incompatible positive templates cannot win and are skipped to keep the
      # local search fast. Decoys are always aligned for explicit noise ranking.
      if (template$template_class[[1]] == "positive" && (!family_ok || !role_ok)) next
      alignment <- .dnmb_rebasefinder_template_alignment(sequence, template_sequences[[template_index]])
      if (base::is.null(alignment)) next
      motif_family <- if (family_unspecified && template$template_class[[1]] == "positive") {
        template$family[[1]]
      } else {
        family
      }
      template_motif_hits <- .dnmb_rebasefinder_homology_motif_hits(
        query, sequence, motif_family, role
      )
      motif_mapping <- .dnmb_rebasefinder_homology_motif_mapping(template_motif_hits, alignment)
      rank_score <- alignment$identity * alignment$template_coverage * alignment$query_span_coverage
      one_rows[[base::length(one_rows) + 1L]] <- base::data.frame(
        query = query,
        family_id = family,
        enzyme_role = role,
        curation_tier = candidate$curation_tier[[1]] %||% NA_character_,
        curation_score = candidate$curation_score[[1]] %||% NA_real_,
        template_id = template$template_id[[1]],
        template_class = template$template_class[[1]],
        template_family = template$family[[1]],
        template_role = template$enzyme_role[[1]],
        template_description = if ("description" %in% base::names(template)) template$description[[1]] else NA_character_,
        motif_family_used = motif_family,
        family_compatible = family_ok,
        role_compatible = role_ok,
        alignment_identity = alignment$identity,
        query_coverage = alignment$query_coverage,
        query_span_coverage = alignment$query_span_coverage,
        template_coverage = alignment$template_coverage,
        aligned_length = alignment$aligned_length,
        query_start = alignment$query_start,
        query_end = alignment$query_end,
        template_start = alignment$template_start,
        template_end = alignment$template_end,
        alignment_score = alignment$score,
        rank_score = rank_score,
        mapped_motifs = motif_mapping$mapped_motifs,
        mapped_motif_pairs = motif_mapping$pairs,
        motif_mapping_complete = motif_mapping$complete,
        template_pdb_path = template$pdb_path[[1]],
        template_fasta_path = template$fasta_path[[1]],
        stringsAsFactors = FALSE
      )
      one_alignments[[template$template_id[[1]]]] <- alignment
      one_motif_tables[[template$template_id[[1]]]] <- template_motif_hits
    }
    if (!base::length(one_rows)) next
    one <- base::do.call(base::rbind, one_rows)
    positive <- one[
      one$template_class == "positive" & one$family_compatible & one$role_compatible,
      , drop = FALSE
    ]
    decoy <- one[one$template_class == "decoy", , drop = FALSE]
    positive <- positive[base::order(-positive$rank_score, -positive$alignment_score), , drop = FALSE]
    decoy <- decoy[base::order(-decoy$rank_score, -decoy$alignment_score), , drop = FALSE]
    if (!base::nrow(positive)) {
      best <- one[base::which.max(one$rank_score), , drop = FALSE]
      best$selection_status <- "no_compatible_positive_template"
    } else {
      threshold_pass <- positive$alignment_identity >= min_identity &
        positive$template_coverage >= min_template_coverage &
        positive$aligned_length >= min_aligned_length
      threshold_pass[base::is.na(threshold_pass)] <- FALSE
      best <- if (base::any(threshold_pass)) {
        positive[base::which(threshold_pass)[[1]], , drop = FALSE]
      } else {
        positive[1, , drop = FALSE]
      }
      best_decoy_score <- if (base::nrow(decoy)) decoy$rank_score[[1]] else NA_real_
      best$best_decoy_id <- if (base::nrow(decoy)) decoy$template_id[[1]] else NA_character_
      best$best_decoy_score <- best_decoy_score
      best$decoy_margin <- if (base::is.finite(best_decoy_score)) best$rank_score - best_decoy_score else NA_real_
      best$decoy_relative_margin <- if (base::is.finite(best_decoy_score)) {
        (best$rank_score - best_decoy_score) / base::max(best_decoy_score, 0.05)
      } else {
        NA_real_
      }
      best$selection_status <- if (best$alignment_identity < min_identity) {
        "below_identity_threshold"
      } else if (best$template_coverage < min_template_coverage) {
        "below_template_coverage_threshold"
      } else if (best$aligned_length < min_aligned_length) {
        "alignment_too_short"
      } else if (base::is.finite(best_decoy_score) &&
                 (best$decoy_margin < min_decoy_margin ||
                  best$decoy_relative_margin < min_decoy_relative_margin)) {
        "decoy_not_beaten"
      } else {
        "template_mapped"
      }
    }
    for (column in c("best_decoy_id", "best_decoy_score", "decoy_margin", "decoy_relative_margin")) {
      if (!column %in% base::names(best)) best[[column]] <- if (column == "best_decoy_id") NA_character_ else NA_real_
    }
    best$model_eligible <- best$selection_status == "template_mapped"
    best$model_status <- if (best$model_eligible) "pending" else "not_eligible"
    best$model_path <- NA_character_
    best$model_expected_path <- NA_character_
    best$model_cache_hit <- FALSE
    best$model_supported <- FALSE
    best$geometry_status <- NA_character_
    audit_rows[[base::length(audit_rows) + 1L]] <- best
    alignment_objects[[query]] <- one_alignments[[best$template_id[[1]]]]
    motif_tables[[query]] <- one_motif_tables[[best$template_id[[1]]]] %||% base::data.frame()
  }
  if (!base::length(audit_rows)) return(empty)
  audit <- base::do.call(base::rbind, audit_rows)
  eligible <- audit$model_eligible %in% TRUE
  work_dir <- base::file.path(output_dir, "promod3_work")
  audit_alignment_dir <- base::file.path(work_dir, "alignments")
  audit_template_dir <- base::file.path(work_dir, "templates")
  base::dir.create(audit_alignment_dir, recursive = TRUE, showWarnings = FALSE)
  base::dir.create(audit_template_dir, recursive = TRUE, showWarnings = FALSE)
  jobs <- list()
  job_audit_rows <- list()
  job_cache_keys <- character()
  layout <- .dnmb_rebasefinder_promod3_layout(cache_root = cache_root, create = TRUE)
  base::dir.create(layout$models, recursive = TRUE, showWarnings = FALSE)
  temp_candidates <- base::unique(c(base::tempdir(), base::Sys.getenv("TMPDIR", unset = ""), "/tmp"))
  temp_candidates <- temp_candidates[
    base::nzchar(temp_candidates) & base::dir.exists(temp_candidates) &
      !base::grepl("[[:space:]]", temp_candidates) & base::file.access(temp_candidates, 2L) == 0L
  ]
  if (!base::length(temp_candidates)) {
    audit$model_status[eligible] <- "backend_unavailable:no_space_free_tempdir"
    files <- .dnmb_rebasefinder_write_homology_audit(output_dir, audit)
    return(list(
      hits = .dnmb_rebasefinder_merge_homology(hits, audit), audit = audit,
      status = "partial", detail = "No writable whitespace-free temporary directory for ProMod3",
      files = files
    ))
  }
  stage_dir <- base::file.path(
    temp_candidates[[1]], "dnmb-promod3",
    base::paste0(.dnmb_rebasefinder_hash_text(base::normalizePath(output_dir, winslash = "/", mustWork = FALSE)), "-", base::Sys.getpid())
  )
  base::unlink(stage_dir, recursive = TRUE, force = TRUE)
  base::dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  base::on.exit(base::unlink(stage_dir, recursive = TRUE, force = TRUE), add = TRUE)
  alignment_dir <- base::file.path(stage_dir, "alignments")
  template_dir <- base::file.path(stage_dir, "templates")
  staged_model_dir <- base::file.path(stage_dir, "models")
  base::dir.create(alignment_dir, recursive = TRUE, showWarnings = FALSE)
  base::dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)
  base::dir.create(staged_model_dir, recursive = TRUE, showWarnings = FALSE)
  wrapper <- .dnmb_rebasefinder_promod3_wrapper_path()
  wrapper_hash <- if (!base::is.na(wrapper) && base::file.exists(wrapper)) {
    base::unname(tools::md5sum(wrapper)[[1]])
  } else {
    "wrapper-unavailable"
  }
  safe_queries <- base::gsub("[^A-Za-z0-9_.-]+", "_", audit$query)
  safe_collision <- base::duplicated(safe_queries) | base::duplicated(safe_queries, fromLast = TRUE)
  if (base::any(safe_collision)) {
    collision_rows <- base::which(safe_collision)
    safe_queries[collision_rows] <- base::vapply(collision_rows, function(i) {
      base::paste0(
        safe_queries[[i]], "__",
        base::substr(.dnmb_rebasefinder_hash_text(audit$query[[i]]), 1L, 8L), "_", i
      )
    }, character(1))
  }
  for (i in base::which(eligible)) {
    safe_query <- safe_queries[[i]]
    file_stem <- base::paste0(safe_query, "__", audit$template_id[[i]])
    audit_alignment_path <- base::file.path(audit_alignment_dir, base::paste0(file_stem, ".fasta"))
    audit_template_path <- base::file.path(audit_template_dir, base::paste0(audit$template_id[[i]], ".pdb"))
    alignment_path <- base::file.path(alignment_dir, base::paste0(file_stem, ".fasta"))
    template_path <- base::file.path(template_dir, base::paste0(audit$template_id[[i]], ".pdb"))
    .dnmb_rebasefinder_write_promod3_alignment(audit_alignment_path, alignment_objects[[audit$query[[i]]]])
    audit_template_path <- .dnmb_rebasefinder_prepare_template_pdb(audit$template_pdb_path[[i]], audit_template_path)
    alignment_ok <- base::file.copy(audit_alignment_path, alignment_path, overwrite = TRUE)
    template_ok <- !base::is.na(audit_template_path) && base::file.copy(audit_template_path, template_path, overwrite = TRUE)
    if (!base::isTRUE(alignment_ok) || !base::isTRUE(template_ok)) template_path <- NA_character_
    if (base::is.na(template_path)) {
      audit$model_status[[i]] <- "failed_template_prepare"
      next
    }
    template_hash <- base::unname(tools::md5sum(audit$template_pdb_path[[i]])[[1]])
    cache_key <- .dnmb_rebasefinder_hash_text(base::paste(
      audit$template_id[[i]], template_hash, wrapper_hash,
      alignment_objects[[audit$query[[i]]]]$aligned_query,
      alignment_objects[[audit$query[[i]]]]$aligned_template,
      "promod3-3.6.0", sep = "\n"
    ))
    cache_model <- base::file.path(layout$models, base::paste0(cache_key, ".pdb"))
    final_model <- base::file.path(output_model_dir, base::paste0(safe_query, ".pdb"))
    audit$model_expected_path[[i]] <- base::normalizePath(final_model, winslash = "/", mustWork = FALSE)
    if (base::file.exists(cache_model) && base::nrow(.dnmb_rebasefinder_read_pdb_ca(cache_model))) {
      copied <- base::file.copy(cache_model, final_model, overwrite = TRUE)
      if (base::isTRUE(copied) && base::nrow(.dnmb_rebasefinder_read_pdb_ca(final_model))) {
        audit$model_status[[i]] <- "template_model_built"
        audit$model_path[[i]] <- base::normalizePath(final_model, winslash = "/", mustWork = FALSE)
        audit$model_cache_hit[[i]] <- TRUE
      } else {
        audit$model_status[[i]] <- "failed_cached_model_copy"
      }
      next
    }
    staged_model <- base::file.path(staged_model_dir, base::paste0(cache_key, ".pdb"))
    existing_job <- base::match(cache_key, job_cache_keys)
    if (!base::is.na(existing_job)) {
      job_audit_rows[[existing_job]] <- c(job_audit_rows[[existing_job]], i)
      next
    }
    jobs[[base::length(jobs) + 1L]] <- base::data.frame(
      query = audit$query[[i]], alignment = alignment_path,
      template = template_path, output = staged_model,
      cache_output = cache_model, final_output = final_model,
      alignment_source = audit_alignment_path, template_source = audit_template_path,
      stringsAsFactors = FALSE
    )
    job_cache_keys <- c(job_cache_keys, cache_key)
    job_audit_rows[[base::length(jobs)]] <- i
  }

  backend <- list(path = "", status = "not_needed", detail = "No eligible uncached model")
  model_manifest <- NA_character_
  if (base::length(jobs)) {
    backend <- .dnmb_rebasefinder_promod3_cli(
      cache_root = cache_root, install = install, verbose = verbose
    )
    pending_rows <- base::unique(base::unlist(job_audit_rows, use.names = FALSE))
    if (!base::nzchar(backend$path)) {
      audit$model_status[pending_rows] <- "backend_unavailable"
    } else {
      jobs_tbl <- base::do.call(base::rbind, jobs)
      jobs_path <- base::file.path(work_dir, "promod3_jobs.tsv")
      model_manifest <- base::file.path(work_dir, "promod3_models.tsv")
      utils::write.table(jobs_tbl, jobs_path, sep = "\t", quote = FALSE, row.names = FALSE)
      if (base::is.na(wrapper)) {
        audit$model_status[pending_rows] <- "wrapper_unavailable"
      } else {
        run <- dnmb_run_external(
          base::file.path(R.home("bin"), "Rscript"),
          args = c(
            wrapper, "--jobs", jobs_path, "--manifest", model_manifest,
            "--pm", backend$path, "--threads", base::as.character(base::max(1L, base::as.integer(cpu))),
            "--timeout", "1800", "--overwrite"
          ),
          required = FALSE
        )
        model_results <- if (base::file.exists(model_manifest)) {
          tryCatch(utils::read.delim(model_manifest, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) base::data.frame())
        } else {
          base::data.frame()
        }
        for (j in base::seq_along(job_audit_rows)) {
          audit_indices <- job_audit_rows[[j]]
          representative <- audit_indices[[1]]
          result_index <- if (base::nrow(model_results)) {
            base::match(audit$query[[representative]], model_results$query)
          } else {
            NA_integer_
          }
          staged_model <- jobs_tbl$output[[j]]
          cache_model <- jobs_tbl$cache_output[[j]]
          valid <- !base::is.na(result_index) &&
            model_results$status[[result_index]] %in% c("ok", "exists") &&
            base::file.exists(staged_model) && base::nrow(.dnmb_rebasefinder_read_pdb_ca(staged_model))
          if (valid) {
            cache_copied <- base::file.copy(staged_model, cache_model, overwrite = TRUE)
            for (i in audit_indices) {
              final_model <- audit$model_expected_path[[i]]
              final_copied <- base::isTRUE(cache_copied) &&
                base::file.copy(cache_model, final_model, overwrite = TRUE)
              if (base::isTRUE(final_copied) && base::nrow(.dnmb_rebasefinder_read_pdb_ca(final_model))) {
                audit$model_status[[i]] <- "template_model_built"
                audit$model_path[[i]] <- base::normalizePath(final_model, winslash = "/", mustWork = FALSE)
              } else {
                audit$model_status[[i]] <- "model_failed:validated_copy"
              }
            }
          } else {
            failure_status <- if (!base::is.na(result_index)) {
              base::paste0("model_failed:", model_results$status[[result_index]])
            } else if (!base::isTRUE(run$ok)) {
              "model_failed"
            } else {
              "model_missing"
            }
            audit$model_status[audit_indices] <- failure_status
          }
        }
      }
    }
  }

  model_rows <- base::which(audit$model_status == "template_model_built" & base::file.exists(audit$model_path))
  for (i in model_rows) {
    motif_hits <- motif_tables[[audit$query[[i]]]]
    if (base::is.null(motif_hits) || !base::nrow(motif_hits)) {
      audit$geometry_status[[i]] <- "motif_pair_missing"
      next
    }
    gene_index <- base::match(audit$query[[i]], genes$locus_tag)
    protein <- genes[gene_index, , drop = FALSE]
    protein$locus_tag <- audit$query[[i]]
    geometry <- .dnmb_rebasefinder_verify_motif_geometry(
      motif_hits = motif_hits,
      protein_table = protein,
      structure_dirs = output_model_dir,
      structure_paths = stats::setNames(audit$model_path[[i]], audit$query[[i]])
    )
    summary <- geometry$summary
    if (base::nrow(summary)) {
      geometry_status <- summary$structural_adjacency_status[[1]]
      role <- base::as.character(audit$enzyme_role[[i]])
      family <- if ("family_id" %in% base::names(motif_hits) && base::nrow(motif_hits)) {
        base::as.character(motif_hits$family_id[[1]])
      } else {
        base::as.character(audit$family_id[[i]])
      }
      canonical_family <- .dnmb_rebasefinder_canonical_structure_family(family)
      restriction_motor <- role %in% c("R", "RM") &&
        canonical_family %in% c("Type I", "Type III")
      nuclease_motifs <- if (base::identical(canonical_family, "Type I")) {
        "HsdR-PD"
      } else if (base::identical(canonical_family, "Type III")) {
        "ResIII-PD"
      } else {
        c("PD-ExK", "HNH", "GIY-YIG", "PLD-HKD")
      }
      has_nuclease <- base::any(motif_hits$motif %in% nuclease_motifs)
      if (restriction_motor && !has_nuclease && geometry_status %in%
          c("homology_model_supported", "homology_model_partial")) {
        geometry_status <- "homology_model_motor_only"
      }
      audit$geometry_status[[i]] <- geometry_status
      audit$model_supported[[i]] <- geometry_status == "homology_model_supported" &&
        (!restriction_motor || has_nuclease)
    } else {
      audit$geometry_status[[i]] <- "motif_pair_missing"
    }
  }
  files <- .dnmb_rebasefinder_write_homology_audit(output_dir, audit)
  files$model_manifest <- if (!base::is.na(model_manifest) && base::file.exists(model_manifest)) {
    base::normalizePath(model_manifest, winslash = "/", mustWork = FALSE)
  } else {
    NA_character_
  }
  hits <- .dnmb_rebasefinder_merge_homology(hits, audit)
  n_mapped <- base::sum(audit$selection_status == "template_mapped", na.rm = TRUE)
  n_built <- base::sum(audit$model_status == "template_model_built", na.rm = TRUE)
  n_supported <- base::sum(audit$model_supported %in% TRUE, na.rm = TRUE)
  list(
    hits = hits,
    audit = audit,
    status = if (n_mapped && n_built < n_mapped) "partial" else "ok",
    detail = base::sprintf(
      "%d/%d selected candidates mapped to a positive template; %d models built; %d motif geometries supported (%s)",
      n_mapped, base::nrow(audit), n_built, n_supported, backend$detail
    ),
    files = files
  )
}

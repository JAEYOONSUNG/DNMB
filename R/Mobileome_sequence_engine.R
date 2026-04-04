.dnmb_empty_sequence_engine <- function() {
  list(
    status = tibble::tibble(
      component = character(),
      status = character(),
      detail = character()
    ),
    orfs = tibble::tibble(),
    hmmer_hits = tibble::tibble(),
    sequence_elements = tibble::tibble(),
    reference_hits = tibble::tibble(),
    files = list()
  )
}

.dnmb_run_sequence_engine <- function(
  parsed,
  output_dir,
  prodigal_mode = "single",
  isescan_dir = NULL,
  isfinder_db = NULL,
  isfinder_fasta = NULL,
  use_genbank_proteins = TRUE,
  verbose = TRUE
) {
  result <- .dnmb_empty_sequence_engine()
  resources <- .dnmb_detect_sequence_resources(isescan_dir = isescan_dir)
  result$status <- resources$status

  # Only hmmsearch + HMM DB required (prodigal optional when using GenBank proteins)
  required_ok <- nzchar(resources$hmmsearch) && file.exists(resources$hmm_db)
  if (!required_ok) {
    result$status <- dplyr::bind_rows(
      result$status,
      tibble::tibble(
        component = "sequence_engine",
        status = "skipped",
        detail = "hmmsearch or ISEScan HMM database is missing"
      )
    )
    return(result)
  }

  engine_dir <- file.path(output_dir, "sequence_engine")
  dir.create(engine_dir, recursive = TRUE, showWarnings = FALSE)

  contig_fasta <- file.path(engine_dir, "genome_contigs.fna")
  .dnmb_mobileome_message(verbose, "Writing contig FASTA for sequence engine")
  .dnmb_write_contig_fasta(parsed$metadata, contig_fasta)

  # --- Use GenBank CDS proteins directly (skip Prodigal) ---
  if (isTRUE(use_genbank_proteins)) {
    .dnmb_mobileome_message(verbose, "Using GenBank CDS proteins (skipping Prodigal)")
    gene_table <- .dnmb_build_gene_table(parsed$features)
    proteins <- gene_table[!is.na(gene_table$translation) & nzchar(gene_table$translation), , drop = FALSE]

    if (!nrow(proteins)) {
      result$status <- dplyr::bind_rows(result$status,
        tibble::tibble(component = "genbank_proteins", status = "empty",
                       detail = "No CDS translations in GenBank"))
      return(result)
    }

    # Write GenBank proteins as FASTA (clean sequences: amino acids only)
    genbank_faa <- file.path(engine_dir, "genbank_proteins.faa")
    con <- file(genbank_faa, "w")
    for (i in seq_len(nrow(proteins))) {
      token <- proteins$locus_tag[i]
      if (is.na(token) || !nzchar(token)) token <- proteins$gene_id[i]
      # Clean: keep only valid amino acid characters
      aa_seq <- gsub("[^A-Za-z*]", "", proteins$translation[i])
      if (nchar(aa_seq) < 10) next
      writeLines(c(paste0(">", token), aa_seq), con)
    }
    close(con)

    # Build ORF table from GenBank CDS (same format as Prodigal output)
    result$orfs <- tibble::tibble(
      protein_token = dplyr::coalesce(proteins$locus_tag, proteins$gene_id),
      contig = proteins$contig,
      start = proteins$start,
      end = proteins$end,
      strand = proteins$strand,
      orf_length_bp = proteins$end - proteins$start + 1L,
      aa_length = nchar(proteins$translation)
    )

    result$files$contig_fasta <- contig_fasta
    result$files$prodigal_faa <- genbank_faa
    result$status <- dplyr::bind_rows(result$status,
      tibble::tibble(component = "genbank_proteins", status = "ok",
                     detail = paste0(nrow(proteins), " CDS proteins from GenBank")))

    hmmer_run <- .dnmb_run_hmmer_pipeline(
      proteins_faa = genbank_faa,
      resources = resources,
      output_dir = engine_dir,
      verbose = verbose
    )
  } else {
    # --- Original Prodigal pipeline ---
    if (!nzchar(resources$prodigal)) {
      result$status <- dplyr::bind_rows(result$status,
        tibble::tibble(component = "sequence_engine", status = "skipped",
                       detail = "prodigal not found and use_genbank_proteins=FALSE"))
      return(result)
    }

    prodigal_run <- tryCatch(
      .dnmb_run_prodigal(
        contig_fasta = contig_fasta,
        output_dir = engine_dir,
        prodigal = resources$prodigal,
        prodigal_mode = prodigal_mode,
        verbose = verbose
      ),
      error = function(e) e
    )
    if (inherits(prodigal_run, "error")) {
      result$status <- dplyr::bind_rows(
        result$status,
        tibble::tibble(component = "prodigal", status = "failed",
                       detail = conditionMessage(prodigal_run)))
      return(result)
    }

    result$files$contig_fasta <- contig_fasta
    result$files$prodigal_gff <- prodigal_run$gff
    result$files$prodigal_faa <- prodigal_run$faa
    result$files$prodigal_fna <- prodigal_run$fna
    result$status <- dplyr::bind_rows(
      result$status,
      tibble::tibble(component = "prodigal", status = "ok", detail = prodigal_run$gff))

    result$orfs <- .dnmb_parse_prodigal_outputs(
      gff_path = prodigal_run$gff,
      faa_path = prodigal_run$faa
    )

    hmmer_run <- .dnmb_run_hmmer_pipeline(
      proteins_faa = prodigal_run$faa,
      resources = resources,
      output_dir = engine_dir,
      verbose = verbose
    )
  }  # end if/else use_genbank_proteins

  result$status <- dplyr::bind_rows(result$status, hmmer_run$status)
  result$files <- c(result$files, hmmer_run$files)
  result$hmmer_hits <- hmmer_run$hits

  if (!nrow(result$orfs) || !nrow(result$hmmer_hits)) {
    result$status <- dplyr::bind_rows(
      result$status,
      tibble::tibble(
        component = "sequence_candidates",
        status = "empty",
        detail = "no ORFs or no HMM hits passed filtering"
      )
    )
    return(result)
  }

  family_rules <- .dnmb_load_isescan_family_rules(resources$constants_path)
  result$sequence_elements <- .dnmb_build_sequence_candidates(
    orfs = result$orfs,
    hmmer_hits = result$hmmer_hits,
    metadata = parsed$metadata,
    family_rules = family_rules
  )
  result$status <- dplyr::bind_rows(
    result$status,
    tibble::tibble(
      component = "sequence_candidates",
      status = if (nrow(result$sequence_elements)) "ok" else "empty",
      detail = paste("candidate elements:", nrow(result$sequence_elements))
    )
  )

  reference_run <- .dnmb_run_reference_search(
    sequence_elements = result$sequence_elements,
    metadata = parsed$metadata,
    output_dir = engine_dir,
    resources = resources,
    isfinder_db = isfinder_db,
    isfinder_fasta = isfinder_fasta,
    verbose = verbose
  )
  result$status <- dplyr::bind_rows(result$status, reference_run$status)
  result$reference_hits <- reference_run$hits
  result$files <- c(result$files, reference_run$files)

  if (nrow(result$sequence_elements) && nrow(result$reference_hits)) {
    top_ref <- result$reference_hits %>%
      dplyr::arrange(.data$query_id, .data$evalue, dplyr::desc(.data$pident), dplyr::desc(.data$qcov)) %>%
      dplyr::group_by(.data$query_id) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(
        sequence_element_id = .data$query_id,
        reference_hit_id = .data$subject_id,
        reference_hit_family = .data$subject_family,
        reference_hit_identity = .data$pident,
        reference_hit_qcov = .data$qcov,
        reference_hit_evalue = .data$evalue
      )

    result$sequence_elements <- result$sequence_elements %>%
      dplyr::left_join(top_ref, by = "sequence_element_id") %>%
      dplyr::mutate(
        element_family = dplyr::coalesce(.data$reference_hit_family, .data$element_family),
        confidence_score = .data$confidence_score + dplyr::if_else(!is.na(.data$reference_hit_id), 1L, 0L),
        confidence = dplyr::case_when(
          .data$confidence_score >= 5L ~ "high",
          .data$confidence_score >= 3L ~ "medium",
          TRUE ~ "low"
        ),
        evidence_summary = dplyr::case_when(
          !is.na(.data$reference_hit_id) ~ stringr::str_squish(paste(
            .data$evidence_summary,
            paste0("reference=", .data$reference_hit_id),
            paste0("pident=", sprintf("%.1f", .data$reference_hit_identity))
          )),
          TRUE ~ .data$evidence_summary
        )
      )
  }

  result
}

.dnmb_detect_sequence_resources <- function(isescan_dir = NULL) {
  isescan_exe <- Sys.which("isescan.py")
  inferred_dir <- isescan_dir
  if (is.null(inferred_dir) || !nzchar(inferred_dir)) {
    if (nzchar(isescan_exe)) {
      inferred_dir <- dirname(isescan_exe)
    } else {
      inferred_dir <- ""
    }
  }

  constants_path <- if (nzchar(inferred_dir)) file.path(inferred_dir, "constants.py") else ""
  hmm_db <- if (nzchar(inferred_dir)) file.path(inferred_dir, "pHMMs", "clusters.faa.hmm") else ""
  single_db <- if (nzchar(inferred_dir)) file.path(inferred_dir, "pHMMs", "clusters.single.faa") else ""

  prodigal <- Sys.which("prodigal")
  hmmsearch <- Sys.which("hmmsearch")
  phmmer <- Sys.which("phmmer")
  blastn <- Sys.which("blastn")
  makeblastdb <- Sys.which("makeblastdb")

  status <- tibble::tibble(
    component = c(
      "isescan_dir",
      "isescan_constants",
      "isescan_hmm_db",
      "isescan_singleton_db",
      "prodigal",
      "hmmsearch",
      "phmmer",
      "blastn",
      "makeblastdb"
    ),
    status = c(
      if (nzchar(inferred_dir) && dir.exists(inferred_dir)) "found" else "missing",
      if (nzchar(constants_path) && file.exists(constants_path)) "found" else "missing",
      if (nzchar(hmm_db) && file.exists(hmm_db)) "found" else "missing",
      if (nzchar(single_db) && file.exists(single_db)) "found" else "missing",
      if (nzchar(prodigal)) "found" else "missing",
      if (nzchar(hmmsearch)) "found" else "missing",
      if (nzchar(phmmer)) "found" else "missing",
      if (nzchar(blastn)) "found" else "missing",
      if (nzchar(makeblastdb)) "found" else "missing"
    ),
    detail = c(
      inferred_dir,
      constants_path,
      hmm_db,
      single_db,
      prodigal,
      hmmsearch,
      phmmer,
      blastn,
      makeblastdb
    )
  )

  list(
    isescan_dir = inferred_dir,
    constants_path = constants_path,
    hmm_db = hmm_db,
    single_db = single_db,
    prodigal = prodigal,
    hmmsearch = hmmsearch,
    phmmer = phmmer,
    blastn = blastn,
    makeblastdb = makeblastdb,
    status = status
  )
}

.dnmb_write_contig_fasta <- function(metadata, path) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(nrow(metadata))) {
    writeLines(paste0(">", metadata$contig[[i]]), con = con)
    writeLines(metadata$sequence[[i]], con = con)
  }
  invisible(path)
}

.dnmb_run_prodigal <- function(contig_fasta, output_dir, prodigal, prodigal_mode = "single", verbose = TRUE) {
  gff_path <- file.path(output_dir, "prodigal.gff")
  faa_path <- file.path(output_dir, "prodigal.faa")
  fna_path <- file.path(output_dir, "prodigal.fna")

  args <- c(
    "-i", shQuote(contig_fasta),
    "-a", shQuote(faa_path),
    "-d", shQuote(fna_path),
    "-o", shQuote(gff_path),
    "-f", "gff",
    "-p", prodigal_mode,
    "-q"
  )

  .dnmb_mobileome_message(verbose, "Running Prodigal")
  .dnmb_run_system_command(prodigal, args)

  if (!file.exists(gff_path) || !file.exists(faa_path)) {
    stop("Prodigal did not produce the expected output files.")
  }

  list(gff = gff_path, faa = faa_path, fna = fna_path)
}

.dnmb_run_system_command <- function(cmd, args) {
  status <- system2(cmd, args = args, stdout = TRUE, stderr = TRUE)
  exit_status <- attr(status, "status")
  if (!is.null(exit_status) && exit_status != 0L) {
    stop(paste(c(cmd, args), collapse = " "), "\n", paste(status, collapse = "\n"))
  }
  invisible(status)
}

.dnmb_parse_prodigal_outputs <- function(gff_path, faa_path) {
  orfs <- .dnmb_parse_prodigal_gff(gff_path)
  if (!nrow(orfs)) {
    return(tibble::tibble())
  }

  aa_set <- Biostrings::readAAStringSet(faa_path)
  aa_tbl <- tibble::tibble(
    aa_header = names(aa_set),
    protein_token = sub("\\s+.*$", "", names(aa_set)),
    aa_sequence = as.character(aa_set),
    aa_length = nchar(as.character(aa_set))
  )

  orfs <- orfs %>%
    dplyr::left_join(aa_tbl, by = c("protein_token" = "protein_token"))

  orfs
}

.dnmb_parse_prodigal_gff <- function(gff_path) {
  lines <- readLines(gff_path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  if (!length(lines)) {
    return(tibble::tibble())
  }

  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[vapply(fields, length, integer(1)) >= 9L]
  if (!length(fields)) {
    return(tibble::tibble())
  }

  gff_tbl <- tibble::tibble(
    contig = vapply(fields, `[[`, character(1), 1L),
    source = vapply(fields, `[[`, character(1), 2L),
    feature_type = vapply(fields, `[[`, character(1), 3L),
    start = as.integer(vapply(fields, `[[`, character(1), 4L)),
    end = as.integer(vapply(fields, `[[`, character(1), 5L)),
    score = vapply(fields, `[[`, character(1), 6L),
    strand = vapply(fields, `[[`, character(1), 7L),
    phase = vapply(fields, `[[`, character(1), 8L),
    attributes = vapply(fields, `[[`, character(1), 9L)
  ) %>%
    dplyr::mutate(
      attr_id = vapply(.data$attributes, .dnmb_extract_gff_attribute, character(1), key = "ID"),
      attr_partial = vapply(.data$attributes, .dnmb_extract_gff_attribute, character(1), key = "partial"),
      protein_token = vapply(
        seq_len(dplyr::n()),
        function(i) .dnmb_guess_prodigal_token(.data$contig[[i]], .data$attr_id[[i]]),
        character(1)
      ),
      orf_length_bp = .data$end - .data$start + 1L
    )

  gff_tbl
}

.dnmb_extract_gff_attribute <- function(attribute_string, key) {
  pattern <- paste0("(^|;)", key, "=([^;]+)")
  match <- regmatches(attribute_string, regexec(pattern, attribute_string, perl = TRUE))[[1]]
  if (length(match) >= 3L) {
    match[3L]
  } else {
    ""
  }
}

.dnmb_guess_prodigal_token <- function(contig, attr_id) {
  if (is.na(attr_id) || !nzchar(attr_id)) {
    return(NA_character_)
  }
  suffix <- sub("^.*_", "", attr_id)
  if (!nzchar(suffix)) {
    suffix <- attr_id
  }
  paste(contig, suffix, sep = "_")
}

.dnmb_run_hmmer_pipeline <- function(proteins_faa, resources, output_dir, verbose = TRUE) {
  status <- tibble::tibble(component = character(), status = character(), detail = character())
  files <- list()
  hits <- tibble::tibble()

  hmm_tbl <- file.path(output_dir, "isescan_hmmsearch.tbl")
  hmm_log <- file.path(output_dir, "isescan_hmmsearch.log")
  .dnmb_mobileome_message(verbose, "Running hmmsearch against ISEScan profile HMMs")
  .dnmb_run_system_command(
    resources$hmmsearch,
    c(
      "--tblout", shQuote(hmm_tbl),
      "--noali",
      "--cpu", "1",
      shQuote(resources$hmm_db),
      shQuote(proteins_faa)
    )
  )
  files$hmm_tbl <- hmm_tbl
  files$hmm_log <- hmm_log
  status <- dplyr::bind_rows(
    status,
    tibble::tibble(component = "hmmsearch", status = "ok", detail = hmm_tbl)
  )
  hits <- dplyr::bind_rows(
    hits,
    .dnmb_parse_hmmer_tblout(hmm_tbl, search_type = "hmmsearch")
  )

  if (nzchar(resources$phmmer) && file.exists(resources$single_db)) {
    single_tbl <- file.path(output_dir, "isescan_singleton_phmmer.tbl")
    single_log <- file.path(output_dir, "isescan_singleton_phmmer.log")
    .dnmb_mobileome_message(verbose, "Running phmmer against ISEScan singleton references")
    .dnmb_run_system_command(
      resources$phmmer,
      c(
        "--tblout", shQuote(single_tbl),
        "--noali",
        "--cpu", "1",
        shQuote(resources$single_db),
        shQuote(proteins_faa)
      )
    )
    files$single_tbl <- single_tbl
    files$single_log <- single_log
    status <- dplyr::bind_rows(
      status,
      tibble::tibble(component = "phmmer", status = "ok", detail = single_tbl)
    )
    hits <- dplyr::bind_rows(
      hits,
      .dnmb_parse_hmmer_tblout(single_tbl, search_type = "phmmer")
    )
  } else {
    status <- dplyr::bind_rows(
      status,
      tibble::tibble(component = "phmmer", status = "skipped", detail = "phmmer or singleton database missing")
    )
  }

  if (nrow(hits)) {
    hits <- hits %>%
      dplyr::mutate(
        family = vapply(.data$query_name, .dnmb_extract_family_from_isescan_query, character(1)),
        cluster_id = vapply(.data$query_name, .dnmb_extract_cluster_id, character(1)),
        pass_default_filter = (.data$best_domain_evalue <= 1e-5 | .data$full_evalue <= 1e-5) &
          (.data$family != "new" | .data$best_domain_evalue <= 1e-20)
      ) %>%
      dplyr::filter(.data$pass_default_filter) %>%
      dplyr::arrange(.data$target_name, .data$best_domain_evalue, .data$full_evalue)
  }

  list(status = status, files = files, hits = hits)
}

.dnmb_parse_hmmer_tblout <- function(path, search_type) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  if (!length(lines)) {
    return(tibble::tibble(
      target_name = character(),
      query_name = character(),
      full_evalue = numeric(),
      full_score = numeric(),
      best_domain_evalue = numeric(),
      best_domain_score = numeric(),
      overlap_number = integer(),
      search_type = character()
    ))
  }

  parsed <- lapply(lines, function(line) {
    fields <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(fields) < 18L) {
      return(NULL)
    }
    tibble::tibble(
      target_name = fields[[1]],
      query_name = fields[[3]],
      full_evalue = suppressWarnings(as.numeric(fields[[5]])),
      full_score = suppressWarnings(as.numeric(fields[[6]])),
      best_domain_evalue = suppressWarnings(as.numeric(fields[[8]])),
      best_domain_score = suppressWarnings(as.numeric(fields[[9]])),
      overlap_number = suppressWarnings(as.integer(fields[[14]])),
      search_type = search_type
    )
  })
  dplyr::bind_rows(parsed)
}

.dnmb_extract_family_from_isescan_query <- function(query_name) {
  clean <- sub("\\.faa$", "", query_name)
  clean <- sub("\\|.*$", "", clean)
  family <- sub("_.*$", "", clean)
  if (family == "IS200") {
    return("IS200/IS605")
  }
  if (grepl("^IS200_IS605", clean)) {
    return("IS200/IS605")
  }
  family
}

.dnmb_extract_cluster_id <- function(query_name) {
  clean <- sub("\\|.*$", "", query_name)
  sub("\\.faa$", "", clean)
}

.dnmb_load_isescan_family_rules <- function(constants_path) {
  if (!nzchar(constants_path) || !file.exists(constants_path)) {
    return(tibble::tibble())
  }
  lines <- readLines(constants_path, warn = FALSE)

  is_len <- .dnmb_parse_python_tuple_dict(lines, "minMaxLen4is", c("is_min_len_bp", "is_max_len_bp"))
  tpase <- .dnmb_parse_python_tuple_dict(lines, "minMax4tpase", c("tpase_min_bp", "tpase_max_bp", "tpase_min_pep_aa"))
  tir <- .dnmb_parse_python_tuple_dict(lines, "minMax4tir", c("tir_min_len_bp", "tir_max_len_bp", "tir_opt_len_bp", "tir_presence_flag"))

  outer_dist <- .dnmb_parse_python_numeric_tuple(lines, "outerDist4ter2tpase", c("outer_near_bp", "outer_far_bp"))
  min_ir_identity <- .dnmb_parse_python_scalar(lines, "minIrIdentity")
  opt_ir_identity <- .dnmb_parse_python_scalar(lines, "optIrIdentity")
  min_dist4ter2orf <- .dnmb_parse_python_scalar(lines, "minDist4ter2orf")

  rules <- is_len %>%
    dplyr::full_join(tpase, by = "family") %>%
    dplyr::full_join(tir, by = "family") %>%
    dplyr::mutate(
      outer_near_bp = outer_dist$outer_near_bp[[1]],
      outer_far_bp = outer_dist$outer_far_bp[[1]],
      min_ir_identity = min_ir_identity,
      opt_ir_identity = opt_ir_identity,
      min_dist4ter2orf = min_dist4ter2orf
    )

  rules
}

.dnmb_parse_python_tuple_dict <- function(lines, var_name, field_names) {
  start_idx <- grep(paste0("^", var_name, "\\s*=\\s*\\{"), lines)
  if (!length(start_idx)) {
    return(tibble::tibble(family = character()))
  }

  block <- character()
  depth <- 0L
  for (line in lines[start_idx[1]:length(lines)]) {
    depth <- depth + stringr::str_count(line, stringr::fixed("{")) - stringr::str_count(line, stringr::fixed("}"))
    block <- c(block, line)
    if (depth == 0L) {
      break
    }
  }

  entry_lines <- block[grepl("^\\s*'[^']+'\\s*:\\s*\\(", block)]
  rows <- lapply(entry_lines, function(line) {
    family <- sub("^\\s*'([^']+)'.*$", "\\1", line)
    tuple_txt <- sub("^.*\\(([^)]*)\\).*$", "\\1", line)
    values <- trimws(strsplit(tuple_txt, ",", fixed = TRUE)[[1]])
    values <- values[nzchar(values)]
    out <- as.list(rep(NA_real_, length(field_names)))
    for (i in seq_len(min(length(field_names), length(values)))) {
      out[[i]] <- suppressWarnings(as.numeric(values[[i]]))
    }
    tibble::as_tibble(stats::setNames(c(list(family = family), out), c("family", field_names)))
  })
  dplyr::bind_rows(rows)
}

.dnmb_parse_python_numeric_tuple <- function(lines, var_name, field_names) {
  idx <- grep(paste0("^", var_name, "\\s*=\\s*\\("), lines)
  if (!length(idx)) {
    return(tibble::as_tibble(stats::setNames(as.list(rep(NA_real_, length(field_names))), field_names)))
  }
  tuple_txt <- sub("^.*\\(([^)]*)\\).*$", "\\1", lines[idx[1]])
  values <- trimws(strsplit(tuple_txt, ",", fixed = TRUE)[[1]])
  values <- values[nzchar(values)]
  out <- as.list(rep(NA_real_, length(field_names)))
  for (i in seq_len(min(length(field_names), length(values)))) {
    out[[i]] <- suppressWarnings(as.numeric(values[[i]]))
  }
  tibble::as_tibble(stats::setNames(out, field_names))
}

.dnmb_parse_python_scalar <- function(lines, var_name) {
  idx <- grep(paste0("^", var_name, "\\s*="), lines)
  if (!length(idx)) {
    return(NA_real_)
  }
  value_txt <- sub("^.*=\\s*", "", lines[idx[1]])
  value_txt <- sub("\\s+#.*$", "", value_txt)
  suppressWarnings(as.numeric(trimws(value_txt)))
}

.dnmb_build_sequence_candidates <- function(orfs, hmmer_hits, metadata, family_rules) {
  if (!nrow(orfs) || !nrow(hmmer_hits)) {
    return(tibble::tibble())
  }

  hit_summary <- hmmer_hits %>%
    dplyr::arrange(.data$target_name, .data$best_domain_evalue, .data$full_evalue) %>%
    dplyr::group_by(.data$target_name) %>%
    dplyr::summarise(
      query_name = dplyr::first(.data$query_name),
      element_family = dplyr::first(.data$family),
      cluster_id = dplyr::first(.data$cluster_id),
      full_evalue = dplyr::first(.data$full_evalue),
      full_score = dplyr::first(.data$full_score),
      best_domain_evalue = dplyr::first(.data$best_domain_evalue),
      best_domain_score = dplyr::first(.data$best_domain_score),
      hit_count = dplyr::n(),
      search_types = paste(unique(.data$search_type), collapse = "; "),
      family_support = paste(unique(.data$family), collapse = "; "),
      .groups = "drop"
    )

  candidates <- orfs %>%
    dplyr::inner_join(hit_summary, by = c("protein_token" = "target_name")) %>%
    dplyr::left_join(metadata %>% dplyr::select(contig, contig_number, definition, sequence_length_bp, sequence), by = "contig") %>%
    dplyr::filter(!is.na(.data$sequence_length_bp), .data$orf_length_bp >= 300L)

  if (!nrow(candidates)) {
    return(tibble::tibble())
  }

  rows <- lapply(seq_len(nrow(candidates)), function(i) {
    row <- candidates[i, , drop = FALSE]
    rule <- .dnmb_get_family_rule(family_rules, row$element_family[[1]])
    tpase_length_ok <- !is.na(rule$tpase_min_bp) &&
      !is.na(rule$tpase_max_bp) &&
      row$orf_length_bp[[1]] >= rule$tpase_min_bp &&
      row$orf_length_bp[[1]] <= rule$tpase_max_bp
    tir <- .dnmb_empty_tir_result()
    if (tpase_length_ok && row$element_family[[1]] != "new" && row$best_domain_evalue[[1]] <= 1e-20) {
      tir <- .dnmb_find_best_tir(
        seq = row$sequence[[1]],
        orf_start = row$start[[1]],
        orf_end = row$end[[1]],
        family_rule = rule
      )
    }

    element_start <- if (isTRUE(tir$found)) tir$start1 else row$start[[1]]
    element_end <- if (isTRUE(tir$found)) tir$end2 else row$end[[1]]
    tsd <- .dnmb_find_best_tsd(
      seq = row$sequence[[1]],
      element_start = element_start,
      element_end = element_end
    )

    size_bp <- element_end - element_start + 1L
    is_length_ok <- !is.na(rule$is_min_len_bp) &&
      !is.na(rule$is_max_len_bp) &&
      size_bp >= rule$is_min_len_bp &&
      size_bp <= rule$is_max_len_bp

    confidence_score <- 1L +
      as.integer(row$best_domain_evalue[[1]] <= 1e-20) +
      as.integer(tpase_length_ok) +
      as.integer(isTRUE(tir$found)) +
      as.integer(isTRUE(tsd$found)) +
      as.integer(row$element_family[[1]] != "new")

    confidence <- dplyr::case_when(
      confidence_score >= 5L ~ "high",
      confidence_score >= 3L ~ "medium",
      TRUE ~ "low"
    )

    tibble::tibble(
      sequence_element_id = paste0("SEQIS_", sprintf("%04d", i)),
      contig = row$contig[[1]],
      contig_number = row$contig_number[[1]],
      definition = row$definition[[1]],
      start = element_start,
      end = element_end,
      size_bp = size_bp,
      strand = row$strand[[1]],
      element_family = row$element_family[[1]],
      cluster_id = row$cluster_id[[1]],
      query_name = row$query_name[[1]],
      query_sources = row$search_types[[1]],
      family_support = row$family_support[[1]],
      hit_count = row$hit_count[[1]],
      full_evalue = row$full_evalue[[1]],
      best_domain_evalue = row$best_domain_evalue[[1]],
      orf_start = row$start[[1]],
      orf_end = row$end[[1]],
      orf_length_bp = row$orf_length_bp[[1]],
      aa_length = row$aa_length[[1]],
      tpase_length_ok = tpase_length_ok,
      is_length_ok = is_length_ok,
      tir_found = isTRUE(tir$found),
      tir_start1 = tir$start1,
      tir_end1 = tir$end1,
      tir_start2 = tir$start2,
      tir_end2 = tir$end2,
      tir_len_bp = tir$tir_len_bp,
      tir_identity = tir$tir_identity,
      tir_score = tir$tir_score,
      tsd_found = isTRUE(tsd$found),
      tsd_len_bp = tsd$tsd_len_bp,
      tsd_seq = tsd$tsd_seq,
      confidence = confidence,
      confidence_score = confidence_score,
      evidence_label = "ISEScan_like_sequence_engine",
      evidence_summary = stringr::str_squish(paste(
        row$query_name[[1]],
        paste0("best_domain_evalue=", signif(row$best_domain_evalue[[1]], 3)),
        paste0("sources=", row$search_types[[1]])
      )),
      reference_hit_id = NA_character_,
      reference_hit_family = NA_character_,
      reference_hit_identity = NA_real_
    )
  })

  dplyr::bind_rows(rows)
}

.dnmb_get_family_rule <- function(family_rules, family) {
  if (!nrow(family_rules)) {
    return(tibble::as_tibble(list(
      family = family,
      is_min_len_bp = 400,
      is_max_len_bp = 10000,
      tpase_min_bp = 300,
      tpase_max_bp = 2100,
      tpase_min_pep_aa = 50,
      tir_min_len_bp = 10,
      tir_max_len_bp = 50,
      tir_opt_len_bp = 20,
      tir_presence_flag = -1,
      outer_near_bp = 150,
      outer_far_bp = 500,
      min_ir_identity = 0.4,
      opt_ir_identity = 0.6,
      min_dist4ter2orf = -150
    )))
  }
  family_name <- family
  hit <- family_rules %>% dplyr::filter(.data$family == .env$family_name)
  if (!nrow(hit)) {
    hit <- family_rules %>% dplyr::filter(.data$family == "new")
  }
  hit[1, , drop = FALSE]
}

.dnmb_find_best_tir <- function(seq, orf_start, orf_end, family_rule) {
  empty <- .dnmb_empty_tir_result()

  if (is.na(seq) || !nzchar(seq)) {
    return(empty)
  }
  if (!is.na(family_rule$tir_presence_flag[[1]]) && family_rule$tir_presence_flag[[1]] == 0) {
    return(empty)
  }

  seq_len <- nchar(seq)
  min_len <- as.integer(dplyr::coalesce(family_rule$tir_min_len_bp[[1]], 10))
  max_len <- as.integer(dplyr::coalesce(family_rule$tir_max_len_bp[[1]], 50))
  opt_len <- as.integer(dplyr::coalesce(family_rule$tir_opt_len_bp[[1]], min_len))
  seed_len <- max(6L, min(12L, opt_len))
  search_far_bp <- min(250L, as.integer(dplyr::coalesce(family_rule$outer_far_bp[[1]], 500)))
  left_from <- max(1L, orf_start - search_far_bp)
  left_to <- min(seq_len - seed_len + 1L, orf_start + abs(as.integer(dplyr::coalesce(family_rule$min_dist4ter2orf[[1]], -150))))
  right_end_from <- max(seed_len, orf_end - abs(as.integer(dplyr::coalesce(family_rule$min_dist4ter2orf[[1]], -150))))
  right_end_to <- min(seq_len, orf_end + search_far_bp)
  if (left_from > left_to || right_end_from > right_end_to) {
    return(empty)
  }

  left_positions <- left_from:left_to
  left_seeds <- substring(seq, left_positions, left_positions + seed_len - 1L)
  seed_map <- split(left_positions, left_seeds)

  best <- empty
  best_score <- -Inf

  for (right_end in right_end_from:right_end_to) {
    right_start_seed <- right_end - seed_len + 1L
    if (right_start_seed < 1L) {
      next
    }
    right_seed <- substring(seq, right_start_seed, right_end)
    key <- .dnmb_revcomp_string(right_seed)
    matched_left <- seed_map[[key]]
    if (is.null(matched_left)) {
      next
    }

    for (left_start in matched_left) {
      max_extend <- min(
        max_len,
        seq_len - left_start + 1L,
        right_end
      )
      mismatches <- 0L
      best_local <- NULL
      for (tir_len in seed_len:max_extend) {
        if (tir_len > seed_len) {
          left_base <- substring(seq, left_start + tir_len - 1L, left_start + tir_len - 1L)
          right_base <- substring(seq, right_end - tir_len + 1L, right_end - tir_len + 1L)
          if (!identical(left_base, .dnmb_complement_base(right_base))) {
            mismatches <- mismatches + 1L
          }
        }
        identity <- (tir_len - mismatches) / tir_len
        if (tir_len < min_len || identity < dplyr::coalesce(family_rule$min_ir_identity[[1]], 0.4)) {
          next
        }

        candidate_start2 <- right_end - tir_len + 1L
        candidate_is_len <- right_end - left_start + 1L
        if (
          !is.na(family_rule$is_min_len_bp[[1]]) &&
          !is.na(family_rule$is_max_len_bp[[1]]) &&
          (candidate_is_len < family_rule$is_min_len_bp[[1]] || candidate_is_len > family_rule$is_max_len_bp[[1]])
        ) {
          next
        }
        if (left_start > orf_start || right_end < orf_end) {
          next
        }

        score <- identity * tir_len - abs(tir_len - opt_len) / 10
        best_local <- list(
          found = TRUE,
          start1 = left_start,
          end1 = left_start + tir_len - 1L,
          start2 = candidate_start2,
          end2 = right_end,
          tir_len_bp = tir_len,
          tir_identity = identity,
          tir_score = score
        )
      }

      if (!is.null(best_local) && best_local$tir_score > best_score) {
        best <- best_local
        best_score <- best_local$tir_score
      }
    }
  }

  best
}

.dnmb_empty_tir_result <- function() {
  list(
    found = FALSE,
    start1 = NA_integer_,
    end1 = NA_integer_,
    start2 = NA_integer_,
    end2 = NA_integer_,
    tir_len_bp = NA_integer_,
    tir_identity = NA_real_,
    tir_score = NA_real_
  )
}

.dnmb_complement_base <- function(base) {
  switch(
    toupper(base),
    A = "T",
    T = "A",
    C = "G",
    G = "C",
    N = "N",
    base
  )
}

.dnmb_revcomp_string <- function(seq) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq)))
}

.dnmb_find_best_tsd <- function(seq, element_start, element_end, min_len = 2L, max_len = 15L) {
  empty <- list(found = FALSE, tsd_len_bp = NA_integer_, tsd_seq = NA_character_)
  if (is.na(seq) || !nzchar(seq) || element_start <= 1L || element_end >= nchar(seq)) {
    return(empty)
  }

  best_len <- NA_integer_
  best_seq <- NA_character_
  for (len in seq(max_len, min_len)) {
    left_start <- element_start - len
    left_end <- element_start - 1L
    right_start <- element_end + 1L
    right_end <- element_end + len
    if (left_start < 1L || right_end > nchar(seq)) {
      next
    }
    left_seq <- substring(seq, left_start, left_end)
    right_seq <- substring(seq, right_start, right_end)
    if (identical(left_seq, right_seq)) {
      best_len <- len
      best_seq <- left_seq
      break
    }
  }

  if (is.na(best_len)) {
    return(empty)
  }
  list(found = TRUE, tsd_len_bp = best_len, tsd_seq = best_seq)
}

.dnmb_run_reference_search <- function(
  sequence_elements,
  metadata,
  output_dir,
  resources,
  isfinder_db = NULL,
  isfinder_fasta = NULL,
  verbose = TRUE
) {
  empty <- list(
    status = tibble::tibble(
      component = "reference_search",
      status = "skipped",
      detail = "no reference database provided"
    ),
    hits = tibble::tibble(),
    files = list()
  )

  if (!nrow(sequence_elements)) {
    empty$status$detail <- "no sequence elements available"
    return(empty)
  }

  if (is.null(isfinder_db) && is.null(isfinder_fasta)) {
    return(empty)
  }

  if (!nzchar(resources$blastn)) {
    empty$status$status <- "failed"
    empty$status$detail <- "blastn is missing"
    return(empty)
  }

  db_prefix <- isfinder_db
  files <- list()

  if (!is.null(isfinder_fasta)) {
    if (!file.exists(isfinder_fasta)) {
      empty$status$status <- "failed"
      empty$status$detail <- paste("reference FASTA not found:", isfinder_fasta)
      return(empty)
    }
    if (!nzchar(resources$makeblastdb)) {
      empty$status$status <- "failed"
      empty$status$detail <- "makeblastdb is missing"
      return(empty)
    }
    db_prefix <- file.path(output_dir, "reference", "isfinder_reference")
    dir.create(dirname(db_prefix), recursive = TRUE, showWarnings = FALSE)
    .dnmb_mobileome_message(verbose, "Building BLAST database for IS reference search")
    .dnmb_run_system_command(
      resources$makeblastdb,
      c("-in", shQuote(normalizePath(isfinder_fasta, mustWork = TRUE)), "-dbtype", "nucl", "-out", shQuote(db_prefix))
    )
    files$reference_fasta <- normalizePath(isfinder_fasta, mustWork = TRUE)
    files$reference_db <- db_prefix
  } else {
    files$reference_db <- db_prefix
  }

  candidate_fasta <- file.path(output_dir, "reference_candidates.fna")
  .dnmb_write_candidate_element_fasta(sequence_elements, metadata, candidate_fasta)
  blast_out <- file.path(output_dir, "isfinder_blast.tsv")
  .dnmb_mobileome_message(verbose, "Running reference BLAST search for IS candidates")
  .dnmb_run_system_command(
    resources$blastn,
    c(
      "-task", "blastn",
      "-query", shQuote(candidate_fasta),
      "-db", shQuote(db_prefix),
      "-outfmt", shQuote("6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"),
      "-max_target_seqs", "5",
      "-num_threads", "1",
      "-out", shQuote(blast_out)
    )
  )
  files$candidate_fasta <- candidate_fasta
  files$blast_tsv <- blast_out

  hits <- .dnmb_parse_blast_tabular(blast_out)
  if (nrow(hits)) {
    hits <- hits %>%
      dplyr::mutate(
        qcov = 100 * .data$length / .data$qlen,
        scov = 100 * .data$length / .data$slen,
        subject_family = vapply(.data$subject_id, .dnmb_extract_is_family, character(1))
      ) %>%
      dplyr::filter(.data$pident >= 70, .data$qcov >= 50)
  }

  status <- tibble::tibble(
    component = "reference_search",
    status = if (nrow(hits)) "ok" else "empty",
    detail = if (nrow(hits)) blast_out else "no BLAST hits passed filtering"
  )

  list(status = status, hits = hits, files = files)
}

.dnmb_write_candidate_element_fasta <- function(sequence_elements, metadata, path) {
  seq_map <- stats::setNames(metadata$sequence, metadata$contig)
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(nrow(sequence_elements))) {
    seq <- seq_map[[sequence_elements$contig[[i]]]]
    if (is.null(seq) || !nzchar(seq)) {
      next
    }
    subseq <- substring(seq, sequence_elements$start[[i]], sequence_elements$end[[i]])
    if (sequence_elements$strand[[i]] == "-") {
      subseq <- .dnmb_revcomp_string(subseq)
    }
    writeLines(paste0(">", sequence_elements$sequence_element_id[[i]]), con = con)
    writeLines(subseq, con = con)
  }
  invisible(path)
}

.dnmb_parse_blast_tabular <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(tibble::tibble())
  }
  tbl <- utils::read.table(path, sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  if (!ncol(tbl)) {
    return(tibble::tibble())
  }
  colnames(tbl) <- c(
    "query_id",
    "subject_id",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"
  )
  tibble::as_tibble(tbl)
}

.dnmb_merge_element_calls <- function(annotation_elements, sequence_elements, genes, merge_gap_bp = 300L) {
  ann <- .dnmb_standardize_annotation_elements(annotation_elements)
  seq <- .dnmb_standardize_sequence_elements(sequence_elements)
  combined <- dplyr::bind_rows(ann, seq) %>%
    dplyr::arrange(.data$contig, .data$start, .data$end)

  if (!nrow(combined)) {
    return(.dnmb_empty_elements())
  }

  merged_rows <- list()
  current <- combined[1, , drop = FALSE]
  for (i in 2:nrow(combined)) {
    next_row <- combined[i, , drop = FALSE]
    if (.dnmb_calls_are_mergeable(current, next_row, merge_gap_bp)) {
      current <- .dnmb_merge_call_pair(current, next_row)
    } else {
      merged_rows[[length(merged_rows) + 1L]] <- current
      current <- next_row
    }
  }
  merged_rows[[length(merged_rows) + 1L]] <- current

  merged <- dplyr::bind_rows(merged_rows) %>%
    dplyr::mutate(
      element_id = paste0("IS_", sprintf("%04d", dplyr::row_number())),
      size_bp = .data$end - .data$start + 1L,
      confidence = vapply(.data$confidence_score, .dnmb_score_to_confidence, character(1))
    ) %>%
    .dnmb_add_gene_support(genes) %>%
    dplyr::select(
      element_id,
      contig,
      contig_number,
      definition,
      start,
      end,
      size_bp,
      strand,
      element_family,
      mobile_element_type,
      confidence,
      confidence_score,
      feature_type,
      locus_tag,
      product,
      evidence_label,
      evidence_summary,
      supporting_gene_count,
      supporting_genes,
      supporting_products,
      source_call_types,
      annotation_element_ids,
      sequence_element_ids,
      tir_found,
      tir_len_bp,
      tir_identity,
      tsd_found,
      tsd_len_bp,
      tsd_seq,
      reference_hit_id,
      reference_hit_family,
      reference_hit_identity
    )

  merged
}

.dnmb_standardize_annotation_elements <- function(annotation_elements) {
  if (!nrow(annotation_elements)) {
    return(tibble::tibble())
  }
  annotation_elements %>%
    dplyr::transmute(
      contig = .data$contig,
      contig_number = .data$contig_number,
      definition = .data$definition,
      start = .data$start,
      end = .data$end,
      strand = .data$strand,
      element_family = .data$element_family,
      mobile_element_type = .data$mobile_element_type,
      confidence_score = dplyr::case_when(
        .data$confidence == "high" ~ 4L,
        .data$confidence == "medium" ~ 3L,
        TRUE ~ 2L
      ),
      feature_type = .data$feature_type,
      locus_tag = .data$locus_tag,
      product = .data$product,
      evidence_label = .data$evidence_label,
      evidence_summary = .data$evidence_summary,
      source_call_types = "annotation",
      annotation_element_ids = .data$element_id,
      sequence_element_ids = NA_character_,
      tir_found = NA,
      tir_len_bp = NA_integer_,
      tir_identity = NA_real_,
      tsd_found = NA,
      tsd_len_bp = NA_integer_,
      tsd_seq = NA_character_,
      reference_hit_id = NA_character_,
      reference_hit_family = NA_character_,
      reference_hit_identity = NA_real_
    )
}

.dnmb_standardize_sequence_elements <- function(sequence_elements) {
  if (!nrow(sequence_elements)) {
    return(tibble::tibble())
  }
  sequence_elements %>%
    dplyr::transmute(
      contig = .data$contig,
      contig_number = .data$contig_number,
      definition = .data$definition,
      start = .data$start,
      end = .data$end,
      strand = .data$strand,
      element_family = .data$element_family,
      mobile_element_type = "insertion_sequence_candidate",
      confidence_score = .data$confidence_score,
      feature_type = "sequence_candidate",
      locus_tag = NA_character_,
      product = .data$query_name,
      evidence_label = .data$evidence_label,
      evidence_summary = .data$evidence_summary,
      source_call_types = "sequence",
      annotation_element_ids = NA_character_,
      sequence_element_ids = .data$sequence_element_id,
      tir_found = .data$tir_found,
      tir_len_bp = .data$tir_len_bp,
      tir_identity = .data$tir_identity,
      tsd_found = .data$tsd_found,
      tsd_len_bp = .data$tsd_len_bp,
      tsd_seq = .data$tsd_seq,
      reference_hit_id = dplyr::coalesce(.data$reference_hit_id, NA_character_),
      reference_hit_family = dplyr::coalesce(.data$reference_hit_family, NA_character_),
      reference_hit_identity = dplyr::coalesce(.data$reference_hit_identity, NA_real_)
    )
}

.dnmb_calls_are_mergeable <- function(current, next_row, merge_gap_bp) {
  if (!identical(current$contig[[1]], next_row$contig[[1]])) {
    return(FALSE)
  }

  cur_fam <- current$element_family[[1]]
  nxt_fam <- next_row$element_family[[1]]
  family_ok <- .dnmb_family_compatible(cur_fam, nxt_fam) ||
    isTRUE(cur_fam == "new") ||
    isTRUE(nxt_fam == "new")

  overlap_or_close <- next_row$start[[1]] <= current$end[[1]] + merge_gap_bp &&
    current$start[[1]] <= next_row$end[[1]] + merge_gap_bp

  source_combo <- current$source_call_types[[1]] != next_row$source_call_types[[1]] ||
    min(current$end[[1]], next_row$end[[1]]) >= max(current$start[[1]], next_row$start[[1]])

  isTRUE(family_ok && overlap_or_close && source_combo)
}

.dnmb_merge_call_pair <- function(current, next_row) {
  current$start[[1]] <- min(current$start[[1]], next_row$start[[1]], na.rm = TRUE)
  current$end[[1]] <- max(current$end[[1]], next_row$end[[1]], na.rm = TRUE)
  current$strand[[1]] <- .dnmb_merge_strand(current$strand[[1]], next_row$strand[[1]])
  current$element_family[[1]] <- .dnmb_merge_family(current$element_family[[1]], next_row$element_family[[1]])
  current$mobile_element_type[[1]] <- .dnmb_merge_unique_values(current$mobile_element_type[[1]], next_row$mobile_element_type[[1]])
  current$confidence_score[[1]] <- max(
    current$confidence_score[[1]],
    next_row$confidence_score[[1]],
    na.rm = TRUE
  ) + as.integer(current$source_call_types[[1]] != next_row$source_call_types[[1]])
  current$feature_type[[1]] <- .dnmb_merge_unique_values(current$feature_type[[1]], next_row$feature_type[[1]])
  current$locus_tag[[1]] <- .dnmb_merge_unique_values(current$locus_tag[[1]], next_row$locus_tag[[1]])
  current$product[[1]] <- .dnmb_merge_unique_values(current$product[[1]], next_row$product[[1]])
  current$evidence_label[[1]] <- .dnmb_merge_unique_values(current$evidence_label[[1]], next_row$evidence_label[[1]])
  current$evidence_summary[[1]] <- .dnmb_merge_unique_values(current$evidence_summary[[1]], next_row$evidence_summary[[1]])
  current$source_call_types[[1]] <- .dnmb_merge_unique_values(current$source_call_types[[1]], next_row$source_call_types[[1]])
  current$annotation_element_ids[[1]] <- .dnmb_merge_unique_values(current$annotation_element_ids[[1]], next_row$annotation_element_ids[[1]])
  current$sequence_element_ids[[1]] <- .dnmb_merge_unique_values(current$sequence_element_ids[[1]], next_row$sequence_element_ids[[1]])

  current$tir_found[[1]] <- isTRUE(current$tir_found[[1]]) || isTRUE(next_row$tir_found[[1]])
  current$tir_len_bp[[1]] <- suppressWarnings(max(current$tir_len_bp[[1]], next_row$tir_len_bp[[1]], na.rm = TRUE))
  if (!is.finite(current$tir_len_bp[[1]])) current$tir_len_bp[[1]] <- NA_integer_
  current$tir_identity[[1]] <- suppressWarnings(max(current$tir_identity[[1]], next_row$tir_identity[[1]], na.rm = TRUE))
  if (!is.finite(current$tir_identity[[1]])) current$tir_identity[[1]] <- NA_real_
  current$tsd_found[[1]] <- isTRUE(current$tsd_found[[1]]) || isTRUE(next_row$tsd_found[[1]])
  current$tsd_len_bp[[1]] <- suppressWarnings(max(current$tsd_len_bp[[1]], next_row$tsd_len_bp[[1]], na.rm = TRUE))
  if (!is.finite(current$tsd_len_bp[[1]])) current$tsd_len_bp[[1]] <- NA_integer_
  current$tsd_seq[[1]] <- .dnmb_merge_unique_values(current$tsd_seq[[1]], next_row$tsd_seq[[1]])
  current$reference_hit_id[[1]] <- .dnmb_merge_unique_values(current$reference_hit_id[[1]], next_row$reference_hit_id[[1]])
  current$reference_hit_family[[1]] <- .dnmb_merge_family(current$reference_hit_family[[1]], next_row$reference_hit_family[[1]])
  current$reference_hit_identity[[1]] <- suppressWarnings(max(current$reference_hit_identity[[1]], next_row$reference_hit_identity[[1]], na.rm = TRUE))
  if (!is.finite(current$reference_hit_identity[[1]])) current$reference_hit_identity[[1]] <- NA_real_

  current
}

.dnmb_merge_unique_values <- function(a, b) {
  vals <- c(a, b)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) {
    return(NA_character_)
  }
  paste(unique(unlist(strsplit(vals, ";\\s*"))), collapse = "; ")
}

.dnmb_add_gene_support <- function(elements, genes) {
  if (!nrow(elements)) {
    return(elements)
  }
  if (!nrow(genes)) {
    elements$supporting_gene_count <- 0L
    elements$supporting_genes <- NA_character_
    elements$supporting_products <- NA_character_
    return(elements)
  }

  support <- lapply(seq_len(nrow(elements)), function(i) {
    overlaps <- genes %>%
      dplyr::filter(
        .data$contig == elements$contig[[i]],
        .data$start <= elements$end[[i]],
        .data$end >= elements$start[[i]]
      )
    tibble::tibble(
      supporting_gene_count = nrow(overlaps),
      supporting_genes = if (nrow(overlaps)) paste(unique(stats::na.omit(overlaps$gene_id)), collapse = "; ") else NA_character_,
      supporting_products = if (nrow(overlaps)) paste(unique(stats::na.omit(overlaps$product)), collapse = "; ") else NA_character_
    )
  })

  dplyr::bind_cols(elements, dplyr::bind_rows(support))
}

.dnmb_score_to_confidence <- function(score) {
  if (is.na(score)) {
    return("low")
  }
  if (score >= 5L) {
    return("high")
  }
  if (score >= 3L) {
    return("medium")
  }
  "low"
}
#' Internal helpers for sequence-based mobileome detection
#'
#' Sequence-engine routines used by the mobileome pipeline to identify and
#' refine mobile element candidates from genome sequence evidence.
#'
#' @name dnmb_internal_mobileome_sequence_engine
#' @keywords internal
#' @noRd
NULL

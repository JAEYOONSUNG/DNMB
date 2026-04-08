.dnmb_write_exec_wrapper <- function(path, command, args_prefix = character()) {
  path <- path.expand(as.character(path)[1])
  command <- path.expand(as.character(command)[1])
  args_prefix <- as.character(args_prefix)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  quoted_args <- if (length(args_prefix)) {
    paste(vapply(args_prefix, shQuote, character(1), type = "sh"), collapse = " ")
  } else {
    ""
  }
  exec_line <- paste(
    "exec",
    shQuote(command, type = "sh"),
    quoted_args,
    "\"$@\""
  )
  script_lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    exec_line
  )
  writeLines(script_lines, con = path)
  Sys.chmod(path, mode = "0755")
  invisible(path)
}

.dnmb_exec_wrapper_target <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  exec_line <- lines[grepl("^exec ", lines)][1]
  if (is.na(exec_line) || !nzchar(exec_line)) {
    return(path.expand(path))
  }
  if (grepl("^exec '[^']+'", exec_line)) {
    return(path.expand(sub("^exec '([^']+)'.*$", "\\1", exec_line)))
  }
  if (grepl('^exec "[^"]+"', exec_line)) {
    return(path.expand(sub('^exec "([^"]+)".*$', "\\1", exec_line)))
  }
  pieces <- strsplit(sub("^exec ", "", exec_line), "[[:space:]]+")[[1]]
  if (!length(pieces)) {
    return(path.expand(path))
  }
  path.expand(gsub("^['\"]|['\"]$", "", pieces[[1]]))
}

.dnmb_exec_wrapper_ok <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  target <- .dnmb_exec_wrapper_target(path)
  if (is.na(target) || !nzchar(target)) {
    return(FALSE)
  }
  if (!file.exists(target)) {
    return(FALSE)
  }
  file.access(target, mode = 1L) == 0L
}

.dnmb_defense_escape_gff_value <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("%", "%25", x, fixed = TRUE)
  x <- gsub(";", "%3B", x, fixed = TRUE)
  x <- gsub("=", "%3D", x, fixed = TRUE)
  x <- gsub(",", "%2C", x, fixed = TRUE)
  x <- gsub("\t|\r|\n", " ", x)
  trimws(x)
}

.dnmb_defense_normalize_nt <- function(x) {
  x <- as.character(x)
  x <- gsub("[[:space:]]+", "", x)
  x <- toupper(x)
  x <- gsub("[^ACGTN]", "N", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

.dnmb_defense_direction_symbol <- function(x) {
  x <- trimws(as.character(x))
  out <- ifelse(
    x %in% c("-", "-1", "reverse", "minus"),
    "-",
    "+"
  )
  out[is.na(x) | !nzchar(x)] <- "+"
  out
}

.dnmb_defense_contig_ids <- function(contig, contig_number = NULL) {
  contig <- as.character(contig)
  contig_number <- as.character(contig_number %||% rep(NA_character_, length(contig)))
  if (length(contig_number) != length(contig)) {
    contig_number <- rep(contig_number, length.out = length(contig))
  }
  key_tbl <- unique(data.frame(contig = contig, contig_number = contig_number, stringsAsFactors = FALSE))
  safe <- ifelse(
    !is.na(key_tbl$contig) & grepl("^[A-Za-z0-9_.-]+$", key_tbl$contig),
    key_tbl$contig,
    NA_character_
  )
  missing_safe <- is.na(safe) | !nzchar(safe)
  if (any(missing_safe)) {
    contig_idx <- suppressWarnings(as.integer(key_tbl$contig_number[missing_safe]))
    safe[missing_safe] <- ifelse(
      !is.na(contig_idx),
      sprintf("DNMB_CONTIG_%03d", contig_idx),
      gsub("[^A-Za-z0-9_.-]+", "_", key_tbl$contig[missing_safe])
    )
  }
  safe[is.na(safe) | !nzchar(safe)] <- sprintf("DNMB_CONTIG_%03d", seq_len(sum(is.na(safe) | !nzchar(safe))))
  if (any(duplicated(safe))) {
    dup_index <- split(seq_along(safe), safe)
    for (nm in names(dup_index)) {
      idx <- dup_index[[nm]]
      if (length(idx) > 1L) {
        safe[idx] <- paste0(safe[idx], "_", seq_along(idx))
      }
    }
  }
  safe[match(paste(contig, contig_number, sep = "\r"), paste(key_tbl$contig, key_tbl$contig_number, sep = "\r"))]
}

.dnmb_defense_query_ids <- function(locus_tag, protein_id = NULL) {
  locus_tag <- .dnmb_module_clean_annotation_key(locus_tag)
  protein_id <- protein_id %||% rep(NA_character_, length(locus_tag))
  protein_id <- as.character(protein_id)
  if (length(protein_id) != length(locus_tag)) {
    protein_id <- rep(protein_id, length.out = length(locus_tag))
  }
  protein_id <- .dnmb_module_clean_annotation_key(protein_id)
  out <- protein_id
  invalid <- is.na(out) | !nzchar(out)
  out[invalid] <- locus_tag[invalid]
  still_invalid <- is.na(out) | !nzchar(out)
  if (any(still_invalid)) {
    out[still_invalid] <- sprintf("DNMB_PROT_%06d", seq_len(sum(still_invalid)))
  }
  duplicated_id <- duplicated(out) | duplicated(out, fromLast = TRUE)
  out[duplicated_id] <- paste0(out[duplicated_id], "__", locus_tag[duplicated_id])
  if (any(duplicated(out))) {
    out <- ave(out, out, FUN = function(values) {
      if (length(values) <= 1L) {
        values
      } else {
        paste0(values, "_", seq_along(values))
      }
    })
  }
  out
}

.dnmb_defense_prepare_cds_genes <- function(genes, require_nt = FALSE) {
  genes <- .dnmb_prepare_query_proteins(genes)
  required_cols <- c("contig", "start", "end", "direction", "product")
  missing_cols <- setdiff(required_cols, names(genes))
  if (length(missing_cols)) {
    stop(
      "Defense modules require genbank_table columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  genes <- as.data.frame(genes, stringsAsFactors = FALSE)
  genes$contig <- as.character(genes$contig)
  genes$contig_id <- .dnmb_defense_contig_ids(genes$contig, genes$contig_number %||% NA_character_)
  genes$product <- as.character(genes$product)
  genes$gene <- as.character(genes$gene %||% NA_character_)
  genes$protein_id <- as.character(genes$protein_id %||% NA_character_)
  genes$start <- suppressWarnings(as.numeric(genes$start))
  genes$end <- suppressWarnings(as.numeric(genes$end))
  genes$direction <- .dnmb_defense_direction_symbol(genes$direction)
  genes$query_id <- .dnmb_defense_query_ids(genes$locus_tag, genes$protein_id)
  genes$translation <- .dnmb_normalize_translation(genes$translation)

  nt_source <- NULL
  if ("rearranged_nt_seq" %in% names(genes)) {
    nt_source <- genes$rearranged_nt_seq
  } else if ("nt_seq" %in% names(genes)) {
    nt_source <- genes$nt_seq
  }
  genes$cds_seq <- .dnmb_defense_normalize_nt(nt_source %||% rep(NA_character_, nrow(genes)))

  keep <- !is.na(genes$contig) &
    nzchar(genes$contig) &
    !is.na(genes$start) &
    !is.na(genes$end) &
    !is.na(genes$translation) &
    nzchar(genes$translation)
  if (isTRUE(require_nt)) {
    keep <- keep & !is.na(genes$cds_seq) & nzchar(genes$cds_seq)
  }
  genes <- genes[keep, , drop = FALSE]
  order_contig <- if ("contig_number" %in% names(genes)) suppressWarnings(as.numeric(genes$contig_number)) else genes$contig_id
  genes <- genes[order(order_contig, genes$start, genes$end, genes$locus_tag), , drop = FALSE]
  rownames(genes) <- NULL
  genes$record_index <- seq_len(nrow(genes))
  genes
}

.dnmb_defense_infer_assembly_label <- function(genes, assembly_id = NULL) {
  assembly_id <- as.character(assembly_id %||% "")[1]
  if (nzchar(trimws(assembly_id))) {
    return(trimws(assembly_id))
  }
  gb_path <- .dnmb_module_detect_genbank(getwd())
  if (!is.null(gb_path) && nzchar(gb_path)) {
    return(tools::file_path_sans_ext(basename(gb_path)))
  }
  contigs <- unique(as.character(genes$contig))
  contigs <- contigs[!is.na(contigs) & nzchar(contigs)]
  if (length(contigs) == 1L) {
    return(contigs[[1L]])
  }
  "DNMB_genbank_table"
}

.dnmb_write_padloc_input <- function(genes,
                                     output_dir,
                                     prefix = "padloc_query") {
  cds <- .dnmb_defense_prepare_cds_genes(genes, require_nt = FALSE)
  faa_path <- file.path(output_dir, paste0(prefix, ".faa"))
  gff_path <- file.path(output_dir, paste0(prefix, ".gff"))
  map_path <- file.path(output_dir, paste0(prefix, "_id_map.tsv"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!nrow(cds)) {
    writeLines(character(), con = faa_path)
    writeLines("##gff-version 3", con = gff_path)
    utils::write.table(cds, file = map_path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(list(faa = faa_path, gff = gff_path, map = cds, map_path = map_path))
  }

  faa_tbl <- data.frame(
    protein_label = cds$query_id,
    protein_seq = cds$translation,
    stringsAsFactors = FALSE
  )
  .dnmb_write_protein_fasta(faa_tbl, faa_path)

  gff_lines <- c("##gff-version 3")
  for (i in seq_len(nrow(cds))) {
    attrs <- c(
      paste0("ID=cds-", .dnmb_defense_escape_gff_value(cds$query_id[[i]])),
      paste0("Name=", .dnmb_defense_escape_gff_value(cds$query_id[[i]])),
      paste0("locus_tag=", .dnmb_defense_escape_gff_value(cds$locus_tag[[i]])),
      paste0("protein_id=", .dnmb_defense_escape_gff_value(cds$query_id[[i]])),
      paste0("product=", .dnmb_defense_escape_gff_value(cds$product[[i]]))
    )
    gff_lines <- c(
      gff_lines,
      paste(
        cds$contig_id[[i]],
        "DNMB",
        "CDS",
        as.integer(cds$start[[i]]),
        as.integer(cds$end[[i]]),
        ".",
        cds$direction[[i]],
        "0",
        paste(attrs, collapse = ";"),
        sep = "\t"
      )
    )
  }
  writeLines(gff_lines, con = gff_path)

  map_tbl <- cds[, c("query_id", "locus_tag", "protein_id", "contig", "contig_id", "start", "end", "direction", "product"), drop = FALSE]
  utils::write.table(map_tbl, file = map_path, sep = "\t", row.names = FALSE, quote = FALSE)
  list(faa = faa_path, gff = gff_path, map = map_tbl, map_path = map_path)
}

.dnmb_write_defensepredictor_input <- function(genes,
                                               output_dir,
                                               prefix = "defense_predictor_query",
                                               assembly_id = NULL) {
  cds <- .dnmb_defense_prepare_cds_genes(genes, require_nt = TRUE)
  ft_path <- file.path(output_dir, paste0(prefix, "_feature_table.tsv"))
  fna_path <- file.path(output_dir, paste0(prefix, "_cds_from_genomic.fna"))
  faa_path <- file.path(output_dir, paste0(prefix, "_protein.faa"))
  map_path <- file.path(output_dir, paste0(prefix, "_id_map.tsv"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!nrow(cds)) {
    writeLines(character(), con = ft_path)
    writeLines(character(), con = fna_path)
    writeLines(character(), con = faa_path)
    utils::write.table(cds, file = map_path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(list(feature_table = ft_path, cds_fna = fna_path, protein_faa = faa_path, map = cds, map_path = map_path))
  }

  assembly_label <- .dnmb_defense_infer_assembly_label(cds, assembly_id = assembly_id)
  feature_tbl <- data.frame(
    "# feature" = "CDS",
    class = "with_protein",
    assembly = assembly_label,
    assembly_unit = if (length(unique(cds$contig_id)) == 1L) "Primary Assembly" else "GenBank Assembly",
    seq_type = if (length(unique(cds$contig_id)) == 1L) "chromosome" else "contig",
    chromosome = if (length(unique(cds$contig)) == 1L) unique(cds$contig)[1] else "",
    genomic_accession = cds$contig_id,
    start = as.integer(cds$start),
    end = as.integer(cds$end),
    strand = cds$direction,
    product_accession = cds$query_id,
    `non-redundant_refseq` = cds$query_id,
    related_accession = "",
    name = cds$product,
    symbol = if ("gene" %in% names(cds)) cds$gene else NA_character_,
    GeneID = "",
    locus_tag = cds$locus_tag,
    feature_interval_length = abs(as.integer(cds$end) - as.integer(cds$start)) + 1L,
    product_length = nchar(cds$translation),
    attributes = paste0(
      "protein_id=", cds$query_id,
      ";locus_tag=", cds$locus_tag,
      ";product=", cds$product
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    feature_tbl,
    file = ft_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  faa_tbl <- data.frame(
    protein_label = cds$query_id,
    protein_seq = cds$translation,
    stringsAsFactors = FALSE
  )
  .dnmb_write_protein_fasta(faa_tbl, faa_path)

  fna_con <- file(fna_path, open = "w")
  on.exit(close(fna_con), add = TRUE)
  for (i in seq_len(nrow(cds))) {
    loc <- if (identical(cds$direction[[i]], "-")) {
      sprintf("complement(%s..%s)", as.integer(cds$start[[i]]), as.integer(cds$end[[i]]))
    } else {
      sprintf("%s..%s", as.integer(cds$start[[i]]), as.integer(cds$end[[i]]))
    }
    header <- paste(
      paste0("lcl|", cds$contig_id[[i]], "_cds_", cds$query_id[[i]], "_", cds$record_index[[i]]),
      paste0("[protein_id=", cds$query_id[[i]], "]"),
      paste0("[locus_tag=", cds$locus_tag[[i]], "]"),
      paste0("[location=", loc, "]"),
      paste0("[gbkey=CDS]"),
      sep = " "
    )
    writeLines(paste0(">", header), con = fna_con)
    writeLines(cds$cds_seq[[i]], con = fna_con)
  }

  map_tbl <- cds[, c("query_id", "locus_tag", "protein_id", "contig", "contig_id", "start", "end", "direction", "product"), drop = FALSE]
  utils::write.table(map_tbl, file = map_path, sep = "\t", row.names = FALSE, quote = FALSE)
  list(
    feature_table = ft_path,
    cds_fna = fna_path,
    protein_faa = faa_path,
    map = map_tbl,
    map_path = map_path
  )
}

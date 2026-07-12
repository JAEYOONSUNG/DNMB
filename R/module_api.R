#' Merge module run outputs by locus tag
#'
#' Combines multiple `run_module()` / `run_module_set()` outputs into a single
#' locus-level table. Shared backbone columns are kept once, while module-
#' specific columns are prefixed with the module or database name so multiple
#' module outputs can coexist in one table.
#'
#' @param ... `dnmb_module_run` objects to merge, or a single named list of
#'   module results such as the output from `run_module_set()`.
#'
#' @return A data frame keyed by `locus_tag`.
#' @export
merge_module_results <- function(...) {
  module_results <- dnmb_collect_module_runs(...)

  if (length(module_results) == 0L) {
    stop("No module results supplied.", call. = FALSE)
  }

  module_tables <- Map(
    dnmb_module_table_with_name,
    module_result = module_results,
    fallback_name = names(module_results)
  )
  merged <- module_tables[[1L]]$table

  if (length(module_tables) > 1L) {
    for (index in seq.int(2L, length(module_tables))) {
      next_table <- module_tables[[index]]$table
      merged <- dnmb_merge_two_module_tables(merged, next_table)
    }
  }

  dnmb_finalize_merged_module_table(merged)
}

#' Append merged module results onto a DNMB table
#'
#' Joins wide module results onto an existing DNMB or GenBank table by
#' `locus_tag`, preserving the original row order and keeping existing DNMB
#' columns unchanged.
#'
#' @param main_table A data frame containing the main DNMB or GenBank table.
#' @param ... `dnmb_module_run` objects to merge before appending, or a single
#'   named list of module results such as the output from `run_module_set()`.
#' @param merged_module_table An optional pre-merged module-wide table keyed by
#'   `locus_tag`.
#'
#' @return `main_table` with new module-prefixed columns appended.
#' @export
append_module_results <- function(main_table, ..., merged_module_table = NULL) {
  dnmb_validate_append_main_table(main_table)

  module_results <- dnmb_collect_module_runs(...)

  if (!is.null(merged_module_table) && length(module_results) > 0L) {
    stop(
      "Supply module runs or merged_module_table, not both.",
      call. = FALSE
    )
  }

  if (is.null(merged_module_table)) {
    if (length(module_results) == 0L) {
      return(main_table)
    }
    merged_module_table <- do.call(merge_module_results, module_results)
  }

  if (!is.data.frame(merged_module_table)) {
    stop("merged_module_table must be a data frame.", call. = FALSE)
  }

  merged_module_table <- dnmb_normalize_prophage_prefix(merged_module_table)

  if (!"locus_tag" %in% names(merged_module_table)) {
    stop("merged_module_table must contain a 'locus_tag' column.", call. = FALSE)
  }

  module_columns <- setdiff(names(merged_module_table), dnmb_backbone_columns())

  if (length(module_columns) == 0L) {
    return(main_table)
  }

  overlapping_columns <- intersect(module_columns, names(main_table))

  if (length(overlapping_columns) > 0L) {
    # Drop overlapping columns from main_table so module results take precedence
    # (common in re-runs where modules were already appended)
    warning(
      sprintf(
        "Overwriting existing module columns: %s",
        paste(overlapping_columns, collapse = ", ")
      ),
      call. = FALSE
    )
    main_table[overlapping_columns] <- NULL
  }

  matched_rows <- match(main_table$locus_tag, merged_module_table$locus_tag)
  appended_columns <- merged_module_table[matched_rows, module_columns, drop = FALSE]

  dbcan_row_column <- grep(
    "(^|_)dbcan_gene_row$",
    names(merged_module_table),
    value = TRUE,
    ignore.case = TRUE
  )
  dbcan_columns <- grep("^dbCAN_", module_columns, value = TRUE, ignore.case = TRUE)
  if (length(dbcan_row_column) && length(dbcan_columns)) {
    dbcan_rows <- suppressWarnings(as.integer(merged_module_table[[dbcan_row_column[[1]]]]))
    row_matched <- match(seq_len(nrow(main_table)), dbcan_rows)
    if (!dnmb_append_row_mapping_is_compatible(main_table, merged_module_table, row_matched)) {
      row_matched <- dnmb_stable_row_match(main_table, merged_module_table)
    }
    unresolved <- is.na(row_matched)
    contig_column <- grep(
      "^dbCAN_dbcan_contig_key$",
      names(merged_module_table),
      value = TRUE,
      ignore.case = TRUE
    )
    if (any(unresolved) && length(contig_column) &&
        all(c("start", "end", "direction") %in% names(main_table)) &&
        all(c("start", "end", "direction") %in% names(merged_module_table))) {
      main_contigs <- .dnmb_dbcan_contig_keys(main_table)$key
      module_contigs <- as.character(merged_module_table[[contig_column[[1]]]])
      main_key <- paste(main_contigs, main_table$start, main_table$end, main_table$direction, sep = "|")
      module_key <- paste(module_contigs, merged_module_table$start, merged_module_table$end, merged_module_table$direction, sep = "|")
      coordinate_match <- dnmb_unique_value_match(main_key, module_key)
      row_matched[unresolved] <- coordinate_match[unresolved]
    }
    include <- !is.na(row_matched)
    for (column_name in dbcan_columns) {
      value <- .dnmb_na_vector_like(merged_module_table[[column_name]], nrow(main_table))
      value[include] <- merged_module_table[[column_name]][row_matched[include]]
      appended_columns[[column_name]] <- value
    }
  }
  rownames(appended_columns) <- NULL

  cbind(main_table, appended_columns, stringsAsFactors = FALSE)
}

dnmb_append_row_mapping_is_compatible <- function(main_table, module_table, matched_rows) {
  mapped <- !is.na(matched_rows)
  if (!any(mapped)) return(FALSE)
  shared <- intersect(names(main_table), names(module_table))
  compare_field <- function(column_name) {
    left <- as.character(main_table[[column_name]][mapped])
    right <- as.character(module_table[[column_name]][matched_rows[mapped]])
    comparable <- !is.na(left) & nzchar(left) & !is.na(right) & nzchar(right)
    list(n = sum(comparable), same = !any(left[comparable] != right[comparable]))
  }
  stable_fields <- intersect(c("protein_id", "contig", "start", "end", "direction"), shared)
  if (length(stable_fields) >= 2L) {
    status <- lapply(stable_fields, compare_field)
    comparable <- sum(vapply(status, `[[`, integer(1), "n"))
    if (comparable > 0L) return(all(vapply(status, `[[`, logical(1), "same")))
  }
  if ("locus_tag" %in% shared) return(compare_field("locus_tag")$same)
  !length(shared)
}

dnmb_validate_append_main_table <- function(main_table) {
  if (!is.data.frame(main_table)) {
    stop("main_table must be a data frame.", call. = FALSE)
  }

  if (!"locus_tag" %in% names(main_table)) {
    stop("main_table must contain a 'locus_tag' column.", call. = FALSE)
  }
}

dnmb_merge_two_module_tables <- function(primary_table, secondary_table) {
  if (dnmb_module_rows_are_aligned(primary_table, secondary_table)) {
    primary_aligned <- primary_table
    secondary_aligned <- secondary_table
  } else {
    stable_rows <- dnmb_stable_row_match(primary_table, secondary_table)
    if (nrow(primary_table) == nrow(secondary_table) && all(!is.na(stable_rows))) {
      primary_aligned <- primary_table
      secondary_aligned <- secondary_table[stable_rows, , drop = FALSE]
      rownames(secondary_aligned) <- NULL
    } else {
      combined_locus_tags <- c(
        primary_table$locus_tag,
        secondary_table$locus_tag[!secondary_table$locus_tag %in% primary_table$locus_tag]
      )
      primary_aligned <- dnmb_align_module_table(primary_table, combined_locus_tags)
      secondary_aligned <- dnmb_align_module_table(secondary_table, combined_locus_tags)
    }
  }
  backbone_columns <- setdiff(dnmb_backbone_columns(), "locus_tag")

  for (column_name in backbone_columns) {
    if (!column_name %in% names(secondary_aligned)) {
      next
    }

    if (!column_name %in% names(primary_aligned)) {
      primary_aligned[[column_name]] <- secondary_aligned[[column_name]]
      next
    }

    primary_aligned[[column_name]] <- dnmb_coalesce_vectors(
      primary_aligned[[column_name]],
      secondary_aligned[[column_name]]
    )
  }

  secondary_specific <- setdiff(names(secondary_aligned), dnmb_backbone_columns())

  for (column_name in secondary_specific) {
    primary_aligned[[column_name]] <- secondary_aligned[[column_name]]
  }

  primary_aligned
}

dnmb_module_rows_are_aligned <- function(primary_table, secondary_table) {
  if (nrow(primary_table) != nrow(secondary_table) || !nrow(primary_table)) {
    return(FALSE)
  }
  shared <- intersect(
    dnmb_backbone_columns(),
    intersect(names(primary_table), names(secondary_table))
  )
  if (!length(shared)) return(FALSE)
  same_column <- function(column_name) {
    identical(
      as.character(primary_table[[column_name]]),
      as.character(secondary_table[[column_name]])
    )
  }
  stable_fields <- intersect(
    c("protein_id", "contig", "start", "end", "direction"),
    shared
  )
  if (length(stable_fields) >= 2L) {
    return(all(vapply(stable_fields, same_column, logical(1))))
  }
  "locus_tag" %in% shared && same_column("locus_tag")
}

dnmb_unique_value_match <- function(target, source) {
  target <- as.character(target)
  source <- as.character(source)
  target_missing <- is.na(target) | !nzchar(target)
  source_missing <- is.na(source) | !nzchar(source)
  target_unique <- !duplicated(target) & !duplicated(target, fromLast = TRUE) & !target_missing
  source_unique <- !duplicated(source) & !duplicated(source, fromLast = TRUE) & !source_missing
  out <- match(target, ifelse(source_unique, source, NA_character_))
  out[!target_unique] <- NA_integer_
  out
}

dnmb_stable_row_match <- function(target_table, source_table) {
  out <- rep(NA_integer_, nrow(target_table))
  if (!nrow(target_table) || !nrow(source_table)) return(out)

  if ("protein_id" %in% names(target_table) && "protein_id" %in% names(source_table)) {
    out <- dnmb_unique_value_match(target_table$protein_id, source_table$protein_id)
  }

  coordinate_fields <- intersect(
    c("contig", "start", "end", "direction"),
    intersect(names(target_table), names(source_table))
  )
  if (all(c("start", "end") %in% coordinate_fields)) {
    target_key <- do.call(paste, c(target_table[coordinate_fields], sep = "|"))
    source_key <- do.call(paste, c(source_table[coordinate_fields], sep = "|"))
    coordinate_match <- dnmb_unique_value_match(target_key, source_key)
    unresolved <- is.na(out)
    out[unresolved] <- coordinate_match[unresolved]
  }

  if ("locus_tag" %in% names(target_table) && "locus_tag" %in% names(source_table)) {
    locus_match <- dnmb_unique_value_match(target_table$locus_tag, source_table$locus_tag)
    unresolved <- is.na(out)
    out[unresolved] <- locus_match[unresolved]
  }
  out
}

dnmb_align_module_table <- function(module_table, locus_tags) {
  matched_rows <- match(locus_tags, module_table$locus_tag)
  aligned_table <- module_table[matched_rows, , drop = FALSE]
  aligned_table$locus_tag <- locus_tags
  rownames(aligned_table) <- NULL
  aligned_table
}

dnmb_collect_module_runs <- function(...) {
  module_results <- list(...)

  if (length(module_results) == 1L &&
      is.list(module_results[[1L]]) &&
      !inherits(module_results[[1L]], "dnmb_module_run") &&
      !is.data.frame(module_results[[1L]])) {
    module_results <- module_results[[1L]]
  }

  if (is.null(names(module_results))) {
    names(module_results) <- rep.int("", length(module_results))
  }

  module_results
}

dnmb_module_table_with_name <- function(module_result, fallback_name = "") {
  module_name <- dnmb_module_name(module_result, fallback_name)
  module_table <- dnmb_module_table(module_result)

  if (!"locus_tag" %in% names(module_table)) {
    stop(
      sprintf("Module '%s' result must contain a 'locus_tag' column.", module_name),
      call. = FALSE
    )
  }

  list(
    name = module_name,
    table = dnmb_prefix_module_columns(module_table, module_name)
  )
}

dnmb_module_name <- function(module_result, fallback_name = "") {
  candidates <- c(
    fallback_name,
    attr(module_result, "database", exact = TRUE),
    attr(module_result, "module", exact = TRUE),
    attr(module_result, "name", exact = TRUE),
    if (is.list(module_result)) module_result[["database"]] else NULL,
    if (is.list(module_result)) module_result[["module"]] else NULL,
    if (is.list(module_result)) module_result[["name"]] else NULL
  )

  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  candidates <- unlist(candidates, use.names = FALSE)
  candidates <- candidates[nzchar(candidates)]

  if (length(candidates) == 0L) {
    stop(
      "Each module result must include a module/database name.",
      call. = FALSE
    )
  }

  dnmb_safe_module_name(candidates[[1L]])
}

dnmb_module_table <- function(module_result) {
  if (is.data.frame(module_result)) {
    return(module_result)
  }

  if (!is.list(module_result)) {
    stop("Module results must be data frames or dnmb_module_run objects.", call. = FALSE)
  }

  table_candidates <- c("output_table", "table", "result", "results")

  for (candidate in table_candidates) {
    module_table <- module_result[[candidate]]
    if (is.data.frame(module_table)) {
      return(module_table)
    }
  }

  stop("Could not locate an output table in a module result.", call. = FALSE)
}

dnmb_prefix_module_columns <- function(module_table, module_name) {
  backbone_columns <- dnmb_backbone_columns()
  module_specific <- setdiff(names(module_table), backbone_columns)

  if (length(module_specific) > 0L) {
    names(module_table)[match(module_specific, names(module_table))] <- paste0(
      module_name,
      "_",
      module_specific
    )
  }

  module_table
}

dnmb_backbone_columns <- function() {
  c("locus_tag", "gene", "product", "protein_id", "contig", "start", "end", "direction")
}

dnmb_finalize_merged_module_table <- function(merged_table) {
  backbone_columns <- dnmb_backbone_columns()
  ordered_backbone <- intersect(backbone_columns, names(merged_table))
  ordered_specific <- setdiff(names(merged_table), ordered_backbone)
  merged_table[c(ordered_backbone, ordered_specific)]
}

dnmb_coalesce_vectors <- function(primary, secondary) {
  if (length(primary) == 0L) {
    return(secondary)
  }

  missing_primary <- is.na(primary)

  if (is.character(primary)) {
    missing_primary <- missing_primary | primary == ""
  }

  primary[missing_primary] <- secondary[missing_primary]
  primary
}

dnmb_safe_module_name <- function(module_name) {
  gsub("[^[:alnum:]_]+", "_", module_name)
}

dnmb_normalize_prophage_prefix <- function(merged_table) {
  if (!is.data.frame(merged_table) || !ncol(merged_table)) {
    return(merged_table)
  }
  nm <- names(merged_table)
  phispy_cols <- grep("^PhiSpy_", nm, value = TRUE)
  if (!length(phispy_cols)) {
    return(merged_table)
  }
  renamed <- sub("^PhiSpy_", "Prophage_", phispy_cols)
  for (idx in seq_along(phispy_cols)) {
    old <- phispy_cols[[idx]]
    new <- renamed[[idx]]
    if (new %in% names(merged_table)) {
      merged_table[[new]] <- NULL
    }
    names(merged_table)[names(merged_table) == old] <- new
  }
  merged_table
}

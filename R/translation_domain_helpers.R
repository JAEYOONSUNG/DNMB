.dnmb_detect_genbank_path <- function(path = getwd()) {
  candidates <- list.files(
    path,
    pattern = "\\.gbk$|\\.gb$|\\.gbff$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  candidates <- candidates[!grepl("(^|/|\\\\)~\\$", candidates)]
  if (!length(candidates)) {
    return(NULL)
  }
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

.dnmb_read_genbank_lines <- function(gb_path = NULL) {
  gb_path <- gb_path %||% .dnmb_detect_genbank_path()
  gb_path <- trimws(as.character(gb_path)[1])
  if (is.na(gb_path) || !nzchar(gb_path) || !file.exists(gb_path)) {
    return(character())
  }
  tryCatch(readLines(gb_path, warn = FALSE), error = function(e) character())
}

.dnmb_extract_genbank_lineage <- function(gb_path = NULL) {
  lines <- .dnmb_read_genbank_lines(gb_path)
  if (!length(lines)) {
    return(NA_character_)
  }
  org_idx <- grep("^  ORGANISM", lines)
  if (!length(org_idx)) {
    return(NA_character_)
  }
  i <- org_idx[[1]] + 1L
  lineage_lines <- character()
  while (i <= length(lines) && grepl("^            ", lines[[i]])) {
    lineage_lines <- c(lineage_lines, trimws(lines[[i]]))
    i <- i + 1L
  }
  lineage <- paste(lineage_lines, collapse = " ")
  lineage <- gsub("\\s+", " ", lineage)
  lineage <- trimws(lineage)
  if (!nzchar(lineage)) {
    return(NA_character_)
  }
  lineage
}

.dnmb_detect_translation_domain <- function(target = NULL,
                                            gb_path = NULL,
                                            fallback = "bacteria") {
  lineage <- .dnmb_extract_genbank_lineage(gb_path)
  lineage_low <- tolower(lineage %||% "")
  if (nzchar(lineage_low)) {
    if (grepl("archaea", lineage_low, fixed = TRUE)) {
      return("archaea")
    }
    if (grepl("bacteria", lineage_low, fixed = TRUE)) {
      return("bacteria")
    }
    if (grepl("viruses", lineage_low, fixed = TRUE)) {
      return("phage")
    }
  }

  product_col <- if (!is.null(target) && is.data.frame(target) && "product" %in% names(target)) {
    tolower(paste(target$product[!is.na(target$product)], collapse = " "))
  } else {
    ""
  }
  if (grepl("16s ribosomal rna|small subunit ribosomal rna|ssu rrna", product_col)) {
    return("bacteria")
  }
  fallback
}

.dnmb_detect_transl_table <- function(target = NULL,
                                      gb_path = NULL,
                                      translation_domain = NULL) {
  translation_domain <- translation_domain %||% .dnmb_detect_translation_domain(target = target, gb_path = gb_path)
  lines <- .dnmb_read_genbank_lines(gb_path)
  transl_table_vals <- character()
  if (length(lines)) {
    transl_lines <- grep("/transl_table=", lines, value = TRUE)
    if (length(transl_lines)) {
      transl_table_vals <- stringr::str_extract(transl_lines, "[0-9]+")
      transl_table_vals <- transl_table_vals[!is.na(transl_table_vals) & nzchar(transl_table_vals)]
    }
  }
  if (length(transl_table_vals)) {
    tt <- suppressWarnings(as.integer(names(sort(table(transl_table_vals), decreasing = TRUE))[1]))
    if (!is.na(tt)) {
      return(tt)
    }
  }
  if (identical(translation_domain, "archaea")) {
    return(11L)
  }
  11L
}

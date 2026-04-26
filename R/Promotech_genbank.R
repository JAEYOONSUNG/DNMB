#' Add Promotech promoter predictions to GenBank
#'
#' @param promotech Path to a Promotech `genome_predictions.csv`,
#'   `promotech_predictions.tsv`, or legacy `promotech.csv` file. If `NULL`,
#'   the current directory is searched.
#' @param genbank Path to a GenBank file. If `NULL`, the current directory is
#'   searched for one `.gb`, `.gbk`, or `.gbff` file.
#' @param output_file Optional GenBank path for the annotated copy.
#' @param snippet_file Optional path for the paste-ready GenBank FEATURES block.
#' @param threshold Optional minimum Promotech score to keep.
#' @return Invisibly returns a list with the written annotation file paths.
#' @export
#'

Promoter_to_genbank <- function(promotech = NULL,
                                genbank = NULL,
                                output_file = NULL,
                                snippet_file = NULL,
                                threshold = NULL) {
  current_dir <- getwd()

  if (is.null(promotech)) {
    promotech_files <- list.files(
      current_dir,
      pattern = "(genome_predictions\\.csv|promotech_predictions\\.tsv|promotech\\.csv)$",
      full.names = TRUE
    )
    if (length(promotech_files) == 0) {
      stop("No Promotech prediction file found in the current directory.", call. = FALSE)
    } else if (length(promotech_files) > 1) {
      stop("Multiple Promotech prediction files found. Please specify `promotech`.", call. = FALSE)
    } else {
      promotech <- promotech_files[1]
    }
  } else {
    if (!file.exists(promotech)) {
      stop(paste("The specified promotech file does not exist:", promotech), call. = FALSE)
    }
  }

  if (is.null(genbank)) {
    gb_files <- list.files(current_dir, pattern = "\\.(gb|gbk|gbff)$", full.names = TRUE)
    if (length(gb_files) == 0) {
      stop("No GenBank files ending with .gb, .gbk, or .gbff found in the current directory.", call. = FALSE)
    } else if (length(gb_files) > 1) {
      stop("Multiple GenBank files found. Please specify `genbank`.", call. = FALSE)
    } else {
      genbank <- gb_files[1]
    }
  } else {
    if (!file.exists(genbank)) {
      stop(paste("The specified GenBank file does not exist:", genbank), call. = FALSE)
    }
  }

  if (is.null(output_file)) {
    base_name <- tools::file_path_sans_ext(genbank)
    extension <- tools::file_ext(genbank)
    output_file <- paste0(base_name, "_promotech.", extension)
  }
  if (is.null(snippet_file)) {
    snippet_file <- file.path(dirname(output_file), "promoter_feature_for_gb")
  }

  predictions <- .dnmb_promotech_read_predictions(promotech, threshold = threshold)
  artifacts <- .dnmb_promotech_write_genbank_artifacts(
    predictions = predictions,
    hits = NULL,
    genbank = genbank,
    output_dir = dirname(output_file),
    snippet_path = snippet_file,
    annotated_genbank_path = output_file
  )

  if (!file.exists(output_file)) {
    stop(.dnmb_module_status_detail(artifacts$status) %||% "Failed to write Promotech GenBank annotations.", call. = FALSE)
  }
  message(paste("Modified GenBank file has been saved as", output_file))
  invisible(artifacts$files)
}

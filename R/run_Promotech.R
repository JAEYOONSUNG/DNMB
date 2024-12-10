#' Save Contig Sequences and Run External Tool
#'
#' This script identifies objects in the R environment that match the pattern `contig_[number]_seq`,
#' saves them as FASTA files, and runs an external Python tool on each file using a conda environment.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Identifies all objects in the R environment matching the pattern `contig_[number]_seq`.
#'   \item Saves each matching object as a `.fna` file in FASTA format.
#'   \item Executes the external Python tool `promotech.py` on each `.fna` file using a conda environment.
#'   \item Outputs the results into a designated `results` folder.
#' }
#'
#' @param conda_env_path A character string specifying the path to the conda environment
#'                       (e.g., `/Users/JaeYoon/miniforge3/envs/promotech_mac_env`).
#' @param promotech_script_path A character string specifying the path to the `promotech.py` script.
#' @param output_dir A character string specifying the directory to store the results. Default is `results`.
#'
#' @return This function does not return a value but saves `.fna` files and runs `promotech.py`.
#'         Outputs from the external tool are stored in the specified results folder.
#'
#' @examples
#' # Example usage:
#' save_and_run_promotech(
#'   conda_env_path = "/Users/JaeYoon/miniforge3/envs/promotech_mac_env",
#'   promotech_script_path = "~/Desktop/Bioinformatics/promotech-master/promotech.py"
#' )
#'
#' @export

save_and_run_promotech <- function(conda_env_path, promotech_script_path, output_dir = "results") {
  # Ensure required packages are available
  if (!requireNamespace("tools", quietly = TRUE)) {
    stop("Package 'tools' is required but not installed.")
  }

  # Get the list of all objects in the environment
  object_names <- ls()

  # Filter for contig_[number]_seq objects
  contig_objects <- grep("^contig_[0-9]+_seq$", object_names, value = TRUE)

  # Loop over each contig object and save it in a .fna file
  for (contig in contig_objects) {
    # Get the sequence data from the object
    seq_data <- get(contig)

    # Create the filename
    file_name <- paste0(contig, ".fna")

    # Write to file in FASTA format
    writeLines(c(paste0(">", contig), seq_data), file_name)
    message(paste("Saved", file_name))
  }

  # Get the list of all .fna files
  fna_files <- list.files(pattern = "\\.fna$")

  # Run promotech.py using the conda environment
  for (file in fna_files) {
    command <- paste(
      "conda run -p", normalizePath(conda_env_path),
      "python", normalizePath(promotech_script_path),
      "-pg -m RF-HOT -f", shQuote(normalizePath(file)),
      "-g -o", shQuote(output_dir)
    )
    message(paste("Running command:", command))
    system(command)
  }

  message("All files processed.")
}

## Get the list of all objects in the environment
#object_names <- ls()
#
## Filter for contig_[number]_seq objects
#contig_objects <- grep("^contig_[0-9]+_seq$", object_names, value = TRUE)
## Loop over each contig object and save it in a .fna file
#for (contig in contig_objects) {
#  # Get the sequence data from the object
#  seq_data <- get(contig)
#  # Create the filename
#  file_name <- paste0(contig, ".fna")
#  # Write to file in FASTA format
#  writeLines(c(paste0(">", contig), seq_data), file_name)
#  message(paste("Saved", file_name))
#}
#
#
#fna_files <- list.files(pattern = "\\.fna$")
#
## Loop over the files
#for (i in seq_along(fna_files)) {
#  command <- paste("conda run -n promotech_mac_env python promotech.py -pg -m RF-HOT -f '", normalizePath(fna_files[i]), "' -g -o results", sep = "")
#  system(command)
#}
#
## Loop over the files
#for (i in seq_along(fna_files)) {
#  command <- paste("conda run -p /Users/JaeYoon/miniforge3/envs/promotech_mac_env python ~/Desktop/Bioinformatics/promotech-master/promotech.py -pg -m RF-HOT -f '", normalizePath(fna_files[i]), "' -g -o results", sep = "")
#  system(command)
#}
#

#' Genbank fna extractor
#'
#' @param df A data frame containing the data from the genbank.
#' @param gb_dir A string specifying the directory where the genbank output files are located. If NULL, the current working directory is used.
#' @return A data frame containing the baseline table.
#' @export
#'

genbank_fna_extractor <- function(gb_dir = NULL) {

  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Biostrings)
  library(plyr)

  # If gb_dir is NULL, use the current working directory
  if (is.null(gb_dir)) {
    gb_dir <- getwd()
  }

  # Get the list of .tsv files in the directory
  gb_files <- list.files(
    gb_dir,
    pattern = "\\.gbk$|\\.gb$|\\.gbff$",
    full.names = FALSE
  )
  # Exclude temporary files starting with '~$'
  gb_files <- gb_files[!grepl("^~\\$", gb_files)]

  # Read the GenBank file
  origin <- read.csv(gb_files, header = FALSE, sep = "\n")

  # Identify the rows where contigs and their lengths start and end
  contig_row <- grep("^LOCUS", origin$V1)
  pre_contig_start <- grep("^LOCUS", origin$V1)
  pre_contig_finish <- c(pre_contig_start[-1] - 1, nrow(origin))

  # Extract contigs based on row positions
  for (i in 1:length(contig_row)) {
    assign(paste("contig", i, sep = "_"), as.data.frame(origin[pre_contig_start[i]:pre_contig_finish[i], ]))
  }

  # Calculate contig lengths
  contig_length <- as.data.frame(grep("source + 1..", origin$V1, ignore.case = TRUE, value = TRUE))
  colnames(contig_length) <- "contig_length"
  contig_length <- separate(contig_length, col = `contig_length`, into = c("start", "end"), sep = "\\.\\.")
  contig_length <- as.data.frame(contig_length$end)
  colnames(contig_length) <- "length"
  genome_size <- sum(as.numeric(as.character(contig_length$length)))

  assign("contig_length", contig_length, envir = .GlobalEnv)
  assign("genome_size", genome_size, envir = .GlobalEnv)

  # Extract contig definitions
  for (i in 1:nrow(contig_length)) {
    assign(paste("contig", i, "definition", sep = "_"),
           gsub("^DEFINITION  |\\.$", "", get(paste("contig", i, sep = "_"))[grep("DEFINITION", get(paste("contig", i, sep = "_"))[,1], value = TRUE), ]))
  }

  # Sequence extraction
  for (i in 1:nrow(contig_length)) {
    assign(paste("contig", i, "seq_start", sep = "_"), grep("ORIGIN", get(paste("contig", i, sep = "_"))[, 1], value = FALSE))
    assign(paste("contig", i, "seq_finish", sep = "_"), grep("^//", get(paste("contig", i, sep = "_"))[, 1], value = FALSE))

    assign(paste("contig", i, "nt_seq", sep = "_"),
           paste(get(paste("contig", i, sep = "_"))[get(paste("contig", i, "seq_start", sep = "_")):get(paste("contig", i, "seq_finish", sep = "_")), 1], sep = "", collapse = ""))

    assign(paste("contig", i, "seq_trim", sep = "_"), gsub("[0-9]{1,}", "", gsub("ORIGIN| 1 +.{1,}", "", get(paste("contig", i, "nt_seq", sep = "_")))))
    assign(paste("contig", i, "seq_final", sep = "_"), gsub(" ", "", get(paste("contig", i, "seq_trim", sep = "_"))))
  }

  # Convert nucleotide sequences to uppercase
  for (i in 1:nrow(contig_length)) {
    assign(paste("contig", i, "seq", "final", sep = "_"), toupper(get(paste("contig", i, "seq", "final", sep = "_"))))
  }

  # Create final sequence data frame for each contig
  for (i in 1:nrow(contig_length)) {
    assign(paste("contig", i, "seq", "final", "temp", sep = "_"),
           data.frame("definition" = get(paste("contig", i, "definition", sep = "_")),
                      "top_strand" = get(paste("contig", i, "seq", "final", sep = "_"))))
  }

  # Combine contig sequences into one final table
  contig_seq_final <- contig_1_seq_final_temp
  for (i in (1:nrow(contig_length))[-1]) {
    contig_temp <- get(paste("contig", i, "seq", "final", "temp", sep = "_"))
    if (!is.na(contig_temp)[1]) {
      contig_seq_final <- rbind.fill(contig_seq_final, contig_temp)
    }
  }

  #Concatenated nucleotide
  #Need to fix system bug for Large character# # FIXME
  for(i in 1:nrow(contig_length)){
    assign(paste("contig",i,"nt_seq",sep="_"),
           paste(((get(paste("contig",i,sep = "_"))[,1])[get(paste("contig",i,"seq_start",sep = "_")):get(paste("contig",i,"seq_finish",sep = "_"))]),sep="",collapse=""))
  }
  for(i in 1:nrow(contig_length)){
    assign(paste("contig",i,"seq_pre",sep="_"), gsub("ORIGIN","", get(paste("contig",i,"nt_seq",sep="_"))))
  }
  for(i in 1:nrow(contig_length)){
    assign(paste("contig",i,"seq_trim",sep="_"), gsub("[0-9]{1,}","", get(paste("contig",i,"seq_pre",sep="_"))))
  }
  for(i in 1:nrow(contig_length)){
    assign(paste("contig",i,"seq_final",sep="_"), gsub(" ", "", get(paste("contig",i,"seq_trim",sep="_"))))
  }

  #Large character concatenate
  {
    for(i in c(1:nrow(contig_length))){
      assign(paste("contig",i,"nt_seq_data_frame",sep="_"),
             na.omit(as.data.frame((get(paste("contig",i,sep = "_"))[,1])[get(paste("contig",i,"seq_start",sep = "_")):get(paste("contig",i,"seq_finish",sep = "_"))])))
    }
    for(i in c(1:nrow(contig_length))){
      assign(paste("contig",i,"nt_seq_data_frame",sep="_"),
             as.data.frame(get(paste("contig",i,"nt_seq_data_frame",sep="_"))[-c(1,as.numeric(nrow(get(paste("contig",i,"nt_seq_data_frame",sep="_"))))),]))
    }
    for(i in c(1:nrow(contig_length))){
      assign(paste("contig",i,"seq_pre",sep="_"), gsub("ORIGIN ","", get(paste("contig",i,"nt_seq",sep="_"))))
      #str_extract(get(paste("contig",i,"nt_seq",sep="_")), " 1 +.{1,}"))
    }
    for(i in c(1:nrow(contig_length))){
      d=get(paste("contig",i,"nt_seq_data_frame",sep="_"))
      colnames(d)="V1"
      assign(paste("contig",i,"nt_seq_data_frame",sep="_"),d)
      rm(d)
    }
    for(i in c(1:nrow(contig_length))){
      assign(paste("contig",i,"seq_trim",sep="_"),data.frame(lapply(get(paste("contig",i,"nt_seq_data_frame",sep="_")),function(x) gsub("[0-9]{1,}","", x))))
    }
    for(i in c(1:nrow(contig_length))){
      assign(paste("contig",i,"seq_trim",sep="_"),data.frame(lapply(get(paste("contig",i,"seq_trim",sep="_")),function(x) gsub(" ","", x))))
    }
    for(i in c(1:nrow(contig_length))){
      assign(paste("contig",i,"seq","final",sep="_"),
             paste(((get(paste("contig",i,"seq","trim",sep = "_"))[,1])[1:nrow(get(paste("contig",i,"seq","trim",sep="_")))]),sep="",collapse=""))
    }
  }

  #To upper character
  for(i in 1:nrow(contig_length)){
    assign(paste("contig",i,"seq",sep="_"), toupper(get(paste("contig",i,"seq","final",sep="_"))), envir = .GlobalEnv)
  }
}


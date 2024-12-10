#' GenBank organizer
#'
#' @param df A data frame containing the data from the genbank.
#' @param gb_dir A string specifying the directory where the genbank output files are located. If NULL, the current working directory is used.
#' @return A data frame containing the baseline table.
#' @export
#'

Genbank_organizer <- function(gb_dir = NULL, save_output = FALSE) {

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

# Loop through each file in gb_files
for (genbank_file in gb_files) {

  # Step 1: Read and prepare the data
  origin <- read.csv(genbank_file, header = FALSE, sep = "\n")
  contig_row <- grep("^LOCUS", origin$V1)
  contig_length_row <- grep("source + 1..", origin$V1, ignore.case = TRUE)
  pre_contig_start <- grep("^LOCUS", origin$V1)
  pre_contig_finish <- c(pre_contig_start[-1] - 1, nrow(origin))

  #Need to check contig numbering
  k <- as.numeric(order(1:nrow(as.data.frame(contig_row))))

  # Extract and assign contig sequences
  for (i in 1:nrow(as.data.frame(contig_row))) {
    for (j in k[i]) {
      assign(paste("contig", j, sep = "_"), as.data.frame(origin[pre_contig_start[i]:pre_contig_finish[i], ]))
    }
  }

      #Calculate contig length-----------------------------------------------------------------------------------------------------------
      {
        contig_length<-as.data.frame(grep("source + 1..",origin$V1, ignore.case = TRUE, value=TRUE))
        colnames(contig_length)<-"contig_length"
        contig_length<-tidyr::separate(contig_length, col = `contig_length`, into=c("start","end"),sep="\\.\\.")
        contig_length<-as.data.frame(contig_length$end)
        contig_length<-cbind(contig_length, as.data.frame(order(as.character(1:1))))
        colnames(contig_length)<-c("length","order")
        contig_length<-contig_length[order(contig_length$order),]
        rownames(contig_length)<-NULL
        contig_length<-as.data.frame(contig_length$length)
        colnames(contig_length)<-"length"
        rm(i,j,k,contig_row,contig_length_row,origin)
        genomesize<-sum(as.numeric(as.character(contig_length$length)))
      }

      #contig definition
      for(i in 1:nrow(contig_length)){
        assign(paste("contig",i,"definition",sep="_"), get(paste("contig",i,sep="_"))[grep("DEFINITION",get(paste("contig",i,sep="_"))[,1], value="TURE"),])
      }
      for(i in 1:nrow(contig_length)){
        assign(paste("contig",i,"definition",sep="_"), gsub("^DEFINITION  |\\.$","", get(paste("contig",i,"definition",sep="_"))))
        #str_extract(get(paste("contig",i,"nt_seq",sep="_")), " 1 +.{1,}"))
      }
      definition <- contig_1_definition

      #sequence extract in gbff
      #contig_seq concatenated (CONTIG join + origin...^//)
      {
        for(i in 1:nrow(contig_length)){
          assign(paste("contig",i,"seq_start",sep="_"), grep("ORIGIN",get(paste("contig",i,sep="_"))[,1], value="FALSE"))
        }
        for(i in 1:nrow(contig_length)){
          assign(paste("contig",i,"seq_finish",sep="_"), grep("^//",get(paste("contig",i,sep="_"))[,1], value="FALSE"))
        }

        #concatenated nucleotide
        #nucleotide---####need to fix Large character########## # FIXME
        for(i in 1:nrow(contig_length)){
          assign(paste("contig",i,"nt_seq",sep="_"),
                 paste(((get(paste("contig",i,sep = "_"))[,1])[get(paste("contig",i,"seq_start",sep = "_")):get(paste("contig",i,"seq_finish",sep = "_"))]),sep="",collapse=""))
        }
        for(i in 1:nrow(contig_length)){
          assign(paste("contig",i,"seq_pre",sep="_"), gsub("ORIGIN| 1 +.{1,}","", get(paste("contig",i,"nt_seq",sep="_"))))
          #str_extract(get(paste("contig",i,"nt_seq",sep="_")), " 1 +.{1,}"))
        }
        for(i in 1:nrow(contig_length)){
          assign(paste("contig",i,"seq_trim",sep="_"), gsub("[0-9]{1,}","", get(paste("contig",i,"seq_pre",sep="_"))))
        }
        for(i in 1:nrow(contig_length)){
          assign(paste("contig",i,"seq_final",sep="_"), gsub(" ", "", get(paste("contig",i,"seq_trim",sep="_"))))
        }
        rm(list=ls(pattern="contig_[0-9]{1,}_seq_trim|pre"))
      }

      #large character concatenate-------------------------------------------------------------------------------------------------------------
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
        rm(list=ls(pattern="contig_[0-9]{1,}_seq_start|finish|trim|pre"))
        rm(list=ls(pattern="contig_[0-9]{1,}_nt_seq"))


        #to upper character
        for(i in 1:nrow(contig_length)){
          assign(paste("contig",i,"seq","final",sep="_"), toupper(get(paste("contig",i,"seq","final",sep="_"))))
        }
      }

      # TODO: sequence save # FIXME
      for(i in 1:nrow(contig_length)){
        assign(paste("contig",i,"seq","final","temp",sep="_"), data.frame("definition"=get(paste("contig",i,"definition",sep="_")),
                                                                          "top_strand"=get(paste("contig",i,"seq","final",sep="_"))
        )
        )
      }


      contig_seq_final <- contig_1_seq_final_temp
      for(i in (1:nrow(contig_length))[-c(1)]){
        if(!is.na(get(paste("contig",i,"seq","final","temp",sep="_")))[1]){
          contig_seq_final <- plyr::rbind.fill(get(paste("contig","seq","final",sep="_")),
                                         get(paste("contig",i,"seq","final","temp",sep="_"))
          )
        }else{}
      }
      rm(list=ls(pattern="contig_[0-9]{1,}_seq_final_temp"))

      contig_seq_final<- contig_seq_final %>% mutate("bottom_strand"=sapply(contig_seq_final$top_strand, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))),
                                                     "length"=nchar(top_strand))



  #contig gbff trimming
  for(i in 1:nrow(contig_length)){
    assign(paste("contig_start",i,sep = "_"),grep("gene + |[\\ ]{1}complement + \\(\\d{1,}[\\.]{2}\\d{1,})",get(paste("contig",i,sep = "_"))[,1], ignore.case = TRUE)[1])
  }
  for(i in 1:nrow(contig_length)){
    assign(paste("contig_finish",i,sep = "_"),grep("^ORIGIN",get(paste("contig",i,sep = "_"))[,1], ignore.case = FALSE))
  }


  #delete gene-vacant contig
  {
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("contig_gbff",i,sep="_"),as.data.frame(get(paste("contig",i,sep="_"))[get(paste("contig","start",i,sep="_")):(get(paste("contig","finish",i,sep="_"))-1),]))}
      else{
      }
    }

    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("gene",i,sep = "_"),grep("gene + |[\\ ]{1}complement + \\(\\d{1,}[\\.]{2}\\d{1,})",get(paste("contig","gbff",i,sep = "_"))[,1], ignore.case = TRUE))}
      else{
      }
    }
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("gene",i,sep = "_"),as.data.frame(get(paste("contig","gbff",i,sep="_"))[get(paste("gene",i,sep="_")),]))}
      else{
      }
    }
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("location",i,sep="_"),as.data.frame(stringr::str_extract((get(paste("gene",i,sep = "_")))[,1], "[0-9]{1,}[\\.]{2}>?[0-9]{1,}")))}
      else{
      }
    }
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        d=get(paste("contig","gbff",i,sep="_"))
        colnames(d)="V1"
        assign(paste("contig","gbff",i,sep="_"),d)}
      else{
      }
    }
    for (i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        d=get(paste("location",i,sep="_"))
        colnames(d)="GENE"
        assign(paste("location",i,sep="_"),d)}
      else{
      }
    }
    rm(d)

    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("CDS",i,sep="_"), tidyr::separate(get(paste("location",i,sep="_")), col="GENE", into=c("start","end"), sep="\\..>?", fill = "right"))}
      else{
      }
    }

    #find out gene location in contig_gbff-------------------------------------------------------------------------------------
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("gene","start",i,sep="_"), grep("gene + |[\\ ]{1}complement + \\(\\d{1,}[\\.]{2}\\d{1,})",get(paste("contig","gbff",i,sep="_"))[,1], ignore.case = TRUE))}
      else{
      }
    }
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("gene","finish",i,sep="_"), gene_finish<-c(get(paste("gene","start",i,sep="_"))[-1]-1,nrow(get(paste("contig","gbff",i,sep="_")))))}
      else{
      }
    }
  }
  rm(list=ls(pattern="location_[0-9]{1,}"))
  rm(list=ls(pattern="contig_finish_[0-9]{1,}"))

  #direction_gff
  for(i in (1:nrow(contig_length))){
    if(!is.na(get(paste0("contig_start_",i)))){
      assign(paste("gene","direction",i,sep="_"), stats::setNames(as.data.frame(grepl("complement",get(paste("gene",i,sep="_"))[,1], ignore.case = TRUE)),"direction"))}
    else{
    }
  }
  for(i in (1:nrow(contig_length))){
    if(!is.na(get(paste0("contig_start_",i)))){
      assign(paste("gene","direction",i,sep="_"), stats::setNames(as.data.frame(ifelse(get(paste("gene","direction",i,sep="_"))==TRUE, "-", "+")),"direction"))}
    else{
    }
  }
  for(i in (1:nrow(contig_length))){
    if(!is.na(get(paste0("contig_start_",i)))){
      assign(paste("gene","direction",i,sep="_"), data.frame(lapply(get(paste("gene","direction",i,sep="_")), as.character), stringsAsFactors=FALSE))}
    else{
    }
  }

  #pre_assign gene number-----------------------------------------------------------------------------------------------------------------
  for(i in (1:nrow(contig_length))){
    if(!is.na(get(paste0("contig_start_",i)))){
      assign(paste("number",i,sep="_"),1:nrow(get(paste("gene",i,sep="_"))))}
    else{
    }
  }
  for(j in (1:nrow(contig_length))){
    if(!is.na(get(paste0("contig_start_",j)))){
      sink(paste("genenumber",j,sep="_"));for(i in 1:nrow(get(paste("gene",j,sep="_")))){
        cat(assign(paste("genenumber",j,sep="_"),
                   paste(rep((get(paste("number",j,sep="_"))[i]),
                             length(as.integer(get(paste("gene","start",j,sep = "_"))[i]):as.integer(get(paste("gene","finish",j,sep="_"))[i])))))
            ,"\n")}
      sink()}
    else{}
  }


  #delete processing---------------------------------------------------------------------------

  rm(list=ls(pattern="genenumber|number_[0-9]{1,}"))
  rm(list=ls(pattern="gene_start_[0-9]{1,}"))
  rm(gene_finish,i,j)

  #for assign gene number read.csv-----------------------------
  for(i in (1:nrow(contig_length))){
    if(!is.na(get(paste0("contig_start_",i)))){
      assign(paste("genenumber",i,sep="_"), utils::read.csv(paste("genenumber",i,sep="_"),header=FALSE, stringsAsFactors = FALSE))
      assign(paste("a",i,sep="_"), as.data.frame(paste0("",(get(paste("genenumber",i,sep="_"))[,1])[1:nrow(get(paste("genenumber",i,sep="_")))],collapse="")))
      d=get(paste("a",i,sep="_"))
      colnames(d)="a"
      assign(paste("a",i,sep="_"),d)
      assign(paste("genenumber",i,sep="_"), tidyr::separate_rows(get(paste("a",i,sep="_")), a, sep=" "))
      assign(paste("genenumber",i,sep="_"), as.data.frame((get(paste("genenumber",i,sep="_")))[-(nrow(get(paste("contig","gbff",i,sep="_")))+1),])) #[Largest+1]
      d=get(paste("genenumber",i,sep="_"))
      colnames(d)="gene_number"
      assign(paste("genenumber",i,sep="_"),d)
    }
    else{}
  }

  rm(list=ls(pattern="a_[0-9]{1,}"))
  rm(d)

  #assign gene start to end---------------------------------------------------------------------------------
  {
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("separated",i,sep="_"), tidyr::separate(get(paste("contig_gbff",i,sep="_")), col="V1", into=c("name","value"), sep="=", fill = "right"))
        assign(paste("phase1",i,sep="_"), cbind(get(paste("genenumber",i,sep="_")),get(paste("separated",i,sep="_"))))
        assign(paste("phase2",i,sep="_"), stats::na.omit(get(paste("phase1",i,sep="_"))))
        assign(paste("phase3",i,sep="_"), unique(get(paste("phase2",i,sep="_"))))
        assign(paste("molten",i,sep="_"), reshape2::melt(get(paste("phase3",i,sep="_")), id=c("gene_number","name")))
      }
      else{}
    }

    rm(list=ls(pattern="phase1|phase2|phase3_[0-9]{1,}"))

    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("protein_id",i,sep="_"), grep("/protein_id",get(paste("molten",i,sep="_"))$name))
        assign(paste("product",i,sep="_"), grep("/product",get(paste("molten",i,sep="_"))$name))
        assign(paste("gene",i,sep="_"), grep("/gene",get(paste("molten",i,sep="_"))$name))
        assign(paste("locus","tag",i,sep="_"), grep("/locus_tag",get(paste("molten",i,sep="_"))$name))
        assign(paste("old_locus","tag",i,sep="_"), grep("/old_locus_tag",get(paste("molten",i,sep="_"))$name))
        assign(paste("translation","tag",i,sep="_"), grep("/translation",get(paste("molten",i,sep="_"))$name))
        assign(paste("note","tag",i,sep="_"), grep("/note",get(paste("molten",i,sep="_"))$name))
        assign(paste("EC","number",i,sep="_"), grep("/EC_number",get(paste("molten",i,sep="_"))$name))
        assign(paste("anticodon",i,sep="_"), grep("/anticodon",get(paste("molten",i,sep="_"))$name))
        #assign(paste("gene","direction",i,sep="_"), grep("/direction",get(paste("molten",i,sep="_"))$name))
      }
      else{}
    }

    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("remolten",i,sep="_"),get(paste("molten",i,sep="_"))[get(paste("gene",i,sep="_")),])
        assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("protein_id",i,sep="_")),]))
        assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("product",i,sep="_")),]))
        assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("locus_tag",i,sep="_")),]))
        assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("old_locus_tag",i,sep="_")),]))
        assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("translation_tag",i,sep="_")),]))
        assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("note_tag",i,sep="_")),]))
        assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("EC_number",i,sep="_")),]))
        assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("anticodon",i,sep="_")),]))
        #assign(paste("remolten",i,sep="_"),rbind(get(paste("remolten",i,sep="_")),get(paste("molten",i,sep="_"))[get(paste("gene_direction",i,sep="_")),]))
      }
      else{}
    }

    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("trimmed",i,sep="_"),reshape2::dcast(get(paste("remolten",i,sep="_")), formula = gene_number ~ name, fun.aggregate = toString))
      }
      else{}
    }

    #order gene number;low to high
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("trimmed",i,sep="_"),get(paste("trimmed",i,sep="_"))[gtools::mixedorder(as.character(get(paste("trimmed",i,sep = "_"))$gene_number)),]) #mixedorder function in gtools package
        d=get(paste("trimmed",i,sep="_"))
        rownames(d)<-NULL
        assign(paste("trimmed",i,sep="_"),d)

      }
      else{}
    }
    rm(d)
  }



  #paste CDS
  {
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("trimmed",i,sep="_"),cbind(get(paste("trimmed",i,sep="_")), get(paste("CDS",i,sep="_"))))
      }
      else{}
    }
    rm(list=ls(pattern="protein_id_|product_|locus_tag_|old_locus_tag|gene_[0-9]{1,}"))
    rm(list=ls(pattern="separated|remolten|molten|location|CDS|genenumber|contig_gbff_[0-9]{1,}"))
  }

  #sink nt_seq
  for(j in (1:nrow(contig_length))){
    if(!is.na(get(paste0("contig_start_",j)))){
      sink(paste("gene","nt",j,sep="_"));
      for(i in 1:nrow(get(paste("trimmed",j,sep="_")))){
        w=substr(get(paste("contig",j,"seq","final",sep="_")),as.integer(get(paste("trimmed",j,sep="_"))$start[i]),as.integer(get(paste("trimmed",j,sep="_"))$end[i]))
        cat(w,"\n")}
      sink()
    }
    else{}
  }

  # check---------------------
  {
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("gene","nt",i,sep="_"), utils::read.csv(paste("gene","nt",i,sep="_"),header=FALSE, stringsAsFactors = FALSE))
      }
      else{}
    }
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        d=get(paste("gene","nt",i,sep="_"))
        colnames(d)="nt_seq"
        assign(paste("gene","nt",i,sep="_"),d)
      }
      else{}
    }
    rm(d)

    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("contig_trim",i,sep="_"), cbind(get(paste("trimmed",i,sep="_")), get(paste("gene","nt",i,sep="_"))))
      }
      else{}
    }

    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("contig","trim",i,sep="_"), get(paste("contig","trim",i,sep="_")) %>% dplyr::rename_all(~ stringr::str_replace_all(., " |\\/", "")))
      }
      else{}
    }

    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("contig_final",i,sep="_"), as.data.frame(cbind(contig=get(paste("contig",i,"definition",sep="_")),
                                                                    get(paste("contig","trim",i,sep="_")),
                                                                    get(paste("gene","direction",i,sep="_")))))
      }
      else{}
    }
  }
  rm(list=ls(pattern="CDS|contig_gbff|gene_nt|genenumber|location|molten|remolten|separated|trimmed_[0-9]{1,}"))


  #Combine all contigs...Should be start contig2
  contig_final<-contig_final_1
  for(i in (1:nrow(contig_length))[-c(1)]){
    if(!is.na(get(paste0("contig_start_",i)))){
      contig_final <- plyr::rbind.fill(get(paste("contig","final",sep="_")), get(paste("contig","final",i,sep="_")))
    }
  }
  rm(list=ls(pattern="contig_trim|gene_finish|gene_direction_[0-9]{1,}"))
  rm(list=ls(pattern="contig_start|finish_[0-9]{1,}"))


  #if needed, factor values change as character
  {
    contig_final$nt_seq<-gsub(" ","",contig_final$nt_seq)
    contig_final$translation<-as.character(contig_final$translation)
    contig_final$translation<-gsub(" ","",contig_final$translation)
    contig_final$translation<-gsub("\n","",contig_final$translation)
    contig_final$amino_acid<-nchar(contig_final$translation)
  }

  contig_final$product <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", contig_final$product, perl=TRUE)
  contig_final$note <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", contig_final$note, perl=TRUE)
  rm(list=ls(pattern="contig_[0-9]{1,}"))

  undefined_gene_seq <- as.data.frame(lengths(base::regmatches(contig_final$nt_seq, base::gregexpr("n", contig_final$nt_seq, ignore.case = TRUE))))
  colnames(undefined_gene_seq)<-"N_number"

  #reverse complementation
  contig_final$nt_seq_rev_comp <- sapply(contig_final$nt_seq, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))))

  contig_final <- contig_final %>% mutate(contig_final, "rearranged_nt_seq"=
                                            case_when(grepl("\\-", `contig_final`$direction) ~ sapply(`contig_final`$nt_seq, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))),
                                                      grepl("\\+", `contig_final`$direction) ~ sapply(`contig_final`$nt_seq, function(x) as.character(Biostrings::DNAString(x))))
  )

  contig_final <- contig_final %>% mutate(contig_final, "Mw"= Peptides::mw(.$translation))
  contig_final <- contig_final %>% mutate(contig_final, "pI"= Peptides::pI(.$translation))

  rm(list=ls(pattern="note_tag_[0-9]{1,}"))
  rm(list=ls(pattern="translation_tag_[0-9]{1,}"))
  rm(list=ls(pattern="EC_number_{1,}"))
  rm(list=ls(pattern="contig_final_{1,}"))

  contig_final <- arrange.vars(contig_final, c("contig"=1,
                                               "gene_number"=2,
                                               "locus_tag"=3,
                                               "protein_id"=ncol(contig_final)-11,
                                               "product"=ncol(contig_final)-10,
                                               "start"=ncol(contig_final)-9,
                                               "end"=ncol(contig_final)-8,
                                               "amino_acid"=ncol(contig_final)-7,
                                               "Mw"=ncol(contig_final)-6,
                                               "pI"=ncol(contig_final)-5,
                                               "translation"=ncol(contig_final)-4,
                                               "rearranged_nt_seq"=ncol(contig_final)-3,
                                               "nt_seq"=ncol(contig_final)-2,
                                               "nt_seq_rev_comp"=ncol(contig_final)-1,
                                               "note"=ncol(contig_final)))
  contig_final <- contig_final %>%
    dplyr::mutate(contig_number = base::cumsum(gene_number == 1))

  contig_final <- contig_final %>%
    dplyr::select(contig, contig_number, dplyr::everything())

  contig_final <- contig_final %>%
    dplyr::mutate(
      "product" = gsub("\r\n|\n|\r", " ", product), # Replace line with a space
      "product" = stringr::str_squish(product)  # Apply str_squish to remove extra spaces
    )

  # export faa
  contig_final %>% dplyr::select(locus_tag, translation) %>% dplyr::filter(!is.na(translation) & stringr::str_trim(translation) != "") %>% dplyr::mutate("locus_tag"=paste0(">",.$locus_tag)) %>% utils::write.table(.,paste0(qdap::beg2char(genbank_file, "."),".faa"), row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)

  delfiles <- dir(path=getwd(), pattern="genenumber_[0-9]{1,}|gene_nt_[0-9]{1,}")
  file.remove(file.path(getwd(), delfiles))


  # Save the output as 'genbank_table.xlsx' if save_output is TRUE
  if (save_output) {
    output_file <- file.path(gb_dir, paste0(gsub("\\:|\\/"," ", definition),".xlsx"))
    openxlsx::write.xlsx(genbank_file, output_file, rowNames = FALSE)
    message(paste("Results have been saved to:", output_file))
  }


  # Assign the result to 'genbank_table' in the global environment
  assign("genbank_table", contig_final, envir = .GlobalEnv)
  message("The result has been saved to the R environment variable 'genbank_table'.")

  assign(paste(definition, genomesize, sep = "_"), contig_final)
}
}

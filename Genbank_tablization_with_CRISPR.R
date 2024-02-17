# Load required R packages dynamically, installing them if they are not already installed.
# This ensures the script has all the dependencies it needs to run.
pacman::p_load(readr, Peptides, tidyr, plyr, dplyr, data.table, stringr, seqinr, magrittr, reutils, reshape2, xlsx, openxlsx, gtools, Biostrings, varhandle, splitstackshape)

arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}


# Set the working directory and identify all GenBank files to be processed.
gb_dir<-getwd()
gb_file <- list.files(gb_dir, pattern = "\\.gbk$|gb$|gbff$") # list

# Iterate through each GenBank file for processing.
# Start
for(genbank in (1:length(gb_file))){
  {
    origin<-read.csv(paste0(gb_file[genbank]), header=FALSE, sep="\n")  #name
    contig_row<-grep("^LOCUS",origin$V1)
    contig_length_row<-grep("source + 1..",origin$V1, ignore.case = TRUE)
    pre_contig_start<-grep("^LOCUS",origin$V1)
    pre_contig_finish<-c(pre_contig_start[-1]-1,nrow(origin))
  }

  #Need to check contig numbering
  k<-as.numeric(order((1:nrow(as.data.frame(contig_row))))) #contig number
  for(i in (1:nrow(as.data.frame(contig_row)))){
    for(j in k[i]){
      assign(paste("contig",j,sep="_"),as.data.frame(origin[pre_contig_start[i]:pre_contig_finish[i],]))
    }}
  
  
  #Calculate contig length-----------------------------------------------------------------------------------------------------------
  {
    contig_length<-as.data.frame(grep("source + 1..",origin$V1, ignore.case = TRUE, value=TRUE))
    colnames(contig_length)<-"contig_length"
    contig_length<-separate(contig_length, col = `contig_length`, into=c("start","end"),sep="\\.\\.")
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
      contig_seq_final <- rbind.fill(get(paste("contig","seq","final",sep="_")), 
                                     get(paste("contig",i,"seq","final","temp",sep="_"))
                                     )
    }else{}
  }
  rm(list=ls(pattern="contig_[0-9]{1,}_seq_final_temp"))
  
  contig_seq_final<- contig_seq_final %>% mutate("bottom_strand"=sapply(contig_seq_final$top_strand, function(x) as.character(Biostrings::reverseComplement(DNAString(x)))),
                                                 "length"=nchar(top_strand))
  
  
  # Feature extraction
  {
    for(genbank in (1:length(gb_file))){
      origin<- readr::read_fwf(paste0(gb_file[genbank]), readr::fwf_widths(c(12,Inf)))
      
      # delete nucleotide sequences
      origin <- origin %>% dplyr::slice(., 1:as.numeric(grep("ORIGIN",origin$X1,ignore.case = FALSE)[1])-1)
      origin <- origin %>% sapply(function(x) str_squish(x)) %>% as.data.frame()
      #origin %>% group_by(X1) %>% summarise()
      
      # report table
      report <- origin %>% dplyr::slice(., 1:as.numeric(grep("FEATURES",origin$X1,ignore.case = FALSE)[1])-1)
      report <- report %>% sapply(function(x) str_squish(x)) %>% as.data.frame()
      
      # feature table 
      feature <- origin %>% dplyr::slice(., as.numeric(grep("FEATURES",origin$X1,ignore.case = FALSE)[1]):nrow(.))
      feature <- feature %>% dplyr::slice(., as.numeric(grep("gene",feature$X1,ignore.case = FALSE)[1]):nrow(.))
      feature <- feature %>% fill(X1)
      #feature <- feature %>% dplyr::filter(!X1 %in% c("gene","CDS","rRNA","tRNA","tmRNA","ncRNA"))
      
      feature_info <- rle(feature$X1)
      feature_info <- data.frame(unclass(feature_info))
      feature_info <- feature_info %>% mutate("end"=cumsum(feature_info$lengths))
      feature_info <- feature_info %>% mutate("start"=end-lengths+1)
      feature_info <- arrange.vars(feature_info, c("values"=1,
                                                   "start"=2,
                                                   "end"=3,
                                                   "lengths"=4)
      )
      
      CRISPR <- feature_info %>% dplyr::filter(str_detect(values,"repeat"))
    if(nrow(CRISPR)==!0){
      CRISPR <- CRISPR %>% mutate("discriminator"=1:nrow(CRISPR))
      temp_rows <- as.numeric()
      
      for(i in(1:nrow(CRISPR))){
        temp_rows <- c(temp_rows, CRISPR$start[i]:CRISPR$end[i])
      }
      CRISPR_table <- feature[temp_rows,] 
      rm(temp_rows)
      CRISPR_table <- CRISPR_table %>% mutate("discriminator"=rep(CRISPR$discriminator, CRISPR$lengths))
      CRISPR_table <- CRISPR_table %>% separate(., col="X2", into=c("category","value"), sep="\\=")
      CRISPR_table <- CRISPR_table %>% separate(., col="category", into=c("category","value2"), sep=" ", fill="right", extra="merge")
      CRISPR_table <- CRISPR_table %>% mutate("category"=gsub("/", "", category)) %>% dplyr::select(-X1)
      CRISPR_table <- CRISPR_table %>% replace(is.na(CRISPR_table), "")
      
      CRISPR_table <- reshape::melt(CRISPR_table, id=c("discriminator","category"),na.rm=TRUE) %>% reshape2::dcast(formula = discriminator ~ category, fun.aggregate = toString)
      
      CRISPR_table <- CRISPR_table %>% mutate_all(funs(str_replace_all(., " ", "")))
      CRISPR_table <- CRISPR_table %>% mutate_all(funs(str_replace_all(., "^,{1,}|,{1,}$", "")))
      CRISPR_table <- CRISPR_table %>% mutate_all(funs(str_replace_all(., "\"", "")))
      CRISPR_table <- CRISPR_table %>% tidyr::separate(region, into=c("repeat_start","repeat_end"), sep="\\.\\.")
    
    for(i in 1:nrow(contig_length)){
      assign(paste("CRISPR_table", i, sep="_"), CRISPR_table %>% mutate("sequence" = str_sub(get(paste("contig",i,"seq_final",sep="_")),repeat_start,repeat_end)))
      assign(paste("CRISPR_table", i, sep="_"), get(paste("CRISPR_table", i, sep="_")) %>% mutate("sequence_length" = nchar(sequence)))
      }
    }
      else{
      }
    }
  }
  
  #Combine all contigs...Should be start contig2
  CRISPR_table<-CRISPR_table_1
  for(i in (1:nrow(contig_length))[-c(1)]){
    #if(!is.na(get(paste0("CRISPR_table_",i)))){
    if(exists(paste0("CRISPR_table_",i))){
      CRISPR_table<-rbind.fill(get(paste("CRISPR_table",sep="_")),get(paste("CRISPR_table",i,sep="_")))
    }
  }
  rm(list=ls(pattern="CRISPR_table_[0-9]{1,}"))
  
  
  library(splitstackshape)
  
  CRISPR_table <- CRISPR_table %>% mutate(rpt_unit_seq= toupper(rpt_unit_seq))
  CRISPR_table <- CRISPR_table %>% mutate_at("sequence", ~na_if(., '')) %>% dplyr::filter(!is.na(sequence))
  #complete.cases(CRISPR_table)
  ####################################################################################################################################
  
  
  CRISPR_spacer <- list()
  for(i in 1:nrow(CRISPR_table)){
  CRISPR_spacer <- append(CRISPR_spacer, strsplit(CRISPR_table$sequence[i], split = CRISPR_table$rpt_unit_seq[i]) %>% setNames(CRISPR_table$rpt_unit_seq[i]))
  }
  
  # FIXME
  CRISPR_by_spacer <- data.frame()
  for (i in 1:length(CRISPR_spacer)){
    "discriminator"=names(CRISPR_spacer[i])
    by_spacer <- cbind(discriminator, CRISPR_spacer[i] %>% as.data.frame() %>% setNames("spacer"))
    CRISPR_by_spacer <- rbind(CRISPR_by_spacer, by_spacer)
    rm(discriminator)
  }
  
  CRISPR_by_spacer <- CRISPR_by_spacer %>% mutate_all(na_if,"")
  CRISPR_by_spacer <- CRISPR_by_spacer %>% drop_na(spacer)
  CRISPR_by_spacer <- CRISPR_by_spacer %>% mutate("length"=nchar(spacer))
  CRISPR_by_spacer_dist <- CRISPR_by_spacer %>% distinct(spacer, .keep_all = TRUE)
  CRISPR_by_spacer_dist_N20 <- CRISPR_by_spacer_dist %>% mutate("N20"=str_sub(spacer, start= -20),
                                                                "length_N20"=nchar(N20))

write.xlsx(x = CRISPR_by_spacer_dist_N20, "CRISPR_by_spacer.xlsx", row.names = FALSE)
# fasta with original spacer
CRISPR_by_spacer_dist_N20 %>% dplyr::select(discriminator, spacer) %>% mutate("discriminator"=paste0(">",.$discriminator)) %>% write.table(.,"CIRSPR_spacer.fasta", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)
# # fasta with N20 from end of spacer
CRISPR_by_spacer_dist_N20 %>% dplyr::select(discriminator, N20) %>% mutate("discriminator"=paste0(">",.$discriminator)) %>% write.table(.,"CIRSPR_spacer_N20.fasta", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)


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
        assign(paste("location",i,sep="_"),as.data.frame(str_extract((get(paste("gene",i,sep = "_")))[,1], "[0-9]{1,}[\\.]{2}>?[0-9]{1,}")))}
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
        assign(paste("CDS",i,sep="_"), separate(get(paste("location",i,sep="_")), col="GENE", into=c("start","end"), sep="\\..>?", fill = "right"))}
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
      assign(paste("gene","direction",i,sep="_"), setNames(as.data.frame(grepl("complement",get(paste("gene",i,sep="_"))[,1], ignore.case = TRUE)),"direction"))}
    else{
    }
  }
  for(i in (1:nrow(contig_length))){
    if(!is.na(get(paste0("contig_start_",i)))){
      assign(paste("gene","direction",i,sep="_"), setNames(as.data.frame(ifelse(get(paste("gene","direction",i,sep="_"))==TRUE, "-", "+")),"direction"))}
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
      assign(paste("genenumber",i,sep="_"),read.csv(paste("genenumber",i,sep="_"),header=FALSE, stringsAsFactors = FALSE))
      assign(paste("a",i,sep="_"),as.data.frame(paste0("",(get(paste("genenumber",i,sep="_"))[,1])[1:nrow(get(paste("genenumber",i,sep="_")))],collapse="")))
      d=get(paste("a",i,sep="_"))
      colnames(d)="a"
      assign(paste("a",i,sep="_"),d)
      assign(paste("genenumber",i,sep="_"),separate_rows(get(paste("a",i,sep="_")), a, sep=" "))
      assign(paste("genenumber",i,sep="_"),as.data.frame((get(paste("genenumber",i,sep="_")))[-(nrow(get(paste("contig","gbff",i,sep="_")))+1),])) #[Largest+1]
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
        assign(paste("separated",i,sep="_"),separate(get(paste("contig_gbff",i,sep="_")), col="V1", into=c("name","value"), sep="=", fill = "right"))
        assign(paste("phase1",i,sep="_"),cbind(get(paste("genenumber",i,sep="_")),get(paste("separated",i,sep="_"))))
        assign(paste("phase2",i,sep="_"),na.omit(get(paste("phase1",i,sep="_"))))
        assign(paste("phase3",i,sep="_"),unique(get(paste("phase2",i,sep="_")))) 
        assign(paste("molten",i,sep="_"),reshape2::melt(get(paste("phase3",i,sep="_")), id=c("gene_number","name")))
      }
      else{}
    }
    
    rm(list=ls(pattern="phase1|phase2|phase3_[0-9]{1,}"))
    
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("protein_id",i,sep="_"),grep("/protein_id",get(paste("molten",i,sep="_"))$name))
        assign(paste("product",i,sep="_"),grep("/product",get(paste("molten",i,sep="_"))$name))
        assign(paste("gene",i,sep="_"),grep("/gene",get(paste("molten",i,sep="_"))$name))
        assign(paste("locus","tag",i,sep="_"),grep("/locus_tag",get(paste("molten",i,sep="_"))$name))
        assign(paste("old_locus","tag",i,sep="_"),grep("/old_locus_tag",get(paste("molten",i,sep="_"))$name))
        assign(paste("translation","tag",i,sep="_"),grep("/translation",get(paste("molten",i,sep="_"))$name))
        assign(paste("note","tag",i,sep="_"),grep("/note",get(paste("molten",i,sep="_"))$name))
        assign(paste("EC","number",i,sep="_"),grep("/EC_number",get(paste("molten",i,sep="_"))$name))
        #assign(paste("gene","direction",i,sep="_"),grep("/direction",get(paste("molten",i,sep="_"))$name))
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
        assign(paste("trimmed",i,sep="_"),get(paste("trimmed",i,sep="_"))[mixedorder(as.character(get(paste("trimmed",i,sep = "_"))$gene_number)),]) #mixedorder function in gtools package
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
        assign(paste("gene","nt",i,sep="_"),read.csv(paste("gene","nt",i,sep="_"),header=FALSE, stringsAsFactors = FALSE))
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
        assign(paste("contig_trim",i,sep="_"),cbind(get(paste("trimmed",i,sep="_")),get(paste("gene","nt",i,sep="_"))))
      }
      else{}
    }
    
    for(i in (1:nrow(contig_length))){
      if(!is.na(get(paste0("contig_start_",i)))){
        assign(paste("contig","trim",i,sep="_"), get(paste("contig","trim",i,sep="_")) %>% rename_all(funs(stringr::str_replace_all(., " |\\/", ""))))
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
      contig_final<-rbind.fill(get(paste("contig","final",sep="_")),get(paste("contig","final",i,sep="_")))
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
  
  undefined_gene_seq <- as.data.frame(lengths(regmatches(contig_final$nt_seq, gregexpr("n", contig_final$nt_seq, ignore.case = TRUE))))
  colnames(undefined_gene_seq)<-"N_number"
  
  #reverse complementation
  contig_final$nt_seq_rev_comp <- sapply(contig_final$nt_seq, function(x) as.character(Biostrings::reverseComplement(DNAString(x))))
  
  contig_final <- contig_final %>% mutate(contig_final, "rearranged_nt_seq"= 
                                            case_when(grepl("\\-", `contig_final`$direction) ~ sapply(`contig_final`$nt_seq, function(x) as.character(Biostrings::reverseComplement(DNAString(x)))),
                                                      grepl("\\+", `contig_final`$direction) ~ sapply(`contig_final`$nt_seq, function(x) as.character(DNAString(x))))
  )
  
  contig_final <- contig_final %>% mutate(contig_final, "Mw"= Peptides::mw(.$translation))
  contig_final <- contig_final %>% mutate(contig_final, "pI"= Peptides::pI(.$translation))
  
  #have to check
  #'3to5'<-contig_final[(grep("\\-", `contig_final`$direction, ignore.case = TRUE)),]
  #`3to5`$`nt_seq`<-sapply(`3to5`$`nt_seq`,function(x) as.character(reverseComplement(DNAString(x))))
  #'5to3'<-contig_final[(grep("\\+", `contig_final`$direction, ignore.case = TRUE)),]
  #rearrange<-rbind(`3to5`,`5to3`)
  
  #nt_seq<-`rearrange`[(mixedorder(`rearrange`$`locus_tag`)),]
  #nt_seq<-setNames(as.data.frame(nt_seq$nt_seq),"rearrange_nt_seq")
  #contig_final<-cbind(contig_final,nt_seq_rev_comp)
  #rm(list=c("3to5","5to3","rearrange","nt_seq","i","j","w"))
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
  
  xlsx::write.xlsx(contig_final,paste0(gsub("\\:|\\/"," ", definition),".xlsx"),row.names = FALSE)
  assign(paste(definition, genomesize, sep = "_"), contig_final)
  
  #faa with simple locus_tag
  #contig_final %>% select(locus_tag, translation) %>% na_if("") %>% drop_na(translation) %>%mutate("locus_tag"=paste0(">",.$locus_tag)) %>% write.table(.,"CDS_fasta.faa", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)
  
  #non-protein deletion/ protein extraction for COG annotation
  annot_protein<-contig_final
  #annot_protein<-annot_protein %>% mutate_all(na_if,"")
  
  empty_as_na <- function(x){
    if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
    ifelse(as.character(x)!="", x, NA)
  }
  annot_protein <- annot_protein %>% mutate_each(funs(empty_as_na)) 
  annot_protein_fasta<-data.frame(lapply(annot_protein, as.character), stringsAsFactors=FALSE)
  #delete NA containing colomn's row
  annot_protein_fasta<-annot_protein_fasta[complete.cases(annot_protein_fasta$translation),]
  annot_protein_fasta<-as.data.frame(cbind(annot_protein_fasta$locus_tag,annot_protein_fasta$translation))
  
  annot_protein_fasta<-data.frame(lapply(annot_protein_fasta, as.character), stringsAsFactors=FALSE)
  annot_protein_fasta<-as.data.frame(cbind(paste0(">",annot_protein_fasta$V1),annot_protein_fasta$V2))
  
  #for COGsoft fasta format gene, genome multi fasta format
  write.table(annot_protein_fasta,"CDS_fasta.faa", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)
  rm(annot_protein_fasta, annot_protein)
  
  
  delfiles <- dir(path=getwd(), pattern="genenumber_[0-9]{1,}|gene_nt_[0-9]{1,}")
  file.remove(file.path(getwd(), delfiles))
  
  rm.all.but(c(paste(definition, genomesize, sep = "_"), 
               'gb_file', 
               'genomesize', 
               #'contig_length', 
               'contig_final',
               'CRISPR_table'))
}




#mylist <- mget(ls())
#contig_final <- rbind.fill(mylist)

pacman::p_load(RFLPtools,ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, rhmmer, stringr, seqinr, ggtree, magrittr, sjmisc, ggseqlogo, genoPlotR, gridExtra, seqinr, ape, rentrez, reutils, stringr, reshape2, xlsx, openxlsx, gtools, Biostrings,varhandle)


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
  
  
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
} 

annotation_list<- names(Filter(isTRUE, eapply(.GlobalEnv, is.data.frame)))


'MEROPS-MPRO' <- read.csv(file="merops_scan_lib.txt", header=FALSE, sep = "\t", quote = "", fill=TRUE)

name <- setNames(as.data.frame(grep("^>", `MEROPS-MPRO`$V1, value = TRUE)),"name")
name <- as.data.frame(name)
name$name <- as.character(name$name)

#IF not working separate with "Error in nchar(u, itype) : invalid multibyte string, element 1"
#Sys.setlocale(category ="LC_ALL", locale="US")

name<- tidyr::separate(name, col="name", sep = " \\- ", into=c("MEROPS_ACCESSION","V1"), extra="merge", fill = "left")
name<- tidyr::separate(name, col="V1", sep = "\\[", into=c("Protease","V1"), extra="merge", fill = "left")
name<- tidyr::separate(name, col="V1", sep = "\\#", into=c("MEROPS_ID","V1"), extra="merge", fill = "left")
name<- tidyr::separate(name, col="V1", sep = "\\#", into=c("SUBFAMILY","V1"), extra="merge", fill = "left")
name<- tidyr::separate(name, col="V1", sep = "\\~", into=c("unit","source"), extra="merge", fill = "left")
name$MEROPS_ID <- gsub("\\]","",name$MEROPS_ID)
name$unit <- gsub("\\{|}","",name$unit)
name$source <- gsub("source ","",name$source)
name<- tidyr::separate(name, col="Protease", sep = "(?>\\()", into=c("protease","organism"), extra="merge", fill = "left")
name<- tidyr::separate(name, col="organism", sep = "(?>\\()", into=c("organism","suborganism"), extra="merge", fill = "left")
name$organism <- gsub("\\)","",name$organism)
name$suborganism <- gsub("\\)","",name$suborganism)
name$MEROPS_ACCESSION <- gsub(">","",name$MEROPS_ACCESSION)

name<- name %>% mutate_all(na_if,"")
name[is.na(name)]  <- ""
name<- as.data.frame(apply(name,2,function(x)gsub('\\s+$', '',x)))




seq_full<-grep("^>",`MEROPS-MPRO`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,as.numeric(nrow(`MEROPS-MPRO`)))
sink("MEROPS-MPRO_seq.csv");for(i in 1:as.numeric(nrow(name))){
  Sequence=paste(`MEROPS-MPRO`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(Sequence,"\n")
}

`MEROPS-MPRO_sequence`<-read.csv("MEROPS-MPRO_seq.csv",header=FALSE,sep="\n")
`MEROPS-MPRO_sequence`$V1 <- gsub("\\s+", "",`MEROPS-MPRO_sequence`$V1)
`MEROPS-MPRO_sequence`$V1 <- gsub("\\W", "",`MEROPS-MPRO_sequence`$V1)
names(`MEROPS-MPRO_sequence`)[1] <- "AA_sequence"

`MEROPS-MPRO`<-cbind(name, `MEROPS-MPRO_sequence`)

rm(`MEROPS-MPRO_sequence`)
rm(name)


write.xlsx(`MEROPS-MPRO`,"MEROPS-MPRO_DB.xlsx")  

`MEROPS-MPRO_fasta` <- `MEROPS-MPRO` %>% select(MEROPS_ACCESSION, AA_sequence)
`MEROPS-MPRO_fasta`$MEROPS_ACCESSION<- paste0(">",`MEROPS-MPRO_fasta`$MEROPS_ACCESSION)
write.table(x = `MEROPS-MPRO_fasta`, file = "MEROPS-MPRO.fasta", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\n")






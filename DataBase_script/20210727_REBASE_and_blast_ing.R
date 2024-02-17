pacman::p_load(RFLPtools, ggplot2, ggnewscale, GGally, ggtree, ggseqlogo, genoPlotR, gridExtra, 
               tidyr, plyr, dplyr, rhmmer, stringr, seqinr, ggtree, magrittr, sjmisc, 
               rentrez, reutils, reshape2, xlsx, openxlsx, gtools, Biostrings, varhandle)


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



#For eggnog-mapper or blastp or hmmer
#non-protein deletion/ protein extraction for COG annotation
annotation_list<- names(Filter(isTRUE, eapply(.GlobalEnv, is.data.frame)))


REBASE_DB <- read.csv(file="REBASE_DB_protein_gold_seqs.txt", header=FALSE, sep = "\t", quote = "", fill=TRUE)
REBASE_DB <- read.csv(file="REBASE_DB_protein_gold_seqs.txt", header=TRUE, quote = "", fill=TRUE)

REBASE_DB <- setNames(as.data.frame(REBASE_DB[!grepl("^<>$",REBASE_DB$REBASE.protein..gold..sequences),]),"REBASE")

name <- setNames(as.data.frame(grep("^>", `REBASE_DB`$REBASE, value = TRUE)),"name")
name <- as.data.frame(name)

test<- separate(name,sep = "\t",col="name", into=paste0("V",1:12))
test<- separate(test,sep = ">REBASE\\:",col="V1", into=c("V1","REBASE"), remove = TRUE) %>% select(!"V1")
test<- separate(test,sep = "EnzType\\:",col="V2", c("V2","EnzType")) %>% select(!"V2")
test<- separate(test,sep = "RecSeq\\:",col="V3", c("Org#","RecSeq"))

test<- separate(test,sep = "OrgName\\: ",col="V4", c("Org#2","OrgName"))
test<- unite(test, col="Org#",`Org#`, `Org#2`, na.rm = TRUE, sep="")
test<- separate(test,sep = "Org\\#\\: ",col="Org#", c("V4","Org#")) %>% select(!"V4")


test<- separate(test,sep = "Source\\: ",col="V5", c("OrgName2","Source"))
test<- separate(test,sep = "OrgName\\: ",col="OrgName2", c("V5","OrgName2")) %>% select(!"V5")
test<- unite(test, col="OrgName",`OrgName`, `OrgName2`, na.rm = TRUE, sep = "")

test<- separate(test,sep = "GenBank\\:",col="V6", c("Source2","GenBank"))
test<- separate(test,sep = "Source\\: ",col="Source2", c("V6","Source2")) %>% select(!"V6")
test<- unite(test, col="Source",`Source`, `Source2`, na.rm = TRUE, sep = "")


test<- separate(test,sep = "SeqLength\\:",col="V7", c("GenBank2","SeqLength"))
test<- separate(test,sep = "GenBank\\:",col="GenBank2", c("V7","GenBank2")) %>% select(!"V7")
test<- unite(test, col="GenBank",`GenBank`, `GenBank2`, na.rm = TRUE, sep = "")


#Locus Protein Id UniProt GI
test<- separate(test,sep = "Locus\\:",col="V8", c("SeqLength2","Locus"))
test<- separate(test,sep = "ProteinId\\:",col="SeqLength2", c("SeqLength2","ProteinId"))
test<- separate(test,sep = "SeqLength\\:",col="SeqLength2", c("V8","SeqLength2")) %>% select(!"V8")
test<- unite(test, col="SeqLength",`SeqLength`, `SeqLength2`, na.rm = TRUE, sep = "")


test<- separate(test,sep = "Locus\\:",col="V9", c("ProteinId2","Locus"))
test<- separate(test,sep = "ProteinId\\:",col="ProteinId2", c("UniProt","ProteinId2"))

test<- separate(test,sep = "ProteinId\\:",col="V10", c("UniProt2","ProteinId3"))
test<- unite(test, col="ProteinId",`ProteinId`, `ProteinId2`,`ProteinId3`, na.rm = TRUE, sep = "")


test<- separate(test,sep = "GI\\:",col="UniProt2", c("UniProt2","GI"))
test<- separate(test,sep = "GI\\:",col="V11", c("UniProt3","GI2"))
test<- unite(test, col="UniProt",`UniProt`, `UniProt2`,`UniProt3`, na.rm = TRUE, sep = "")
test<- separate(test,sep = "UniProt\\:",col="UniProt", c("V9","UniProt")) %>% select(!"V9")


test<- separate(test,sep = "GI\\:",col="V12", c("V10","GI3")) %>% select(!"V10")
test<- unite(test, col="GI",`GI`, `GI2`,`GI3`, na.rm = TRUE, sep = "")

test<- test %>% mutate_all(na_if,"")
test[is.na(test)]  <- ""

REBASE_DB_protein_gold_seqs <- test

#z<- test %>%dplyr::select(-contains("EnzType"))
"\t$"

seq_full<-grep("^>",`REBASE_DB`$REBASE) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,as.numeric(nrow(`REBASE_DB`)))
sink("REBASE_seq.csv");for(i in 1:as.numeric(nrow(name))){
  Sequence=paste(`REBASE_DB`$REBASE[seq_start[i]:seq_finish[i]],collapse = "")
  cat(Sequence,"\n")
}

`REBASE_sequence`<-read.csv("REBASE_seq.csv",header=FALSE,sep="\n")
REBASE_sequence$V1 <- gsub("\\s+", "",REBASE_sequence$V1)
REBASE_sequence$V1 <- gsub("\\W", "",REBASE_sequence$V1)
names(REBASE_sequence)[1] <- "AA_sequence"

`REBASE_DB_protein_gold_seqs`<-cbind(REBASE_DB_protein_gold_seqs, REBASE_sequence)

rm(REBASE_sequence)
rm(name)
rm(REBASE_DB)

write.xlsx(REBASE_DB_protein_gold_seqs,"REBASE_DB.xlsx")  

REBASE_DB_fasta <- REBASE_DB_protein_gold_seqs %>% select(REBASE, AA_sequence)
REBASE_DB_fasta$REBASE<- paste0(">",REBASE_DB_fasta$REBASE)
write.table(x = REBASE_DB_fasta, file = "All_REBASE_Gold_Standards_Protein.fasta", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\n")

######################################################################################################################################################
#Extract .. test...
extract(col = "document",
        into = c("is_score", "is_movie"),
        regex = "(??Á¡).*?(??È­)",
        remove=TRUE) 

z <- name %>% tidyr::extract(col=name, into="REBASE", regex="^>REBASE:+")
z<- as.data.frame(str_extract(name$name, "\\bRecSeq:.{1,}\t{1}$"))


REBASE <- setNames(as.data.frame(substr(name$name, regexpr("^>REBASE:",name$name), regexpr("\t",name$name))),"REBASE")
z <- setNames(as.data.frame(str_rtr_replace_alle$name, substr(name$name, regexpr("^>REBASE:",name$name), regexpr("\t",\ame$name)),"")),"name")

Recstr_replace(name$name[53],">REBASE:Aco12261II \\(RM.Aco12261II\\)\t","")
Seq <- setNames(as.data.frame(substr(stre$name, regexpr("^RecSeq:",name$name), regexpr("\t",name$name))),"REBASE")

as.vector(regexpr("\t",name$name)[[1]])
as.vector(attr(regexpr("ing",x),"match.length"))

name<- separate(name, col=name, into=c("V1","REBASE"),sep=">REBASE:",extra="merge", fill="right")
name<- separate(name, col=REBASE, into=c("REBASE","EnzType"),sep="EnzType:",extra="merge", fill="right")
name<- separate(name, col=EnzType, into=c("EnzType","RecSeq"),sep="RecSeq:",extra="merge", fill="right")

######################################################################################################################################################









#Blast with my data
REBASE_results <- read.csv(file="Results", header=FALSE)
REBASE_results <- REBASE_results %>% filter(!grepl("^#",REBASE_results$V1))
REBASE_results <- separate(REBASE_results,col="V1",into=c("subject","acc.ver","% identity","alignment length","mismatches","gqp opens","q. start","q. end","s. start","s. end","evalue","bit score"),sep="\t")
REBASE_results <- REBASE_results %>% filter(!is.na(acc.ver))







#pfam hmmer
blast_dir<-getwd()
blast_file <- list.files(blast_dir, pattern = "Results") # list


REBASE_result<- read.csv(file="Results",header=FALSE)
for(REBASE_list in 1:length(blast_file)){
  assign(paste(gsub("_REBASE\\.txt","",blast_file[REBASE_list]),"pfam", sep="_"), read.blast(paste0(blast_file[REBASE_list])))  #name
}
read.blast(file = "Results",sep = "\t")

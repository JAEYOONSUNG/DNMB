{
  `NDARO`<-read.csv("AMRProt", header=FALSE, sep="\n")
  name <- grep("^>", `NDARO`$V1, value = TRUE)
  name <- as.data.frame(name)
  seq_full<-grep("^>",`NDARO`$V1) 
  seq_start<-seq_full+1
  seq_finish<-c(seq_full[-1]-1,nrow(`NDARO`))
  sink("NDARO_seq.csv");for(i in 1:nrow(`name`)){
    w=paste(`NDARO`$V1[seq_start[i]:seq_finish[i]],collapse = "")
    cat(w,"\n")
  }
  `NDARO_sequence`<-read.csv("NDARO_seq.csv",header=FALSE,sep="\n")
  `NDARO_final`<-cbind(name, `NDARO_sequence`)
  rm(i,seq_finish,seq_full,seq_start,w,`NDARO_sequence`,name,`NDARO`)
  sink()
}

write.table(NDARO_final,"NDARO_protein_fasta.txt",sep="\n",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")



colnames(NDARO_final)<-c("name","name2")
NDARO_sep <- separate(NDARO_final, col=name, into=c("V1","V2","V3","V4","V5","V6"), sep=">", fill = "right")
NDARO_sepp<-unite(NDARO_sep, V1, V2, V3, remove=TRUE, sep="",na.rm=TRUE)
NDARO_sepp<-NDARO_sepp[,c(2,5)]
NDARO_seppp <- separate(NDARO_sepp, col=name2, into=c("V2","V3","V4","V5","V6","V7"), sep=">", fill = "right")
NDARO_sepppp<-unite(NDARO_seppp, A1, V1, V2, remove=TRUE, sep=" ",na.rm=TRUE)
write.table(NDARO_sepppp,"NDARO_protein_fasta_t.txt",sep="\n",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")

NDARO_seppppp<-read.csv("NDARO_protein_fasta_t.txt", header=FALSE, sep="\n")
NDARO<-setNames(as.data.frame(gsub("^gi","\\>gi",NDARO_seppppp$V1)),"V1")
NDARO<-NDARO %>% mutate_if(is.factor, as.character)

`NDARO_right`<-separate(NDARO, col=V1, into=c("V2","V3","V4","V5"), sep="\\|", fill = "right")
NDARO_gi<-setNames(as.data.frame(NDARO_right$V5),"V1")
NDARO_gi$V1<-as.character(NDARO_gi$V1)
NDARO_gi2<-na.omit(NDARO_gi)
uniq_NDARO_gi2<-unique(NDARO_gi2)

#download...
for(i in 1:nrow(NDARO_gi2)){
  ids<-esearch(NDARO_gi2$V1[i], db = "protein", rettype = "uilist", retmode = "xml",retstart = 0, retmax = 99999999, usehistory = TRUE)
  efetch(ids, db = "protein", rettype = "fasta", retmax = 99999999, retmode="text",outfile=paste0(i,"_protein",".txt"))
}

for(i in 1){
  assign(paste("fasta",i,sep="_"),read.table(paste0(i,"_protein",".txt"),header=FALSE, stringsAsFactors = FALSE, col.names = "V1", sep = "\n"))
}
fasta<-fasta_1

for(i in (1:(nrow(NDARO_gi2)))[-1]){
  assign(paste("fasta",i,sep="_"),read.table(paste0(i,"_protein",".txt"),header=FALSE, stringsAsFactors = FALSE, col.names = "V1", sep = "\n"))
  fasta<-rbind(fasta,get(paste("fasta",i,sep="_")))
  rm(list=ls(pattern="fasta_[0-9]{1,}"))
}
###
name <- grep("^>", `fasta`$V1, value = TRUE)
name <- as.data.frame(name)
seq_full<-grep("^>",`fasta`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,nrow(`fasta`))
sink("fasta.txt");for(i in 1:nrow(`name`)){
  w=paste(`fasta`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(w,"\n")
}
`fasta_trim`<-read.table("fasta.txt",header=FALSE,sep="\n")
`fasta_final`<-cbind(name, fasta_trim)
fasta_final2<-unique(fasta_final)

write.table(fasta_final2,"NDARO_db.txt",sep="\n",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")


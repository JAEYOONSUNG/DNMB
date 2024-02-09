{
  `victor`<-read.csv("victor_gen_downloads_protein.php.txt", header=FALSE, sep="\n")
  name <- grep("^>", `victor`$V1, value = TRUE)
  name <- as.data.frame(name)
  seq_full<-grep("^>",`victor`$V1) 
  seq_start<-seq_full+1
  seq_finish<-c(seq_full[-1]-1,nrow(`victor`))
  sink("victor_seq.csv");for(i in 1:nrow(`name`)){
    w=paste(`victor`$V1[seq_start[i]:seq_finish[i]],collapse = "")
    cat(w,"\n")
  }
  `victor_sequence`<-read.csv("victor_seq.csv",header=FALSE,sep="\n")
  `victor_final`<-cbind(name, `victor_sequence`)
  rm(i,seq_finish,seq_full,seq_start,w,`victor_sequence`,name,`victor`)
  sink()
}

write.table(victor_sepp,"victor_protein_fasta.txt",sep="\n",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")
colnames(victor_final)<-c("name","name2")
victor_sep <- separate(victor_final, col=name, into=c("V1","V2","V3","V4","V5","V6"), sep=">", fill = "right")
victor_sepp<-unite(victor_sep, V1, V2, V3, remove=TRUE, sep="",na.rm=TRUE)
victor_sepp<-victor_sepp[,c(2,5)]
victor_seppp <- separate(victor_sepp, col=name2, into=c("V2","V3","V4","V5","V6","V7"), sep=">", fill = "right")
victor_sepppp<-unite(victor_seppp, A1, V1, V2, remove=TRUE, sep=" ",na.rm=TRUE)
write.table(victor_sepppp,"victor_protein_fasta_t.txt",sep="\n",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")

victor_seppppp<-read.csv("victor_protein_fasta_t.txt", header=FALSE, sep="\n")
Victor<-setNames(as.data.frame(gsub("^gi","\\>gi",victor_seppppp$V1)),"V1")
Victor<-Victor %>% mutate_if(is.factor, as.character)

`Victor_right`<-separate(Victor, col=V1, into=c("V2","V3","V4","V5"), sep="\\|", fill = "right")
Victor_gi<-setNames(as.data.frame(Victor_right$V5),"V1")
Victor_gi$V1<-as.character(Victor_gi$V1)
Victor_gi2<-na.omit(Victor_gi)
uniq_Victor_gi2<-unique(Victor_gi2)

#download...
for(i in 1:nrow(Victor_gi2)){
  ids<-esearch(Victor_gi2$V1[i], db = "protein", rettype = "uilist", retmode = "xml",retstart = 0, retmax = 99999999, usehistory = TRUE)
  efetch(ids, db = "protein", rettype = "fasta", retmax = 99999999, retmode="text",outfile=paste0(i,"_protein",".txt"))
}

for(i in 1){
  assign(paste("fasta",i,sep="_"),read.table(paste0(i,"_protein",".txt"),header=FALSE, stringsAsFactors = FALSE, col.names = "V1", sep = "\n"))
}
fasta<-fasta_1

for(i in (1:(nrow(Victor_gi2)))[-1]){
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

write.table(fasta_final2,"victor_db.txt",sep="\n",row.names=FALSE, col.names=FALSE, quote=FALSE, na="")


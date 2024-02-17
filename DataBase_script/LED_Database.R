# Load necessary libraries for data manipulation, plotting, and bioinformatics analysis.
# pacman::p_load ensures that missing packages are installed and then loaded.
pacman::p_load(RFLPtools, ggplot2, ggnewscale, GGally, ggtree, ggseqlogo, genoPlotR, gridExtra,
               tidyr, plyr, dplyr, rhmmer, stringr, seqinr, 
               magrittr, reutils, reshape2, xlsx, openxlsx, gtools, Biostrings, varhandle, data.table)


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

##############################################################################################################
# Lipase engineering database (LED DB) v.4.1 # 
##############################################################################################################
#The Lipase Engineering Database (LED) version4.1.0
#part of the BioCatNet
#This Internet database integrates information on sequence, structure and function of lipases and related proteins sharing the same a/b hydrolase fold.

#The BioCatNet Database System aims to collect and present comprehensive information about biocatalysts: sequence, structure, educts, products, environmental conditions and kinetics. Moreover, it will reveal structure-function relationships to boost development of novel biocatalysts.

#To contribute to BioCatNet by submission of functional information, please refer to the information in the BioCatNet Wiki.

#The LED was build based on the different architectures of α/β???hydrolases. All α/β???hydrolases consist of a catalytically active core domain (the α/β???hydrolase fold), containing the catalytic triad and the oxyanion hole. They might also contain additional structural modules, such as a lid, a cap, and N??? or C???terminal domains. The lid can emerge from five different positions: between β???strands β+1/β+2, β???1/β0, β???4/β???3, β+3/β+4, or between the N???terminus and β???3. Further, they can contain single caps, double caps, or N???terminal caps, and a N-terminal or a C-terminal domain. The presence of these structural modules determine the architecture of the α/β???hydrolase, resulting in 12 superfamilies. In addition, α/β???hydrolases can be distinguished by their oxyanion hole signature (GX???, GGGX???, or Y???type), which is indicated for each homologous family. For all homologous families with available structure information, the catalytic triad, the oxyanion hole residues, and the structural modules are annotated.
#The LED release 3.0 (2009) is still accessible at www.led.uni-stuttgart.de.

#https://wiki.biocatnet.de/?page=10
#[T#xxxxx] Taxon
#[S#xxxxx] Sequence
#[P#xxxxx] Protein
#[HF#cxxx] Homologous Family
#[HFG#xxx] Homologous Family Group
#[SF#xxxx] Superfamily
#[SFG#xxx] Superfamily Group

#[SF#1] Core #This superfamily contains α/β-hydrolases consisting only of the core domain.
#[SF#2] Lid (beta +1 & beta +2) #This superfamily contains α/β-hydrolases with an additional lid domain between β-strand +1 and +2. The lid is mobile and can open and close the entrance to the active site.
#[SF#3] Lid (beta -1 & beta 0) #This superfamily contains α/β-hydrolases with an additional lid domain between β-strand -1 and 0. The lid is mobile and can open and close the entrance to the active site.
#[SF#4] Lid (beta -4 & beta -3) #This superfamily contains α/β-hydrolases with an additional lid domain between β-strand -4 and -3. The lid is mobile and can open and close the entrance to the active site.
#[SF#5] Single cap #This superfamily contains α/β-hydrolases with an additional single cap domain that covers the active site and emerges between β-strand +1 and +2.
#[SF#6] Double cap #This superfamily contains α/β-hydrolases with two additional caps. One cap is covering the active site like the single cap and emerges between β-strand +1 and +2, while the second cap is covering the first cap from the N-terminal end of the protein.
#[SF#7] N-terminal cap #This superfamily contains α/β-hydrolases with an additional N-terminal cap domain that folds back over the protein and thereby covers the active site.
#[SF#8] N-terminal domain #This superfamily contains α/β-hydrolases with an additional N-terminal β-propeller domain. The propeller is formed by β-sheets that form a tunnel like structure on top of the protein.
#[SF#9] C-terminal domain #This superfamily contains α/β-hydrolases with an additional C-terminal β-sandwich domain. The C-terminal domain does not cover or interfere with the active site.
#[SF#10] Lid (beta +3 & beta +4) and C-terminal domain #This superfamily contains α/β-hydrolases with two additional domains: an additional C-terminal β-sandwich domain and an additional lid domain between β-strand +3 and +4 that can open and close the entrance to the active site.
#[SF#11] Lid (N-terminal) and C-terminal domain #This superfamily contains α/β-hydrolases with two additional domains: an additional C-terminal β-sandwich domain and an additional N-terminal lid domain that can open and close the entrance to the active site.
#[SF#12] Cap and N-terminal domain #This superfamily contains α/β-hydrolases with two additional domains: an additional N-terminal Rossmann fold domain and an additional single cap domain emerging between β-strand +1 and +2.
#[SF#13] Carboxylesterase-like #This superfamily contains α/β-hydrolases that are remote homologues to carboxylesterases

# Load LED DB data, providing sequence, structure, and functional information of lipases and related proteins.
'LED_DB' <- read.csv(file="led.SFID1+SFID2+SFID3+SFID4+SFID5+SFID6+SFID7+SFID8+SFID9+SFID10+SFID11+SFID12+SFID13.fasta", header=FALSE, sep = "\t", quote = "", fill=TRUE)

# Extract header lines from the LED_DB FASTA file and create a dataframe.
name <- setNames(as.data.frame(grep("^>", `LED_DB`$V1, value = TRUE)),"name")
name <- as.data.frame(name)
name$name <- as.character(name$name)

# Use regular expressions to extract specific annotations from the FASTA headers.
name <- name %>% mutate(`S# [sequence]`=str_extract(name,"sid\\|[0-9]{1,}\\|"))
name <- name %>% mutate(`P# [protein]`=str_extract(name,"pid\\|[0-9]{1,}\\|"))
name <- name %>% mutate(`HF# [homologous family]`=str_extract(name,"hfid\\|[0-9]{1,}\\|"))
name <- name %>% mutate(`SF# [superfamily]`=str_extract(name,"sfid\\|[0-9]{1,}\\|"))
name <- name %>% mutate(`taxonID`=str_extract(name,"gb|gi|emb|dbj\\|[0-9]{1,}\\|"))
name <- name %>% mutate(`taxonID`=str_extract(name,"taxonID\\|[0-9]{1,}\\|"))

# Separate the 'name' column into more specific annotations, removing unnecessary parts.
name<- tidyr::separate(name, col="name", sep = "sid\\|[0-9]{1,}\\|", into=c("vacant","name"), extra="merge", fill = "left") %>% select(-vacant)
name<- tidyr::separate(name, col="name", sep = "pid\\|[0-9]{1,}\\|", into=c("vacant","name"), extra="merge", fill = "left") %>% select(-vacant)
name<- tidyr::separate(name, col="name", sep = "hfid\\|[0-9]{1,}\\|", into=c("vacant","name"), extra="merge", fill = "left") %>% select(-vacant)
name<- tidyr::separate(name, col="name", sep = "sfid\\|[0-9]{1,}\\|", into=c("vacant","name"), extra="merge", fill = "left") %>% select(-vacant)
name<- tidyr::separate(name, col="name", sep = "taxonID\\|[0-9]{1,}\\|", into=c("protein_accession","taxonomy"), extra="merge", fill = "left")

# Identify sequence start and end positions in the LED_DB data frame.
# make concatenated fasta
seq_full<-grep("^>",`LED_DB`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,as.numeric(nrow(`LED_DB`)))

# Concatenate sequences and write to a CSV file.
sink("LED_DB_seq.csv");for(i in 1:as.numeric(nrow(name))){
  Sequence=paste(`LED_DB`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(Sequence,"\n")
}

# Read the concatenated sequences back into R, clean up, and add to the main data frame.
`LED_DB_sequence`<-read.csv("LED_DB_seq.csv",header=FALSE,sep="\n")
`LED_DB_sequence`$V1 <- gsub("\\s+", "",`LED_DB_sequence`$V1)
`LED_DB_sequence`$V1 <- gsub("\\W", "",`LED_DB_sequence`$V1)
names(`LED_DB_sequence`)[1] <- "sequence"

`LED_DB`<-cbind(name, `LED_DB_sequence`)
rm.all.but("LED_DB")

# Rearrange variables for consistency and clarity.
`LED_DB` <- arrange.vars(`LED_DB`, c("protein_accession"=1,
                                     "taxonomy"=2,
                                     "taxonID"=3,
                                     "SF# [superfamily]"=4,
                                     "HF# [homologous family]"=5,
                                     "P# [protein]"=6,
                                     "S# [sequence]"=7))

`LED_DB`$V1 <- gsub("\\W", "",`LED_DB`$V1)

# Clean up annotations by removing extraneous characters.
LED_DB$taxonID <- gsub("\\||taxonID", "",LED_DB$taxonID)
LED_DB$`SF# [superfamily]` <- gsub("\\||sfid", "",LED_DB$`SF# [superfamily]`)
LED_DB$`HF# [homologous family]` <- gsub("\\||hfid", "",LED_DB$`HF# [homologous family]`)
LED_DB$`P# [protein]` <- gsub("\\||pid", "",LED_DB$`P# [protein]`)
LED_DB$`S# [sequence]` <- gsub("\\||sid", "",LED_DB$`S# [sequence]`)

getwd()
Database_list <- list.files(pattern = "\\.tsv$")
sFam <- Database_list %>%  setNames(nm = gsub("\\.tsv","",.)) %>% map_df(~read.csv(.,sep="\t") %>% select(-1), .id = "sFam") 
names(sFam)[2] <- "#"
names(sFam)[3] <- "HF [homologous family]"
names(sFam)[4] <- "[HFG#xxx] Homologous Family Group"
sFam <- sFam %>% select(-data)

LED_DB <- merge(x=LED_DB, y=sFam, by.x="HF# [homologous family]", by.y="#", all.x=TRUE)
LED_DB_fasta <- LED_DB %>% dplyr::select(`S# [sequence]`,sequence) %>% mutate("S# [sequence]"=paste0(">",.$`S# [sequence]`))
write.table(x = LED_DB_fasta, file = "LED_DB_fasta.fasta", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\n")

getwd()
Database_list <- list.files(pattern = "gdbstructures\\.tsv$")
LED_structures <- Database_list %>%  setNames(nm = gsub("\\.tsv","",.)) %>% map_df(~read.csv(.,sep="\t"), .id = "Structure") 
names(LED_structures)[2] <- "SF# [superfamily]"
names(LED_structures)[3] <- "SF [superfamily]"
names(LED_structures)[4] <- "HF# [homologous family]"
names(LED_structures)[5] <- "HF [homologous family]"








#
for (counter in Database_list)
{
  assign(x = strsplit(x = counter,
                      split = "[-.]")[[1]][3], # for the pattern in your filenames
         value = readxl::read_xlsx(file = counter))
}

gg<-df_list_rbindlist %>% filter (locus_tag %in% F)

locustag<- select(a, contains('locus')) ##locus ??? ???????????? ??? 추출

##추출??? ????????? ??????????????? 각각 ????????? 분리
F<- F %>% unlist() %>% unique() 




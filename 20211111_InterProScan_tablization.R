pacman::p_load(ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr, seqinr, ggtree, magrittr, sjmisc, ggseqlogo, genoPlotR, gridExtra, seqinr, ape, rentrez, reutils, stringr, reshape2, xlsx, openxlsx, gtools, Biostrings, data.table, rhmmer, RFLPtools)

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

################################################################################################################################
                                                        # Bracket start #
################################################################################################################################
{

################################################################################################################################
                                                        # EggNOG emapper #
################################################################################################################################
  #http://eggnog-mapper.embl.de (eggNOG-mapper)

  #Load EggNOG mapper 
{ 
  EggNOG_dir<-getwd()
  EggNOG_file <- list.files(EggNOG_dir, pattern = "\\.annotations\\.xlsx$") # list
  
  for(EggNOG_list in 1:length(EggNOG_file)){
    Egg_temp <-  openxlsx::read.xlsx(paste0(EggNOG_file[EggNOG_list]), startRow = 3)  #name
    assign(paste(gsub("\\.annotations\\.xlsx","",EggNOG_file[EggNOG_list]),"eggNOG", sep="_"), Egg_temp[-((nrow(Egg_temp)-2):nrow(Egg_temp)),])
    rm(Egg_temp)
  }
  
  #  for(EggNOG_list in 3){
  #    assign(paste(gsub("\\.annotations\\.xlsx","",EggNOG_file[EggNOG_list]),"EggNOG_merged", sep="_"), 
  #           merge(get(paste(gsub("\\.annotations\\.xlsx","",EggNOG_file[EggNOG_list]))), get(paste(gsub("\\.annotations\\.xlsx","",EggNOG_file[EggNOG_list]),"eggNOG", sep="_")), 
  #                 by.x="locus_tag",by.y="Query", sort=FALSE))
  #  }
  
  annotation_list <- grep("emapper_eggNOG$",names(Filter(isTRUE, eapply(.GlobalEnv, is.data.frame))),value = TRUE)
  eggNOG <- get(annotation_list) %>% mutate(query=gsub("gnl\\|","",.$query))
  eggNOG <- eggNOG %>% mutate(query=gsub("\\|","\\:",.$query))
  eggNOG <- eggNOG %>% mutate(query=gsub("extdb:","",.$query))
  eggNOG <- eggNOG %>% mutate(query=gsub("lcl\\:(AC|NC|BG|NT|NW|NZ)\\_([a-zA-Z]{1,})?[0-9]{1,}\\.[0-9]{1,}_prot_","",.$query))
  
  #have to check after 2 rows when occurring errors
    #eggNOG <- eggNOG %>% mutate(query=gsub("_[0-9]{1,}$","",.$query))
    #eggNOG <- eggNOG %>% filter(!str_detect(query,"^[0-9]{1,}$"))
  
  
  #AC_	Genomic	Complete genomic molecule, usually alternate assembly
  #NC_	Genomic	Complete genomic molecule, usually reference assembly
  #NG_	Genomic	Incomplete genomic region
  #NT_	Genomic	Contig or scaffold, clone-based or WGSa
  #NW_	Genomic	Contig or scaffold, primarily WGSa
  #NZ_ Complete genomes and unfinished WGS data
  
  # Check #######################################################################################################################
  #Merge EggNOG mapper
  #contig_final <- merge(contig_final, eggNOG,by.x="protein_id",by.y="protein_id", sort=FALSE, all.x=TRUE)
  {if(grepl("(.)*_(RS)?[0-9]{1,}$",eggNOG$`query`[1])){
    assign("contig_final",merge(x=contig_final, y=eggNOG, by.x="locus_tag", by.y="query", sort=FALSE, all.x=TRUE))}
    else{assign("contig_final",merge(x=contig_final, y=eggNOG, by.x="protein_id", by.y="query", sort=FALSE, all.x=TRUE))
    }}
  contig_final <- contig_final %>% distinct(locus_tag, .keep_all = TRUE)
  rm(list=c("EggNOG_dir","EggNOG_file","EggNOG_list","annotation_list"))
}
################################################################################################################################
                                                           # InterProScan #
################################################################################################################################
  # Load InterProScan
{  
src_dir<-getwd()
  InterPro_search <- list.files(src_dir, pattern = "\\.tsv$") # list
  InterPro_search<- read.csv(InterPro_search, header = FALSE, quote = "", sep = "\t")
  InterPro_site <- list.files(src_dir, pattern = "\\.tsv.sites$") # list
  InterPro_site<- read.csv(InterPro_site, header = FALSE, quote = "", sep = "\t")
  #InterPro_site<- InterPro_site %>% arrange(., query)
  
  assign(paste("InterPro_site",sep="_"), setNames(as.data.frame(InterPro_site, ignore.case = TRUE),
                                                    c("query", #1
                                                      "Sequence MD5 digest", #2
                                                      "Sequence length", #3
                                                      "Identifier", #4
                                                      "Signature ac", #5
                                                      "rpsblast-site start", #6
                                                      "rpsblast-site end", #7
                                                      "numLocations", #8
                                                      "site-location residue", #9
                                                      "site-location residue start", #10
                                                      "site-location residue end", #11
                                                      "rpsblast-site description" #12
                                                      ))
         )
  
  
  
  for(i in (1:length(src_dir))){
    if(ncol(InterPro_search)==15){
      assign(paste("InterPro_search",sep="_"), setNames(as.data.frame(InterPro_search, ignore.case = TRUE),
                                                        c("query", "Sequence MD5 digest", "Sequence length", "Analysis", "Signature accession", "Signature description", "Start location", "Stop location", "Score", "Status", "Date", "InterPro annotations - accession", "InterPro annotations - description", "GO annotations", "Pathways annotations")))}
    else{assign(paste("InterPro_search",sep="_"), setNames(as.data.frame(InterPro_search, ignore.case = TRUE),
                                                           c("query", "Sequence MD5 digest", "Sequence length", "Analysis", "Signature accession", "Signature description", "Start location", "Stop location", "Score", "Status", "Date", "InterPro annotations - accession", "InterPro annotations - description")))}
  }
  
  # prepare output
  InterPro <- reshape2::melt(InterPro_search, id.vars=c("query","Analysis"), measure.vars=c("Signature accession","Signature description","Score"))
  InterPro <- reshape2::dcast(InterPro, formula = `query`  ~ variable + Analysis, fun.aggregate = toString)
  InterPro <- setcolorder(InterPro, order(sub('.*_', '', names(InterPro))))
  InterPro <- arrange.vars(InterPro, c("query"=1))
  InterPro$`query` <- gsub("gnl\\|","",InterPro$`query`)
  InterPro$`query` <- gsub("\\|","\\:",InterPro$`query`)
  InterPro$`query` <- gsub("extdb:","",InterPro$`query`)
  # Check #######################################################################################################################
#  contig_final <- merge(x=contig_final, y=InterPro, by.x="protein_id", by.y="query", sort=FALSE, all.x=TRUE)
#  contig_final <- merge(x=contig_final, y=InterPro, by.x="locus_tag", by.y="query", sort=FALSE, all.x=TRUE)
  {if(grepl("(.)*_(RS)?[0-9]{1,}$",InterPro$`query`[1])){
        assign("contig_final",merge(x=contig_final, y=InterPro, by.x="locus_tag", by.y="query", sort=FALSE, all.x=TRUE))}
    else{assign("contig_final",merge(x=contig_final, y=InterPro, by.x="protein_id", by.y="query", sort=FALSE, all.x=TRUE))
  }}
}
  
#[TODO] GO and Reactome aggregation
  ################################################################################################################################
    # save xlsx #
  ################################################################################################################################

  {
    wb <- createWorkbook()
    addWorksheet(wb, "1.GenBank_table")
    addWorksheet(wb, "2.InterPro_site")
    addWorksheet(wb, "3.InterPro_search")
    addWorksheet(wb, "4.Codon_usage")
    
    writeData(wb, "1.GenBank_table", contig_final, startRow = 1, startCol = 1)
    writeData(wb, "2.InterPro_site", InterPro_site, startRow = 1, startCol = 1)
    writeData(wb, "3.InterPro_search", InterPro_search, startRow = 1, startCol = 1)
    writeData(wb, "4.Codon_usage", codon, startRow = 1, startCol = 1)
    
    saveWorkbook(wb, file = "total.xlsx", overwrite = TRUE)  
  }
  
  
  
  
  
  
  
  
  ################################################################################################################################
  # Load Database #
  ################################################################################################################################
  {
    working_dir<-getwd()
    #go to the database directory
    setwd("../../DataBase")
    #file_names=as.list(dir(pattern="*\\.RData$"))
    #Database_list <- list.files(src_dir, pattern = "\\.RData$") # list
    Database_list <- list.files(getwd(), pattern = "\\.RData$") # list
    lapply(Database_list, load, environment())
    rm(list=c("Database_list","src_dir","sFam","MEROPS-MPEP","MEROPS-MPRO",ls(pattern="fasta$")))
    #back to the working directory
    setwd(working_dir)
  }
  
  
  ################################################################################################################################
                                                              # PSORTb #
  ################################################################################################################################
  # Load PSORTb
{  
src_dir <- getwd()
  PSORTb <- list.files(src_dir, pattern = "PSORTb\\.txt$") # list
  PSORTb <- read.csv(PSORTb, header = TRUE, quote = "", sep = "\t", row.names = NULL)  
  names(PSORTb)[1] <- "query"
  PSORTb$query <- substr(PSORTb$query, regexpr("^gln",PSORTb$query), regexpr(" ",PSORTb$query)-1)
  PSORTb$query <- gsub("gnl\\|","", PSORTb$query)
  PSORTb$query <- gsub("\\|","\\:", PSORTb$query)
  PSORTb$query <- gsub("extdb:","", PSORTb$query)
  PSORTb <- PSORTb %>% dplyr::select("query", "ModHMM._Localization", "SCL.BLAST._Localization", "Signal._Localization", "Extracellular_Score", "Final_Localization_Details")
  
  # Check #######################################################################################################################
  # Load PSORTb
  #contig_final <- merge(x=contig_final, y=PSORTb, by.x="protein_id", by.y="query", sort=FALSE, all.x=TRUE)
  #contig_final <- merge(x=contig_final, y=PSORTb, by.x="locus_tag", by.y="query", sort=FALSE, all.x=TRUE)
  {if(grepl("(.)*_(RS)?[0-9]{1,}$",PSORTb$`query`[1])){
    assign("contig_final",merge(x=contig_final, y=PSORTb, by.x="locus_tag", by.y="query", sort=FALSE, all.x=TRUE))}
    else{assign("contig_final",merge(x=contig_final, y=PSORTb, by.x="protein_id", by.y="query", sort=FALSE, all.x=TRUE))
    }}
}
  ################################################################################################################################
                                                              # MEROPS #
  ################################################################################################################################
  # Load MEROPS-MP
{ 
 src_dir <- getwd()
  `MEROPS-MP_hit` <- list.files(src_dir, pattern = "MEROPS_MP\\.txt$") # list
  `MEROPS-MP_hit` <- read.csv(`MEROPS-MP_hit`, header = FALSE, quote = "", sep = "\t", row.names = NULL)  
  names(`MEROPS-MP_hit`) <- c("subject","acc.ver","% identity","alignment length","mismatches","gqp opens","q. start","q. end","s. start","s. end","evalue","bit score")
  
  `MEROPS-MP_hit`$subject <- gsub("gnl\\|","", `MEROPS-MP_hit`$subject)
  `MEROPS-MP_hit`$subject <- gsub("\\|","\\:", `MEROPS-MP_hit`$subject)
  `MEROPS-MP_hit`$subject <- gsub("extdb:","", `MEROPS-MP_hit`$subject)
  # Check #######################################################################################################################
  `MEROPS-MP_hit` <- `MEROPS-MP_hit` %>% mutate(subject=gsub("lcl\\:(AC|NC|BG|NT|NW|NZ)\\_([a-zA-Z])?{1,}[0-9]{1,}\\.[0-9]{1,}_prot_","",.$subject))
  `MEROPS-MP_hit` <- `MEROPS-MP_hit` %>% mutate(subject=gsub("_[0-9]{1,}$","",.$subject))
  `MEROPS-MP_hit` <- `MEROPS-MP_hit` %>% filter(!str_detect(subject,"^[0-9]{1,}$")) %>% dplyr::rename(., "subject" = subject)
  
  
  
  # Extract the best hit
  `MEROPS-MP_best_hit` <- `MEROPS-MP_hit` %>% group_by(subject) %>% 
    arrange(evalue, desc(`bit score`), desc(`alignment length`), desc(`% identity`)) %>% filter(row_number() <= 1) %>% ungroup()
  `MEROPS-MP_best_hit` <- `MEROPS-MP_best_hit` %>% filter(`% identity` > 50.0)
  
  `MEROPS-MP_best_hit` <- merge(x=`MEROPS-MP_best_hit`, y=`MEROPS-MP`, by.x="acc.ver", by.y="MEROPS_ACCESSION", suffixes = c("_merops", ""), sort=FALSE, all.x=TRUE)
  # Check #######################################################################################################################
  #`MEROPS-MP_best_hit` <- merge(x=`MEROPS-MP_best_hit`, y=contig_final, by.x="subject", by.y="protein_id", sort=FALSE, all.x=TRUE)
  {if(grepl("(.)*_(RS)?[0-9]{1,}$",`MEROPS-MP_best_hit`$`subject`[1])){
    assign("MEROPS-MP_best_hit",merge(x=`MEROPS-MP_best_hit`, y=contig_final, by.x="subject", by.y="locus_tag", suffixes = c("_merops", ""), sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "locus_tag" = subject))}
    else{assign("MEROPS-MP_best_hit",merge(x=`MEROPS-MP_best_hit`, y=contig_final, by.x="subject", by.y="protein_id", suffixes = c("_merops", ""), sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "protein_id" = subject))
    }}
  
  `MEROPS-MP_best_hit` <- arrange.vars(`MEROPS-MP_best_hit`, c("locus_tag"=1,
                                                               "protein_id"=2,
                                                               "product"=3,
                                                               "Protease"=4,
                                                               "organism"=5,
                                                               "suborganism"=6,
                                                               "MEROPS_ID"=7,
                                                               "% identity"=8,
                                                               "alignment length"=9,
                                                               "mismatches"=10,
                                                               "gqp opens"=11,
                                                               "q. start"=12,
                                                               "q. end"=13,
                                                               "s. start"=14,
                                                               "s. end"=15,
                                                               "evalue_merops"=16,
                                                               "bit score"=17,
                                                               "AA_sequence"=18))
}
  
  # Need to filter manually
  # Save as xlsx file to curate manually
  write.xlsx(`MEROPS-MP_best_hit`, file="MEROPS-MP_best_hit_before_filter.xlsx", sheetName="before_filter", row.names=FALSE, overwrite = TRUE)
  #read.table()
  
  
  ################################################################################################################################
                                                           # dbCAN (CAZy) #
  ################################################################################################################################
  # Load dbCAN
  {
  src_dir <- getwd()
  dbCAN <- list.files(src_dir, pattern = "dbCAN\\.txt$") # list
  dbCAN <- read.csv(dbCAN, header = TRUE, quote = "", sep = "\t", row.names = NULL)  
  names(dbCAN) <- paste("dbCAN", names(dbCAN))
  names(dbCAN)[1] <- "query"
  dbCAN$query <- gsub("gnl\\|","", dbCAN$query)
  dbCAN$query <- gsub("\\|","\\:", dbCAN$query)
  dbCAN$query <- gsub("extdb:","", dbCAN$query)
  # Merge dbCAN
  # Check #######################################################################################################################
  #dbCAN <- merge(x=dbCAN, y=contig_final, by.x="query", by.y="protein_id", sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "protein_id" = query)
  #dbCAN <- merge(x=dbCAN, y=contig_final, by.x="query", by.y="locus_tag", sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "locus_tag" = query)
  {if(grepl("(.)*_(RS)?[0-9]{1,}$",dbCAN$query[1])){
    assign("dbCAN",merge(x=dbCAN, y=contig_final, by.x="query", by.y="locus_tag", suffixes = c("_dbCAN", ""), sort=FALSE, all.x=TRUE)%>% dplyr::rename(., "locus_tag" = query))}
    else{assign("dbCAN",merge(x=dbCAN, y=contig_final, by.x="query", by.y="protein_id", suffixes = c("_dbCAN", ""), sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "protein_id" = query))
    }}

  
  {
    dbCAN <- dbCAN %>% mutate(`dbCAN HMMER`=gsub("\\([0-9]{1,}\\-[0-9]{1,}\\)","", .$`dbCAN HMMER`))
    dbCAN <- dbCAN %>% mutate_at(vars(`CAZy`, `dbCAN HMMER`, `dbCAN Hotpep`, `dbCAN DIAMOND`, `dbCAN Signalp`), funs(gsub("^N|\\-$", "", .)))
    dbCAN <- dbCAN %>% mutate_at(vars(`CAZy`, `dbCAN HMMER`, `dbCAN Hotpep`, `dbCAN DIAMOND`, `dbCAN Signalp`), funs(gsub(" ", "", .)))
    dbCAN <- dbCAN %>% mutate_at(vars(`CAZy`, `dbCAN HMMER`, `dbCAN Hotpep`, `dbCAN DIAMOND`, `dbCAN Signalp`), funs(gsub(",", "+", .)))
    dbCAN <- dbCAN %>% unite("CAZy Family", c(`CAZy`,`dbCAN HMMER`, `dbCAN Hotpep`, `dbCAN DIAMOND`), sep="+", na.rm = TRUE, remove = FALSE)
    dbCAN <- dbCAN %>% mutate_at(vars(`CAZy Family`), funs(gsub("^\\+{1,}|\\+{1,}$", "", .))) # delete additionally appended prefix and suffix
    dbCAN <- dbCAN %>% mutate_at(vars(`CAZy Family`), funs(gsub("([\\+])\\1+", "\\1", ., perl=TRUE))) # repeated character replacement
    dbCAN <- dbCAN %>% mutate(`CAZy Family` = strsplit(.$`CAZy Family`,"\\+") %>% Map(unique,.)) %>% unnest(cols = `CAZy Family`) %>% mutate(`CAZy Family` = qdap::beg2char(.$`CAZy Family`, "_")) %>% distinct(locus_tag, `CAZy Family`, .keep_all = TRUE)
    dbCAN <- dbCAN %>% mutate(., "CAZy Class" = case_when(grepl("^GH", .$`CAZy Family`) ~ "GH",
                                                          grepl("^GT", .$`CAZy Family`) ~ "GT",
                                                          grepl("^PL", .$`CAZy Family`) ~ "PL",
                                                          grepl("^CE", .$`CAZy Family`) ~ "CE",
                                                          grepl("^AA", .$`CAZy Family`) ~ "AA",
                                                          grepl("^CBM", .$`CAZy Family`) ~ "CBM")
    )
  }
  dbCAN <- arrange.vars(dbCAN, c("locus_tag"=1,
                                 "protein_id"=2,
                                 "product"=3,
                                 "CAZy Class"=4,
                                 "CAZy Family"=5,
                                 "CAZy"=6,
                                 "dbCAN HMMER"=7,
                                 "dbCAN Hotpep"=8,
                                 "dbCAN DIAMOND"=9,
                                 "dbCAN Signalp"=10,
                                 "dbCAN X.ofTools"=11)
                        )
  
  write.xlsx(`dbCAN`, file="dbCAN_best_hit.xlsx", sheetName="dbCAN", row.names=FALSE)
  
  
  # Load minor dbCAN database
    # hmmer
  src_dir <- getwd()
  dbCAN_hmmer <- list.files(src_dir, pattern = "hmmer\\.out\\.txt$") # list
  dbCAN_hmmer <- read.table(file=dbCAN_hmmer, sep="\t", header = TRUE)
  dbCAN_hmmer <- arrange.vars(dbCAN_hmmer, c("Gene.ID"=1))
  names(dbCAN_hmmer) <- paste("dbCAN", names(dbCAN_hmmer))
  names(dbCAN_hmmer)[1] <- "query"
  dbCAN_hmmer$query <- gsub("gnl\\|","", dbCAN_hmmer$query)
  dbCAN_hmmer$query <- gsub("\\|","\\:", dbCAN_hmmer$query)
  dbCAN_hmmer$query <- gsub("extdb:","", dbCAN_hmmer$query)
  
    # Hotpep
  src_dir <- getwd()
  dbCAN_Hotpep <- list.files(src_dir, pattern = "Hotpep\\.out\\.txt$") # list
  dbCAN_Hotpep <- read.csv(dbCAN_Hotpep, header = TRUE, quote = "", sep = "\t", row.names = NULL)  
  dbCAN_Hotpep <- arrange.vars(dbCAN_Hotpep, c("Gene.ID"=1))
  names(dbCAN_Hotpep) <- paste("dbCAN_Hotpep", names(dbCAN_Hotpep))
  names(dbCAN_Hotpep)[1] <- "query"
  dbCAN_Hotpep$query <- gsub("gnl\\|","", dbCAN_Hotpep$query)
  dbCAN_Hotpep$query <- gsub("\\|","\\:", dbCAN_Hotpep$query)
  dbCAN_Hotpep$query <- gsub("extdb:","", dbCAN_Hotpep$query)
  
    # diamond
  src_dir<-getwd()
  dbCAN_diamond <- list.files(src_dir, pattern = "diamond\\.out\\.txt$") # list
  dbCAN_diamond <- read.csv(dbCAN_diamond, header = TRUE, quote = "", sep = "\t", row.names = NULL)  
  dbCAN_diamond <- arrange.vars(dbCAN_diamond, c("Gene.ID"=1))
  names(dbCAN_diamond) <- paste("dbCAN_diamond", names(dbCAN_diamond))
  names(dbCAN_diamond)[1] <- "query"
  dbCAN_diamond$query <- gsub("gnl\\|","", dbCAN_diamond$query)
  dbCAN_diamond$query <- gsub("\\|","\\:", dbCAN_diamond$query)
  dbCAN_diamond$query <- gsub("extdb:","", dbCAN_diamond$query)
}
  
  ################################################################################################################################
                                                          # LED #
  ################################################################################################################################
  # Blastp condition (evalue 1e-3, max_target_seqs 30)
{
  src_dir <- getwd()
  LED <- list.files(src_dir, pattern = "LED\\.txt$") # list
  #LED <- RFLPtools::read.blast(LED, sep = "\t")
  LED <- read.csv(LED, header = FALSE, quote = "", sep = "\t", row.names = NULL)  
  names(LED) <- c("subject","acc.ver","% identity","alignment length","mismatches","gqp opens","q. start","q. end","s. start","s. end","evalue","bit score")
  
  LED$subject <- gsub("gnl\\|","", LED$subject)
  LED$subject <- gsub("\\|","\\:", LED$subject)
  LED$subject <- gsub("extdb:","", LED$subject)
  
  # Extract the best hit
  LED_best_hit <- LED %>% group_by(subject) %>% top_n(1)
  LED_best_hit <- LED_best_hit %>% group_by(subject) %>% 
    arrange(evalue, desc(`bit score`), desc(`alignment length`), desc(`% identity`)) %>% filter(row_number() <= 1) %>% ungroup()
  LED_best_hit %>% filter(`% identity` > 40.0)
  
  LED_best_hit <- merge(x=LED_best_hit, y=LED_DB, by.x="acc.ver", by.y="S# [sequence]", sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "S# [sequence]" = 'acc.ver')
  
  # Check #######################################################################################################################
  #LED_best_hit <- merge(x=LED_best_hit, y=contig_final, by.x="subject", by.y="protein_id", sort=FALSE, all.x=TRUE)
  #LED_best_hit <- merge(x=LED_best_hit, y=contig_final, by.x="subject", by.y="locus_tag", sort=FALSE, all.x=TRUE)
  {if(grepl("(.)*_(RS)?[0-9]{1,}$",LED_best_hit$`subject`[1])){
    assign("LED_best_hit",merge(x=LED_best_hit, y=contig_final, by.x="subject",by.y="locus_tag", suffixes = c("_led", ""), sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "locus_tag" = subject))}
    else{assign("LED_best_hit",merge(x=LED_best_hit, y=contig_final,by.x="subject", by.y="protein_id", suffixes = c("_led", ""), sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "protein_id" = subject))
    }}  
  LED_best_hit <- arrange.vars(LED_best_hit, c("protein_id"=1,
                                               "product"=2))
  LED_best_hit <- LED_best_hit %>% arrange(evalue_led, desc(`bit score`), desc(`alignment length`), desc(`% identity`))
  # Need to filter manually
  # Save as xlsx file to curate manually
  write.xlsx(`LED_best_hit`, file="LED_best_hit_before_filter.xlsx", sheetName="before_filter", row.names=FALSE)
}
  ################################################################################################################################
                                                              # REBASE #
  ################################################################################################################################
  # Blastp condition (evalue 1e-3, max_target_seqs 30)
{
  src_dir <- getwd()
  REBASE <- list.files(src_dir, pattern = "REBASE\\.txt$") # list
  #REBASE <- RFLPtools::read.blast(REBASE, sep = "\t")
  REBASE <- read.csv(REBASE, header = FALSE, quote = "", sep = "\t", row.names = NULL)  
  names(REBASE) <- c("subject","acc.ver","% identity","alignment length","mismatches","gqp opens","q. start","q. end","s. start","s. end","evalue","bit score")
  
  REBASE$subject <- gsub("gnl\\|","", REBASE$subject)
  REBASE$subject <- gsub("\\|","\\:", REBASE$subject)
  REBASE$subject <- gsub("extdb:","", REBASE$subject)
  
  # Extract the best hit
  REBASE_best_hit <- REBASE %>% group_by(subject) %>% top_n(1)
  REBASE_best_hit <- REBASE_best_hit %>% group_by(subject) %>% 
    arrange(evalue, desc(`bit score`), desc(`alignment length`), desc(`% identity`)) %>% filter(row_number() <= 1) %>% ungroup()
  REBASE_best_hit %>% filter(`% identity` > 40.0)
  
  REBASE_best_hit <- merge(x=REBASE_best_hit, y=REBASE_DB_protein_gold_seqs, by.x="acc.ver", by.y="REBASE", sort=FALSE, all.x=TRUE)
  REBASE_best_hit <- arrange.vars(REBASE_best_hit, c("subject"=1,
                                                     "acc.ver"=2,
                                                     "EnzType"=3,
                                                     "RecSeq"=4,
                                                     "% identity"=5,
                                                     "alignment length"=6,
                                                     "mismatches"=7,
                                                     "gqp opens"=8,
                                                     "q. start"=9,
                                                     "q. end"=10,
                                                     "s. start"=11,
                                                     "s. end"=12,
                                                     "evalue"=13,
                                                     "bit score"=14))
  # Check #######################################################################################################################
  #REBASE_best_hit <- merge(x=REBASE_best_hit, y=contig_final, by.x="subject", by.y="protein_id", sort=FALSE, all.x=TRUE)
  #REBASE_best_hit <- merge(x=REBASE_best_hit, y=contig_final, by.x="subject", by.y="locus_tag", sort=FALSE, all.x=TRUE)
  {if(grepl("(.)*_(RS)?[0-9]{1,}$",REBASE_best_hit$`subject`[1])){
    assign("REBASE_best_hit",merge(x=REBASE_best_hit, y=contig_final, by.x="subject", by.y="locus_tag", suffixes = c("_rebase", ""), sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "locus_tag"="subject"))}
    else{assign("REBASE_best_hit",merge(x=REBASE_best_hit, y=contig_final, by.x="subject", by.y="protein_id", suffixes = c("_rebase", ""), sort=FALSE, all.x=TRUE) %>% dplyr::rename(., "protein_id"="subject"))
    }}  
  REBASE_best_hit <- REBASE_best_hit %>% arrange(evalue_rebase, desc(`bit score`), desc(`alignment length`), desc(`% identity`))
  # Need to filter manually
  # Save as xlsx file to curate manually
  REBASE_best_hit <- arrange.vars(REBASE_best_hit, c("locus_tag"=1,
                                                     "protein_id"=2,
                                                     "product"=3)
                                  )
  write.xlsx(`REBASE_best_hit`, file="REBASE_best_hit_before_filter.xlsx", sheetName="before_filter", row.names=FALSE)
}
  
  ################################################################################################################################
                                                            # Bracket end #
  ################################################################################################################################
  rm(list=c("eggNOG","EggNOG_file","EggNOG_list",
            "PSORTb",
            "LED_DB","LED_structures",
            "REBASE_DB_protein_gold_seqs",
            "DoriC",
            "MEROPS-MP",
            "counter",
            "i",
            ls(pattern="dir$")))
  }
  
  
  ################################################################################################################################
                                                        # Write xlsx format #
##################################################################################################################################
{
wb <- createWorkbook()
addWorksheet(wb, "1.GenBank_table")
addWorksheet(wb, "2.InterPro_site")
addWorksheet(wb, "3.InterPro_search")
addWorksheet(wb, "4.Codon_usage")

writeData(wb, "1.GenBank_table", contig_final, startRow = 1, startCol = 1)
writeData(wb, "2.InterPro_site", InterPro_site, startRow = 1, startCol = 1)
writeData(wb, "3.InterPro_search", InterPro_search, startRow = 1, startCol = 1)
writeData(wb, "4.Codon_usage", codon, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "total.xlsx", overwrite = TRUE)  
}

#openxlsx::write.xlsx(contig_final, file="filename.xlsx", sheetName="1. InterPro", row.names=FALSE, append = TRUE)
#openxlsx::write.xlsx(InterPro_site, file="filename.xlsx", sheetName="2. InterPro", row.names=FALSE, append = TRUE,overwrite = TRUE)
#openxlsx::write.xlsx(InterPro_search, file="filename.xlsx", sheetName="3. InterPro", row.names=FALSE, append = TRUE)
#xlsx::write.xlsx(`MEROPS-MP_best_hit_filtered`, file="filename.xlsx", sheetName="2.MEROPS", append=TRUE, row.names=FALSE)
#xlsx::write.xlsx(LED_best_hit_filtered, file="filename.xlsx", sheetName="3.LED", append=TRUE, row.names=FALSE)
#xlsx::write.xlsx(dbCAN, file="filename.xlsx", sheetName="4. CAZy", append=TRUE, row.names=FALSE)
#xlsx::write.xlsx(dbCAN_diamond, file="filename.xlsx", sheetName="5. CAZy_diamond", append=TRUE, row.names=FALSE)
#xlsx::write.xlsx(dbCAN_hmmer, file="filename.xlsx", sheetName="6. CAZy_hmmer", append=TRUE, row.names=FALSE)
#xlsx::write.xlsx(dbCAN_Hotpep, file="filename.xlsx", sheetName="7. CAZy_Hotpep", append=TRUE, row.names=FALSE)
#xlsx::write.xlsx(REBASE_best_hit_filtered, file="filename.xlsx", sheetName="8. REBASE", append=TRUE, row.names=FALSE)
  
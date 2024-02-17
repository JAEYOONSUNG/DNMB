pacman::p_load(qdap, qdapRegex, data.table, tidyr, plyr, dplyr, rhmmer, stringr, seqinr,magrittr, sjmisc, purrr,
               reutils, reshape2, xlsx, openxlsx, gtools, Biostrings, varhandle,
               ggpubr, gtable, RFLPtools, ggplot2, ggnewscale, GGally, ggtree, ggseqlogo, genoPlotR, 
               gridExtra, cowplot, ComplexHeatmap, ggdendro, patchwork, forcats)

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
                                              # MEROPS-MPRO DB #
##############################################################################################################
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


# make concatenated fasta
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

################################################################################################################################################
#Blast with my data
src_dir<-getwd()
`MEROPS-MPRO_results` <- list.files(src_dir, pattern = "MEROPS\\-MPRO\\.txt$") # list
`MEROPS-MPRO_results` <- read.csv(file=`MEROPS-MPRO_results`, header=FALSE)
`MEROPS-MPRO_results` <- `MEROPS-MPRO_results` %>% filter(!grepl("^#", `MEROPS-MPRO_results`$V1))
`MEROPS-MPRO_results` <- separate(`MEROPS-MPRO_results`, col="V1",into=c("subject","acc.ver","% identity","alignment length","mismatches","gqp opens","q. start","q. end","s. start","s. end","evalue","bit score"),sep="\t")
`MEROPS-MPRO_results` <- `MEROPS-MPRO_results` %>% filter(!is.na(acc.ver))
`MEROPS-MPRO_results` <- merge(x=`MEROPS-MPRO_results`, y=`MEROPS-MPRO`, by.x="acc.ver", by.y="MEROPS_ACCESSION", all.x=TRUE)

write.xlsx(x = `MEROPS-MPRO_results`, file = "MEROPS-MPRO_results.xlsx", quote = FALSE, row.names = FALSE, col.names = TRUE)

##############################################################################################################
                                                  # MEROPS-MPEP DB #
##############################################################################################################
'MEROPS-MPEP' <- read.csv(file="pepunit.lib.txt", header=FALSE, sep = "\t", quote = "", fill=TRUE, fileEncoding="latin1")

name <- setNames(as.data.frame(grep("^>", `MEROPS-MPEP`$V1, value = TRUE)),"name")
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
name<- tidyr::separate(name, col="Protease", sep = "\\((?=[^\\(]*$)", into=c("protease","suborganism"), extra="merge", fill = "left")
name<- tidyr::separate(name, col="protease", sep = "(?>\\()", into=c("protease","organism"), extra="merge", fill = "right")
name$organism <- gsub("\\)","",name$organism)
name$suborganism <- gsub("\\)","",name$suborganism)
name$MEROPS_ACCESSION <- gsub(">","",name$MEROPS_ACCESSION)
name<- name %>% mutate_all(na_if,"")
name[is.na(name)]  <- ""
name<- as.data.frame(apply(name,2,function(x)gsub('\\s+$', '',x)))


# make concatenated fasta
seq_full<-grep("^>",`MEROPS-MPEP`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,as.numeric(nrow(`MEROPS-MPEP`)))
sink("MEROPS-MPEP_seq.csv");for(i in 1:as.numeric(nrow(name))){
  Sequence=paste(`MEROPS-MPEP`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(Sequence,"\n")
}

`MEROPS-MPEP_sequence`<-read.csv("MEROPS-MPEP_seq.csv",header=FALSE,sep="\n")
`MEROPS-MPEP_sequence`$V1 <- gsub("\\s+", "",`MEROPS-MPEP_sequence`$V1)
`MEROPS-MPEP_sequence`$V1 <- gsub("\\W", "",`MEROPS-MPEP_sequence`$V1)
names(`MEROPS-MPEP_sequence`)[1] <- "AA_sequence"

`MEROPS-MPEP`<-cbind(name, `MEROPS-MPEP_sequence`)

rm(`MEROPS-MPEP_sequence`)
rm(name)

write.xlsx(`MEROPS-MPEP`,"MEROPS-MPEP_DB.xlsx")  


`MEROPS-MPEP_fasta` <- `MEROPS-MPEP` %>% select(MEROPS_ACCESSION, AA_sequence)
`MEROPS-MPEP_fasta`$MEROPS_ACCESSION<- paste0(">",`MEROPS-MPEP_fasta`$MEROPS_ACCESSION)
write.table(x = `MEROPS-MPEP_fasta`, file = "MEROPS-MPEP.fasta", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\n")

#Blast with my data
src_dir<-getwd()
`MEROPS-MPEP_results` <- list.files(src_dir, pattern = "MEROPS\\-MPEP\\.txt$") # list
`MEROPS-MPEP_results` <- read.csv(file=`MEROPS-MPEP_results`, header=FALSE)
`MEROPS-MPEP_results` <- `MEROPS-MPEP_results` %>% filter(!grepl("^#", `MEROPS-MPEP_results`$V1))
`MEROPS-MPEP_results` <- separate(`MEROPS-MPEP_results`, col="V1",into=c("subject","acc.ver","% identity","alignment length","mismatches","gqp opens","q. start","q. end","s. start","s. end","evalue","bit score"),sep="\t")
`MEROPS-MPEP_results` <- `MEROPS-MPEP_results` %>% filter(!is.na(acc.ver))
`MEROPS-MPEP_results` <- merge(x=`MEROPS-MPEP_results`, y=`MEROPS-MPEP`, by.x="acc.ver", by.y="MEROPS_ACCESSION", all.x=TRUE)

write.xlsx(x = `MEROPS-MPEP_results`, file = "MEROPS-MPEP_results.xlsx", quote = FALSE, row.names = FALSE, col.names = TRUE)



##############################################################################################################
                                                  # MEROPS-MP DB #
##############################################################################################################
'MEROPS-MP' <- read.csv(file="protease.lib.txt", header=FALSE, sep = "\t", quote = "", fill=TRUE, fileEncoding="latin1")
'MEROPS-MP' <- read.table(file="protease.lib.txt", header=FALSE, sep = "\t", quote = "", fill=TRUE, fileEncoding="latin1")


name <- setNames(as.data.frame(grep("^>", `MEROPS-MP`$V1, value = TRUE)),"name")
name <- as.data.frame(name)
name$name <- as.character(name$name)

#IF not working separate with "Error in nchar(u, itype) : invalid multibyte string, element 1"
#Sys.setlocale(category ="LC_ALL", locale="US")

name<- tidyr::separate(name, col="name", sep = " \\- ", into=c("MEROPS_ACCESSION","V1"), extra="merge", fill = "left")
name<- tidyr::separate(name, col="V1", sep = "\\[", into=c("Protease","MEROPS_ID"), extra="merge", fill = "left")
name$MEROPS_ID <- gsub("\\]","",name$MEROPS_ID)

name<- tidyr::separate(name, col="Protease", sep = "\\((?=[^\\(]*$)", into=c("Protease","suborganism"), extra="merge", fill = "left")
name<- tidyr::separate(name, col="Protease", sep = "(?>\\()", into=c("Protease","organism"), extra="merge", fill = "right")
name$organism <- gsub("\\)","",name$organism)
name$suborganism <- gsub("\\)","",name$suborganism)
name$MEROPS_ACCESSION <- gsub(">","",name$MEROPS_ACCESSION)
name<- name %>% mutate_all(na_if,"")
name[is.na(name)]  <- ""
name<- as.data.frame(apply(name,2,function(x)gsub('\\s+$', '',x)))


# make concatenated fasta
seq_full<-grep("^>",`MEROPS-MP`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,as.numeric(nrow(`MEROPS-MP`)))
sink("MEROPS-MP_seq.csv");for(i in 1:as.numeric(nrow(name))){
  Sequence=paste(`MEROPS-MP`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(Sequence,"\n")
}

`MEROPS-MP_sequence`<-read.csv("MEROPS-MP_seq.csv",header=FALSE,sep="\n")
`MEROPS-MP_sequence`$V1 <- gsub("\\s+", "",`MEROPS-MP_sequence`$V1)
`MEROPS-MP_sequence`$V1 <- gsub("\\W", "",`MEROPS-MP_sequence`$V1)
names(`MEROPS-MP_sequence`)[1] <- "AA_sequence"

`MEROPS-MP`<-cbind(name, `MEROPS-MP_sequence`)

rm(`MEROPS-MP_sequence`)
rm(name)

write.xlsx(`MEROPS-MP`,"MEROPS-MP_DB.xlsx")  

`MEROPS-MP_fasta` <- `MEROPS-MP` %>% select(MEROPS_ACCESSION, AA_sequence)
`MEROPS-MP_fasta`$MEROPS_ACCESSION<- paste0(">",`MEROPS-MP_fasta`$MEROPS_ACCESSION)
write.table(x = `MEROPS-MP_fasta`, file = "MEROPS-MP.fasta", quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\n")



##############################################################################################################################


src_dir<-getwd()
`MEROPS-MP_best_hit_filtered` <- list.files(src_dir, pattern = "MEROPS\\-MP_best_hit_after_filter\\.txt$") # list
`MEROPS-MP_best_hit_filtered` <- read.csv(file=`MEROPS-MP_best_hit_filtered`, header=FALSE)

#`MEROPS-MP_results` <- `MEROPS-MP_results` %>% filter(!grepl("^#", `MEROPS-MP_results`$V1))
#`MEROPS-MP_results` <- separate(`MEROPS-MP_results`, col="V1",into=c("subject","acc.ver","% identity","alignment length","mismatches","gqp opens","q. start","q. end","s. start","s. end","evalue","bit score"),sep="\t")
#`MEROPS-MP_results` <- `MEROPS-MP_results` %>% filter(!is.na(acc.ver))
#`MEROPS-MP_results` <- merge(x=`MEROPS-MP_results`, y=`MEROPS-MP`, by.x="acc.ver", by.y="MEROPS_ACCESSION", all.x=TRUE)
#`MEROPS-MP_best_hit` <- `MEROPS-MP_results` %>% group_by(subject) %>% filter(`% identity` == max(`% identity`)) %>% arrange(subject) %>% ungroup()

MEROPS-MP_best_hit_filtered %>% dplyr::select(c("locus_tag",names(`MEROPS-MP_hit`)))

`MEROPS-MP_best_hit_filtered` <- merge(`MEROPS-MP_best_hit_filtered`, contig_final, by="protein_id", all.x=TRUE)
`MEROPS-MP_best_hit_filtered` <- arrange.vars(`MEROPS-MP_best_hit_filtered`, c("gene"=15,
                                                             "product"=16,
                                                             "amino_acid"=17,
                                                             "Description"=18))



##############################################################################################################################
                                                    ######Visualization#######
##############################################################################################################################

# data processing to plot MEROPS (SUB)TYPE
`MEROPS-MP_best_hit` <- `MEROPS-MP_best_hit` %>% mutate(MEROPS_TYPE = substr(MEROPS_ID, 1, 1))
`MEROPS-MP_best_hit` <- `MEROPS-MP_best_hit` %>% mutate(MEROPS_SUBTYPE = str_extract(MEROPS_ID,".{1,}.+?(?=\\.)"))
`MEROPS_TYPE_plot` <- `MEROPS-MP_best_hit` %>% group_by(MEROPS_TYPE) %>% summarise(count=n()) %>% as.data.frame()
`MEROPS_SUBTYPE_plot` <- `MEROPS-MP_best_hit` %>% group_by(MEROPS_SUBTYPE) %>% summarise(count=n()) %>% as.data.frame()
`MEROPS_SUBTYPE_plot` <- `MEROPS_SUBTYPE_plot` %>% mutate(MEROPS_TYPE = substr(MEROPS_SUBTYPE, 1, 1))

theme_bw <- function(base_size = 12) {
  structure(list(
    axis.line =         theme_blank(),
    axis.text.x =       theme_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 1),
    axis.text.y =       theme_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
    axis.ticks =        theme_segment(colour = "black", size = 0.2),
    axis.title.x =      theme_text(size = base_size, vjust = 1),
    axis.title.y =      theme_text(size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.3, "lines"),
    axis.ticks.margin = unit(0.5, "lines"),
    
    legend.background = theme_rect(colour=NA), 
    legend.key =        theme_rect(colour = "grey80"),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       theme_text(size = base_size * 0.8),
    legend.title =      theme_text(size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",
    
    panel.background =  theme_rect(fill = "white", colour = NA), 
    panel.border =      theme_rect(fill = NA, colour="grey50"), 
    panel.grid.major =  theme_line(colour = "grey90", size = 0.2),
    panel.grid.minor =  theme_line(colour = "grey98", size = 0.5),
    panel.margin =      unit(0.25, "lines"),
    
    strip.background =  theme_rect(fill = "grey80", colour = "grey50"), 
    strip.text.x =      theme_text(size = base_size * 0.8),
    strip.text.y =      theme_text(size = base_size * 0.8, angle = -90),
    
    plot.background =   theme_rect(colour = NA),
    plot.title =        theme_text(size = base_size * 1.2),
    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
  ), class = "options")
}


#histogram
p1 <- ggplot(`MEROPS_SUBTYPE_plot`, aes(x = MEROPS_SUBTYPE, y = count,fill = MEROPS_TYPE))+    # x??: ž???? ??, ž???? ???? ???? ?????? ?׷��??ϴ?.
  #ggtitle("Protease") +
  geom_bar(stat="identity", color="black", width = 0.8)+
  #geom_text(aes(label=`count`), vjust=1.6, color="white")+
  scale_fill_manual(values = c(A="#ABCE30", #Aspartic
                               C="#003215", #Cysteine
                               G="#E3E3E3", #Glutamic
                               I="#0E8BC4", #Inhibitor
                               M="#58B4EE", #Metallo
                               N="#015352", #Asparagine
                               P="#FA1A0D", #Mixed
                               S="#fd9d00", #Serine
                               T="#f8d42e", #Threonine
                               U="#f8d45e"  #Unknown
                               )
                    )+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color='grey'),
        panel.grid.minor.y = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
         axis.text.y = element_text(size = 16), 
         axis.title.x = element_blank(),
         axis.title.y = element_text(size=14, face="bold"))+
  theme(panel.background = element_rect(fill = NA, 
                                        colour = "black",
                                        linetype = 'solid',
                                        size = 1))+
  theme(legend.position = "bottom", 
        #legend.title = element_blank()
        ) +
  #guides(fill=guide_legend(nrow=1))+
  labs(y = "Count", x= "Type", fill= "Type")+
  scale_x_discrete(expand = c(0.02, 0.02))+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, ceiling(max(`MEROPS_SUBTYPE_plot`$count)/10)*10), 
                     breaks=seq(0,ceiling(max(`MEROPS_SUBTYPE_plot`$count)/10)*10,5))

#`MEROPS_SUBTYPE_plot`$MEROPS_SUBTYPE <- factor(`MEROPS_SUBTYPE_plot`$MEROPS_SUBTYPE,                                   # Factor levels in decreasing order
#                                               levels = `MEROPS_SUBTYPE_plot`$MEROPS_SUBTYPE[order(`MEROPS_SUBTYPE_plot`$count, decreasing = TRUE)])

# Stacked bar plot
# Look factors
levels(as.factor(`MEROPS_TYPE_plot`$MEROPS_TYPE))
# change the sequence

`MEROPS_TYPE_plot`$MEROPS_TYPE <- factor(`MEROPS_TYPE_plot`$MEROPS_TYPE, 
                                         levels=`MEROPS_TYPE_plot` %>% arrange(desc(count)) %>% dplyr::select(MEROPS_TYPE) %>% pull())

`MEROPS_TYPE_plot` <- `MEROPS_TYPE_plot` %>%  mutate(TYPE = "MEROPS")
`MEROPS_TYPE_plot` <- `MEROPS_TYPE_plot` %>%  mutate(Proportion = round(count/sum(.$count)*100,1))
p2<- ggplot(`MEROPS_TYPE_plot`, aes(x=TYPE, y = Proportion, fill = MEROPS_TYPE)) +
  geom_col() +
  geom_text(aes(label = paste0(Proportion, "%")),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(A="#ABCE30",
                               C="#003215",
                               G="#E3E3E3",
                               I="#0E8BC4",
                               M="#58B4EE",
                               N="#015352",
                               P="#FA1A0D",
                               S="#fd9d00",
                               T="#f8d42e",
                               U="#f8d45e")) +
  theme(panel.background = element_rect(fill = NA, 
                                        colour = "black",
                                        linetype = 'solid',
                                        size = 1))+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color='grey'),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, face="bold"))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0), position = "right")+ 
  theme(legend.position = "none", 
        legend.title = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(r = 20))) +
  ylab(NULL) +
  xlab("Proportion")


# plot combine
prow <- plot_grid( p1 + theme(legend.position="none"),
                   p2 + theme(legend.position="none"),
                   rel_widths = c(8.5, 1.5),
                   align = 'vh',
                   labels = c("A", "B"),
                   hjust = -1,
                   nrow = 1
)
legend_common <- get_legend(p1 + theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)))

p <- plot_grid(prow, legend_common, ncol = 1, rel_heights = c(1, 0.2)) 

p

ggsave("MEROPS_bar_plot.pdf", dpi = 300) 


##############################################################################################################################
pacman::p_load(ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr, seqinr, ggtree, magrittr, sjmisc, ggseqlogo, genoPlotR, gridExtra, seqinr, ape, rentrez, reutils, stringr, reshape2, xlsx, openxlsx, gtools, Biostrings, data.table, rhmmer, RFLPtools)

################################################################################################################################
                                              # Load MEROPS_TYPE_corrplot or dotplot #
################################################################################################################################
read.table("MEROPS_PepList.xlsx",header=TRUE, quote = "",na.strings = "") # only at first time

{
  working_dir<-getwd()
  #go to the database directory
  #setwd("../../DataBase")
  #file_names=as.list(dir(pattern="*\\.RData$")) 
  #Database_list <- list.files(src_dir, pattern = "\\.RData$") # list
  Database_list <- list.files(getwd(), pattern = "MEROPS-MP_best_hit_before_filter\\.xlsx$") # list
  df.list <- lapply(Database_list, read.xlsx)
  names(df.list) <- Database_list %>% gsub("_MEROPS-MP_best_hit_before_filter.xlsx","",.)
  
  for(i in 1:length(df.list)){
    assign(names(df.list[i]),df.list[i] %>% as.data.frame() %>% dplyr::select(contains("MEROPS_ID")) %>% setNames(.,"MEROPS_ID") %>% 
             mutate(MEROPS_TYPE = substr(MEROPS_ID, 1, 1)) %>% mutate(MEROPS_SUBTYPE = str_extract(MEROPS_ID,".{1,}.+?(?=\\.)")) %>% 
             group_by(MEROPS_TYPE) %>% summarise(count = n()) %>% as.data.frame() %>% setnames(., old = "count", new = names(df.list[i]))
    )
  }
  # merge all the datas
  df.MEROPS_TYPE <- mget(ls(pattern = names(df.list) %>% paste0(collapse = "|")))
  df.MEROPS_TYPE <- Reduce(full_join, df.MEROPS_TYPE)
  df.MEROPS_TYPE<- df.MEROPS_TYPE %>%
    arrange(match(MEROPS_TYPE,  mixedsort(df.MEROPS_TYPE$MEROPS_TYPE)))
  #list(df.MEROPS_TYPE) %>% reduce(full_join, by = "MEROPS_TYPE")
  
  # data processing to plot MEROPS (SUB)TYPE
  for(i in 1:length(df.list)){
    assign(names(df.list[i]),df.list[i] %>% as.data.frame() %>% dplyr::select(contains("MEROPS_ID")) %>% setNames(.,"MEROPS_ID") %>% 
             mutate(MEROPS_SUBTYPE = substr(MEROPS_ID, 1, 1)) %>% mutate(MEROPS_SUBTYPE = str_extract(MEROPS_ID,".{1,}.+?(?=\\.)")) %>% 
             group_by(MEROPS_SUBTYPE) %>% summarise(count = n()) %>% as.data.frame() %>% setnames(., old = "count", new = names(df.list[i]))
    )
  }
  df.MEROPS_SUBTYPE <- mget(ls(pattern = names(df.list) %>% paste0(collapse = "|")))
  df.MEROPS_SUBTYPE <- Reduce(full_join, df.MEROPS_SUBTYPE)
  df.MEROPS_SUBTYPE<- df.MEROPS_SUBTYPE %>%
    arrange(match(MEROPS_SUBTYPE,  mixedsort(df.MEROPS_SUBTYPE$MEROPS_SUBTYPE)))
  
  list2env(df.list,.GlobalEnv) # set them into the global environment
  
  rm(list=c("Database_list","src_dir","sFam","MEROPS-MPEP","MEROPS-MPRO",ls(pattern="fasta$")))
  #back to the working directory
  setwd(working_dir)
}



# make data square to calculate euclidean distance
df.MEROPS_TYPE_melt = melt(df.MEROPS_TYPE)

mat <- df.MEROPS_TYPE_melt %>%  
  #select(-cell_ct, -cell_exp_ct, -Group) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = MEROPS_TYPE, values_from = value) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$variable  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
v_clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
#h_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
############ NOTICE THE t() above)

ddgram_col <- as.dendrogram(v_clust)
#ddgram_col <- as.dendrogram(h_clust)
ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()
ggtree_plot_col



# create plot using ggplot function. geom_point() is used for creating dot plots

dotplot <- df.MEROPS_TYPE_melt %>% 
  #mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
  #       Gene = factor(Gene, levels = gene_pos_table$gene)) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=factor(variable, levels = v_clust$labels[v_clust$order]), 
             y=MEROPS_TYPE,
             color = value, 
             size = value
             )) + 
  scale_x_discrete(labels=function(x) gsub("_", " ",x))+
  scale_size(range = c(0, 10))+

  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.0, hjust=1.0),
        axis.title.x = element_text(size = 16, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(panel.background = element_rect(fill = NA, 
                                        colour = "black",
                                        linetype = 'solid',
                                        size = 1))+
  
  labs(x="Strains",
       y="Merops Type") +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,80), oob = scales::squish, name = 'Count')

dotplot

ggtree_plot <- ggtree_plot_col + xlim2(dotplot)

#ggtree_plot


labels <- ggplot(df.MEROPS_TYPE_melt %>% 
                   mutate(`Genus` = qdap::beg2char(df.MEROPS_TYPE_melt$variable, "_") %>% as.factor(),
                          `Species` = qdapRegex::ex_between(.$variable, "_", "_") %>% map(.,1) %>% unlist() %>% as.factor(),
                          MEROPS_TYPE = factor(MEROPS_TYPE, levels = v_clust$labels[v_clust$order]),
                          variable = factor(variable, levels = v_clust$labels[v_clust$order])), 
                 aes(x = df.MEROPS_TYPE_melt$variable, y = 1, fill = Species)) + 
  geom_tile() + 
  #scale_fill_brewer(palette = 'Paired') + 
  scale_fill_manual(values =c("stearothermophilus"="#CA0020",
                              "zalihae"="#F4A582",
                              "caldolyticus"="#4DAC26",
                              "thermocatenulatus"="#B8E186",
                              "subterraneus"="#0571B0",
                              "thermodenitrificans"="#92C5DE",
                              "thermoleovorans"="#B2ABD2",
                              "kaustophilus"="#5E3C99",
                              "caldoxylosilyticus"="#80CDC1",
                              "thermoglucosidasius"="#018571",
                              "toebii"="#DFC27D",
                              "sp."="#BABABA",
                              "genomosp."="#404040"))+
  theme_nothing() +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  aplot::xlim2(dotplot)

labels_genus <- ggplot(df.MEROPS_TYPE_melt %>% 
                   mutate(`Genus` = qdap::beg2char(df.MEROPS_TYPE_melt$variable, "_") %>% as.factor(),
                          `Species` = qdapRegex::ex_between(.$variable, "_", "_") %>% map(.,1) %>% unlist() %>% as.factor(),
                          MEROPS_TYPE = factor(MEROPS_TYPE, levels = v_clust$labels[v_clust$order]),
                          variable = factor(variable, levels = v_clust$labels[v_clust$order])), 
                 aes(x = df.MEROPS_TYPE_melt$variable, y = 1, fill = Genus)) + 
  geom_tile() + 
  #scale_fill_brewer(palette = 'Spectral') + 
  scale_fill_manual(values=c(`Geobacillus`="#E66101",
                             `Parageobacillus`="#FDB863"))+
  theme_nothing() +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  aplot::xlim2(dotplot)

legend <- plot_grid(get_legend(labels_genus + theme(legend.position=c("0.3","0.5"))),
                    get_legend(labels + theme(legend.position=c("0.1","0.5")))
)
  ggtree_plot + 
  labels_genus +
  labels + 
  dotplot + 
  legend + 
  plot_layout(ncol = 1, heights = c(0.5, 0.25, 0.25, 
                                    4.0,
                                    1.0)
  )


ggsave("Merops_type_by_strains_dot_plot.pdf", 
       #dpi = 300, 
       device = "pdf", 
       width = 500,
       height = 250,
       units = "mm",
       limitsize = FALSE,
       #scale=1.5
       )
################################################################################################################################



################################################################################################################################
                                  # Load MEROPS_SUBTYPE_corrplot or dotplot #
################################################################################################################################

# make data square to calculate euclidean distance
df.MEROPS_SUBTYPE_melt = melt(df.MEROPS_SUBTYPE)

mat <- df.MEROPS_SUBTYPE_melt %>%  
  #select(-cell_ct, -cell_exp_ct, -Group) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = MEROPS_SUBTYPE, values_from = value) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$variable  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
v_clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
#h_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
############ NOTICE THE t() above)

ddgram_col <- as.dendrogram(v_clust)
#ddgram_col <- as.dendrogram(h_clust)
ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()
ggtree_plot_col



# create plot using ggplot function. geom_point() is used for creating dot plots

dotplot <- df.MEROPS_SUBTYPE_melt %>% 
  #mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
  #       Gene = factor(Gene, levels = gene_pos_table$gene)) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=factor(variable, levels = v_clust$labels[v_clust$order]), 
             y=MEROPS_SUBTYPE,
             color = value, 
             size = value
  )) + 
  scale_x_discrete(labels=function(x) gsub("_", " ",x))+
  scale_size(range = c(0, 10))+
  
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.0, hjust=1.0),
        axis.title.x = element_text(size = 16, face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(panel.background = element_rect(fill = NA, 
                                        colour = "black",
                                        linetype = 'solid',
                                        size = 1))+
  
  labs(x="Strains",
       y="Merops Type") +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,16), oob = scales::squish, name = 'Count')

dotplot

ggtree_plot <- ggtree_plot_col + xlim2(dotplot)

#ggtree_plot


labels <- ggplot(df.MEROPS_SUBTYPE_melt %>% 
                   mutate(`Genus` = qdap::beg2char(df.MEROPS_SUBTYPE_melt$variable, "_") %>% as.factor(),
                          `Species` = qdapRegex::ex_between(.$variable, "_", "_") %>% map(.,1) %>% unlist() %>% as.factor(),
                          MEROPS_SUBTYPE = factor(MEROPS_SUBTYPE, levels = v_clust$labels[v_clust$order]),
                          variable = factor(variable, levels = v_clust$labels[v_clust$order])), 
                 aes(x = df.MEROPS_SUBTYPE_melt$variable, y = 1, fill = Species)) + 
  geom_tile() + 
  #scale_fill_brewer(palette = 'Paired') + 
  scale_fill_manual(values =c("stearothermophilus"="#CA0020",
                              "zalihae"="#F4A582",
                              "caldolyticus"="#4DAC26",
                              "thermocatenulatus"="#B8E186",
                              "subterraneus"="#0571B0",
                              "thermodenitrificans"="#92C5DE",
                              "thermoleovorans"="#B2ABD2",
                              "kaustophilus"="#5E3C99",
                              "caldoxylosilyticus"="#80CDC1",
                              "thermoglucosidasius"="#018571",
                              "toebii"="#DFC27D",
                              "sp."="#BABABA",
                              "genomosp."="#404040"))+
  theme_nothing() +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  aplot::xlim2(dotplot)

labels_genus <- ggplot(df.MEROPS_SUBTYPE_melt %>% 
                         mutate(`Genus` = qdap::beg2char(df.MEROPS_SUBTYPE_melt$variable, "_") %>% as.factor(),
                                `Species` = qdapRegex::ex_between(.$variable, "_", "_") %>% map(.,1) %>% unlist() %>% as.factor(),
                                MEROPS_SUBTYPE = factor(MEROPS_SUBTYPE, levels = v_clust$labels[v_clust$order]),
                                variable = factor(variable, levels = v_clust$labels[v_clust$order])), 
                       aes(x = df.MEROPS_SUBTYPE_melt$variable, y = 1, fill = Genus)) + 
  geom_tile() + 
  #scale_fill_brewer(palette = 'Spectral') + 
  scale_fill_manual(values=c(`Geobacillus`="#E66101",
                             `Parageobacillus`="#FDB863"))+
  theme_nothing() +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  aplot::xlim2(dotplot)

legend <- plot_grid(get_legend(labels_genus + theme(legend.position=c("0.3","0.5"))),
                    get_legend(labels + theme(legend.position=c("0.1","0.5")))
)
ggtree_plot + 
  labels_genus +
  labels + 
  dotplot + 
  legend + 
  plot_layout(ncol = 1, heights = c(0.5, 0.25, 0.25, 
                                    12.0,
                                    1.0)
  )


ggsave("MEROPS_SUBTYPE_by_strains_dot_plot.pdf", 
       #dpi = 300, 
       device = "pdf", 
       width = 500,
       height = 700,
       units = "mm",
       limitsize = FALSE,
       #scale=1.5
)




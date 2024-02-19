pacman::p_load(Peptides, ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr, seqinr, ggtree, magrittr, sjmisc, ggseqlogo, genoPlotR, gridExtra, seqinr, ape, rentrez, reutils, stringr, reshape2, xlsx, openxlsx, gtools, Biostrings, data.table, varhandle)

todor::todor_file("20230607_RBS_region_extraction.R")

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

read.table("\\.fna")


RBS <- contig_final %>% dplyr::filter(product=="16S ribosomal RNA") %>% dplyr::select(rearranged_nt_seq) %>% mutate(RBS=str_sub(.$rearranged_nt_seq, start = -9))
RBS$RBS <- sapply(RBS$RBS,function(x) as.character(Biostrings::reverseComplement(DNAString(x))))

# RBS Sequence library with Start codon
max(RBS$RBS) #sort(table(RBS$RBS), decreasing = TRUE)

# TODO: change nucleotide to contig_seq_final$top_strand and contig_seq_final$bottom_strand
#for(i in nrow(contig_seq_final)){}

strand_plus_RBS <- vmatchPattern(toupper(max(RBS$RBS)), contig_seq_final$top_strand %>% toupper(),
                                 max.mismatch=2, min.mismatch=0,
                                 with.indels=FALSE, fixed=TRUE,
                                 algorithm="auto") %>% as.data.frame() %>% #dplyr::select(c(start,end,width)) %>% 
                                 mutate(strand="+")
temp_df <- data.frame()
for(i in 1:nrow(contig_seq_final)){
  temp_RBS_df <- strand_plus_RBS %>% dplyr::filter(group==i) %>% mutate(RBS= str_sub(contig_seq_final$top_strand[i], start, end),
                                                                        TR= str_sub(contig_seq_final$top_strand[i], start, end+100))
  temp_df <- rbind.fill(temp_df, temp_RBS_df)
}             
strand_plus_RBS <- temp_df
rm(temp_df)

range <- "{1,10}" # set the range
strand_plus_RBS <- strand_plus_RBS %>% mutate(RBS_hit= stringr::str_extract(strand_plus_RBS$TR, paste0(strand_plus_RBS$RBS,".",range,"ATG")))
strand_plus_RBS <- strand_plus_RBS %>% mutate(spacer= nchar(RBS_hit)-nchar("ATG")-nchar(RBS))
strand_plus_RBS <- strand_plus_RBS %>% mutate(ORFstart= start + nchar(RBS) + spacer)
strand_plus_RBS <- strand_plus_RBS %>% mutate(`ORFmatch`=case_when(ORFstart %in% contig_final$start ~ "match",
                                                                   TRUE~ "unmatch"))

################################################################################################################################################
# TODO
strand_minus_RBS <- vmatchPattern(toupper(max(RBS$RBS)), contig_seq_final$bottom_strand %>% toupper(),
                                 max.mismatch=2, min.mismatch=0,
                                 with.indels=FALSE, fixed=TRUE,
                                 algorithm="auto") %>% as.data.frame() %>% #dplyr::select(c(start,end,width)) %>% 
  mutate(strand="-")
temp_df <- data.frame()
for(i in 1:nrow(contig_seq_final)){
  temp_RBS_df <- strand_minus_RBS %>% dplyr::filter(group==i) %>% mutate(RBS= str_sub(contig_seq_final$bottom_strand[i], start, end),
                                                                        TR= str_sub(contig_seq_final$bottom_strand[i], start, end+100))
  temp_df <- rbind.fill(temp_df, temp_RBS_df)
}             
strand_minus_RBS <- temp_df
rm(temp_df)

range <- "{1,10}" # set the range
strand_minus_RBS <- strand_minus_RBS %>% mutate(RBS_hit= stringr::str_extract(strand_minus_RBS$TR, paste0(strand_minus_RBS$RBS,".",range,"ATG")))
strand_minus_RBS <- strand_minus_RBS %>% mutate(spacer= nchar(RBS_hit)-nchar("ATG")-nchar(RBS))

strand_minus_RBS <- strand_minus_RBS %>% mutate(ORFstart= start + nchar(RBS) + spacer)

temp_df <- data.frame()
for(i in 1:nrow(contig_seq_final)){
  temp_RBS_df <- strand_minus_RBS %>% dplyr::filter(group==i) %>% mutate(ORFstart= nchar(contig_seq_final$bottom_strand[i])-as.numeric(start + nchar(RBS) + spacer)+1)
  temp_df <- rbind.fill(temp_df, temp_RBS_df)
}             
strand_minus_RBS <- temp_df
rm(temp_df, temp_RBS_df)

strand_minus_RBS <- strand_minus_RBS %>% mutate(`ORFmatch`=case_when(ORFstart %in% (contig_final %>% dplyr::filter(direction=="-") %>% dplyr::select(end) %>% pull) ~ "match",
                                                                   TRUE~ "unmatch"))

# filtering
putative_RBS <- rbind.fill(strand_plus_RBS, strand_minus_RBS)
putative_matched_RBS<- putative_RBS %>% filter(ORFmatch=="match")


# plotting
# seqlogo
# Create custom colour scheme
cs1 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), groups=c('A', 'T', 'G', 'C'), 
                      cols=c('#76B56A', '#D06461', '#3B81A1', '#EFCA70'))

# Generate sequence logo
ggseqlogo(putative_matched_RBS$RBS, col_scheme=cs1)

plot <- ggplot(putative_matched_RBS, aes(x = spacer)) 

plot <- plot + geom_histogram(aes(y=..density..), color="black", fill = "#EBE5DE", binwidth = 1, alpha = 1.0) +
  scale_x_continuous(breaks = seq(1, 10, by=1))
#plot <- plot + geom_density()
plot+theme_bw()+
  theme(panel.background = element_blank(),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(angle = 90, size=30),
        axis.text.x = element_text(color = "grey20", size = 30),
        axis.text.y = element_text(color = "grey20", size = 30),  
        plot.title = element_text(size=30)
  )+
  xlab("Spacer") + 
  ylab("Density") +
  ggtitle('RBS to AUG')



pacman::p_load(stringi, rlist, tidyverse, purrr, readr, Peptides, ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr, seqinr, ggtree, magrittr, sjmisc, ggseqlogo, genoPlotR, gridExtra, seqinr, ape, rentrez, reutils, stringr, reshape2, xlsx, openxlsx, gtools, Biostrings, data.table, varhandle)

# attach promotech promoter prediction result to genbank for the snapgene visualization

# processing my promotech_output
xlsx_dir <- getwd()
xlsx_file <- list.files(xlsx_dir, pattern = "genome_predictions\\.csv$") # list
genome_predictions <- read.csv(xlsx_file, sep="\t")

genome_predictions <- genome_predictions %>% mutate(position= case_when(strand=="-" ~ paste0("complement","(",start,"..",end,")"),
                                                                        strand=="+" ~ paste0(start,"..",end)
                                                                        )
                                                    )
# for snapgene format
#primer_bind     complement(3617..3634)
#/note="M13 rev"
#/note="common sequencing primer, one of multiple similar variants"
#/note="color: #a020f0; direction: LEFT"
#misc_feature    2591..2806
#/note="pRplsWT"
#/note="color: #ffcc99; direction: RIGHT"

# for primers
#primer_bind     2611..2627
#/label=James
#/note="'hello'"
#/note="color: black; sequence: gtttttgcgccgcccgg"

#primer_bind     complement(2611..2627)
#/label=James
#/note="'hello'"
#/note="color: black; sequence: gtttttgcgccgcccgg"


Feature<- paste0(paste0(rep(" ",5), collapse=""), 
                 "Promoter",
                 paste0(rep(" ", (21-5-nchar("Promoter"))) ,collapse=""))

paste0(rep(" ",21),collapse="")

Feature<- paste0(Feature, genome_predictions$position,
           "\n",
           paste0(rep(" ",21),collapse=""), paste0("/gene=",genome_predictions$score))
write.table(Feature, "promoter_feature_for_gb",col.names = FALSE,row.names = FALSE, quote = FALSE)

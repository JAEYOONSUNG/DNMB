pacman::p_load(qdap, qdapRegex, ComplexHeatmap, ggdendro, patchwork, forcats, data.table, purrr, cowplot, ggpubr, gtable, RFLPtools,ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, rhmmer, stringr, seqinr, ggtree, magrittr, sjmisc, ggseqlogo, genoPlotR, gridExtra, seqinr, ape, rentrez, reutils, stringr, reshape2, xlsx, openxlsx, gtools, Biostrings,varhandle)

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
{
  working_dir<-getwd()
  #go to the database directory
  #setwd("../../DataBase")
  #file_names=as.list(dir(pattern="*\\.RData$")) 
  #Database_list <- list.files(src_dir, pattern = "\\.RData$") # list
  Database_list <- list.files(getwd(), pattern = "_dbCAN_best_hit\\.xlsx$") # list
  df.list <- lapply(Database_list, read.xlsx)
  names(df.list) <- Database_list %>% gsub("_dbCAN_best_hit.xlsx","",.)
  
  list2env(df.list,.GlobalEnv) # set them into the global environment
  
  #rm(list=c(,ls(pattern="")))
  #back to the working directory
  setwd(working_dir)
}

##################################################################################################################
# data processing to plot CAZy Class
for(i in 1:length(df.list)){
  assign(names(df.list[i]),df.list[i] %>% as.data.frame() %>% dplyr::select(contains("CAZy.Class")) %>% setNames(.,"CAZy Class") %>% 
           group_by(`CAZy Class`) %>% summarise(count = n()) %>% as.data.frame() %>% setnames(., old = "count", new = names(df.list[i]))
  )
}
df.CAZy_Class <- mget(ls(pattern = names(df.list) %>% paste0(collapse = "|")))
df.CAZy_Class <- Reduce(full_join, df.CAZy_Class)
df.CAZy_Class<- df.CAZy_Class %>%
  arrange(match(`CAZy Class`,  mixedsort(df.CAZy_Class$`CAZy Class`)))


# make data square to calculate euclidean distance
df.CAZy_Class_melt = melt(df.CAZy_Class)

mat <- df.CAZy_Class_melt %>%  
  #select(-cell_ct, -cell_exp_ct, -Group) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = `CAZy Class`, values_from = value) %>% 
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

dotplot <- df.CAZy_Class_melt %>% 
  #mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
  #       Gene = factor(Gene, levels = gene_pos_table$gene)) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=factor(variable, levels = v_clust$labels[v_clust$order]), 
             y=`CAZy Class`,
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
       y="CAZy Class") +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::inferno(75), limits = c(0,75), oob = scales::squish, name = 'Count')
  #max(df.CAZy_Class_melt$value, na.rm=TRUE)
dotplot

ggtree_plot <- ggtree_plot_col + aplot::xlim2(dotplot)

#ggtree_plot


labels <- ggplot(df.CAZy_Class_melt %>% 
                   mutate(`Genus` = qdap::beg2char(df.CAZy_Class_melt$variable, "_") %>% as.factor(),
                          `Species` = qdapRegex::ex_between(.$variable, "_", "_") %>% map(.,1) %>% unlist() %>% as.factor(),
                          `CAZy Class` = factor(`CAZy Class`, levels = v_clust$labels[v_clust$order]),
                          variable = factor(variable, levels = v_clust$labels[v_clust$order])), 
                 aes(x = df.CAZy_Class_melt$variable, y = 1, fill = Species)) + 
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


labels_genus <- ggplot(df.CAZy_Class_melt %>% 
                         mutate(`Genus` = qdap::beg2char(df.CAZy_Class_melt$variable, "_") %>% as.factor(),
                                `Species` = qdapRegex::ex_between(.$variable, "_", "_") %>% map(.,1) %>% unlist() %>% as.factor(),
                                `CAZy Class` = factor(`CAZy Class`, levels = v_clust$labels[v_clust$order]),
                                variable = factor(variable, levels = v_clust$labels[v_clust$order])), 
                       aes(x = df.CAZy_Class_melt$variable, y = 1, fill = Genus)) + 
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


ggsave("CAZy_Class_by_strains_dot_plot.pdf", 
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
# Load CAZy_Family_corrplot or dotplot #
################################################################################################################################
# data processing to plot CAZy Family
for(i in 1:length(df.list)){
  assign(names(df.list[i]),df.list[i] %>% as.data.frame() %>% dplyr::select(contains("CAZy.Family")) %>% setNames(.,"CAZy Family") %>% 
           group_by(`CAZy Family`) %>% summarise(count = n()) %>% as.data.frame() %>% setnames(., old = "count", new = names(df.list[i]))
  )
}
df.CAZy_Family <- mget(ls(pattern = names(df.list) %>% paste0(collapse = "|")))
df.CAZy_Family <- Reduce(full_join, df.CAZy_Family)
df.CAZy_Family<- df.CAZy_Family %>%
  arrange(match(`CAZy Family`,  mixedsort(df.CAZy_Family$`CAZy Family`)))


# make data square to calculate euclidean distance
df.CAZy_Family_melt = melt(df.CAZy_Family)

mat <- df.CAZy_Family_melt %>%  
  #select(-cell_ct, -cell_exp_ct, -Group) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = 'CAZy Family', values_from = value) %>% 
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

dotplot <- df.CAZy_Family_melt %>% 
  #mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
  #       Gene = factor(Gene, levels = gene_pos_table$gene)) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=factor(variable, levels = v_clust$labels[v_clust$order]), 
             y=`CAZy Family`,
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
       y="CAZy Family") +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::inferno(25), limits = c(0,25), oob = scales::squish, name = 'Count')

dotplot

ggtree_plot <- ggtree_plot_col + aplot::xlim2(dotplot)

#ggtree_plot


labels <- ggplot(df.CAZy_Family_melt %>% 
                   mutate(`Genus` = qdap::beg2char(df.CAZy_Family_melt$variable, "_") %>% as.factor(),
                          `Species` = qdapRegex::ex_between(.$variable, "_", "_") %>% map(.,1) %>% unlist() %>% as.factor(),
                          `CAZy Family` = factor(`CAZy Family`, levels = v_clust$labels[v_clust$order]),
                          variable = factor(variable, levels = v_clust$labels[v_clust$order])), 
                 aes(x = df.CAZy_Family_melt$variable, y = 1, fill = Species)) + 
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

labels_genus <- ggplot(df.CAZy_Family_melt %>% 
                         mutate(`Genus` = qdap::beg2char(df.CAZy_Family_melt$variable, "_") %>% as.factor(),
                                `Species` = qdapRegex::ex_between(.$variable, "_", "_") %>% map(.,1) %>% unlist() %>% as.factor(),
                                `CAZy Family` = factor(`CAZy Family`, levels = v_clust$labels[v_clust$order]),
                                variable = factor(variable, levels = v_clust$labels[v_clust$order])), 
                       aes(x = df.CAZy_Family_melt$variable, y = 1, fill = Genus)) + 
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


ggsave("CAZy_Family_by_strains_dot_plot.pdf", 
       #dpi = 300, 
       device = "pdf", 
       width = 500,
       height = 700,
       units = "mm",
       limitsize = FALSE,
       #scale=1.5
)





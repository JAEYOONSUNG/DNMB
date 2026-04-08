suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(tidyr)
  library(cowplot); library(circlize); library(ComplexHeatmap)
})
pkg_dir <- "/Users/JaeYoon/Dropbox/0.Personal folder/5. Bioinformatics/DNMB"
for (f in c("R/db_modules.R","R/module_api.R","R/Mobileome_pipeline.R",
            "R/Mobileome_sequence_engine.R","R/Mobileome_variant_engine.R",
            "R/Mobileome_comparative.R","R/Mobileome_landing_pads.R",
            "R/Mobileome_auto_comparative.R","R/db_module_iselement.R"))
  source(file.path(pkg_dir, f), local = FALSE)

parsed <- .dnmb_parse_genbank_features("GCF_030369615.1.gbff")
genes <- .dnmb_predict_gene_essentiality(.dnmb_build_gene_table(parsed$features))
result <- readRDS("iselement_full_result.rds")
elements <- result$elements; census <- result$census; target_models <- result$target_models
landing_pads <- read.delim("dnmb_module_iselement/iselement_landing_pads.tsv", check.names=FALSE)
genome_len <- max(genes$end, na.rm=TRUE)
ctg <- "NZ_CP128494"
plot_dir <- file.path(getwd(), "visualizations")

ess_track <- genes %>% filter(contig == ctg)
is_track <- elements %>% filter(contig == ctg) %>%
  mutate(fam = ifelse(is.na(element_family)|!nzchar(element_family), "unknown", element_family))
lp_track <- landing_pads %>% filter(contig == ctg)
fam_cols <- c(IS110="#7E57C2",IS982="#26A69A",IS701="#EF5350",IS630="#42A5F5",
              IS66="#FFA726",IS4="#AB47BC",IS21="#26C6DA",IS3="#EC407A",
              IS481="#5C6BC0",ISL3="#66BB6A","IS200/IS605"="#8D6E63",
              IS256="#78909C",IS5="#D4E157",IS1595="#FF7043",unknown="#9E9E9E")

# ======================================================================
# Page 1: Circlize genome map with IS links (publication quality)
# ======================================================================
pdf(file.path(plot_dir, "ISelement_overview.pdf"), width=14, height=14)

par(mar=c(1,1,2,1))
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 2,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.005, 0.005),
  points.overflow.warning = FALSE
)

# Initialize
genome_bed <- data.frame(chr=ctg, start=0, end=genome_len)
circos.genomicInitialize(genome_bed, plotType=NULL, major.by=200000, axis.labels.cex=0.55)

# Title
title(paste0("IS Element Distribution & Landing Pad Map\n",
             ctg, " (", round(genome_len/1e6,2), " Mb) | ",
             nrow(is_track), " IS elements across ", length(unique(is_track$fam)), " families"),
      cex.main=1.2, font.main=2, line=0.5)

# --- Track 0: Genomic axis with kb labels ---
circos.track(ylim=c(0,1), track.height=0.04, bg.border=NA, panel.fun=function(x,y) {
  circos.genomicAxis(h="top", major.by=200000, 
                     labels.cex=0.45, labels.facing="clockwise",
                     )
})

# --- Track 1: Gene Essentiality (genomicTrack with BED) ---
ess_bed <- data.frame(
  chr = ctg,
  start = ess_track$start,
  end = ess_track$end,
  value = ess_track$essentiality_score
)
ess_col_fun <- colorRamp2(c(0, 0.2, 0.45, 0.7, 1.0),
                           c("#43A047","#A5D6A7","#FFFDE7","#FFAB91","#D32F2F"))

circos.genomicTrack(ess_bed, track.height=0.06, bg.border="grey50", bg.lwd=0.5,
                    panel.fun=function(region, value, ...) {
  circos.genomicRect(region, value, ytop.column=1, ybottom=0,
                     col=ess_col_fun(value[[1]]), border=NA)
})

# --- Track 2: IS Elements (genomicTrack, family-colored) ---
is_bed <- data.frame(
  chr = ctg,
  start = is_track$start,
  end = is_track$end,
  family = is_track$fam
)
circos.genomicTrack(is_bed, ylim=c(0,1), track.height=0.08, bg.border="grey50", bg.lwd=0.5,
                    panel.fun=function(region, value, ...) {
  for (i in seq_len(nrow(region))) {
    col <- fam_cols[value$family[i]]
    if (is.na(col)) col <- "#9E9E9E"
    circos.rect(region$start[i], 0, region$end[i], 1, col=col, border="grey20", lwd=0.2)
  }
})

# --- Track 3: Landing Pad Score (bar/heatmap) ---
lp_bed <- data.frame(
  chr = ctg,
  start = lp_track$region_start,
  end = lp_track$region_end,
  score = lp_track$landing_pad_score
)
lp_col_fun <- colorRamp2(c(0.55, 0.65, 0.72, 0.78, 0.85),
                          c("#E53935","#FF8F00","#FDD835","#66BB6A","#1B5E20"))

circos.genomicTrack(lp_bed, track.height=0.10, bg.border="grey50", bg.lwd=0.5,
                    panel.fun=function(region, value, ...) {
  circos.genomicRect(region, value, ytop.column=1, ybottom=0,
                     col=lp_col_fun(value[[1]]), border=NA)
  # Reference line at 0.7
  circos.lines(CELL_META$cell.xlim, c(0.7, 0.7), lty=3, col="grey40", lwd=0.5)
})

# --- Track 4: IS Density (genomicDensity) ---
is_density_bed <- data.frame(chr=ctg, start=is_track$start, end=is_track$end)
circos.genomicDensity(is_density_bed, track.height=0.07, bg.border="grey50", bg.lwd=0.5,
                      col="#9C27B0",
                      window.size=30000, overlap=TRUE)

# --- Track 5: GC skew-like track (essential gene density for contrast) ---
ess_high_bed <- data.frame(
  chr=ctg,
  start=ess_track$start[ess_track$essentiality_score >= 0.45],
  end=ess_track$end[ess_track$essentiality_score >= 0.45]
)
if (nrow(ess_high_bed) > 2) {
  circos.genomicDensity(ess_high_bed, track.height=0.05, bg.border="grey50", bg.lwd=0.5,
                        col="#EF5350",
                        window.size=50000, overlap=TRUE)
}

# ============================================================
# IS family LINKS — connect same-family IS elements across genome
# ============================================================
# Filter families with >= 2 elements for linking
link_families <- is_track %>%
  count(fam, sort=TRUE) %>%
  filter(n >= 2, fam != "unknown") %>%
  pull(fam)

for (fam_name in link_families) {
  fam_elements <- is_track %>% filter(fam == fam_name) %>% arrange(start)
  if (nrow(fam_elements) < 2) next
  
  col <- fam_cols[fam_name]
  if (is.na(col)) col <- "#9E9E9E"
  link_col <- adjustcolor(col, alpha.f=0.25)
  border_col <- adjustcolor(col, alpha.f=0.5)
  
  # Connect consecutive pairs
  for (i in seq_len(nrow(fam_elements) - 1)) {
    mid1_start <- fam_elements$start[i]
    mid1_end <- fam_elements$end[i]
    mid2_start <- fam_elements$start[i+1]
    mid2_end <- fam_elements$end[i+1]
    
    # Only link if distance > 20kb (avoid cluttering nearby elements)
    if ((mid2_start - mid1_end) > 20000) {
      circos.link(ctg, c(mid1_start, mid1_end),
                  ctg, c(mid2_start, mid2_end),
                  col=link_col, border=border_col, lwd=0.5,
                  h.ratio=0.3)
    }
  }
}

# ============================================================
# Center legends (clean layout)
# ============================================================
# IS Family legend
fam_counts <- sort(table(is_track$fam), decreasing=TRUE)
legend_fams <- names(fam_counts)[seq_len(min(15, length(fam_counts)))]

text(0, 0.42, "IS Element Families", cex=0.85, font=2)

y_start <- 0.36
for (i in seq_along(legend_fams)) {
  y <- y_start - (i-1) * 0.042
  rect(-0.38, y-0.014, -0.32, y+0.014, col=fam_cols[legend_fams[i]], border="grey50", lwd=0.3)
  text(-0.30, y, paste0(legend_fams[i], " (", fam_counts[legend_fams[i]], ")"),
       cex=0.5, adj=c(0, 0.5), family="sans")
}

# Track legend
text(0.15, 0.36, "Tracks (outer to inner)", cex=0.7, font=2)
track_labels <- c("1. Gene Essentiality", "2. IS Elements (family-colored)",
                  "3. Landing Pad Score", "4. IS Element Density",
                  "5. Essential Gene Density")
for (i in seq_along(track_labels)) {
  text(0.15, 0.30 - (i-1)*0.035, track_labels[i], cex=0.5, adj=c(0,0.5))
}

# Essentiality scale
text(0.15, 0.08, "Gene Essentiality", cex=0.6, font=2)
ess_scale_cols <- ess_col_fun(seq(0, 1, length.out=20))
for (i in seq_along(ess_scale_cols)) {
  x_pos <- 0.15 + (i-1) * 0.012
  rect(x_pos, 0.03, x_pos+0.012, 0.06, col=ess_scale_cols[i], border=NA)
}
text(0.15, 0.015, "Low", cex=0.45, adj=c(0,0.5))
text(0.15 + 20*0.012, 0.015, "High", cex=0.45, adj=c(1,0.5))

# LP Score scale
text(0.15, -0.05, "Landing Pad Score", cex=0.6, font=2)
lp_scale_cols <- lp_col_fun(seq(0.55, 0.85, length.out=20))
for (i in seq_along(lp_scale_cols)) {
  x_pos <- 0.15 + (i-1) * 0.012
  rect(x_pos, -0.10, x_pos+0.012, -0.07, col=lp_scale_cols[i], border=NA)
}
text(0.15, -0.115, "0.55", cex=0.45, adj=c(0,0.5))
text(0.15 + 20*0.012, -0.115, "0.85", cex=0.45, adj=c(1,0.5))

# Link explanation
text(0, -0.20, "Lines connect same-family IS elements", cex=0.55, font=3, col="grey30")
text(0, -0.24, "(link color = IS family color)", cex=0.5, col="grey40")

circos.clear()

# ======================================================================
# Page 2: ggplot2 panels (Census + Recognition + Landing Pads)
# ======================================================================
# Census
census_plot <- census %>% arrange(desc(element_count)) %>%
  mutate(label=ifelse(!is.na(recognition_motif)&nzchar(recognition_motif),
                      paste0(family,"  [",recognition_motif,"]"),family),
         label=factor(label,levels=rev(label)))
cl <- census_plot %>%
  select(label,high_confidence_count,medium_confidence_count,low_confidence_count) %>%
  pivot_longer(-label,names_to="lev",values_to="cnt") %>%
  mutate(lev=factor(gsub("_confidence_count","",lev),levels=c("low","medium","high")))
p_census <- ggplot(cl,aes(x=label,y=cnt,fill=lev))+geom_col(width=0.7)+
  geom_text(data=census_plot,aes(x=label,y=element_count+1,label=element_count),inherit.aes=FALSE,hjust=0,size=3.5)+
  scale_fill_manual(values=c(low="#FFC107",medium="#42A5F5",high="#66BB6A"),name="Confidence")+
  coord_flip()+labs(title="IS Element Census",subtitle=paste0(sum(census$element_count)," elements, ",nrow(census)," families"),
                    x=NULL,y="Count")+
  theme_minimal(base_size=11)+
  theme(plot.title=element_text(face="bold",size=14),plot.subtitle=element_text(size=10,color="grey40"),
        axis.text.y=element_text(size=10),panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
        legend.position=c(0.82,0.18),legend.background=element_rect(fill="white",color="grey80",linewidth=0.3))

# Recognition
tm <- target_models %>% filter(n_elements >= 2) %>% arrange(desc(n_elements)) %>%
  mutate(flab=factor(paste0(family," (n=",n_elements,")"),levels=rev(paste0(family," (n=",n_elements,")"))),
         motif_d=ifelse(!is.na(model_motif),model_motif,"-"),
         tsd_d=ifelse(!is.na(dominant_tsd_len),paste0(dominant_tsd_len,"bp"),"-"))
p_recog <- ggplot(tm,aes(y=flab))+
  geom_tile(aes(x=0.5,fill=model_confidence),width=0.22,height=0.85,color="white")+
  geom_text(aes(x=0.78,label=motif_d),hjust=0,size=5,fontface="bold")+
  geom_text(aes(x=1.6,label=paste0("TSD: ",tsd_d)),hjust=0,size=4,color="grey30")+
  scale_fill_manual(values=c(high="#66BB6A",medium="#42A5F5",low="#FFC107"),name="Model Conf.",na.value="#E0E0E0")+
  scale_x_continuous(limits=c(0.3,2.5))+
  labs(title="IS Recognition Sequences",subtitle="Target site duplication patterns",y=NULL)+
  theme_void(base_size=11)+
  theme(plot.title=element_text(face="bold",size=14,hjust=0),
        plot.subtitle=element_text(size=10,color="grey40",hjust=0),
        axis.text.y=element_text(size=11,hjust=1),plot.margin=margin(10,20,10,10),
        legend.position="bottom",legend.text=element_text(size=9),legend.key.size=unit(0.4,"cm"))

# Landing pad heatmap
top_n <- min(15,nrow(lp_track))
top_lp <- lp_track %>% head(top_n) %>%
  mutate(ls=ifelse(nzchar(left_gene_product),substr(left_gene_product,1,22),""),
         rs=ifelse(nzchar(right_gene_product),substr(right_gene_product,1,22),""),
         label=paste0(landing_pad_id," (",round(region_start/1000),"kb, ",region_size_bp,"bp)\n",ls," | ",rs),
         label=factor(label,levels=rev(label)))
dl <- top_lp %>%
  select(label,comparative_conservation,essentiality_safety,is_recognition_avoidance,
         region_size_score,is_distance_score,redundancy_buffer,landing_pad_score) %>%
  pivot_longer(-label,names_to="sub",values_to="val") %>%
  mutate(sub_label=factor(case_when(
    sub=="comparative_conservation"~"Comp.\nConserv.\n(0.25)",sub=="essentiality_safety"~"Ess.\nSafety\n(0.20)",
    sub=="is_recognition_avoidance"~"IS Recog.\nAvoid.\n(0.15)",sub=="region_size_score"~"Region\nSize\n(0.15)",
    sub=="is_distance_score"~"IS\nDist.\n(0.15)",sub=="redundancy_buffer"~"Gene\nRedund.\n(0.10)",
    sub=="landing_pad_score"~"TOTAL"
  ),levels=c("Comp.\nConserv.\n(0.25)","Ess.\nSafety\n(0.20)","IS Recog.\nAvoid.\n(0.15)",
             "Region\nSize\n(0.15)","IS\nDist.\n(0.15)","Gene\nRedund.\n(0.10)","TOTAL")))

p_heatmap <- ggplot(dl,aes(x=sub_label,y=label,fill=val))+
  geom_tile(color="white",linewidth=0.5)+
  geom_text(aes(label=sprintf("%.2f",val)),size=3.2,
            fontface=ifelse(dl$sub=="landing_pad_score","bold","plain"))+
  scale_fill_gradientn(colors=c("#E53935","#FF8F00","#FDD835","#9CCC65","#2E7D32"),
                       values=scales::rescale(c(0,0.3,0.5,0.7,1)),limits=c(0,1),name="Score")+
  labs(title="Top Landing Pads - Sub-score Breakdown",
       subtitle="Ranked by composite score | Genome-wide spacing",x=NULL,y=NULL)+
  theme_minimal(base_size=11)+
  theme(plot.title=element_text(face="bold",size=14),plot.subtitle=element_text(size=10,color="grey40"),
        panel.grid=element_blank(),axis.text.x=element_text(size=9),axis.text.y=element_text(size=8),
        legend.position="right",legend.key.size=unit(0.4,"cm"))

# Assemble page 2
left_col <- plot_grid(p_census, p_recog, ncol=1, rel_heights=c(1.0,0.7))
page2 <- plot_grid(left_col, p_heatmap, ncol=2, rel_widths=c(0.42,0.58))

print(page2)
dev.off()

cat("Saved:", round(file.size(file.path(plot_dir,"ISelement_overview.pdf"))/1024,1), "KB\n")

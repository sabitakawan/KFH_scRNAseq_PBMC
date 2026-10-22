library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

####Plot for Gene-set enrichment analysis (GSEA) of GO-BP for DEgenes in monocytes from patients with KFH compared to healthy control in Supplementary Table/ Table S7 [Figure 6C]
####Analysis done in GSEA_MonoMacroDCs.R

GSEA_GO_Untreated_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                        sheet="Table S7", skip=6)

####Keeping only the relevant GO-terms for visualization
GSEA_GO_Untreated_enriched_terms<-GSEA_GO_Untreated_enriched_terms[c(-1,-4,-8,-11,-13:-17,-21:-25,-30),]
GSEA_GO_Untreated_enriched_terms_df<- data.frame(
  pathways = GSEA_GO_Untreated_enriched_terms$Description,
  NES = GSEA_GO_Untreated_enriched_terms$NES,
  qvalue = GSEA_GO_Untreated_enriched_terms$qvalues,
  Gene_Ratio = (unlist(lapply(GSEA_GO_Untreated_enriched_terms$core_enrichment, function(x) length(unlist(strsplit(x, split= "\\/" )))))/GSEA_GO_Untreated_enriched_terms$setSize)
)
GSEA_GO_Untreated_enriched_terms_df$type = "upregulated"
GSEA_GO_Untreated_enriched_terms_df$type[GSEA_GO_Untreated_enriched_terms_df$NES < 0] = "downregulated"
GSEA_GO_Untreated_enriched_terms_df$type<-factor(GSEA_GO_Untreated_enriched_terms_df$type, levels=c("upregulated","downregulated"))
GSEA_GO_Untreated_enriched_terms_df$log_qvalue <- -log10(GSEA_GO_Untreated_enriched_terms_df$qvalue)
GSEA_GO_Untreated_enriched_terms_df$Gene_Ratio <- round(GSEA_GO_Untreated_enriched_terms_df$Gene_Ratio*100)

GSEA_GO_Untreated_enriched_terms_df <- GSEA_GO_Untreated_enriched_terms_df %>% arrange(type, Gene_Ratio)
GSEA_GO_Untreated_enriched_terms_df$pathways <- factor(GSEA_GO_Untreated_enriched_terms_df$pathways, levels=GSEA_GO_Untreated_enriched_terms_df$pathways)

####PLOT
GSEA_GO_Untreated_enriched_terms_plot<- ggplot(GSEA_GO_Untreated_enriched_terms_df,
            aes(x     = Gene_Ratio,
                y     = pathways,
                size  = log_qvalue,
                color = NES)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(type ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(type = c(upregulated   = "Upregulated",
                                          downregulated = "Downregulated"))) +
  
  scale_colour_gradient2(
    midpoint = 0,
    low      = "darkslateblue",
    mid      = "white",
    high     = "darkred",
    name     = "NES",
    guide    = guide_colorbar(
      barwidth  = 0.8,
      barheight = 6,
      ticks     = TRUE)
  ) +
  
  scale_size_continuous(
    name   = expression(-log[10](italic(q))),
    range  = c(3, 9),
    breaks = pretty(GSEA_GO_Untreated_enriched_terms_df$log_qvalue, n = 4)
  ) +
  
  scale_y_discrete(labels = label_wrap(50)) +
  
  labs(x = "Gene Ratio",
       y = NULL) +
  
  theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    # Axes
    axis.text.x      = element_text(size = 9,  color = "black"),
    axis.text.y      = element_text(size = 9,  color = "black"),
    axis.title.x     = element_text(size = 10, color = "black", margin = margin(t = 6)),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    
    # No gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Facet strip styling
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.4),
    strip.text       = element_text(size = 10, color = "black", face = "bold"),
    
    # Legend
    legend.title     = element_text(size = 9,  color = "black"),
    legend.text      = element_text(size = 8,  color = "black"),
    legend.key.size  = unit(0.4, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.box       = "vertical",
    legend.position  = "right",
    legend.margin    = margin(l = 4),
    
    # Margins
    plot.margin      = margin(t = 8, r = 6, b = 8, l = 6)
  )
GSEA_GO_Untreated_enriched_terms_plot

ggsave("Figure6C.pdf",
       plot   = GSEA_GO_Untreated_enriched_terms_plot,
       width  = 7,
       height = 5)


############## Gene-set enrichment analysis (GSEA) of KEGG pathway for DEgenes in monocytes from patients with KFH compared to healthy control presented in Supplementary Table/Table S9 [Figure 6E]
############# Analysis done in GSEA_MonoMacroDCs.R
GSEA_KEGG_Untreated_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                       sheet="Table S9", skip=5)

GSEA_KEGG_Untreated_enriched_terms<-GSEA_KEGG_allSignificantTerms_Patient_vs_Control_withinUntreated
GSEA_KEGG_Untreated_enriched_terms_df<- data.frame(
  pathways = GSEA_KEGG_Untreated_enriched_terms$Description,
  NES = GSEA_KEGG_Untreated_enriched_terms$NES,
  qvalue = GSEA_KEGG_Untreated_enriched_terms$qvalues,
  Gene_Ratio = (unlist(lapply(GSEA_KEGG_Untreated_enriched_terms$core_enrichment, function(x) length(unlist(strsplit(x, split= "\\/" )))))/GSEA_KEGG_Untreated_enriched_terms$setSize)
)

GSEA_KEGG_Untreated_enriched_terms_df$type = "upregulated"
GSEA_KEGG_Untreated_enriched_terms_df$type[GSEA_KEGG_Untreated_enriched_terms_df$NES < 0] = "downregulated"
GSEA_KEGG_Untreated_enriched_terms_df$type<-factor(GSEA_KEGG_Untreated_enriched_terms_df$type, levels=c("upregulated","downregulated"))
GSEA_KEGG_Untreated_enriched_terms_df$log_qvalue <- -log10(GSEA_KEGG_Untreated_enriched_terms_df$qvalue)
GSEA_KEGG_Untreated_enriched_terms_df$Gene_Ratio <- round(GSEA_KEGG_Untreated_enriched_terms_df$Gene_Ratio*100)

GSEA_KEGG_Untreated_enriched_terms_df <- GSEA_KEGG_Untreated_enriched_terms_df %>% arrange(type, Gene_Ratio)
GSEA_KEGG_Untreated_enriched_terms_df$pathways <- factor(GSEA_KEGG_Untreated_enriched_terms_df$pathways, levels=GSEA_KEGG_Untreated_enriched_terms_df$pathways)

####Plot
GSEA_KEGG_Untreated_enriched_terms_plot<- ggplot(GSEA_KEGG_Untreated_enriched_terms_df,
                                   aes(x     = Gene_Ratio,
                                       y     = pathways,
                                       size  = log_qvalue,
                                       color = NES)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(type ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(type = c(upregulated   = "Upregulated",
                                          downregulated = "Downregulated"))) +
  
  scale_colour_gradient2(
    midpoint = 0,
    low      = "darkslateblue",
    mid      = "white",
    high     = "darkred",
    name     = "NES",
    guide    = guide_colorbar(
      barwidth  = 0.8,
      barheight = 6,
      ticks     = TRUE)
  ) +
  
  scale_size_continuous(
    name   = expression(-log[10](italic(q))),
    range  = c(3, 9),
    breaks = pretty(GSEA_KEGG_Untreated_enriched_terms_df$log_qvalue, n = 4)
  ) +
  
  scale_y_discrete(labels = label_wrap(50)) +
  
  labs(x = "Gene Ratio",
       y = NULL) +
  
  theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    # Axes
    axis.text.x      = element_text(size = 9,  color = "black"),
    axis.text.y      = element_text(size = 9,  color = "black"),
    axis.title.x     = element_text(size = 10, color = "black", margin = margin(t = 6)),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    
    # No gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Facet strip styling
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.4),
    strip.text       = element_text(size = 10, color = "black", face = "bold"),
    
    # Legend
    legend.title     = element_text(size = 9,  color = "black"),
    legend.text      = element_text(size = 8,  color = "black"),
    legend.key.size  = unit(0.4, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.box       = "vertical",
    legend.position  = "right",
    legend.margin    = margin(l = 4),
    
    # Margins
    plot.margin      = margin(t = 8, r = 6, b = 8, l = 6)
  )
GSEA_KEGG_Untreated_enriched_terms_plot

ggsave("Figure6E.pdf",
       plot   = GSEA_KEGG_Untreated_enriched_terms_plot,
       width  = 6,
       height = 4)



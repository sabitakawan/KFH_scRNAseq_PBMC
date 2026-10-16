library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

####Gene-set enrichment analysis (GSEA) of GO-BP for DEgenes in monocytes from patients with KFH compared to healthy control [Figure 6D]

GSEA_Untreated_GO<-GSEA_GO_ReducedTerms_similarityScore_Patient_vs_Control_withinUntreated
GSEA_Untreated_GO<-GSEA_Untreated_GO[c(-1,-4,-6:-10,-14:-18,-20,-23,-30),]
GSEA_Untreated_GO_df<- data.frame(
  pathways = GSEA_Untreated_GO$Description,
  NES = GSEA_Untreated_GO$NES,
  qvalue = GSEA_Untreated_GO$qvalues,
  Gene_Ratio = (unlist(lapply(GSEA_Untreated_GO$core_enrichment, function(x) length(unlist(strsplit(x, split= "\\/" )))))/GSEA_Untreated_GO$setSize)
)
GSEA_Untreated_GO_df$type = "upregulated"
GSEA_Untreated_GO_df$type[GSEA_Untreated_GO_df$NES < 0] = "downregulated"
GSEA_Untreated_GO_df$type<-factor(GSEA_Untreated_GO_df$type, levels=c("upregulated","downregulated"))
GSEA_Untreated_GO_df$log_qvalue <- -log10(GSEA_Untreated_GO_df$qvalue)
GSEA_Untreated_GO_df$Gene_Ratio <- round(GSEA_Untreated_GO_df$Gene_Ratio*100)

GSEA_Untreated_GO_df <- GSEA_Untreated_GO_df %>% arrange(type, Gene_Ratio)
GSEA_Untreated_GO_df$pathways <- factor(GSEA_Untreated_GO_df$pathways, levels=GSEA_Untreated_GO_df$pathways)

####PLOT
GSEA_Untreated_GO_BP_plot<- ggplot(GSEA_Untreated_GO_df,
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
    breaks = pretty(GSEA_Untreated_GO_df$log_qvalue, n = 4)
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
GSEA_Untreated_GO_BP_plot

ggsave("GSEA_Untreated_GO_BP_plot_Figure6C.pdf",
       plot   = GSEA_Untreated_GO_BP_plot,
       width  = 7,
       height = 5)


############## Gene-set enrichment analysis (GSEA) of KEGG pathway for DEgenes in monocytes from patients with KFH compared to healthy control [Figure 6E]
GSEA_Untreated_KEGG<-GSEA_KEGG_allSignificantTerms_Patient_vs_Control_withinUntreated
GSEA_Untreated_KEGG_df<- data.frame(
  pathways = GSEA_Untreated_KEGG$Description,
  NES = GSEA_Untreated_KEGG$NES,
  qvalue = GSEA_Untreated_KEGG$qvalues,
  Gene_Ratio = (unlist(lapply(GSEA_Untreated_KEGG$core_enrichment, function(x) length(unlist(strsplit(x, split= "\\/" )))))/GSEA_Untreated_KEGG$setSize)
)

GSEA_Untreated_KEGG_df$type = "upregulated"
GSEA_Untreated_KEGG_df$type[GSEA_Untreated_KEGG_df$NES < 0] = "downregulated"
GSEA_Untreated_KEGG_df$type<-factor(GSEA_Untreated_KEGG_df$type, levels=c("upregulated","downregulated"))
GSEA_Untreated_KEGG_df$log_qvalue <- -log10(GSEA_Untreated_KEGG_df$qvalue)
GSEA_Untreated_KEGG_df$Gene_Ratio <- round(GSEA_Untreated_KEGG_df$Gene_Ratio*100)

GSEA_Untreated_KEGG_df <- GSEA_Untreated_KEGG_df %>% arrange(type, Gene_Ratio)
GSEA_Untreated_KEGG_df$pathways <- factor(GSEA_Untreated_KEGG_df$pathways, levels=GSEA_Untreated_KEGG_df$pathways)

####Plot
GSEA_Untreated_KEGG_plot<- ggplot(GSEA_Untreated_KEGG_df,
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
    breaks = pretty(GSEA_Untreated_KEGG_df$log_qvalue, n = 4)
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
GSEA_Untreated_KEGG_plot

ggsave("GSEA_Untreated_KEGGplot_Figure6E.pdf",
       plot   = GSEA_Untreated_GO_BP_plot,
       width  = 7,
       height = 5)



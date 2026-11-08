
##### Overrepresentation analysis of KEGG pathway for upregulated and downregulated DEGs in LPS stimulated monocytes in patient with KFH compared to healthy control [Figure 7B]
library(Seurat)
library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
library(clusterProfiler)
library("org.Hs.eg.db")
library(ggh4x)
OrgDb <- org.Hs.eg.db
Symbol_col <- 'SYMBOL'
padj.thr <- 1e-03

############Semantic analysis GO-BP for upregulated DEGs in unstimulated monocytes compared between patients with KFH and healthy controls [FIGURE 5B]
## Load Table S3 with GO terms for upregulated DEGs

KFH_LPS_markergenes <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                        sheet="Table S3", skip=3)

if (sum(KFH_LPS_markergenes$significant) >= 100) {
  up_genes_to_use_LPS <- KFH_LPS_markergenes$category=="upregulated" & KFH_LPS_markergenes$avg_log2FC > log2(1.5)
  do_genes_to_use_LPS <- KFH_LPS_markergenes$category=="downregulated" & KFH_LPS_markergenes$avg_log2FC < -log2(1.5)
} else {
  up_genes_to_use_LPS <- KFH_LPS_markergenes$category=="upregulated"
  do_genes_to_use_LPS <- KFH_LPS_markergenes$category=="downregulated"
}
genes_upLPS   <- KFH_LPS_markergenes$gene[up_genes_to_use_LPS]
genes_downLPS <- KFH_LPS_markergenes$gene[do_genes_to_use_LPS]


# Upregulated
df_up_LPS <- bitr(genes_upLPS, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
genes_up_LPS_entrez <- unique(df_up_LPS$ENTREZID)

# Downregulated
df_down_LPS <- bitr(genes_downLPS, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
genes_down_LPS_entrez <- unique(df_down_LPS$ENTREZID)

####Load Uiverse gene set from Supplementary_Table, Table S21
universe_geneset<- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                          sheet="Table S21", skip=2)

LPS_universe <- bitr(universe_geneset$`Universe gene set`, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
universeLPS_entrez <- unique(LPS_universe$ENTREZID)


# Upregulated pathways
if(length(genes_up_LPS_entrez) >= 10){
  ekegg_upLPS <- enrichKEGG(
    gene = genes_up_LPS_entrez,
    universe = universeLPS_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
}

# Downregulated pathways
if(length(genes_down_LPS_entrez) >= 10){
  ekegg_downLPS <- enrichKEGG(
    gene = genes_down_LPS_entrez,
    universe = universeLPS_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
}

ekegg_upLPS   <- setReadable(ekegg_upLPS, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
ekegg_downLPS <- setReadable(ekegg_downLPS, OrgDb=org.Hs.eg.db, keyType="ENTREZID")

# Upregulated
ekegg_upLPS@result[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue")]

# Downregulated
ekegg_downLPS@result[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue")]

sig_upLPS   <- ekegg_upLPS@result[ekegg_upLPS@result$p.adjust < 0.05, ]
sig_downLPS <- ekegg_downLPS@result[ekegg_downLPS@result$p.adjust < 0.05, ]

ekegg_upLPS@result$direction <- "Up"
ekegg_downLPS@result$direction <- "Down"

sig_upLPS <- ekegg_upLPS@result[ekegg_upLPS@result$p.adjust < 0.05, ]
sig_downLPS <- ekegg_downLPS@result[ekegg_downLPS@result$p.adjust < 0.05, ]

combinedLPS_raw <- rbind(ekegg_upLPS@result, ekegg_downLPS@result)

cols <- c("ID", "Description", "GeneRatio", "BgRatio",
          "pvalue", "p.adjust", "qvalue", "geneID", "Count", "direction")

KFH_KEGG_LPS <- combinedLPS_raw[, cols]
####write.csv(combinedLPS[, cols],
###         "KEGG_LPS_mono_clean.csv",
###         row.names = FALSE)


# Check your column name for up/down labels
colnames(KFH_KEGG_LPS)
unique(KFH_KEGG_LPS$direction)    # adjust "regulation" to your actual column name

KFH_KEGG_LPS_Filtered<- KFH_KEGG_LPS |>
  mutate(GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"),
                                    function(x) as.numeric(x[1]) / as.numeric(x[2]))) |>
  filter(p.adjust < 0.05,
         category != "Human Diseases") |>
  mutate(
    Description  = str_wrap(Description, width = 45),
    direction   = factor(direction, levels = c("Up", "Down"))  # adjust to your labels
  ) |>
  arrange(direction, GeneRatio_numeric)

KFH_KEGG_LPS_Filtered<-KFH_KEGG_LPS_Filtered[c(-1:-3,-9,-12, -15,-20),]
KFH_KEGG_LPS_Filtered_plot <- ggplot(KFH_KEGG_LPS_Filtered,
                      aes(x     = GeneRatio_numeric,
                          y     = reorder(Description, GeneRatio_numeric),
                          size  = Count,
                          color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  # Facet by direction — two clean panels
  facet_grid(direction ~ .,
             scales    = "free_y",
             space     = "free_y",
             labeller  = labeller(direction = c(Up   = "Upregulated",
                                                Down = "Downregulated"))) +
  
  # Colour scale
  scale_color_gradient(low  = "#c0392b",
                       high = "#2980b9",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(
                         barwidth  = 0.8,
                         barheight = 6,
                         ticks     = TRUE)) +
  
  # Dot size
  scale_size_continuous(name   = "Gene\ncount",
                        range  = c(3, 9),
                        breaks = pretty(KFH_KEGG_LPS_Filtered$Count, n = 4)) +
  
  labs(x = "Gene ratio",
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

KFH_KEGG_LPS_Filtered_plot

ggsave("Figure5C.pdf",
       plot   = KFH_KEGG_LPS_Filtered_plot,
       width  = 7,
       height = 5)   



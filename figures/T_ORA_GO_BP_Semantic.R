###Libraries needed for the analysis
library(tidyr)
library(dplyr)
library(ggplot2)
library(rrvgo)
library(GOSemSim)
library(stringr)

############Semantic analysis GO-BP for upregulated DEGs in unstimulated T cells compared between patients with KFH and healthy controls. Data is presented in Supplementary Table/Table S8  [FIGURE 6C]
## Load Table S8B with GO terms for upregulated DEGs analysed [analysis done in diffexp_and_ORA_T.R]

Untreated_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                        sheet="Table S8B", skip=3)
simMatrix_Untreated_upreg<- calculateSimMatrix(Untreated_GO_upreg_enriched_terms$ID,
                                               orgdb="org.Hs.eg.db",
                                               ont="BP",
                                               method="Rel")
scores_Untreated_upreg <- setNames(-log10(Untreated_GO_upreg_enriched_terms$qvalue), Untreated_GO_upreg_enriched_terms$ID)
reducedTerms_Untreated_upreg <- reduceSimMatrix(simMatrix_Untreated_upreg,
                                                scores_Untreated_upreg,
                                                threshold=0.7,
                                                orgdb="org.Hs.eg.db")
colnames(reducedTerms_Untreated_upreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_Untreated_upreg <- reducedTerms_Untreated_upreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_Untreated_upreg <- Untreated_GO_upreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_Untreated_GO_upreg <- merge(rrvgo_mapping_Untreated_upreg, original_data_subset_Untreated_upreg, by = "ID", all.x = TRUE)

# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
parents_Untreated_upreg <- final_publication_supplement_Untreated_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")

top_parents_Untreated <- parents_Untreated_upreg %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot

GO_term_Untreated_T<- ggplot(top_parents_Untreated,
                           aes(x     = ratio,
                               y     = reorder(parentTerm, ratio),
                               size  = Count,
                               color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(Regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(Regulation = "Upregulated"
                                                )) +
  
  scale_color_gradient(low  = "#c0392b",
                       high = "#2980b9",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(
                         barwidth  = 0.8,
                         barheight = 6,
                         ticks     = TRUE)) +
  
  scale_size_continuous(name   = "Gene\ncount",
                        range  = c(3, 9),
                        breaks = pretty(top_parents_Untreated$Count, n = 4)) +
  
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
GO_term_Untreated_T
ggsave("Figure6D.pdf",
       plot   = GO_term_Untreated_T,
       width  = 6.5,
       height = 4.5)


################Semantic analysis GO-BP for upregulated DEGs in LPS-stimulated T cells compared between patients with KFH and healthy controls. Data is presented in Supplementary Table/Table S13  [FIGURE 7C]
## Load Table S13B with GO terms for upregulated DEGs analysed [analysis done in diffexp_and_ORA_T.R]

LPS_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                        sheet="Table S13B", skip=3)
simMatrix_LPS_upreg<- calculateSimMatrix(LPS_GO_upreg_enriched_terms$ID,
                                               orgdb="org.Hs.eg.db",
                                               ont="BP",
                                               method="Rel")
scores_LPS_upreg <- setNames(-log10(LPS_GO_upreg_enriched_terms$qvalue), LPS_GO_upreg_enriched_terms$ID)
reducedTerms_LPS_ATP_upreg <- reduceSimMatrix(simMatrix_LPS_upreg,
                                                scores_LPS_upreg,
                                                threshold=0.7,
                                                orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_upreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_LPS_upreg <- reducedTerms_LPS_upreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_LPS_upreg <- LPS_GO_upreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_LPS_GO_upreg <- merge(rrvgo_mapping_LPS_upreg, original_data_subset_LPS_upreg, by = "ID", all.x = TRUE)


# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
parents_LPS_upreg <- final_publication_supplement_LPS_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")


top_parents_LPS <- parents_LPS_upreg %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot GO terms in a single plot

GO_term_LPS_T<- ggplot(top_parents_LPS,
                             aes(x     = ratio,
                                 y     = reorder(parentTerm, ratio),
                                 size  = Count,
                                 color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(Regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(Regulation = "Upregulated"
             )) +
  
  scale_color_gradient(low  = "#c0392b",
                       high = "#2980b9",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(
                         barwidth  = 0.8,
                         barheight = 6,
                         ticks     = TRUE)) +
  
  scale_size_continuous(name   = "Gene\ncount",
                        range  = c(3, 9),
                        breaks = pretty(top_parents_LPS$Count, n = 4)) +
  
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
GO_term_LPS_T
ggsave("Figure7D.pdf",
       plot   = GO_term_LPS_T,
       width  = 7.5,
       height = 4)

################Semantic analysis GO-BP for upregulated DEGs in LPS+ATP-stimulated T cells compared between patients with KFH and healthy controls. Data is presented in Supplementary Table/Table S18  [FIGURE 8C]
## Load Table S18B with GO terms for upregulated DEGs analysed [analysis done in diffexp_and_ORA_T.R]


LPS_ATP_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                  sheet="Table S18B", skip=3)
simMatrix_LPS_ATP_upreg<- calculateSimMatrix(LPS_ATP_GO_upreg_enriched_terms$ID,
                                         orgdb="org.Hs.eg.db",
                                         ont="BP",
                                         method="Rel")
scores_LPS_ATP_upreg <- setNames(-log10(LPS_ATP_GO_upreg_enriched_terms$qvalue), LPS_ATP_GO_upreg_enriched_terms$ID)
reducedTerms_LPS_ATP_upreg <- reduceSimMatrix(simMatrix_LPS_ATP_upreg,
                                          scores_LPS_ATP_upreg,
                                          threshold=0.7,
                                          orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_ATP_upreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_LPS_ATP_upreg <- reducedTerms_LPS_ATP_upreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_LPS_ATP_upreg <- LPS_ATP_GO_upreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_LPS_ATP_GO_upreg <- merge(rrvgo_mapping_LPS_ATP_upreg, original_data_subset_LPS_ATP_upreg, by = "ID", all.x = TRUE)


# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
# (This extracts the clean Gene Ratio and Score for just the parent terms)

parents_LPS_ATP_upreg <- final_publication_supplement_LPS_ATP_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")


top_parents_LPS_ATP <- parents_LPS_ATP_upreg %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot

GO_term_LPS_ATP_T<- ggplot(top_parents_LPS_ATP,
                       aes(x     = ratio,
                           y     = reorder(parentTerm, ratio),
                           size  = Count,
                           color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(Regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(Regulation = "Upregulated"
             )) +
  
  scale_color_gradient(low  = "#c0392b",
                       high = "#2980b9",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(
                         barwidth  = 0.8,
                         barheight = 6,
                         ticks     = TRUE)) +
  
  scale_size_continuous(name   = "Gene\ncount",
                        range  = c(3, 9),
                        breaks = pretty(top_parents_LPS_ATP$Count, n = 4)) +
  
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
GO_term_LPS_ATP_T
ggsave("Figure8C.pdf",
       plot   = GO_term_LPS_ATP_T,
       width  = 6,
       height = 4)

library(tidyr)
library(dplyr)
library(ggplot2)
library(rrvgo)
library(GOSemSim)
library(stringr)

#####ORA_GO_BP enrichment analysis for upregulated and downregulated DEGs within LPS and LPS+ATP stimulated monocytes compared to their unstimulated condition in patients with KFH and healthy controls separately. 
####Analysis pipeline is presented in diffexp_and_ORA.R

###########Semantic analysis GO-BP for downregulated DEGs within LPS stimulated monocytes compared to unstimulated monocytes from patients with KFH [FIGURES4B]
## Load Table S14B for reproducing Figure S4B
WithinKFH_LPS_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                        sheet="Table S14B", skip=3)
simMatrix_withinKFH_LPS_upreg<- calculateSimMatrix(WithinKFH_LPS_GO_upreg_enriched_terms$ID,
                                               orgdb="org.Hs.eg.db",
                                               ont="BP",
                                               method="Rel")
scores_withinKFH_LPS_upreg <- setNames(-log10(WithinKFH_LPS_GO_upreg_enriched_terms$qvalue), WithinKFH_LPS_GO_upreg_enriched_terms$ID)
reducedTerms_withinKFH_LPS_upreg <- reduceSimMatrix(simMatrix_withinKFH_LPS_upreg,
                                                   scores_withinKFH_LPS_upreg,
                                                   threshold=0.7,
                                                   orgdb="org.Hs.eg.db")
colnames(reducedTerms_withinKFH_LPS_upreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_reducedTerms_withinKFH_LPS_upreg <- reducedTerms_withinKFH_LPS_upreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_withinKFH_LPS_upreg <- WithinKFH_LPS_GO_upreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_WithinKFH_LPS_GO_upreg_enriched_terms_GO_upreg <- merge(rrvgo_mapping_reducedTerms_withinKFH_LPS_upreg, original_data_subset_withinKFH_LPS_upreg, by = "ID", all.x = TRUE)

# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
parents_WithinKFH_LPS_GO_upreg <- final_publication_supplement_WithinKFH_LPS_GO_upreg_enriched_terms_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")

top_parents_WithinKFH_LPS_GO_upreg <- parents_WithinKFH_LPS_GO_upreg %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot

top_parents_WithinKFH_LPS_GO_upreg_plot<- ggplot(top_parents_WithinKFH_LPS_GO_upreg,
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
                        breaks = pretty(top_parents_WithinKFH_LPS_GO_upreg$Count, n = 4)) +
  
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
top_parents_WithinKFH_LPS_GO_upreg_plot
ggsave("FigureS4B.pdf",
       plot   = top_parents_WithinKFH_LPS_GO_upreg_plot,
       width  = 7,
       height = 4)

###########Semantic analysis GO-BP for downregulated DEGs within LPS stimulated monocytes compared to unstimulated monocytes from Healthy control [FIGURES4C]
## Load Table S15B for reproducing Figure S4C
WithinHC_LPS_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                            sheet="Table S15B", skip=3)
simMatrix_withinHC_LPS_upreg<- calculateSimMatrix(WithinHC_LPS_GO_upreg_enriched_terms$ID,
                                                   orgdb="org.Hs.eg.db",
                                                   ont="BP",
                                                   method="Rel")
scores_withinHC_LPS_upreg <- setNames(-log10(WithinHC_LPS_GO_upreg_enriched_terms$qvalue), WithinHC_LPS_GO_upreg_enriched_terms$ID)
reducedTerms_withinHC_LPS_upreg <- reduceSimMatrix(simMatrix_withinHC_LPS_upreg,
                                                    scores_withinHC_LPS_upreg,
                                                    threshold=0.7,
                                                    orgdb="org.Hs.eg.db")
colnames(reducedTerms_withinHC_LPS_upreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_reducedTerms_withinHC_LPS_upreg <- reducedTerms_withinHC_LPS_upreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_withinHC_LPS_upreg <- WithinHC_LPS_GO_upreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_WithinHC_LPS_GO_upreg_enriched_terms_GO_upreg <- merge(rrvgo_mapping_reducedTerms_withinHC_LPS_upreg, original_data_subset_withinHC_LPS_upreg, by = "ID", all.x = TRUE)

# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
parents_WithinHC_LPS_GO_upreg <- final_publication_supplement_WithinHC_LPS_GO_upreg_enriched_terms_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")

top_parents_WithinHC_LPS_GO_upreg <- parents_WithinHC_LPS_GO_upreg %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot

top_parents_WithinHC_LPS_GO_upreg_plot<- ggplot(top_parents_WithinHC_LPS_GO_upreg,
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
                        breaks = pretty(top_parents_WithinHC_LPS_GO_upreg$Count, n = 4)) +
  
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
top_parents_WithinHC_LPS_GO_upreg_plot
ggsave("FigureS4C.pdf",
       plot   = top_parents_WithinHC_LPS_GO_upreg_plot,
       width  = 7,
       height = 4)



###########Semantic analysis GO-BP for downregulated DEGs within LPS+ATP stimulated monocytes compared to unstimulated monocytes from patient with KFH [FIGURES4D]
## Load Table S20B for reproducing Figure S4D
WithinKFH_LPS_ATP_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                            sheet="Table S20B", skip=3)
simMatrix_withinKFH_LPS_ATP_upreg<- calculateSimMatrix(WithinKFH_LPS_ATP_GO_upreg_enriched_terms$ID,
                                                   orgdb="org.Hs.eg.db",
                                                   ont="BP",
                                                   method="Rel")
scores_withinKFH_LPS_ATP_upreg <- setNames(-log10(WithinKFH_LPS_ATP_GO_upreg_enriched_terms$qvalue), WithinKFH_LPS_ATP_GO_upreg_enriched_terms$ID)
reducedTerms_withinKFH_LPS_ATP_upreg <- reduceSimMatrix(simMatrix_withinKFH_LPS_ATP_upreg,
                                                    scores_withinKFH_LPS_ATP_upreg,
                                                    threshold=0.7,
                                                    orgdb="org.Hs.eg.db")
colnames(reducedTerms_withinKFH_LPS_ATP_upreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_reducedTerms_withinKFH_LPS_ATP_upreg <- reducedTerms_withinKFH_LPS_ATP_upreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_withinKFH_LPS_ATP_upreg <- WithinKFH_LPS_ATP_GO_upreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_WithinKFH_LPS_ATP_GO_upreg_enriched_terms_GO_upreg <- merge(rrvgo_mapping_reducedTerms_withinKFH_LPS_ATP_upreg, original_data_subset_withinKFH_LPS_ATP_upreg, by = "ID", all.x = TRUE)

# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
# (This extracts the clean Gene Ratio and Score for just the parent terms)
parents_WithinKFH_LPS_ATP_GO_upreg <- final_publication_supplement_WithinKFH_LPS_ATP_GO_upreg_enriched_terms_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")

top_parents_WithinKFH_LPS_ATP_GO_upreg <- parents_WithinKFH_LPS_ATP_GO_upreg %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot

top_parents_WithinKFH_LPS_ATP_GO_upreg_plot<- ggplot(top_parents_WithinKFH_LPS_ATP_GO_upreg,
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
                        breaks = pretty(top_parents_WithinKFH_LPS_ATP_GO_upreg$Count, n = 4)) +
  
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
top_parents_WithinKFH_LPS_ATP_GO_upreg_plot
ggsave("FigureS4D.pdf",
       plot   = top_parents_WithinKFH_LPS_ATP_GO_upreg_plot,
       width  = 7,
       height = 4)

###########Semantic analysis GO-BP for downregulated DEGs within LPS+ATP stimulated monocytes compared to unstimulated monocytes from Healthy control [FIGURES4E]
## Load Table S19B for reproducing Figure S4E
WithinHC_LPS_ATP_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                           sheet="Table S19B", skip=3)
simMatrix_withinHC_LPS_ATP_upreg<- calculateSimMatrix(WithinHC_LPS_ATP_GO_upreg_enriched_terms$ID,
                                                  orgdb="org.Hs.eg.db",
                                                  ont="BP",
                                                  method="Rel")
scores_withinHC_LPS_ATP_upreg <- setNames(-log10(WithinHC_LPS_ATP_GO_upreg_enriched_terms$qvalue), WithinHC_LPS_ATP_GO_upreg_enriched_terms$ID)
reducedTerms_withinHC_LPS_ATP_upreg <- reduceSimMatrix(simMatrix_withinHC_LPS_ATP_upreg,
                                                   scores_withinHC_LPS_ATP_upreg,
                                                   threshold=0.7,
                                                   orgdb="org.Hs.eg.db")
colnames(reducedTerms_withinHC_LPS_ATP_upreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_reducedTerms_withinHC_LPS_ATP_upreg <- reducedTerms_withinHC_LPS_ATP_upreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_withinHC_LPS_ATP_upreg <- WithinHC_LPS_ATP_GO_upreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_WithinHC_LPS_ATP_GO_upreg_enriched_terms_GO_upreg <- merge(rrvgo_mapping_reducedTerms_withinHC_LPS_ATP_upreg, original_data_subset_withinHC_LPS_ATP_upreg, by = "ID", all.x = TRUE)

# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
# (This extracts the clean Gene Ratio and Score for just the parent terms)
parents_WithinHC_LPS_ATP_GO_upreg <- final_publication_supplement_WithinHC_LPS_ATP_GO_upreg_enriched_terms_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")

top_parents_WithinHC_LPS_ATP_GO_upreg <- parents_WithinHC_LPS_ATP_GO_upreg %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot

top_parents_WithinHC_LPS_ATP_GO_upreg_plot<- ggplot(top_parents_WithinHC_LPS_ATP_GO_upreg,
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
                        breaks = pretty(top_parents_WithinHC_LPS_ATP_GO_upreg$Count, n = 4)) +
  
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
top_parents_WithinHC_LPS_ATP_GO_upreg_plot
ggsave("FigureS4E.pdf",
       plot   = top_parents_WithinHC_LPS_ATP_GO_upreg_plot,
       width  = 7,
       height = 4)

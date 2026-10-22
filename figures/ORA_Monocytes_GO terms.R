###Libraries needed for the analysis
library(tidyr)
library(dplyr)
library(ggplot2)
library(rrvgo)
library(GOSemSim)
library(stringr)

############Semantic analysis GO-BP for upregulated DEGs in unstimulated monocytes compared between patients with KFH and healthy controls. Data is presented in Supplementary Table/Table S5 [FIGURE 6B]
## Load Table S5B with GO terms for upregulated DEGs analysed [analysis done in diffexp_and_ORA.R]

Untreated_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                                              sheet="Table S5B", skip=3)
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


#####Semantic analysis GO-BP for upregulated DEGs in unstimulated monocytes compared between patients with KFH and healthy controls Data is presented in Supplementary Table/Table S6 [FIGURE 6B]
## Load Table S6B with GO terms for downregulated DEGs analysed [analysis done in analysis done in diffexp_and_ORA.R]

Untreated_GO_downreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                        sheet="Table S6B", skip=3)
simMatrix_Untreated_downreg<- calculateSimMatrix(Untreated_GO_downreg_enriched_terms$ID,
                                               orgdb="org.Hs.eg.db",
                                               ont="BP",
                                               method="Rel")
scores_Untreated_downreg <- setNames(-log10(Untreated_GO_downreg_enriched_terms$qvalue), Untreated_GO_downreg_enriched_terms$ID)
reducedTerms_Untreated_downreg <- reduceSimMatrix(simMatrix_Untreated_downreg,
                                                scores_Untreated_downreg,
                                                threshold=0.7,
                                                orgdb="org.Hs.eg.db")
colnames(reducedTerms_Untreated_downreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_Untreated_downreg <- reducedTerms_Untreated_downreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_Untreated_downreg <- Untreated_GO_downreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_Untreated_GO_downreg <- merge(rrvgo_mapping_Untreated_downreg, original_data_subset_Untreated_downreg, by = "ID", all.x = TRUE)

###ggplot

# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
parents_Untreated_downreg <- final_publication_supplement_Untreated_GO_downreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Downregulated")

parents_Untreated_upreg <- final_publication_supplement_Untreated_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")

###Keeping the GO-terms that are relevant for visualization 
parents_up<-parents_up[c(-1:-3,-8:-11,-13,-14,-16,-21,-21,-23,-24),]

# Combine Up and Down into a single dataframe
combined_parents_Untreated <- bind_rows(parents_Untreated_upreg, parents_Untreated_downreg)

top_parents_Untreated <- combined_parents %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot

GO_term_Untreated<- ggplot(top_parents_Untreated,
                    aes(x     = ratio,
                        y     = reorder(parentTerm, ratio),
                        size  = Count,
                        color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(Regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(Regulation = c(Upregulated   = "Upregulated",
                                                Downregulated = "Downregulated"))) +
  
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
GO_term_Untreated
ggsave("Figure6B.pdf",
       plot   = GO_term_Untreated,
       width  = 7,
       height = 8)


############Semantic analysis GO-BP for upregulated DEGs in LPS-stimulated monocytes compared between patients with KFH and healthy controls. Data is presented in Supplementary Table/Table S12 [FIGURE 6B]
## Load Table S12B with GO terms for upregulated DEGs analysed [analysis done in diffexp_and_ORA.R]

LPS_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                        sheet="Table S12B", skip=3)
simMatrix_LPS_upreg<- calculateSimMatrix(LPS_GO_upreg_enriched_terms$ID,
                                               orgdb="org.Hs.eg.db",
                                               ont="BP",
                                               method="Rel")
scores_LPS_upreg <- setNames(-log10(LPS_GO_upreg_enriched_terms$qvalue), LPS_GO_upreg_enriched_terms$ID)
reducedTerms_LPS_upreg <- reduceSimMatrix(simMatrix_LPS_upreg,
                                                scores_LPS_upreg,
                                                threshold=0.7,
                                                orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_upreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_LPS_upreg <- reducedTerms_LPS_upreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_LPS_upreg <- LPS_GO_upreg_enriched_terms[, c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_LPS_GO_upreg <- merge(rrvgo_mapping_LPS_upreg, original_data_subset_LPS_upreg, by = "ID", all.x = TRUE)


#####Semantic analysis GO-BP for upregulated DEGs in unstimulated monocytes compared between patients with KFH and healthy controls. Data is presented in Supplementary Table/ Table S10 [FIGURE 6B]
## Load Table S10B with GO terms for downregulated DEGs analysed [analysis done in diffexp_and_ORA.R]

LPS_GO_downreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                          sheet="Table S10B", skip=3)
simMatrix_LPS_downreg<- calculateSimMatrix(LPS_GO_downreg_enriched_terms$ID,
                                                 orgdb="org.Hs.eg.db",
                                                 ont="BP",
                                                 method="Rel")
scores_LPS_downreg <- setNames(-log10(LPS_GO_downreg_enriched_terms$qvalue), LPS_GO_downreg_enriched_terms$ID)
reducedTerms_LPS_downreg <- reduceSimMatrix(simMatrix_LPS_downreg,
                                                  scores_LPS_downreg,
                                                  threshold=0.7,
                                                  orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_downreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_LPS_downreg <- reducedTerms_LPS_downreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_LPS_downreg <- LPS_GO_downreg_enriched_terms[, c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_LPS_GO_downreg <- merge(rrvgo_mapping_LPS_downreg, original_data_subset_LPS_downreg, by = "ID", all.x = TRUE)

###ggplot

# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
parents_LPS_downreg <- final_publication_supplement_LPS_GO_downreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Downregulated")

###Keeping the GO-terms that are relevant for visualization 
parents_LPS_downreg <- parents_LPS_downreg[c(-2,-5:-7,-9:-14,-17,-19,-21,-23,-24,-27:-30),]

parents_LPS_upreg <- final_publication_supplement_LPS_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")

###Keeping the GO-terms that are relevant for visualization 
parents_LPS_upreg<-parents_LPS_upreg[c(-3,-7,-9),]

# Combine Up and Down into a single dataframe
combined_parents_LPS <- bind_rows(parents_LPS_upreg, parents_LPS_downreg)

top_parents_LPS <- combined_parents_LPS %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 20) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot
# Convert Regulation to a factor with Upregulated first
top_parents_LPS$Regulation <- factor(top_parents_LPS$Regulation, 
                                     levels = c("Upregulated", "Downregulated"))

GO_term_LPS<- ggplot(top_parents_LPS,
                           aes(x     = ratio,
                               y     = reorder(parentTerm, ratio),
                               size  = Count,
                               color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(Regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(Regulation = c(Upregulated   = "Upregulated",
                                                Downregulated = "Downregulated"))) +
  
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
GO_term_LPS
ggsave("Figure7B.pdf",
       plot   = GO_term_LPS,
       width  = 9,
       height = 8)


#####
############Semantic analysis GO-BP for upregulated DEGs in LPS+ATP-stimulated monocytes compared between patients with KFH and healthy controls. Data is presented in Supplementary Table/ Table S17[FIGURE 6B]
## Load Table S17B with GO terms for upregulated DEGs analysed [analysis done in diffexp_and ORA.R]

LPS_ATP_GO_upreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                  sheet="Table S17B", skip=3)
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
original_data_subset_LPS_ATP_upreg <- LPS_ATP_GO_upreg_enriched_terms[, c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_LPS_ATP_GO_upreg <- merge(rrvgo_mapping_LPS_ATP_upreg, original_data_subset_LPS_ATP_upreg, by = "ID", all.x = TRUE)


#####Semantic analysis GO-BP for upregulated DEGs in unstimulated monocytes compared between patients with KFH and healthy controls. Data is presented in Supplementary Table/ Table S17 [FIGURE 6B]
## Load Table S16B with GO terms for downregulated DEGs analysed [analysis done in diffexp_and ORA.R]

LPS_ATP_GO_downreg_enriched_terms <- readxl::read_excel("~/Desktop/Supplementary_Table.xlsx",
                                                    sheet="Table S16B", skip=3)
simMatrix_LPS_ATP_downreg<- calculateSimMatrix(LPS_ATP_GO_downreg_enriched_terms$ID,
                                           orgdb="org.Hs.eg.db",
                                           ont="BP",
                                           method="Rel")
scores_LPS_ATP_downreg <- setNames(-log10(LPS_ATP_GO_downreg_enriched_terms$qvalue), LPS_ATP_GO_downreg_enriched_terms$ID)
reducedTerms_LPS_ATP_downreg <- reduceSimMatrix(simMatrix_LPS_ATP_downreg,
                                            scores_LPS_ATP_downreg,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_ATP_downreg)[1]<-"ID"

# 2. Extract the core clustering columns from your rrvgo output
rrvgo_mapping_LPS_ATP_downreg <- reducedTerms_LPS_ATP_downreg[, c("ID", "term", "parent", "parentTerm", "score")]

# 3. Extract the statistical values and the gene list from your original data frame
original_data_subset_LPS_ATP_downreg <- LPS_ATP_GO_downreg_enriched_terms[, c("ID","Description", "GeneRatio", "BgRatio", "p.adjust", "qvalue", "geneID", "Count", "ratio")]

# 4. Merge them together using the GO ID as the unique key
final_publication_supplement_LPS_ATP_GO_downreg <- merge(rrvgo_mapping_LPS_ATP_downreg, original_data_subset_LPS_ATP_downreg, by = "ID", all.x = TRUE)

###ggplot

# Filter your rrvgo outputs to keep ONLY the representative Parent Rows
# (This extracts the clean Gene Ratio and Score for just the parent terms)
parents_LPS_ATP_downreg <- final_publication_supplement_LPS_ATP_GO_downreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Downregulated")

###Keeping the GO-terms that are relevant for visualization 
parents_LPS_ATP_downreg<-parents_LPS_ATP_downreg[c(-1,-7,-8,-14,-16,-17),]

parents_LPS_ATP_upreg <- final_publication_supplement_LPS_ATP_GO_upreg %>% 
  filter(ID == parent) %>% 
  mutate(Regulation = "Upregulated")

###Keeping the GO-terms that are relevant for visualization 
parents_LPS_ATP_upreg<-parents_LPS_ATP_upreg[c(-1:-3,-9,-12,-13,-15:-18,-21,-26:-33,),]

# Combine Up and Down into a single dataframe
combined_parents_LPS_ATP <- bind_rows(parents_LPS_ATP_upreg, parents_LPS_ATP_downreg)

top_parents_LPS_ATP <- combined_parents_LPS_ATP %>%
  group_by(Regulation) %>%
  slice_max(order_by = score, n = 10) %>% # Keeps the top 10 highest-scoring terms per group
  ungroup()

#######Plot upregulated and downregulated GO terms in a single plot

GO_term_LPS_ATP<- ggplot(top_parents_LPS_ATP,
                     aes(x     = ratio,
                         y     = reorder(parentTerm, ratio),
                         size  = Count,
                         color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(Regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(Regulation = c(Upregulated   = "Upregulated",
                                                Downregulated = "Downregulated"))) +
  
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
GO_term_LPS_ATP
ggsave("Figure8B.pdf",
       plot   = GO_term_LPS_ATP,
       width  = 7,
       height = 8)
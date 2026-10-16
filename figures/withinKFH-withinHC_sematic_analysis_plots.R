library(tidyr)
library(dplyr)
library(ggplot2)
library(rrvgo)
library(GOSemSim)
library(stringr)



###########Semantic analysis GO-BP for downregulated DEGs within LPS stimulated monocytes compared to unstimulated monocytes from patients with KFH [FIGURES4B]

## Load Table S15 for reproducing Figure S4B
turunen_scRNAseq_within_Patient_LPS_GO_downreg_enriched_terms_26072023 <- readxl::read_excel("~/Desktop/Manuscript_scRNA/Manucript_SK_coauthors/Supplementary Data/Supplementary_Table.xlsx",
                                                                                             sheet="Table S15", skip=3)

simMatrix_LPS_down<- calculateSimMatrix(turunen_scRNAseq_within_Patient_LPS_GO_downreg_enriched_terms_26072023$ID,
                                              orgdb="org.Hs.eg.db",
                                              ont="BP",
                                              method="Rel")

scores_LPS_down <- setNames(-log10(turunen_scRNAseq_within_Patient_LPS_GO_downreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_within_Patient_LPS_GO_downreg_enriched_terms_26072023$ID)
reducedTerms_LPS_down <- reduceSimMatrix(simMatrix_LPS_down,
                                               scores_LPS_down,
                                               threshold=0.7,
                                               orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_down)[1]<-"ID"
parent_term_LPS_down<-unique(reducedTerms_LPS_down$parentTerm)

#write.csv(reducedTerms_LPS_down, file = "reducedTerms_withinKFH_LPS_down.csv", row.names = FALSE)

withinKFH_LPS_down <- turunen_scRNAseq_within_Patient_LPS_GO_downreg_enriched_terms_26072023[which(turunen_scRNAseq_within_Patient_LPS_GO_downreg_enriched_terms_26072023$Description %in% parent_term_LPS_down),]
#write.csv(withinKFH_LPS_down, file = "sematic_withinKFH_LPS_down.csv", row.names = FALSE)


###########Semantic analysis GO-BP for upregulated DEGs within LPS stimulated monocytes compared to unstimulated monocytes from patients with KFH [FIGURES4B]

turunen_scRNAseq_within_Patient_LPS_GO_upreg_enriched_terms_26072023 <- readxl::read_excel("~/Desktop/Manuscript_scRNA/Manucript_SK_coauthors/Supplementary Data/Supplementary_Table.xlsx",
                                                                                             sheet="Table S15", skip=3)

simMatrix_LPS_up<- calculateSimMatrix(turunen_scRNAseq_within_Patient_LPS_GO_upreg_enriched_terms_26072023$ID,
                                            orgdb="org.Hs.eg.db",
                                            ont="BP",
                                            method="Rel")
scores_LPS_up <- setNames(-log10(turunen_scRNAseq_within_Patient_LPS_GO_upreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_within_Patient_LPS_GO_upreg_enriched_terms_26072023$ID)
reducedTerms_LPS_up <- reduceSimMatrix(simMatrix_LPS_up,
                                             scores_LPS_up,
                                             threshold=0.7,
                                             orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_up)[1]<-"ID"
parent_term_LPS_up<-unique(reducedTerms_LPS_up$parentTerm)

write.csv(reducedTerms_LPS_up, file = "reducedTerms_withinLPS_up.csv", row.names = FALSE)

withinKFH_LPS_Up <- turunen_scRNAseq_within_Patient_LPS_GO_upreg_enriched_terms_26072023[which(turunen_scRNAseq_within_Patient_LPS_GO_upreg_enriched_terms_26072023$Description %in% parent_term_LPS_up),]
write.csv(withinKFH_LPS_Up, file = "sematic_withinKFH_LPS_Up.csv", row.names = FALSE)

#####Plot the GO-terms (upregulated and downregulated in the same plot)

# Add regulation label to each
withinKFH_LPS_down$regulation   <- "down"
withinKFH_LPS_Up$regulation   <- "up"

withinKFH_LPS_Up<-withinKFH_LPS_Up[c(-11:-79),]
withinKFH_LPS_down<-withinKFH_LPS_down[c(-11:-56),]
# Combine
withinKFH_GO_LPS_df <- rbind(withinKFH_LPS_Up, withinKFH_LPS_down)

# Verify
unique(withinKFH_GO_LPS_df$regulation)  # should show both "up" and "down"


withinKFH_GO_LPS_df_filtered <- withinKFH_GO_LPS_df |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = c("up", "down"))  # adjust to your labels
  )

####Dotplot
withinKFH_GO_LPS<- ggplot(withinKFH_GO_LPS_df_filtered,
                          aes(x     = ratio,
                              y     = reorder(Description, ratio),
                              size  = Count,
                              color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(regulation = c(up   = "Upregulated",
                                                down = "Downregulated"))) +
  
  scale_color_gradient(low  = "#c0392b",
                       high = "#2980b9",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(
                         barwidth  = 0.8,
                         barheight = 6,
                         ticks     = TRUE)) +
  
  scale_size_continuous(name   = "Gene\ncount",
                        range  = c(3, 9),
                        breaks = pretty(withinKFH_GO_LPS_df_filtered$Count, n = 4)) +
  
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
withinKFH_GO_LPS

## Figure S4B ("KFH_GO_LPS)

ggsave("FigureS4B.pdf",
       plot   = withinKFH_GO_LPS,
       width  = 7,
       height = 8)


#######

###########Semantic analysis GO-BP for downregulated DEGs within LPS stimulated monocytes compared to unstimulated monocytes from patients with KFH [FIGURES4B]

simMatrix_LPS_ATP_down<- calculateSimMatrix(turunen_scRNAseq_within_Patient_LPS_ATP_GO_downreg_enriched_terms_26072023$ID,
                                        orgdb="org.Hs.eg.db",
                                        ont="BP",
                                        method="Rel")
scores_LPS_ATP_down <- setNames(-log10(turunen_scRNAseq_within_Patient_LPS_ATP_GO_downreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_within_Patient_LPS_ATP_GO_downreg_enriched_terms_26072023$ID)
reducedTerms_LPS_ATP_down <- reduceSimMatrix(simMatrix_LPS_ATP_down,
                                         scores_LPS_ATP_down,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_ATP_down)[1]<-"ID"
parent_term_LPS_ATP_down<-unique(reducedTerms_LPS_ATP_down$parentTerm)

write.csv(reducedTerms_LPS_ATP_down, file = "reducedTerms_withinKFH_LPS_ATP_down.csv", row.names = FALSE)

withinKFH_LPS_ATP_down <- turunen_scRNAseq_within_Patient_LPS_ATP_GO_downreg_enriched_terms_26072023[which(turunen_scRNAseq_within_Patient_LPS_ATP_GO_downreg_enriched_terms_26072023$Description %in% parent_term_LPS_ATP_down),]
write.csv(withinKFH_LPS_ATP_down, file = "sematic_withinKFH_LPS_ATP_down.csv", row.names = FALSE)


###########Semantic analysis GO-BP for upregulated DEGs within LPS+ATP stimulated monocytes compared to unstimulated monocytes from patients with KFH [FIGURES5B]


simMatrix_LPS_ATP_up<- calculateSimMatrix(turunen_scRNAseq_within_Patient_LPS_ATP_GO_upreg_enriched_terms_26072023$ID,
                                      orgdb="org.Hs.eg.db",
                                      ont="BP",
                                      method="Rel")
scores_LPS_ATP_up <- setNames(-log10(turunen_scRNAseq_within_Patient_LPS_ATP_GO_upreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_within_Patient_LPS_ATP_GO_upreg_enriched_terms_26072023$ID)
reducedTerms_LPS_ATP_up <- reduceSimMatrix(simMatrix_LPS_ATP_up,
                                       scores_LPS_ATP_up,
                                       threshold=0.7,
                                       orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_ATP_up)[1]<-"ID"
parent_term_LPS_ATP_up<-unique(reducedTerms_LPS_ATP_up$parentTerm)

write.csv(reducedTerms_LPS_ATP_up, file = "reducedTerms_withinLPS_ATP_up.csv", row.names = FALSE)

withinKFH_LPS_ATP_Up <- turunen_scRNAseq_within_Patient_LPS_ATP_GO_upreg_enriched_terms_26072023[which(turunen_scRNAseq_within_Patient_LPS_ATP_GO_upreg_enriched_terms_26072023$Description %in% parent_term_LPS_ATP_up),]
write.csv(withinKFH_LPS_ATP_Up, file = "sematic_withinKFH_LPS_ATP_Up.csv", row.names = FALSE)

#####Plot the GO-terms (upregulated and downregulated in the same plot)

# Add regulation label to each
withinKFH_LPS_ATP_down$regulation   <- "down"
withinKFH_LPS_ATP_Up$regulation   <- "up"

withinKFH_LPS_ATP_Up<-withinKFH_LPS_ATP_Up[c(-11:-82),]
withinKFH_LPS_ATP_down<-withinKFH_LPS_ATP_down[c(-11:-57),]
# Combine
withinKFH_GO_LPS_ATP_df <- rbind(withinKFH_LPS_ATP_Up, withinKFH_LPS_ATP_down)

# Verify
unique(withinKFH_GO_LPS_ATP_df$regulation)  # should show both "up" and "down"


withinKFH_GO_LPS_ATP_df_filtered <- withinKFH_GO_LPS_ATP_df |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = c("up", "down"))  # adjust to your labels
  )

####Dotplot
withinKFH_GO_LPS_ATP<- ggplot(withinKFH_GO_LPS_ATP_df_filtered,
                          aes(x     = ratio,
                              y     = reorder(Description, ratio),
                              size  = Count,
                              color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(regulation = c(up   = "Upregulated",
                                                down = "Downregulated"))) +
  
  scale_color_gradient(low  = "#c0392b",
                       high = "#2980b9",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(
                         barwidth  = 0.8,
                         barheight = 6,
                         ticks     = TRUE)) +
  
  scale_size_continuous(name   = "Gene\ncount",
                        range  = c(3, 9),
                        breaks = pretty(withinKFH_GO_LPS_ATP_df_filtered$Count, n = 4)) +
  
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
withinKFH_GO_LPS_ATP

ggsave("KFH_GO_LPS_ATP_FigureS5B.pdf",
       plot   = withinKFH_GO_LPS_ATP,
       width  = 7,
       height = 8)


####################Healthy Control #########################
###########Semantic analysis GO-BP for downregulated DEGs within LPS stimulated monocytes compared to unstimulated monocytes from Healthy control [FIGURES4C]

simMatrix_withinHC_LPS_down<- calculateSimMatrix(turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023$ID,
                                        orgdb="org.Hs.eg.db",
                                        ont="BP",
                                        method="Rel")
scores_withinHC_LPS_down <- setNames(-log10(turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023$ID)
reducedTerms__withinHC_LPS_down <- reduceSimMatrix(simMatrix_withinHC_LPS_down,
                                                   scores_withinHC_LPS_down,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")
colnames(reducedTerms__withinHC_LPS_down)[1]<-"ID"
parent_term_withinHC_LPS_down<-unique(reducedTerms__withinHC_LPS_down$parentTerm)

write.csv(reducedTerms__withinHC_LPS_down, file = "reducedTerms_withinHC_LPS_down.csv", row.names = FALSE)

withinHC_LPS_down <- turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023[which(turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023$Description %in% parent_term_withinHC_LPS_down),]
write.csv(withinHC_LPS_down, file = "sematic_withinHC_LPS_down.csv", row.names = FALSE)


###########Semantic analysis GO-BP for upregulated DEGs within LPS stimulated monocytes compared to unstimulated monocytes from patients with KFH [FIGURES4B]


simMatrix_withinHC_LPS_up<- calculateSimMatrix(turunen_scRNAseq_within_Control_LPS_GO_upreg_enriched_terms_26072023$ID,
                                      orgdb="org.Hs.eg.db",
                                      ont="BP",
                                      method="Rel")
scores_withinHC_LPS_up <- setNames(-log10(turunen_scRNAseq_within_Control_LPS_GO_upreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_within_Control_LPS_GO_upreg_enriched_terms_26072023$ID)
reducedTerms_withinHC_LPS_up <- reduceSimMatrix(simMatrix_withinHC_LPS_up,
                                                scores_withinHC_LPS_up,
                                       threshold=0.7,
                                       orgdb="org.Hs.eg.db")
colnames(reducedTerms_withinHC_LPS_up)[1]<-"ID"
parent_term_withinHC_LPS_up<-unique(reducedTerms_withinHC_LPS_up$parentTerm)

write.csv(reducedTerms_withinHC_LPS_up, file = "reducedTerms_withinHC_LPS_up.csv", row.names = FALSE)

withinHC_LPS_Up <- turunen_scRNAseq_within_Control_LPS_GO_upreg_enriched_terms_26072023[which(turunen_scRNAseq_within_Control_LPS_GO_upreg_enriched_terms_26072023$Description %in% parent_term_withinHC_LPS_up),]
write.csv(withinHC_LPS_Up, file = "sematic_withinHC_LPS_Up.csv", row.names = FALSE)

#####Plot the GO-terms (upregulated and downregulated in the same plot)

# Add regulation label to each
withinHC_LPS_down$regulation   <- "down"
withinHC_LPS_Up$regulation   <- "up"

withinHC_LPS_Up<-withinHC_LPS_Up[c(-11:-73),]
withinHC_LPS_down<-withinHC_LPS_down[c(-11:-54),]
# Combine
withinHC_GO_LPS_df <- rbind(withinHC_LPS_Up, withinHC_LPS_down)

# Verify
unique(withinHC_GO_LPS_df$regulation)  # should show both "up" and "down"


withinHC_GO_LPS_df_filtered <- withinHC_GO_LPS_df |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = c("up", "down"))  # adjust to your labels
  )

####Dotplot
withinHC_GO_LPS<- ggplot(withinHC_GO_LPS_df_filtered,
                          aes(x     = ratio,
                              y     = reorder(Description, ratio),
                              size  = Count,
                              color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(regulation = c(up   = "Upregulated",
                                                down = "Downregulated"))) +
  
  scale_color_gradient(low  = "#c0392b",
                       high = "#2980b9",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(
                         barwidth  = 0.8,
                         barheight = 6,
                         ticks     = TRUE)) +
  
  scale_size_continuous(name   = "Gene\ncount",
                        range  = c(3, 9),
                        breaks = pretty(withinHC_GO_LPS_df_filtered$Count, n = 4)) +
  
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
withinHC_GO_LPS

ggsave("HC_GO_LPS_FigureS4C.pdf",
       plot   = withinHC_GO_LPS,
       width  = 7,
       height = 8)


#######

###########Semantic analysis GO-BP for downregulated DEGs within LPS stimulated monocytes compared to unstimulated monocytes from patients with HC [FIGURES4B]

simMatrix_withinHC_LPS_ATP_down<- calculateSimMatrix(turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023$ID,
                                            orgdb="org.Hs.eg.db",
                                            ont="BP",
                                            method="Rel")
scores_withinHC_LPS_ATP_down <- setNames(-log10(turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023$ID)
reducedTerms_withinHC_LPS_ATP_down <- reduceSimMatrix(simMatrix_withinHC_LPS_ATP_down,
                                             scores_withinHC_LPS_ATP_down,
                                             threshold=0.7,
                                             orgdb="org.Hs.eg.db")
colnames(reducedTerms_withinHC_LPS_ATP_down)[1]<-"ID"
parent_term_withinHC_LPS_ATP_down<-unique(reducedTerms_withinHC_LPS_ATP_down$parentTerm)

write.csv(reducedTerms_withinHC_LPS_ATP_down, file = "reducedTerms_withinHC_LPS_ATP_down.csv", row.names = FALSE)

withinHC_LPS_ATP_down <- turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023[which(turunen_scRNAseq_within_Control_LPS_ATP_GO_downreg_enriched_terms_26072023$Description %in% parent_term_withinHC_LPS_ATP_down),]
write.csv(withinHC_LPS_ATP_down, file = "sematic_withinHC_LPS_ATP_down.csv", row.names = FALSE)


###########Semantic analysis GO-BP for upregulated DEGs within LPS+ATP stimulated monocytes compared to unstimulated monocytes from patients with HC [FIGURES5B]


simMatrix_withinHC_LPS_ATP_up<- calculateSimMatrix(turunen_scRNAseq_within_Control_LPS_ATP_GO_upreg_enriched_terms_26072023$ID,
                                          orgdb="org.Hs.eg.db",
                                          ont="BP",
                                          method="Rel")
scores_withinHC_LPS_ATP_up <- setNames(-log10(turunen_scRNAseq_within_Control_LPS_ATP_GO_upreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_within_Control_LPS_ATP_GO_upreg_enriched_terms_26072023$ID)
reducedTerms_withinHC_LPS_ATP_up <- reduceSimMatrix(simMatrix_withinHC_LPS_ATP_up,
                                           scores_withinHC_LPS_ATP_up,
                                           threshold=0.7,
                                           orgdb="org.Hs.eg.db")
colnames(reducedTerms_withinHC_LPS_ATP_up)[1]<-"ID"
parent_term_withinHC_LPS_ATP_up<-unique(reducedTerms_withinHC_LPS_ATP_up$parentTerm)

write.csv(reducedTerms_withinHC_LPS_ATP_up, file = "reducedTerms_withinHC_LPS_ATP_up.csv", row.names = FALSE)

withinHC_LPS_ATP_Up <- turunen_scRNAseq_within_Control_LPS_ATP_GO_upreg_enriched_terms_26072023[which(turunen_scRNAseq_within_Control_LPS_ATP_GO_upreg_enriched_terms_26072023$Description %in% parent_term_withinHC_LPS_ATP_up),]
write.csv(withinHC_LPS_ATP_Up, file = "sematic_withinHC_LPS_ATP_Up.csv", row.names = FALSE)

#####Plot the GO-terms (upregulated and downregulated in the same plot)

# Add regulation label to each
withinHC_LPS_ATP_down$regulation   <- "down"
withinHC_LPS_ATP_Up$regulation   <- "up"

withinHC_LPS_ATP_Up<-withinHC_LPS_ATP_Up[c(-11:-88),]
withinHC_LPS_ATP_down<-withinHC_LPS_ATP_down[c(-11:-54),]
# Combine
withinHC_GO_LPS_ATP_df <- rbind(withinHC_LPS_ATP_Up, withinHC_LPS_ATP_down)


withinHC_GO_LPS_ATP_df_filtered <- withinHC_GO_LPS_ATP_df |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = c("up", "down"))  # adjust to your labels
  )

####Dotplot
withinHC_GO_LPS_ATP<- ggplot(withinHC_GO_LPS_ATP_df_filtered,
                              aes(x     = ratio,
                                  y     = reorder(Description, ratio),
                                  size  = Count,
                                  color = p.adjust)) +
  geom_point(alpha = 0.85) +
  
  facet_grid(regulation ~ .,
             scales   = "free_y",
             space    = "free_y",
             labeller = labeller(regulation = c(up   = "Upregulated",
                                                down = "Downregulated"))) +
  
  scale_color_gradient(low  = "#c0392b",
                       high = "#2980b9",
                       name = "Adjusted\np-value",
                       guide = guide_colorbar(
                         barwidth  = 0.8,
                         barheight = 6,
                         ticks     = TRUE)) +
  
  scale_size_continuous(name   = "Gene\ncount",
                        range  = c(3, 9),
                        breaks = pretty(withinHC_GO_LPS_ATP_df_filtered$Count, n = 4)) +
  
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
withinHC_GO_LPS_ATP

ggsave("HC_GO_LPS_ATP_FigureS5D.pdf",
       plot   = withinHC_GO_LPS_ATP,
       width  = 7,
       height = 8)



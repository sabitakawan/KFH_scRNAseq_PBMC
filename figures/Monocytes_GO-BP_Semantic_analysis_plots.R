library(tidyr)
library(dplyr)
library(ggplot2)
library(rrvgo)
library(GOSemSim)
library(stringr)



###########Semantic analysis Untreated_GO-BP for upregulated DEGs in patients with KFH compared to healthy controls [FIGURE6B]

simMatrix_Untreated_down<- calculateSimMatrix(turunen_scRNAseq_Untreated_disease_specific_GO_downreg_enriched_terms_26072023$ID,
                                        orgdb="org.Hs.eg.db",
                                        ont="BP",
                                        method="Rel")
scores_Untreated_down <- setNames(-log10(turunen_scRNAseq_Untreated_disease_specific_GO_downreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_Untreated_disease_specific_GO_downreg_enriched_terms_26072023$ID)
reducedTerms_Untreated_down <- reduceSimMatrix(simMatrix_Untreated_down,
                                         scores_Untreated_down,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")
colnames(reducedTerms_Untreated_down)[1]<-"ID"
parent_term_Untreated_down<-unique(reducedTerms_Untreated_down$parentTerm)

write.csv(reducedTerms_Untreated_down, file = "reducedTerms_Untreated_down.csv", row.names = FALSE)

KFH_Untreated_down <- turunen_scRNAseq_Untreated_disease_specific_GO_downreg_enriched_terms_26072023[which(turunen_scRNAseq_Untreated_disease_specific_GO_downreg_enriched_terms_26072023$Description %in% parent_term_Untreated_down),]
write.csv(KFH_Untreated_down, file = "sematic_KFH_Untreated_down.csv", row.names = FALSE)


###########Untreated_GO-BP for upregulated DEGs in patients with KFH compared to healthy controls

simMatrix_Untreated_up<- calculateSimMatrix(turunen_scRNAseq_Untreated_disease_specific_GO_upreg_enriched_terms_26072023$ID,
                                      orgdb="org.Hs.eg.db",
                                      ont="BP",
                                      method="Rel")
scores_Untreated_up <- setNames(-log10(turunen_scRNAseq_Untreated_disease_specific_GO_upreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_Untreated_disease_specific_GO_upreg_enriched_terms_26072023$ID)
reducedTerms_Untreated_up <- reduceSimMatrix(simMatrix_Untreated_up,
                                       scores_Untreated_up,
                                       threshold=0.7,
                                       orgdb="org.Hs.eg.db")
colnames(reducedTerms_Untreated_up)[1]<-"ID"
parent_term_Untreated_up<-unique(reducedTerms_Untreated_up$parentTerm)

write.csv(reducedTerms_Untreated_up, file = "reducedTerms_Untreated_up.csv", row.names = FALSE)

KFH_Untreated_Up <- turunen_scRNAseq_Untreated_disease_specific_GO_upreg_enriched_terms_26072023[which(turunen_scRNAseq_Untreated_disease_specific_GO_upreg_enriched_terms_26072023$Description %in% parent_term_Untreated_up),]
write.csv(KFH_Untreated_Up, file = "sematic_KFH_Untreated_Up.csv", row.names = FALSE)

#####Plot the GO-terms (upregulated and downregulated in the same plot)

# Add regulation label to each
KFH_Untreated_down$regulation   <- "down"

KFH_Untreated_Up$regulation   <- "up"
KFH_Untreated_Up<-KFH_Untreated_Up[c(-2,-3,-5,-7,-9,-10:-15,-22:-24),]

# Combine
KFH_GO_Untreated_df <- rbind(KFH_Untreated_Up, KFH_Untreated_down)

# Verify
unique(KFH_GO_Untreated_df$regulation)  # should show both "up" and "down"


KFH_GO_Untreated_df_filtered <- KFH_GO_Untreated_df |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = c("up", "down"))  # adjust to your labels
  )

####Dotplot
KFH_GO_Untreated<- ggplot(KFH_GO_Untreated_df_filtered,
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
                        breaks = pretty(KFH_GO_Untreated_df_filtered$Count, n = 4)) +
  
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
KFH_GO_Untreated

ggsave("KFH_GO_Untreated_Figure6B.pdf",
       plot   = KFH_GO_Untreated,
       width  = 7,
       height = 8)


###########Semantic analysis LPS_GO-BP for downregulated DEGs in patients with KFH compared to healthy controls [FIGURE7B]

simMatrix_LPS_down<- calculateSimMatrix(turunen_scRNAseq_LPS_disease_specific_GO_downreg_enriched_terms_26072023$ID,
                                       orgdb="org.Hs.eg.db",
                                       ont="BP",
                                       method="Rel")
scores_LPS_down <- setNames(-log10(turunen_scRNAseq_LPS_disease_specific_GO_downreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_LPS_disease_specific_GO_downreg_enriched_terms_26072023$ID)
reducedTerms_LPS_down <- reduceSimMatrix(simMatrix_LPS_down,
                                        scores_LPS_down,
                                        threshold=0.7,
                                        orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_down)[1]<-"ID"
parent_term_LPS_down<-unique(reducedTerms_LPS_down$parentTerm)

write.csv(reducedTerms_LPS_down, file = "reducedTerms_LPS_down.csv", row.names = FALSE)

KFH_LPS_down <- turunen_scRNAseq_LPS_disease_specific_GO_downreg_enriched_terms_26072023[which(turunen_scRNAseq_LPS_disease_specific_GO_downreg_enriched_terms_26072023$Description %in% parent_term_LPS_down),]
write.csv(KFH_LPS_down, file = "sematic_KFH_LPS_down.csv", row.names = FALSE)


###########LPS_GO-BP for upregulated DEGs in patients with KFH compared to healthy controls

simMatrix_LPS_up<- calculateSimMatrix(turunen_scRNAseq_LPS_disease_specific_GO_upreg_enriched_terms_26072023$ID,
                                     orgdb="org.Hs.eg.db",
                                     ont="BP",
                                     method="Rel")
scores_LPS_up <- setNames(-log10(turunen_scRNAseq_LPS_disease_specific_GO_upreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_LPS_disease_specific_GO_upreg_enriched_terms_26072023$ID)
reducedTerms_LPS_up <- reduceSimMatrix(simMatrix_LPS_up,
                                      scores_LPS_up,
                                      threshold=0.7,
                                      orgdb="org.Hs.eg.db")
colnames(reducedTerms_LPS_up)[1]<-"ID"
parent_term_LPS_up<-unique(reducedTerms_LPS_up$parentTerm)

write.csv(reducedTerms_LPS_up, file = "reducedTerms_LPS_up.csv", row.names = FALSE)

KFH_LPS_Up <- turunen_scRNAseq_LPS_disease_specific_GO_upreg_enriched_terms_26072023[which(turunen_scRNAseq_LPS_disease_specific_GO_upreg_enriched_terms_26072023$Description %in% parent_term_LPS_up),]
write.csv(KFH_LPS_Up, file = "sematic_KFH_LPS_Up.csv", row.names = FALSE)

#####Plot the GO-terms (upregulated and downregulated in the same plot)

# Add regulation label to each
KFH_LPS_Up<-sematic_KFH_LPS_Up
KFH_LPS_down<-sematic_KFH_LPS_down
colnames(KFH_LPS_down)[1]<-"ONTOLOGY"

KFH_LPS_Up$regulation   <- "up"
KFH_LPS_Up<-KFH_LPS_Up[c(-6,-7,-9),]
KFH_LPS_down$regulation   <- "down"
KFH_LPS_down<-KFH_LPS_down[c(-1,-4,-5,-7,-10,-11:-15,-17,-18:-21,-24,-25,-27,-28,-30),]

# Combine
KFH_GO_LPS_df <- rbind(KFH_LPS_Up, KFH_LPS_down)

# Verify
unique(KFH_GO_LPS_df$regulation)  # should show both "up" and "down"


KFH_GO_LPS_df_filtered <- KFH_GO_LPS_df |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = c("up", "down"))  # adjust to your labels
  )

####Dotplot
KFH_GO_LPS<- ggplot(KFH_GO_LPS_df_filtered,
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
                        breaks = pretty(KFH_GO_LPS_df_filtered$Count, n = 4)) +
  
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
KFH_GO_LPS

ggsave("KFH_GO_LPS_Figure7B.pdf",
       plot   = KFH_GO_LPS,
       width  = 7,
       height = 8)


###########Semantic analysis LPS+ATP_GO-BP for upregulated DEGs in patients with KFH compared to healthy controls [Figure8B]

simMatrix_LA_down<- calculateSimMatrix(turunen_scRNAseq_LPS_ATP_disease_specific_GO_downreg_enriched_terms_26072023$ID,
                                                  orgdb="org.Hs.eg.db",
                                                  ont="BP",
                                                  method="Rel")
scores_LA_down <- setNames(-log10(turunen_scRNAseq_LPS_ATP_disease_specific_GO_downreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_LPS_ATP_disease_specific_GO_downreg_enriched_terms_26072023$ID)
reducedTerms_LA_down <- reduceSimMatrix(simMatrix_LA_down,
                                            scores_LA_down,
                                            threshold=0.7,
                                            orgdb="org.Hs.eg.db")
colnames(reducedTerms_LA_down)[1]<-"ID"
parent_term_LA_down<-unique(reducedTerms_LA_down$parentTerm)

write.csv(reducedTerms_LA_down, file = "reducedTerms_LA_down.csv", row.names = FALSE)

KFH_LA_down <- turunen_scRNAseq_LPS_ATP_disease_specific_GO_downreg_enriched_terms_26072023[which(turunen_scRNAseq_LPS_ATP_disease_specific_GO_downreg_enriched_terms_26072023$Description %in% parent_term_LA_down),]
write.csv(KFH_LA_down, file = "sematic_KFH_LA_down.csv", row.names = FALSE)


###########LPS+ATP_GO-BP for upregulated DEGs in patients with KFH compared to healthy controls_FIGURE8B

simMatrix_LA_up<- calculateSimMatrix(turunen_scRNAseq_LPS_ATP_disease_specific_GO_upreg_enriched_terms_26072023$ID,
                                       orgdb="org.Hs.eg.db",
                                       ont="BP",
                                       method="Rel")
scores_LA_up <- setNames(-log10(turunen_scRNAseq_LPS_ATP_disease_specific_GO_upreg_enriched_terms_26072023$qvalue), turunen_scRNAseq_LPS_ATP_disease_specific_GO_upreg_enriched_terms_26072023$ID)
reducedTerms_LA_up <- reduceSimMatrix(simMatrix_LA_up,
                                        scores_LA_up,
                                        threshold=0.7,
                                        orgdb="org.Hs.eg.db")
colnames(reducedTerms_LA_up)[1]<-"ID"
parent_term_LA_up<-unique(reducedTerms_LA_up$parentTerm)

write.csv(reducedTerms_LA_up, file = "reducedTerms_LA_up.csv", row.names = FALSE)

KFH_LA_Up <- turunen_scRNAseq_LPS_ATP_disease_specific_GO_upreg_enriched_terms_26072023[which(turunen_scRNAseq_LPS_ATP_disease_specific_GO_upreg_enriched_terms_26072023$Description %in% parent_term_LA_up),]
write.csv(KFH_LA_Up, file = "sematic_KFH_LA_Up.csv", row.names = FALSE)

#####Plot the GO-terms (upregulated and downregulated in the same plot)

# Add regulation label to each
KFH_LA_down <-sematic_KFH_LA_down
KFH_LA_Up <- sematic_KFH_LA_Up


KFH_LA_down$regulation   <- "down"
KFH_LA_down<-KFH_LA_down[c(-1,-7,-8,-10,-16,-18),]
KFH_LA_Up$regulation   <- "up"
KFH_LA_Up<-KFH_LA_Up[c(-4,-6,-11,-14,-16,-18,-19:-32),]

# Combine
KFH_GO_LA_df <- rbind(KFH_LA_Up, KFH_LA_down)

# Verify
unique(KFH_GO_LA_df$regulation)  # should show both "up" and "down"


KFH_GO_df_filtered <- KFH_GO_LA_df |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = c("up", "down"))  # adjust to your labels
  )

####Dotplot KFH_LPS+ATP GO terms
KFH_GO_LA<- ggplot(KFH_GO_df_filtered,
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
                        breaks = pretty(KFH_GO_df_filtered$Count, n = 4)) +
  
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
KFH_GO_LA

ggsave("KFH_GO_LA_Figure8B.pdf",
       plot   = KFH_GO_LA,
       width  = 7,
       height = 8)


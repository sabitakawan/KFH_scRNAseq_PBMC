library(stringr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(rrvgo)
library(GOSemSim)
library(clusterProfiler)
library("org.Hs.eg.db")
library(ggh4x)
OrgDb <- org.Hs.eg.db
Symbol_col <- 'SYMBOL'
padj.thr <- 1e-03

####Marker genes and GO analysis for T cells
seurat_default_options <- list(
  Seurat.memsafe = TRUE,
  Seurat.warn.umap.uwot = TRUE,
  Seurat.checkdots = "warn",
  Seurat.limma.wilcox.msg = TRUE,
  Seurat.Rfast2.msg = TRUE,
  Seurat.warn.vlnplot.split = TRUE
)

process_ego <- function(ego, max.char = 50, stitle) { # ego = ego_do; max.char=50; stitle = "Downregulated"
  ego@result$ratio <- sapply(ego@result$GeneRatio, FUN = function(x) {
    splt <- as.numeric(unlist(strsplit(x, split = '[/]')))
    return(splt[1]/splt[2])
  })
  ego@result$Count <- sapply(ego@result$GeneRatio, FUN = function(x) {
    splt <- as.integer(unlist(strsplit(x, split = '[/]')))
    return(splt[1])
  })
  to.plt <- ego@result %>%
    slice_max(n = 10, order_by = -(p.adjust), with_ties = FALSE) %>%
    mutate(term = paste0(ID,': ', Description)) %>%
    mutate(term = ifelse(nchar(term) > max.char, paste0(substr(term, 1, max.char),'...'), term))
  # Plot
  # Create unique levels for each term within ontology
  to.plt$unique_term <- paste(to.plt$ONTOLOGY, to.plt$term, sep = ": ")
  p <- to.plt %>%
    ggplot(aes(x = ratio, y = reorder(unique_term, ratio), size = factor(x = Count, levels = min(Count):max(Count)), fill = p.adjust)) +
    geom_point(shape = 21) +
    scale_fill_distiller(palette = 'PuBuGn', direction = -1) +
    scale_y_discrete(labels = function(x) gsub(".*? (GO:\\d+): ", "\\1 ", x)) +  # remove ontology prefix
    facet_nested(ONTOLOGY ~ ., scales = "free_y", space="free_y") +
    xlab('Gene ratio') +
    ylab('Ontology term') +
    labs(subtitle = stitle, size="Count") +
    guides(size = guide_legend(order = 1), fill = guide_legend(order = 2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(list(ego, p))
}

xm <- readRDS("/projappl/project_2011471/Seurat_object_harmony_annotated.RDS")
xm$broad_annot <- xm$cluster_harmony
levels(xm$broad_annot) <- c("T", "NK", "T", "T", "T", "Myeloid", "T", "T", "B", "T", 
                            "T", "B", "T", "T", "T", "NK", "T? Metamyelocyte?", "T", 
                            "B", "NK", "T", "T", "Myeloid", "T?", "T?", "Megakaryocyte")
levels(xm$broad_annot)[grep("\\?", levels(xm$broad_annot))] <- "?"

xm$narrow_annot <- xm$cluster_harmony
levels(xm$narrow_annot) <- c("CD4+ Central Memory", "NK", "CD8+", "CD4+ Naive", "CD4+ Naive", 
                             "Monocytes, DCs, Macrophages", "CD4+ Naive", "CD4+ Naive", "B", 
                             "CD4+ Naive", "CD4+ Effector Memory", "B", "CD8+", "CD4+ Naive", 
                             "CD4+ Naive", "NK", "CD8+?", "CD4+ Naive", "B", "NK", "CD4+ Naive", 
                             "CD4+ Central Memory", "Neutrophil", "CD4+?", "CD4+ Central Memory?", 
                             "Megakaryocyte")
levels(xm$narrow_annot)[grep("\\?", levels(xm$narrow_annot))] <- "?"

for (ct in c("B", "T","NK")) {  ct = "B"; ct = "T"; ct = "NK"
if (!dir.exists(paste0("/projappl/project_2011471/ct_diffexp/", ct))) {
  dir.create(paste0("/projappl/project_2011471/ct_diffexp/", ct), recursive = TRUE)
}
}

so <- subset(x = xm, subset = broad_annot==ct) # so <- subset(x = xm, subset = narrow_annot=="Neutrophil")
#openxlsx::write.xlsx(x=table(so$disease_status, so$condition), file = paste0("/workstation/project_results/turunen_project1/ct_diffexp/", ct, "/turunen_scRNAseq_how_many_cells_each_26072023.xlsx"))

# Now, make the contrasts between patient and control cells, within each condition
table_of_DEGs <- data.frame(contrast=character(0), base_filter=numeric(0), DEG=numeric(0), upreg=numeric(0), downreg=numeric(0), strong_upreg=numeric(0), strong_downreg=numeric(0))
for (cond in unique(so$condition)) {cond = unique(so$condition)[1]; cond = "LPS-ATP"
path_of_interest <- paste0("/projappl/project_2011471/ct_diffexp//", ct, "/withinCondition_diseaseStatusComparison/turunen_scRNAseq_", cond, "_disease-specific_DEG")
if (!dir.exists(dirname(path_of_interest))) {
  dir.create(dirname(path_of_interest))
}
}
cat(cond, "\n")
sso <- subset(so, subset = condition == cond)
sso <- SetIdent(sso, value = "disease_status")
marker_genes <- FindMarkers(sso, ident.1 = "Patient", ident.2 = "Control")
marker_genes$significant <- marker_genes$p_val_adj<0.05
marker_genes$gene <- rownames(marker_genes)
marker_genes <- marker_genes %>% mutate(category = case_when(avg_log2FC > 0 & p_val_adj < 0.05 ~ "upregulated",
                                                             avg_log2FC < 0 & p_val_adj < 0.05 ~ "downregulated",
                                                             TRUE ~ "below threshold"))
marker_genes$NFKB_target_genes <- marker_genes$gene %in% c("MYD88", "TICAM1", "FADD")
marker_genes$NLRP3_inflammasome_target_genes <- marker_genes$gene %in% c("CASP1", "NLRP3", "PYCARD")
so$samp <- substr(so$sample, 1, 2)
for (k_elim in 1:3) { # k_elim=1
  k_eli <- paste0("P", k_elim) # This will hold the sample name
  k_el <- paste0("without_P", k_elim) # This will hold the column name
  print(k_eli)
  sso <- subset(so, subset = condition == cond & samp != k_eli)
  sso <- SetIdent(sso, value = "disease_status")
  marker_genes2 <- FindMarkers(sso, ident.1 = "Patient", ident.2 = "Control")
  marker_genes2$significant <- marker_genes2$p_val_adj<0.05
  marker_genes2$gene <- rownames(marker_genes2)
  marker_genes2 <- marker_genes2 %>% mutate(category = case_when(avg_log2FC > 0 & p_val_adj < 0.05 ~ "upregulated",
                                                                 avg_log2FC < 0 & p_val_adj < 0.05 ~ "downregulated",
                                                                 TRUE ~ "below threshold"))
  marker_genes[,k_el] <- "below threshold"
  marker_genes[,k_el] <- marker_genes2$category[match(marker_genes$gene, marker_genes2$gene)]
  marker_genes[is.na(marker_genes[,k_el]),k_el] <- "below threshold"
}

# Compute unity score
for (k_elim in 1:3) { # k_elim=1
  k_eli <- paste0("P", k_elim) # This will hold the sample name
  k_el <- paste0("without_P", k_elim) # This will hold the column name
  marker_genes[,paste0(k_el, "_match")] <- as.integer(marker_genes$category == marker_genes[,k_el])
}
marker_genes$unity_score <- marker_genes[,"without_P1_match"] + marker_genes[,"without_P2_match"] + marker_genes[,"without_P3_match"]
marker_genes[,"without_P1_match"] <- marker_genes[,"without_P2_match"] <- marker_genes[,"without_P3_match"] <- NULL
marker_genes <- marker_genes %>% dplyr::select(dplyr::everything(), "without_P1", "without_P2", "without_P3")

table_of_DEGs[nrow(table_of_DEGs)+1,] <- list(
  paste0("within_", cond, "_Patient_Vs_Control"),
  nrow(marker_genes),
  sum(marker_genes$p_val_adj < 0.05),
  sum(marker_genes$category == "upregulated"),
  sum(marker_genes$category == "downregulated"),
  sum(marker_genes$p_val_adj < 0.05 & marker_genes$avg_log2FC > log2(1.5)),
  sum(marker_genes$p_val_adj < 0.05 & marker_genes$avg_log2FC < -log2(1.5))
)
#save markers selecting also for Unity score
openxlsx::write.xlsx(x = marker_genes, file = paste0(path_of_interest, "s_26072023.xlsx"))
if (nrow(
  marker_genes[abs(marker_genes$avg_log2FC)>log2(1.5) & marker_genes$p_val_adj < 0.05 & marker_genes$unity_score == 3,]
)>0) {
  openxlsx::write.xlsx(x = marker_genes[abs(marker_genes$avg_log2FC)>log2(1.5) & marker_genes$p_val_adj < 0.05 & marker_genes$unity_score == 3,],
                       file = paste0(path_of_interest, "s_strong_26072023.xlsx"))
}


###########Untreated_GO-BP for upregulated DEGs in patients with KFH compared to healthy controls [Figure 6D]

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

#####Plot the GO-terms

# Add regulation label to each
KFH_Untreated_Up$regulation   <- "up"
KFH_Untreated_Up<-KFH_Untreated_Up[c(-4),]

KFH_GO_Untreated_df <- KFH_Untreated_Up |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = "up")  # adjust to your labels
  )

####Dotplot
KFH_GO_Untreated<- ggplot(KFH_GO_Untreated_df,
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
                        breaks = pretty(KFH_GO_Untreated_df$Count, n = 4)) +
  
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

ggsave("KFH_GO_Untreated_Figure6D.pdf",
       plot   = KFH_GO_Untreated,
       width  = 7,
       height = 8)


###########Semantic analysis LPS_GO-BP for upregulated DEGs in patients with KFH compared to healthy controls [FIGURE7D]

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

####### Semantic analysis LPS_GO-BP for downregulated DEGs in patients with KFH compared to healthy controls [FIGURE7B]
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


#####Plot the GO-terms (upregulated and downregulated in the same plot)

# Add regulation label to each
KFH_LPS_Up$regulation   <- "up"
KFH_LPS_Up<-KFH_LPS_Up[c(-8),]
KFH_LPS_down$regulation   <- "down"
KFH_LPS_down<-KFH_LPS_down[c(-2,-3,-5,-6,-7,-8, -9,-11,-12,-14),]

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

ggsave("KFH_GO_LPS_Figure7D.pdf",
       plot   = KFH_GO_LPS,
       width  = 7,
       height = 8)


###########Semantic analysis LPS+ATP_GO-BP for upregulated DEGs in patients with KFH compared to healthy controls [Figure8D]

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

KFH_LA_Up$regulation   <- "up"
KFH_LA_Up<-KFH_LA_Up[c(-2),]

KFH_GO_LA_df<- KFH_LA_Up |>
  filter(p.adjust < 0.05, ONTOLOGY=="BP") |>
  group_by(regulation) |>
  slice_min(p.adjust, n = 20) |>          # top 15 per direction — adjust as needed
  ungroup() |>
  mutate(
    Description = str_wrap(Description, width = 45),
    regulation  = factor(regulation, levels = "up")  # adjust to your labels
  )

####Dotplot KFH_LPS+ATP GO terms
KFH_GO_LA<- ggplot(KFH_GO_LA_df,
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
                        breaks = pretty(KFH_GO_LA_df$Count, n = 4)) +
  
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

ggsave("KFH_GO_LA_Figure8D.pdf",
       plot   = KFH_GO_LA,
       width  = 7,
       height = 8)




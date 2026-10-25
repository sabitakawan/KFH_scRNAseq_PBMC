

## run with R version 4.1.3, Seurat Version 4.0.3
## other packages

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggh4x)

OrgDb <- org.Hs.eg.db
Symbol_col <- 'SYMBOL'
padj.thr <- 1e-03

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

xm <- readRDS("saved/ZenodoObjs/Seurat_object_harmony_annotated.RDS")
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


ct <- "Monocytes, DCs, Macrophages"  #(narrow_annot)

so <- subset(x = xm, subset = narrow_annot==ct)



## within condition

pathToFolder2 <- "diffExpResults/withinCondition_diseaseStatusComparison/monocytes/"

for (cond in unique(so$disease_status)) { 

  cat(cond, "\n")
  sso <- subset(so, subset = disease_status == cond)
  
  comparisons = list(c("LPS","Untreated"),
                     c("LPS-ATP","Untreated"),
                     c("LPS-ATP","LPS"))
  
  sapply(comparisons, function(y){
    
    
    sso_sub <- subset(sso, subset = condition %in% y)
    sso_sub <- SetIdent(sso_sub, value = "condition")
    
    marker_genes <- FindMarkers(
      sso_sub,
      ident.1 = y[1],
      ident.2 = y[2],
      logfc.threshold = 0,
      min.pct = 0,
      min.diff.pct = -Inf
    )
    
    marker_genes$significant <- marker_genes$p_val_adj<0.05
    marker_genes$gene <- rownames(marker_genes)
    marker_genes <- marker_genes %>% mutate(category = case_when(avg_log2FC > 0 & p_val_adj < 0.05 ~ "upregulated",
                                                                 avg_log2FC < 0 & p_val_adj < 0.05 ~ "downregulated",
                                                                 TRUE ~ "below threshold"))
    
    write.table(marker_genes, file=paste0(pathToFolder2,cond,"_",y[1],"_",y[2],"_DEGs.tsv"),
                row.names=F,
                col.names=T,
                sep="\t",
                quote=F)
    

  }, simplify=F)
  
  
}
  
  


## within treatment

pathToDirectory <- "diffExpResults/withinTreatment_diseaseStatusComparison/monocytes/"


for (cond in unique(so$condition)) { 
  path_of_interest <- paste0(pathToDirectory, cond, "_", gsub(", ","_", ct), "_disease-specific_DEG")

  cat(cond, "\n")
  sso <- subset(so, subset = condition == cond)
  sso <- SetIdent(sso, value = "disease_status")
  
  marker_genes <- FindMarkers(
    sso,
    ident.1 = "Patient",
    ident.2 = "Control",
    logfc.threshold = 0,
    min.pct = 0,
    min.diff.pct = -Inf
  )
  
  marker_genes$significant <- marker_genes$p_val_adj<0.05
  marker_genes$gene <- rownames(marker_genes)
  marker_genes <- marker_genes %>% mutate(category = case_when(avg_log2FC > 0 & p_val_adj < 0.05 ~ "upregulated",
                                                               avg_log2FC < 0 & p_val_adj < 0.05 ~ "downregulated",
                                                               TRUE ~ "below threshold"))
  marker_genes$NFKB_target_genes <- marker_genes$gene %in% c("MYD88", "TICAM1", "FADD")
  marker_genes$NLRP3_inflammasome_target_genes <- marker_genes$gene %in% c("CASP1", "NLRP3", "PYCARD")
  
  write.table(marker_genes, file=paste0(path_of_interest,"/all_results.txt"))
  
  so$samp <- substr(so$sample, 1, 2)
  for (k_elim in 1:3) { # k_elim=1
    k_eli <- paste0("P", k_elim) # This will hold the sample name
    k_el <- paste0("without_P", k_elim) # This will hold the column name
    print(k_eli)
    sso <- subset(so, subset = condition == cond & samp != k_eli)
    sso <- SetIdent(sso, value = "disease_status")
    marker_genes2 <- FindMarkers(sso, ident.1 = "Patient", ident.2 = "Control",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 min.diff.pct = -Inf)
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
  
  write.table(marker_genes, file=paste0(path_of_interest, "s_29062026.xlsx"))

  if (nrow(
    marker_genes[abs(marker_genes$avg_log2FC)>log2(1.5) & marker_genes$p_val_adj < 0.05 & marker_genes$unity_score == 3,]
  )>0) {
    write.table(marker_genes[abs(marker_genes$avg_log2FC)>log2(1.5) & marker_genes$p_val_adj < 0.05 & marker_genes$unity_score == 3,],
                file=paste0(path_of_interest, "s_strong_29062026.xlsx"))
  }
  
}






























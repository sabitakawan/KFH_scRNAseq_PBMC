

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


ct <- "T"                            #(broad_annot)

so <- subset(x = xm, subset = broad_annot==ct)



## within condition
pathToFolder2 <- "diffExpResults/withinCondition_diseaseStatusComparison/tcells/"


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
pathToDirectory <- "diffExpResults/withinTreatment_diseaseStatusComparison/tcells/"


for (cond in unique(so$condition)) { 
  path_of_interest <- paste0(pathToDirectory, cond, "_tcells_disease-specific_DEG")

  cat(cond, "\n")
  sso <- subset(so, subset = condition == cond)
  sso <- SetIdent(sso, value = "disease_status")
  
  totalcells <-  dim(sso)[2]
  
  set.seed(123)
  sampled <- sample(1:totalcells, size = 0.25 * length(1:totalcells))
  
  sso_sub <- sso[,sampled]

  marker_genes <- FindMarkers(
    sso_sub,
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
  
  
  write.table(marker_genes, file=paste0(path_of_interest, "_19072026.xlsx"))

  
}






























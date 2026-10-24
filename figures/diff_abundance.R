library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(speckle)
library(limma)
library(tidyr)
library(dplyr)
library(edgeR)

obj <- readRDS("saved/ZenodoObjs/Seurat_object_harmony_annotated.RDS")
obj$sampleDonor <- sapply(strsplit(obj$orig.ident,"_"), function(x) x[2])

##Strategy: glmQLFit  (Aaron Lun)
#https://bioconductor.org/books/3.16/OSCA.multisample/differential-abundance.html


matrix_cts_broad <- as.data.frame(table(obj$sampleDonor, obj$condition, obj$broad_annot))
colnames(matrix_cts_broad) <- c("Sample", "Treatment", "CellType", "Counts")
matrix_cts_broad$CellType <- as.character(matrix_cts_broad$CellType)
matrix_cts_broad$CellType[matrix_cts_broad$CellType == "?"] <- "Unknown"
matrix_cts_broad$Condition <- ifelse(grepl("C", matrix_cts_broad$Sample), "Control", "Patient")

counts_wide_broad <- matrix_cts_broad %>%
  pivot_wider(id_cols = c(Sample, Treatment, Condition),
              names_from = CellType, values_from = Counts, values_fill = 0)

cell_cols_broad <- setdiff(colnames(counts_wide_broad), c("Sample", "Treatment", "Condition"))


matrix_cts <- as.data.frame(table(obj$sampleDonor, obj$condition, obj$narrow_annot))
colnames(matrix_cts) <- c("Sample", "Treatment", "CellType", "Counts")
matrix_cts$CellType <- as.character(matrix_cts$CellType)
matrix_cts$CellType[matrix_cts$CellType == "?"] <- "Unknown"
matrix_cts$Condition <- ifelse(grepl("C", matrix_cts$Sample), "Control", "Patient")


counts_wide <- matrix_cts %>%
  pivot_wider(id_cols = c(Sample, Treatment, Condition),
              names_from = CellType, values_from = Counts, values_fill = 0)

cell_cols <- setdiff(colnames(counts_wide), c("Sample", "Treatment", "Condition"))





run_edgeR_DA <- function(sub_counts, cell_cols) {
  
  # celltype x sample count matrix
  count_mat <- t(as.matrix(sub_counts[, cell_cols]))
  colnames(count_mat) <- as.character(sub_counts$Sample)
  storage.mode(count_mat) <- "integer"
  
  meta <- sub_counts %>% select(Sample, Condition)
  meta$Condition <- factor(meta$Condition, levels = c("Control", "Patient"))
  
  y.ab <- DGEList(counts = count_mat, samples = meta)
  
  keep <- filterByExpr(y.ab, group = y.ab$samples$Condition)
  cat("Cell types kept:", sum(keep), "of", length(keep), "\n")
  print(names(keep)[!keep])  # show what got dropped
  
  y.ab <- y.ab[keep, , keep.lib.sizes = FALSE]
  
  y.ab <- calcNormFactors(y.ab)  # TMM normalization, handles depth here
  
  design <- model.matrix(~Condition, data = y.ab$samples)
  
  y.ab <- estimateDisp(y.ab, design, trend = "none")
  
  fit.ab <- glmQLFit(y.ab, design, robust = TRUE)
  res <- glmQLFTest(fit.ab, coef = "ConditionPatient")
  
  resdf <- as.data.frame(topTags(res, n = Inf))
  resdf$CellType <- rownames(resdf)
  rownames(resdf) <- NULL
  
  resdf
}

treatments <- unique(counts_wide$Treatment)

edgeR_per_treatment <- lapply(treatments, function(trt) {
  cat("\n===", trt, "===\n")
  sub <- counts_wide %>% filter(Treatment == trt)
  res <- run_edgeR_DA(sub, cell_cols)
  res$Treatment <- trt
  res
})
names(edgeR_per_treatment) <- treatments

edgeR_treatment_df <- bind_rows(edgeR_per_treatment)


## ---- 3. global test (sum counts across treatments per donor) ----
global_counts <- counts_wide %>%
  group_by(Sample, Condition) %>%
  summarise(across(all_of(cell_cols), sum), .groups = "drop")

cat("\n=== Global ===\n")
edgeR_global_df <- run_edgeR_DA(global_counts, cell_cols)
edgeR_global_df$Treatment <- "Global"

## ---- 4. combine everything ----
edgeR_all <- bind_rows(edgeR_global_df, edgeR_treatment_df) %>%
  select(CellType, Treatment, logFC, logCPM, F, PValue, FDR) %>%
  arrange(Treatment, PValue)

edgeR_all


edgeR_all <- edgeR_all %>%
  mutate(sig_label = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01  ~ "**",
    FDR < 0.05  ~ "*",
    FDR < 0.1   ~ ".",
    TRUE ~ "ns"
  ))

# edgeR_all treatments logically, order cell types by global estimate for readability
edgeR_all$Treatment <- factor(edgeR_all$Treatment, 
                            levels = c("Global", "Untreated", "LPS", "LPS-ATP"))


da_test_edgeR <- ggplot(edgeR_all, aes(x = Treatment, y = CellType, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), color = "grey25", size = 6, vjust = 0.7) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name = "KFH effect\n(log2FC, edgeR)",
                       guide = guide_colorbar(barheight = unit(6, "cm"),   # <- taller bar
                                              barwidth  = unit(0.75, "cm"))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, size=14, hjust=1, vjust=0.5, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        legend.text = element_text(size = 13), 
        plot.title=element_text(size=16, hjust=0.5))+
  ggtitle("Differential abundance (KFH vs controls)")


## Main Figure 2D

pdf(file="figs/mainFigures/Figure2D.pdf", width=8, height =8)
plot(da_test_edgeR)
dev.off()



## broad


treatments <- unique(counts_wide_broad$Treatment)

edgeR_per_treatment_broad <- lapply(treatments, function(trt) {
  cat("\n===", trt, "===\n")
  sub <- counts_wide_broad %>% filter(Treatment == trt)
  res <- run_edgeR_DA(sub, cell_cols_broad)
  res$Treatment <- trt
  res
})
names(edgeR_per_treatment_broad) <- treatments

edgeR_treatment_df <- bind_rows(edgeR_per_treatment_broad)


## ---- 3. global test (sum counts across treatments per donor) ----
global_counts <- counts_wide_broad %>%
  group_by(Sample, Condition) %>%
  summarise(across(all_of(cell_cols_broad), sum), .groups = "drop")

cat("\n=== Global ===\n")
edgeR_global_df <- run_edgeR_DA(global_counts, cell_cols_broad)
edgeR_global_df$Treatment <- "Global"

## ---- 4. combine everything ----
edgeR_all <- bind_rows(edgeR_global_df, edgeR_treatment_df) %>%
  select(CellType, Treatment, logFC, logCPM, F, PValue, FDR) %>%
  arrange(Treatment, PValue)

edgeR_all


edgeR_all <- edgeR_all %>%
  mutate(sig_label = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01  ~ "**",
    FDR < 0.05  ~ "*",
    FDR < 0.1   ~ ".",
    TRUE ~ "ns"
  ))

# edgeR_all treatments logically, order cell types by global estimate for readability
edgeR_all$Treatment <- factor(edgeR_all$Treatment, 
                              levels = c("Global", "Untreated", "LPS", "LPS-ATP"))


da_test_edgeR_2 <- ggplot(edgeR_all, aes(x = Treatment, y = CellType, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), color = "grey25", size = 6, vjust = 0.7) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, name = "KFH effect\n(log2FC, edgeR)",
                       guide = guide_colorbar(barheight = unit(6, "cm"),   # <- taller bar
                                              barwidth  = unit(0.75, "cm"))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, size=14, hjust=1, vjust=0.5, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        legend.text = element_text(size = 13), 
        plot.title=element_text(size=16, hjust=0.5))+
  ggtitle("Differential abundance (KFH vs controls)")


## Figure S1C ##

pdf(file="figs/suppFigures/FigureS1C.pdf", width=8, height =8)
plot(da_test_edgeR_2)
dev.off()












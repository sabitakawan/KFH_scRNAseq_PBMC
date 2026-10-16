library(RColorBrewer)
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(KEGGREST)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(Matrix)
library(circlize)
library(tibble)
library(tidyr)

####Load seurat object
Seurat_object_harmony_annotated <- readRDS("/scratch/project_2011471/Seurat_object_harmony_annotated.RDS")
Seurat_object_harmony_annotated[["harmony"]]@assay.used <- "SCT"


####UMAP
celltype_colors <- c(
  "B" = "#63B3E1",
  "CD4+ Central Memory" = "#D160FF",
  "CD4+ Effector Memory" = "#FFA2B0",
  "CD4+ Naive" = "#FF87F7",
  "CD8+" = "#FF9DA3",
  "Megakaryocyte" = "#2961B1",
  "Monocytes, DCs, Macrophages" = "#81e686",
  "Neutrophil" = "#42584B",
  "NK" = "#1D8E32",
  "?" = "#D90021"
)
UMAP_narrow_annot<-DimPlot(Seurat_object_harmony_annotated, group.by = "narrow_annot", pt.size = 0.5, cols=celltype_colors)
UMAP_narrow_annot


########
celltype_colors_broad <- c(
  "B" = "#63B3E1",
  "T" = "#FF9DA3",
  "Megakaryocyte" = "#2961B1",
  "Myeloid" = "#81e686",
  "Neutrophil" = "#42584B",
  "NK" = "#1D8E32",
  "?" = "#D90021"
)
UMAP_broad_annot<-DimPlot(Seurat_object_harmony_annotated, group.by = "broad_annot", pt.size = 0.5, cols=celltype_colors_broad)
UMAP_broad_annot

ggsave("UMAP_narrow_annotFigure2B.pdf", UMAP_narrow_annot, width=8, height=6, dpi=300)
ggsave("UMAP_broad_annot.pdf_FigureS1A", UMAP_narrow_annot, width=8, height=6, dpi=300)


####NLRP3 expression

UMAP_NLRP3<-FeaturePlot(
  Seurat_object_harmony_annotated,
  features = "NLRP3",
  reduction = "umap",
  max.cutoff = "q95", 
  cols = c("grey","darkred"),  # low to high expression
  pt.size = 0.5
)

UMAP_NLRP3

####IL1B expression

UMAP_IL1B<-FeaturePlot(
  Seurat_object_harmony_annotated,
  features = "IL1B",
  reduction = "umap",
  max.cutoff = "q95", 
  cols = c("grey","darkred"),  # low to high expression
  pt.size = 0.5
)

UMAP_IL1B

ggsave("UMAP_IL1B_Figure4B.pdf", UMAP_IL1B, width=8, height=6, dpi=300)
ggsave("UMAP_NLRP3_Figure4A.pdf", UMAP_IL1B, width=8, height=6, dpi=300)


#####Cell number
# Check metadata columns
colnames(Seurat_object_harmony_annotated@meta.data)

# Count retained cells per sample x condition
df <- Seurat_object_harmony_annotated@meta.data %>%
  count(sample, disease_status, condition, name = "retained_cells")

# Order: group -> condition -> sample
c

cell_number<-ggplot(df, aes(x = sample, y = retained_cells)) +
  geom_text(aes(label = retained_cells), vjust = -0.5, size = 4) +
  geom_col(width = 0.9) +
  labs(x = NULL, y = "Retained cell number", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
cell_number

ggsave("Number of cells_Figure2A.pdf", cell_number, width=8, height=6, dpi=300)


#####cell proportion

# Count cells per sample x cell type [Broad_annotation Figure2C]
cell_proportion_df <- Seurat_object_harmony_annotated@meta.data %>%
  count(sample, disease_status, condition, narrow_annot, name = "cellnumber")

# Sample order: group -> condition -> sample
sample_order <- cell_proportion_df %>%
  distinct(sample, disease_status, condition) %>%
  mutate(
    group     = factor(disease_status, levels = c("Control", "Patient")),
    condition = factor(condition, levels = c("Untreated", "LPS", "LPS-ATP"))
  ) %>%
  arrange(group, condition, sample) %>%
  pull(sample)


# Cell type order (top of stack -> bottom, same as legend)
ct_levels <- c("B",
               "CD4+ Central Memory",
               "CD4+ Effector Memory",
               "CD4+ Naive",
               "CD8+",
               "Megakaryocyte",
               "Monocytes, DCs, Macrophages",
               "Neutrophil",
               "NK",
               "?")

# Matched colors (RColorBrewer "Paired", as in your image)
ct_colors <- c("B"                                     = "#A6CEE3",
               "CD4+ Central Memory"                   = "#1F78B4",
               "CD4+ Effector Memory"                  = "#B2DF8A",
               "CD4+ Naive"                            = "#33A02C",
               "CD8+"                                  = "#FB9A99",
               "Megakaryocyte"                         = "#E31A1C",
               "Monocytes, DCs, Macrophages"           = "#FDBF6F",
               "Neutrophil"                            = "#FF7F00",
               "NK"                                    = "#CAB2D6",
               "?"                                     = "#6A3D9A")


cell_proportion_df <- cell_proportion_df %>% 
  mutate(sample = factor(sample, levels = sample_order),
         CellTypes  = factor(narrow_annot, levels = ct_levels))


cell_proportion_narrow_annot<-ggplot(cell_proportion_df, aes(x = sample, y = cellnumber, fill = CellTypes)) +
  geom_col(position = "fill", width = 1) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = ct_colors, drop = FALSE) +
  labs(x = NULL, y = "Percentage", fill = "CellTypes") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cell_proportion_narrow_annot


# Count cells per sample x cell type [Broad_annotation Figure S1B]
cell_proportion_broad_df <- Seurat_object_harmony_annotated@meta.data %>%
  count(sample, disease_status, condition, broad_annot, name = "cellnumber")

# Sample order: group -> condition -> sample
broad_sample_order <- cell_proportion_broad_df %>%
  distinct(sample, disease_status, condition) %>%
  mutate(
    group     = factor(disease_status, levels = c("Control", "Patient")),
    condition = factor(condition, levels = c("Untreated", "LPS", "LPS-ATP"))
  ) %>%
  arrange(group, condition, sample) %>%
  pull(sample)


ct_levels_broad <- c("B",
                     "T",
                     "Megakaryocyte",
                     "Myeloid",
                     "Neutrophil",
                     "NK",
                     "?")

# Matched colors (RColorBrewer "Paired")
ct_colors_broad <- c("B"                               = "#A6CEE3",
                     "Myeloid"                               = "#B2DF8A",
                     "NK"                                    = "#33A02C",
                     "T"                                     = "#FB9A99",
                     "?"                                     = "#E31A1C",
                     "Megakaryocyte"                         = "#1F78B4")


cell_proportion_broad_df <- cell_proportion_broad_df %>% 
  mutate(sample = factor(sample, levels = sample_order),
         CellTypes  = factor(broad_annot, levels = ct_levels_broad))


cell_proportion_broad_annot<-ggplot(cell_proportion_broad_df, aes(x = sample, y = cellnumber, fill = CellTypes)) +
  geom_col(position = "fill", width = 1) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = ct_colors_broad, drop = FALSE) +
  labs(x = NULL, y = "Percentage", fill = "CellTypes") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cell_proportion_broad_annot


ggsave("cell_proportion_narrow_annot_Figure2C.pdf", cell_proportion_narrow_annot, width=8, height=6, dpi=300)
ggsave("cell_proportion_broad_annot_FigureS1B.pdf", cell_proportion_broad_annot, width=8, height=6, dpi=300)

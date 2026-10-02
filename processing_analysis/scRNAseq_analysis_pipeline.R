#!/usr/bin/env Rscript
# ==============================================================================
# scRNA-seq processing and integration pipeline
#
# Overview
# --------
# This script implements a two-stage scRNA-seq analysis pipeline:
#
#   Stage 1 (per sample)   Ambient-RNA correction (SoupX), doublet detection
#                           (DoubletFinder), and standard per-sample QC /
#                           clustering / cell-type scoring (SingleR).
#
#   Stage 2 (integration)  Merging of all per-sample objects, QC filtering,
#                           SCTransform normalization, batch integration with
#                           Harmony, cluster annotation (SingleR + manual
#                           curation from marker genes), and a per-cell-type
#                           differential expression analysis across conditions.
#
# Stage 1 is designed to be run once per sample. 
# Stage 2 is run once, after every sample has been processed by Stage 1.
#
# Usage
# -----
#   # Stage 1: process a single sample (1-based index into `SAMPLE_DIRS`)
#   Rscript scRNAseq_analysis_pipeline.R preprocess <sample_index>
#
#   # Stage 2: integrate all processed samples
#   Rscript scRNAseq_analysis_pipeline.R integrate
#
# All paths and analysis parameters are collected in the "Configuration"
# section below; nothing else in the script should need to be edited to
# rerun the pipeline on a new machine or a new set of samples.
# ==============================================================================

## =============================================================================
## 1. Configuration
## =============================================================================
# Edit this block only. Everything below derives from these values.

# --- Directory layout ---------------------------------------------------------
BASE_DIR        <- "path_to_dir/KFH_scRNAseq/"
CELLRANGER_DIR  <- file.path(BASE_DIR, "cellranger")            # one subfolder per sample, each with outs/{raw,filtered}_feature_bc_matrix
SAMPLE_RDS_DIR  <- file.path(BASE_DIR, "sample_seurats")        # Stage 1 output: one .rds per sample
CHECKPOINT_DIR  <- file.path(BASE_DIR, "checkpoints")           # Stage 2 intermediate objects (resumable)
PLOTS_DIR       <- file.path(BASE_DIR, "plots")                 # Stage 2 QC/diagnostic figures
MARKERS_DIR     <- file.path(BASE_DIR, "cluster_markers")       # per-cluster marker gene tables
DE_DIR          <- file.path(BASE_DIR, "differential_expression") # per-cell-type DE tables

for (d in c(SAMPLE_RDS_DIR, CHECKPOINT_DIR, PLOTS_DIR, MARKERS_DIR, DE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Sample sheet --------------------------------------------------------------
# One directory per sample is expected under CELLRANGER_DIR; sample order here
# defines the numeric sample index used by Stage 1 and by SAMPLE_RDS_DIR file
# names (<index>.rds).
SAMPLE_DIRS  <- list.files(CELLRANGER_DIR, full.names = TRUE)
SAMPLE_NAMES <- basename(SAMPLE_DIRS)
N_SAMPLES    <- length(SAMPLE_DIRS)

# Prefix stripped from `orig.ident` to obtain a readable sample label, e.g.
# "SC022_C1_LPS-ATP" -> "C1_LPS-ATP".
SAMPLE_NAME_PREFIX <- "SC022_"

# Maps the subject/batch code embedded in each sample name to a numeric batch
# ID (samples processed together on the same day/chip are treated as one
# Harmony batch). Update this if the sample naming scheme changes.
BATCH_MAP <- c(P1 = 1, P2 = 2, C1 = 3, C2 = 4, C3 = 5, P3 = 6)

# --- Analysis parameters --------------------------------------------------------
N_PCS         <- 50    # principal components used throughout (PCA, neighbors, UMAP, Harmony)
CLUSTER_RES   <- 0.8   # Louvain clustering resolution (per-sample and integrated)
RANDOM_SEED   <- 123

# Per-sample doublet rate model: fitted linear relationship between the
# number of cells recovered and the expected multiplet rate, derived from
# 10x Genomics' published cell-recovery / multiplet-rate table. 
DOUBLET_RATE_INTERCEPT <- 5.355917e-16 # residual/rounding artifact, not a meaningful term
DOUBLET_RATE_SLOPE     <- 8.000000e-04

# Post-merge cell QC thresholds (applied once, after merging all samples)
QC_MIN_FEATURES <- 500
QC_MAX_FEATURES <- 7500
QC_MIN_COUNTS   <- 1000
QC_MAX_COUNTS   <- 1.5e5
QC_MAX_PCT_MT   <- 10

# Classic marker genes (used for visual QC of the integrated clustering)
MARKER_GENES <- c(
  "CD3D", "CD3E", "CD3G", "CCR7", "SELL", "LTB", "IL7R", "S100A4",
  "CD8A", "CD8B", "GZMB", "PRF1", "FCGR3A", "NKG7", "GNLY", "KLRC1",
  "CD14", "LYZ", "S100A9", "CSF3R", "IFITM1", "IFITM2", "IFITM3",
  "CD79A", "MS4A1", "CD1C", "ITGAX", "CLEC4C", "GP9", "CD34", "THY1",
  "ENG", "KIT", "PROM1", "CCL5", "IL32", "PTPRCAP", "PF4"
)

# Gene of interest used for the final expression-overlay figure
GENE_OF_INTEREST <- "NLRP3"

# If TRUE, Stage 2 will reuse any existing checkpoint .rds files instead of
# recomputing expensive steps (SCTransform, Harmony). Set to FALSE to force
# a full recompute.
USE_CHECKPOINTS <- TRUE

## =============================================================================
## 2. Libraries
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat) ## v.4.0.8
  library(SoupX)
  library(DoubletFinder)
  library(SingleR)
  library(celldex)
  library(glmGamPoi)
  library(harmony)
  library(dittoSeq)
  library(dplyr)
  library(ggplot2)
  library(BiocParallel)
})

set.seed(RANDOM_SEED)

## =============================================================================
## 3. Shared helper functions
## =============================================================================

#' Save an object to `path`, or load it from `path` if it already exists.
#'
#' Used as a lightweight checkpointing mechanism so that expensive steps
#' (SCTransform, Harmony) do not need to be recomputed when re-running or
#' resuming the integration stage.
#'
#' @param path Path to the checkpoint .rds file.
#' @param compute A zero-argument function that computes the object if no
#'   checkpoint is found.
#' @param use_checkpoint If FALSE, always recompute and overwrite the file.
run_or_load <- function(path, compute, use_checkpoint = USE_CHECKPOINTS) {
  if (use_checkpoint && file.exists(path)) {
    message("Loading checkpoint: ", path)
    return(readRDS(path))
  }
  message("Computing (no checkpoint found): ", path)
  result <- compute()
  saveRDS(result, path)
  result
}

#' Save a ggplot/plot object and immediately free it.
save_plot <- function(plot_obj, filename, width, height, dpi = 600) {
  ggsave(file.path(PLOTS_DIR, filename), plot_obj, width = width, height = height, units = "in", dpi = dpi)
  invisible(NULL)
}

#' pK selection helper for DoubletFinder.
#'
#' Reimplementation of `DoubletFinder::find.pK()` that does not attempt to
#' plot anything. The original `find.pK()` calls base-R plotting functions
#' unconditionally, which fails when the pipeline is run non-interactively
#' (e.g. inside a batch job with no graphics device).
#'
#' @param sweep_stats Output of `DoubletFinder::summarizeSweep()`.
#' @return A data frame of per-pK BCmetric summary statistics, identical in
#'   content to `find.pK()`'s return value.
find_pK_no_plot <- function(sweep_stats) {
  has_auc <- "AUC" %in% colnames(sweep_stats)
  pk_values <- unique(sweep_stats$pK)

  ncols <- if (has_auc) 6 else 5
  bc_mvn <- as.data.frame(matrix(0L, nrow = length(pk_values), ncol = ncols))
  colnames(bc_mvn) <- if (has_auc) {
    c("ParamID", "pK", "MeanAUC", "MeanBC", "VarBC", "BCmetric")
  } else {
    c("ParamID", "pK", "MeanBC", "VarBC", "BCmetric")
  }
  bc_mvn$pK <- pk_values
  bc_mvn$ParamID <- seq_along(pk_values)
 
  for (x in seq_along(pk_values)) {
    idx <- which(sweep_stats$pK == pk_values[x])
    if (has_auc) bc_mvn$MeanAUC[x] <- mean(sweep_stats[idx, "AUC"])
    bc_mvn$MeanBC[x]   <- mean(sweep_stats[idx, "BCreal"])
    bc_mvn$VarBC[x]    <- sd(sweep_stats[idx, "BCreal"])^2
    bc_mvn$BCmetric[x] <- bc_mvn$MeanBC[x] / bc_mvn$VarBC[x]
  }

  return(bc_mvn)
}

#' Run DoubletFinder on a clustered Seurat object and attach a `doublet`
#' metadata column ("Singlet" / "Doublet").
#'
#' Two passes of DoubletFinder are run: the first estimates the raw expected
#' number of doublets from the assumed formation rate, and the second
#' re-estimates it after excluding homotypic (same cell type) doublets. 
#' The adjusted (second-pass) call is what gets stored.
#'
#' @param seurat_obj A normalized, PCA-reduced, and clustered Seurat object.
#' @param n_pcs Number of PCs to use (must match the PCs used for clustering).
find_doublets <- function(seurat_obj, n_pcs = N_PCS) {
  sweep_res    <- paramSweep_v3(seurat_obj, PCs = 1:n_pcs, sct = FALSE)
  sweep_stats  <- summarizeSweep(sweep_res, GT = FALSE)
  bc_mvn       <- find_pK_no_plot(sweep_stats)

  best_pK <- as.numeric(as.character(
    bc_mvn$pK[which(bc_mvn$BCmetric == max(bc_mvn$BCmetric))]
  ))

  homotypic_prop <- modelHomotypic(Idents(seurat_obj))
  n_cells <- nrow(seurat_obj@meta.data)
  expected_doublet_rate <- (DOUBLET_RATE_INTERCEPT + n_cells * DOUBLET_RATE_SLOPE) / 100
  n_exp_raw      <- round(expected_doublet_rate * n_cells)
  n_exp_adjusted <- round(n_exp_raw * (1 - homotypic_prop))

  seurat_obj <- doubletFinder_v3(
    seurat_obj, PCs = 1:n_pcs, pN = 0.25, pK = best_pK,
    nExp = n_exp_raw, reuse.pANN = FALSE, sct = FALSE
  )
  pann_col <- paste0("pANN_0.25_", best_pK, "_", n_exp_raw)
  seurat_obj <- doubletFinder_v3(
    seurat_obj, PCs = 1:n_pcs, pN = 0.25, pK = best_pK,
    nExp = n_exp_adjusted, reuse.pANN = pann_col, sct = FALSE
  )

  classification_col <- paste0("DF.classifications_0.25_", best_pK, "_", n_exp_adjusted)
  seurat_obj$doublet <- seurat_obj@meta.data[[classification_col]]
  return(seurat_obj)
}

#' Annotate cluster identities of a Seurat object with SingleR, using the
#' Human Primary Cell Atlas reference, and return both the fine-grained and
#' the broad ("main") label for each cluster.
#'
#' @param seurat_obj Seurat object with an active cluster identity column.
#' @param cluster_col Name of the metadata column holding cluster labels.
#' @param reference Reference dataset (e.g. from `celldex::HumanPrimaryCellAtlasData()`).
#' @return A named list with `fine` and `broad` character vectors, in cell
#'   order matching `colnames(seurat_obj)`.
annotate_clusters_singleR <- function(seurat_obj, cluster_col, reference) {
  test_assay <- GetAssayData(seurat_obj)
  clusters <- seurat_obj[[cluster_col, drop = TRUE]]

  prediction <- SingleR(test = test_assay, ref = reference, clusters = clusters, labels = reference$label.fine)
  pred_df <- as.data.frame(prediction)
  pred_df <- data.frame(label = pred_df$pruned.labels, cluster = rownames(pred_df))

  fine_by_cluster <- pred_df$label[match(clusters, pred_df$cluster)]

  label_map <- unique(as.data.frame(reference@colData@listData)[, c("label.main", "label.fine")])
  broad_by_cluster <- label_map$label.main[match(fine_by_cluster, label_map$label.fine)]

  list(fine = fine_by_cluster, broad = broad_by_cluster)
}

## =============================================================================
## 4. Stage 1 -- per-sample preprocessing
## =============================================================================
#
# For each sample: read the raw and filtered Cell Ranger matrices, build a
# preliminary Seurat object to (a) get cluster labels for SoupX and (b) call
# doublets, correct counts for ambient RNA with SoupX, then rebuild and
# re-analyze the corrected object. Cell-type scores from SingleR are added
# for a fast, cluster-level annotation; a finer, integrated annotation is
# produced later in Stage 2.

preprocess_sample <- function(sample_index) {
  stopifnot(sample_index >= 1, sample_index <= N_SAMPLES)
  sample_name <- SAMPLE_NAMES[sample_index]
  message(sprintf("[Stage 1] Sample %d/%d: %s", sample_index, N_SAMPLES, sample_name))

  raw_dir      <- file.path(SAMPLE_DIRS[sample_index], "outs", "raw_feature_bc_matrix")
  filtered_dir <- file.path(SAMPLE_DIRS[sample_index], "outs", "filtered_feature_bc_matrix")

  # Barcodes are suffixed with the sample index (instead of Cell Ranger's
  # default "-1") so that they stay unique once all samples are merged.
  message("Reading raw and filtered count matrices")
  raw_counts <- Read10X(raw_dir)
  colnames(raw_counts) <- sub("-1", paste0("-", sample_index), colnames(raw_counts), fixed = TRUE)

  filtered_counts <- Read10X(filtered_dir)
  colnames(filtered_counts) <- sub("-1", paste0("-", sample_index), colnames(filtered_counts), fixed = TRUE)

  # --- Preliminary object: used only to get cluster labels for SoupX and to
  # call doublets prior to ambient-RNA correction. ---
  message("Building preliminary Seurat object")
  seu <- CreateSeuratObject(counts = filtered_counts, project = sample_name, min.cells = 0, min.features = 0)
  seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL][[:digit:]]")
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu, vars.to.regress = "percent.mt")
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:N_PCS)
  seu <- FindClusters(seu, resolution = CLUSTER_RES)
  seu <- suppressMessages(suppressWarnings(find_doublets(seu)))
  doublet_calls <- seu@meta.data[, "doublet", drop = FALSE]

  # --- Ambient RNA correction (SoupX) ---
  message("Correcting ambient RNA with SoupX")
  soup_channel <- SoupChannel(raw_counts, filtered_counts)
  soup_channel <- setClusters(soup_channel, seu$seurat_clusters)
  soup_channel <- autoEstCont(soup_channel)
  corrected_counts <- adjustCounts(soup_channel)

  # --- Final per-sample object, built from SoupX-corrected counts ---
  message("Re-analyzing SoupX-corrected counts")
  seu_adj <- CreateSeuratObject(counts = corrected_counts, project = sample_name, min.cells = 0, min.features = 0)
  seu_adj[["percent.mt"]]   <- PercentageFeatureSet(seu_adj, pattern = "^MT-")
  seu_adj[["percent.ribo"]] <- PercentageFeatureSet(seu_adj, pattern = "^RP[SL][[:digit:]]")
  seu_adj <- NormalizeData(seu_adj)
  seu_adj <- FindVariableFeatures(seu_adj, selection.method = "vst", nfeatures = 2000)
  seu_adj <- ScaleData(seu_adj, vars.to.regress = "percent.mt")
  seu_adj <- RunPCA(seu_adj, verbose = FALSE)
  seu_adj <- RunUMAP(seu_adj, reduction = "pca", dims = 1:N_PCS)
  seu_adj <- RunTSNE(seu_adj, reduction = "pca", dims = 1:N_PCS)
  seu_adj <- FindNeighbors(seu_adj, reduction = "pca", dims = 1:N_PCS)
  seu_adj <- FindClusters(seu_adj, resolution = CLUSTER_RES)

  # Doublet calls were made on the pre-correction object; carry them over by
  # cell barcode (SoupX does not add or remove cells, only adjusts counts).
  seu_adj$doublet <- doublet_calls$doublet[match(rownames(seu_adj@meta.data), rownames(doublet_calls))]

  # --- Cluster-level cell type annotation (SingleR, HPCA reference) ---
  message("Annotating clusters with SingleR (Human Primary Cell Atlas)")
  hpca_ref <- HumanPrimaryCellAtlasData()
  singleR_labels <- annotate_clusters_singleR(seu_adj, "seurat_clusters", hpca_ref)
  seu_adj$SingleR_celltype_cluster <- singleR_labels$fine

  out_path <- file.path(SAMPLE_RDS_DIR, paste0(sample_index, ".rds"))
  message("Writing ", out_path)
  saveRDS(seu_adj, out_path)
  message("Sample ", sample_index, " finished successfully")
  invisible(seu_adj)
}

## =============================================================================
## 5. Stage 2 -- integration and downstream analysis
## =============================================================================

#' Parse subject/batch, treatment condition, and disease status out of a
#' sample's `orig.ident`, and add a readable `sample` label.
#'
#' Expects `orig.ident` values of the form "<prefix><subject>_<condition>",
#' e.g. "SC022_C1_LPS-ATP", where <subject> is one of the codes in
#' BATCH_MAP and <condition> is one of "Untreated", "LPS", "LPS-ATP".
add_sample_metadata <- function(seurat_obj) {
  ident <- seurat_obj$orig.ident

  seurat_obj$sample <- sub(SAMPLE_NAME_PREFIX, "", ident, fixed = TRUE)

  seurat_obj$batch <- NA_integer_
  for (code in names(BATCH_MAP)) {
    seurat_obj$batch[grepl(paste0("_", code, "_"), ident, fixed = TRUE)] <- BATCH_MAP[[code]]
  }

  condition_token <- vapply(
    strsplit(gsub("LPS_ATP", "LPS-ATP", ident, fixed = TRUE), "_", fixed = TRUE),
    function(parts) parts[[3]],
    character(1)
  )
  seurat_obj$condition <- factor(condition_token, levels = c("Untreated", "LPS", "LPS-ATP"))

  seurat_obj$disease_status <- ifelse(grepl("_C", ident, fixed = TRUE), "Control", "Patient")

  return(seurat_obj)
}

run_integration <- function() {
  ## --- 5.1 Load per-sample objects and add cohort metadata -------------------
  sample_files <- file.path(SAMPLE_RDS_DIR, paste0(seq_len(N_SAMPLES), ".rds"))
  missing <- !file.exists(sample_files)
  if (any(missing)) {
    stop("Missing Stage 1 output for sample(s): ", paste(which(missing), collapse = ", "),
         ". Run `preprocess_sample()` for these first.")
  }

  message("[Stage 2] Loading ", N_SAMPLES, " per-sample objects")
  sample_list <- pbapply::pblapply(sample_files, function(f) add_sample_metadata(readRDS(f)))

  ## --- 5.2 Merge samples -------------------------------------------------------
  # `merge.dr = "umap"` keeps each sample's own (pre-integration) UMAP
  # embedding around for later QC comparison against the integrated UMAP.
  merged <- run_or_load(file.path(CHECKPOINT_DIR, "merged.rds"), function() {
    obj <- sample_list[[1]]
    for (i in 2:length(sample_list)) {
      message("Merging sample ", i, "/", length(sample_list))
      obj <- merge(x = obj, y = sample_list[[i]], merge.dr = "umap")
    }
    return(obj)
  })
  rm(sample_list); gc()

  ## --- 5.3 QC filtering --------------------------------------------------------
  p_pre <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                    ncol = 3, group.by = "sample", pt.size = 0, split.by = "condition")
  save_plot(p_pre, "QC_pre_filtering.png", width = 16.6, height = 11.6)

  merged <- subset(
    merged,
    subset = nFeature_RNA > QC_MIN_FEATURES & nFeature_RNA < QC_MAX_FEATURES &
      nCount_RNA > QC_MIN_COUNTS & nCount_RNA < QC_MAX_COUNTS &
      percent.mt < QC_MAX_PCT_MT & doublet == "Singlet"
  )

  p_post <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                     ncol = 3, group.by = "sample", pt.size = 0, split.by = "condition")
  save_plot(p_post, "QC_post_filtering.png", width = 16.6, height = 11.6)
  rm(p_pre, p_post); gc()

  ## --- 5.4 Normalization (SCTransform) -----------------------------------------
  merged <- run_or_load(file.path(CHECKPOINT_DIR, "merged_sct.rds"), function() {
    SCTransform(merged, method = "glmGamPoi", verbose = TRUE, vars.to.regress = "percent.mt")
  })

  ## --- 5.5 Pre-integration diagnostics ------------------------------------------
  merged <- RunPCA(merged, verbose = FALSE)

  # Keep the per-sample UMAP embeddings computed in Stage 1 for comparison,
  # then clear the "umap" slot so RunUMAP() below writes the integrated one.
  merged[["ind_samp_umap"]] <- CreateDimReducObject(
    embeddings = merged@reductions[["umap"]]@cell.embeddings,
    key = "UMAPindsamp_", assay = DefaultAssay(merged)
  )
  merged[["umap"]] <- NULL

  p <- dittoDimPlot(merged, "sample", reduction.use = "ind_samp_umap", split.by = "sample", show.others = FALSE)
  save_plot(p, "UMAP_per_sample_preintegration.png", width = 8.3, height = 5.8)
  p <- dittoDimPlot(merged, "sample", reduction.use = "ind_samp_umap")
  save_plot(p, "UMAP_per_sample_preintegration_overlay.png", width = 8.3, height = 5.8)
  rm(p); gc()

  set.seed(RANDOM_SEED) # RunUMAP is stochastic; fixed here for reproducibility
  merged <- RunUMAP(merged, reduction = "pca", dims = 1:N_PCS)
  p <- dittoDimPlot(merged, "sample", reduction.use = "umap", split.by = "sample")
  save_plot(p, "UMAP_merged_preharmony.png", width = 8.3, height = 5.8)
  p <- UMAPPlot(merged, group.by = c("batch", "disease_status"), split.by = "condition")
  save_plot(p, "UMAP_batch_disease_preharmony.png", width = 8.3, height = 5.8)
  rm(p); gc()

  ## --- 5.6 Batch integration (Harmony) -------------------------------------------
  merged <- run_or_load(file.path(CHECKPOINT_DIR, "merged_harmony.rds"), function() {
    set.seed(RANDOM_SEED)
    meta <- merged@meta.data
    # coerce to character.
    meta[, c("batch", "disease_status", "condition", "sample")] <-
      lapply(meta[, c("batch", "disease_status", "condition", "sample")], as.character)

    harmony_embeddings <- HarmonyMatrix(
      data_mat = Embeddings(merged, "pca"),
      meta_data = meta,
      vars_use = "sample",
      theta = 2,
      lambda = 1,
      do_pca = FALSE
    )
    merged[["harmony"]] <- CreateDimReducObject(embeddings = harmony_embeddings, key = "harmony_", assay = c("SCT", "RNA"))
    merged <- RunUMAP(merged, reduction = "harmony", dims = 1:N_PCS)
    merged
  })

  p <- dittoDimPlot(merged, "sample", reduction.use = "umap", split.by = "sample")
  save_plot(p, "UMAP_merged_postharmony.png", width = 8.3, height = 5.8)
  p <- UMAPPlot(merged, group.by = c("batch", "disease_status"), split.by = "condition")
  save_plot(p, "UMAP_batch_disease_postharmony.png", width = 8.3, height = 5.8)
  rm(p); gc()

  ## --- 5.7 Clustering on the Harmony-corrected space --------------------------
  merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:N_PCS)
  merged <- FindClusters(merged, graph.name = "SCT_nn", resolution = CLUSTER_RES)

  cluster_levels <- as.character(sort(as.numeric(as.character(unique(merged$SCT_nn_res.0.8)))))
  merged$cluster_harmony <- factor(merged$SCT_nn_res.0.8, levels = cluster_levels)

  p <- UMAPPlot(merged, group.by = "cluster_harmony", label = TRUE, cols = dittoColors())
  save_plot(p, "UMAP_clusters_harmony.png", width = 8.3, height = 5.8)
  rm(p); gc()

  p <- ggplot(merged@meta.data, aes(x = cluster_harmony, fill = sample)) +
    geom_bar(position = "fill") +
    labs(y = "Proportion", x = "Cluster") +
    scale_fill_manual(values = dittoColors()) +
    geom_text(aes(label = after_stat(count)), stat = "count", position = "fill", size = 1.5, alpha = 0.4)
  save_plot(p, "sample_proportions_by_cluster.png", width = 8.3, height = 5.8)
  rm(p); gc()

  ## --- 5.8 Cell type annotation (SingleR, HPCA reference) ----------------------
  hpca_ref <- HumanPrimaryCellAtlasData()

  message("SingleR: cluster-level annotation")
  cluster_labels <- annotate_clusters_singleR(merged, "cluster_harmony", hpca_ref)
  merged$SingleR_cluster       <- cluster_labels$fine
  merged$SingleR_cluster_broad <- cluster_labels$broad

  p <- UMAPPlot(merged, group.by = "SingleR_cluster")
  save_plot(p, "UMAP_SingleR_cluster.png", width = 10, height = 5.8)
  p <- UMAPPlot(merged, group.by = "SingleR_cluster_broad")
  save_plot(p, "UMAP_SingleR_cluster_broad.png", width = 10, height = 5.8)
  rm(p); gc()

  message("SingleR: single-cell-level annotation (chunked)")
  chunk_size <- 5000
  n_cells <- ncol(merged)
  cell_chunks <- split(seq_len(n_cells), (seq_len(n_cells) - 1) %/% chunk_size)

  bp_param <- MulticoreParam(workers = 8)
  sct_data <- GetAssayData(merged, assay = "SCT", slot = "data")
  chunk_predictions <- pbapply::pblapply(cell_chunks, function(idx) {
    SingleR(test = sct_data[, idx], ref = hpca_ref, assay.type.test = 1,
            labels = hpca_ref$label.fine, num.threads = 8, BPPARAM = bp_param)
  })
  predictions <- data.table::rbindlist(lapply(chunk_predictions, function(p) {
    df <- as.data.frame(p)
    df$cell_id <- rownames(p)
    df[, c("pruned.labels", "cell_id")]
  }))

  merged$SingleR_singleCell <- predictions$pruned.labels[match(colnames(merged), predictions$cell_id)]
  label_map <- unique(as.data.frame(hpca_ref@colData@listData)[, c("label.main", "label.fine")])
  merged$SingleR_singleCell_broad <- label_map$label.main[match(merged$SingleR_singleCell, label_map$label.fine)]
  rm(hpca_ref, chunk_predictions, predictions, cell_chunks, bp_param, sct_data); gc()

  p <- UMAPPlot(merged, group.by = "SingleR_singleCell")
  save_plot(p, "UMAP_SingleR_singlecell.png", width = 20, height = 11.6)
  p <- UMAPPlot(merged, group.by = "SingleR_singleCell_broad", label = TRUE, label.size = 2)
  save_plot(p, "UMAP_SingleR_singlecell_broad.png", width = 10, height = 5.8)
  rm(p); gc()

  merged <- run_or_load(file.path(CHECKPOINT_DIR, "merged_annotated.rds"), function() merged)

  ## --- 5.9 Marker gene visualization --------------------------------------------
  present_markers <- intersect(MARKER_GENES, rownames(merged))
  if (length(present_markers) < length(MARKER_GENES)) {
    warning(length(MARKER_GENES) - length(present_markers), " marker genes not found in the assay and were skipped")
  }

  p <- VlnPlot(merged, features = present_markers, pt.size = 0)
  save_plot(p, "markers_violin.png", width = 33.2, height = 23.2, dpi = 500)
  rm(p); gc()

  p <- FeaturePlot(merged, features = present_markers, reduction = "umap")
  save_plot(p, "markers_umap.png", width = 33.2, height = 34.8, dpi = 400)
  rm(p); gc()

  ## --- 5.10 Cluster marker genes (differential expression vs. rest) -------------
  merged <- SetIdent(merged, value = "cluster_harmony")
  cluster_markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  saveRDS(cluster_markers, file.path(MARKERS_DIR, "cluster_markers_harmony.rds"))

  for (cl in unique(cluster_markers$cluster)) {
    write.csv(
      cluster_markers[cluster_markers$cluster == cl, ],
      file.path(MARKERS_DIR, paste0("cluster_", cl, "_markers.csv")),
      row.names = FALSE, quote = FALSE
    )
  }

  top_markers <- cluster_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
  p <- DotPlot(merged, features = unique(top_markers$gene)) + RotatedAxis() +
    scale_colour_gradient2(low = "red", mid = "white", high = "green") +
    theme(plot.background = element_rect(fill = "white"))
  save_plot(p, "markers_dotplot_top2_per_cluster.png", width = 12, height = 8)
  rm(p); gc()

  ## --- 5.11 Manual coarse/fine cell-type annotation -----------------------------
  # NOTE: these labels were assigned by inspecting `cluster_markers_harmony.rds`
  # / the dot plot above against known marker genes, for the specific cluster
  # numbering obtained on the reference run of this pipeline. If the input
  # data, QC thresholds, or package versions change, cluster numbers and
  # membership may shift -- re-derive these mappings from the marker table
  # before trusting them on a rerun.
  broad_annotation_by_cluster <- c(
    "T", "NK", "T", "T", "T", "Myeloid", "T", "T", "B", "T",
    "T", "B", "T", "T", "T", "NK", "T? Metamyelocyte?", "T",
    "B", "NK", "T", "T", "Myeloid", "T?", "T?", "Megakaryocyte"
  )
  narrow_annotation_by_cluster <- c(
    "CD4+ Central Memory", "NK", "CD8+", "CD4+ Naive", "CD4+ Naive",
    "Monocytes, DCs, Macrophages", "CD4+ Naive", "CD4+ Naive", "B",
    "CD4+ Naive", "CD4+ Effector Memory", "B", "CD8+", "CD4+ Naive",
    "CD4+ Naive", "NK", "CD8+?", "CD4+ Naive", "B", "NK", "CD4+ Naive",
    "CD4+ Central Memory", "Neutrophil", "CD4+?", "CD4+ Central Memory?",
    "Megakaryocyte"
  )
  stopifnot(length(broad_annotation_by_cluster) == length(levels(merged$cluster_harmony)))
  stopifnot(length(narrow_annotation_by_cluster) == length(levels(merged$cluster_harmony)))

  merged$broad_annot <- merged$cluster_harmony
  levels(merged$broad_annot) <- broad_annotation_by_cluster
  levels(merged$broad_annot)[grepl("\\?", levels(merged$broad_annot))] <- "?"

  merged$narrow_annot <- merged$cluster_harmony
  levels(merged$narrow_annot) <- narrow_annotation_by_cluster
  levels(merged$narrow_annot)[grepl("\\?", levels(merged$narrow_annot))] <- "?"

  ## --- 5.12 Gene-of-interest overlay on annotated clusters ----------------------
  plot_df <- cbind(merged@meta.data, merged@reductions$umap@cell.embeddings)
  plot_df$gene_expr <- GetAssayData(merged, assay = "SCT", slot = "data")[GENE_OF_INTEREST, ]
  plot_df <- plot_df[plot_df$broad_annot != "?", ]

  p <- ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = broad_annot), plot_df, size = 10) +
    geom_point(aes(x = UMAP_1, y = UMAP_2), plot_df, size = 8, color = "white") +
    ggnewscale::new_scale_color() +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = gene_expr), plot_df[order(plot_df$gene_expr), ], size = 2) +
    labs(color = paste0(GENE_OF_INTEREST, " expression"))
  save_plot(p, paste0("clusters_and_", GENE_OF_INTEREST, ".png"), width = 16.6, height = 9.6)
  rm(p); gc()

  narrow_order <- c("CD4+ Naive", "CD4+ Central Memory", "CD4+ Effector Memory", "CD8+", "B",
                     "?", "NK", "Monocytes, DCs, Macrophages", "Neutrophil", "Megakaryocyte")
  plot_df$narrow_annot <- factor(as.character(plot_df$narrow_annot), levels = narrow_order)
  p <- ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2, color = narrow_annot), plot_df, size = 2)
  save_plot(p, "clusters_annotations.png", width = 16.6, height = 9.6)
  rm(p, plot_df); gc()

  ## --- 5.13 Final metadata cleanup ------------------------------------------------
  # Drop intermediate/redundant columns produced along the way; keep
  # everything analysis-relevant.
  columns_to_drop <- c("RNA_snn_res.0.8", "SCT_nn_res.0.8", "SingleR_celltype_cluster")
  for (col in columns_to_drop) merged[[col]] <- NULL

  final_path <- file.path(CHECKPOINT_DIR, "seurat_final_annotated.rds")
  saveRDS(merged, final_path)
  message("Saved final annotated object to ", final_path)

  ## --- 5.14 Differential expression: treatment vs. untreated, per cell type ------
  message("Running per-cell-type differential expression (LPS vs. Untreated, LPS-ATP vs. Untreated)")
  merged <- SetIdent(merged, value = "condition")
  cell_types <- setdiff(unique(as.character(merged$narrow_annot)), "?")

  de_results <- setNames(vector("list", length(cell_types)), cell_types)
  for (cell_type in cell_types) {
    message(" - ", cell_type)
    subset_obj <- subset(merged, subset = narrow_annot == cell_type)

    markers_lps <- tryCatch(
      FindMarkers(subset_obj, ident.1 = "LPS", ident.2 = "Untreated", group.by = "condition"),
      error = function(e) { warning("LPS vs Untreated failed for ", cell_type, ": ", conditionMessage(e)); NULL }
    )
    markers_lps_atp <- tryCatch(
      FindMarkers(subset_obj, ident.1 = "LPS-ATP", ident.2 = "Untreated", group.by = "condition"),
      error = function(e) { warning("LPS-ATP vs Untreated failed for ", cell_type, ": ", conditionMessage(e)); NULL }
    )

    de_results[[cell_type]] <- list(LPS = markers_lps, `LPS-ATP` = markers_lps_atp)

    safe_name <- gsub("[^A-Za-z0-9_+-]", "_", cell_type)
    if (!is.null(markers_lps)) {
      write.csv(markers_lps, file.path(DE_DIR, paste0(safe_name, "_LPS_vs_Untreated.csv")))
    }
    if (!is.null(markers_lps_atp)) {
      write.csv(markers_lps_atp, file.path(DE_DIR, paste0(safe_name, "_LPS-ATP_vs_Untreated.csv")))
    }
  }
  saveRDS(de_results, file.path(DE_DIR, "differential_expression_by_celltype.rds"))

  message("[Stage 2] Integration and downstream analysis complete.")
  invisible(merged)
}

## =============================================================================
## 6. Entry point
## =============================================================================
# Allows the script to be run non-interactively for either stage:
#   Rscript scRNAseq_analysis_pipeline.R preprocess <sample_index>
#   Rscript scRNAseq_analysis_pipeline.R integrate
# When sourced interactively (no command-line arguments), neither stage runs
# automatically -- call `preprocess_sample()` / `run_integration()` directly.

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[1] == "preprocess") {
    if (length(args) < 2) stop("Usage: Rscript scRNAseq_analysis_pipeline.R preprocess <sample_index>")
    preprocess_sample(as.integer(args[2]))
  } else if (length(args) >= 1 && args[1] == "integrate") {
    run_integration()
  } else {
    message(
      "No stage specified. Usage:\n",
      "  Rscript scRNAseq_analysis_pipeline.R preprocess <sample_index>\n",
      "  Rscript scRNAseq_analysis_pipeline.R integrate"
    )
  }
}


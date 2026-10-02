#!/usr/bin/env Rscript
# ==============================================================================
# Per-cell-type differential expression analysis
#
# For each cell type of interest, this script runs two families of contrasts
# on the final, annotated, integrated Seurat object:
#
#   (A) Patient vs. Control, within each treatment condition
#         (Untreated / LPS / LPS-ATP)
#   (B) LPS vs. Untreated and LPS-ATP vs. Untreated, within each disease
#         status (Patient / Control)
#
# For every contrast it:
#   - runs Seurat::FindMarkers() and categorizes genes as up/down/below
#     threshold from adjusted p-value and log2 fold-change,
#   - computes a leave-one-donor-out "unity score": how many of the
#     relevant donors' single-subject exclusions still agree with the
#     full-cohort call, used to flag genes whose significance is not
#     being driven by any one subject,
#   - writes the full marker table and a "robust" subset to .xlsx,
#   - draws two volcano plots (top hits; NF-kB / NLRP3-inflammasome
#     target genes highlighted),
#   - runs GO enrichment (up- and down-regulated genes separately) when
#     there are enough significant genes.
#
# This script is downstream of, and expects as input, the final annotated
# object produced by the integration pipeline (`seurat_final_annotated.rds`),
# which already carries `narrow_annot`, `condition`, `disease_status`, and
# `sample` metadata columns. 
#
# Usage:
#   Rscript diffexp_and_ORA.R
# ==============================================================================

## =============================================================================
## 1. Configuration
## =============================================================================

BASE_DIR            <- "path_to_dir/KFH_scRNAseq/"
SEURAT_OBJECT_PATH  <- file.path(BASE_DIR, "checkpoints", "seurat_final_annotated.rds")
DE_OUTPUT_DIR       <- file.path(BASE_DIR, "ct_diffexp")

# Cell types (as labeled in `narrow_annot`) to run this analysis on.
CELL_TYPES_OF_INTEREST <- c("Neutrophil", "Monocytes, DCs, Macrophages")

# Statistical thresholds
PADJ_THRESHOLD           <- 0.05
STRONG_LOG2FC_THRESHOLD  <- log2(1.5)   # "strong"/"robust" effect-size cutoff
MANY_SIG_GENES_THRESHOLD <- 100         # above this many significant genes, tighten the GO input set
MIN_GENES_FOR_GO         <- 10          # clusterProfiler::enrichGO() needs a reasonably sized gene set

# Gene sets highlighted on volcano plots
NFKB_TARGET_GENES  <- c("MYD88", "TICAM1", "FADD")
NLRP3_TARGET_GENES <- c("CASP1", "NLRP3", "PYCARD")

# GO enrichment settings
GO_ORG_DB     <- "org.Hs.eg.db"
GO_KEY_TYPE   <- "SYMBOL"
GO_ONTOLOGY   <- "ALL"
GO_MIN_PCT_DETECTED <- 1   # minimum % of cells a gene must be detected in to enter the GO background/universe

## =============================================================================
## 2. Libraries
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggh4x)
  library(openxlsx)
})

options(SEURAT_OPTIONS)
OrgDb <- get(GO_ORG_DB)

dir.create(DE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## 3. Helper functions
## =============================================================================

#' Add `gene`, `significant`, and `category` columns to a FindMarkers() result.
#'
#' @param marker_genes Output of `Seurat::FindMarkers()`.
categorize_markers <- function(marker_genes) {
  marker_genes$gene <- rownames(marker_genes)
  marker_genes$significant <- marker_genes$p_val_adj < PADJ_THRESHOLD
  marker_genes$category <- dplyr::case_when(
    marker_genes$avg_log2FC > 0 & marker_genes$p_val_adj < PADJ_THRESHOLD ~ "upregulated",
    marker_genes$avg_log2FC < 0 & marker_genes$p_val_adj < PADJ_THRESHOLD ~ "downregulated",
    TRUE ~ "below threshold"
  )
  marker_genes
}

#' Flag genes belonging to the study's pathways of interest.
flag_pathway_genes <- function(marker_genes) {
  marker_genes$NFKB_target_gene <- marker_genes$gene %in% NFKB_TARGET_GENES
  marker_genes$NLRP3_inflammasome_target_gene <- marker_genes$gene %in% NLRP3_TARGET_GENES
  marker_genes
}

#' Leave-one-donor-out robustness ("unity") score.
#'
#' For each donor in `donors_to_test`, excludes that donor's cells from
#' `subset_obj`, reruns the same two-group contrast, and checks whether each
#' gene's regulation call (up / down / below threshold) still agrees with the
#' full-cohort call in `marker_genes`. The number of agreeing reruns is
#' returned as `unity_score`, alongside `n_donors_tested` (its maximum
#' possible value) -- a gene is only as "robust" as `unity_score ==
#' n_donors_tested`.
#'
#' Genes absent from a given leave-one-out run (e.g. filtered out for low
#' expression once a donor's cells are removed) are counted as disagreeing,
#' which is the conservative choice.
#'
#' @param marker_genes Output of `categorize_markers()` for the full-cohort contrast.
#' @param subset_obj The Seurat object subset used for the full-cohort contrast.
#' @param donors_to_test Character vector of donor codes to leave out, one at a time.
#' @param donor_col Metadata column identifying each donor (e.g. "donor").
#' @param group_col Metadata column defining the two groups being compared.
#' @param ident_1,ident_2 The two group labels being compared (as passed to `FindMarkers()`).
compute_unity_score <- function(marker_genes, subset_obj, donors_to_test, donor_col, group_col, ident_1, ident_2) {
  agreement <- integer(nrow(marker_genes))

  for (donor in donors_to_test) {
    keep_cells <- colnames(subset_obj)[subset_obj[[donor_col, drop = TRUE]] != donor]
    loo_obj <- subset(subset_obj, cells = keep_cells)
    loo_obj <- SetIdent(loo_obj, value = group_col)

    loo_markers <- tryCatch(
      FindMarkers(loo_obj, ident.1 = ident_1, ident.2 = ident_2, group.by = group_col),
      error = function(e) {
        message("  Leave-one-out (excluding ", donor, ") failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(loo_markers)) next

    loo_markers <- categorize_markers(loo_markers)
    matched_category <- loo_markers$category[match(marker_genes$gene, loo_markers$gene)]
    # A gene missing from the leave-one-out result (e.g. filtered out once
    # this donor's cells are removed) does not count as agreement.
    agreement <- agreement + (!is.na(matched_category) & matched_category == marker_genes$category)
  }

  marker_genes$unity_score <- agreement
  marker_genes$n_donors_tested <- length(donors_to_test)
  marker_genes
}

#' A gene is a "robust" DEG if it clears the strong effect-size and
#' significance thresholds, AND its regulation call agreed across every
#' leave-one-donor-out rerun.
is_robust_deg <- function(marker_genes) {
  abs(marker_genes$avg_log2FC) > STRONG_LOG2FC_THRESHOLD &
    marker_genes$p_val_adj < PADJ_THRESHOLD &
    marker_genes$unity_score == marker_genes$n_donors_tested
}

#' Write the full marker table and, if non-empty, the "robust" subset to .xlsx.
write_deg_outputs <- function(marker_genes, output_prefix) {
  openxlsx::write.xlsx(marker_genes, paste0(output_prefix, "_all.xlsx"))
  robust <- marker_genes[is_robust_deg(marker_genes), ]
  if (nrow(robust) > 0) {
    openxlsx::write.xlsx(robust, paste0(output_prefix, "_robust.xlsx"))
  }
}

#' Volcano plot: log2 fold-change vs. -log10(adjusted p-value), colored by
#' regulation category, with either the top-N hits or a specific gene set
#' labeled.
#'
#' @param marker_genes Categorized marker table.
#' @param highlight_genes If given, label exactly these genes (when present).
#'   Otherwise, label the `label_top_n` genes with the largest
#'   |log2FC| x -log10(padj) product.
#' @param label_top_n Number of top hits to label when `highlight_genes` is NULL.
#' @param title Plot title.
make_volcano_plot <- function(marker_genes, highlight_genes = NULL, label_top_n = 5, title = "") {
  category_colors <- c(upregulated = "#66c2a5", downregulated = "#fc8d62", "below threshold" = "grey")
  category_labels <- c("upregulated", "downregulated", "below\nthreshold")

  if (is.null(highlight_genes)) {
    ranking <- order(abs(marker_genes$avg_log2FC) * -log10(marker_genes$p_val_adj), decreasing = TRUE)
    label_data <- marker_genes[ranking[seq_len(min(label_top_n, nrow(marker_genes)))], ]
  } else {
    label_data <- marker_genes[marker_genes$gene %in% highlight_genes, ]
  }

  ggplot(marker_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), fill = category)) +
    scale_color_manual(values = category_colors, breaks = names(category_colors), labels = category_labels) +
    scale_fill_manual(values = category_colors, breaks = names(category_colors), labels = category_labels) +
    geom_point(shape = 21, alpha = 0.5) +
    geom_hline(yintercept = -log10(PADJ_THRESHOLD), linetype = "dashed", alpha = 0.2) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.2) +
    labs(x = "log2(fold-change)", y = "-log10(adjusted P-value)", title = title) +
    guides(color = "none", alpha = "none", fill = "none") +
    geom_label_repel(data = label_data, mapping = aes(alpha = p_val_adj < PADJ_THRESHOLD, label = gene))
}

#' Format an `enrichGO()` result into a plottable data frame and dot plot.
#'
#' @param ego An `enrichResult` from `clusterProfiler::enrichGO()`.
#' @param max_char Maximum term-label length before truncating with "...".
#' @param stitle Plot subtitle (e.g. "Upregulated" / "Downregulated").
#' @return A list of `list(ego, plot)`: the (annotated) `ego` object and the ggplot.
process_ego <- function(ego, max_char = 50, stitle) {
  ego@result$ratio <- vapply(ego@result$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/", fixed = TRUE)[[1]])
    parts[1] / parts[2]
  }, numeric(1))
  ego@result$Count <- vapply(ego@result$GeneRatio, function(x) {
    as.integer(strsplit(x, "/", fixed = TRUE)[[1]][1])
  }, integer(1))

  to_plot <- ego@result %>%
    slice_max(n = 10, order_by = -(p.adjust), with_ties = FALSE) %>%
    mutate(term = paste0(ID, ": ", Description)) %>%
    mutate(term = ifelse(nchar(term) > max_char, paste0(substr(term, 1, max_char), "..."), term))
  to_plot$unique_term <- paste(to_plot$ONTOLOGY, to_plot$term, sep = ": ")

  p <- to_plot %>%
    ggplot(aes(x = ratio, y = reorder(unique_term, ratio),
               size = factor(Count, levels = min(Count):max(Count)), fill = p.adjust)) +
    geom_point(shape = 21) +
    scale_fill_distiller(palette = "PuBuGn", direction = -1) +
    scale_y_discrete(labels = function(x) gsub(".*? (GO:\\d+): ", "\\1 ", x)) + # drop ontology prefix
    facet_nested(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    xlab("Gene ratio") + ylab("Ontology term") +
    labs(subtitle = stitle, size = "Count") +
    guides(size = guide_legend(order = 1), fill = guide_legend(order = 2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list(ego, p)
}

#' Background gene universe for GO enrichment: genes detected in at least
#' `GO_MIN_PCT_DETECTED`% of cells in the given Seurat object.
gene_universe <- function(seurat_obj) {
  counts <- seurat_obj@assays$RNA@counts
  pct_detected <- rowSums(counts > 0) / ncol(counts) * 100
  rownames(counts)[pct_detected >= GO_MIN_PCT_DETECTED]
}

#' Save one direction's (up- or down-regulated) enrichGO() result: an .xlsx
#' of the enrichment table and a dot plot, via `process_ego()`. No-op if
#' `ego` is NULL or has no enriched terms.
save_ego_result <- function(ego, direction_label, output_prefix) {
  if (is.null(ego) || nrow(ego) == 0) return(invisible(NULL))
  pathways_dir <- file.path(dirname(output_prefix), "pathways")
  dir.create(pathways_dir, recursive = TRUE, showWarnings = FALSE)
  result <- process_ego(ego, stitle = direction_label)
  base <- file.path(pathways_dir, paste0(basename(output_prefix), "_GO_", tolower(direction_label), "_enriched_terms"))
  openxlsx::write.xlsx(result[[1]]@result, paste0(base, ".xlsx"))
  ggsave(paste0(base, ".png"), result[[2]], width = 7, height = 7, units = "in", dpi = 600)
}

# --- DEG summary table (one row per contrast, one table per cell type) -------

DEG_SUMMARY_COLNAMES <- c(
  "Contrast",
  "Base filter (|logFC| > 0.25)",
  "Differentially expressed (Padj < 0.05)",
  "Upregulated (Padj < 0.05)",
  "Downregulated (Padj < 0.05)",
  sprintf("Strongly upregulated (Padj < 0.05 & logFC > %.2f)", STRONG_LOG2FC_THRESHOLD),
  sprintf("Strongly downregulated (Padj < 0.05 & logFC < -%.2f)", STRONG_LOG2FC_THRESHOLD)
)

make_empty_deg_summary <- function() {
  data.frame(contrast = character(0), base_filter = numeric(0), DEG = numeric(0),
             upreg = numeric(0), downreg = numeric(0), strong_upreg = numeric(0), strong_downreg = numeric(0))
}

add_deg_summary_row <- function(table_of_degs, contrast_name, marker_genes) {
  table_of_degs[nrow(table_of_degs) + 1, ] <- list(
    contrast_name,
    nrow(marker_genes),
    sum(marker_genes$significant),
    sum(marker_genes$category == "upregulated"),
    sum(marker_genes$category == "downregulated"),
    sum(marker_genes$significant & marker_genes$avg_log2FC >  STRONG_LOG2FC_THRESHOLD),
    sum(marker_genes$significant & marker_genes$avg_log2FC < -STRONG_LOG2FC_THRESHOLD)
  )
  table_of_degs
}

## =============================================================================
## 4. Load the annotated, integrated Seurat object
## =============================================================================

xm <- readRDS(SEURAT_OBJECT_PATH)

required_columns <- c("narrow_annot", "condition", "disease_status", "sample")
missing_columns <- setdiff(required_columns, colnames(xm@meta.data))
if (length(missing_columns) > 0) {
  stop(
    "The Seurat object at ", SEURAT_OBJECT_PATH, " is missing required metadata column(s): ",
    paste(missing_columns, collapse = ", "),
    ". Run the integration pipeline (which derives `narrow_annot` from `cluster_harmony`) first."
  )
}

## =============================================================================
## 5. Per-cell-type differential expression
## =============================================================================

for (ct in CELL_TYPES_OF_INTEREST) {
  message("== Cell type: ", ct, " ==")
  ct_dir <- file.path(DE_OUTPUT_DIR, ct)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  so <- subset(xm, subset = narrow_annot == ct)
  so$donor <- substr(so$sample, 1, 2) # subject code, e.g. "P1", "C2"

  openxlsx::write.xlsx(
    as.data.frame.matrix(table(so$disease_status, so$condition)),
    file.path(ct_dir, "cell_counts_by_group.xlsx"),
    rowNames = TRUE
  )

  table_of_degs <- make_empty_deg_summary()

  ## --- Contrast A: Patient vs. Control, within each treatment condition -----
  contrast_a_dir <- file.path(ct_dir, "withinCondition_diseaseStatusComparison")
  dir.create(contrast_a_dir, showWarnings = FALSE)

  for (cond in unique(so$condition)) {
    message(" - ", cond, ": Patient vs Control")
    output_prefix <- file.path(contrast_a_dir, paste0(cond, "_disease_specific_DEG"))

    sso <- subset(so, subset = condition == cond)
    sso <- SetIdent(sso, value = "disease_status")
    marker_genes <- FindMarkers(sso, ident.1 = "Patient", ident.2 = "Control")
    marker_genes <- categorize_markers(marker_genes)
    marker_genes <- flag_pathway_genes(marker_genes)

    # Leave out one patient donor at a time.
    patient_donors <- unique(sso$donor[sso$disease_status == "Patient"])
    marker_genes <- compute_unity_score(
      marker_genes, sso, donors_to_test = patient_donors,
      donor_col = "donor", group_col = "disease_status", ident_1 = "Patient", ident_2 = "Control"
    )

    table_of_degs <- add_deg_summary_row(table_of_degs, paste0("within_", cond, "_Patient_vs_Control"), marker_genes)
    write_deg_outputs(marker_genes, output_prefix)

    p_top <- make_volcano_plot(marker_genes, title = paste0(ct, ": ", cond, ", Patient vs Control"))
    ggsave(paste0(output_prefix, "_volcano.png"), p_top, width = 7, height = 5, units = "in", dpi = 600)

    p_pathways <- make_volcano_plot(
      marker_genes, highlight_genes = c(NFKB_TARGET_GENES, NLRP3_TARGET_GENES),
      title = paste0(ct, ": ", cond, ", Patient vs Control")
    )
    ggsave(paste0(output_prefix, "_volcano_NFKB_NLRP3.png"), p_pathways, width = 7, height = 5, units = "in", dpi = 600)

    if (sum(marker_genes$significant) >= MIN_GENES_FOR_GO) {
      universe_gs <- gene_universe(sso)

      # Above MANY_SIG_GENES_THRESHOLD significant genes, tighten to the
      # strong effect-size cutoff so enrichment stays focused on the most
      # confident hits and computationally tractable.
      if (sum(marker_genes$significant) >= MANY_SIG_GENES_THRESHOLD) {
        up_genes_to_use <- marker_genes$category == "upregulated"   & marker_genes$avg_log2FC >  STRONG_LOG2FC_THRESHOLD
        do_genes_to_use <- marker_genes$category == "downregulated" & marker_genes$avg_log2FC < -STRONG_LOG2FC_THRESHOLD
      } else {
        up_genes_to_use <- marker_genes$category == "upregulated"
        do_genes_to_use <- marker_genes$category == "downregulated"
      }

      rm(ego_up, ego_do); gc()
      if (sum(up_genes_to_use) >= MIN_GENES_FOR_GO) {
        ego_up <- enrichGO(gene = marker_genes$gene[up_genes_to_use],
                            universe = universe_gs, keyType = GO_KEY_TYPE, OrgDb = OrgDb,
                            ont = GO_ONTOLOGY, pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = FALSE)
      }
      if (sum(do_genes_to_use) >= MIN_GENES_FOR_GO) {
        ego_do <- enrichGO(gene = marker_genes$gene[do_genes_to_use],
                            universe = universe_gs, keyType = GO_KEY_TYPE, OrgDb = OrgDb,
                            ont = GO_ONTOLOGY, pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = FALSE)
      }

      if (!exists("ego_up") && !exists("ego_do")) {
        cat("Not enough genes to derive significant pathway terms for celltype", ct, "within", cond, "between Patient and Controls\n")
      } else {
        if (exists("ego_do")) save_ego_result(ego_do, "Downregulated", output_prefix)
        if (exists("ego_up")) save_ego_result(ego_up, "Upregulated", output_prefix)
      }
    }
  }

  ## --- Contrast B: LPS / LPS-ATP vs. Untreated, within each disease status --
  contrast_b_dir <- file.path(ct_dir, "withinDiseaseStatus_conditionComparison")
  dir.create(contrast_b_dir, showWarnings = FALSE)

  for (dis in unique(so$disease_status)) {
    for (cond in c("LPS", "LPS-ATP")) {
      message(" - ", dis, " | ", cond, " vs Untreated")
      output_prefix <- file.path(contrast_b_dir, paste0("within_", dis, "_", cond, "_DEG"))

      sso <- subset(so, subset = disease_status == dis & condition %in% c(cond, "Untreated"))
      sso <- SetIdent(sso, value = "condition")
      marker_genes <- FindMarkers(sso, ident.1 = cond, ident.2 = "Untreated")
      marker_genes <- categorize_markers(marker_genes)
      marker_genes <- flag_pathway_genes(marker_genes)

      # Leave out one donor at a time 
      donors_in_group <- unique(sso$donor)
      marker_genes <- compute_unity_score(
        marker_genes, sso, donors_to_test = donors_in_group,
        donor_col = "donor", group_col = "condition", ident_1 = cond, ident_2 = "Untreated"
      )

      table_of_degs <- add_deg_summary_row(table_of_degs, paste0("within_", dis, "_", cond, "_vs_Untreated"), marker_genes)
      write_deg_outputs(marker_genes, output_prefix)

      p_top <- make_volcano_plot(marker_genes, title = paste0(ct, ": ", dis, ", ", cond, " vs Untreated"))
      ggsave(paste0(output_prefix, "_volcano.png"), p_top, width = 7, height = 5, units = "in", dpi = 600)

      p_pathways <- make_volcano_plot(
        marker_genes, highlight_genes = c(NFKB_TARGET_GENES, NLRP3_TARGET_GENES),
        title = paste0(ct, ": ", dis, ", ", cond, " vs Untreated")
      )
      ggsave(paste0(output_prefix, "_volcano_NFKB_NLRP3.png"), p_pathways, width = 7, height = 5, units = "in", dpi = 600)

      if (sum(marker_genes$significant) >= MIN_GENES_FOR_GO) {
        universe_gs <- gene_universe(sso)

        if (sum(marker_genes$significant) >= MANY_SIG_GENES_THRESHOLD) {
          up_genes_to_use <- marker_genes$category == "upregulated"   & marker_genes$avg_log2FC >  STRONG_LOG2FC_THRESHOLD
          do_genes_to_use <- marker_genes$category == "downregulated" & marker_genes$avg_log2FC < -STRONG_LOG2FC_THRESHOLD
        } else {
          up_genes_to_use <- marker_genes$category == "upregulated"
          do_genes_to_use <- marker_genes$category == "downregulated"
        }

        rm(ego_up, ego_do); gc()
        if (sum(up_genes_to_use) >= MIN_GENES_FOR_GO) {
          ego_up <- enrichGO(gene = marker_genes$gene[up_genes_to_use],
                              universe = universe_gs, keyType = GO_KEY_TYPE, OrgDb = OrgDb,
                              ont = GO_ONTOLOGY, pAdjustMethod = "BH",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = FALSE)
        }
        if (sum(do_genes_to_use) >= MIN_GENES_FOR_GO) {
          ego_do <- enrichGO(gene = marker_genes$gene[do_genes_to_use],
                              universe = universe_gs, keyType = GO_KEY_TYPE, OrgDb = OrgDb,
                              ont = GO_ONTOLOGY, pAdjustMethod = "BH",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = FALSE)
        }

        if (!exists("ego_up") && !exists("ego_do")) {
          cat("Not enough genes to derive significant pathway terms for celltype", ct, "within", dis, "for", cond, "vs Untreated\n")
        } else {
          if (exists("ego_do")) save_ego_result(ego_do, "Downregulated", output_prefix)
          if (exists("ego_up")) save_ego_result(ego_up, "Upregulated", output_prefix)
        }
      }
    }
  }

  colnames(table_of_degs) <- DEG_SUMMARY_COLNAMES
  openxlsx::write.xlsx(table_of_degs, file.path(ct_dir, "table_of_DEGs_summary.xlsx"))
  message("Finished cell type: ", ct)
}

message("Per-cell-type differential expression analysis complete.")

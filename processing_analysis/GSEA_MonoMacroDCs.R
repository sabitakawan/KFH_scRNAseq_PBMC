#!/usr/bin/env Rscript
# ==============================================================================
# GSEA (Gene Set Enrichment Analysis) -- Monocytes, DCs, Macrophages
#
# For the Monocytes/DCs/Macrophages compartment, this script runs pre-ranked
# GSEA (GO, KEGG, Reactome) for the same two contrast families used in the
# per-cell-type differential expression analysis:
#
#   (A) Patient vs. Control, within each treatment condition
#   (B) LPS / LPS-ATP vs. Untreated, within each disease status
#
# Unlike the (threshold-based) differential expression script, GSEA here is
# run on the *entire* ranked gene list (no log-fold-change or detection-rate
# filtering), since GSEA's power comes from using the full ranking rather
# than a pre-filtered gene set.
#
# For each contrast, it:
#   - ranks genes by log2 fold-change (Entrez IDs, via a cached biomaRt
#     symbol<->Entrez mapping table),
#   - runs gseGO(), gseKEGG(), and gsePathway() (Reactome), saving both the
#     raw result and a qvalue-filtered, symbol-remapped summary table,
#   - reduces redundant GO terms by semantic similarity (rrvgo) per
#     ontology (BP/MF/CC) and draws a similarity heatmap,
#   - draws "parental term" dot plots (GO parent terms; full filtered
#     tables for KEGG/Reactome, which have no redundancy-reduction step).
#
# It then produces a handful of specific cross-contrast comparison plots
# (e.g. "does this LPS-ATP response differ between Patient and Control")
# and gene-level violin plots for three marker panels tied to GO terms that
# came up as enriched in this analysis.
#
#
# Usage:
#   Rscript GSEA_MonoMacroDCs.R
# ==============================================================================

## =============================================================================
## 1. Configuration
## =============================================================================


BASE_DIR           <- "path_to_dir/KFH_scRNAseq/"
GSEA_OUTPUT_DIR     <- file.path(BASE_DIR, "tgi", "GSEA", "MonoMacroDC")
BIOM_CACHE_PATH     <- file.path(BASE_DIR, "tgi", "GSEA", "bioM_genesDataset.Robj")

SEURAT_OBJECT_PATH <- "path_to_dir/KFH_scRNAseq/checkpoints/seurat_final_annotated.rds"

CELL_TYPE_OF_INTEREST <- "Monocytes, DCs, Macrophages"

# GSEA parameters (shared by gseGO / gseKEGG / gsePathway)
GSEA_MIN_SET_SIZE <- 10
GSEA_MAX_SET_SIZE <- 500
GSEA_N_PERM       <- 10000
QVALUE_THRESHOLD  <- 0.05

# GSEA ranks the *entire* gene list, so FindMarkers is run without an
# effect-size or detection-rate pre-filter (unlike the DE script).
GSEA_FINDMARKERS_ARGS <- list(logfc.threshold = 0, min.pct = 0.1, only.pos = FALSE)

# rrvgo redundant-GO-term reduction
RRVGO_SIM_THRESHOLD <- 0.65
RRVGO_SIM_METHOD    <- "Rel"

# Gene panels for the violin-plot figures, one subdirectory per panel,
# named after the GO term that motivated it.
GENE_PANELS <- list(
  "cellular_response_to_interleukin-1" = c(
    "ZC3H12A", "CCL2", "CCL7", "CEBPB", "IRAK3", "PLCB1", "GBP3", "EDN1",
    "GBP1", "CCL4", "CCL24", "IL1B", "CCL3", "IL6", "CCL20", "CXCL8"
  ),
  "antigen_processing_and_presentation_of_peptide_antigen_via_MHC_class_II" = c(
    "HLA-DRB5", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "LGMN", "HLA-DPA1",
    "CD74", "HLA-DMA", "HLA-DRB1", "HLA-DQA1", "HLA-DRA"
  ),
  "response_to_interferon-alpha" = c(
    "IFITM3", "IFIT2", "GAS6", "LAMP3", "IFIT3", "IFNAR2", "EIF2AK2", "IFITM2"
  )
)


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
  library(biomaRt)
  library(ReactomePA)
  library(xlsx)
  library(rrvgo)
  library(GOSemSim)
})
# Not used by this script (dropped from the working draft): msigdbr,
# forcats, dittoSeq -- none of their functions are called anywhere below.

options(SEURAT_OPTIONS)
OrgDb <- org.Hs.eg.db

dir.create(GSEA_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## 3. Helper functions -- gene ID mapping
## =============================================================================

#' Rank genes by log2 fold-change for GSEA, keyed by Entrez ID.
#'
#' clusterProfiler's GSEA functions require a sorted, named numeric vector.
#' Genes with no Entrez mapping in `bioM` are dropped. When a symbol maps to
#' more than one Entrez ID, the first mapping is used (matching the original
#' script's `ifelse()`-based lookup, which implicitly does the same).
#'
#' @param marker_genes Output of `Seurat::FindMarkers()` (rownames = gene symbols).
#' @param bioM Data frame with `hgnc_symbol` and `entrezgene_id` columns.
rank_genes_for_gsea <- function(marker_genes, bioM) {
  logFC <- marker_genes$avg_log2FC
  gene_id <- vapply(rownames(marker_genes), function(sym) {
    if (sym %in% bioM$hgnc_symbol) as.character(bioM$entrezgene_id[bioM$hgnc_symbol == sym][1]) else NA_character_
  }, character(1))

  logFC <- logFC[!is.na(gene_id)]
  names(logFC) <- gene_id[!is.na(gene_id)]
  sort(logFC, decreasing = TRUE)
}

#' Map a "/"-separated string of Entrez IDs (as stored in an enrichResult's
#' `core_enrichment` column) back to gene symbols.
#'
#' @param core_enrichment Character vector of "/"-separated Entrez ID strings.
#' @param bioM Data frame with `hgnc_symbol` and `entrezgene_id` columns.
remap_core_enrichment_to_symbols <- function(core_enrichment, bioM) {
  Reduce("rbind", lapply(strsplit(core_enrichment, "/", fixed = TRUE), function(entrez_ids) {
    symbols <- vapply(entrez_ids, function(id) {
      if (id %in% bioM$entrezgene_id) as.character(bioM$hgnc_symbol[bioM$entrezgene_id == id][1]) else NA_character_
    }, character(1))
    paste(symbols, collapse = "/")
  }))
}

## =============================================================================
## 4. Helper functions -- running and saving GSEA
## =============================================================================

GSEA_METHOD_LABELS <- list(
  GO       = list(rds = "GO",       xlsx = "GO"),
  kegg     = list(rds = "kegg",     xlsx = "KEGG"),
  reactome = list(rds = "reactome", xlsx = "Reactome")
)

#' Build a clean summary table from a raw `gseaResult`, filter to
#' qvalue < `QVALUE_THRESHOLD`, remap `core_enrichment` Entrez IDs back to
#' gene symbols, sort by NES, and save (filtered RDS + .xlsx).
#'
#' @param gsea_result A `gseaResult` (from gseGO/gseKEGG/gsePathway).
#' @param bioM Entrez<->symbol mapping table.
#' @param method One of "GO", "kegg", "reactome" (selects output naming).
#' @param contrast_label Used to build output filenames (e.g. "Patient_vs_Control_LPS").
#' @param has_ontology If TRUE (GO only), also writes separate BP/MF/CC .xlsx sheets.
#' @return The filtered, remapped summary table (invisibly).
process_and_save_gsea <- function(gsea_result, bioM, method, contrast_label, has_ontology) {
  labels <- GSEA_METHOD_LABELS[[method]]

  cols <- c("ID", "Description", "setSize", "NES", "pvalue", "qvalue", "rank", "core_enrichment")
  if (has_ontology) cols <- append(cols, "ONTOLOGY", after = 1)
  tbl <- as.data.frame(gsea_result@result)[, cols]
  colnames(tbl)[colnames(tbl) == "qvalue"] <- "qvalues"

  tbl <- tbl[tbl$qvalues < QVALUE_THRESHOLD, ]
  tbl$core_enrichment <- remap_core_enrichment_to_symbols(tbl$core_enrichment, bioM)
  tbl <- tbl[order(tbl$NES, decreasing = TRUE), ]

  saveRDS(tbl, file.path(GSEA_OUTPUT_DIR, paste0("GSEA_", labels$rds, "_tableFiltered_", contrast_label, ".RDS")))

  xlsx_path <- file.path(GSEA_OUTPUT_DIR, paste0("GSEA_", labels$xlsx, "_tableFiltered_", contrast_label, ".xlsx"))
  if (has_ontology) {
    xlsx::write.xlsx(tbl[tbl$ONTOLOGY == "BP", ], xlsx_path, sheetName = "Biological Process")
    xlsx::write.xlsx(tbl[tbl$ONTOLOGY == "MF", ], xlsx_path, sheetName = "Molecular Function", append = TRUE)
    xlsx::write.xlsx(tbl[tbl$ONTOLOGY == "CC", ], xlsx_path, sheetName = "Cellular Component", append = TRUE)
  } else {
    xlsx::write.xlsx(tbl, xlsx_path)
  }
  invisible(tbl)
}

#' Run the full GSEA suite (GO, KEGG, Reactome) for one ranked gene list and
#' save all raw ("allTerms") and filtered results.
#'
#' @param logFC Named, sorted numeric vector (Entrez IDs -> log2FC), from `rank_genes_for_gsea()`.
#' @param contrast_label Used to build every output filename for this contrast.
#' @param bioM Entrez<->symbol mapping table.
run_gsea_for_contrast <- function(logFC, contrast_label, bioM) {
  message("  GSEA (GO): ", contrast_label)
  gsea_go <- gseGO(geneList = logFC, OrgDb = OrgDb, ont = "ALL",
                    minGSSize = GSEA_MIN_SET_SIZE, maxGSSize = GSEA_MAX_SET_SIZE,
                    pvalueCutoff = 1, verbose = FALSE, nPermSimple = GSEA_N_PERM, eps = 0)
  saveRDS(gsea_go, file.path(GSEA_OUTPUT_DIR, paste0("GSEA_GO_allTerms_", contrast_label, ".RDS")))
  process_and_save_gsea(gsea_go, bioM, "GO", contrast_label, has_ontology = TRUE)

  message("  GSEA (KEGG): ", contrast_label)
  gsea_kegg <- gseKEGG(geneList = logFC, organism = "hsa", keyType = "kegg",
                        minGSSize = GSEA_MIN_SET_SIZE, maxGSSize = GSEA_MAX_SET_SIZE,
                        pvalueCutoff = 1, verbose = FALSE, nPermSimple = GSEA_N_PERM, eps = 0)
  saveRDS(gsea_kegg, file.path(GSEA_OUTPUT_DIR, paste0("GSEA_kegg_allTerms_", contrast_label, ".RDS")))
  process_and_save_gsea(gsea_kegg, bioM, "kegg", contrast_label, has_ontology = FALSE)

  message("  GSEA (Reactome): ", contrast_label)
  gsea_reactome <- ReactomePA::gsePathway(geneList = logFC, organism = "human",
                        minGSSize = GSEA_MIN_SET_SIZE, maxGSSize = GSEA_MAX_SET_SIZE,
                        pvalueCutoff = 1, verbose = FALSE, nPermSimple = GSEA_N_PERM, eps = 0)
  saveRDS(gsea_reactome, file.path(GSEA_OUTPUT_DIR, paste0("GSEA_reactome_allTerms_", contrast_label, ".RDS")))
  process_and_save_gsea(gsea_reactome, bioM, "reactome", contrast_label, has_ontology = FALSE)

  invisible(NULL)
}

## =============================================================================
## 5. Helper functions -- redundant GO term reduction and dot plots
## =============================================================================

#' Group a set of enriched GO terms (one ontology) by semantic similarity,
#' save a similarity heatmap and the reduced "parent term" table, and return
#' the unique parent-term descriptions.
#'
#' Positive- and negative-NES terms are reduced separately, since up- and
#' down-regulated enrichment are usually driven by different biological
#' processes and mixing them would distort the similarity clustering.
#'
#' @param go_table Filtered GSEA GO table (qvalue < threshold), for ONE ontology.
#' @param semdata Precomputed `GOSemSim::godata()` object for this ontology.
#' @param ont Ontology code ("BP", "MF", or "CC").
#' @param name Label used to build this contrast/ontology's output filenames.
#' @return Character vector of unique parent-term GO descriptions (possibly length 0).
reduce_go_terms <- function(go_table, semdata, ont, name) {
  rownames(go_table) <- go_table$ID
  reduced_all <- NULL
  parent_terms <- character(0)

  reduce_one_direction <- function(direction_table, main_label) {
    if (nrow(direction_table) <= 1) return(NULL)
    sim_matrix <- calculateSimMatrix(direction_table$ID, semdata = semdata, orgdb = "org.Hs.eg.db",
                                      ont = ont, method = RRVGO_SIM_METHOD)
    scores <- setNames(abs(direction_table$NES), direction_table$ID)
    reduced <- reduceSimMatrix(sim_matrix, scores, threshold = RRVGO_SIM_THRESHOLD, orgdb = "org.Hs.eg.db")
    reduced <- reduced[order(reduced$cluster), ]
    reduced$parentTerm <- factor(reduced$parentTerm, levels = unique(reduced$parentTerm))
    print(heatmapPlot(sim_matrix[rownames(reduced), rownames(reduced)], reduced,
                       annotateParent = TRUE, annotationLabel = "parentTerm", fontsize = 6,
                       cluster_rows = FALSE, cluster_cols = FALSE, labels_row = NULL,
                       labels_col = reduced$term, main = main_label))
    reduced
  }

  pdf(file.path(GSEA_OUTPUT_DIR, paste0("Heatmap_SemSim_scores", name, ".pdf")), width = 15, height = 15)
  pos_reduced <- reduce_one_direction(go_table[go_table$NES > 0, ], "NES > 0")
  if (!is.null(pos_reduced)) {
    reduced_all <- rbind(reduced_all, pos_reduced)
    parent_terms <- c(parent_terms, unique(pos_reduced$parentTerm))
  }
  neg_reduced <- reduce_one_direction(go_table[go_table$NES < 0, ], "NES < 0")
  if (!is.null(neg_reduced)) {
    reduced_all <- rbind(reduced_all, neg_reduced)
    parent_terms <- c(parent_terms, unique(neg_reduced$parentTerm))
  }
  dev.off()

  if (!is.null(reduced_all)) {
    xlsx::write.xlsx(reduced_all[, -6], file.path(GSEA_OUTPUT_DIR, paste0("reducedTerms_table", name, ".xlsx")))
  }

  unique(parent_terms)
}

#' "Parental term" GSEA dot plot: gene ratio (x) vs. pathway (y), colored by
#' NES, sized by -log10(qvalue), faceted by regulation direction.
#'
#' @param go_subset A filtered GSEA table (GO/KEGG/Reactome), restricted to
#'   the terms to plot (e.g. parent terms for GO; the full filtered table
#'   for KEGG/Reactome, which have no redundancy-reduction step).
#' @param out_path PDF path to save to.
make_parental_terms_dotplot <- function(go_subset, out_path) {
  if (nrow(go_subset) == 0) return(invisible(NULL))

  df <- data.frame(
    pathways = go_subset$Description,
    NES = go_subset$NES,
    qvalue = go_subset$qvalues,
    Gene_Ratio = vapply(go_subset$core_enrichment, function(x) length(strsplit(x, "/", fixed = TRUE)[[1]]), integer(1)) / go_subset$setSize
  )
  df$type <- ifelse(df$NES < 0, "downregulated", "upregulated")
  df$log_qvalue <- -log10(df$qvalue)
  df$Gene_Ratio <- round(df$Gene_Ratio * 100)
  if (nrow(df) == 0) return(invisible(NULL))

  df <- df %>% arrange(type, Gene_Ratio)
  df$pathways <- factor(df$pathways, levels = df$pathways)

  n_types <- length(unique(df$type))
  plot_width <- if (n_types == 2) {
    (8 / 50) * max(nchar(as.character(df$pathways))) + 4
  } else {
    (4 / 50) * max(nchar(as.character(df$pathways))) + 4
  }
  plot_height <- (5 / 20) * nrow(df) + 4

  pdf(out_path, width = plot_width, height = plot_height)
  print(
    ggplot(df, aes(x = Gene_Ratio, y = pathways)) +
      geom_point(aes(size = log_qvalue, color = NES)) +
      theme_bw(base_size = 14) +
      scale_size_continuous(limits = c(1.3, max(df$log_qvalue))) +
      scale_colour_gradient2(midpoint = 0, low = "darkslateblue", mid = "white", high = "darkred") +
      ylab(NULL) + ggtitle("GSEA") + facet_grid(. ~ type)
  )
  dev.off()
}

#' For one contrast: reduce redundant GO terms and draw a parental-term dot
#' plot per ontology (BP/MF/CC), then draw plain parental-term dot plots for
#' KEGG and Reactome (no redundancy-reduction step for these).
#'
#' @param contrast_label Identifies which saved GSEA result files to read back.
#' @param go_semdata Named list of precomputed `GOSemSim::godata()` objects, one per ontology.
summarize_gsea_terms <- function(contrast_label, go_semdata) {
  go_table <- readRDS(file.path(GSEA_OUTPUT_DIR, paste0("GSEA_GO_tableFiltered_", contrast_label, ".RDS")))

  for (ont in c("BP", "MF", "CC")) {
    ont_table <- go_table[go_table$ONTOLOGY == ont, ]
    if (nrow(ont_table) == 0) next
    name <- paste0("GSEA_GO_", ont, "_", contrast_label)

    parent_terms <- reduce_go_terms(ont_table, go_semdata[[ont]], ont, name)
    dotplot_table <- ont_table[ont_table$Description %in% parent_terms, ]
    make_parental_terms_dotplot(dotplot_table, file.path(GSEA_OUTPUT_DIR, paste0("DotPlot_ParentalTerms_", name, ".pdf")))
  }

  kegg_name <- paste0("GSEA_KEGG_", contrast_label)
  kegg_table <- readRDS(file.path(GSEA_OUTPUT_DIR, paste0("GSEA_kegg_tableFiltered_", contrast_label, ".RDS")))
  make_parental_terms_dotplot(kegg_table, file.path(GSEA_OUTPUT_DIR, paste0("DotPlot_ParentalTerms_", kegg_name, ".pdf")))

  reactome_name <- paste0("GSEA_Reactome_", contrast_label)
  reactome_table <- readRDS(file.path(GSEA_OUTPUT_DIR, paste0("GSEA_reactome_tableFiltered_", contrast_label, ".RDS")))
  make_parental_terms_dotplot(reactome_table, file.path(GSEA_OUTPUT_DIR, paste0("DotPlot_ParentalTerms_", reactome_name, ".pdf")))

  invisible(NULL)
}

## =============================================================================
## 6. Helper functions -- cross-contrast comparison plots
## =============================================================================

#' Read a saved "reduced GO terms" table and return its unique parent GO IDs.
load_parent_term_ids <- function(path) {
  unique(xlsx::read.xlsx(path, sheetIndex = 1)$parent)
}

#' Load a saved GSEA "all terms" result and pull NES / qvalue / Description
#' for a specific set of Biological Process GO term IDs.
load_go_bp_terms <- function(rds_path, term_ids) {
  result <- readRDS(rds_path)@result
  result <- result[result$ONTOLOGY == "BP", ]
  result[term_ids, c("NES", "qvalue", "Description")]
}

#' Grouped bar plot of NES per pathway, one bar per `disease_condition`.
barplot_NES_grouped <- function(df, out_path, width = 15, height = 9) {
  df$disease_condition <- factor(df$disease_condition, levels = c("Patient", "Control"))
  pdf(out_path, width = width, height = height)
  print(
    ggplot(data = df, aes(x = NES, y = Description, fill = disease_condition)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(panel.background = element_rect(fill = "white", colour = "black"))
  )
  dev.off()
}

#' Bar plot of NES per pathway (single series, no grouping).
barplot_NES <- function(df, out_path, width = 15, height = 9) {
  pdf(out_path, width = width, height = height)
  print(
    ggplot(data = df, aes(x = NES, y = Description)) +
      geom_bar(stat = "identity") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(panel.background = element_rect(fill = "white", colour = "black"))
  )
  dev.off()
}

## =============================================================================
## 7. Load data, map gene IDs, and prepare the cell-type subset
## =============================================================================

xm <- readRDS(SEURAT_OBJECT_PATH)

required_columns <- c("narrow_annot", "condition", "disease_status")
missing_columns <- setdiff(required_columns, colnames(xm@meta.data))
if (length(missing_columns) > 0) {
  stop(
    "The Seurat object at ", SEURAT_OBJECT_PATH, " is missing required metadata column(s): ",
    paste(missing_columns, collapse = ", "),
    ". Run the integration pipeline first."
  )
}

# Entrez <-> HGNC symbol mapping, used to rank genes for GSEA and to remap
# `core_enrichment` gene lists back to symbols. Re-fetched (and the cache
# file overwritten) on every run, matching the working draft.
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
bioM <- getBM(filters = "hgnc_symbol", values = rownames(xm),
               attributes = c("entrezgene_id", "hgnc_symbol"), mart = mart)
save(bioM, file = BIOM_CACHE_PATH)

# Remove ribosomal and mitochondrial genes before GSEA (their strong,
# non-specific signal otherwise dominates enrichment results).
mt_genes <- grep("^MT-", rownames(xm@assays$SCT@counts), value = TRUE)
ribo_genes <- grep("^RP[SL]", rownames(xm@assays$SCT@counts), value = TRUE)
xm <- subset(xm, features = setdiff(rownames(xm@assays$SCT@counts), c(mt_genes, ribo_genes)))

so <- subset(xm, subset = narrow_annot == CELL_TYPE_OF_INTEREST)
# Re-normalize after removing ribosomal/mitochondrial genes -- SCTransform's
# regularized regression should be refit on the trimmed feature set rather
# than reused from before gene removal.
so <- SCTransform(so, method = "glmGamPoi", verbose = FALSE, vars.to.regress = "percent.mt")

## =============================================================================
## 8. Phase 1 -- run and save GSEA for every contrast
## =============================================================================

message("== Phase 1: running GSEA ==")

## --- Contrast A: Patient vs. Control, within each treatment condition -----
for (cond in unique(so$condition)) {
  message("Patient vs Control within ", cond)
  sso <- subset(so, subset = condition == cond)
  sso <- SetIdent(sso, value = "disease_status")
  marker_genes <- do.call(FindMarkers, c(list(object = sso, ident.1 = "Patient", ident.2 = "Control"), GSEA_FINDMARKERS_ARGS))

  logFC <- rank_genes_for_gsea(marker_genes, bioM)
  run_gsea_for_contrast(logFC, paste0("Patient_vs_Control_", cond), bioM)
}

## --- Contrast B: LPS / LPS-ATP vs. Untreated, within each disease status --
for (dis in unique(so$disease_status)) {
  for (cond in c("LPS", "LPS-ATP")) {
    message(cond, " vs Untreated within ", dis)
    sso <- subset(so, subset = disease_status == dis & condition %in% c(cond, "Untreated"))
    sso <- SetIdent(sso, value = "condition")
    marker_genes <- do.call(FindMarkers, c(list(object = sso, ident.1 = cond, ident.2 = "Untreated"), GSEA_FINDMARKERS_ARGS))

    logFC <- rank_genes_for_gsea(marker_genes, bioM)
    run_gsea_for_contrast(logFC, paste0(cond, "_vs_Untreated_", dis), bioM)
  }
}

## =============================================================================
## 9. Phase 2 -- redundant-term reduction and parental-term dot plots
## =============================================================================

message("== Phase 2: reducing GO terms and building dot plots ==")

# Computed once and reused across every contrast (it depends only on the
# ontology and organism database, not on any contrast's results).
GO_SEMDATA <- list(
  BP = GOSemSim::godata(org.Hs.eg.db, ont = "BP"),
  MF = GOSemSim::godata(org.Hs.eg.db, ont = "MF"),
  CC = GOSemSim::godata(org.Hs.eg.db, ont = "CC")
)

for (cond in unique(so$condition)) {
  summarize_gsea_terms(paste0("Patient_vs_Control_", cond), GO_SEMDATA)
}
for (dis in unique(so$disease_status)) {
  for (cond in c("LPS", "LPS-ATP")) {
    summarize_gsea_terms(paste0(cond, "_vs_Untreated_", dis), GO_SEMDATA)
  }
}

## =============================================================================
## 10. Cross-contrast comparison plots
## =============================================================================
# These reproduce a small set of specific, hand-picked comparisons (not a
# systematic sweep over every possible pairing) 

message("== Cross-contrast comparison plots ==")

## --- Does the LPS-ATP response (vs. Untreated) within Controls differ from
##     the same response within Patients? -------------------------------------
{
  term_ids <- load_parent_term_ids(file.path(GSEA_OUTPUT_DIR, "reducedTerms_tableGSEA_GO_BP_LPS-ATP_vs_Untreated_Control.xlsx"))
  control_rt <- load_go_bp_terms(file.path(GSEA_OUTPUT_DIR, "GSEA_GO_allTerms_LPS-ATP_vs_Untreated_Control.RDS"), term_ids)
  control_rt$disease_condition <- "Control"
  patient_rt <- load_go_bp_terms(file.path(GSEA_OUTPUT_DIR, "GSEA_GO_allTerms_LPS-ATP_vs_Untreated_Patient.RDS"), term_ids)
  patient_rt$disease_condition <- "Patient"

  rt_pos <- rbind(control_rt[control_rt$NES > 0, ], patient_rt[patient_rt$NES > 0, ])
  barplot_NES_grouped(rt_pos, file.path(GSEA_OUTPUT_DIR, "Barplot_reducedTerms_LPS-ATP_vs_Untreated_positive_withinControl_comparison_with_Patient.pdf"))

  rt_neg <- rbind(control_rt[control_rt$NES < 0, ], patient_rt[patient_rt$NES < 0, ])
  barplot_NES_grouped(rt_neg, file.path(GSEA_OUTPUT_DIR, "Barplot_reducedTerms_LPS-ATP_vs_Untreated_negative_withinControl_comparison_with_Patient.pdf"))
}

## --- Are the Control-driven treatment-response terms also seen in the
##     Patient-vs-Control comparison (within LPS / LPS-ATP)? -----------------
for (cond in c("LPS", "LPS-ATP")) {
  term_ids <- load_parent_term_ids(file.path(GSEA_OUTPUT_DIR, paste0("reducedTerms_tableGSEA_GO_BP_", cond, "_vs_Untreated_Control.xlsx")))
  control_result <- readRDS(file.path(GSEA_OUTPUT_DIR, paste0("GSEA_GO_allTerms_", cond, "_vs_Untreated_Control.RDS")))@result
  tmp <- control_result[term_ids, ]
  pos_ids <- rownames(tmp[tmp$NES > 0, ])
  neg_ids <- rownames(tmp[tmp$NES < 0, ])

  go_table <- readRDS(file.path(GSEA_OUTPUT_DIR, paste0("GSEA_GO_tableFiltered_Patient_vs_Control_", cond, ".RDS")))
  rownames(go_table) <- go_table$ID

  pos_data <- go_table[pos_ids[pos_ids %in% rownames(go_table)], ]
  barplot_NES(pos_data, file.path(GSEA_OUTPUT_DIR, paste0("Barplot_reducedTerms_", cond, "_vs_Untreated_positive_in_Patient_vs_Control_GSEA.pdf")), width = 15, height = 9)

  neg_data <- go_table[neg_ids[neg_ids %in% rownames(go_table)], ]
  barplot_NES(neg_data, file.path(GSEA_OUTPUT_DIR, paste0("Barplot_reducedTerms_", cond, "_vs_Untreated_negative_in_Patient_vs_Control_GSEA.pdf")), width = 10, height = 3)
}

## --- Do the terms significant in the Patient-vs-Control comparison (within
##     LPS / LPS-ATP) show a differential treatment response between Patient
##     and Control? ------------------------------------------------------------
for (cond in c("LPS", "LPS-ATP")) {
  term_ids <- load_parent_term_ids(file.path(GSEA_OUTPUT_DIR, paste0("reducedTerms_tableGSEA_GO_BP_Patient_vs_Control_", cond, ".xlsx")))
  control_rt <- load_go_bp_terms(file.path(GSEA_OUTPUT_DIR, paste0("GSEA_GO_allTerms_", cond, "_vs_Untreated_Control.RDS")), term_ids)
  control_rt$disease_condition <- "Control"
  patient_rt <- load_go_bp_terms(file.path(GSEA_OUTPUT_DIR, paste0("GSEA_GO_allTerms_", cond, "_vs_Untreated_Patient.RDS")), term_ids)
  patient_rt$disease_condition <- "Patient"

  rt <- rbind(control_rt, patient_rt)
  rt <- rt[!is.na(rt$Description), ]
  barplot_NES_grouped(rt, file.path(GSEA_OUTPUT_DIR, paste0("Barplot_reducedTerms_Patient_vs_Control_within", cond, "_", cond, "induction.pdf")))
}

## =============================================================================
## 11. Gene-panel violin plots
## =============================================================================
# Marker panels tied to GO terms that came up as enriched above.

message("== Gene-panel violin plots ==")

so$condition <- factor(so$condition, levels = c("Untreated", "LPS", "LPS-ATP"))

for (panel_name in names(GENE_PANELS)) {
  panel_dir <- file.path(GSEA_OUTPUT_DIR, paste0("Gene_plot_", panel_name))
  dir.create(panel_dir, showWarnings = FALSE)
  for (gene in GENE_PANELS[[panel_name]]) {
    pdf(file.path(panel_dir, paste0("vln_plot_MonoMacroDC_condition_", gene, ".pdf")))
    print(
      VlnPlot(so, gene, pt.size = 0.5, group.by = "condition", split.by = "disease_status", raster = FALSE, assay = "SCT") +
        ggtitle(paste0(gene, " - Monocytes, Macrophages and DCs"))
    )
    dev.off()
  }
}

message("GSEA analysis complete.")

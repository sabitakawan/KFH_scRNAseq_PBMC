The dataset deposited here corresponds to the matrix of counts (raw counts and SCT-transformed) from the single-cell RNA-seq profiling of peripheral blood mononuclear cells (PBMCs) derived from healthy individuals (n=3) and patients diagnosed with Keratitis Fugax Hereditaria (KFH, n=3). Cells derived from these six individuals were subjected to three different conditions, leading to a total of 18 samples in the dataset.

Untreated ("unstimulated")
Primed with lipopolysaccharides ("LPS")
Activated with lipopolysaccharides and adenosinetriphosphate ("LPS+ATP").

The object is a Seurat object saved with a ".rds" extension in Zenodo as "Seurat_object_harmony_annotated.RDS"
md5code: f25adcd4301b58e91f93ea7b979a3468

Title of the manuscript: Single-cell RNA sequencing Reveals a Tolerance-like Monocyte Signature in NLRP3-asscoiated Keratitis Fugax Hereditaria


| R Script                                    | Figures                             |
| :------------------------------------------ | :---------------------------------- |
| `UMAP_and_Cellnumber.R`                     | Figure 1A, 1B, 1C, 2A, 2B, S1A, S1B |
| `Differential abundance.R`                  | Figure 1D, Figure S1C               |
| `Plot_inflammasome.R`                       | Figure 2D                           |
| `rrho_plots.R`                              | Figure 3A, 3B, 3C (left panel)      |
| `Signatures_withinTreatments_log2fc_plots.R`| Figure 3A, 3B, 3C (right panel)     |
| `ORA_monocytes_GO_terms.R`                  | Figure 4B, 4C, 5B, 6B               |
| `GSEA_GO_BP_KEGG_MONO_Unstimulated.R`       | Figure 4C, 4E                       |
| `LPS_monocytes_KEGG_ORA.R`                  | Figure 5C                           |
| `T_ORA_GO_BO_Semantic.R`                    | Figure 4D, 5D, 6C                   |
| `Within_HC_KFH_GO_terms_Semantic.R`         | Figures S4B, S4C, S4D, S4E          |
| `Pca_vca_gostats_oob_monocytes.R`           | Figures S2A, S2B, S2C               |


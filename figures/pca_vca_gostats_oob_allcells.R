library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggrepel)


obj <- readRDS("saved/ZenodoObjs/Seurat_object_harmony_annotated.RDS")

tmp <- as.data.frame(table(obj$sample))
colnames(tmp) <- c("Sample", "Ncells")
tmp$Sample <- gsub(" ","_", tmp$Sample)

tmp$Sample <- factor(tmp$Sample, levels=c(tmp$Sample[grepl("Untreated",tmp$Sample) & grepl("^C",tmp$Sample)],
                                          tmp$Sample[grepl("_LPS$",tmp$Sample) & grepl("^C",tmp$Sample)],
                                          tmp$Sample[grepl("_ATP$",tmp$Sample) & grepl("^C",tmp$Sample)],
                                          tmp$Sample[grepl("Untreated",tmp$Sample) & grepl("^P",tmp$Sample)],
                                          tmp$Sample[grepl("_LPS$",tmp$Sample) & grepl("^P",tmp$Sample)],
                                          tmp$Sample[grepl("_ATP$",tmp$Sample) & grepl("^P",tmp$Sample)]))

library(scales)

cellnum <- ggplot(tmp, aes(x = Sample, y = Ncells)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Ncells), vjust = 0.5, hjust=-0.25, angle=90, size=5) +
  theme_classic() +
  xlab("") +
  ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title.y=element_text(size=16))+
  scale_y_continuous(labels = label_comma(), limits=c(0,10000)) 

plot(cellnum)



## filter-out lowly expressed genes
selected_f <- names(which(rowSums(obj@assays$SCT["counts"])>0))
obj_all <- subset(obj, features = selected_f)

obj_all$sampleId <- gsub(" ","_", obj_all$sample)
samples_to_discard <- names(which(table(obj_all$sampleId)<10))

obj_all <- obj_all[,!obj_all$sampleId %in% samples_to_discard]

## Pseudobulk by sample

pseudo_obj_all <- AggregateExpression(
  obj_all,
  group.by = "sampleId",  # Replace with your grouping variable (e.g., sample, condition)
  assays = "SCT",
  slot = "counts",  # Raw counts
  return.seurat=T
)

## Find highly variable genes
pseudo_obj_all <- FindVariableFeatures(pseudo_obj_all)

## Make sure to have the sample ids in the metadata
pseudo_obj_all$disease_status <- ifelse(gsub("[0-9]-.+","",pseudo_obj_all@meta.data$orig.ident)=="C", "Control", "Patient")
pseudo_obj_all$treatment <- sapply(strsplit(pseudo_obj_all@meta.data$orig.ident,"-"), function(x) paste(x[-1], collapse="_"))

matrixForPca <- prcomp(pseudo_obj_all@assays$SCT["scale.data"][VariableFeatures(pseudo_obj_all),], scale = FALSE)
dfPca <- as.data.frame(matrixForPca$rotation)
stopifnot(all(stopifnot(rownames(pseudo_obj_all@meta.data)==rownames(dfPca))))

dfPca <- cbind(pseudo_obj_all@meta.data, dfPca)
dfPca$sampleId <- factor(dfPca$sampleId, levels=unique(sort(dfPca$sampleId)))

pseudo_obj_all$donor <- gsub("-.+","", pseudo_obj_all@meta.data$orig.ident)

## Just to know how much the principle components of each PC explain (since we will plot just the first two)
explained_variance <- matrixForPca$sdev^2 / sum(matrixForPca$sdev^2)
explained_variance_non_scientific <- format(round(explained_variance*100,2), scientific = FALSE)
names(explained_variance_non_scientific) <- paste0("PC", 1:length(explained_variance_non_scientific))
print(explained_variance_non_scientific)


## Get colors for the samples (I think in your case they are 12)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(pseudo_obj_all$treatment))
colVec <- setNames(getPalette(colourCount),
                   unique(pseudo_obj_all$treatment))

dfPca$ind <- gsub("-.+","",dfPca$sampleId)


##Plot just with two principal components: PC1-PC2
pseudoPCA_all <- ggplot(dfPca, aes(x=PC1, y=PC2, col=treatment, shape=disease_status, label=ind))+
  geom_point(size=5)+
  theme_bw()+
  scale_color_manual(name="Samples", values=colVec)+
  theme(plot.title=element_text(hjust=0.5, face="bold", size=14),
        legend.position="right")+
  ggtitle(paste0("OT Pseudobulk scaled HVG before gsva"))+
  xlab(paste0("PC1, var.explained = ", explained_variance_non_scientific["PC1"], " %"))+
  ylab(paste0("PC2, var.explained = ", explained_variance_non_scientific["PC2"], " %"))+
  geom_text_repel()

plot(pseudoPCA_all)

## Main Figure 3A ##


## To save plot as PDF
pdf(file="figs/mainFigures/Figure3A.pdf")
plot(pseudoPCA_all)
dev.off()


## grid with PCA correlation

pc_pairs <- combn(paste0("PC", 1:6), 2, simplify = FALSE)

library(dplyr)
library(purrr)

plot_df <- map_dfr(pc_pairs, function(pcs) {
  dfPca %>%
    transmute(
      PCx = .data[[pcs[1]]],
      PCy = .data[[pcs[2]]],
      PCx_name = pcs[1],
      PCy_name = pcs[2]
    )
})

plot_df$sampleId <- gsub("\\..+","",rownames(plot_df)) 
plot_df$disease_status <- dfPca[match(plot_df$sampleId, dfPca$sampleId),]$disease_status
plot_df$treatment <- dfPca[match(plot_df$sampleId, dfPca$sampleId),]$treatment


gridpca_first6 <- ggplot(plot_df, aes(PCx, PCy, col=treatment, shape=disease_status)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_grid(
    PCy_name ~ PCx_name,
    scales = "free",
    drop = TRUE
  ) +
  theme_bw()+
  theme(axis.text.x=element_text(hjust=1, vjust=0.5, size=8, angle=90),
        axis.text.y=element_text(size=8))

plot(gridpca_first6)


######################
## variancePartiton ##
######################

library(variancePartition)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(sva)

#### differential expression

se <- SummarizedExperiment(assays=list(counts=as.matrix(pseudo_obj_all@assays$SCT["counts"])),
                           colData=pseudo_obj_all@meta.data)

dge <- DGEList(counts=assays(se)$counts, genes=as.data.frame(mcols(se)))
dim(dge)

assays(se)$logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=0.5)

## Check the distribution

avgexp <- rowMeans(assays(se)$logCPM)

mask <- avgexp > -2
dim(se)
se.filt <- se[mask, ]
dim(se.filt)

dim(dge)
dge.filt <- dge[mask, ]
dim(dge.filt)

## Copy without normalisation of dge.filt counts to use later for voomQualityWeights
dge.filt.noTMM <- dge.filt

## Renormalise again
assays(se.filt)$logCPM <- edgeR::cpm(dge.filt, normalized.lib.sizes=TRUE,
                                     log=TRUE, prior.count=0.5)

colData(se.filt)$disease_status <- factor(colData(se.filt)$disease_status, levels=c("Control","Patient"))
colData(se.filt)$treatment <- factor(colData(se.filt)$treatment, levels=c("Untreated","LPS","LPS_ATP"))
se.filt$treatment <- relevel(se.filt$treatment, ref="Untreated")


################
### with SVA ###
################

mod <- model.matrix(~ disease_status + treatment, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))

IQRs <- apply(assays(se.filt)$logCPM, 1, IQR)
sv <- sva(assays(se.filt)$logCPM[IQRs > quantile(IQRs, prob=0.9), ],
          mod=mod, mod0=mod0)

cnames <- c(colnames(mod), paste0("SV", 1:sv$n))
mod <- cbind(mod, sv$sv)
colnames(mod) <- cnames
head(mod, n=3)

dge.filt <- estimateDisp(dge.filt, mod, robust=TRUE)
prdf <- cut(dge.filt$prior.df, breaks=c(0, 1, 2, 3, 4, 5))
table(prdf, useNA="always")

v <- voomWithQualityWeights(dge.filt.noTMM, mod, normalize="quantile")


## Run with voomWeights (recommended)
form <- ~  treatment + disease_status + SV1 + SV2 + SV3

colData(se.filt)$SV1 <- mod[, "SV1"]
colData(se.filt)$SV2 <- mod[, "SV2"]
colData(se.filt)$SV3 <- mod[, "SV3"]

varPart2 <- fitExtractVarPartModel(v, form, colData(se.filt))

## Check colinearity
res <- fitVarPartModel(v, form, colData(se.filt))
col_all <- sapply(1:length(res), function(x) colinearityScore(res[[x]]))
summary(col_all)

vp2 <- sortCols(varPart2)
pdf(file=paste0(dirQCplots,"variancePartition_voom_withSVA_all.pdf"))
plotVarPart(vp2)
dev.off()

mat_var2 <- as.data.frame(varPart2)
result <- data.frame(
  top_var = max_value <- apply(mat_var2, 1, max),
  top_factor = colnames(mat_var2)[max.col(mat_var2)]
)

disease_most_explained <- subset(result, top_factor=="disease_status")
disease_most_explained_ordered <- disease_most_explained[order(disease_most_explained$top_var, decreasing=T),]


tmp_var2 <- as.data.frame(table(result$top_factor))
colnames(tmp_var2) <- c("var","nGenes")

tmp_var2$var <- factor(tmp_var2$var, levels=tmp_var2$var[order(tmp_var2$nGenes, decreasing=F)])

## Main Figure 3D ##

pdf(file="figs/mainFigures/Figure3D.pdf", width=6.25, height=6.25)
ggplot(data=tmp_var2, aes(y=var, x=nGenes, fill=var))+
  geom_bar(stat="identity", color="black")+
  theme_bw()+
  scale_x_continuous(labels = scales::comma_format(), limits=c(0,11000), breaks=seq(0,10000,2500))+
  ylab("")+
  xlab("Genes with more VE")+
  geom_text(aes(x=nGenes, label=nGenes), hjust=-0.15, 
            color="black", size=4)+
  theme(axis.text=element_text(size=12.5, color="black"),
        axis.title=element_text(size=15),
        legend.position="none")+
  scale_fill_manual(values=c("SV3"="#8DD3C7",
                             "SV2"="#FFFFB3",
                             "SV1"="#BEBADA",
                             "disease_status"="#FB8072",
                             "Residuals"="#80B1D3",
                             "treatment"="#FDB462"))
dev.off()


### Run GO enrichment

library(GOstats)


GOenrichmentAndReport <- function(geneIds, universeGeneIds,
                                  minSize=3, maxSize=300, minCount=3,
                                  minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                  label="allDE"){
  
  
  GOparams <- new("GOHyperGParams", geneIds=geneIds,
                  universeGeneIds=universeGeneIds,
                  annotation="org.Hs.eg.db",
                  ontology="BP", pvalueCutoff=0.05, conditional=TRUE,
                  minSizeCutoff=3, maxSizeCutoff=300, orCutoff=1.5,
                  testDirection="over")
  
  hgOverGOBP <- hyperGTest(GOparams)
  
  report <- data.frame(GOBPID=as.vector(names(geneIdUniverse(hgOverGOBP))),
                       Pvalue=pvalues(hgOverGOBP),
                       OddsRatio=oddsRatios(hgOverGOBP),
                       ExpCount=expectedCounts(hgOverGOBP),
                       Count=geneCounts(hgOverGOBP),
                       Size=universeCounts(hgOverGOBP),
                       stringsAsFactors=FALSE)
  
  ## discard gene sets that do not meet a minimum and maximum number of genes
  report <- report[report$Size >= minSize & report$Size <= maxSize, , drop=FALSE]
  
  ## discard gene sets that show a p.value>0.05
  report <- report[report$Pvalue < p.value, , drop=FALSE]
  
  ## discard gene sets that do not satisfy the OR cutoff
  report <- report[report$OddsRatio >= minOddsRatio & report$Count >= minCount, , drop=FALSE]
  
  ## apply the maximum-number-of-GO-terms-reported cutoff and sort by odds ratio
  maxReported <- min(nrow(report))
  report <- report[sort(report$OddsRatio, decreasing=TRUE, index.return=TRUE)$ix[1:maxReported], ]
  
  if (dim(report[complete.cases(report),])[1]==0){
    message <- "No GO terms enriched"
    print(message)
    return(message)
  }
  
  ## add the symbol and GO term description, information
  reportGenes <- geneIdsByCategory(hgOverGOBP, report$GOBPID)
  reportGeneSyms <- lapply(reportGenes, annotate::getSYMBOL, "org.Hs.eg.db")
  highlightGeneMask <- lapply(reportGenes, function(x, hgenes, fmt) x %in% hgenes, highlightGenes)
  reportGeneSyms <- mapply(function(genes, mask, fmt) ifelse(mask, sprintf(fmt, genes), genes),
                           reportGeneSyms, highlightGeneMask, MoreArgs=list(fmt=highlightStr), SIMPLIFY=FALSE)
  reportGeneSyms <- sapply(reportGeneSyms, paste, collapse=", ")
  reportGenes <- sapply(reportGenes, paste, collapse=", ")
  report <- cbind(report,
                  Term=sapply(mget(report$GOBPID, GOstats:::GOenv("TERM")), Term),
                  GeneSyms=reportGeneSyms,
                  label=label)
  rownames(report) <- NULL
  
  return(report)
  
}

geneUniverse <- data.frame(symbol=rownames(se.filt))


library(org.Hs.eg.db)

tabCorr <- AnnotationDbi::select(org.Hs.eg.db, geneUniverse$symbol, "ENTREZID", "SYMBOL")
tabCorr <- tabCorr[!is.na(tabCorr$ENTREZID),]
tabCorr <- tabCorr[!duplicated(tabCorr$SYMBOL),]

geneUniverse$entrezid <- tabCorr[match(geneUniverse$symbol,tabCorr$SYMBOL),]$ENTREZID
geneUniverse <- subset(geneUniverse, !is.na(entrezid))



covariates <- sort(names(table(result$top_factor))[-match(c("Residuals","SV1","SV2","SV3"),names(table(result$top_factor)))])

cov_goenrich <- sapply(covariates, function(x){
  
  tmp <- subset(result, top_factor==x)
  tmp <- tmp[order(tmp$top_var, decreasing=T),]
  
  prop_entrezid <- geneUniverse[match(rownames(tmp), geneUniverse$symbol),]$entrezid
  prop_entrezid <- prop_entrezid[!is.na(prop_entrezid)]
  
  if (length(prop_entrezid)>200){
    
    prop_go <- GOenrichmentAndReport(prop_entrezid, geneUniverse$entrezid,
                                     minSize=5, maxSize=300, minCount=5,
                                     minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                     label=x)
  } else {
    
    prop_go <- GOenrichmentAndReport(prop_entrezid, geneUniverse$entrezid,
                                     minSize=3, maxSize=300, minCount=3,
                                     minOddsRatio=1.5, p.value=0.05, highlightGenes=NULL, highlightStr="*%s*",
                                     label=x)
  }
  
  
  prop_go$ranking <- 1:dim(prop_go)[1]
  #allres_intersect <- rbind(upreg_go[1:10,], downreg_go[1:10,])
  prop_go$short <- prop_go$Term
  mask_long<- nchar(prop_go$Term)>30
  
  
  prop_go$short[mask_long] <- sapply(prop_go[mask_long,]$Term, function(x){
    mask <- unlist(sapply(strsplit(x, " "), function(y) cumsum(nchar(y))>30, simplify=F))
    sapply(strsplit(x, " "), function(y) paste0(paste0(y[!mask], collapse=" "),"\n",paste0(y[mask], collapse=" ")))
  })
  
  return(prop_go)
  
}, simplify=F)

cov_goenrich_high <- cov_goenrich
cov_goenrich_high <- sapply(cov_goenrich_high, function(x) x[1:10,], simplify=F)

cov_goenrich_high <- do.call("rbind", cov_goenrich_high)
rownames(cov_goenrich_high) <- NULL

cov_goenrich_high <- cov_goenrich_high %>%
  group_by(label) %>%
  mutate(short = fct_reorder(short, ranking, .desc = TRUE)) %>% 
  ungroup()

gostats_plot <- ggplot(data=subset(cov_goenrich_high, label=="disease_status"),
                       aes(x=OddsRatio, y=short, size=Count, col=-log10(Pvalue)))+
  theme_bw()+
  geom_point()+
  facet_wrap(~label, nrow=2, scales="free_x")+
  ylab("")+xlab("Odds ratio")+
  scale_color_viridis_c()+
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        axis.text=element_text(size=10))+
  facet_wrap(~label, scales="free_y")+
  theme(axis.title.y=element_text(size=6))+
  ggtitle("GO enrichment per covariate")+
  geom_text(aes(y=short, x=-3, label=ranking), size=4, col="black")


## Main Figure 3G ##

pdf(file="figs/mainFigures/Figure3G.pdf", width=12, height=4)
plot(gostats_plot)
dev.off()




library(devtools)
library(oob)
library(ggbeeswarm)

PCAres <- PCA(assays(se.filt)$logCPM)
varExp <- barplotPercentVar(PCAres, returnGraph = TRUE)



resPC_aov <-
  PCaov(
    PCAres, colData = as.data.frame(colData(se.filt)[, c("disease_status", "treatment","SV1","SV2","SV3")])
  )

resPC_aov$PC <- factor(resPC_aov$PC, levels=unique(resPC_aov$PC))


## Supplementary Figure 2A ##

pdf(file="figs/suppFigures/FigureS2A.pdf"
    , width=10, height=6)
multiplot(
  rows = 2, cols=1,
  barplotPercentVar(PCAres, returnGraph = TRUE
  ),
  ggBorderedFactors(
    ggplot(resPC_aov, aes(
      x = PC, y = sumSqPercent, fill = feature
    )) +
      geom_beeswarm(pch = 21, size = 4, cex = 3) +
      xlab("Principal component") + ylab("% Sum of Squares") +
      scale_fill_manual(values = oobColors(nlevels(resPC_aov$feature))) +
      theme(
        panel.grid.major.y = element_line(colour = "grey75"),
        panel.grid.minor.y = element_line(colour = "grey75"),
        panel.background = element_rect(fill = NA, colour = "black")
      )
  ))
dev.off()


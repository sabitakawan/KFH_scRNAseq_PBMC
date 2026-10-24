

library(ggplot2)
library(RColorBrewer)
library(grid)

inflammasome_monocytes <- c("NLRP3","CASP1","PYCARD","NEK7","IL1B","IL18","GSDMD")
inflammasome_tcells <- c("IFNG")


###############
## MONOCYTES ##
###############

## within treatment
pathToFolder <- "diffExpResults/withinTreatment_diseaseStatusComparison/monocytes/"
mono_within_lps_kfh_vs_ctrl <- read.table(paste0(pathToFolder, "LPS_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))
mono_within_lpsatp_kfh_vs_ctrl <- read.table(paste0(pathToFolder, "LPS-ATP_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))
mono_within_untr_kfh_vs_ctrl <- read.table(paste0(pathToFolder, "Untreated_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))

## within condition (disease)
pathToFolder2 <- "diffExpResults/withinCondition_diseaseStatusComparison/monocytes/"
mono_within_kfh_lps_vs_untr <- read.table(paste0(pathToFolder2, "Patient_LPS_Untreated_DEGs.tsv"),sep="\t", header=T)
mono_within_kfh_lpsatp_vs_untr <- read.table(paste0(pathToFolder2, "Patient_LPS-ATP_Untreated_DEGs.tsv"),sep="\t", header=T)
mono_within_kfh_lpsatp_vs_lps <- read.table(paste0(pathToFolder2, "Patient_LPS-ATP_LPS_DEGs.tsv"),sep="\t", header=T)

## within condition (control)
mono_within_ctrl_lps_vs_untr <- read.table(paste0(pathToFolder2, "Control_LPS_Untreated_DEGs.tsv"),sep="\t", header=T)
mono_within_ctrl_lpsatp_vs_untr <- read.table(paste0(pathToFolder2, "Control_LPS-ATP_Untreated_DEGs.tsv"),sep="\t", header=T)
mono_within_ctrl_lpsatp_vs_lps <- read.table(paste0(pathToFolder2, "Control_LPS-ATP_LPS_DEGs.tsv"),sep="\t", header=T)


vars <- ls()[grepl("mono_", ls())]

stopifnot(sapply(vars, function(x){
  df <- get(x)
  all(inflammasome_monocytes %in% df$gene)
}))



alldfs <- sapply(vars, function(x){
  
  df <- get(x)
  df <- df[df$gene %in% inflammasome_monocytes,]
  df <- df[order(df$gene),]
  
  
  to_return <- data.frame(gene=df$gene,
                          avg_log2FC=df$avg_log2FC,
                          p_val_adj=df$p_val_adj,
                          test=x,
                          ctype="Monocytes")
  
  return(to_return)
  
}, simplify=F)


alldfs <- do.call("rbind", alldfs)
rownames(alldfs) <- NULL


## pval_signif

pvalConverter <- function(vectorPval){
  symbolVec <- rep("ns", length(vectorPval))
  
  if (any(vectorPval<0.05)){
    symbolVec[vectorPval<0.05] <- "*"
  }
  
  if (any(vectorPval<0.01)){
    symbolVec[vectorPval<0.01] <- "**"
  }
  
  if (any(vectorPval<0.001)){
    symbolVec[vectorPval<0.001] <- "***"
  }                       
  return(symbolVec)
}

alldfs$signifLevel <- pvalConverter(alldfs$p_val_adj)



myPalette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
maxOdds <- ceiling(max(alldfs$avg_log2FC, na.rm=T))
minOdds <- -maxOdds


size_values=c("ns"=NA,
              "*"=1,
              "**"=3,
              "***"=4)
my_breaks <- seq(minOdds, maxOdds, 2)



alldfs$test <- factor(
  alldfs$test,
  levels = c("mono_within_untr_kfh_vs_ctrl","mono_within_lps_kfh_vs_ctrl","mono_within_lpsatp_kfh_vs_ctrl",
             "mono_within_kfh_lps_vs_untr","mono_within_kfh_lpsatp_vs_untr","mono_within_kfh_lpsatp_vs_lps",
             "mono_within_ctrl_lps_vs_untr","mono_within_ctrl_lpsatp_vs_untr","mono_within_ctrl_lpsatp_vs_lps"))

alldfs$gene <- factor(
  alldfs$gene,
  levels = c("CASP1","GSDMD","IL18","IL1B","NEK7","NLRP3","PYCARD"))



monoPlot <- ggplot(data=alldfs, aes(x=test, y=gene, fill=avg_log2FC))+
  geom_tile()+
  theme_classic()+
  theme(axis.text.x=element_text(vjust=0.5, hjust=1, angle=90),
        plot.title=element_text(hjust=0.5, face="bold", size=15),
        axis.title=element_text(size=14, face="bold"),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  scale_fill_gradient2(
    low = brewer.pal(9, "RdBu")[9],
    mid = "white",
    high = brewer.pal(9, "RdBu")[1],
    midpoint = 0,
    limits = c(my_breaks[1],  my_breaks[length(my_breaks)]),
    breaks = my_breaks,
    oob = scales::squish
  )+
  scale_size_manual(
    name = "Adj.Pval",
    values = size_values,
    breaks = rev(c("***", "**", "*")),
    labels = rev(c("<0.001", "<0.01", "<0.05"))
  )+
  #ggtitle("Number of DE genes / Gene set enrichment")+
  geom_point(data=alldfs, aes(x=test, y=gene, size=signifLevel), inherit.aes = F)+
  theme(legend.title=element_text(face="bold"))+
  guides(
    fill = guide_colourbar(
      barheight = unit(5.5, "cm"),
      barwidth = unit(0.75, "cm")
    ),
    size = guide_legend(order = 2)
  )+
  scale_x_discrete(
    labels = c("KFH vs Ctrl (Untr)", "KFH vs Ctrl (LPS)", "KFH vs Ctrl (LPS-ATP)",
               "LPS vs Untr (KFH)","LPS-ATP vs Untr (KFH)","LPS-ATP vs LPS (KFH)",
               "LPS vs Untr (Ctrl)","LPS-ATP vs Untr (Ctrl)","LPS-ATP vs LPS (Ctrl)")
  )+
  ggtitle("Monocytes: NLRP3 inflammasome")+
  xlab("")+ylab("")


## Main Figure 4D ##
pdf(file="figs/mainFigures/Figure4D.pdf")
plot(monoPlot)
dev.off()




#############################
## T-CELLS (subsample 25%) ##
#############################


## within treatment
pathToFolder3 <- "diffExpResults/withinTreatment_diseaseStatusComparison/tcells/"
tcells_within_lps_kfh_vs_ctrl <- read.table(paste0(pathToFolder3, "LPS_tcells_disease-specific_DEG_19072026.xlsx"))
tcells_within_lpsatp_kfh_vs_ctrl <- read.table(paste0(pathToFolder3, "LPS-ATP_tcells_disease-specific_DEG_19072026.xlsx"))
tcells_within_untr_kfh_vs_ctrl <- read.table(paste0(pathToFolder3, "Untreated_tcells_disease-specific_DEG_19072026.xlsx"))

## within condition (disease)
pathToFolder4 <- "diffExpResults/withinCondition_diseaseStatusComparison/tcells/"
tcells_within_kfh_lps_vs_untr <- read.table(paste0(pathToFolder4, "Patient_LPS_Untreated_DEGs.tsv"),sep="\t", header=T)
tcells_within_kfh_lpsatp_vs_untr <- read.table(paste0(pathToFolder4, "Patient_LPS-ATP_Untreated_DEGs.tsv"),sep="\t", header=T)
tcells_within_kfh_lpsatp_vs_lps <- read.table(paste0(pathToFolder4, "Patient_LPS-ATP_LPS_DEGs.tsv"),sep="\t", header=T)

## within condition (control)
tcells_within_ctrl_lps_vs_untr <- read.table(paste0(pathToFolder4, "Control_LPS_Untreated_DEGs.tsv"),sep="\t", header=T)
tcells_within_ctrl_lpsatp_vs_untr <- read.table(paste0(pathToFolder4, "Control_LPS-ATP_Untreated_DEGs.tsv"),sep="\t", header=T)
tcells_within_ctrl_lpsatp_vs_lps <- read.table(paste0(pathToFolder4, "Control_LPS-ATP_LPS_DEGs.tsv"),sep="\t", header=T)


vars <- ls()[grepl("tcells_", ls())]

stopifnot(sapply(vars, function(x){
  df <- get(x)
  all(inflammasome_tcells %in% df$gene)
}))



alldfs <- sapply(vars, function(x){
  
  df <- get(x)
  df <- df[df$gene %in% inflammasome_tcells,]
  df <- df[order(df$gene),]
  
  
  to_return <- data.frame(gene=df$gene,
                          avg_log2FC=df$avg_log2FC,
                          p_val_adj=df$p_val_adj,
                          test=x,
                          ctype="T-cells")
  
  return(to_return)
  
}, simplify=F)


alldfs <- do.call("rbind", alldfs)
rownames(alldfs) <- NULL


## pval_signif

pvalConverter <- function(vectorPval){
  symbolVec <- rep("ns", length(vectorPval))
  
  if (any(vectorPval<0.05)){
    symbolVec[vectorPval<0.05] <- "*"
  }
  
  if (any(vectorPval<0.01)){
    symbolVec[vectorPval<0.01] <- "**"
  }
  
  if (any(vectorPval<0.001)){
    symbolVec[vectorPval<0.001] <- "***"
  }                       
  return(symbolVec)
}

alldfs$signifLevel <- pvalConverter(alldfs$p_val_adj)



myPalette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
maxOdds <- ceiling(max(alldfs$avg_log2FC, na.rm=T))
minOdds <- -maxOdds


size_values=c("ns"=NA,
              "*"=0.5,
              "**"=1,
              "***"=4)
my_breaks <- seq(minOdds, maxOdds, 2)



alldfs$test <- factor(
  alldfs$test,
  levels = c("tcells_within_untr_kfh_vs_ctrl","tcells_within_lps_kfh_vs_ctrl","tcells_within_lpsatp_kfh_vs_ctrl",
             "tcells_within_kfh_lps_vs_untr","tcells_within_kfh_lpsatp_vs_untr","tcells_within_kfh_lpsatp_vs_lps",
             "tcells_within_ctrl_lps_vs_untr","tcells_within_ctrl_lpsatp_vs_untr","tcells_within_ctrl_lpsatp_vs_lps"))


tcellsPlot <- ggplot(data=alldfs, aes(x=test, y=gene, fill=avg_log2FC))+
  geom_tile()+
  theme_classic()+
  theme(axis.text.x=element_text(vjust=0.5, hjust=1, angle=90),
        plot.title=element_text(hjust=0.5, face="bold", size=15),
        axis.title=element_text(size=14, face="bold"),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  scale_fill_gradient2(
    low = brewer.pal(9, "RdBu")[9],
    mid = "white",
    high = brewer.pal(9, "RdBu")[1],
    midpoint = 0,
    limits = c(my_breaks[1],  my_breaks[length(my_breaks)]),
    breaks = c(-1,-0.5,0.,0.5,1),
    oob = scales::squish
  )+
  scale_size_manual(
    name = "Adj.Pval",
    values = size_values,
    breaks = rev(c("***", "**", "*")),
    labels = rev(c("<0.001", "<0.01", "<0.05"))
  )+
  #ggtitle("Number of DE genes / Gene set enrichment")+
  geom_point(data=alldfs, aes(x=test, y=gene, size=signifLevel), inherit.aes = F)+
  theme(legend.title=element_text(face="bold"))+
  guides(
    fill = guide_colourbar(
      barheight = unit(5.5, "cm"),
      barwidth = unit(0.75, "cm")
    ),
    size = guide_legend(order = 2)
  )+
  scale_x_discrete(
    labels = c("KFH vs Ctrl (Untr)", "KFH vs Ctrl (LPS)", "KFH vs Ctrl (LPS-ATP)",
               "LPS vs Untr (KFH)","LPS-ATP vs Untr (KFH)","LPS-ATP vs LPS (KFH)",
               "LPS vs Untr (Ctrl)","LPS-ATP vs Untr (Ctrl)","LPS-ATP vs LPS (Ctrl)")
  )+
  ggtitle("T-cells: IFNG")+
  xlab("")+ylab("")



## Main Figure 7E ##
pdf(file="figs/mainFigures/Figure7E.pdf", width=8, height=3.5)
plot(tcellsPlot)
dev.off()




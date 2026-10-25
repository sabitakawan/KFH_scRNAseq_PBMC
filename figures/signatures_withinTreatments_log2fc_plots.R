
library(ggplot2)
library(RColorBrewer)
library(grid)
library(ggrepel)


#################
### Monocytes ###
#################

## within treatment

pathToFolder <- "diffExpResults/withinTreatment_diseaseStatusComparison/monocytes/"
mono_within_lps_kfh_vs_ctrl <- read.table(paste0(pathToFolder, "LPS_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))
mono_within_lpsatp_kfh_vs_ctrl <- read.table(paste0(pathToFolder, "LPS-ATP_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))
mono_within_untr_kfh_vs_ctrl <- read.table(paste0(pathToFolder, "Untreated_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))



## create function

plotterFunction <- function(obj1=mono_within_lps_kfh_vs_ctrl, obj2=mono_within_untr_kfh_vs_ctrl,
                            def.cond1= "KFH vs ctrl (LPS) [in log2FC]", def.cond2= "KFH vs ctrl (Untreated) [in log2FC]",
                            title_text="Effect of priming on disease",
                            figureName="log2fc_cor_lps_vs_untr_on_kfh_vs_ctrl.pdf"){
  
  genes_group1 <- obj1[!(obj1$pct.1==0 & obj1$pct.2==0),]$gene
  genes_group2 <- obj2[!(obj2$pct.1==0 & obj2$pct.2==0),]$gene
  union_genes <- union(genes_group1, genes_group2)
  
  obj1 <- obj1[obj1$gene %in% union_genes,]
  obj2 <- obj2[obj2$gene %in% union_genes,]
  
  obj2 <- obj2[match(obj1$gene, obj2$gene),]
  stopifnot(all(obj2$gene==obj1$gene))
  
  
  toPlot <- data.frame(gene=obj1$gene,
                       cond.1=def.cond1,
                       cond.2=def.cond2,
                       cond.1.log2fc=obj1$avg_log2FC,
                       cond.2.log2fc=obj2$avg_log2FC,
                       sign.1=ifelse(obj1$p_val_adj>0.05 | is.na(obj1$p_val_adj), FALSE, TRUE),
                       sign.2=ifelse(obj2$p_val_adj>0.05 | is.na(obj2$p_val_adj), FALSE, TRUE),
                       dir.1=ifelse(abs(obj1$avg_log2FC)>log2(1.5), "de","not_de"),
                       dir.2=ifelse(abs(obj2$avg_log2FC)>log2(1.5), "de","not_de"))
  
  
  mask1 <- toPlot$dir.1=="de"
  mask2 <- toPlot$dir.2=="de"
  
  toPlot[mask1,]$dir.1 <- ifelse(toPlot[mask1,]$cond.1.log2fc>0, "upreg","downreg")
  toPlot[mask2,]$dir.2 <- ifelse(toPlot[mask2,]$cond.2.log2fc>0, "upreg","downreg")
  
  toPlot$sign.1.strong <- ifelse(toPlot$sign.1==TRUE & toPlot$dir.1!="not_de", TRUE, FALSE)
  toPlot$sign.2.strong <- ifelse(toPlot$sign.2==TRUE & toPlot$dir.2!="not_de", TRUE, FALSE)
  
  ### colors
  
  toPlot$group <- "NS !(FDR<0.05 & abs(FC)>1.5)"
  
  # both up
  toPlot$group[
    toPlot$sign.1.strong &
      toPlot$sign.2.strong &
      toPlot$dir.1=="upreg" &
      toPlot$dir.2=="upreg"
  ] <- "Both up"
  
  # one up
  toPlot$group[
    xor(toPlot$sign.1.strong & toPlot$dir.1=="upreg",
        toPlot$sign.2.strong & toPlot$dir.2=="upreg")
  ] <- "One up"
  
  # both down
  toPlot$group[
    toPlot$sign.1.strong &
      toPlot$sign.2.strong &
      toPlot$dir.1=="downreg" &
      toPlot$dir.2=="downreg"
  ] <- "Both down"
  
  # one down
  toPlot$group[
    xor(toPlot$sign.1.strong & toPlot$dir.1=="downreg",
        toPlot$sign.2.strong & toPlot$dir.2=="downreg")
  ] <- "One down"
  
  # opposite directions
  toPlot$group[
    toPlot$sign.1.strong &
      toPlot$sign.2.strong &
      toPlot$dir.1 != toPlot$dir.2
  ] <- "Opposite"
  
  
  cols <- c(
    "Both up"   = "#B22222",
    "One up"    = "#FB6A4A",
    "Both down" = "#08519C",
    "One down"  = "#6BAED6",
    "Opposite"  = "#7B68EE",
    "NS !(FDR<0.05 & abs(FC)>1.5)"        = "grey85"
  )
  
  toPlot$opposite <- with(toPlot,
                          sign.1.strong &
                            sign.2.strong &
                            dir.1 %in% c("upreg","downreg") &
                            dir.2 %in% c("upreg","downreg") &
                            dir.1 != dir.2)
  
  toPlot <- ggplot(subset(toPlot, !(sign.1==FALSE & sign.2==FALSE)),
         aes(x=cond.1.log2fc, y=cond.2.log2fc, color=group, shape = opposite))+
    geom_point(alpha=0.8, size=1.5)+
    theme_bw()+
    xlab(def.cond1)+
    ylab(def.cond2)+
    geom_hline(yintercept=log2(1.5), linetype="dashed")+
    geom_vline(xintercept=log2(1.5), linetype="dashed")+
    geom_hline(yintercept=-log2(1.5), linetype="dashed")+
    geom_vline(xintercept=-log2(1.5), linetype="dashed")+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    ggtitle(title_text)+
    theme(plot.title=element_text(size=15, hjust=0.5),
          axis.text=element_text(size=13))+
    geom_text_repel(data=subset(toPlot, sign.1.strong==TRUE | sign.2.strong==TRUE), aes(label=gene),
                    show.legend = FALSE)+
    scale_colour_manual(name="Direction DE",values = cols)+
    scale_shape_manual(name="Opposite DE", values = c(`FALSE` = 16, `TRUE` = 18))
  
  
  pdf(file=paste0("figs/mainFigures/",figureName), width=10, height=8)
  plot(toPlot)
  dev.off()
  
}


## KFH vs ctrl (LPS vs untr)
## Main Figure 5A (right) ##

plotterFunction(obj1=mono_within_lps_kfh_vs_ctrl, obj2=mono_within_untr_kfh_vs_ctrl,
                def.cond1= "KFH vs ctrl (LPS) [in log2FC]", def.cond2= "KFH vs ctrl (Untreated) [in log2FC]",
                title_text="Effect of priming on disease",
                figureName="Figure5A.pdf")
  
## KFH vs ctrl (LPS-ATP vs untr)
## Main Figure 5B (right) ##

plotterFunction(obj1=mono_within_lpsatp_kfh_vs_ctrl, obj2=mono_within_untr_kfh_vs_ctrl,
                def.cond1= "KFH vs ctrl (LPS-ATP) [in log2FC]", def.cond2= "KFH vs ctrl (Untreated) [in log2FC]",
                title_text="Effect of priming+activation on disease",
                figureName="Figure5B.pdf")

## KFH vs ctrl (LPS-ATP vs LPS)
## Main Figure 5C (right) ##

plotterFunction(obj1=mono_within_lpsatp_kfh_vs_ctrl, obj2=mono_within_lps_kfh_vs_ctrl,
                def.cond1= "KFH vs ctrl (LPS-ATP) [in log2FC]", def.cond2= "KFH vs ctrl (LPS) [in log2FC]",
                title_text="Effect of activation on disease",
                figureName="Figure5C.pdf")


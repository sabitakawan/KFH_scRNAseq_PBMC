
library(readxl)
library(RRHO)

## Code for RRHO plots is adapted from: Qi et al. 2026, Neuron (DOI: 10.1016/j.neuron.2026.09.022)

pathToDirectory <- "diffExpResults/withinTreatment_diseaseStatusComparison/monocytes/"

#Untreated_Monocytes_DCs_Macrophages_disease-specific_DEG
untreated_mono <- read.table(paste0(pathToDirectory,"Untreated_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))
untreated_mono <- untreated_mono[!(untreated_mono$pct.1==0 & untreated_mono$pct.2==0),] ## keep only genes with non-zero expression in monocytes

#LPS_Monocytes_DCs_Macrophages_disease-specific_DE
lps_treated_mono <- read.table(paste0(pathToDirectory,"LPS_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))
lps_treated_mono <- lps_treated_mono[!(lps_treated_mono$pct.1==0 & lps_treated_mono$pct.2==0),]

#LPS_Monocytes_DCs_Macrophages_disease-specific_DE
lpsatp_treated_mono <- read.table(paste0(pathToDirectory,"LPS-ATP_Monocytes_DCs_Macrophages_disease-specific_DEGs_29062026.xlsx"))
lpsatp_treated_mono <- lpsatp_treated_mono[!(lpsatp_treated_mono$pct.1==0 & lpsatp_treated_mono$pct.2==0),]

length(intersect(rownames(untreated_mono), rownames(lps_treated_mono)))
# [1] 16920

length(intersect(rownames(untreated_mono), rownames(lpsatp_treated_mono)))
# [1] 16181

length(intersect(rownames(lps_treated_mono), rownames(lpsatp_treated_mono)))
# [1] 16354


RRHO_plot_custom<-function (list1, list2, stepsize = defaultStepSize(list1, list2), 
                            labels, alternative, plots = FALSE, outputdir = NULL, BY = FALSE, 
                            log10.ind = FALSE) 
{
  if (length(list1[, 1]) != length(unique(list1[, 1]))) 
    stop("Non-unique gene identifier found in list1")
  if (length(list2[, 1]) != length(unique(list2[, 1]))) 
    stop("Non-unique gene identifier found in list2")
  if (plots && (missing(outputdir) || missing(labels))) 
    stop("When plots=TRUE, outputdir and labels are required.")
  if (!(alternative == "two.sided" || alternative == "enrichment")) 
    stop("Wrong alternative specified.")
  result <- list(hypermat = NA, hypermat.counts = NA, hypermat.signs = NA, 
                 hypermat.by = NA, n.items = nrow(list1), stepsize = stepsize, 
                 log10.ind = log10.ind, call = match.call())
  list1 <- list1[order(list1[, 2], decreasing = TRUE), ]
  list2 <- list2[order(list2[, 2], decreasing = TRUE), ]
  nlist1 <- length(list1[, 1])
  nlist2 <- length(list2[, 1])
  N <- max(nlist1, nlist2)
  .hypermat <- RRHO:::numericListOverlap(list1[, 1], list2[, 1], 
                                         stepsize, alternative)
  hypermat <- .hypermat$log.pval
  if (log10.ind) 
    hypermat <- hypermat * log10(exp(1))
  if (BY) {
    hypermatvec <- matrix(hypermat, nrow = nrow(hypermat) * 
                            ncol(hypermat), ncol = 1)
    hypermat.byvec <- p.adjust(exp(-hypermatvec), method = "BY")
    hypermat.by <- matrix(-log(hypermat.byvec), nrow = nrow(hypermat), 
                          ncol = ncol(hypermat))
    if (log10.ind) 
      hypermat.by <- hypermat.by * log10(exp(1))
    result$hypermat.by <- hypermat.by
  }
  if (plots) {
    hypermat.signed <- hypermat * .hypermat$signs
    
    # --- clip infinite cells so image() and color.bar() agree on range ---
    finite.ind <- is.finite(hypermat.signed)
    posmax <- max(hypermat.signed[finite.ind], na.rm = TRUE)
    negmin <- min(hypermat.signed[finite.ind], na.rm = TRUE)
    hypermat.signed[hypermat.signed ==  Inf] <- posmax
    hypermat.signed[hypermat.signed == -Inf] <- negmin
    # ----------------------------------------------------------------
    
    color.bar <- function(lut, min, max = -min, nticks = 11, 
                          ticks = seq(min, max, len = nticks), title = "") {
      scale <- (length(lut) - 1)/(max - min)
      plot(c(0, 10), c(min, max), type = "n", bty = "n", 
           xaxt = "n", xlab = "", yaxt = "n", ylab = "")
      mtext(title, 2, 2.3, cex = 0.8)
      axis(2, round(ticks, 0), las = 1, cex.lab = 0.8)
      for (i in 1:(length(lut) - 1)) {
        y <- (i - 1)/scale + min
        rect(0, y, 10, y + 1/scale, col = lut[i], 
             border = NA)
      }
    }
    
    .filename <- paste("RRHOMap", labels[1], "_VS_", 
                       labels[2], ".pdf", sep = "")
    .filename <- paste("RRHOMap", labels[1], "_VS_", 
                       labels[2], ".pdf", sep = "")
    pdf(file = paste(outputdir, .filename, sep = "/"),width = 12, height = 12)
    jet.colors <- colorRampPalette(c("#00007F", "blue", 
                                     "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
                                     "red", "#7F0000"))
    layout(matrix(c(rep(1, 5), 2), 1, 6, byrow = TRUE))
    image(hypermat.signed, xlab = "", ylab = "", col = jet.colors(100), 
          axes = FALSE, main = "",
          zlim = c(-35, 750))   # <-- add this so image() bins match color.bar
    mtext(labels[2], 2, 0.5)
    mtext(labels[1], 1, 0.5)
    finite.ind <- is.finite(hypermat.signed)   # now redundant but harmless — everything is finite after clipping
    color.bar(jet.colors(100), min = -35, max = 750, nticks = 6, title = "signed -log(P-value)")
    dev.off()
    }
    return(NA)}


## Comparison1 (Untreated vs LPS)

## Main Figure 5A (left) ##

name_i <- "monocytes_LPS"
name_j <- "monocytes_Untreated"
untreated_mono_c1 <- untreated_mono[intersect(rownames(untreated_mono),rownames(lps_treated_mono)),c("gene","avg_log2FC")]
lps_treated_mono_c1 <- lps_treated_mono[intersect(rownames(lps_treated_mono),rownames(untreated_mono)),c("gene","avg_log2FC")]
RRHO_plot_custom(lps_treated_mono_c1,untreated_mono_c1,stepsize = 100, labels = c(name_i,name_j),plots = T,
                 alternative = "two.sided",
                 outputdir = "figs/mainFigures/Figure5A.pdf") #the figure is the heatmap in this folder





## Comparison2 (Untreated vs LPS-ATP)

## Main Figure 5B (left) ##


name_i <- "monocytes_LPS_ATP"
name_j <- "monocytes_Untreated"
untreated_mono_c2 <- untreated_mono[intersect(rownames(untreated_mono),rownames(lpsatp_treated_mono)),c("gene","avg_log2FC")]
lpsatp_treated_mono_c2 <- lpsatp_treated_mono[intersect(rownames(lpsatp_treated_mono),rownames(untreated_mono)),c("gene","avg_log2FC")]
RRHO_plot_custom(lpsatp_treated_mono_c2,untreated_mono_c2,stepsize = 100, labels = c(name_i,name_j),plots = T,
                 alternative = "two.sided",
                 outputdir = "figs/mainFigures/Figure5B.pdf") #the figure is the heatmap in this folder


## Comparison3 (LPS vs LPS-ATP)

## Main Figure 5C (left) ##

name_i <- "monocytes_LPS_ATP"
name_j <- "monocytes_LPS"
lps_treated_mono_c3 <- lps_treated_mono[intersect(rownames(lps_treated_mono),rownames(lpsatp_treated_mono)),c("gene","avg_log2FC")]
lpsatp_treated_mono_c3 <- lpsatp_treated_mono[intersect(rownames(lpsatp_treated_mono),rownames(lps_treated_mono)),c("gene","avg_log2FC")]


RRHO_plot_custom(lpsatp_treated_mono_c3, lps_treated_mono_c3, stepsize = 100, labels = c(name_i,name_j),plots = T,
                 alternative = "two.sided",
                 outputdir = "figs/mainFigures/Figure5C.pdf") #the figure is the heatmap in this folder



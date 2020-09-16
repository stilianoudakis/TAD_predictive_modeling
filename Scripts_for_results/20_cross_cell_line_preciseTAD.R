#20_cross_cell_line_preciseTAD

library(S4Vectors)
library(GenomicRanges)
library(ggplot2)
library(caret)
library(randomForest)
library(glmnet)
library(gbm)
library(DMwR)
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(cluster)
library(pROC)
library(SuperLearner)
library(data.table)
library(PRROC)
library(bigmemory)
library(ranger)
library(dplyr)
library(parallel)
library(doSNOW)
library(foreach)
library(pbapply)
library(ggpubr)
library(scales)
library(GGally)
library(VennDiagram)
library(network)
library(sna)
library(gtable)
library(grid)
library(magrittr) 
library(DT)
library(tidyverse)
library(reshape)
library(ROSE)
library(lattice)
library(corrplot)
library(cluster)
library(RColorBrewer)
#library(GenometriCorr)
library(ggsignif)
library(Vennerable)
library(VennDiagram)
library(ChIPpeakAnno)
library(EnrichedHeatmap)

source("Z:/TAD_data_analysis/functions_for_R_package/createTADdata.R")
source("Z:/TAD_data_analysis/functions_for_R_package/preciseTAD.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADrandomForest.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADrfe2.R")
source("Z:/TAD_data_analysis/functions_for_R_package/distance_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/percent_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/count_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/binary_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/signal_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/annots_to_granges_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/extract_boundaries_func.R")

venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Venn(SetNames=SetNames, Weight=Weight)
}

genomicCorr.jaccard = function(query, reference, restrict = NULL) {
  if(is.null(restrict)) {
    res = sum(width(GenomicRanges::intersect(query, reference))) / sum(as.numeric(width(GenomicRanges::union(query, reference))))
  } else {
    gr1 = GenomicRanges::intersect(query, reference)
    gr1 = GenomicRanges::intersect(gr1, restrict)
    
    gr2 = GenomicRanges::union(query, reference)
    gr2 = GenomicRanges::intersect(gr2, restrict)
    res = sum(width(gr1)) / sum(width(gr2))
  }
  return(res)
}

#arrowhead

##gm12878 on gm12878 vs k562 on gm12878

g_on_g <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  g_on_g[[i]] <- called_and_pred$PTBP
  
}
g_on_g <- do.call("c", g_on_g)
g_on_g <- flank(g_on_g, 5000, both=TRUE)

k_on_g <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_k562_on_gm12878.rds"))
  
  k_on_g[[i]] <- called_and_pred$PTBP
  
}
k_on_g <- do.call("c", k_on_g)
k_on_g <- flank(k_on_g, 5000, both=TRUE)


res <- makeVennDiagram(Peaks=list(g_on_g,
                                  k_on_g),
                       NameOfPeaks=c("GM12878 on GM12878", "K562 on GM12878"))

venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "purple"
gp[["Face"]][["01"]]$fill <-  "lightblue"
gp[["Face"]][["10"]]$fill <-  "red"
#gp$Set$Set1$col="green"
#gp$Set$Set2$col="green"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .01
gp[["SetText"]][["Set2"]]$cex <- .01

gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(g_on_g,k_on_g)

##k562 on k562 vs gm12878 on k562

k_on_k <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  k_on_k[[i]] <- called_and_pred$PTBP
  
}
k_on_k <- do.call("c", k_on_k)
k_on_k <- flank(k_on_k, 5000, both=TRUE)

g_on_k <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_gm12878_on_k562.rds"))
  
  g_on_k[[i]] <- called_and_pred$PTBP
  
}
g_on_k <- do.call("c", g_on_k)
g_on_k <- flank(g_on_k, 5000, both=TRUE)

res <- makeVennDiagram(Peaks=list(k_on_k,
                                  g_on_k),
                       NameOfPeaks=c("K562 on K562", "GM12878 on K562"))

venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "purple"
gp[["Face"]][["01"]]$fill <-  "lightblue"
gp[["Face"]][["10"]]$fill <-  "red"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .01
gp[["SetText"]][["Set2"]]$cex <- .01

gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(k_on_k,g_on_k)


#peakachu

##gm12878 on gm12878 vs k562 on gm12878

g_on_g <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  g_on_g[[i]] <- called_and_pred$PTBP
  
}
g_on_g <- do.call("c", g_on_g)
g_on_g <- flank(g_on_g, 10000, both=TRUE)

k_on_g <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_k562_on_gm12878_peakachu.rds"))
  
  k_on_g[[i]] <- called_and_pred$PTBP
  
}
k_on_g <- do.call("c", k_on_g)
k_on_g <- flank(k_on_g, 10000, both=TRUE)

res <- makeVennDiagram(Peaks=list(g_on_g,
                                  k_on_g),
                       NameOfPeaks=c("GM12878 on GM12878", "K562 on GM12878"))

venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "purple"
gp[["Face"]][["01"]]$fill <-  "lightblue"
gp[["Face"]][["10"]]$fill <-  "red"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .01
gp[["SetText"]][["Set2"]]$cex <- .01

gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(g_on_g,k_on_g)

##k562 on k562 vs gm12878 on k562

k_on_k <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  k_on_k[[i]] <- called_and_pred$PTBP
  
}
k_on_k <- do.call("c", k_on_k)
k_on_k <- flank(k_on_k, 10000, both=TRUE)

g_on_k <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_gm12878_on_k562_peakachu.rds"))
  
  g_on_k[[i]] <- called_and_pred$PTBP
  
}
g_on_k <- do.call("c", g_on_k)
g_on_k <- flank(g_on_k, 10000, both=TRUE)


res <- makeVennDiagram(Peaks=list(k_on_k,
                                  g_on_k),
                       NameOfPeaks=c("K562 on K562", "GM12878 on K562"))

venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "purple"
gp[["Face"]][["01"]]$fill <-  "lightblue"
gp[["Face"]][["10"]]$fill <-  "red"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .01
gp[["SetText"]][["Set2"]]$cex <- .01

gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(k_on_k,g_on_k)

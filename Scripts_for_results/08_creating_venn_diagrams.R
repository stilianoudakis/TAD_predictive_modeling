#08_creating_venn_diagrams

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
library(GenometriCorr)
library(ggsignif)
library(Vennerable)
library(VennDiagram)
library(ChIPpeakAnno)
library(EnrichedHeatmap)
library(ggsci)

scales::show_col(pal_lancet("lanonc")(8))
mycols = pal_lancet("lanonc")(8)

w <- wes_palette(n=3, "GrandBudapest2")
w
str(w)

source("Z:/TAD_data_analysis/functions_for_R_package/preciseTAD.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADRF.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADrfe.R")
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


#nested within cell line

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  true_bound_list[[i]] <- called_and_pred$CTBP
  pred_bound_list[[i]] <- called_and_pred$PTBP
  
}
all_true_bounds_gm12878_a <- do.call("c", true_bound_list)
all_pred_bounds_gm12878_a <- do.call("c", pred_bound_list)
all_true_bounds_gm12878_a <- flank(all_true_bounds_gm12878_a, 5000, both=TRUE)
all_pred_bounds_gm12878_a <- flank(all_pred_bounds_gm12878_a, 5000, both=TRUE)


pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  true_bound_list[[i]] <- called_and_pred$CTBP
  pred_bound_list[[i]] <- called_and_pred$PTBP
  
}
all_true_bounds_gm12878_p <- do.call("c", true_bound_list)
all_pred_bounds_gm12878_p <- do.call("c", pred_bound_list)
all_true_bounds_gm12878_p <- flank(all_true_bounds_gm12878_p, 10000, both=TRUE)
all_pred_bounds_gm12878_p <- flank(all_pred_bounds_gm12878_p, 10000, both=TRUE)


pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  true_bound_list[[i]] <- called_and_pred$CTBP
  pred_bound_list[[i]] <- called_and_pred$PTBP
}
all_true_bounds_k562_a <- do.call("c", true_bound_list)
all_pred_bounds_k562_a <- do.call("c", pred_bound_list)
all_true_bounds_k562_a <- flank(all_true_bounds_k562_a, 5000, both=TRUE)
all_pred_bounds_k562_a <- flank(all_pred_bounds_k562_a, 5000, both=TRUE)


pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  true_bound_list[[i]] <- called_and_pred$CTBP
  pred_bound_list[[i]] <- called_and_pred$PTBP
}
all_true_bounds_k562_p <- do.call("c", true_bound_list)
all_pred_bounds_k562_p <- do.call("c", pred_bound_list)
all_true_bounds_k562_p <- flank(all_true_bounds_k562_p, 10000, both=TRUE)
all_pred_bounds_k562_p <- flank(all_pred_bounds_k562_p, 10000, both=TRUE)

## arrowhead vs peakachu

res <- makeVennDiagram(Peaks=list(all_pred_bounds_gm12878_a,
                                  all_pred_bounds_gm12878_p),
                       NameOfPeaks=c("Arrowhead", "PEAKACHU"))
venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "#D8A499"
gp[["Face"]][["01"]]$fill <-  "#E6A0C4"
gp[["Face"]][["10"]]$fill <-  "#C6CDF7"
gp$Set$Set1$col <- "#C6CDF7"
gp$Set$Set2$col <- "#E6A0C4"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .1
gp[["SetText"]][["Set2"]]$cex <- .1
gp$SetText$Set1$col <- "white"
gp$SetText$Set2$col <- "white"
gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(all_pred_bounds_gm12878_a,
                    all_pred_bounds_gm12878_p)

res <- makeVennDiagram(Peaks=list(all_pred_bounds_k562_a,
                                  all_pred_bounds_k562_p),
                       NameOfPeaks=c("Arrowhead", "PEAKACHU"))
venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "#D8A499"
gp[["Face"]][["01"]]$fill <-  "#E6A0C4"
gp[["Face"]][["10"]]$fill <-  "#C6CDF7"
gp$Set$Set1$col <- "#C6CDF7"
gp$Set$Set2$col <- "#E6A0C4"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .1
gp[["SetText"]][["Set2"]]$cex <- .1
gp$SetText$Set1$col <- "white"
gp$SetText$Set2$col <- "white"
gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(all_pred_bounds_k562_a,
                    all_pred_bounds_k562_p)

#across cell line

res <- makeVennDiagram(Peaks=list(all_true_bounds_gm12878_a,
                                  all_true_bounds_k562_a),
                       NameOfPeaks=c("ArrowheadGM", "ArrowheadK"))
venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "#D8A499"
gp[["Face"]][["01"]]$fill <-  "#E6A0C4"
gp[["Face"]][["10"]]$fill <-  "#C6CDF7"
gp$Set$Set1$col <- "#C6CDF7"
gp$Set$Set2$col <- "#E6A0C4"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .1
gp[["SetText"]][["Set2"]]$cex <- .1
gp$SetText$Set1$col <- "white"
gp$SetText$Set2$col <- "white"
gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(all_true_bounds_gm12878_a,
                    all_true_bounds_k562_a)



res <- makeVennDiagram(Peaks=list(all_pred_bounds_gm12878_a,
                                  all_pred_bounds_k562_a),
                       NameOfPeaks=c("preciseTADGM", "preciseTADK"))
venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "#D8A499"
gp[["Face"]][["01"]]$fill <-  "#E6A0C4"
gp[["Face"]][["10"]]$fill <-  "#C6CDF7"
gp$Set$Set1$col <- "#C6CDF7"
gp$Set$Set2$col <- "#E6A0C4"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .1
gp[["SetText"]][["Set2"]]$cex <- .1
gp$SetText$Set1$col <- "white"
gp$SetText$Set2$col <- "white"
gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(all_pred_bounds_gm12878_a,
                    all_pred_bounds_k562_a)


res <- makeVennDiagram(Peaks=list(all_pred_bounds_gm12878_p,
                                  all_pred_bounds_k562_p),
                       NameOfPeaks=c("preciseTADGM", "preciseTADK"))
venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "#D8A499"
gp[["Face"]][["01"]]$fill <-  "#E6A0C4"
gp[["Face"]][["10"]]$fill <-  "#C6CDF7"
gp$Set$Set1$col <- "#C6CDF7"
gp$Set$Set2$col <- "#E6A0C4"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .1
gp[["SetText"]][["Set2"]]$cex <- .1
gp$SetText$Set1$col <- "white"
gp$SetText$Set2$col <- "white"
gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(all_pred_bounds_gm12878_p,
                    all_pred_bounds_k562_p)


res <- makeVennDiagram(Peaks=list(all_true_bounds_gm12878_p,
                                  all_true_bounds_k562_p),
                       NameOfPeaks=c("PeakachuGM", "PeakachuK"))
venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "#D8A499"
gp[["Face"]][["01"]]$fill <-  "#E6A0C4"
gp[["Face"]][["10"]]$fill <-  "#C6CDF7"
gp$Set$Set1$col <- "#C6CDF7"
gp$Set$Set2$col <- "#E6A0C4"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- .1
gp[["SetText"]][["Set2"]]$cex <- .1
gp$SetText$Set1$col <- "white"
gp$SetText$Set2$col <- "white"
gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

genomicCorr.jaccard(all_true_bounds_gm12878_p,
                    all_true_bounds_k562_p)


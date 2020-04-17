#08 Creating venn diagrams

## Loading packages and functions

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

### Extracting called and predicted boundaries

#### GM12878

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
chrs <- paste0("CHR", c(2:8,10:15,17:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[2]]
  pred_bound_list[[i]] <- called_and_pred[[3]]
  
}

all_true_bounds_gm12878 <- do.call("c", true_bound_list)
all_pred_bounds_gm12878 <- do.call("c", pred_bound_list)

all_true_bounds_gm12878 <- flank(all_true_bounds_gm12878, 5000, both=TRUE)
all_pred_bounds_gm12878 <- flank(all_pred_bounds_gm12878, 5000, both=TRUE)

#### K562

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
chrs <- paste0("CHR", c(2:8,10:15,17:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[2]]
  pred_bound_list[[i]] <- called_and_pred[[3]]
  
}

all_true_bounds_k562 <- do.call("c", true_bound_list)
all_pred_bounds_k562 <- do.call("c", pred_bound_list)

all_true_bounds_k562 <- flank(all_true_bounds_k562, 5000, both=TRUE)
all_pred_bounds_k562 <- flank(all_pred_bounds_k562, 5000, both=TRUE)

##### Within cell line (preciseTAD vs ARROWHEAD)

venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Venn(SetNames=SetNames, Weight=Weight)
}

###### GM12878

res <- makeVennDiagram(Peaks=list(all_pred_bounds_gm12878,
                                  all_true_bounds_gm12878),
                       NameOfPeaks=c("preciseTAD", "ARROWHEAD"))

venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "purple"
gp[["Face"]][["01"]]$fill <-  "lightblue"
gp[["Face"]][["10"]]$fill <-  "red"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- 2
gp[["SetText"]][["Set2"]]$cex <- 2

gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

###### K562

res <- makeVennDiagram(Peaks=list(all_pred_bounds_k562,
                                  all_true_bounds_k562),
                       NameOfPeaks=c("preciseTAD", "ARROWHEAD"))

venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "purple"
gp[["Face"]][["01"]]$fill <-  "lightblue"
gp[["Face"]][["10"]]$fill <-  "red"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- 2
gp[["SetText"]][["Set2"]]$cex <- 2

gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

##### Across cell line 

###### preciseTAD

res <- makeVennDiagram(Peaks=list(all_pred_bounds_gm12878,
                                  all_pred_bounds_k562),
                       NameOfPeaks=c("GM12878", "K562"))

venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "purple"
gp[["Face"]][["01"]]$fill <-  "lightblue"
gp[["Face"]][["10"]]$fill <-  "red"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- 2
gp[["SetText"]][["Set2"]]$cex <- 2

gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

###### ARROWHEAD

res <- makeVennDiagram(Peaks=list(all_true_bounds_gm12878,
                                  all_true_bounds_k562),
                       NameOfPeaks=c("GM12878", "K562"))

venn.pt <- venn_cnt2venn(res$vennCounts)
venn.pt <- compute.Venn(venn.pt)
gp <- VennThemes(venn.pt)
gp[["Face"]][["11"]]$fill <-  "purple"
gp[["Face"]][["01"]]$fill <-  "lightblue"
gp[["Face"]][["10"]]$fill <-  "red"
gp[["FaceText"]][["10"]]$cex <- 2
gp[["FaceText"]][["11"]]$cex <- 2
gp[["FaceText"]][["01"]]$cex <- 2
gp[["SetText"]][["Set1"]]$cex <- 2
gp[["SetText"]][["Set2"]]$cex <- 2

gridExtra::grid.arrange(grid::grid.grabExpr(plot(venn.pt, gp = gp, show=list(Universe=FALSE))), top=textGrob("", gp=gpar(fontsize=50)))

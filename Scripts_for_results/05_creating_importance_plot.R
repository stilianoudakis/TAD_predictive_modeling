#05_creating_importance_plot

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

chrs <- paste0("CHR",c(1:8,10:22))
rfeperchr <- list()
for(i in c("GM12878")){
  for(k in 1:length(chrs)){
    #rfeModelResultsList <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/",
    #                                      i,
    #                                      "/",
    #                                      "5kb/results_by_chr/",
    #                                      chrs[k],
    #                                      "/rus/distance/TADrfe_holdout.rds"))
    
    #/home/stilianoudakisc/TAD_data_analysis/miscellaneous/testing_rfe/k562/
    rfeModelResultsList <- readRDS(paste0("Z:/TAD_data_analysis/miscellaneous/testing_rfe/TADrfe_holdout_",
                                          chrs[k],
                                          ".rds"))
    
    set.size = 8
    rfeModelResults <- rfeModelResultsList[[2]]
    rfe.set <- rfeModelResults[rfeModelResults$Variables==set.size & rfeModelResults$Resample=="Fold1",]
    rfe.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)
    rfe.order <- order(rfe.set[, c("x")], decreasing = TRUE)
    rfe.set <- rfe.set[rfe.order, ]
    rfe.set$CHR <- chrs[k]
    rfeperchr[[k]] <- rfe.set
  }
}
rfeperchr <- do.call(rbind.data.frame, rfeperchr)
rfeperchr <- reshape(rfeperchr, idvar = "Group.1", timevar = "CHR", direction = "wide")
names(rfeperchr) <- c("TFBS", chrs)
rfeperchr[is.na(rfeperchr)] <- 0
rfeperchr <- rfeperchr[order(rowMeans(rfeperchr[,-1], na.rm = TRUE),decreasing = TRUE),]
rfe.set.mat <- as.matrix(rfeperchr[,-1])
rownames(rfe.set.mat) <- rfeperchr$TFBS
rfe.set.mat <- t(apply(rfe.set.mat,1,rev))

rownames(rfe.set.mat)[2] <- "Znf143"
rfe.set.mat <- rfe.set.mat[1:8,]

distance.col = dist(rfe.set.mat, method = "euclidean")
my_palette <- colorRampPalette(c("white", "#00468BFF"))(ncol(rfe.set.mat))

Heatmap(rfe.set.mat,
        col=my_palette,
        border=TRUE,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        show_row_names = TRUE,
        row_names_side="left",
        show_column_names = FALSE,
        column_names_gp=gpar(fontsize = 15),
        row_names_gp=gpar(fontsize = 15),
        column_names_rot = 45,
        column_names_centered =FALSE,
        show_heatmap_legend = FALSE,
        heatmap_legend_param = list(title = "",
                                    legend_direction = "vertical", 
                                    legend_height = unit(3, "in"),
                                    labels_gp = gpar(fontsize = 15)),
        rect_gp = gpar(col = "black", lwd = 1),
        column_dend_side = "top",column_names_side = "top")



chrs <- paste0("CHR",c(1:8,10:22))
rfeperchr <- list()
for(i in c("GM12878")){
  for(k in 1:length(chrs)){
    rfeModelResultsList <- readRDS(paste0("Z:/TAD_data_analysis/",
                                          i,
                                          "/",
                                          "10kb/results_by_chr/",
                                          chrs[k],
                                          "/rus/distance/TADrfe_holdout_peakachu.rds"))
    
    set.size = 8
    rfeModelResults <- rfeModelResultsList[[2]]
    rfe.set <- rfeModelResults[rfeModelResults$Variables==set.size & rfeModelResults$Resample=="Fold1",]
    rfe.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)
    rfe.order <- order(rfe.set[, c("x")], decreasing = TRUE)
    rfe.set <- rfe.set[rfe.order, ]
    rfe.set$CHR <- chrs[k]
    rfeperchr[[k]] <- rfe.set
  }
}
rfeperchr <- do.call(rbind.data.frame, rfeperchr)
rfeperchr <- reshape(rfeperchr, idvar = "Group.1", timevar = "CHR", direction = "wide")
names(rfeperchr) <- c("TFBS", chrs)
rfeperchr[is.na(rfeperchr)] <- 0
rfeperchr <- rfeperchr[order(rowMeans(rfeperchr[,-1], na.rm = TRUE),decreasing = TRUE),]
rfe.set.mat <- as.matrix(rfeperchr[,-1])
rownames(rfe.set.mat) <- rfeperchr$TFBS
rfe.set.mat <- t(apply(rfe.set.mat,1,rev))

rfe.set.mat <- rfe.set.mat[1:8,]

distance.col = dist(rfe.set.mat, method = "euclidean")
my_palette <- colorRampPalette(c("white", "#AD002AFF"))(ncol(rfe.set.mat))

Heatmap(rfe.set.mat,
        col=my_palette,
        border=TRUE,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        show_row_names = TRUE,
        row_names_side="left",
        show_column_names = FALSE,
        column_names_gp=gpar(fontsize = 15),
        row_names_gp=gpar(fontsize = 15),
        column_names_rot = 45,
        column_names_centered =FALSE,
        show_heatmap_legend = FALSE,
        heatmap_legend_param = list(title = "",
                                    legend_direction = "vertical", 
                                    legend_height = unit(3, "in"),
                                    labels_gp = gpar(fontsize = 15)),
        rect_gp = gpar(col = "black", lwd = 1),
        column_dend_side = "top",column_names_side = "top")


###############################################################################################

chrs <- paste0("CHR",c(1:8,10:22))
rfeperchr <- list()
for(i in c("K562")){
  for(k in 1:length(chrs)){
    rfeModelResultsList <- readRDS(paste0("z:/TAD_data_analysis/",
                                          i,
                                          "/",
                                          "5kb/results_by_chr/",
                                          chrs[k],
                                          "/rus/distance/TADrfe_holdout.rds"))
    
    #/home/stilianoudakisc/TAD_data_analysis/miscellaneous/testing_rfe/k562/
    #rfeModelResultsList <- readRDS(paste0("Z:/TAD_data_analysis/miscellaneous/testing_rfe/TADrfe_holdout_",
    #                                      chrs[k],
    #                                      ".rds"))
    
    set.size = 8
    rfeModelResults <- rfeModelResultsList[[2]]
    rfe.set <- rfeModelResults[rfeModelResults$Variables==set.size & rfeModelResults$Resample=="Fold1",]
    rfe.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)
    rfe.order <- order(rfe.set[, c("x")], decreasing = TRUE)
    rfe.set <- rfe.set[rfe.order, ]
    rfe.set$CHR <- chrs[k]
    rfeperchr[[k]] <- rfe.set
  }
}
rfeperchr <- do.call(rbind.data.frame, rfeperchr)
rfeperchr <- reshape(rfeperchr, idvar = "Group.1", timevar = "CHR", direction = "wide")
names(rfeperchr) <- c("TFBS", chrs)
rfeperchr[is.na(rfeperchr)] <- 0
rfeperchr <- rfeperchr[order(rowMeans(rfeperchr[,-1], na.rm = TRUE),decreasing = TRUE),]
rfe.set.mat <- as.matrix(rfeperchr[,-1])
rownames(rfe.set.mat) <- rfeperchr$TFBS
rfe.set.mat <- t(apply(rfe.set.mat,1,rev))

rfe.set.mat <- rfe.set.mat[1:8,]

distance.col = dist(rfe.set.mat, method = "euclidean")
my_palette <- colorRampPalette(c("white", "#00468BFF"))(ncol(rfe.set.mat))

Heatmap(rfe.set.mat,
        col=my_palette,
        border=TRUE,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        show_row_names = TRUE,
        row_names_side="left",
        show_column_names = FALSE,
        column_names_gp=gpar(fontsize = 15),
        row_names_gp=gpar(fontsize = 15),
        column_names_rot = 45,
        column_names_centered =FALSE,
        show_heatmap_legend = FALSE,
        heatmap_legend_param = list(title = "",
                                    legend_direction = "vertical", 
                                    legend_height = unit(3, "in"),
                                    labels_gp = gpar(fontsize = 15)),
        rect_gp = gpar(col = "black", lwd = 1),
        column_dend_side = "top",column_names_side = "top")


chrs <- paste0("CHR",c(1:8,10:22))
rfeperchr <- list()
for(i in c("K562")){
  for(k in 1:length(chrs)){
    rfeModelResultsList <- readRDS(paste0("Z:/TAD_data_analysis/",
                                          i,
                                          "/",
                                          "10kb/results_by_chr/",
                                          chrs[k],
                                          "/rus/distance/TADrfe_holdout_peakachu.rds"))
    
    set.size = 8
    rfeModelResults <- rfeModelResultsList[[2]]
    rfe.set <- rfeModelResults[rfeModelResults$Variables==set.size & rfeModelResults$Resample=="Fold1",]
    rfe.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)
    rfe.order <- order(rfe.set[, c("x")], decreasing = TRUE)
    rfe.set <- rfe.set[rfe.order, ]
    rfe.set$CHR <- chrs[k]
    rfeperchr[[k]] <- rfe.set
  }
}
rfeperchr <- do.call(rbind.data.frame, rfeperchr)
rfeperchr <- reshape(rfeperchr, idvar = "Group.1", timevar = "CHR", direction = "wide")
names(rfeperchr) <- c("TFBS", chrs)
rfeperchr[is.na(rfeperchr)] <- 0
rfeperchr <- rfeperchr[order(rowMeans(rfeperchr[,-1], na.rm = TRUE),decreasing = TRUE),]
rfe.set.mat <- as.matrix(rfeperchr[,-1])
rownames(rfe.set.mat) <- rfeperchr$TFBS
rfe.set.mat <- t(apply(rfe.set.mat,1,rev))

rfe.set.mat <- rfe.set.mat[1:8,]

distance.col = dist(rfe.set.mat, method = "euclidean")
my_palette <- colorRampPalette(c("white", "#AD002AFF"))(ncol(rfe.set.mat))

Heatmap(rfe.set.mat,
        col=my_palette,
        border=TRUE,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        show_row_names = TRUE,
        row_names_side="left",
        show_column_names = FALSE,
        column_names_gp=gpar(fontsize = 15),
        row_names_gp=gpar(fontsize = 15),
        column_names_rot = 45,
        column_names_centered =FALSE,
        show_heatmap_legend = FALSE,
        heatmap_legend_param = list(title = "",
                                    legend_direction = "vertical", 
                                    legend_height = unit(3, "in"),
                                    labels_gp = gpar(fontsize = 15)),
        rect_gp = gpar(col = "black", lwd = 1),
        column_dend_side = "top",column_names_side = "top")

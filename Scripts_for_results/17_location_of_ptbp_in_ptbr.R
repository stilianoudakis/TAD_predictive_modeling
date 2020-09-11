#17_location_of_ptbp_in_ptbr

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

#arrowhead

##gm12878

ptbp_location <- list()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  ptbp_location_chr <- numeric()
  for(j in 1:length(called_and_pred$PTBP)){
    dist_from_left_ptbr <- start(called_and_pred$PTBP[j]) - start(called_and_pred$PTBR[j]) + 1
    w <- width(called_and_pred$PTBR[j])
    
    ptbp_location_chr[j] <- dist_from_left_ptbr/w
  }
  
  ptbp_location[[i]] <- ptbp_location_chr
  
}

ptbp_location_arrow_gm12878 <- unlist(ptbp_location)
summary(ptbp_location_arrow_gm12878)
sd(ptbp_location_arrow_gm12878)


##k562

ptbp_location <- list()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  ptbp_location_chr <- numeric()
  for(j in 1:length(called_and_pred$PTBP)){
    dist_from_left_ptbr <- start(called_and_pred$PTBP[j]) - start(called_and_pred$PTBR[j]) + 1
    w <- width(called_and_pred$PTBR[j])
    
    ptbp_location_chr[j] <- dist_from_left_ptbr/w
  }
  
  ptbp_location[[i]] <- ptbp_location_chr
  
}

ptbp_location_arrow_k562 <- unlist(ptbp_location)
summary(ptbp_location_arrow_k562)
sd(ptbp_location_arrow_k562)

##########################################################

#peakachu

##gm12878

ptbp_location <- list()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  ptbp_location_chr <- numeric()
  for(j in 1:length(called_and_pred$PTBP)){
    dist_from_left_ptbr <- start(called_and_pred$PTBP[j]) - start(called_and_pred$PTBR[j]) + 1
    w <- width(called_and_pred$PTBR[j])
    
    ptbp_location_chr[j] <- dist_from_left_ptbr/w
  }
  
  ptbp_location[[i]] <- ptbp_location_chr
  
}

ptbp_location_peak_gm12878 <- unlist(ptbp_location)
summary(ptbp_location_peak_gm12878)
sd(ptbp_location_peak_gm12878)


##k562

ptbp_location <- list()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  ptbp_location_chr <- numeric()
  for(j in 1:length(called_and_pred$PTBP)){
    dist_from_left_ptbr <- start(called_and_pred$PTBP[j]) - start(called_and_pred$PTBR[j]) + 1
    w <- width(called_and_pred$PTBR[j])
    
    ptbp_location_chr[j] <- dist_from_left_ptbr/w
  }
  
  ptbp_location[[i]] <- ptbp_location_chr
  
}

ptbp_location_peak_k562 <- unlist(ptbp_location)
summary(ptbp_location_peak_k562)
sd(ptbp_location_peak_k562)

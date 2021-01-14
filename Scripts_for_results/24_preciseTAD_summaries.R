#24_preciseTAD_summaries

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


#gm12878

##arrowhead

PTBRWidth_mean_a <- numeric()
PTBRCoverage_mean_a <- numeric()
DistanceBetweenPTBR_mean_a <- numeric()
NumSubRegions_mean_a <- numeric()
SubRegionWidth_mean_a <- numeric()
DistBetweenSubRegions_mean_a <- numeric()

PTBRWidth_sd_a <- numeric()
PTBRCoverage_sd_a <- numeric()
DistanceBetweenPTBR_sd_a <- numeric()
NumSubRegions_sd_a <- numeric()
SubRegionWidth_sd_a <- numeric()
DistBetweenSubRegions_sd_a <- numeric()

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  PTBRWidth_mean_a[i] <- called_and_pred$Summaries$PTBRWidth$mean
  PTBRCoverage_mean_a[i] <- called_and_pred$Summaries$PTBRCoverage$mean
  DistanceBetweenPTBR_mean_a[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$mean
  NumSubRegions_mean_a[i] <- called_and_pred$Summaries$NumSubRegions$mean
  SubRegionWidth_mean_a[i] <- called_and_pred$Summaries$SubRegionWidth$mean
  DistBetweenSubRegions_mean_a[i] <- called_and_pred$Summaries$DistBetweenSubRegions$mean
  
  PTBRWidth_sd_a[i] <- called_and_pred$Summaries$PTBRWidth$sd
  PTBRCoverage_sd_a[i] <- called_and_pred$Summaries$PTBRCoverage$sd
  DistanceBetweenPTBR_sd_a[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$sd
  NumSubRegions_sd_a[i] <- called_and_pred$Summaries$NumSubRegions$sd
  SubRegionWidth_sd_a[i] <- called_and_pred$Summaries$SubRegionWidth$sd
  DistBetweenSubRegions_sd_a[i] <- called_and_pred$Summaries$DistBetweenSubRegions$sd
  
}

## peakachu

PTBRWidth_mean_p <- numeric()
PTBRCoverage_mean_p <- numeric()
DistanceBetweenPTBR_mean_p <- numeric()
NumSubRegions_mean_p <- numeric()
SubRegionWidth_mean_p <- numeric()
DistBetweenSubRegions_mean_p <- numeric()

PTBRWidth_sd_p <- numeric()
PTBRCoverage_sd_p <- numeric()
DistanceBetweenPTBR_sd_p <- numeric()
NumSubRegions_sd_p <- numeric()
SubRegionWidth_sd_p <- numeric()
DistBetweenSubRegions_sd_p <- numeric()

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  PTBRWidth_mean_p[i] <- called_and_pred$Summaries$PTBRWidth$mean
  PTBRCoverage_mean_p[i] <- called_and_pred$Summaries$PTBRCoverage$mean
  DistanceBetweenPTBR_mean_p[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$mean
  NumSubRegions_mean_p[i] <- called_and_pred$Summaries$NumSubRegions$mean
  SubRegionWidth_mean_p[i] <- called_and_pred$Summaries$SubRegionWidth$mean
  DistBetweenSubRegions_mean_p[i] <- called_and_pred$Summaries$DistBetweenSubRegions$mean
  
  PTBRWidth_sd_p[i] <- called_and_pred$Summaries$PTBRWidth$sd
  PTBRCoverage_sd_p[i] <- called_and_pred$Summaries$PTBRCoverage$sd
  DistanceBetweenPTBR_sd_p[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$sd
  NumSubRegions_sd_p[i] <- called_and_pred$Summaries$NumSubRegions$sd
  SubRegionWidth_sd_p[i] <- called_and_pred$Summaries$SubRegionWidth$sd
  DistBetweenSubRegions_sd_p[i] <- called_and_pred$Summaries$DistBetweenSubRegions$sd
  
}

##making table

PTBRWidth_mean_a <- round(mean(PTBRWidth_mean_a),1)
PTBRCoverage_mean_a <- round(mean(PTBRCoverage_mean_a),1)
DistanceBetweenPTBR_mean_a <- round(mean(DistanceBetweenPTBR_mean_a),1)
NumSubRegions_mean_a <- round(mean(NumSubRegions_mean_a),1)
SubRegionWidth_mean_a <- round(mean(SubRegionWidth_mean_a),1)
DistBetweenSubRegions_mean_a <- round(mean(DistBetweenSubRegions_mean_a),1)

PTBRWidth_sd_a <- round(mean(PTBRWidth_sd_a),1)
PTBRCoverage_sd_a <- round(mean(PTBRCoverage_sd_a),1)
DistanceBetweenPTBR_sd_a <- round(mean(DistanceBetweenPTBR_sd_a),1)
NumSubRegions_sd_a <- round(mean(NumSubRegions_sd_a),1)
SubRegionWidth_sd_a <- round(mean(SubRegionWidth_sd_a),1)
DistBetweenSubRegions_sd_a <- round(mean(DistBetweenSubRegions_sd_a),1)
  

PTBRWidth_mean_p <- round(mean(PTBRWidth_mean_p),1)
PTBRCoverage_mean_p <- round(mean(PTBRCoverage_mean_p),1)
DistanceBetweenPTBR_mean_p <- round(mean(DistanceBetweenPTBR_mean_p),1)
NumSubRegions_mean_p <- round(mean(NumSubRegions_mean_p),1)
SubRegionWidth_mean_p <- round(mean(SubRegionWidth_mean_p),1)
DistBetweenSubRegions_mean_p <- round(mean(DistBetweenSubRegions_mean_p),1)

PTBRWidth_sd_p <- round(mean(PTBRWidth_sd_p),1)
PTBRCoverage_sd_p <- round(mean(PTBRCoverage_sd_p),1)
DistanceBetweenPTBR_sd_p <- round(mean(DistanceBetweenPTBR_sd_p),1)
NumSubRegions_sd_p <- round(mean(NumSubRegions_sd_p),1)
SubRegionWidth_sd_p <- round(mean(SubRegionWidth_sd_p),1)
DistBetweenSubRegions_sd_p <- round(mean(DistBetweenSubRegions_sd_p),1)


cbind(c(paste0(PTBRWidth_mean_a, " (", PTBRWidth_sd_a, ")"),
        paste0(PTBRCoverage_mean_a, " (", PTBRCoverage_sd_a, ")"),
        paste0(DistanceBetweenPTBR_mean_a, " (", DistanceBetweenPTBR_sd_a, ")"),
        paste0(NumSubRegions_mean_a, " (", NumSubRegions_sd_a, ")"),
        paste0(SubRegionWidth_mean_a, " (", SubRegionWidth_sd_a, ")"),
        paste0(DistBetweenSubRegions_mean_a, " (", DistBetweenSubRegions_sd_a, ")")),
      c(paste0(PTBRWidth_mean_p, " (", PTBRWidth_sd_p, ")"),
        paste0(PTBRCoverage_mean_p, " (", PTBRCoverage_sd_p, ")"),
        paste0(DistanceBetweenPTBR_mean_p, " (", DistanceBetweenPTBR_sd_p, ")"),
        paste0(NumSubRegions_mean_p, " (", NumSubRegions_sd_p, ")"),
        paste0(SubRegionWidth_mean_p, " (", SubRegionWidth_sd_p, ")"),
        paste0(DistBetweenSubRegions_mean_p, " (", DistBetweenSubRegions_sd_p, ")")))

################################################################################################

#k562

##arrowhead

PTBRWidth_mean_a <- numeric()
PTBRCoverage_mean_a <- numeric()
DistanceBetweenPTBR_mean_a <- numeric()
NumSubRegions_mean_a <- numeric()
SubRegionWidth_mean_a <- numeric()
DistBetweenSubRegions_mean_a <- numeric()

PTBRWidth_sd_a <- numeric()
PTBRCoverage_sd_a <- numeric()
DistanceBetweenPTBR_sd_a <- numeric()
NumSubRegions_sd_a <- numeric()
SubRegionWidth_sd_a <- numeric()
DistBetweenSubRegions_sd_a <- numeric()

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  PTBRWidth_mean_a[i] <- called_and_pred$Summaries$PTBRWidth$mean
  PTBRCoverage_mean_a[i] <- called_and_pred$Summaries$PTBRCoverage$mean
  DistanceBetweenPTBR_mean_a[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$mean
  NumSubRegions_mean_a[i] <- called_and_pred$Summaries$NumSubRegions$mean
  SubRegionWidth_mean_a[i] <- called_and_pred$Summaries$SubRegionWidth$mean
  DistBetweenSubRegions_mean_a[i] <- called_and_pred$Summaries$DistBetweenSubRegions$mean
  
  PTBRWidth_sd_a[i] <- called_and_pred$Summaries$PTBRWidth$sd
  PTBRCoverage_sd_a[i] <- called_and_pred$Summaries$PTBRCoverage$sd
  DistanceBetweenPTBR_sd_a[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$sd
  NumSubRegions_sd_a[i] <- called_and_pred$Summaries$NumSubRegions$sd
  SubRegionWidth_sd_a[i] <- called_and_pred$Summaries$SubRegionWidth$sd
  DistBetweenSubRegions_sd_a[i] <- called_and_pred$Summaries$DistBetweenSubRegions$sd
  
}

## peakachu

PTBRWidth_mean_p <- numeric()
PTBRCoverage_mean_p <- numeric()
DistanceBetweenPTBR_mean_p <- numeric()
NumSubRegions_mean_p <- numeric()
SubRegionWidth_mean_p <- numeric()
DistBetweenSubRegions_mean_p <- numeric()

PTBRWidth_sd_p <- numeric()
PTBRCoverage_sd_p <- numeric()
DistanceBetweenPTBR_sd_p <- numeric()
NumSubRegions_sd_p <- numeric()
SubRegionWidth_sd_p <- numeric()
DistBetweenSubRegions_sd_p <- numeric()

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  PTBRWidth_mean_p[i] <- called_and_pred$Summaries$PTBRWidth$mean
  PTBRCoverage_mean_p[i] <- called_and_pred$Summaries$PTBRCoverage$mean
  DistanceBetweenPTBR_mean_p[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$mean
  NumSubRegions_mean_p[i] <- called_and_pred$Summaries$NumSubRegions$mean
  SubRegionWidth_mean_p[i] <- called_and_pred$Summaries$SubRegionWidth$mean
  DistBetweenSubRegions_mean_p[i] <- called_and_pred$Summaries$DistBetweenSubRegions$mean
  
  PTBRWidth_sd_p[i] <- called_and_pred$Summaries$PTBRWidth$sd
  PTBRCoverage_sd_p[i] <- called_and_pred$Summaries$PTBRCoverage$sd
  DistanceBetweenPTBR_sd_p[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$sd
  NumSubRegions_sd_p[i] <- called_and_pred$Summaries$NumSubRegions$sd
  SubRegionWidth_sd_p[i] <- called_and_pred$Summaries$SubRegionWidth$sd
  DistBetweenSubRegions_sd_p[i] <- called_and_pred$Summaries$DistBetweenSubRegions$sd
  
}

##making table

PTBRWidth_mean_a <- round(mean(PTBRWidth_mean_a),1)
PTBRCoverage_mean_a <- round(mean(PTBRCoverage_mean_a),1)
DistanceBetweenPTBR_mean_a <- round(mean(DistanceBetweenPTBR_mean_a),1)
NumSubRegions_mean_a <- round(mean(NumSubRegions_mean_a),1)
SubRegionWidth_mean_a <- round(mean(SubRegionWidth_mean_a),1)
DistBetweenSubRegions_mean_a <- round(mean(DistBetweenSubRegions_mean_a),1)

PTBRWidth_sd_a <- round(mean(PTBRWidth_sd_a),1)
PTBRCoverage_sd_a <- round(mean(PTBRCoverage_sd_a),1)
DistanceBetweenPTBR_sd_a <- round(mean(DistanceBetweenPTBR_sd_a),1)
NumSubRegions_sd_a <- round(mean(NumSubRegions_sd_a),1)
SubRegionWidth_sd_a <- round(mean(SubRegionWidth_sd_a),1)
DistBetweenSubRegions_sd_a <- round(mean(DistBetweenSubRegions_sd_a),1)


PTBRWidth_mean_p <- round(mean(PTBRWidth_mean_p),1)
PTBRCoverage_mean_p <- round(mean(PTBRCoverage_mean_p),1)
DistanceBetweenPTBR_mean_p <- round(mean(DistanceBetweenPTBR_mean_p),1)
NumSubRegions_mean_p <- round(mean(NumSubRegions_mean_p),1)
SubRegionWidth_mean_p <- round(mean(SubRegionWidth_mean_p),1)
DistBetweenSubRegions_mean_p <- round(mean(DistBetweenSubRegions_mean_p),1)

PTBRWidth_sd_p <- round(mean(PTBRWidth_sd_p),1)
PTBRCoverage_sd_p <- round(mean(PTBRCoverage_sd_p),1)
DistanceBetweenPTBR_sd_p <- round(mean(DistanceBetweenPTBR_sd_p),1)
NumSubRegions_sd_p <- round(mean(NumSubRegions_sd_p),1)
SubRegionWidth_sd_p <- round(mean(SubRegionWidth_sd_p),1)
DistBetweenSubRegions_sd_p <- round(mean(DistBetweenSubRegions_sd_p),1)


cbind(c(paste0(PTBRWidth_mean_a, " (", PTBRWidth_sd_a, ")"),
        paste0(PTBRCoverage_mean_a, " (", PTBRCoverage_sd_a, ")"),
        paste0(DistanceBetweenPTBR_mean_a, " (", DistanceBetweenPTBR_sd_a, ")"),
        paste0(NumSubRegions_mean_a, " (", NumSubRegions_sd_a, ")"),
        paste0(SubRegionWidth_mean_a, " (", SubRegionWidth_sd_a, ")"),
        paste0(DistBetweenSubRegions_mean_a, " (", DistBetweenSubRegions_sd_a, ")")),
      c(paste0(PTBRWidth_mean_p, " (", PTBRWidth_sd_p, ")"),
        paste0(PTBRCoverage_mean_p, " (", PTBRCoverage_sd_p, ")"),
        paste0(DistanceBetweenPTBR_mean_p, " (", DistanceBetweenPTBR_sd_p, ")"),
        paste0(NumSubRegions_mean_p, " (", NumSubRegions_sd_p, ")"),
        paste0(SubRegionWidth_mean_p, " (", SubRegionWidth_sd_p, ")"),
        paste0(DistBetweenSubRegions_mean_p, " (", DistBetweenSubRegions_sd_p, ")")))
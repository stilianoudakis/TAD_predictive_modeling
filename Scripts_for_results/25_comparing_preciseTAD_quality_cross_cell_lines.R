#25_comparing_preciseTAD_quality_cross_cell_lines

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

#arrowhead

## g_on_g

g_on_g_PTBRWidth <- numeric()
g_on_g_PTBRCoverage <- numeric()
g_on_g_DistanceBetweenPTBR <- numeric()
g_on_g_NumSubRegions <- numeric()
g_on_g_SubRegionWidth <- numeric()
g_on_g_DistBetweenSubRegions <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  g_on_g_PTBRWidth[i] <- called_and_pred$Summaries$PTBRWidth$median
  g_on_g_PTBRCoverage[i] <- called_and_pred$Summaries$PTBRCoverage$median
  g_on_g_DistanceBetweenPTBR[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$median
  g_on_g_NumSubRegions[i] <- called_and_pred$Summaries$NumSubRegions$median
  g_on_g_SubRegionWidth[i] <- called_and_pred$Summaries$SubRegionWidth$median
  g_on_g_DistBetweenSubRegions[i] <- called_and_pred$Summaries$DistBetweenSubRegions$median
  
}


## k_on_g

k_on_g_PTBRWidth <- numeric()
k_on_g_PTBRCoverage <- numeric()
k_on_g_DistanceBetweenPTBR <- numeric()
k_on_g_NumSubRegions <- numeric()
k_on_g_SubRegionWidth <- numeric()
k_on_g_DistBetweenSubRegions <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_k562_on_gm12878.rds"))
  
  k_on_g_PTBRWidth[i] <- called_and_pred$Summaries$PTBRWidth$median
  k_on_g_PTBRCoverage[i] <- called_and_pred$Summaries$PTBRCoverage$median
  k_on_g_DistanceBetweenPTBR[i] <- called_and_pred$Summaries$DistanceBetweenPTBR$median
  k_on_g_NumSubRegions[i] <- called_and_pred$Summaries$NumSubRegions$median
  k_on_g_SubRegionWidth[i] <- called_and_pred$Summaries$SubRegionWidth$median
  k_on_g_DistBetweenSubRegions[i] <- called_and_pred$Summaries$DistBetweenSubRegions$median
  
}

## tests

t.test(g_on_g_PTBRWidth,k_on_g_PTBRWidth)
t.test(g_on_g_PTBRCoverage,k_on_g_PTBRCoverage)
t.test(g_on_g_DistanceBetweenPTBR,k_on_g_DistanceBetweenPTBR)
t.test(g_on_g_NumSubRegions,k_on_g_NumSubRegions)
t.test(g_on_g_SubRegionWidth,k_on_g_SubRegionWidth)
t.test(g_on_g_DistBetweenSubRegions,k_on_g_DistBetweenSubRegions)

wilcox.test(g_on_g_PTBRWidth,k_on_g_PTBRWidth)
wilcox.test(g_on_g_PTBRCoverage,k_on_g_PTBRCoverage)
wilcox.test(g_on_g_DistanceBetweenPTBR,k_on_g_DistanceBetweenPTBR)
wilcox.test(g_on_g_NumSubRegions,k_on_g_NumSubRegions)
wilcox.test(g_on_g_SubRegionWidth,k_on_g_SubRegionWidth)
wilcox.test(g_on_g_DistBetweenSubRegions,k_on_g_DistBetweenSubRegions)


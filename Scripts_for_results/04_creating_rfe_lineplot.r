#04 Creating RFE line plot

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

### GM12878

#### Optimal parameters: 5 kb, TFBS, distance type predictors, RUS

rfeModelResultsList <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/WHOLE/rus/distance/TADrfe_3.rds"))
rfeModelResults <- rfeModelResultsList[[1]]
rfeModelResults <- rfeModelResults[,c(1,11)]

ggplot(rfeModelResults, aes(x=Variables, y=Accuracy)) + 
  #geom_errorbar(aes(ymin=AveMCC-AveMCCSD, ymax=AveMCC+AveMCCSD), 
  #              width=2, 
  #              position=position_dodge(2),
  #              size=1) +
  geom_line(size=2, color="red") +
  geom_point(size=4, color="red") +
  ylim(.6,.8) +
  xlab("Number of top ranked annotations") +
  ylab("Cross-Validated Accuracy") +
  scale_color_discrete(name="Cell line")+
  scale_x_continuous(breaks = c(2,4,8,16,32,52),
                     labels = c(2,4,8,16,32,52))+
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")

#######################################################################

### K562

#### Optimal parameters: 5 kb, TFBS, distance type predictors, RUS

rfeModelResultsList <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/WHOLE/rus/distance/TADrfe_3.rds"))
rfeModelResults <- rfeModelResultsList[[1]]
rfeModelResults <- rfeModelResults[,c(1,11)]

ggplot(rfeModelResults, aes(x=Variables, y=Accuracy)) + 
  #geom_errorbar(aes(ymin=AveMCC-AveMCCSD, ymax=AveMCC+AveMCCSD), 
  #              width=2, 
  #              position=position_dodge(2),
  #              size=1) +
  geom_line(size=2, color="red") +
  geom_point(size=4, color="red") +
  ylim(.6,.8) +
  xlab("Number of top ranked annotations") +
  ylab("Cross-Validated Accuracy") +
  scale_color_discrete(name="Cell line")+
  scale_x_continuous(breaks = c(2,4,8,16,32,52),
                     labels = c(2,4,8,16,32,52))+
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


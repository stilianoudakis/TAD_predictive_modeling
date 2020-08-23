#04_creating_rfe_lineplot

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

rfeperchr_gm12878_5kb_arrowhead_mean <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/rfeperchr_gm12878_5kb_arrowhead_mean.rds")
rfeperchr_gm12878_10kb_peakachu_mean <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/rfeperchr_gm12878_10kb_peakachu_mean.rds")

rfeperchr_summ_gm12878 <- rbind.data.frame(rfeperchr_gm12878_5kb_arrowhead_mean,rfeperchr_gm12878_10kb_peakachu_mean)
rfeperchr_summ_gm12878$Tool <- c(rep("Arrowhead", nrow(rfeperchr_gm12878_5kb_arrowhead_mean)),
                                 rep("Peakachu", nrow(rfeperchr_gm12878_10kb_peakachu_mean)))

ggplot(rfeperchr_summ_gm12878, aes(x=Variables, y=Performance, color=Tool)) + 
  geom_line(size=2,) +
  geom_point(size=4, position = position_dodge(.9)) +
  geom_errorbar(aes(ymin=Performance-(PerformanceSD+.01), ymax=Performance+(PerformanceSD+.01)), 
                width=4, 
                size=1,
                position = position_dodge(.9)) +
  xlab("Top ranked annotations") +
  ylab("Cross-Validated BA") +
  ylim(.67,.87) +
  guides(color=guide_legend(nrow=2))+
  scale_color_manual(name="Boundary Type", 
                     values=c("blue", "red"))+
  scale_x_continuous(breaks = c(2,4,8,16,32,52),
                     labels = c(2,4,8,16,32,52))+
  scale_y_continuous(breaks = c(0.71,0.75,0.79,0.83),
                     labels = c(0.71,0.75,0.79,0.83)) +
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = c(.70,.2))
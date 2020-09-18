#22_evaluating_normalized_enrichment

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

pt_ne_grid <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/pt_ne_grid.rds")

pt_ne_grid_peakachu <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/pt_ne_grid_peakachu.rds")

pt_ne_grid_k562 <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/pt_ne_grid_k562.rds")

pt_ne_grid_peakachu_k562 <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/pt_ne_grid_peakachu_k562.rds")

grid_df <- rbind.data.frame(pt_ne_grid$Summaries$ParameterGrid,
                            pt_ne_grid_peakachu$Summaries$ParameterGrid,
                            pt_ne_grid_k562$Summaries$ParameterGrid,
                            pt_ne_grid_peakachu_k562$Summaries$ParameterGrid)
names(grid_df) <- c("Threshold","Epsilon", "NE", "NEsd")
grid_df$CL <- c(rep("GM12878", nrow(pt_ne_grid$Summaries$ParameterGrid)+nrow(pt_ne_grid_peakachu$Summaries$ParameterGrid)),
                rep("K562", nrow(pt_ne_grid_k562$Summaries$ParameterGrid)+nrow(pt_ne_grid_peakachu_k562$Summaries$ParameterGrid)))
grid_df$Tool <- c(rep("Arrowhead", nrow(pt_ne_grid$Summaries$ParameterGrid)),
                  rep("Peakachu", nrow(pt_ne_grid_peakachu$Summaries$ParameterGrid)),
                  rep("Arrowhead", nrow(pt_ne_grid_k562$Summaries$ParameterGrid)),
                  rep("Peakachu", nrow(pt_ne_grid_peakachu_k562$Summaries$ParameterGrid)))
grid_df$Threshold <- as.factor(as.character(grid_df$Threshold))
grid_df$Epsilon <- factor(as.character(grid_df$Epsilon), levels=c("1000",
                                                                     "5000",
                                                                     "10000",
                                                                     "15000",
                                                                     "20000",
                                                                     "25000"))

ggplot(grid_df, aes(x=Epsilon, y=NE, color=Tool, group=Tool)) +
  geom_line(size=1.5, position=position_dodge(0.5))+
  geom_point(size=4, position=position_dodge(0.5))+
  geom_errorbar(aes(ymin=NE-NEsd, ymax=NE+NEsd), width=1,
                position=position_dodge(0.5)) +
  facet_grid(CL ~ Threshold)+
  theme_minimal() +
  theme_bw()+
  xlab("Epsilon") + 
  ylab("Normalized Enrichment") + #ylim(-.1,.5) +
  scale_color_manual(name="Ground Truth",
                       values=c("#00468BFF", "#AD002AFF"))+
  guides(shape=FALSE, linetype=FALSE)+
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  strip.text.x = element_text(size = 15),
  strip.text.y = element_text(size = 15),
  #panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.text=element_text(size=15),
  legend.title=element_text(size=15),
  plot.title = element_text(size=15),
  legend.position="bottom")

#03_creating_annotation_specific_barplots

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

badata_gm12878_5kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/badata_gm12878_5kb_arrowhead_annotations.rds")

all_rfperfs <- data.frame(Perf=badata_gm12878_5kb)
all_rfperfs$Resolution <- rep("5 kb", length(badata_gm12878_5kb))
all_rfperfs$Annotation <- c("BroadHMM","HM","TFBS","ALL")
all_rfperfs$Annotation <- factor(all_rfperfs$Annotation, levels = c("HM", "BroadHMM", "TFBS", "ALL"))
all_rfperfs$Chromosome <- c(rep("CHR1", 4),
                            rep("CHR2", 4),
                            rep("CHR3", 4),
                            rep("CHR4", 4),
                            rep("CHR5", 4),
                            rep("CHR6", 4),
                            rep("CHR7", 4),
                            rep("CHR8", 4),
                            rep("CHR10", 4),
                            rep("CHR11", 4),
                            rep("CHR12", 4),
                            rep("CHR13", 4),
                            rep("CHR14", 4),
                            rep("CHR15", 4),
                            rep("CHR16", 4),
                            rep("CHR17", 4),
                            rep("CHR18", 4),
                            rep("CHR19", 4),
                            rep("CHR20", 4),
                            rep("CHR21", 4),
                            rep("CHR22", 4))

all_rfperfs <- all_rfperfs %>% 
  group_by(Resolution, Annotation) %>% 
  dplyr::summarise(Performance = mean(Perf, na.rm = TRUE),
                   PerformanceSD = sd(Perf, na.rm = TRUE))

all_rfperfs_a <- all_rfperfs


badata_gm12878_10kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/badata_gm12878_10kb_peakachu_annotations.rds")

all_rfperfs <- data.frame(Perf=badata_gm12878_10kb)
all_rfperfs$Resolution <- rep("10 kb", length(badata_gm12878_10kb))
all_rfperfs$Annotation <- c("BroadHMM","HM","TFBS","ALL")
all_rfperfs$Annotation <- factor(all_rfperfs$Annotation, levels = c("HM", "BroadHMM", "TFBS", "ALL"))
all_rfperfs$Chromosome <- c(rep("CHR1", 4),
                            rep("CHR2", 4),
                            rep("CHR3", 4),
                            rep("CHR4", 4),
                            rep("CHR5", 4),
                            rep("CHR6", 4),
                            rep("CHR7", 4),
                            rep("CHR8", 4),
                            rep("CHR10", 4),
                            rep("CHR11", 4),
                            rep("CHR12", 4),
                            rep("CHR13", 4),
                            rep("CHR14", 4),
                            rep("CHR15", 4),
                            rep("CHR16", 4),
                            rep("CHR17", 4),
                            rep("CHR18", 4),
                            rep("CHR19", 4),
                            rep("CHR20", 4),
                            rep("CHR21", 4),
                            rep("CHR22", 4))

all_rfperfs <- all_rfperfs %>% 
  group_by(Resolution, Annotation) %>% 
  dplyr::summarise(Performance = mean(Perf, na.rm = TRUE),
                   PerformanceSD = sd(Perf, na.rm = TRUE))

all_rfperfs_p <- all_rfperfs

all_rfperfs_a_p <- rbind.data.frame(all_rfperfs_a,all_rfperfs_p)
all_rfperfs_a_p$Tool <- c(rep("Arrowhead", nrow(all_rfperfs_a)),
                          rep("Peakachu", nrow(all_rfperfs_p)))

ggplot(all_rfperfs_a_p, aes(x=Annotation, y=Performance, fill=Tool)) +
  facet_grid(. ~ Tool) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Performance-PerformanceSD, ymax=Performance+PerformanceSD), width=.75, size=1) +
  guides(fill=FALSE) +
  scale_fill_manual(values = c("blue", "red"))+
  theme_minimal() +
  theme_bw()+
  xlab("Genomic Class") + 
  ylab("Balanced Accuracy") + 
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  strip.text.x = element_text(size = 15),
  legend.text=element_text(size=15),
  legend.title=element_text(size=15),
  legend.position = "bottom",
  plot.title = element_text(size=15))
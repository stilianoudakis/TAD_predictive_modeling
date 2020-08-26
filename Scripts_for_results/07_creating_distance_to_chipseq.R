#07_creating_distance_to_chipseq

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

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS",
                                             pattern="*.bed",
                                             signal=4)

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[3]]
  pred_bound_list[[i]] <- called_and_pred[[2]]
  
}

all_true_bounds <- do.call("c", true_bound_list)
all_pred_bounds <- do.call("c", pred_bound_list)

pred_v_called_df_a <- data.frame(LogDist = c(#preciseTAD
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1),
  
  #ARROWHEAD
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)
),
Annotation = c(#preciseTAD
  rep("CTCF",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("RAD21",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("SMC3",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ZNF143",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1))),
  
  #ARROWHEAD
  rep("CTCF",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("RAD21",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("SMC3",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ZNF143",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))),
BoundReg = c(#preciseTAD
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1))),
  
  #ARROWHEAD
  rep("Arrowhead", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("Arrowhead", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("Arrowhead", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("Arrowhead", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))))

pred_v_called_df_a$BoundReg <- factor(pred_v_called_df_a$BoundReg, levels=c("Arrowhead", "preciseTAD"))

#######################

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[3]]
  pred_bound_list[[i]] <- called_and_pred[[2]]
  
}

all_true_bounds <- do.call("c", true_bound_list)
all_pred_bounds <- do.call("c", pred_bound_list)

pred_v_called_df_p <- data.frame(LogDist = c(#preciseTAD
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1),
  
  #PEAKACHU
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)
),
Annotation = c(#preciseTAD
  rep("CTCF",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("RAD21",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("SMC3",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ZNF143",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1))),
  
  #PEAKACHU
  rep("CTCF",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("RAD21",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("SMC3",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ZNF143",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))),
BoundReg = c(#preciseTAD
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1))),
  
  #PEAKACHU
  rep("Peakachu", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("Peakachu", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("Peakachu", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("Peakachu", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))))

pred_v_called_df_p$BoundReg <- factor(pred_v_called_df_p$BoundReg, levels=c("Peakachu", "preciseTAD"))

######################

#plotting

a <- ggplot(pred_v_called_df_a[which(pred_v_called_df_a$Annotation=="CTCF"),], aes(x=BoundReg, y = LogDist, fill=BoundReg))  +   
  stat_boxplot(geom ='errorbar', width = 0.2, size=1.2) + 
  geom_boxplot(outlier.shape = NA, color="black", size=1.2) +
  geom_signif(test = "wilcox.test", 
              comparisons = list(c("preciseTAD","Arrowhead")),
              vjust = 0,
              textsize = 5,
              size = .5,
              #step_increase = .5,
              color="black") +
  theme_minimal()+
  theme_bw()+
  ylab("Log2 Distance to CTCF")+
  xlab("") +
  scale_fill_manual(values=c("blue",
                             "forestgreen"))+
  scale_color_manual(values=c("blue",
                              "forestgreen")) +
  guides(color=FALSE, fill=FALSE)+
  theme(axis.text.x = element_text(size=20,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")

b <- ggplot(pred_v_called_df_p[which(pred_v_called_df_p$Annotation=="CTCF"),], aes(x=BoundReg, y = LogDist, fill=BoundReg, color=BoundReg))  +   
  stat_boxplot(geom ='errorbar', width = 0.2, size=1.2, color="black") + 
  geom_boxplot(outlier.shape = NA, color="black", size=1.2) +
  geom_signif(test = "wilcox.test", 
              comparisons = list(c("preciseTAD","Peakachu")),
              vjust = 0,
              textsize = 5,
              size = .5,
              #step_increase = .5,
              color="black") +
  theme_minimal()+
  theme_bw()+
  ylab("")+
  xlab("") +
  scale_fill_manual(values=c("red",
                             "forestgreen"))+
  scale_color_manual(values=c("red",
                              "forestgreen")) +
  guides(color=FALSE, fill=FALSE)+
  theme(axis.text.x = element_text(size=20,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")

ggarrange(a,b,ncol = 2)

################################################################################################

#k562 

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/k562/topTFBS",
                                             pattern="*.bed",
                                             signal=4)

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[3]]
  pred_bound_list[[i]] <- called_and_pred[[2]]
  
}

all_true_bounds <- do.call("c", true_bound_list)
all_pred_bounds <- do.call("c", pred_bound_list)

pred_v_called_df_a <- data.frame(LogDist = c(#preciseTAD
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1),
  
  #ARROWHEAD
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)
),
Annotation = c(#preciseTAD
  rep("CTCF",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("RAD21",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("SMC3",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ZNF143",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1))),
  
  #ARROWHEAD
  rep("CTCF",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("RAD21",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("SMC3",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ZNF143",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))),
BoundReg = c(#preciseTAD
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1))),
  
  #ARROWHEAD
  rep("Arrowhead", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("Arrowhead", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("Arrowhead", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("Arrowhead", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))))

pred_v_called_df_a$BoundReg <- factor(pred_v_called_df_a$BoundReg, levels=c("Arrowhead", "preciseTAD"))

#######################

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[3]]
  pred_bound_list[[i]] <- called_and_pred[[2]]
  
}

all_true_bounds <- do.call("c", true_bound_list)
all_pred_bounds <- do.call("c", pred_bound_list)

pred_v_called_df_p <- data.frame(LogDist = c(#preciseTAD
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1),
  log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1),
  
  #PEAKACHU
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1),
  log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)
),
Annotation = c(#preciseTAD
  rep("CTCF",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("RAD21",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("SMC3",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ZNF143",length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1))),
  
  #PEAKACHU
  rep("CTCF",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("RAD21",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("SMC3",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ZNF143",length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))),
BoundReg = c(#preciseTAD
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("preciseTAD", length(log2(mcols(distanceToNearest(all_pred_bounds, genomicElements.GR[[4]]))$distance+1))),
  
  #PEAKACHU
  rep("Peakachu", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("Peakachu", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("Peakachu", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("Peakachu", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))))

pred_v_called_df_p$BoundReg <- factor(pred_v_called_df_p$BoundReg, levels=c("Peakachu", "preciseTAD"))


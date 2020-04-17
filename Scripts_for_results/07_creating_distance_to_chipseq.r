#07 Creating distance to chip-seq plots

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

#### Create granges from top tfbs

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS",
                                             pattern="*.bed",
											 signal=4)

##### Extracting called and predicted boundaries

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[2]]
  pred_bound_list[[i]] <- called_and_pred[[3]]
  
}

all_true_bounds <- do.call("c", true_bound_list)
all_pred_bounds <- do.call("c", pred_bound_list)

#all_true_bounds <- flank(all_true_bounds, 2500, both=TRUE)
#all_pred_bounds <- flank(all_pred_bounds, 2500, both=TRUE)


###### Distance plots

pred_v_called_df <- data.frame(LogDist = c(#preciseTAD
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
  rep("ARROWHEAD", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("ARROWHEAD", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("ARROWHEAD", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ARROWHEAD", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))))

pred_v_called_df$BoundReg <- factor(pred_v_called_df$BoundReg, levels=c("preciseTAD", "ARROWHEAD"))

ggplot(pred_v_called_df, aes(x=BoundReg, y = LogDist, fill=BoundReg, color=BoundReg))  + 
  facet_grid(~ Annotation) +    
  stat_boxplot(geom ='errorbar', width = 0.2) + 
  geom_boxplot(outlier.shape = NA) +
  geom_signif(test = "wilcox.test", 
              comparisons = list(c("preciseTAD","ARROWHEAD")),
              vjust = 0,
              textsize = 4,
              size = .5,
              step_increase = .08,
              color="black") +
  theme_minimal()+
  theme_bw()+
  ylab("Log2 Distance to Annotation")+
  xlab("") +
  scale_fill_manual(name="Boundary",
                    values=c("forestgreen",
                             "blue"))+
  scale_color_manual(values=c("black",
                              "black")) +
  guides(color=FALSE, linetype=FALSE)+
  theme(axis.text.x = element_text(size=15,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")

###########################################################################

### K562

#### Create granges from top tfbs

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/k562/topTFBS",
                                             pattern="*.bed",
											 signal=4)

##### Extracting called and predicted boundaries

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
chrs <- paste0("CHR", c(2:8,10:15,17:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[2]]
  pred_bound_list[[i]] <- called_and_pred[[3]]
  
}

all_true_bounds <- do.call("c", true_bound_list)
all_pred_bounds <- do.call("c", pred_bound_list)

#all_true_bounds <- flank(all_true_bounds, 5000, both=TRUE)
#all_pred_bounds <- flank(all_pred_bounds, 5000, both=TRUE)


###### Distance plots

pred_v_called_df <- data.frame(LogDist = c(#preciseTAD
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
  rep("ARROWHEAD", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[1]]))$distance+1))),
  rep("ARROWHEAD", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[2]]))$distance+1))),
  rep("ARROWHEAD", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[3]]))$distance+1))),
  rep("ARROWHEAD", length(log2(mcols(distanceToNearest(all_true_bounds, genomicElements.GR[[4]]))$distance+1)))))

pred_v_called_df$BoundReg <- factor(pred_v_called_df$BoundReg, levels=c("preciseTAD", "ARROWHEAD"))

ggplot(pred_v_called_df, aes(x=BoundReg, y = LogDist, fill=BoundReg, color=BoundReg))  + 
  facet_grid(~ Annotation) +    
  stat_boxplot(geom ='errorbar', width = 0.2) + 
  geom_boxplot(outlier.shape = NA) +
  geom_signif(test = "wilcox.test", 
              comparisons = list(c("preciseTAD","ARROWHEAD")),
              vjust = 0,
              textsize = 4,
              size = .5,
              step_increase = .08,
              color="black") +
  theme_minimal()+
  theme_bw()+
  ylab("Log2 Distance to Annotation")+
  xlab("") +
  scale_fill_manual(name="Boundary",
                    values=c("forestgreen",
                             "blue"))+
  scale_color_manual(values=c("black",
                              "black")) +
  guides(color=FALSE, linetype=FALSE)+
  theme(axis.text.x = element_text(size=15,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")

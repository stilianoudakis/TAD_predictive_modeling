library(GenomicRanges)
library(caret)
library(data.table)
library(gbm)
library(randomForest)
library(glmnet)
library(pROC)
library(plyr)
library(dplyr)
library(ggplot2)
library(DMwR)
library(gridExtra)
library(pROC)
library(ROCR)
#library(leaps)
library(cluster)
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

# GM12878; rus; distance types

## 5kb

mccdata_gm12878_5kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
    if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF_",
                                 j, 
                                 ".rds"))
      perfdata <- perfdata[[3]]
    }else{
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF.rds"))
      perfdata <- perfdata[[3]]
    }
    
    mccdata_gm12878_5kb <- rbind(mccdata_gm12878_5kb,
                                 perfdata[15,2])
    
  }
}

## 10 kb

mccdata_gm12878_10kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
    if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF_",
                                 j, 
                                 ".rds"))
      perfdata <- perfdata[[3]]
    }else{
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF.rds"))
      perfdata <- perfdata[[3]]
    }
    
    mccdata_gm12878_10kb <- rbind(mccdata_gm12878_10kb,
                                  perfdata[15,2])
    
  }
}

## 25 kb

mccdata_gm12878_25kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
    if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/25kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF_",
                                 j, 
                                 ".rds"))
      perfdata <- perfdata[[3]]
    }else{
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/25kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF.rds"))
      perfdata <- perfdata[[3]]
    }
    
    mccdata_gm12878_25kb <- rbind(mccdata_gm12878_25kb,
                                  perfdata[15,2])
    
  }
}

## 50 kb

mccdata_gm12878_50kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
    if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/50kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF_",
                                 j, 
                                 ".rds"))
      perfdata <- perfdata[[3]]
    }else{
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/50kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF.rds"))
      perfdata <- perfdata[[3]]
    }
    
    mccdata_gm12878_50kb <- rbind(mccdata_gm12878_50kb,
                                  perfdata[15,2])
    
  }
}

## 100 kb

mccdata_gm12878_100kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
    if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/100kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF_",
                                 j, 
                                 ".rds"))
      perfdata <- perfdata[[3]]
    }else{
      perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/100kb/results_by_chr/",
                                 i, 
                                 "/rus/distance/TADRF.rds"))
      perfdata <- perfdata[[3]]
    }
    
    mccdata_gm12878_100kb <- rbind(mccdata_gm12878_100kb,
                                   perfdata[15,2])
    
  }
}

all_rfperfs <- rbind.data.frame(mccdata_gm12878_5kb,
                                mccdata_gm12878_10kb,
                                mccdata_gm12878_25kb,
                                mccdata_gm12878_50kb,
                                mccdata_gm12878_100kb)
all_rfperfs$Resolution <- c(rep("5 kb", length(mccdata_gm12878_5kb)),
                            rep("10 kb", length(mccdata_gm12878_10kb)),
                            rep("25 kb", length(mccdata_gm12878_25kb)),
                            rep("50 kb", length(mccdata_gm12878_50kb)),
                            rep("100 kb", length(mccdata_gm12878_100kb)))
all_rfperfs$Annotation <- c("CRE","HM","TFBS","ALL")
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
all_rfperfs$Resolution <- factor(all_rfperfs$Resolution, levels = c("5 kb", "10 kb", "25 kb", "50 kb", "100 kb"))
all_rfperfs$Annotation <- factor(all_rfperfs$Annotation, levels = c("HM", "CRE", "TFBS", "ALL"))
all_rfperfs$Chromosome <- factor(all_rfperfs$Chromosome, levels = paste0("CHR", c(1:8,10:22)))

all_rfperfs_mcc <- all_rfperfs %>%
  dplyr::group_by(Resolution, Annotation) %>%
  dplyr::summarise(MCC = mean(V1, na.rm = TRUE),
                   MCCSD = sd(V1, na.rm = TRUE))

ggplot(all_rfperfs_mcc, aes(x=Resolution, y=MCC, group=Annotation, fill=Annotation)) +
  geom_bar(stat="identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=MCC-MCCSD, ymax=MCC+MCCSD), width=.5,
                position=position_dodge(.9), alpha=0.75) +
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + ylab("MCC") + 
  scale_fill_discrete(name="Genomic Class")+
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20),
  legend.position = "bottom",
  plot.title = element_text(size=20))

####################################################################################################

# K562; rus; distance types
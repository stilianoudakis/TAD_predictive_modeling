#03 Creating annotation specific barplots

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

### GM12878 (rus; distance types)

#### 5kb

badata_gm12878_5kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_gm12878_5kb <- rbind(badata_gm12878_5kb,
                               perfdata[10,2])
  
}

#### 10kb

badata_gm12878_10kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/10kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/10kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_gm12878_10kb <- rbind(badata_gm12878_10kb,
                               perfdata[10,2])
  
}

#### 25kb

badata_gm12878_25kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/25kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/25kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_gm12878_25kb <- rbind(badata_gm12878_25kb,
                               perfdata[10,2])
  
}

#### 50kb

badata_gm12878_50kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/50kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/50kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_gm12878_50kb <- rbind(badata_gm12878_50kb,
                               perfdata[10,2])
  
}

#### 100kb

badata_gm12878_100kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/100kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/100kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_gm12878_100kb <- rbind(badata_gm12878_100kb,
                               perfdata[10,2])
  
}

all_rfperfs <- rbind.data.frame(badata_gm12878_5kb,
                                badata_gm12878_10kb,
                                badata_gm12878_25kb,
                                badata_gm12878_50kb,
                                badata_gm12878_100kb)
all_rfperfs$Resolution <- c(rep("5 kb", length(badata_gm12878_5kb)),
                            rep("10 kb", length(badata_gm12878_10kb)),
                            rep("25 kb", length(badata_gm12878_25kb)),
                            rep("50 kb", length(badata_gm12878_50kb)),
                            rep("100 kb", length(badata_gm12878_100kb)))
all_rfperfs$Annotation <- c("CRE","HM","TFBS","ALL")
all_rfperfs$Resolution <- factor(all_rfperfs$Resolution, levels = c("5 kb", "10 kb", "25 kb", "50 kb", "100 kb"))
all_rfperfs$Annotation <- factor(all_rfperfs$Annotation, levels = c("HM", "CRE", "TFBS", "ALL"))

ggplot(all_rfperfs, aes(x=Resolution, y=V1, group=Annotation, fill=Annotation)) +
  geom_bar(stat="identity", position = position_dodge()) +
  #geom_errorbar(aes(ymin=BA-BASD, ymax=BA+BASD), width=.5,
  #              position=position_dodge(.9), alpha=0.75) +
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + ylab("BA") + 
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

# K562 (rus; distance types)

## 5kb

badata_k562_5kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/5kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/5kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_k562_5kb <- rbind(badata_k562_5kb,
                               perfdata[10,2])
  
}

## 10kb

badata_k562_10kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/10kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/10kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_k562_10kb <- rbind(badata_k562_10kb,
                               perfdata[10,2])
  
}

## 25kb

badata_k562_25kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/25kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/25kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_k562_25kb <- rbind(badata_k562_25kb,
                               perfdata[10,2])
  
}

## 50kb

badata_k562_50kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/50kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/50kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_k562_50kb <- rbind(badata_k562_50kb,
                               perfdata[10,2])
  
}

## 100kb

badata_k562_100kb <- numeric()

for(j in c("BroadHMM", "HistoneModifications", "TFBS", "all")){
  if(j %in% c("BroadHMM", "HistoneModifications", "TFBS")){
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/100kb/results_by_chr/WHOLE/rus/distance/TADRF_",
                               j, 
                               ".rds"))
    perfdata <- perfdata[[3]]
  }else{
    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/100kb/results_by_chr/WHOLE/rus/distance/TADRF_3.rds"))
    perfdata <- perfdata[[3]]
  }
  
  badata_k562_100kb <- rbind(badata_k562_100kb,
                               perfdata[10,2])
  
}

all_rfperfs <- rbind.data.frame(badata_k562_5kb,
                                badata_k562_10kb,
                                badata_k562_25kb,
                                badata_k562_50kb,
                                badata_k562_100kb)
all_rfperfs$Resolution <- c(rep("5 kb", length(badata_k562_5kb)),
                            rep("10 kb", length(badata_k562_10kb)),
                            rep("25 kb", length(badata_k562_25kb)),
                            rep("50 kb", length(badata_k562_50kb)),
                            rep("100 kb", length(badata_k562_100kb)))
all_rfperfs$Annotation <- c("CRE","HM","TFBS","ALL")
all_rfperfs$Resolution <- factor(all_rfperfs$Resolution, levels = c("5 kb", "10 kb", "25 kb", "50 kb", "100 kb"))
all_rfperfs$Annotation <- factor(all_rfperfs$Annotation, levels = c("HM", "CRE", "TFBS", "ALL"))

ggplot(all_rfperfs, aes(x=Resolution, y=V1, group=Annotation, fill=Annotation)) +
  geom_bar(stat="identity", position = position_dodge()) +
  #geom_errorbar(aes(ymin=BA-BASD, ymax=BA+BASD), width=.5,
  #              position=position_dodge(.9), alpha=0.75) +
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + ylab("BA") + 
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



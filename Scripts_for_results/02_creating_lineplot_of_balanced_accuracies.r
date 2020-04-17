#02 Creating line plot of balanced accuracies

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

####5kb

badata_gm12878_5kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_gm12878_5kb <- rbind(badata_gm12878_5kb,
                                   perfdata[10,2])
    }
  }
    
####10kb

badata_gm12878_10kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_gm12878_10kb <- rbind(badata_gm12878_10kb,
                                   perfdata[10,2])
    }
  }
  
####25kb

badata_gm12878_25kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/25kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/25kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_gm12878_25kb <- rbind(badata_gm12878_25kb,
                                   perfdata[10,2])
    }
  }
  
####50kb

badata_gm12878_50kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/50kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/50kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_gm12878_50kb <- rbind(badata_gm12878_50kb,
                                   perfdata[10,2])
    }
  }
  
####100kb

badata_gm12878_100kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/100kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/100kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_gm12878_100kb <- rbind(badata_gm12878_100kb,
                                   perfdata[10,2])
    }
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
all_rfperfs$Resampling <- c(rep("None", 4),
                            rep("ROS", 4),
                            rep("RUS", 4),
                            rep("SMOTE", 4))
all_rfperfs$Predictor <- c("Distance", "OC", "OP", "Signal")

ggplot(all_rfperfs, aes(x=Resolution, y=Performance, color=Resampling, group=Resampling)) +
  geom_line(size=2, position=position_dodge(0.5))+
  geom_point(size=5, position=position_dodge(0.5))+
  #geom_errorbar(aes(ymin=Performance-PerformanceSD, ymax=Performance+PerformanceSD), width=1,
  #              position=position_dodge(0.5)) +
  facet_grid(. ~ Predictor)+
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + 
  ylab("Balanced Accuracy") + 
  #ylim(-.1,.5) +
  scale_color_discrete(name="Resampling Technique")+
  guides(shape=FALSE, linetype=FALSE)+
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  strip.text.x = element_text(size = 15),
  strip.text.y = element_text(size = 15),
  #panel.spacing = unit(2, "lines"),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20),
  plot.title = element_text(size=20),
  legend.position="bottom")

########################################################################################

### K562

####5kb

badata_k562_5kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_k562_5kb <- rbind(badata_k562_5kb,
                                   perfdata[10,2])
    }
  }
    
####10kb

badata_k562_10kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_k562_10kb <- rbind(badata_k562_10kb,
                                   perfdata[10,2])
    }
  }
  
####25kb

badata_k562_25kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/25kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/25kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_k562_25kb <- rbind(badata_k562_25kb,
                                   perfdata[10,2])
    }
  }
  
####50kb

badata_k562_50kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/50kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/50kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_k562_50kb <- rbind(badata_k562_50kb,
                                   perfdata[10,2])
    }
  }
  
####100kb

badata_k562_100kb <- numeric()
for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op", "signal")){
      print(c(j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/100kb/results_by_chr/",
                                 "WHOLE", 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF_3.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/100kb/results_by_chr/",
                                   "WHOLE", 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF_3.rds"))
        perfdata <- perfdata[[3]]						   
      }else{perfdata=data.frame(Metric = c("TN",
                                           "FN",
                                           "FP",
                                           "TP",
                                           "Total",
                                           "Sensitivity",
                                           "Specificity",
                                           "Kappa",
                                           "Accuracy",
										   "BalancedAccuracy",
                                           "Precision",
                                           "FPR",
                                           "FNR",
                                           "FOR",
                                           "NPV",
                                           "MCC",
                                           "F1",
                                           "AUC",
                                           "Youden",
                                           "AUPRC"), 
                                Performance=NA)}
      
      badata_k562_100kb <- rbind(badata_k562_100kb,
                                   perfdata[10,2])
    }
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
all_rfperfs$Resampling <- c(rep("None", 4),
                            rep("ROS", 4),
                            rep("RUS", 4),
                            rep("SMOTE", 4))
all_rfperfs$Predictor <- c("Distance", "OC", "OP", "Signal")

ggplot(all_rfperfs, aes(x=Resolution, y=Performance, color=Resampling, group=Resampling)) +
  geom_line(size=2, position=position_dodge(0.5))+
  geom_point(size=5, position=position_dodge(0.5))+
  #geom_errorbar(aes(ymin=Performance-PerformanceSD, ymax=Performance+PerformanceSD), width=1,
  #              position=position_dodge(0.5)) +
  facet_grid(. ~ Predictor)+
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + 
  ylab("Balanced Accuracy") + 
  #ylim(-.1,.5) +
  scale_color_discrete(name="Resampling Technique")+
  guides(shape=FALSE, linetype=FALSE)+
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  strip.text.x = element_text(size = 15),
  strip.text.y = element_text(size = 15),
  #panel.spacing = unit(2, "lines"),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20),
  plot.title = element_text(size=20),
  legend.position="bottom")

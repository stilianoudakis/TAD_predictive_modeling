#19_creating_table_of_model_performances

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

#Metric  Performance
#1                TN 3.291200e+04
#2                FN 2.260000e+02
#3                FP 1.474300e+04
#4                TP 1.369000e+03
#5             Total 4.925000e+04
#6       Sensitivity 8.583072e-01
#7       Specificity 6.906306e-01
#8             Kappa 1.016840e-01
#9          Accuracy 6.960609e-01
#10 BalancedAccuracy 7.744689e-01
#11        Precision 8.496773e-02
#12              FPR 3.093694e-01
#13              FNR 6.819965e-03
#14              FOR 6.819965e-03
#15              NPV 9.931800e-01
#16              MCC 2.071189e-01
#17               F1 1.546281e-01
#18              AUC 8.430922e-01
#19           Youden 5.489378e-01
#20            AUPRC 1.358211e-01



#arrowhead

cl <- c("GM12878","K562")
resolution <- c("5kb","10kb","25kb","50kb","100kb")
resampling <- c("none","ros","rus","smote")
predictor <- c("signal","oc","op","distance")
chromosome <- paste0("CHR",c(1:2))
#metric <- c(9,18,11,17,20)

acc_full <- numeric()
auc_full <- numeric()
prec_full <- numeric()
f1_full <- numeric()
auprc_full <- numeric()

for(j in 1:length(cl)){
  for(k in 1:length(resolution)){
    for(l in 1:length(resampling)){
      for(m in 1:length(predictor)){
        acc <- numeric()
        auc <- numeric()
        prec <- numeric()
        f1 <- numeric()
        auprc <- numeric()
        for(n in 1:length(chromosome)){
          #for(o in metric){
            #"Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/CHR1/rus/distance/TADRF_holdout.rds"
          
            print(paste0("Z:/TAD_data_analysis/",
                         cl[j],
                         "/",
                         resolution[k],
                         "/results_by_chr/",
                         chromosome[n], 
                         "/",
                         resampling[l], 
                         "/", 
                         predictor[m], 
                         "/TADRF_holdout.rds"))
            
            perfdata <- readRDS(paste0("Z:/TAD_data_analysis/",
                                       cl[j],
                                       "/",
                                       resolution[k],
                                       "/results_by_chr/",
                                       chromosome[n], 
                                       "/",
                                       resampling[l], 
                                       "/", 
                                       predictor[m], 
                                       "/TADRF_holdout.rds"))
            
            
            acc[n] <- perfdata[[3]][9,2]
            auc[n] <- perfdata[[3]][18,2]
            prec[n] <- perfdata[[3]][11,2]
            f1[n] <- perfdata[[3]][17,2]
            auprc[n] <- perfdata[[3]][20,2]
          #}
        }
        acc_full <- rbind(acc_full, mean(acc))
        auc_full <- rbind(auc_full, mean(auc))
        prec_full <- rbind(prec_full, mean(prec))
        f1_full <- rbind(f1_full, mean(f1))
        auprc_full <- rbind(auprc_full, mean(auprc))
      }
    }
  }
}

perfdataframe_arrowhead <- data.frame(Tool="Arrowhead",
                            CellLine=c(rep("GM12878",5*4*4),
                                       rep("K562",5*4*4)),
                            Resolution=c(rep("5kb",4*4),
                                         rep("10kb",4*4),
                                         rep("25kb",4*4),
                                         rep("50kb",4*4),
                                         rep("100kb",4*4)),
                            Resampling=c(rep("None",4),
                                         rep("ROS",4),
                                         rep("RUS",4),
                                         rep("SMOTE",4)),
                            Predictor=c("Signal","OC","OP","Distance"),
                            Accuracy=acc_full,
                            AUC=auc_full,
                            Precision=prec_full,
                            F1Score=f1_full,
                            AUPRC=auprc_full
)

#peakachu

cl <- c("GM12878","K562")
resolution <- c("10kb")
resampling <- c("none","ros","rus","smote")
predictor <- c("signal","oc","op","distance")
chromosome <- paste0("CHR",c(1:2))
#metric <- c(9,18,11,17,20)

acc_full <- numeric()
auc_full <- numeric()
prec_full <- numeric()
f1_full <- numeric()
auprc_full <- numeric()

for(j in 1:length(cl)){
  for(k in 1:length(resolution)){
    for(l in 1:length(resampling)){
      for(m in 1:length(predictor)){
        acc <- numeric()
        auc <- numeric()
        prec <- numeric()
        f1 <- numeric()
        auprc <- numeric()
        for(n in 1:length(chromosome)){
          #for(o in metric){
          #"Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/CHR1/rus/distance/TADRF_holdout.rds"
          
          print(paste0("Z:/TAD_data_analysis/",
                       cl[j],
                       "/",
                       resolution[k],
                       "/results_by_chr/",
                       chromosome[n], 
                       "/",
                       resampling[l], 
                       "/", 
                       predictor[m], 
                       "/TADRF_holdout_peakachu.rds"))
          
          perfdata <- readRDS(paste0("Z:/TAD_data_analysis/",
                                     cl[j],
                                     "/",
                                     resolution[k],
                                     "/results_by_chr/",
                                     chromosome[n], 
                                     "/",
                                     resampling[l], 
                                     "/", 
                                     predictor[m], 
                                     "/TADRF_holdout_peakachu.rds"))
          
          
          acc[n] <- perfdata[[3]][9,2]
          auc[n] <- perfdata[[3]][18,2]
          prec[n] <- perfdata[[3]][11,2]
          f1[n] <- perfdata[[3]][17,2]
          auprc[n] <- perfdata[[3]][20,2]
          #}
        }
        acc_full <- rbind(acc_full, mean(acc))
        auc_full <- rbind(auc_full, mean(auc))
        prec_full <- rbind(prec_full, mean(prec))
        f1_full <- rbind(f1_full, mean(f1))
        auprc_full <- rbind(auprc_full, mean(auprc))
      }
    }
  }
}


perfdataframe_peakachu <- data.frame(Tool="Peakachu",
                                      CellLine=c(rep("GM12878",1*4*4),
                                                 rep("K562",1*4*4)),
                                      Resolution=c(rep("10kb",4*4)),
                                      Resampling=c(rep("None",4),
                                                   rep("ROS",4),
                                                   rep("RUS",4),
                                                   rep("SMOTE",4)),
                                      Predictor=c("Signal","OC","OP","Distance"),
                                      Accuracy=acc_full,
                                      AUC=auc_full,
                                      Precision=prec_full,
                                      F1Score=f1_full,
                                      AUPRC=auprc_full
)


perfdataframe_full <- rbind.data.frame(perfdataframe_arrowhead,perfdataframe_peakachu)

perfdataframe_full <- readRDS("Z:/TAD_data_analysis/miscellaneous/perfdataframe_full.rds")

write.csv(perfdataframe_full, "C:/Users/stili/Downloads/perfdataframe_full.csv", quote=FALSE, row.names = FALSE)
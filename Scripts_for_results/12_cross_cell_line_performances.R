#12_cross_cell_line_performances

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
library(ComplexHeatmap)
library(circlize)
library(cotools)
library(ggcorrplot)
library(circlize)
library(egg)

source("Z:/TAD_data_analysis/functions_for_R_package/preciseTAD.R")
source("Z:/TAD_data_analysis/functions_for_R_package/createTADdata.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADRF.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADrfe.R")
source("Z:/TAD_data_analysis/functions_for_R_package/distance_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/percent_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/count_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/binary_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/signal_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/annots_to_granges_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/extract_boundaries_func.R")

genomicElements.GR_gm12878 <- annots_to_granges_func(filepath = "/home/stilianoudakisc/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS_nonames",
                                                     pattern="*.bed",
                                                     signal=4)
genomicElements.GR_k562 <- annots_to_granges_func(filepath = "/home/stilianoudakisc/TAD_data_analysis/annotations/all_common_annotations/k562/topTFBS_nonames",
                                                  pattern="*.bed",
                                                  signal=4)

domains_arrowhead_gm12878 <- read.table("/home/stilianoudakisc/TAD_data_analysis/GSE63525_data/arrowhead_output/GM12878/GM12878_domain_data_5000.b.txt", header=F)
bounds.GR_arrowhead_gm12878 <- extract_boundaries_func(domains.mat=domains_arrowhead_gm12878, 
                                                       preprocess=FALSE, 
                                                       CHR=paste0("CHR", c(1:8,10:22)), 
                                                       resolution=5000)

domains_arrowhead_k562 <- read.table("/home/stilianoudakisc/TAD_data_analysis/GSE63525_data/arrowhead_output/K562/K562_domain_data_5000.b.txt", header=F)
bounds.GR_arrowhead_k562 <- extract_boundaries_func(domains.mat=domains_arrowhead_k562, 
                                                    preprocess=FALSE, 
                                                    CHR=paste0("CHR", c(1:8,10:22)), 
                                                    resolution=5000)

domains_peakachu_gm12878 <- readRDS("/home/stilianoudakisc/TAD_data_analysis/peakachu/peakachu_GM12878.rds")
bounds.GR_peakachu_gm12878 <- extract_boundaries_func(domains.mat=domains_peakachu_gm12878, 
                                                      preprocess=FALSE, 
                                                      CHR=paste0("CHR", c(1:8,10:22)), 
                                                      resolution=5000)

domains_peakachu_k562 <- readRDS("/home/stilianoudakisc/TAD_data_analysis/peakachu/peakachu_K562.rds")
bounds.GR_peakachu_k562 <- extract_boundaries_func(domains.mat=domains_peakachu_k562, 
                                                   preprocess=FALSE, 
                                                   CHR=paste0("CHR", c(1:8,10:22)), 
                                                   resolution=5000)



#arrowhead

##gm12878 on gm12878

arrowhead_g_on_g <- list()
arrowhead_g_on_g_auc <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  TADRF_arrowhead_g_on_g <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                           chrs[i],
                                           "/rus/distance/TADRF_holdout_gm12878_on_gm12878.rds"))
  
  arrowhead_g_on_g_auc[i] <- TADRF_arrowhead_g_on_g[[3]][18,2]
  
  tadData_arrowhead_gm12878 <- createTADdata(bounds.GR=bounds.GR_arrowhead_gm12878,
                                             resolution=5000,
                                             genomicElements.GR=genomicElements.GR_gm12878,
                                             featureType="distance",
                                             resampling="none",
                                             trainCHR=chrs[i],
                                             predictCHR=NULL,
                                             seed=123)
  
  pred_g_on_g <- as.vector(predict(TADRF_arrowhead_g_on_g[[1]], 
                                   newdata=tadData_arrowhead_gm12878[[1]][,-1], 
                                   type="prob")[,"Yes"])
  ROC <- pROC::roc(tadData_arrowhead_gm12878[[1]]$y, pred_g_on_g, quiet = TRUE)
  arrowhead_g_on_g[[i]] <- data.frame(specificity=ROC$specificities,sensitivity=ROC$sensitivities,chr=chrs[i])
  
}
saveRDS(arrowhead_g_on_g, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_g_on_g.rds")
saveRDS(arrowhead_g_on_g_auc, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_g_on_g_auc.rds")

##gm12878 on k562

arrowhead_g_on_k <- list()
arrowhead_g_on_k_auc <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  TADRF_arrowhead_g_on_k <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                           chrs[i],
                                           "/rus/distance/TADRF_holdout_gm12878_on_k562.rds"))
  
  arrowhead_g_on_k_auc[i] <- TADRF_arrowhead_g_on_k[[3]][18,2]
  
  tadData_arrowhead_k562 <- createTADdata(bounds.GR=bounds.GR_arrowhead_k562,
                                          resolution=5000,
                                          genomicElements.GR=genomicElements.GR_k562,
                                          featureType="distance",
                                          resampling="none",
                                          trainCHR=chrs[i],
                                          predictCHR=NULL,
                                          seed=123)
  
  pred_g_on_k <- as.vector(predict(TADRF_arrowhead_g_on_k[[1]], 
                                   newdata=tadData_arrowhead_k562[[1]][,-1], 
                                   type="prob")[,"Yes"])
  ROC <- pROC::roc(tadData_arrowhead_k562[[1]]$y, pred_g_on_k, quiet = TRUE)
  arrowhead_g_on_k[[i]] <- data.frame(specificity=ROC$specificities,sensitivity=ROC$sensitivities,chr=chrs[i])
  
}
saveRDS(arrowhead_g_on_k, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_g_on_k.rds")
saveRDS(arrowhead_g_on_k_auc, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_g_on_k_auc.rds")

#peakachu

##gm12878 on gm12878

peakachu_g_on_g <- list()
peakachu_g_on_g_auc <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  TADRF_peakachu_g_on_g <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                          chrs[i],
                                          "/rus/distance/TADRF_holdout_gm12878_on_gm12878_peakachu.rds"))
  
  peakachu_g_on_g_auc[i] <- TADRF_peakachu_g_on_g[[3]][18,2]
  
  tadData_peakachu_gm12878 <- createTADdata(bounds.GR=bounds.GR_peakachu_gm12878,
                                            resolution=5000,
                                            genomicElements.GR=genomicElements.GR_gm12878,
                                            featureType="distance",
                                            resampling="none",
                                            trainCHR=chrs[i],
                                            predictCHR=NULL,
                                            seed=123)
  
  pred_g_on_g <- as.vector(predict(TADRF_peakachu_g_on_g[[1]], 
                                   newdata=tadData_peakachu_gm12878[[1]][,-1], 
                                   type="prob")[,"Yes"])
  ROC <- pROC::roc(tadData_peakachu_gm12878[[1]]$y, pred_g_on_g, quiet = TRUE)
  peakachu_g_on_g[[i]] <- data.frame(specificity=ROC$specificities,sensitivity=ROC$sensitivities,chr=chrs[i])
  
}
saveRDS(peakachu_g_on_g, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_g_on_g.rds")
saveRDS(peakachu_g_on_g_auc, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_g_on_g_auc.rds")

##gm12878 on k562

peakachu_g_on_k <- list()
peakachu_g_on_k_auc <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  TADRF_peakachu_g_on_k <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                          chrs[i],
                                          "/rus/distance/TADRF_holdout_gm12878_on_k562_peakachu.rds"))
  
  peakachu_g_on_k_auc[i] <- TADRF_peakachu_g_on_k[[3]][18,2]
  
  tadData_peakachu_k562 <- createTADdata(bounds.GR=bounds.GR_peakachu_k562,
                                         resolution=5000,
                                         genomicElements.GR=genomicElements.GR_k562,
                                         featureType="distance",
                                         resampling="none",
                                         trainCHR=chrs[i],
                                         predictCHR=NULL,
                                         seed=123)
  
  pred_g_on_k <- as.vector(predict(TADRF_peakachu_g_on_k[[1]], 
                                   newdata=tadData_peakachu_k562[[1]][,-1], 
                                   type="prob")[,"Yes"])
  ROC <- pROC::roc(tadData_peakachu_k562[[1]]$y, pred_g_on_k, quiet = TRUE)
  peakachu_g_on_k[[i]] <- data.frame(specificity=ROC$specificities,sensitivity=ROC$sensitivities,chr=chrs[i])
  
}
saveRDS(peakachu_g_on_k, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_g_on_k.rds")
saveRDS(peakachu_g_on_k_auc, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_g_on_k_auc.rds")

###################################################################################################

#arrowhead

##k562 on k562

arrowhead_k_on_k <- list()
arrowhead_k_on_k_auc <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  TADRF_arrowhead_k_on_k <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/5kb/results_by_chr/",
                                           chrs[i],
                                           "/rus/distance/TADRF_holdout_k562_on_k562.rds"))
  
  arrowhead_k_on_k_auc[i] <- TADRF_arrowhead_k_on_k[[3]][18,2]
  
  tadData_arrowhead_k562 <- createTADdata(bounds.GR=bounds.GR_arrowhead_k562,
                                          resolution=5000,
                                          genomicElements.GR=genomicElements.GR_k562,
                                          featureType="distance",
                                          resampling="none",
                                          trainCHR=chrs[i],
                                          predictCHR=NULL,
                                          seed=123)
  
  pred_k_on_k <- as.vector(predict(TADRF_arrowhead_k_on_k[[1]], 
                                   newdata=tadData_arrowhead_k562[[1]][,-1], 
                                   type="prob")[,"Yes"])
  ROC <- pROC::roc(tadData_arrowhead_k562[[1]]$y, pred_k_on_k, quiet = TRUE)
  arrowhead_k_on_k[[i]] <- data.frame(specificity=ROC$specificities,sensitivity=ROC$sensitivities,chr=chrs[i])
  
}
saveRDS(arrowhead_k_on_k, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_k_on_k.rds")
saveRDS(arrowhead_k_on_k_auc, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_k_on_k_auc.rds")

##k562 on gm12878

arrowhead_k_on_g <- list()
arrowhead_k_on_g_auc <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  TADRF_arrowhead_k_on_g <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/5kb/results_by_chr/",
                                           chrs[i],
                                           "/rus/distance/TADRF_holdout_k562_on_gm12878.rds"))
  
  arrowhead_k_on_g_auc[i] <- TADRF_arrowhead_k_on_g[[3]][18,2]
  
  tadData_arrowhead_gm12878 <- createTADdata(bounds.GR=bounds.GR_arrowhead_gm12878,
                                             resolution=5000,
                                             genomicElements.GR=genomicElements.GR_gm12878,
                                             featureType="distance",
                                             resampling="none",
                                             trainCHR=chrs[i],
                                             predictCHR=NULL,
                                             seed=123)
  
  pred_k_on_g <- as.vector(predict(TADRF_arrowhead_k_on_g[[1]], 
                                   newdata=tadData_arrowhead_gm12878[[1]][,-1], 
                                   type="prob")[,"Yes"])
  ROC <- pROC::roc(tadData_arrowhead_gm12878[[1]]$y, pred_k_on_g, quiet = TRUE)
  arrowhead_k_on_g[[i]] <- data.frame(specificity=ROC$specificities,sensitivity=ROC$sensitivities,chr=chrs[i])
  
}
saveRDS(arrowhead_k_on_g, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_k_on_g.rds")
saveRDS(arrowhead_k_on_g_auc, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_k_on_g_auc.rds")

#peakachu

##k562 on k562

peakachu_k_on_k <- list()
peakachu_k_on_k_auc <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  TADRF_peakachu_k_on_k <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/10kb/results_by_chr/",
                                          chrs[i],
                                          "/rus/distance/TADRF_holdout_k562_on_k562_peakachu.rds"))
  
  peakachu_k_on_k_auc[i] <- TADRF_peakachu_k_on_k[[3]][18,2]
  
  tadData_peakachu_k562 <- createTADdata(bounds.GR=bounds.GR_peakachu_k562,
                                         resolution=5000,
                                         genomicElements.GR=genomicElements.GR_k562,
                                         featureType="distance",
                                         resampling="none",
                                         trainCHR=chrs[i],
                                         predictCHR=NULL,
                                         seed=123)
  
  pred_k_on_k <- as.vector(predict(TADRF_peakachu_k_on_k[[1]], 
                                   newdata=tadData_peakachu_k562[[1]][,-1], 
                                   type="prob")[,"Yes"])
  ROC <- pROC::roc(tadData_peakachu_k562[[1]]$y, pred_k_on_k, quiet = TRUE)
  peakachu_k_on_k[[i]] <- data.frame(specificity=ROC$specificities,sensitivity=ROC$sensitivities,chr=chrs[i])
  
}
saveRDS(peakachu_k_on_k, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_k_on_k.rds")
saveRDS(peakachu_k_on_k_auc, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_k_on_k_auc.rds")

##k562 on gm12878

peakachu_k_on_g <- list()
peakachu_k_on_g_auc <- numeric()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  print(chrs[i])
  TADRF_peakachu_k_on_g <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/10kb/results_by_chr/",
                                          chrs[i],
                                          "/rus/distance/TADRF_holdout_k562_on_gm12878_peakachu.rds"))
  
  peakachu_k_on_g_auc[i] <- TADRF_peakachu_k_on_g[[3]][18,2]
  
  tadData_peakachu_gm12878 <- createTADdata(bounds.GR=bounds.GR_peakachu_gm12878,
                                            resolution=5000,
                                            genomicElements.GR=genomicElements.GR_gm12878,
                                            featureType="distance",
                                            resampling="none",
                                            trainCHR=chrs[i],
                                            predictCHR=NULL,
                                            seed=123)
  
  pred_k_on_g <- as.vector(predict(TADRF_peakachu_k_on_g[[1]], 
                                   newdata=tadData_peakachu_gm12878[[1]][,-1], 
                                   type="prob")[,"Yes"])
  ROC <- pROC::roc(tadData_peakachu_gm12878[[1]]$y, pred_k_on_g, quiet = TRUE)
  peakachu_k_on_g[[i]] <- data.frame(specificity=ROC$specificities,sensitivity=ROC$sensitivities,chr=chrs[i])
  
}
saveRDS(peakachu_k_on_g, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_k_on_g.rds")
saveRDS(peakachu_k_on_g_auc, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_k_on_g_auc.rds")

#######################################################################################################################

##arrowhead

###g on g

arrowhead_g_on_g <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_g_on_g.rds")

for(i in 1:length(arrowhead_g_on_g)){
  m <- max(unlist(lapply(arrowhead_g_on_g, nrow)))
  if(nrow(arrowhead_g_on_g[[i]]) != m){
    d <- m - nrow(arrowhead_g_on_g[[i]])
    arrowhead_g_on_g[[i]] <- rbind.data.frame(data.frame(specificity=rep(0,d),
                                                         sensitivity=rep(1,d),
                                                         chr=rep(arrowhead_g_on_g[[i]]$chr[1],d)),
                                              arrowhead_g_on_g[[i]])
  }
}

arrowhead_g_on_g <- do.call("rbind.data.frame",arrowhead_g_on_g)
meanroc_g_on_g <- arrowhead_g_on_g %>% group_by(rep(c(1:502),21)) %>% 
  summarize(meanSens = mean(sensitivity), meanSpec = mean(specificity))

meanroc_g_on_g$sdSensunder <- arrowhead_g_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>% 
  mutate(sdSens=meanroc_g_on_g$meanSens-sdSens) %>% 
  select(sdSens)

meanroc_g_on_g$sdSensover <- arrowhead_g_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>%
  mutate(sdSens=meanroc_g_on_g$meanSens+sdSens) %>% 
  select(sdSens)

meanroc_g_on_g$sdSpecunder <- arrowhead_g_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_g_on_g$meanSpec-sdSpec) %>% 
  select(sdSpec)

meanroc_g_on_g$sdSpecover <- arrowhead_g_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_g_on_g$meanSpec+sdSpec) %>% 
  select(sdSpec)


###k on g

arrowhead_k_on_g <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_k_on_g.rds")

for(i in 1:length(arrowhead_k_on_g)){
  m <- max(unlist(lapply(arrowhead_k_on_g, nrow)))
  if(nrow(arrowhead_k_on_g[[i]]) != m){
    d <- m - nrow(arrowhead_k_on_g[[i]])
    arrowhead_k_on_g[[i]] <- rbind.data.frame(data.frame(specificity=rep(0,d),
                                                         sensitivity=rep(1,d),
                                                         chr=rep(arrowhead_k_on_g[[i]]$chr[1],d)),
                                              arrowhead_k_on_g[[i]])
  }
}

arrowhead_k_on_g <- do.call("rbind.data.frame",arrowhead_k_on_g)
meanroc_k_on_g <- arrowhead_k_on_g %>% group_by(rep(c(1:502),21)) %>% 
  summarize(meanSens = mean(sensitivity), meanSpec = mean(specificity))

meanroc_k_on_g$sdSensunder <- arrowhead_k_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>% 
  mutate(sdSens=meanroc_k_on_g$meanSens-sdSens) %>% 
  select(sdSens)

meanroc_k_on_g$sdSensover <- arrowhead_k_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>%
  mutate(sdSens=meanroc_k_on_g$meanSens+sdSens) %>% 
  select(sdSens)

meanroc_k_on_g$sdSpecunder <- arrowhead_k_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_k_on_g$meanSpec-sdSpec) %>% 
  select(sdSpec)

meanroc_k_on_g$sdSpecover <- arrowhead_k_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_k_on_g$meanSpec+sdSpec) %>% 
  select(sdSpec)


###plotting

a1 <- ggplot() + 
  #g on g
  ##under
  #geom_line(data=meanroc_g_on_g,aes(x=1-sdSpecunder$sdSpec, y=sdSensunder$sdSens),
  #          size=.5,
  #          #linetype="solid", 
  #          color="blue") +
  ##over
  #geom_line(data=meanroc_g_on_g,aes(x=1-sdSpecover$sdSpec, y=sdSensover$sdSens),
  #          size=.5,
  #          #linetype="solid", 
  #          color="blue") +
geom_ribbon(data=meanroc_g_on_g, 
            aes(#ymin=sdSensunder$sdSens,
              #ymax=sdSensover$sdSens,
              xmin=1-meanSpec,
              xmax=1-sdSpecover$sdSpec,
              x=1-sdSpecunder$sdSpec,
              y=sdSensover$sdSens), fill="blue", alpha=0.2) +
  geom_ribbon(data=meanroc_g_on_g, 
              aes(#ymin=sdSensunder$sdSens,
                #ymax=sdSensover$sdSens,
                xmin=1-sdSpecunder$sdSpec,
                xmax=1-meanSpec,
                x=1-sdSpecunder$sdSpec,
                y=sdSensover$sdSens), fill="blue", alpha=0.2) +
  geom_ribbon(data=meanroc_g_on_g, 
              aes(ymin=sdSensunder$sdSens,
                  ymax=sdSensover$sdSens,
                  #xmin=1-meanSpec,
                  #xmax=1-sdSpecunder$sdSpec,
                  x=1-sdSpecunder$sdSpec,
                  y=sdSensover$sdSens), fill="blue", alpha=0.2) +
  ##mean
  geom_line(data=meanroc_g_on_g,aes(x=1-meanSpec, y=meanSens),
            size=1.5,
            linetype="solid", 
            color="blue") +
  #k on g
  ##under
  #geom_line(data=meanroc_k_on_g,aes(x=1-sdSpecunder$sdSpec, y=sdSensunder$sdSens),
  #          size=.5,
  #          #linetype="dashed", 
  #          color="blue") +
  ##over
  #geom_line(data=meanroc_k_on_g,aes(x=1-sdSpecover$sdSpec, y=sdSensover$sdSens),
  #          size=.5,
  #          #linetype="dashed", 
  #          color="blue") +
geom_ribbon(data=meanroc_k_on_g, 
            aes(#ymin=sdSensunder$sdSens,
              #ymax=sdSensover$sdSens,
              xmin=1-meanSpec,
              xmax=1-sdSpecover$sdSpec,
              x=1-sdSpecunder$sdSpec,
              y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  geom_ribbon(data=meanroc_k_on_g, 
              aes(#ymin=sdSensunder$sdSens,
                #ymax=sdSensover$sdSens,
                xmin=1-sdSpecunder$sdSpec,
                xmax=1-meanSpec,
                x=1-sdSpecunder$sdSpec,
                y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  geom_ribbon(data=meanroc_k_on_g, 
              aes(ymin=sdSensunder$sdSens,
                  ymax=sdSensover$sdSens,
                  #xmin=1-meanSpec,
                  #xmax=1-sdSpecunder$sdSpec,
                  x=1-sdSpecunder$sdSpec,
                  y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  ##mean
  geom_line(data=meanroc_k_on_g,aes(x=1-meanSpec, y=meanSens),
            size=1.5,
            linetype="dashed", 
            color="black") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  scale_color_manual(name = "Training/Testing",
                     labels = c("GM on GM",
                                "K on GM"),
                     values = c("blue","black")) +
  scale_x_continuous(breaks=c(0, 0.5, 1),
                     labels=c(0, 0.5, 1),
                     limits=c(0,1))+
  scale_y_continuous(breaks=c(0, 0.5, 1),
                     labels=c(0, 0.5, 1),
                     limits=c(0,1))+
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


##peakachu

###g on g

peakachu_g_on_g <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_g_on_g.rds")

for(i in 1:length(peakachu_g_on_g)){
  m <- max(unlist(lapply(peakachu_g_on_g, nrow)))
  if(nrow(peakachu_g_on_g[[i]]) != m){
    d <- m - nrow(peakachu_g_on_g[[i]])
    peakachu_g_on_g[[i]] <- rbind.data.frame(data.frame(specificity=rep(0,d),
                                                        sensitivity=rep(1,d),
                                                        chr=rep(peakachu_g_on_g[[i]]$chr[1],d)),
                                             peakachu_g_on_g[[i]])
  }
}

peakachu_g_on_g <- do.call("rbind.data.frame",peakachu_g_on_g)
meanroc_g_on_g <- peakachu_g_on_g %>% group_by(rep(c(1:502),21)) %>% 
  summarize(meanSens = mean(sensitivity), meanSpec = mean(specificity))

meanroc_g_on_g$sdSensunder <- peakachu_g_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>% 
  mutate(sdSens=meanroc_g_on_g$meanSens-sdSens) %>% 
  select(sdSens)

meanroc_g_on_g$sdSensover <- peakachu_g_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>%
  mutate(sdSens=meanroc_g_on_g$meanSens+sdSens) %>% 
  select(sdSens)

meanroc_g_on_g$sdSpecunder <- peakachu_g_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_g_on_g$meanSpec-sdSpec) %>% 
  select(sdSpec)

meanroc_g_on_g$sdSpecover <- peakachu_g_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_g_on_g$meanSpec+sdSpec) %>% 
  select(sdSpec)


###k on g

peakachu_k_on_g <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_k_on_g.rds")

for(i in 1:length(peakachu_k_on_g)){
  m <- max(unlist(lapply(peakachu_k_on_g, nrow)))
  if(nrow(peakachu_k_on_g[[i]]) != m){
    d <- m - nrow(peakachu_k_on_g[[i]])
    peakachu_k_on_g[[i]] <- rbind.data.frame(data.frame(specificity=rep(0,d),
                                                        sensitivity=rep(1,d),
                                                        chr=rep(peakachu_k_on_g[[i]]$chr[1],d)),
                                             peakachu_k_on_g[[i]])
  }
}

peakachu_k_on_g <- do.call("rbind.data.frame",peakachu_k_on_g)
meanroc_k_on_g <- peakachu_k_on_g %>% group_by(rep(c(1:502),21)) %>% 
  summarize(meanSens = mean(sensitivity), meanSpec = mean(specificity))

meanroc_k_on_g$sdSensunder <- peakachu_k_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>% 
  mutate(sdSens=meanroc_k_on_g$meanSens-sdSens) %>% 
  select(sdSens)

meanroc_k_on_g$sdSensover <- peakachu_k_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>%
  mutate(sdSens=meanroc_k_on_g$meanSens+sdSens) %>% 
  select(sdSens)

meanroc_k_on_g$sdSpecunder <- peakachu_k_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_k_on_g$meanSpec-sdSpec) %>% 
  select(sdSpec)

meanroc_k_on_g$sdSpecover <- peakachu_k_on_g %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_k_on_g$meanSpec+sdSpec) %>% 
  select(sdSpec)


###plotting

b1 <- ggplot() + 
  #g on g
  ##under
  #geom_line(data=meanroc_g_on_g,aes(x=1-sdSpecunder$sdSpec, y=sdSensunder$sdSens),
  #          size=.5,
  #          #linetype="solid", 
  #          color="red") +
  ##over
  #geom_line(data=meanroc_g_on_g,aes(x=1-sdSpecover$sdSpec, y=sdSensover$sdSens),
  #          size=.5,
  #          #linetype="solid", 
  #          color="red") +
geom_ribbon(data=meanroc_g_on_g, 
            aes(#ymin=sdSensunder$sdSens,
              #ymax=sdSensover$sdSens,
              xmin=1-meanSpec,
              xmax=1-sdSpecover$sdSpec,
              x=1-sdSpecunder$sdSpec,
              y=sdSensover$sdSens), fill="red", alpha=0.2) +
  geom_ribbon(data=meanroc_g_on_g, 
              aes(#ymin=sdSensunder$sdSens,
                #ymax=sdSensover$sdSens,
                xmin=1-sdSpecunder$sdSpec,
                xmax=1-meanSpec,
                x=1-sdSpecunder$sdSpec,
                y=sdSensover$sdSens), fill="red", alpha=0.2) +
  geom_ribbon(data=meanroc_g_on_g, 
              aes(ymin=sdSensunder$sdSens,
                  ymax=sdSensover$sdSens,
                  #xmin=1-meanSpec,
                  #xmax=1-sdSpecunder$sdSpec,
                  x=1-sdSpecunder$sdSpec,
                  y=sdSensover$sdSens), fill="red", alpha=0.2) +
  ##mean
  geom_line(data=meanroc_g_on_g,aes(x=1-meanSpec, y=meanSens),
            size=1.5,
            linetype="solid", 
            color="red") +
  #k on g
  ##under
  #geom_line(data=meanroc_k_on_g,aes(x=1-sdSpecunder$sdSpec, y=sdSensunder$sdSens),
  #          size=.5,
  #          #linetype="dashed", 
  #          color="red") +
  ##over
  #geom_line(data=meanroc_k_on_g,aes(x=1-sdSpecover$sdSpec, y=sdSensover$sdSens),
  #          size=.5,
  #          #linetype="dashed", 
  #          color="red") +
geom_ribbon(data=meanroc_k_on_g, 
            aes(#ymin=sdSensunder$sdSens,
              #ymax=sdSensover$sdSens,
              xmin=1-meanSpec,
              xmax=1-sdSpecover$sdSpec,
              x=1-sdSpecunder$sdSpec,
              y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  geom_ribbon(data=meanroc_k_on_g, 
              aes(#ymin=sdSensunder$sdSens,
                #ymax=sdSensover$sdSens,
                xmin=1-sdSpecunder$sdSpec,
                xmax=1-meanSpec,
                x=1-sdSpecunder$sdSpec,
                y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  geom_ribbon(data=meanroc_k_on_g, 
              aes(ymin=sdSensunder$sdSens,
                  ymax=sdSensover$sdSens,
                  #xmin=1-meanSpec,
                  #xmax=1-sdSpecunder$sdSpec,
                  x=1-sdSpecunder$sdSpec,
                  y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  ##mean
  geom_line(data=meanroc_k_on_g,aes(x=1-meanSpec, y=meanSens),
            size=1.5,
            linetype="dashed", 
            color="black") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  scale_color_manual(name = "Training/Testing",
                     labels = c("GM on GM",
                                "K on GM"),
                     values = c("red","black")) +
  scale_x_continuous(breaks=c(0, 0.5, 1),
                     labels=c(0, 0.5, 1),
                     limits=c(0,1))+
  scale_y_continuous(breaks=c(0, 0.5, 1),
                     labels=c(0, 0.5, 1),
                     limits=c(0,1))+
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


f <- ggarrange(a1, b1, ncol = 2)

annotate_figure(f,
                bottom = text_grob("1-Specificity", 
                                   color = "black",
                                   size = 20),
                left = text_grob("Sensitivity", 
                                 color = "black",
                                 size = 20,
                                 rot = 90))

#########################################################################################

##arrowhead

###k on k

arrowhead_k_on_k <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_k_on_k.rds")

for(i in 1:length(arrowhead_k_on_k)){
  m <- max(unlist(lapply(arrowhead_k_on_k, nrow)))
  if(nrow(arrowhead_k_on_k[[i]]) != m){
    d <- m - nrow(arrowhead_k_on_k[[i]])
    arrowhead_k_on_k[[i]] <- rbind.data.frame(data.frame(specificity=rep(0,d),
                                                         sensitivity=rep(1,d),
                                                         chr=rep(arrowhead_k_on_k[[i]]$chr[1],d)),
                                              arrowhead_k_on_k[[i]])
  }
}

arrowhead_k_on_k <- do.call("rbind.data.frame",arrowhead_k_on_k)
meanroc_k_on_k <- arrowhead_k_on_k %>% group_by(rep(c(1:502),21)) %>% 
  summarize(meanSens = mean(sensitivity), meanSpec = mean(specificity))

meanroc_k_on_k$sdSensunder <- arrowhead_k_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>% 
  mutate(sdSens=meanroc_k_on_k$meanSens-sdSens) %>% 
  select(sdSens)

meanroc_k_on_k$sdSensover <- arrowhead_k_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>%
  mutate(sdSens=meanroc_k_on_k$meanSens+sdSens) %>% 
  select(sdSens)

meanroc_k_on_k$sdSpecunder <- arrowhead_k_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_k_on_k$meanSpec-sdSpec) %>% 
  select(sdSpec)

meanroc_k_on_k$sdSpecover <- arrowhead_k_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_k_on_k$meanSpec+sdSpec) %>% 
  select(sdSpec)


###g on k

arrowhead_g_on_k <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_g_on_k.rds")

for(i in 1:length(arrowhead_g_on_k)){
  m <- max(unlist(lapply(arrowhead_g_on_k, nrow)))
  if(nrow(arrowhead_g_on_k[[i]]) != m){
    d <- m - nrow(arrowhead_g_on_k[[i]])
    arrowhead_g_on_k[[i]] <- rbind.data.frame(data.frame(specificity=rep(0,d),
                                                         sensitivity=rep(1,d),
                                                         chr=rep(arrowhead_g_on_k[[i]]$chr[1],d)),
                                              arrowhead_g_on_k[[i]])
  }
}

arrowhead_g_on_k <- do.call("rbind.data.frame",arrowhead_g_on_k)
meanroc_g_on_k <- arrowhead_g_on_k %>% group_by(rep(c(1:502),21)) %>% 
  summarize(meanSens = mean(sensitivity), meanSpec = mean(specificity))

meanroc_g_on_k$sdSensunder <- arrowhead_g_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>% 
  mutate(sdSens=meanroc_g_on_k$meanSens-sdSens) %>% 
  select(sdSens)

meanroc_g_on_k$sdSensover <- arrowhead_g_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>%
  mutate(sdSens=meanroc_g_on_k$meanSens+sdSens) %>% 
  select(sdSens)

meanroc_g_on_k$sdSpecunder <- arrowhead_g_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_g_on_k$meanSpec-sdSpec) %>% 
  select(sdSpec)

meanroc_g_on_k$sdSpecover <- arrowhead_g_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_g_on_k$meanSpec+sdSpec) %>% 
  select(sdSpec)


a2 <- ggplot() + 
  #k on k
  ##under
  #geom_line(data=meanroc_k_on_k,aes(x=1-sdSpecunder$sdSpec, y=sdSensunder$sdSens),
  #          size=.5,
  #          #linetype="solid", 
  #          color="blue") +
  ##over
  #geom_line(data=meanroc_k_on_k,aes(x=1-sdSpecover$sdSpec, y=sdSensover$sdSens),
  #          size=.5,
  #          #linetype="solid", 
  #          color="blue") +
geom_ribbon(data=meanroc_k_on_k, 
            aes(#ymin=sdSensunder$sdSens,
              #ymax=sdSensover$sdSens,
              xmin=1-meanSpec,
              xmax=1-sdSpecover$sdSpec,
              x=1-sdSpecunder$sdSpec,
              y=sdSensover$sdSens), fill="blue", alpha=0.2) +
  geom_ribbon(data=meanroc_k_on_k, 
              aes(#ymin=sdSensunder$sdSens,
                #ymax=sdSensover$sdSens,
                xmin=1-sdSpecunder$sdSpec,
                xmax=1-meanSpec,
                x=1-sdSpecunder$sdSpec,
                y=sdSensover$sdSens), fill="blue", alpha=0.2) +
  geom_ribbon(data=meanroc_k_on_k, 
              aes(ymin=sdSensunder$sdSens,
                  ymax=sdSensover$sdSens,
                  #xmin=1-meanSpec,
                  #xmax=1-sdSpecunder$sdSpec,
                  x=1-sdSpecunder$sdSpec,
                  y=sdSensover$sdSens), fill="blue", alpha=0.2) +
  ##mean
  geom_line(data=meanroc_k_on_k,aes(x=1-meanSpec, y=meanSens),
            size=1.5,
            linetype="solid", 
            color="blue") +
  #g on k
  ##under
  #geom_line(data=meanroc_g_on_k,aes(x=1-sdSpecunder$sdSpec, y=sdSensunder$sdSens),
  #          size=.5,
  #          #linetype="dashed", 
  #          color="blue") +
  ##over
  #geom_line(data=meanroc_g_on_k,aes(x=1-sdSpecover$sdSpec, y=sdSensover$sdSens),
  #          size=.5,
  #          #linetype="dashed", 
  #          color="blue") +
geom_ribbon(data=meanroc_g_on_k, 
            aes(#ymin=sdSensunder$sdSens,
              #ymax=sdSensover$sdSens,
              xmin=1-meanSpec,
              xmax=1-sdSpecover$sdSpec,
              x=1-sdSpecunder$sdSpec,
              y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  geom_ribbon(data=meanroc_g_on_k, 
              aes(#ymin=sdSensunder$sdSens,
                #ymax=sdSensover$sdSens,
                xmin=1-sdSpecunder$sdSpec,
                xmax=1-meanSpec,
                x=1-sdSpecunder$sdSpec,
                y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  geom_ribbon(data=meanroc_g_on_k, 
              aes(ymin=sdSensunder$sdSens,
                  ymax=sdSensover$sdSens,
                  #xmin=1-meanSpec,
                  #xmax=1-sdSpecunder$sdSpec,
                  x=1-sdSpecunder$sdSpec,
                  y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  ##mean
  geom_line(data=meanroc_g_on_k,aes(x=1-meanSpec, y=meanSens),
            size=1.5,
            linetype="dashed", 
            color="black") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  scale_color_manual(name = "Training/Testing",
                     labels = c("K on K",
                                "GM on K"),
                     values = c("blue","black")) +
  scale_x_continuous(breaks=c(0, 0.5, 1),
                     labels=c(0, 0.5, 1),
                     limits=c(0,1))+
  scale_y_continuous(breaks=c(0, 0.5, 1),
                     labels=c(0, 0.5, 1),
                     limits=c(0,1))+
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


##peakachu

###k on k

peakachu_k_on_k <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_k_on_k.rds")

for(i in 1:length(peakachu_k_on_k)){
  m <- max(unlist(lapply(peakachu_k_on_k, nrow)))
  if(nrow(peakachu_k_on_k[[i]]) != m){
    d <- m - nrow(peakachu_k_on_k[[i]])
    peakachu_k_on_k[[i]] <- rbind.data.frame(data.frame(specificity=rep(0,d),
                                                        sensitivity=rep(1,d),
                                                        chr=rep(peakachu_k_on_k[[i]]$chr[1],d)),
                                             peakachu_k_on_k[[i]])
  }
}

peakachu_k_on_k <- do.call("rbind.data.frame",peakachu_k_on_k)
meanroc_k_on_k <- peakachu_k_on_k %>% group_by(rep(c(1:502),21)) %>% 
  summarize(meanSens = mean(sensitivity), meanSpec = mean(specificity))

meanroc_k_on_k$sdSensunder <- peakachu_k_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>% 
  mutate(sdSens=meanroc_k_on_k$meanSens-sdSens) %>% 
  select(sdSens)

meanroc_k_on_k$sdSensover <- peakachu_k_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>%
  mutate(sdSens=meanroc_k_on_k$meanSens+sdSens) %>% 
  select(sdSens)

meanroc_k_on_k$sdSpecunder <- peakachu_k_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_k_on_k$meanSpec-sdSpec) %>% 
  select(sdSpec)

meanroc_k_on_k$sdSpecover <- peakachu_k_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_k_on_k$meanSpec+sdSpec) %>% 
  select(sdSpec)


###g on k

peakachu_g_on_k <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_g_on_k.rds")

for(i in 1:length(peakachu_g_on_k)){
  m <- max(unlist(lapply(peakachu_g_on_k, nrow)))
  if(nrow(peakachu_g_on_k[[i]]) != m){
    d <- m - nrow(peakachu_g_on_k[[i]])
    peakachu_g_on_k[[i]] <- rbind.data.frame(data.frame(specificity=rep(0,d),
                                                        sensitivity=rep(1,d),
                                                        chr=rep(peakachu_g_on_k[[i]]$chr[1],d)),
                                             peakachu_g_on_k[[i]])
  }
}

peakachu_g_on_k <- do.call("rbind.data.frame",peakachu_g_on_k)
meanroc_g_on_k <- peakachu_g_on_k %>% group_by(rep(c(1:502),21)) %>% 
  summarize(meanSens = mean(sensitivity), meanSpec = mean(specificity))

meanroc_g_on_k$sdSensunder <- peakachu_g_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>% 
  mutate(sdSens=meanroc_g_on_k$meanSens-sdSens) %>% 
  select(sdSens)

meanroc_g_on_k$sdSensover <- peakachu_g_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSens = sd(sensitivity)) %>%
  mutate(sdSens=meanroc_g_on_k$meanSens+sdSens) %>% 
  select(sdSens)

meanroc_g_on_k$sdSpecunder <- peakachu_g_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_g_on_k$meanSpec-sdSpec) %>% 
  select(sdSpec)

meanroc_g_on_k$sdSpecover <- peakachu_g_on_k %>% 
  group_by(rep(c(1:502),21)) %>% 
  summarize(sdSpec = sd(specificity)) %>% 
  mutate(sdSpec=meanroc_g_on_k$meanSpec+sdSpec) %>% 
  select(sdSpec)


b2 <- ggplot() + 
  #k on k
  ##under
  #geom_line(data=meanroc_k_on_k,aes(x=1-sdSpecunder$sdSpec, y=sdSensunder$sdSens),
  #          size=.5,
  #          #linetype="solid", 
  #          color="red") +
  ##over
  #geom_line(data=meanroc_k_on_k,aes(x=1-sdSpecover$sdSpec, y=sdSensover$sdSens),
  #          size=.5,
  #          #linetype="solid", 
  #          color="red") +
geom_ribbon(data=meanroc_k_on_k, 
            aes(#ymin=sdSensunder$sdSens,
              #ymax=sdSensover$sdSens,
              xmin=1-meanSpec,
              xmax=1-sdSpecover$sdSpec,
              x=1-sdSpecunder$sdSpec,
              y=sdSensover$sdSens), fill="red", alpha=0.2) +
  geom_ribbon(data=meanroc_k_on_k, 
              aes(#ymin=sdSensunder$sdSens,
                #ymax=sdSensover$sdSens,
                xmin=1-sdSpecunder$sdSpec,
                xmax=1-meanSpec,
                x=1-sdSpecunder$sdSpec,
                y=sdSensover$sdSens), fill="red", alpha=0.2) +
  geom_ribbon(data=meanroc_k_on_k, 
              aes(ymin=sdSensunder$sdSens,
                  ymax=sdSensover$sdSens,
                  #xmin=1-meanSpec,
                  #xmax=1-sdSpecunder$sdSpec,
                  x=1-sdSpecunder$sdSpec,
                  y=sdSensover$sdSens), fill="red", alpha=0.2) +
  ##mean
  geom_line(data=meanroc_k_on_k,aes(x=1-meanSpec, y=meanSens),
            size=1.5,
            linetype="solid", 
            color="red") +
  #g on k
  ##under
  #geom_line(data=meanroc_g_on_k,aes(x=1-sdSpecunder$sdSpec, y=sdSensunder$sdSens),
  #          size=.5,
  #          #linetype="dashed", 
  #          color="red") +
  ##over
  #geom_line(data=meanroc_g_on_k,aes(x=1-sdSpecover$sdSpec, y=sdSensover$sdSens),
  #          size=.5,
  #          #linetype="dashed", 
  #          color="red") +
geom_ribbon(data=meanroc_g_on_k, 
            aes(#ymin=sdSensunder$sdSens,
              #ymax=sdSensover$sdSens,
              xmin=1-meanSpec,
              xmax=1-sdSpecover$sdSpec,
              x=1-sdSpecunder$sdSpec,
              y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  geom_ribbon(data=meanroc_g_on_k, 
              aes(#ymin=sdSensunder$sdSens,
                #ymax=sdSensover$sdSens,
                xmin=1-sdSpecunder$sdSpec,
                xmax=1-meanSpec,
                x=1-sdSpecunder$sdSpec,
                y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  geom_ribbon(data=meanroc_g_on_k, 
              aes(ymin=sdSensunder$sdSens,
                  ymax=sdSensover$sdSens,
                  #xmin=1-meanSpec,
                  #xmax=1-sdSpecunder$sdSpec,
                  x=1-sdSpecunder$sdSpec,
                  y=sdSensover$sdSens), fill="gray", alpha=0.3) +
  ##mean
  geom_line(data=meanroc_g_on_k,aes(x=1-meanSpec, y=meanSens),
            size=1.5,
            linetype="dashed", 
            color="black") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  scale_color_manual(name = "Training/Testing",
                     labels = c("K on K",
                                "GM on K"),
                     values = c("red","black")) +
  scale_x_continuous(breaks=c(0, 0.5, 1),
                     labels=c(0, 0.5, 1),
                     limits=c(0,1))+
  scale_y_continuous(breaks=c(0, 0.5, 1),
                     labels=c(0, 0.5, 1),
                     limits=c(0,1))+
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


f <- ggarrange(a2, b2, ncol = 2)

annotate_figure(f,
                bottom = text_grob("1-Specificity", 
                                   color = "black",
                                   size = 20),
                left = text_grob("Sensitivity", 
                                 color = "black",
                                 size = 20,
                                 rot = 90))

#########################################################################################

arrowhead_g_on_g_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_g_on_g_auc.rds")
mean(arrowhead_g_on_g_auc)
arrowhead_k_on_g_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_k_on_g_auc.rds")
mean(arrowhead_k_on_g_auc)


peakachu_g_on_g_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_g_on_g_auc.rds")
mean(peakachu_g_on_g_auc)
peakachu_k_on_g_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_k_on_g_auc.rds")
mean(peakachu_k_on_g_auc)


arrowhead_k_on_k_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_k_on_k_auc.rds")
mean(arrowhead_k_on_k_auc)
arrowhead_g_on_k_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/arrowhead_g_on_k_auc.rds")
mean(arrowhead_g_on_k_auc)



peakachu_k_on_k_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_k_on_k_auc.rds")
mean(peakachu_k_on_k_auc)
peakachu_g_on_k_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/cross_cell_line_sens_spec/peakachu_g_on_k_auc.rds")
mean(peakachu_g_on_k_auc)


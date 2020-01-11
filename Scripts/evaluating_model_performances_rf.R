library(ggplot2)
library(dplyr)

# GM12878

##5kb

mccdata_gm12878_5kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_gm12878_5kb <- rbind(mccdata_gm12878_5kb,
                                   perfdata[15,2])
    }
  }
}




##10kb


mccdata_gm12878_10kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_gm12878_10kb <- rbind(mccdata_gm12878_10kb,
                                    perfdata[15,2])
    }
  }
}




##25kb


mccdata_gm12878_25kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/25kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/25kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_gm12878_25kb <- rbind(mccdata_gm12878_25kb,
                                    perfdata[15,2])
    }
  }
}



##50kb


mccdata_gm12878_50kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/50kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/50kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_gm12878_50kb <- rbind(mccdata_gm12878_50kb,
                                    perfdata[15,2])
    }
  }
}




##100kb


mccdata_gm12878_100kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/GM12878/100kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/100kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_gm12878_100kb <- rbind(mccdata_gm12878_100kb,
                                     perfdata[15,2])
    }
  }
}


all_rfperfs <- rbind.data.frame(mccdata_gm12878_5kb,
                                mccdata_gm12878_10kb,
                                mccdata_gm12878_25kb,
                                mccdata_gm12878_50kb,
                                mccdata_gm12878_100kb)
all_rfperfs$CL <- c(rep("GM12878", length(mccdata_gm12878_10kb)*5))
all_rfperfs$Resolution <- c(rep("5 kb", length(mccdata_gm12878_5kb)),
                            rep("10 kb", length(mccdata_gm12878_10kb)),
                            rep("25 kb", length(mccdata_gm12878_25kb)),
                            rep("50 kb", length(mccdata_gm12878_50kb)),
                            rep("100 kb", length(mccdata_gm12878_100kb)))
all_rfperfs$Resampling <- c(rep("None", 3),
                            rep("ROS", 3),
                            rep("RUS", 3),
                            rep("SMOTE", 3))
all_rfperfs$Predictor <- c("Distance", "OC", "OP")
all_rfperfs$Chromosome <- c(rep("CHR1", 12),
                            rep("CHR2", 12),
                            rep("CHR3", 12),
                            rep("CHR4", 12),
                            rep("CHR5", 12),
                            rep("CHR6", 12),
                            rep("CHR7", 12),
                            rep("CHR8", 12),
                            rep("CHR10", 12),
                            rep("CHR11", 12),
                            rep("CHR12", 12),
                            rep("CHR13", 12),
                            rep("CHR14", 12),
                            rep("CHR15", 12),
                            rep("CHR16", 12),
                            rep("CHR17", 12),
                            rep("CHR18", 12),
                            rep("CHR19", 12),
                            rep("CHR20", 12),
                            rep("CHR21", 12),
                            rep("CHR22", 12))

all_rfperfs$Resolution <- factor(all_rfperfs$Resolution, levels = c("5 kb", "10 kb", "25 kb", "50 kb", "100 kb"))
all_rfperfs$Resampling <- factor(all_rfperfs$Resampling, levels = c("None", "ROS", "RUS", "SMOTE"))
all_rfperfs$Predictor <- factor(all_rfperfs$Predictor, levels = c("OC", "OP", "Distance"))
all_rfperfs$Chromosome <- factor(all_rfperfs$Chromosome, levels = paste0("CHR", c(1:8,10:22)))

all_rfperfs_mcc <- all_rfperfs %>%
  dplyr::group_by(Resolution, Resampling, Predictor) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))

ggplot(all_rfperfs_mcc, aes(x=Resolution, y=Performance, color=Resampling, group=Resampling)) +
  geom_line(size=2, position=position_dodge(0.5))+
  geom_point(size=5, position=position_dodge(0.5))+
  geom_errorbar(aes(ymin=Performance-PerformanceSD, ymax=Performance+PerformanceSD), width=1,
                position=position_dodge(0.5)) +
  facet_grid(. ~ Predictor)+
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + 
  ylab("") + #ylim(-.1,.5) +
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

ggplot(all_rfperfs_mcc, aes(x=Resolution, y=Performance, color=Predictor, fill=Predictor, group=Predictor)) +
  geom_bar(stat="identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=Performance-PerformanceSD, ymax=Performance+PerformanceSD), width=.75,
                position=position_dodge(.9), alpha=0.75, color="black", size=1) +
  facet_grid(. ~ Resampling)+
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + 
  ylab("MCC") + #ylim(-.1,.5) +
  scale_fill_discrete(name="Predictor Type")+
  guides(color=FALSE)+
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



####################################################################################

#K562 

##5kb

mccdata_k562_5kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_k562_5kb <- rbind(mccdata_k562_5kb,
                                perfdata[15,2])
    }
  }
}




##10kb


mccdata_k562_10kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_k562_10kb <- rbind(mccdata_k562_10kb,
                                 perfdata[15,2])
    }
  }
}




##25kb


mccdata_k562_25kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/25kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/25kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_k562_25kb <- rbind(mccdata_k562_25kb,
                                 perfdata[15,2])
    }
  }
}



##50kb


mccdata_k562_50kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/50kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/50kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_k562_50kb <- rbind(mccdata_k562_50kb,
                                 perfdata[15,2])
    }
  }
}




##100kb


mccdata_k562_100kb <- numeric()

for(i in paste0("CHR", c(1:8,10:22))){
  for(j in c("none","ros","rus","smote")){
    for(k in c("distance", "oc", "op")){
      print(c(i,j,k))
      files <- list.files(paste0("Z:/TAD_data_analysis/K562/100kb/results_by_chr/",
                                 i, 
                                 "/",
                                 j, 
                                 "/", 
                                 k, 
                                 "/"))
      if("TADRF.rds" %in% files){
        perfdata <- readRDS(paste0("Z:/TAD_data_analysis/K562/100kb/results_by_chr/",
                                   i, 
                                   "/",
                                   j, 
                                   "/", 
                                   k, 
                                   "/TADRF.rds"))
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
      
      mccdata_k562_100kb <- rbind(mccdata_k562_100kb,
                                  perfdata[15,2])
    }
  }
}


all_rfperfs <- rbind.data.frame(mccdata_k562_5kb,
                                mccdata_k562_10kb,
                                mccdata_k562_25kb,
                                mccdata_k562_50kb,
                                mccdata_k562_100kb)
all_rfperfs$CL <- c(rep("K562", length(mccdata_k562_10kb)*5))
all_rfperfs$Resolution <- c(rep("5 kb", length(mccdata_k562_5kb)),
                            rep("10 kb", length(mccdata_k562_10kb)),
                            rep("25 kb", length(mccdata_k562_25kb)),
                            rep("50 kb", length(mccdata_k562_50kb)),
                            rep("100 kb", length(mccdata_k562_100kb)))
all_rfperfs$Resampling <- c(rep("None", 3),
                            rep("ROS", 3),
                            rep("RUS", 3),
                            rep("SMOTE", 3))
all_rfperfs$Predictor <- c("Distance", "OC", "OP")
all_rfperfs$Chromosome <- c(rep("CHR1", 12),
                            rep("CHR2", 12),
                            rep("CHR3", 12),
                            rep("CHR4", 12),
                            rep("CHR5", 12),
                            rep("CHR6", 12),
                            rep("CHR7", 12),
                            rep("CHR8", 12),
                            rep("CHR10", 12),
                            rep("CHR11", 12),
                            rep("CHR12", 12),
                            rep("CHR13", 12),
                            rep("CHR14", 12),
                            rep("CHR15", 12),
                            rep("CHR16", 12),
                            rep("CHR17", 12),
                            rep("CHR18", 12),
                            rep("CHR19", 12),
                            rep("CHR20", 12),
                            rep("CHR21", 12),
                            rep("CHR22", 12))

all_rfperfs$Resolution <- factor(all_rfperfs$Resolution, levels = c("5 kb", "10 kb", "25 kb", "50 kb", "100 kb"))
all_rfperfs$Resampling <- factor(all_rfperfs$Resampling, levels = c("None", "ROS", "RUS", "SMOTE"))
all_rfperfs$Predictor <- factor(all_rfperfs$Predictor, levels = c("OC", "OP", "Distance"))
all_rfperfs$Chromosome <- factor(all_rfperfs$Chromosome, levels = paste0("CHR", c(1:8,10:22)))

all_rfperfs_mcc <- all_rfperfs %>%
  dplyr::group_by(Resolution, Resampling, Predictor) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))

ggplot(all_rfperfs_mcc, aes(x=Resolution, y=Performance, color=Resampling, group=Resampling)) +
  geom_line(size=2, position=position_dodge(0.5))+
  geom_point(size=5, position=position_dodge(0.5))+
  geom_errorbar(aes(ymin=Performance-PerformanceSD, ymax=Performance+PerformanceSD), width=1,
                position=position_dodge(0.5)) +
  facet_grid(. ~ Predictor)+
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + 
  ylab("MCC") + #ylim(-.1,.5) +
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

ggplot(all_rfperfs_mcc, aes(x=Resolution, y=Performance, color=Predictor, fill=Predictor, group=Predictor)) +
  geom_bar(stat="identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=Performance-PerformanceSD, ymax=Performance+PerformanceSD), width=.75,
                position=position_dodge(.9), alpha=0.75, color="black", size=1) +
  facet_grid(. ~ Resampling)+
  theme_minimal() +
  theme_bw()+
  xlab("Resolution") + 
  ylab("MCC") + #ylim(-.1,.5) +
  scale_fill_discrete(name="Predictor Type")+
  guides(color=FALSE)+
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

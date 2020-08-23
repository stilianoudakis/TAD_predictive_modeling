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


arrowhead_g_on_g <- readRDS("Z:/TAD_data_analysis/miscellaneous/arrowhead_g_on_g.rds")
arrowhead_g_on_g_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/arrowhead_g_on_g_auc.rds") 
arrowhead_g_on_k <- readRDS("Z:/TAD_data_analysis/miscellaneous/arrowhead_g_on_k.rds")
arrowhead_g_on_k_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/arrowhead_g_on_k_auc.rds") 
arrowhead_k_on_g <- readRDS("Z:/TAD_data_analysis/miscellaneous/arrowhead_k_on_g.rds")
arrowhead_k_on_g_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/arrowhead_k_on_g_auc.rds") 
arrowhead_k_on_k <- readRDS("Z:/TAD_data_analysis/miscellaneous/arrowhead_k_on_k.rds")
arrowhead_k_on_k_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/arrowhead_k_on_k_auc.rds") 


peakachu_g_on_g <- readRDS("Z:/TAD_data_analysis/miscellaneous/peakachu_g_on_g.rds")
peakachu_g_on_g_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/peakachu_g_on_g_auc.rds") 
peakachu_g_on_k <- readRDS("Z:/TAD_data_analysis/miscellaneous/peakachu_g_on_k.rds")
peakachu_g_on_k_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/peakachu_g_on_k_auc.rds") 
peakachu_k_on_k <- readRDS("Z:/TAD_data_analysis/miscellaneous/peakachu_k_on_k.rds")
peakachu_k_on_k_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/peakachu_k_on_k_auc.rds") 
peakachu_k_on_g <- readRDS("Z:/TAD_data_analysis/miscellaneous/peakachu_k_on_g.rds")
peakachu_k_on_g_auc <- readRDS("Z:/TAD_data_analysis/miscellaneous/peakachu_k_on_g_auc.rds") 


arrowhead_g_on_k <- do.call("rbind.data.frame", arrowhead_g_on_k)

arrowhead_g_on_g[[12]] <- rbind.data.frame(arrowhead_g_on_g[[12]],
                                           data.frame(specificity=1,sensitivity=0,chr="CHR13"),
                                           data.frame(specificity=1,sensitivity=0,chr="CHR13"),
                                           data.frame(specificity=1,sensitivity=0,chr="CHR13"))
arrowhead_g_on_g[[13]] <- rbind.data.frame(arrowhead_g_on_g[[13]],
                                           data.frame(specificity=1,sensitivity=0,chr="CHR14"))
arrowhead_g_on_g[[20]] <- rbind.data.frame(arrowhead_g_on_g[[20]],
                                           data.frame(specificity=1,sensitivity=0,chr="CHR21"),
                                           data.frame(specificity=1,sensitivity=0,chr="CHR21"),
                                           data.frame(specificity=1,sensitivity=0,chr="CHR21"),
                                           data.frame(specificity=1,sensitivity=0,chr="CHR21"))
arrowhead_g_on_g <- do.call("rbind.data.frame", arrowhead_g_on_g)

arrowhead_k_on_k[[15]] <- rbind.data.frame(arrowhead_k_on_k[[15]],
                                           data.frame(specificity=1,sensitivity=0,chr="CHR16"))
arrowhead_k_on_k[[17]] <- rbind.data.frame(arrowhead_k_on_k[[17]],
                                           data.frame(specificity=1,sensitivity=0,chr="CHR18"),
                                           data.frame(specificity=1,sensitivity=0,chr="CHR18"))
arrowhead_k_on_k[[20]] <- rbind.data.frame(arrowhead_k_on_k[[20]],
                                           data.frame(specificity=1,sensitivity=0,chr="CHR21"),
                                           data.frame(specificity=1,sensitivity=0,chr="CHR21"))
arrowhead_k_on_k <- do.call("rbind.data.frame", arrowhead_k_on_k)

arrowhead_k_on_g[[12]] <- rbind.data.frame(arrowhead_k_on_g[[12]],
                                           data.frame(specificity=1,sensitivity=0,chr="CHR13"))
arrowhead_k_on_g <- do.call("rbind.data.frame", arrowhead_k_on_g)



peakachu_g_on_g[[20]] <- rbind.data.frame(peakachu_g_on_g[[20]],
                                          data.frame(specificity=1,sensitivity=0,chr="CHR21"))
peakachu_g_on_g <- do.call("rbind.data.frame", peakachu_g_on_g)

peakachu_g_on_k[[20]] <- rbind.data.frame(peakachu_g_on_k[[20]],
                                          data.frame(specificity=1,sensitivity=0,chr="CHR21"))
peakachu_g_on_k <- do.call("rbind.data.frame", peakachu_g_on_k)

peakachu_k_on_g[[20]] <- rbind.data.frame(peakachu_k_on_g[[20]],
                                          data.frame(specificity=1,sensitivity=0,chr="CHR21"),
                                          data.frame(specificity=1,sensitivity=0,chr="CHR21"))
peakachu_k_on_g <- do.call("rbind.data.frame", peakachu_k_on_g)

peakachu_k_on_k <- do.call("rbind.data.frame", peakachu_k_on_k)


# GM12878 
##arrowhead
data <- rbind.data.frame(arrowhead_g_on_g, arrowhead_g_on_k)
data$id <- rep(factor(rep(c(1:502),21)),2)
data$CL <- c(rep("GM12878 on GM12878", nrow(arrowhead_g_on_g)),
             rep("GM12878 on K562", nrow(arrowhead_g_on_k)))

a1 <- reshape(data[which(data$CL=="GM12878 on GM12878"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "specificity")
a1 <- rowMeans(a1[,-c(1:3)])
b1 <- reshape(data[which(data$CL=="GM12878 on GM12878"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "sensitivity")
b1 <- rowMeans(b1[,-c(1:3)])
ave_sens_spec_gm12878_on_gm12878 <- data.frame(specificity=a1, sensitivity=b1)

a2 <- reshape(data[which(data$CL=="GM12878 on K562"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "specificity")
a2 <- rowMeans(a2[,-c(1:3)])
b2 <- reshape(data[which(data$CL=="GM12878 on K562"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "sensitivity")
b2 <- rowMeans(b2[,-c(1:3)])
ave_sens_spec_gm12878_on_k562 <- data.frame(specificity=a2, sensitivity=b2)

ave_sens_spec_arrow <- rbind.data.frame(ave_sens_spec_gm12878_on_gm12878,ave_sens_spec_gm12878_on_k562)
ave_sens_spec_arrow$Training_Testing <- c(rep("GM12878 on GM12878", nrow(ave_sens_spec_gm12878_on_gm12878)),
                                          rep("GM12878 on K562", nrow(ave_sens_spec_gm12878_on_k562)))
ave_sens_spec_arrow$Tool <- "Arrowhead"

##peakachu
data <- rbind.data.frame(peakachu_g_on_g, peakachu_g_on_k)
data$id <- rep(factor(rep(c(1:502),21)),2)
data$CL <- c(rep("GM12878 on GM12878", nrow(peakachu_g_on_g)),
             rep("GM12878 on K562", nrow(peakachu_g_on_k)))

a1 <- reshape(data[which(data$CL=="GM12878 on GM12878"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "specificity")
a1 <- rowMeans(a1[,-c(1:3)])
b1 <- reshape(data[which(data$CL=="GM12878 on GM12878"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "sensitivity")
b1 <- rowMeans(b1[,-c(1:3)])
ave_sens_spec_gm12878_on_gm12878 <- data.frame(specificity=a1, sensitivity=b1)

a2 <- reshape(data[which(data$CL=="GM12878 on K562"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "specificity")
a2 <- rowMeans(a2[,-c(1:3)])
b2 <- reshape(data[which(data$CL=="GM12878 on K562"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "sensitivity")
b2 <- rowMeans(b2[,-c(1:3)])
ave_sens_spec_gm12878_on_k562 <- data.frame(specificity=a2, sensitivity=b2)

ave_sens_spec_peak <- rbind.data.frame(ave_sens_spec_gm12878_on_gm12878,ave_sens_spec_gm12878_on_k562)
ave_sens_spec_peak$Training_Testing <- c(rep("GM12878 on GM12878", nrow(ave_sens_spec_gm12878_on_gm12878)),
                                         rep("GM12878 on K562", nrow(ave_sens_spec_gm12878_on_k562)))
ave_sens_spec_peak$Tool <- "Peakachu"

all_ave_sens_spec_df <- rbind.data.frame(ave_sens_spec_arrow,ave_sens_spec_peak)

ggplot(all_ave_sens_spec_df, aes(x=1-specificity, y=sensitivity, linetype=Training_Testing, color=Tool)) +
  facet_grid(. ~ Tool) +
  geom_text(
    label="AUC: ", 
    x=.70,
    y=.40,
    size = 5,
    color="black"
  ) +
  geom_rect(
    xmin = .55,
    xmax=.95,
    ymin=.3,
    ymax=.5,
    alpha=0,
    #color=Tool,
    fill="white"
  ) +
  geom_text(
    label="AUC: ", 
    x=.70,
    y=.15,
    size = 5,
    color="black"
  ) +
  geom_rect(
    xmin = .55,
    xmax=.95,
    ymin=.05,
    ymax=.25,
    alpha=0,
    #color=Tool,
    fill="white",
    linetype="dashed"
  ) +
  geom_line(size=1) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  scale_linetype_discrete(name="Training/Testing") +
  scale_color_manual(values = c("blue","red")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  guides(color=FALSE, linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")

##################################################################################################################

# K562
##arrowhead
data <- rbind.data.frame(arrowhead_k_on_k, arrowhead_k_on_g)
data$id <- rep(factor(rep(c(1:502),21)),2)
data$CL <- c(rep("K562 on K562", nrow(arrowhead_k_on_k)),
             rep("K562 on GM12878", nrow(arrowhead_k_on_g)))

a1 <- reshape(data[which(data$CL=="K562 on K562"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "specificity")
a1 <- rowMeans(a1[,-c(1:3)])
b1 <- reshape(data[which(data$CL=="K562 on K562"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "sensitivity")
b1 <- rowMeans(b1[,-c(1:3)])
ave_sens_spec_k562_on_k562 <- data.frame(specificity=a1, sensitivity=b1)

a2 <- reshape(data[which(data$CL=="K562 on GM12878"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "specificity")
a2 <- rowMeans(a2[,-c(1:3)])
b2 <- reshape(data[which(data$CL=="K562 on GM12878"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "sensitivity")
b2 <- rowMeans(b2[,-c(1:3)])
ave_sens_spec_k562_on_gm12878 <- data.frame(specificity=a2, sensitivity=b2)

ave_sens_spec_arrow <- rbind.data.frame(ave_sens_spec_k562_on_k562,ave_sens_spec_k562_on_gm12878)
ave_sens_spec_arrow$Training_Testing <- c(rep("K562 on K562", nrow(ave_sens_spec_k562_on_k562)),
                                          rep("K562 on GM12878", nrow(ave_sens_spec_k562_on_gm12878)))
ave_sens_spec_arrow$Tool <- "Arrowhead"

##peakachu
data <- rbind.data.frame(peakachu_k_on_k, peakachu_k_on_g)
data$id <- rep(factor(rep(c(1:502),21)),2)
data$CL <- c(rep("K562 on K562", nrow(peakachu_k_on_k)),
             rep("K562 on GM12878", nrow(peakachu_k_on_g)))

a1 <- reshape(data[which(data$CL=="K562 on K562"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "specificity")
a1 <- rowMeans(a1[,-c(1:3)])
b1 <- reshape(data[which(data$CL=="K562 on K562"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "sensitivity")
b1 <- rowMeans(b1[,-c(1:3)])
ave_sens_spec_k562_on_k562 <- data.frame(specificity=a1, sensitivity=b1)

a2 <- reshape(data[which(data$CL=="K562 on GM12878"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "specificity")
a2 <- rowMeans(a2[,-c(1:3)])
b2 <- reshape(data[which(data$CL=="K562 on GM12878"),], 
              idvar = "id",
              timevar="chr",
              direction="wide",
              v.names = "sensitivity")
b2 <- rowMeans(b2[,-c(1:3)])
ave_sens_spec_k562_on_gm12878 <- data.frame(specificity=a2, sensitivity=b2)

ave_sens_spec_peak <- rbind.data.frame(ave_sens_spec_k562_on_k562,ave_sens_spec_k562_on_gm12878)
ave_sens_spec_peak$Training_Testing <- c(rep("K562 on K562", nrow(ave_sens_spec_k562_on_k562)),
                                         rep("K562 on GM12878", nrow(ave_sens_spec_k562_on_gm12878)))
ave_sens_spec_peak$Tool <- "Peakachu"

all_ave_sens_spec_df <- rbind.data.frame(ave_sens_spec_arrow,ave_sens_spec_peak)

ggplot(all_ave_sens_spec_df, aes(x=1-specificity, y=sensitivity, linetype=Training_Testing, color=Tool)) +
  facet_grid(. ~ Tool) +
  geom_text(
    label="AUC: ", 
    x=.70,
    y=.40,
    size = 5,
    color="black"
  ) +
  geom_rect(
    xmin = .55,
    xmax=.95,
    ymin=.3,
    ymax=.5,
    alpha=0,
    #color=Tool,
    fill="white"
  ) +
  geom_text(
    label="AUC: ", 
    x=.70,
    y=.15,
    size = 5,
    color="black"
  ) +
  geom_rect(
    xmin = .55,
    xmax=.95,
    ymin=.05,
    ymax=.25,
    alpha=0,
    #color=Tool,
    fill="white",
    linetype="dashed"
  ) +
  geom_line(size=1) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  scale_linetype_discrete(name="Training/Testing") +
  scale_color_manual(values = c("blue","red")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  guides(color=FALSE, linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")

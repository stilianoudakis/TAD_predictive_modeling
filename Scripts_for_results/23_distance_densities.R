#23_distance_densities

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
#library(GenometriCorr)
library(ggsignif)
library(Vennerable)
library(VennDiagram)
library(ChIPpeakAnno)
library(EnrichedHeatmap)

source("Z:/TAD_data_analysis/functions_for_R_package/createTADdata.R")
source("Z:/TAD_data_analysis/functions_for_R_package/preciseTAD.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADrandomForest.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADrfe2.R")
source("Z:/TAD_data_analysis/functions_for_R_package/distance_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/percent_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/count_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/binary_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/signal_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/annots_to_granges_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/extract_boundaries_func.R")

domains <- read.table("Z:/TAD_data_analysis/GM12878/5kb/GM12878_domain_data_5000.b.txt")
bounds.GR <- extract_boundaries_func(domains.mat=domains, 
                                     preprocess=FALSE, 
                                     CHR=paste0("CHR", c(1:8,10:22)), 
                                     resolution=5000)

genomicElements.GR.gm12878 <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/TFBS",
                                             pattern="*.bed",
                                             signal=4)

tadData_gm12878_5kb <- createTADdata(bounds.GR=bounds.GR,
                         resolution=5000,
                         genomicElements.GR=genomicElements.GR.gm12878,
                         featureType="distance",
                         resampling="none",
                         trainCHR=paste0("CHR",c(1:8,10:22)),
                         predictCHR=NULL,
                         seed=123)

tadData_gm12878_10kb <- createTADdata(bounds.GR=bounds.GR,
                                     resolution=10000,
                                     genomicElements.GR=genomicElements.GR.gm12878,
                                     featureType="distance",
                                     resampling="none",
                                     trainCHR=paste0("CHR",c(1:8,10:22)),
                                     predictCHR=NULL,
                                     seed=123)

tadData_gm12878_25kb <- createTADdata(bounds.GR=bounds.GR,
                                     resolution=25000,
                                     genomicElements.GR=genomicElements.GR.gm12878,
                                     featureType="distance",
                                     resampling="none",
                                     trainCHR=paste0("CHR",c(1:8,10:22)),
                                     predictCHR=NULL,
                                     seed=123)

tadData_gm12878_50kb <- createTADdata(bounds.GR=bounds.GR,
                                     resolution=50000,
                                     genomicElements.GR=genomicElements.GR.gm12878,
                                     featureType="distance",
                                     resampling="none",
                                     trainCHR=paste0("CHR",c(1:8,10:22)),
                                     predictCHR=NULL,
                                     seed=123)

tadData_gm12878_100kb <- createTADdata(bounds.GR=bounds.GR,
                                     resolution=100000,
                                     genomicElements.GR=genomicElements.GR.gm12878,
                                     featureType="distance",
                                     resampling="none",
                                     trainCHR=paste0("CHR",c(1:8,10:22)),
                                     predictCHR=NULL,
                                     seed=123)

genomicElements.GR.k562 <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/k562/TFBS",
                                                  pattern="*.bed",
                                                  signal=4)

tadData_k562_5kb <- createTADdata(bounds.GR=bounds.GR,
                                  resolution=5000,
                                  genomicElements.GR=genomicElements.GR.k562,
                                  featureType="distance",
                                  resampling="none",
                                  trainCHR=paste0("CHR",c(1:8,10:22)),
                                  predictCHR=NULL,
                                  seed=123)

tadData_k562_10kb <- createTADdata(bounds.GR=bounds.GR,
                                   resolution=10000,
                                   genomicElements.GR=genomicElements.GR.k562,
                                   featureType="distance",
                                   resampling="none",
                                   trainCHR=paste0("CHR",c(1:8,10:22)),
                                   predictCHR=NULL,
                                   seed=123)

tadData_k562_25kb <- createTADdata(bounds.GR=bounds.GR,
                                   resolution=25000,
                                   genomicElements.GR=genomicElements.GR.k562,
                                   featureType="distance",
                                   resampling="none",
                                   trainCHR=paste0("CHR",c(1:8,10:22)),
                                   predictCHR=NULL,
                                   seed=123)

tadData_k562_50kb <- createTADdata(bounds.GR=bounds.GR,
                                   resolution=50000,
                                   genomicElements.GR=genomicElements.GR.k562,
                                   featureType="distance",
                                   resampling="none",
                                   trainCHR=paste0("CHR",c(1:8,10:22)),
                                   predictCHR=NULL,
                                   seed=123)

tadData_k562_100kb <- createTADdata(bounds.GR=bounds.GR,
                                    resolution=100000,
                                    genomicElements.GR=genomicElements.GR.k562,
                                    featureType="distance",
                                    resampling="none",
                                    trainCHR=paste0("CHR",c(1:8,10:22)),
                                    predictCHR=NULL,
                                    seed=123)

tadData_gm12878_5kb <- melt(tadData_gm12878_5kb[[1]][,-1])
tadData_gm12878_10kb <- melt(tadData_gm12878_10kb[[1]][,-1])
tadData_gm12878_25kb <- melt(tadData_gm12878_25kb[[1]][,-1])
tadData_gm12878_50kb <- melt(tadData_gm12878_50kb[[1]][,-1])
tadData_gm12878_100kb <- melt(tadData_gm12878_100kb[[1]][,-1])

tadData_k562_5kb <- melt(tadData_k562_5kb[[1]][,-1])
tadData_k562_10kb <- melt(tadData_k562_10kb[[1]][,-1])
tadData_k562_25kb <- melt(tadData_k562_25kb[[1]][,-1])
tadData_k562_50kb <- melt(tadData_k562_50kb[[1]][,-1])
tadData_k562_100kb <- melt(tadData_k562_100kb[[1]][,-1])

tadData_gm12878_5kb_nolog <- tadData_gm12878_5kb 
tadData_gm12878_10kb_nolog <- tadData_gm12878_10kb 
tadData_gm12878_25kb_nolog <- tadData_gm12878_25kb 
tadData_gm12878_50kb_nolog <- tadData_gm12878_50kb 
tadData_gm12878_100kb_nolog <- tadData_gm12878_100kb 
tadData_k562_5kb_nolog <- tadData_k562_5kb 
tadData_k562_10kb_nolog <- tadData_k562_10kb 
tadData_k562_25kb_nolog <- tadData_k562_25kb 
tadData_k562_50kb_nolog <- tadData_k562_50kb 
tadData_k562_100kb_nolog <- tadData_k562_100kb 

tadData_gm12878_5kb_nolog$value <- exp(tadData_gm12878_5kb_nolog$value) - 1
tadData_gm12878_10kb_nolog$value <- exp(tadData_gm12878_10kb_nolog$value) - 1
tadData_gm12878_25kb_nolog$value <- exp(tadData_gm12878_25kb_nolog$value) - 1
tadData_gm12878_50kb_nolog$value <- exp(tadData_gm12878_50kb_nolog$value) - 1
tadData_gm12878_100kb_nolog$value <- exp(tadData_gm12878_100kb_nolog$value) - 1
tadData_k562_5kb_nolog$value <- exp(tadData_k562_5kb_nolog$value) - 1
tadData_k562_10kb_nolog$value <- exp(tadData_k562_10kb_nolog$value) - 1
tadData_k562_25kb_nolog$value <- exp(tadData_k562_25kb_nolog$value) - 1
tadData_k562_50kb_nolog$value <- exp(tadData_k562_50kb_nolog$value) - 1
tadData_k562_100kb_nolog$value <- exp(tadData_k562_100kb_nolog$value) - 1

#gm12878

##untransformed

###5kb

ggplot(tadData_gm12878_5kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###10kb

ggplot(tadData_gm12878_10kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###25kb

ggplot(tadData_gm12878_25kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###50kb

ggplot(tadData_gm12878_50kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###100kb

ggplot(tadData_gm12878_100kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

##log2 transformed

###5kb

ggplot(tadData_gm12878_5kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###10kb

ggplot(tadData_gm12878_10kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###25kb

ggplot(tadData_gm12878_25kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###50kb

ggplot(tadData_gm12878_50kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###100kb

ggplot(tadData_gm12878_100kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")


#k562

##untransformed

###5kb

ggplot(tadData_k562_5kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###10kb

ggplot(tadData_k562_10kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###25kb

ggplot(tadData_k562_25kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###50kb

ggplot(tadData_k562_50kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###100kb

ggplot(tadData_k562_100kb_nolog, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="red") +
  xlab("Genomic Distance (mb)") + 
  ylab("Density") +
  xlim(0,30000000) +
  scale_x_continuous(limits = c(0,3e+7),
                     breaks = c(0, 1e+07, 2e+07, 3e+07),
                     labels = c(0,1,2,3))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

##log2 transformed

###5kb

ggplot(tadData_k562_5kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###10kb

ggplot(tadData_k562_10kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###25kb

ggplot(tadData_k562_25kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###50kb

ggplot(tadData_k562_50kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")

###100kb

ggplot(tadData_k562_100kb, aes(x=value, group=variable)) +
  geom_density(alpha=.5, color="blue") +
  xlab("Log2 Genomic Distance") + 
  ylab("Density") +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.title = element_text(size=15),
        legend.position = "none")
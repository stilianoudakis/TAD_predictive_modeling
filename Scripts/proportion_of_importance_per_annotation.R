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

#GM12878

annots <- list.files("Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878", pattern = "*.bed")
annots <- gsub(paste(c("Gm12878-",".bed", 
                       "-Haib", 
                       "-Sydh",
                       "-Uta",
                       "-Uw",
                       "-Broad",
                       "-Uchicago"), collapse = "|"), "", annots)
annots <- annots[!duplicated(annots)]

chromstates <- list.files("Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/BroadHMM")
chromstates <- chromstates[grep("Gm12878", chromstates)]
chromstates <- gsub(paste0(c("Gm12878-", ".bed"), collapse = "|"),"",chromstates)

histones <- list.files("Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/HistoneModifications")
histones <- histones[grep("Gm12878", histones)]
histones <- gsub(paste0(c("Gm12878-", ".bed"), collapse = "|"),"",histones)

tfbs <- list.files("Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/TFBS")
tfbs <- gsub(paste0(c("Gm12878-", ".bed", "-Haib", "-Sydh", "-Broad", "-Uw", "-Uta", "-Chicago"), collapse = "|"),"", tfbs)
tfbs <- tfbs[!duplicated(tfbs)]

##5kb

varimp_gm12878_5kb <- matrix(nrow=length(annots), ncol=21)
rownames(varimp_gm12878_5kb) <- annots

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  varimps <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                            chrs[i], 
                            "/rus/distance/TADRF.rds"))
  varimps <- varimps[[2]]
  varimps$Feature <- gsub(paste0(c("`",
                                   "_dist",
                                   "Gm12878-",".bed", 
                                   "-Haib", 
                                   "-Sydh",
                                   "-Uta",
                                   "-Uw",
                                   "-Broad",
                                   "-Uchicago"), collapse = "|"), "", varimps$Feature)	
  varimps <- varimps[order(varimps$Feature),]
  
  varimp_gm12878_5kb[match(varimps$Feature,annots),i] <- varimps$Importance
}

varimp_gm12878_5kb <- data.frame(Feature = rownames(varimp_gm12878_5kb),
                                 aveimp5kb = rowMeans(varimp_gm12878_5kb, na.rm = TRUE),
                                 sdimp5kb = apply(varimp_gm12878_5kb,1, sd, na.rm = TRUE))
rownames(varimp_gm12878_5kb) <- NULL

##10kb

varimp_gm12878_10kb <- matrix(nrow=length(annots), ncol=21)
rownames(varimp_gm12878_10kb) <- annots

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  varimps <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                            chrs[i], 
                            "/rus/distance/TADRF.rds"))
  varimps <- varimps[[2]]
  varimps$Feature <- gsub(paste0(c("`",
                                   "_dist",
                                   "Gm12878-",".bed", 
                                   "-Haib", 
                                   "-Sydh",
                                   "-Uta",
                                   "-Uw",
                                   "-Broad",
                                   "-Uchicago"), collapse = "|"), "", varimps$Feature)	
  varimps <- varimps[order(varimps$Feature),]
  
  varimp_gm12878_10kb[match(varimps$Feature,annots),i] <- varimps$Importance
}

varimp_gm12878_10kb <- data.frame(Feature = rownames(varimp_gm12878_10kb),
                                  aveimp10kb = rowMeans(varimp_gm12878_10kb, na.rm = TRUE),
                                  sdimp10kb = apply(varimp_gm12878_10kb,1, sd, na.rm = TRUE))
rownames(varimp_gm12878_10kb) <- NULL

##25kb

varimp_gm12878_25kb <- matrix(nrow=length(annots), ncol=21)
rownames(varimp_gm12878_25kb) <- annots

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  varimps <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/25kb/results_by_chr/",
                            chrs[i], 
                            "/rus/distance/TADRF.rds"))
  varimps <- varimps[[2]]
  varimps$Feature <- gsub(paste0(c("`",
                                   "_dist",
                                   "Gm12878-",".bed", 
                                   "-Haib", 
                                   "-Sydh",
                                   "-Uta",
                                   "-Uw",
                                   "-Broad",
                                   "-Uchicago"), collapse = "|"), "", varimps$Feature)
  varimps <- varimps[order(varimps$Feature),]
  
  varimp_gm12878_25kb[match(varimps$Feature,annots),i] <- varimps$Importance
}

varimp_gm12878_25kb <- data.frame(Feature = rownames(varimp_gm12878_25kb),
                                  aveimp25kb = rowMeans(varimp_gm12878_25kb, na.rm = TRUE),
                                  sdimp25kb = apply(varimp_gm12878_25kb,1, sd, na.rm = TRUE))
rownames(varimp_gm12878_25kb) <- NULL

##50kb

varimp_gm12878_50kb <- matrix(nrow=length(annots), ncol=21)
rownames(varimp_gm12878_50kb) <- annots

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  varimps <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/50kb/results_by_chr/",
                            chrs[i], 
                            "/rus/distance/TADRF.rds"))
  varimps <- varimps[[2]]
  varimps$Feature <- gsub(paste0(c("`",
                                   "_dist",
                                   "Gm12878-",".bed", 
                                   "-Haib", 
                                   "-Sydh",
                                   "-Uta",
                                   "-Uw",
                                   "-Broad",
                                   "-Uchicago"), collapse = "|"), "", varimps$Feature)
  varimps <- varimps[order(varimps$Feature),]
  
  varimp_gm12878_50kb[match(varimps$Feature,annots),i] <- varimps$Importance
}

varimp_gm12878_50kb <- data.frame(Feature = rownames(varimp_gm12878_50kb),
                                  aveimp50kb = rowMeans(varimp_gm12878_50kb, na.rm = TRUE),
                                  sdimp50kb = apply(varimp_gm12878_50kb,1, sd, na.rm = TRUE))
rownames(varimp_gm12878_50kb) <- NULL

##100kb

varimp_gm12878_100kb <- matrix(nrow=length(annots), ncol=21)
rownames(varimp_gm12878_100kb) <- annots

chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  varimps <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/100kb/results_by_chr/",
                            chrs[i], 
                            "/rus/distance/TADRF.rds"))
  varimps <- varimps[[2]]
  varimps$Feature <- gsub(paste0(c("`",
                                   "_dist",
                                   "Gm12878-",".bed", 
                                   "-Haib", 
                                   "-Sydh",
                                   "-Uta",
                                   "-Uw",
                                   "-Broad",
                                   "-Uchicago"), collapse = "|"), "", varimps$Feature)
  varimps <- varimps[order(varimps$Feature),]
  
  varimp_gm12878_100kb[match(varimps$Feature,annots),i] <- varimps$Importance
}

varimp_gm12878_100kb <- data.frame(Feature = rownames(varimp_gm12878_100kb),
                                   aveimp100kb = rowMeans(varimp_gm12878_100kb, na.rm = TRUE),
                                   sdimp100kb = apply(varimp_gm12878_100kb,1, sd, na.rm = TRUE))
rownames(varimp_gm12878_100kb) <- NULL

all_propimp_gm12878 <- data.frame(Importance=c(varimp_gm12878_5kb$aveimp5kb,
                                               varimp_gm12878_10kb$aveimp10kb,
                                               varimp_gm12878_25kb$aveimp25kb,
                                               varimp_gm12878_50kb$aveimp50kb,
                                               varimp_gm12878_100kb$aveimp100kb),
                                  Feature=c(as.character(varimp_gm12878_10kb$Feature),
                                            as.character(varimp_gm12878_10kb$Feature),
                                            as.character(varimp_gm12878_10kb$Feature),
                                            as.character(varimp_gm12878_10kb$Feature),
                                            as.character(varimp_gm12878_10kb$Feature)),
                                  Resolution=c(rep("5 kb", 78),
                                               rep("10 kb", 78),
                                               rep("25 kb", 78),
                                               rep("50 kb", 78),
                                               rep("100 kb", 78))
)

all_propimp_gm12878$Annotation <- ifelse(all_propimp_gm12878$Feature %in% chromstates, "CRE",
                                         ifelse(all_propimp_gm12878$Feature %in% histones, "HM", "TFBS"))

all_propimp_gm12878_sum <- all_propimp_gm12878 %>%
  dplyr::group_by(Resolution, Annotation) %>%
  dplyr::summarise(Performance = mean(Importance, na.rm = TRUE),
            PerformanceSD = sd(Importance, na.rm = TRUE))

ggplot(all_propimp_gm12878_sum, aes(fill=Annotation, y=Performance, x=Resolution)) + 
  geom_bar( stat="identity", position="fill") +
  xlab("Resolution") +
  ylab("Proportion of Importance") +
  scale_fill_discrete(name="Genomic Class")+
  theme_bw() +
  theme(axis.text.x = element_text(size=15,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20))

#####################################################################################
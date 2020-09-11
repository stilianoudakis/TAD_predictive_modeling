#18_saving_ptbrs_ptbps

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

OUT <- createWorkbook()

#arrowhead

##gm12878

###ptbr

ptbr <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  ptbr[[i]] <- called_and_pred$PTBR
}

ptbr_arrow_gm12878 <- do.call("c", ptbr)
ptbr_arrow_gm12878_df <- data.frame(Chromosome = as.character(seqnames(ptbr_arrow_gm12878)),
                                    Start = start(ptbr_arrow_gm12878),
                                    End = end(ptbr_arrow_gm12878))
addWorksheet(OUT, "arrowhead_gm12878_ptbr")
writeData(OUT, sheet = "arrowhead_gm12878_ptbr", x = ptbr_arrow_gm12878_df)

###ptbp

ptbp <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  ptbp[[i]] <- called_and_pred$PTBP
}

ptbp_arrow_gm12878 <- do.call("c", ptbp)
ptbp_arrow_gm12878_df <- data.frame(Chromosome = as.character(seqnames(ptbr_arrow_gm12878)),
                                    Start = start(ptbp_arrow_gm12878),
                                    End = start(ptbp_arrow_gm12878)+1)
addWorksheet(OUT, "arrowhead_gm12878_ptbp")
writeData(OUT, sheet = "arrowhead_gm12878_ptbp", x = ptbp_arrow_gm12878_df)

##k562

###ptbr

ptbr <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  ptbr[[i]] <- called_and_pred$PTBR
}

ptbr_arrow_k562 <- do.call("c", ptbr)
ptbr_arrow_k562_df <- data.frame(Chromosome = as.character(seqnames(ptbr_arrow_k562)),
                                 Start = start(ptbr_arrow_k562),
                                 End = end(ptbr_arrow_k562))
addWorksheet(OUT, "arrowhead_k562_ptbr")
writeData(OUT, sheet = "arrowhead_k562_ptbr", x = ptbr_arrow_k562_df)

###ptbp

ptbp <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout.rds"))
  ptbp[[i]] <- called_and_pred$PTBP
}

ptbp_arrow_k562 <- do.call("c", ptbp)
ptbp_arrow_k562_df <- data.frame(Chromosome = as.character(seqnames(ptbr_arrow_k562)),
                                 Start = start(ptbp_arrow_k562),
                                 End = start(ptbp_arrow_k562)+1)
addWorksheet(OUT, "arrowhead_k562_ptbp")
writeData(OUT, sheet = "arrowhead_k562_ptbp", x = ptbp_arrow_k562_df)

######################################################################

#peakachu

##gm12878

###ptbr

ptbr <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  ptbr[[i]] <- called_and_pred$PTBR
}

ptbr_peak_gm12878 <- do.call("c", ptbr)
ptbr_peak_gm12878_df <- data.frame(Chromosome = as.character(seqnames(ptbr_peak_gm12878)),
                                   Start = start(ptbr_peak_gm12878),
                                   End = end(ptbr_peak_gm12878))
addWorksheet(OUT, "peakachu_gm12878_ptbr")
writeData(OUT, sheet = "peakachu_gm12878_ptbr", x = ptbr_peak_gm12878_df)

###ptbp

ptbp <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  ptbp[[i]] <- called_and_pred$PTBP
}

ptbp_peak_gm12878 <- do.call("c", ptbp)
ptbp_peak_gm12878_df <- data.frame(Chromosome = as.character(seqnames(ptbr_peak_gm12878)),
                                   Start = start(ptbp_peak_gm12878),
                                   End = start(ptbp_peak_gm12878)+1)
addWorksheet(OUT, "peakachu_gm12878_ptbp")
writeData(OUT, sheet = "peakachu_gm12878_ptbp", x = ptbp_peak_gm12878_df)

##k562

###ptbr

ptbr <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  ptbr[[i]] <- called_and_pred$PTBR
}

ptbr_peak_k562 <- do.call("c", ptbr)
ptbr_peak_k562_df <- data.frame(Chromosome = as.character(seqnames(ptbr_peak_k562)),
                                Start = start(ptbr_peak_k562),
                                End = end(ptbr_peak_k562))
addWorksheet(OUT, "peakachu_k562_ptbr")
writeData(OUT, sheet = "peakachu_k562_ptbr", x = ptbr_peak_k562_df)

###ptbp

ptbp <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_peakachu.rds"))
  ptbp[[i]] <- called_and_pred$PTBP
}

ptbp_peak_k562 <- do.call("c", ptbp)
ptbp_peak_k562_df <- data.frame(Chromosome = as.character(seqnames(ptbr_peak_k562)),
                                Start = start(ptbp_peak_k562),
                                End = start(ptbp_peak_k562)+1)
addWorksheet(OUT, "peakachu_k562_ptbp")
writeData(OUT, sheet = "peakachu_k562_ptbp", x = ptbp_peak_k562_df)



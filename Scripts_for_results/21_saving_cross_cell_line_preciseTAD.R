#21_saving_cross_cell_line_preciseTAD

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

#arrowhead

##g on k

g_on_k <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_gm12878_on_k562.rds"))
  
  g_on_k[[i]] <- called_and_pred$PTBP
  
}
g_on_k <- do.call("c", g_on_k)

g_on_k_df <- data.frame(chr=as.character(seqnames(g_on_k)),
                                 start=start(g_on_k),
                                 end=start(g_on_k)+1)

write.table(g_on_k_df, "C:/Users/stili/Documents/TAD_miscellaneous/data/bedfiles_for_deeptools/arrowhead_preciseTAD_g_on_k.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

##k on g

k_on_g <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_k562_on_gm12878.rds"))
  
  k_on_g[[i]] <- called_and_pred$PTBP
  
}
k_on_g <- do.call("c", k_on_g)

k_on_g_df <- data.frame(chr=as.character(seqnames(k_on_g)),
                        start=start(k_on_g),
                        end=start(k_on_g)+1)

write.table(k_on_g_df, "C:/Users/stili/Documents/TAD_miscellaneous/data/bedfiles_for_deeptools/arrowhead_preciseTAD_k_on_g.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

#peakachu

##g on k

g_on_k <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_gm12878_on_k562_peakachu.rds"))
  
  g_on_k[[i]] <- called_and_pred$PTBP
  
}
g_on_k <- do.call("c", g_on_k)

g_on_k_df <- data.frame(chr=as.character(seqnames(g_on_k)),
                        start=start(g_on_k),
                        end=start(g_on_k)+1)

write.table(g_on_k_df, "C:/Users/stili/Documents/TAD_miscellaneous/data/bedfiles_for_deeptools/peakachu_preciseTAD_g_on_k.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

##k on g

k_on_g <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/10kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD_holdout_k562_on_gm12878_peakachu.rds"))
  
  k_on_g[[i]] <- called_and_pred$PTBP
  
}
k_on_g <- do.call("c", k_on_g)

k_on_g_df <- data.frame(chr=as.character(seqnames(k_on_g)),
                        start=start(k_on_g),
                        end=start(k_on_g)+1)

write.table(k_on_g_df, "C:/Users/stili/Documents/TAD_miscellaneous/data/bedfiles_for_deeptools/peakachu_preciseTAD_k_on_g.bed",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

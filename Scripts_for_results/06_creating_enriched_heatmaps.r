#06 Creating enriched heatmaps

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

#### Create granges from top tfbs

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS",
                                             pattern="*.bed",
											 signal=4)

##### Extracting called and predicted boundaries

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
chrs <- paste0("CHR", c(2:8,10:15,17:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[3]]
  pred_bound_list[[i]] <- called_and_pred[[2]]
  
}

all_true_bounds <- do.call("c", true_bound_list)
all_pred_bounds <- do.call("c", pred_bound_list)

length(all_true_bounds) #
length(all_pred_bounds) #

###### Enriched Heatmaps

enriched_func <- function(boundaries, annotations, e=5, w=50){
  axis_names <- c(paste0("-",e,"kb"), "Center", paste0(e,"kb"))
  
  norm_mat <- normalizeToMatrix(annotations, 
                         boundaries, 
                         value_column = "coverage", 
                         extend = e*1000, 
						 w = 50,
                         mean_mode = "w0",
                         keep = c(0,.999))
  EH <- EnrichedHeatmap(norm_mat, 
                col = c("blue",colorRampPalette(c("white", "pink", "red"))(99)), 
                #heatmap_legend_param = list(show_heatmap_legend = FALSE),
                row_title="ChIP-seq Peaks",
                show_heatmap_legend = FALSE,
                row_title_side="left", 
                column_title="", 
                axis_name = axis_names,
				        column_title_gp = gpar(fontsize = 10),
				        row_title_gp = gpar(fontsize = 10),
				        axis_name_gp = gpar(fontsize = 10),
                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col ="red", 
                                                                                   fontsize = 10,
                                                                                   lwd = 2), 
                                                                         axis_param = list(side = "right",
                                                                                           facing = "inside"))))
  return(EH)
}

####### CTCF

enriched_func(all_pred_bounds, genomicElements.GR[[1]])
enriched_func(all_true_bounds, genomicElements.GR[[1]])

####### RAD21

enriched_func(all_pred_bounds, genomicElements.GR[[2]])
enriched_func(all_true_bounds, genomicElements.GR[[2]])

####### SMC3

enriched_func(all_pred_bounds, genomicElements.GR[[3]])
enriched_func(all_true_bounds, genomicElements.GR[[3]])

####### ZNF143

enriched_func(all_pred_bounds, genomicElements.GR[[4]])
enriched_func(all_true_bounds, genomicElements.GR[[4]])

############################################################################

### K562

#### Create granges from top tfbs

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/k562/topTFBS",
                                             pattern="*.bed",
											 signal=4)

##### Extracting called and predicted boundaries

pred_bound_list <- GRangesList()
true_bound_list <- GRangesList()
chrs <- paste0("CHR", c(1:8,10:22))
chrs <- paste0("CHR", c(2:8,10:15,17:22))
for(i in 1:length(chrs)){
  called_and_pred <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/",
                                    chrs[i],
                                    "/rus/distance/preciseTAD.rds"))
  
  true_bound_list[[i]] <- called_and_pred[[3]]
  pred_bound_list[[i]] <- called_and_pred[[2]]
  
}

all_true_bounds <- do.call("c", true_bound_list)
all_pred_bounds <- do.call("c", pred_bound_list)

length(all_true_bounds) #
length(all_pred_bounds) #

###### Enriched Heatmaps

####### CTCF

enriched_func(all_pred_bounds, genomicElements.GR[[1]])
enriched_func(all_true_bounds, genomicElements.GR[[1]])

####### RAD21

enriched_func(all_pred_bounds, genomicElements.GR[[2]])
enriched_func(all_true_bounds, genomicElements.GR[[2]])

####### SMC3

enriched_func(all_pred_bounds, genomicElements.GR[[3]])
enriched_func(all_true_bounds, genomicElements.GR[[3]])

####### ZNF143

enriched_func(all_pred_bounds, genomicElements.GR[[4]])
enriched_func(all_true_bounds, genomicElements.GR[[4]])


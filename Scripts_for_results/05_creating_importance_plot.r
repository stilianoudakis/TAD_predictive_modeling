#05 Creating importance plot

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

#### Importance barchart/heatmap

rfeModelResultsList <- readRDS(paste0("Z:/TAD_data_analysis/GM12878/5kb/results_by_chr/WHOLE/rus/distance/TADrfe_3.rds"))

set.size = 8
rfeModelResults <- rfeModelResultsList[[2]]
rfe.set <- rfeModelResults[rfeModelResults$Variables==set.size,]
rfe.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)
rfe.order <- order(rfe.set[, c("x")], decreasing = TRUE)
rfe.set <- rfe.set[rfe.order, ]
rfe.set$Group.1 <- factor(rfe.set$Group.1, levels=c("Smc3ab9263", 
                                                    "Znf143",
                                                    "Rad21",
                                                    "Ctcf",
                                                    "Pol2",
                                                    "Elf1sc631",  
                                                    "Pol24h8",    
                                                    "Tbp",        
                                                    "Taf1",      
                                                    "Mef2a",      
                                                    "Mxi1",       
                                                    "Mazab85725"))
ggplot(rfe.set, aes(x=reorder(Group.1,x), y=x)) +
  geom_bar(stat='identity', fill="red") +
  coord_flip() +
  xlab("Transcription Factor Binding Sites") +
  ylab("Predictive Importance") +
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


rfeModelResults <- rfeModelResultsList[[2]]
rfe.set <- rfeModelResults[rfeModelResults$Variables==set.size,]
rfe.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var,rfe.set$Resample), mean)
rfe.set <- reshape(rfe.set, idvar = "Group.1", timevar = "Group.2", direction = "wide")
rfe.set.mat <- as.matrix(rfe.set[,-1])
rownames(rfe.set.mat) <- rfe.set$Group.1
rfe.set.mat[is.na(rfe.set.mat)] <- 0
rfe.set.mat <- rfe.set.mat[ order(rowMeans(rfe.set.mat)), ]

rfe.set.mat <- apply(rfe.set.mat,2,rev)

distance.col = dist(rfe.set.mat, method = "euclidean")
cluster.col = hclust(distance.col, method = "average")
my_palette <- colorRampPalette(c("white", "pink", "red"))(ncol(rfe.set.mat))

heatmap.2(t(rfe.set.mat),
          #labRow=rownames(rfe.set.mat),
          #labCol=colnames(rfe.set.mat),
          #cellnote = ,  # same data set for cell labels
          #main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins=c(12,8),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          na.color = "gray",
          #breaks=col_breaks,    # enable color transition at specified limits
          key=FALSE,
          #keysize = 1.5,
          #key.title = "Elastic-Net Coefficient",
          dendrogram="column",     # only draw a row dendrogram,
          Rowv = FALSE,
          Colv = reorder(as.dendrogram(cluster.col), 37:1),
          cexRow=.5,
          cexCol=1,
          #labRow = FALSE,
          #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), 
          #lhei=c(.5,6,0),
          rowsep=0:nrow(t(rfe.set.mat)),
          colsep=0:ncol(t(rfe.set.mat)),
          sepwidth=c(0.00001,0.00001),
          sepcolor="black"
) 

##########################################################################

### K562

#### Importance barchart/heatmap

rfeModelResultsList <- readRDS(paste0("Z:/TAD_data_analysis/K562/5kb/results_by_chr/WHOLE/rus/distance/TADrfe_3.rds"))

set.size = 8
rfeModelResults <- rfeModelResultsList[[2]]
rfe.set <- rfeModelResults[rfeModelResults$Variables==set.size,]
rfe.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var), mean)
rfe.order <- order(rfe.set[, c("x")], decreasing = TRUE)
rfe.set <- rfe.set[rfe.order, ]
rfe.set$Group.1 <- factor(rfe.set$Group.1, levels=c("Smc3ab9263", 
                                                    "Znf143",
                                                    "Rad21",
                                                    "Ctcf",
                                                    "Pol2",
                                                    "Elf1sc631",  
                                                    "Pol24h8",    
                                                    "Tbp",        
                                                    "Taf1",      
                                                    "Mef2a",      
                                                    "Mxi1",       
                                                    "Mazab85725"))
ggplot(rfe.set, aes(x=reorder(Group.1,x), y=x)) +
  geom_bar(stat='identity', fill="red") +
  coord_flip() +
  xlab("Transcription Factor Binding Sites") +
  ylab("Predictive Importance") +
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


rfeModelResults <- rfeModelResultsList[[2]]
rfe.set <- rfeModelResults[rfeModelResults$Variables==set.size,]
rfe.set <- aggregate(rfe.set[, c("Overall")], list(rfe.set$var,rfe.set$Resample), mean)
rfe.set <- reshape(rfe.set, idvar = "Group.1", timevar = "Group.2", direction = "wide")
rfe.set.mat <- as.matrix(rfe.set[,-1])
rownames(rfe.set.mat) <- rfe.set$Group.1
rfe.set.mat[is.na(rfe.set.mat)] <- 0
rfe.set.mat <- rfe.set.mat[ order(rowMeans(rfe.set.mat)), ]

rfe.set.mat <- apply(rfe.set.mat,2,rev)

distance.col = dist(rfe.set.mat, method = "euclidean")
cluster.col = hclust(distance.col, method = "average")
my_palette <- colorRampPalette(c("white", "pink", "red"))(ncol(rfe.set.mat))

heatmap.2(t(rfe.set.mat),
          #labRow=rownames(rfe.set.mat),
          #labCol=colnames(rfe.set.mat),
          #cellnote = ,  # same data set for cell labels
          #main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins=c(12,8),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          na.color = "gray",
          #breaks=col_breaks,    # enable color transition at specified limits
          key=FALSE,
          #keysize = 1.5,
          #key.title = "Elastic-Net Coefficient",
          dendrogram="column",     # only draw a row dendrogram,
          Rowv = FALSE,
          Colv = reorder(as.dendrogram(cluster.col), 37:1),
          cexRow=.5,
          cexCol=1,
          #labRow = FALSE,
          #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), 
          #lhei=c(.5,6,0),
          rowsep=0:nrow(t(rfe.set.mat)),
          colsep=0:ncol(t(rfe.set.mat)),
          sepwidth=c(0.00001,0.00001),
          sepcolor="black"
) 

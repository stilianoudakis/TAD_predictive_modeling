#09 visualizing example of predicted TAD boundaries

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

## Extracting TAD boundaries using the extract_boundaries_func function for ARROWHEAD defined TADs at 10kb resolution

domains <- read.table("Z:/TAD_data_analysis/GM12878/5kb/GM12878_domain_data_5000.b.txt")
bounds.GR <- extract_boundaries_func(domains.mat=domains, 
                                     preprocess=FALSE, 
                                     CHR=paste0("CHR", c(1:8,10:22)), 
                                     resolution=5000)

## Obtaining list of GRanges objects of ChIP-seq data used for modelling

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS",
                                             pattern="*.bed",
                                             signal=4)

## running TADRF

tadModel <- TADRF(bounds.GR=bounds.GR,
                  resolution=5000,
                  genomicElements.GR=genomicElements.GR,
                  featureType="distance",
                  resampling="rus",
                  trainCHR=paste0("CHR", c(1:8,10:22)),
                  predictCHR=NULL,
                  cvFolds=3,
                  ntrees=100,
                  cvMetric="Accuracy",
                  verbose=TRUE,
                  seed=123,
                  model=TRUE,
                  importances=TRUE,
                  performances=FALSE)

bounds.GR <- extract_boundaries_func(domains.mat=domains, 
                                     preprocess=FALSE, 
                                     CHR="CHR22", 
                                     resolution=5000)

pt <- preciseTAD(bounds.GR=bounds.GR,
                      genomicElements.GR=genomicElements.GR,
                      featureType="distance",
                      CHR="CHR22",
                      chromCoords=list(17389000,18011000),
                      tadModel=tadModel[[1]],
                      threshold=1.0,
                      flank=NULL,
                      verbose=TRUE,
                      seed=123,
                      parallel=FALSE,
                      cores=NULL,
                      splits=NULL,
                      method.Dist="euclidean",
                      DBSCAN=TRUE,
                      DBSCAN_params=list(5000,3),
                      method.Clust=NULL,
                      PTBR=TRUE,
                      CLARA=TRUE,
                      samples=100)

getChrLength <- function(genome = "BSgenome.Hsapiens.UCSC.hg19"){
  g <- getBSgenome(genome, masked=FALSE)
  data.frame(chrom=1:24, length=seqlengths(g)[1:24])
}
.chrAsNum <- function(tbl){
  tbl$chrom <- gsub("chr", "", tbl$chrom)
  tbl$chrom[tbl$chrom=="X"] <- 23
  tbl$chrom[tbl$chrom=="Y"] <- 24
  tbl$chrom <- as.numeric(tbl$chrom)
  tbl[order(tbl$chrom),]
}
getCentromeres <- function( genome="hg19" ){
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case of failure, try another mirror
  if(inherits(mySession, "try-error"))
    mySession <- browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  obj <- ucscTableQuery(mySession, table="gap")
  tbl <- getTable(obj)
  tbl <- tbl[tbl$type=="centromere", c("chrom", "chromStart", "chromEnd")]
  colnames(tbl)[2:3] <- c("centromerStart", "centromerEnd")
  .chrAsNum(tbl)
}
makeHg19 <- function(){
  tbl <- merge(getChrLength(), getCentromeres(), by="chrom")
  cumlen <- c(0, cumsum(as.numeric(tbl$length))[-nrow(tbl)])
  cbind.data.frame(tbl, cumlen=cumlen)    
}
hg19 <- makeHg19()
hg19$chrom <- paste0("CHR",hg19$chrom)

seqDataTest <- c(17389000:18011000)

CHR="CHR22"
  
  if("TRUE" %in% table(seqDataTest %in% c(hg19$centromerStart[hg19$chrom==CHR]:hg19$centromerEnd[hg19$chrom==CHR]))){
    centromereTestStart <- hg19$centromerStart[hg19$chrom==CHR]
    centromereTestEnd <- hg19$centromerEnd[hg19$chrom==CHR]
    seqDataTest <- seqDataTest[-which(seqDataTest %in% c(centromereTestStart:centromereTestEnd))]
  }
  
test_data <- matrix(nrow=length(seqDataTest),
                    ncol=length(genomicElements.GR))
p <- lapply(genomicElements.GR, function(x){distance_func(GRanges(seqnames=tolower(CHR),
                                                                    IRanges(start=seqDataTest,
                                                                            end=seqDataTest)), 
                                                            x)})
for(i in 1:length(genomicElements.GR)){
  test_data[,i] <- p[[i]]
}
colnames(test_data) <- gsub(paste0(c("K562", "Gm12878", "-", "Broad", "Haib", "Sydh", "Uta", "Uw", "Uchicago"), collapse = "|"),"",names(genomicElements.GR))
test_data <- apply(test_data,2,function(x){log(x + 1, base=2)})
predictions <- predict(tadModel[[1]],newdata=test_data,type="prob")[,"Yes"]
															

predprobdata <- data.frame(basenum=c(17389000:18011000),
                           prob=predictions)
predprobdata$predicted <- ifelse(predprobdata$basenum %in% start(pt[[2]]), "Yes", "No")
predprobdata$called <- ifelse(predprobdata$basenum %in% start(pt[[3]]), "Yes", "No")

calleddata <- predprobdata[which(predprobdata$called=="Yes"),]

predicteddata <- predprobdata[which(predprobdata$prob>=1.0),]

predicteddata2 <- predprobdata[which(predprobdata$predicted=="Yes"),]
predicteddata2 <- predicteddata2[-c(3),]

ggplot() + 
  geom_line(data=predprobdata, aes(x=basenum, y=prob), 
            color='black',
            size=1)+
  #geom_vline(data=predicteddata, aes(xintercept=basenum, 
  #                                   color="red"),
  #           size=.5,
  #           show.legend = TRUE)+
  geom_vline(data=calleddata, aes(xintercept=basenum, 
                                  color="blue"),
             size=.75,
             show.legend = TRUE)+
  geom_vline(data=predicteddata2, aes(xintercept=basenum, 
                                      color="green"),
             size=.75,
             show.legend = TRUE)+
  xlab("Base Pair Coordinate") +
  ylab("Probability") +
  scale_color_manual(name = '', 
                     values =c("blue", "green"),
                     labels = c("ARROWHEAD",
                                "preciseTAD"))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")




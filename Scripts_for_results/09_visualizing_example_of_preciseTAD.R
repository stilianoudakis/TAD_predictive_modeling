#09_visualizing_example_of_preciseTAD

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

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS",
                                             pattern="*.bed",
                                             signal=4)

tadData <- createTADdata(bounds.GR=bounds.GR,
                         resolution=5000,
                         genomicElements.GR=genomicElements.GR,
                         featureType="distance",
                         resampling="rus",
                         trainCHR=paste0("CHR",c(1:8,10:22))[-which(paste0("CHR",c(1:8,10:22))=="CHR14")],
                         predictCHR="CHR14",
                         seed=123)

tadModel <- TADrandomForest(trainData=tadData[[1]],
                            testData=NULL,
                            tuneParams=list(mtry=2,
                                            ntree=500,
                                            nodesize=1),
                            cvFolds=3,
                            cvMetric="Accuracy",
                            verbose=TRUE,
                            seed=123,
                            model=TRUE,
                            importances=TRUE,
                            impMeasure="MDA",
                            performances=FALSE)

bounds.GR <- extract_boundaries_func(domains.mat=domains, 
                                     preprocess=FALSE, 
                                     CHR="CHR14", 
                                     resolution=5000)

pt <- preciseTAD(bounds.GR=bounds.GR,
                 genomicElements.GR=genomicElements.GR,
                 featureType="distance",
                 CHR="CHR14",
                 chromCoords=list(50085000,50800000),
                 tadModel=tadModel[[1]],
                 threshold=1,
                 flank=NULL,
                 verbose=TRUE,
                 seed=123,
                 parallel=FALSE,
                 cores=NULL,
                 splits=NULL,
                 method.Dist="euclidean",
                 DBSCAN=TRUE,
                 DBSCAN_params=list(10000,3),
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

#seqDataTest <- c(17389000:18011000)
#seqDataTest <- c(25745000:26176000)
#seqDataTest <- c(15150000:17270000)
seqDataTest <- c(50085000:50800000)

#CHR="CHR22"
#CHR="CHR21"
CHR="CHR14"

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

predprobdata <- data.frame(#basenum=c(17389000:18011000),
  #basenum=c(15150000:17270000),
  basenum=c(50085000:50800000),
  prob=predictions)

predprobdata$predicted <- ifelse(predprobdata$basenum %in% start(pt$PTBP), "Yes", "No")

predprobdata$called <- ifelse(predprobdata$basenum %in% start(pt$CTBP), "Yes", "No")

calleddata <- predprobdata[which(predprobdata$called=="Yes"),]

predicteddata <- predprobdata[which(predprobdata$prob>=1.0),]

predicteddata2 <- predprobdata[which(predprobdata$predicted=="Yes"),]

ptbrdf <- data.frame(x=c(1,2),y=c(50315887,50334548))

ggplot() + 
  #probabilities
  geom_line(data=predprobdata, aes(x=basenum, y=prob), 
            color='black',
            size=1)+
  #ctbp
  geom_vline(data=calleddata, aes(xintercept=basenum, 
                                  color="blue"),
             size=1.5,
             show.legend = TRUE)+
  #ptbp
  geom_vline(data=predicteddata2, aes(xintercept=basenum, 
                                      color="forestgreen"),
             size=1.5,
             show.legend = TRUE)+
  #ptbr
  geom_vline(data=ptbrdf, aes(xintercept=y, 
                                      color="yellow"),
             size=1,
             show.legend = TRUE)+
  xlab("Base Coordinate (mb)") +
  ylab("") +
  scale_color_manual(name = '', 
                     values =c("blue", "forestgreen", "yellow"),
                     labels = c("ARROWHEAD",
                                "preciseTAD",
                                "PTBR"))+
  scale_y_continuous(breaks=c(0,.5,1), limits=c(0,1))+
  scale_x_continuous(breaks=c(50200000,50400000,50600000,50800000), 
                     labels=c("50.2","50.4","50.6","50.8"))+
  guides(color=guide_legend(nrow=3,byrow=TRUE)) +
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




ggplot() + 
  #probabilities
  geom_line(data=predprobdata, aes(x=basenum, y=prob), 
            color='black',
            size=1)+
  #ctbp
  geom_vline(data=calleddata, aes(xintercept=basenum, 
                                  color="blue"),
             size=1.5,
             show.legend = TRUE)+
  #ptbp
  geom_vline(data=predicteddata2, aes(xintercept=basenum, 
                                      color="forestgreen"),
             size=1.5,
             show.legend = TRUE)+
  #ptbr
  geom_vline(data=ptbrdf, aes(xintercept=y, 
                              color="yellow"),
             size=1,
             show.legend = TRUE)+
  xlab("Base Coordinate (mb)") +
  ylab("") +
  xlim(50314887,50335548) +
  scale_color_manual(name = '', 
                     values =c("blue", "forestgreen", "yellow"),
                     labels = c("ARROWHEAD",
                                "preciseTAD",
                                "PTBR"))+
  scale_y_continuous(breaks=c(0,.5,1), limits=c(0,1))+
  #scale_x_continuous(breaks=c(50200000,50400000,50600000,50800000), 
  #                   labels=c("50.2","50.4","50.6","50.8"))+
  guides(color=guide_legend(nrow=3,byrow=TRUE)) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.ticks.y = element_blank(),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "none")


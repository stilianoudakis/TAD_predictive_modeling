#Model filtering

#Loading Libraries

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
library(leaps)

# Creating Response Vector Y

## Reading in TAD data

domains <- read.table("/home/stilianoudakisc/TAD_data_analysis/GM12878/arrowhead_data.txt", header=T)
domains <- domains[,1:3]
head(domains)
dim(domains) #9274    3

#removing chrX from boundary data
domains <- domains[-which(domains$Chromosome=="chrX"),]

#creating granges object out of tad boundary data

#concatenating boundary coordinates into one long vector
coords <- domains
colnames(coords)[2:3] <- c("coordinate", "coordinate")
coords <- rbind.data.frame(coords[,c(1,2)],coords[,c(1,3)])
coords$Chromosome <- as.character(coords$Chromosome)

#sorting by chromosome and coordinate
coords <- coords[order(as.numeric(substr(coords$Chromosome,4,5)), coords$coordinate, decreasing = FALSE),]

#remove duplicates for coordinates that are conjoined
coords <- coords[!duplicated(coords),]
dim(coords) #16770     2

coords$lagchrom <- lag(coords$Chromosome)
coords$lagchrom[which(coords$Chromosome!=coords$lagchrom)] <- NA
coords$lagcoord <- lag(coords$coordinate)
coords$lagcoord[which(is.na(coords$lagchrom))] <- NA

bounds <- GRanges(seqnames=coords$Chromosome, ranges=IRanges(start=coords$coordinate, width=1))
bounds_f <- flank(bounds, width=5000, both=TRUE)

#saving data
saveRDS(bounds, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/bounds.rds")
saveRDS(bounds_f, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/bounds_f.rds")


## Reading in binned genome in the form of contact matrix at 10kb resolution

binslist5 <- read.table("/home/stilianoudakisc/TAD_data_analysis/GSE63525_data/contact_matrices/binned_data/5kb/binslist5.txt")

dim(binslist5) #536396      2

#ordering the bins according to left endpoint
binslist5 <- binslist5[order(as.numeric(substr(binslist5$V1,4,5)), binslist5$V2, decreasing=FALSE),]

#extracting and renaming first 2 columns
binslist5 <- binslist5[,1:2]
colnames(binslist5) <- c("Chromosome", "Coordinate")


#creating a granges object from binned genome
binslist5 <- GRanges(seqnames = binslist5$Chromosome, ranges = IRanges(start = binslist5$Coordinate,
                                                                         width = 5000))
																		 

#saving the binned genome data
saveRDS(binslist5, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/binslist5.rds")


## Creating Response Vector Y (1 if tad boundary is in bin; 0 if not)

#Finding where the TAD boundaries are inside the genomic bins
findOverlaps(binslist5, bounds_f)

y <- countOverlaps(binslist5, bounds_f)
length(y) #536396
table(y)
#     0      1      2
#504576  30126   1694



#investigating where there are multiple overlaps
head(which(y==2))  #324 1456 1479 1516 1604 1698

#looking at row 169 of binned genome
findOverlaps(binslist5[324], bounds_f)
#      queryHits subjectHits
#      <integer>   <integer>
#  [1]         1           9
#  [2]         1          10

bounds_f[9:10]
#      seqnames             ranges strand
#         <Rle>          <IRanges>  <Rle>
#  [1]     chr1 [1855000, 1864999]      *
#  [2]     chr1 [1860000, 1869999]      *


#recode multiple overlaps as 1
y <- ifelse(y>=1,1,0)
prop.table(table(y))  
#         0          1
#0.94067816 0.05932184

mcols(binslist5)$y <- y


## Creating the data frame with resonse vector for modeling

gm12878_5kb <- data.frame(y = y)


saveRDS(gm12878_5kb, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/gm12878_5kb.rds")

#creating a granges object from the center of every range in the binned genome granges object
#this will be used to calculate distances from genomic feature to center of genomic bin
binslist5_center <- resize(binslist5, width = 1, fix = "center")

saveRDS(binslist5_center, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/binslist5_center.rds")

# Creating Feature Space

## 3D subcompartments

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/3D_subcompartments/")

temp = list.files()

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}


##DGV

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DGV/")

temp = list.files()

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}
  
## GERP

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/GERP/")

temp = list.files()

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}


## super_enhancers

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/super_enhancers/")

temp = list.files()
temp <- temp[grep("GM12878",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}

## VMR

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/VMRs/")

temp = list.files()

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}
  
## BroadHMM

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/BroadHMM/")

temp <- list.files()
temp <- temp[grep("Gm12878",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}
  
## Combined

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/Combined/")

temp <- list.files()
temp <- temp[grep("Gm12878",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}

## DNase I


setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DNaseI/")

temp <- list.files()
temp <- temp[grep("Gm12878",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}

## Histone Modifications

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/HistoneModifications/")

temp <- list.files()
temp <- temp[grep("Gm12878",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}

## Adding Chromosome information to the data

#gm12878_5kb$CHR <- seqnames(binslist5)

#gm12878_5kb$CHR <- as.character(gm12878_5kb$CHR)


## Saving the data

dim(gm12878_5kb)

head(gm12878_5kb)

saveRDS(gm12878_5kb, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/gm12878_5kb.rds")

## Distance Summaries for each feature from region to TAD boundary

#grobjects <- ls()
#grobjects <- grobjects[grep("center", grobjects)]

#grlist <- rep( list(GRangesList()), length(grobjects) )

#for(i in 1:length(grobjects)){
  
#  x <- get(grobjects[i])
  
#  grlist[[i]] <- x
  
  #grlist[[i]] <- grlist
  
#}

#grlist <- setNames(grlist, grobjects)

#grlist <- grlist[-which(names(grlist)=="binslist5_center")]

#length(grlist)

#for(i in 1:length(x)){
#  
#  if(is.element("chrX", as.character(seqnames(x[[i]])))){
#    x[[i]] <- x[[i]][-which(as.character(seqnames(x[[i]]))=="chrX")]
#  }
#  
#  if(is.element("chrY", as.character(seqnames(x[[i]])))){
#    x[[i]] <- x[[i]][-which(as.character(seqnames(x[[i]]))=="chrY")]
#  }
#  
#  mcols(x[[i]])$distance <- mcols(distanceToNearest(x[[i]], bounds))$distance
#  mcols(x[[i]])$logdistance <- log(mcols(distanceToNearest(x[[i]], bounds))$distance+1,
#                                          base = 2)
#  
#}


#saveRDS(grlist, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/grlist.rds")

#Summaries

#meandist <- unlist(lapply(grlist, function(x){mean(mcols(x)$distance)}))
#mediandist <- unlist(lapply(grlist, function(x){median(mcols(x)$distance)}))
#rangedist <- unlist(lapply(grlist, function(x){max(mcols(x)$distance)})) - unlist(lapply(grlist, function(x){min(mcols(x)$distance)}))

#meanlogdist <- unlist(lapply(grlist, function(x){mean(mcols(x)$logdistance)}))
#medianlogdist <- unlist(lapply(grlist, function(x){median(mcols(x)$logdistance)}))
#rangelogdist <- unlist(lapply(grlist, function(x){max(mcols(x)$logdistance, base = 2)})) - unlist(lapply(grlist, function(x){min(mcols(x)$logdistance, base = 2)}))

#distancetab <- data.frame(Feature = names(unlist(lapply(grlist, function(x){mean(mcols(x)$distance)}))),
#                          Mean = meandist,
#                          Median = mediandist,
#                          Range = rangedist,
#                          MeanLog = meanlogdist,
#                          MedianLog = medianlogdist,
#                          RangeLog = rangelogdist)

#rownames(distancetab) <- NULL

#distancetab[,2:7] <- round(distancetab[,2:7], 1)

#saveRDS(distancetab, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/distancetab.rds")


#Plots
#d <- density(mcols(grlist[[1]])$distance)
#png(filename="/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/distance_plot.png")
#plot(d, xlab="Distance", 
#     main="Distance from Region to TAD Boundary", 
#     xlim=c(0,1000000),
#     ylim=c(0,3.5e-05))
#for(i in 1:(length(grlist)-1)){
#  lines(density(mcols(grlist[[i+1]])$distance), col=i+1)
#}
#dev.off()


#d <- density(mcols(grlist[[1]])$logdistance)
#png(filename="/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/log_distance_plot.png")
#plot(d, xlab="Distance", 
#     main="Log Distance from Region to TAD Boundary")
#for(i in 1:(length(grlist)-1)){
#  lines(density(mcols(grlist[[i+1]])$logdistance), col=i+1)
#}
#dev.off()


# Model Filtering

##Taking log2 transform of distance variables

cols <- grep("dist",colnames(gm12878_5kb))
gm12878_5kb[,cols] <- apply(gm12878_5kb[,cols],2,function(x){log(x + 1, base=2)})


##Changing binary variables to factors

#cols <- c(intersect(grep("score",colnames(gm12878_5kb), invert = TRUE),
#          grep("dist",colnames(gm12878_5kb), invert = TRUE)))
gm12878_5kb$y <- factor(gm12878_5kb$y)


##Changing levels of response (y) to yes no

levels(gm12878_5kb$y) <- c("No", "Yes")

##Removing zero variance predictors

nzv <- nearZeroVar(gm12878_5kb[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])

nzvar

gm12878_5kb_f <- gm12878_5kb[, -which(colnames(gm12878_5kb) %in% nzvar)]

dim(gm12878_5kb) 
dim(gm12878_5kb_f) 

saveRDS(gm12878_5kb_f, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/gm12878_5kb_f.rds")


setwd("/home/stilianoudakisc/TAD_data_analysis/model_filtering/")

cols <- names(gm12878_5kb_f)[grep("dist", names(gm12878_5kb_f))]
featdistdf <- gm12878_5kb_f[,c(1, which(names(gm12878_5kb_f) %in% cols))]

set.seed(123)
sampids <- sample(which(featdistdf$y=="No"),
                        length(which(featdistdf$y=="Yes")),
                        replace = FALSE)  
  
featdistdf2 <- rbind.data.frame(featdistdf[which(featdistdf$y=="Yes"),],
                           featdistdf[sampids,])
						   
saveRDS(featdistdf2,"/home/stilianoudakisc/TAD_data_analysis/model_filtering/featdistdf2.rds")						   
						   
for(i in 1:(dim(featdistdf2)[2]-1)){
  p <- ggplot(featdistdf2,aes(x=featdistdf2[,i+1],group=y,fill=y))+
	labs(x = colnames(featdistdf2)[i+1], y = "Count") +
	scale_fill_manual(name="TAD in Bin", values=c("red","blue")) +
    geom_histogram(position="identity",alpha=0.5,binwidth=0.25)+
    theme_bw()
  ggsave(p, file=paste0("plot_", colnames(featdistdf2)[i+1],".png"))
}


disttest <- lapply(featdistdf2[,-1], function(i) t.test(i ~ featdistdf2$y))
distancetab <- matrix(nrow=(dim(featdistdf2)[2]-1),ncol=3)
rownames(distancetab) <- names(featdistdf2)[-1]
for(i in 1:(dim(featdistdf2)[2]-1)){
	distancetab[i,1] <- as.numeric(disttest[[i]][[5]][1])
	distancetab[i,2] <- as.numeric(disttest[[i]][[5]][2])
	distancetab[i,3] <- disttest[[i]][[3]]
}
colnames(distancetab) <- c("TAD not in Bin", "TAD in Bin", "P-Value")

distancetab <- as.data.frame(distancetab)

saveRDS(distancetab, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/distancetab.rds")





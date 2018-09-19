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

domains <- read.table("/home/stilianoudakisc/TAD_data_analysis/K562/arrowhead_k562_data.txt", header=T)
head(domains)
dim(domains) # 5975    3

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
dim(coords) #10293     2

coords$lagchrom <- lag(coords$Chromosome)
coords$lagchrom[which(coords$Chromosome!=coords$lagchrom)] <- NA
coords$lagcoord <- lag(coords$coordinate)
coords$lagcoord[which(is.na(coords$lagchrom))] <- NA

bounds <- GRanges(seqnames=coords$Chromosome, ranges=IRanges(start=coords$coordinate, width=1))
bounds_f <- flank(bounds, width=10000, both=TRUE)

#saving data
saveRDS(bounds, "/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/bounds.rds")
saveRDS(bounds_f, "/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/bounds_f.rds")


## Reading in binned genome in the form of contact matrix at 10kb resolution

binslist10 <- read.table("/home/stilianoudakisc/TAD_data_analysis/K562/binslist10.txt")
dim(binslist10) #268097      2
head(binslist10)


#ordering the bins according to left endpoint
binslist10 <- binslist10[order(as.numeric(substr(binslist10$V1,4,5)), binslist10$V2, decreasing=FALSE),]
head(binslist10)

#extracting and renaming first 2 columns
colnames(binslist10) <- c("Chromosome", "Coordinate")


#creating a granges object from binned genome
binslist10 <- GRanges(seqnames = binslist10$Chromosome, ranges = IRanges(start = binslist10$Coordinate,
                                                                         width = 10000))
																		 

#saving the binned genome data
saveRDS(binslist10, "/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/binslist10.rds")


## Creating Response Vector Y (1 if tad boundary is in bin; 0 if not)

#Finding where the TAD boundaries are inside the genomic bins
findOverlaps(binslist10, bounds_f)

y <- countOverlaps(binslist10, bounds_f)
length(y) #268097
table(y)
#     0      1      2
#249016  17593   1488



#investigating where there are multiple overlaps
head(which(y==2))  #107 215 314 638 736 737

#looking at row 169 of binned genome
findOverlaps(binslist10[107], bounds_f)
#      queryHits subjectHits
#      <integer>   <integer>
#  [1]         1           4
#  [2]         1           5


bounds_f[4:5]
#      seqnames             ranges strand
#         <Rle>          <IRanges>  <Rle>
#  [1]     chr1 [1230000, 1249999]      *
#  [2]     chr1 [1240000, 1259999]      *



#recode multiple overlaps as 1
y <- ifelse(y>=1,1,0)
prop.table(table(y))  
#       0        1
#0.928828 0.071172

mcols(binslist10)$y <- y


## Creating the data frame with resonse vector for modeling

k562_10kb <- data.frame(y = y)


saveRDS(k562_10kb, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/k562_10kb.rds")
dim(k562_10kb) #268097      1

#creating a granges object from the center of every range in the binned genome granges object
#this will be used to calculate distances from genomic feature to center of genomic bin
binslist10_center <- resize(binslist10, width = 1, fix = "center")

saveRDS(binslist10_center, "/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/binslist10_center.rds")


# Creating Feature Space


##DGV

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DGV/")

temp = list.files()

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist10_center, dat_gr_center))$distance
	
	k562_10kb <- cbind.data.frame(k562_10kb,d)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist10, dat_gr)
	
	k562_10kb <- cbind.data.frame(k562_10kb,c)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
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
	
	d <- mcols(distanceToNearest(binslist10_center, dat_gr_center))$distance
	
	k562_10kb <- cbind.data.frame(k562_10kb,d)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist10, dat_gr)
	
	k562_10kb <- cbind.data.frame(k562_10kb,c)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}


## super_enhancers

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/super_enhancers/")

temp = list.files()
temp <- temp[grep("K562",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist10_center, dat_gr_center))$distance
	
	k562_10kb <- cbind.data.frame(k562_10kb,d)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist10, dat_gr)
	
	k562_10kb <- cbind.data.frame(k562_10kb,c)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
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
	
	d <- mcols(distanceToNearest(binslist10_center, dat_gr_center))$distance
	
	k562_10kb <- cbind.data.frame(k562_10kb,d)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist10, dat_gr)
	
	k562_10kb <- cbind.data.frame(k562_10kb,c)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}
  
## BroadHMM

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/BroadHMM/")

temp <- list.files()
temp <- temp[grep("K562",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist10_center, dat_gr_center))$distance
	
	k562_10kb <- cbind.data.frame(k562_10kb,d)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist10, dat_gr)
	
	k562_10kb <- cbind.data.frame(k562_10kb,c)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}
  
## Combined

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/Combined/")

temp <- list.files()
temp <- temp[grep("K562",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist10_center, dat_gr_center))$distance
	
	k562_10kb <- cbind.data.frame(k562_10kb,d)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist10, dat_gr)
	
	k562_10kb <- cbind.data.frame(k562_10kb,c)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}

## DNase I


setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DNaseI/")

temp <- list.files()
temp <- temp[grep("K562",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist10_center, dat_gr_center))$distance
	
	k562_10kb <- cbind.data.frame(k562_10kb,d)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist10, dat_gr)
	
	k562_10kb <- cbind.data.frame(k562_10kb,c)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}

## Histone Modifications

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/HistoneModifications/")

temp <- list.files()
temp <- temp[grep("K562",temp)]

for(i in 1:length(temp)){
	dat <- read.table(temp[i],header = FALSE, sep="\t")
	dat_gr <- GRanges(seqnames=dat$V1,
						IRanges(start=dat$V2, 
								end=dat$V3))
	dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
	
	d <- mcols(distanceToNearest(binslist10_center, dat_gr_center))$distance
	
	k562_10kb <- cbind.data.frame(k562_10kb,d)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist10, dat_gr)
	
	k562_10kb <- cbind.data.frame(k562_10kb,c)
	names(k562_10kb)[length(names(k562_10kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}

## Adding Chromosome information to the data

#k562_10kb$CHR <- seqnames(binslist10)

#k562_10kb$CHR <- as.character(k562_10kb$CHR)


## Saving the data

dim(k562_10kb)

head(k562_10kb)

saveRDS(k562_10kb, "/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/k562_10kb.rds")


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

cols <- grep("dist",colnames(k562_10kb))
k562_10kb[,cols] <- apply(k562_10kb[,cols],2,function(x){log(x + 1, base=2)})


##Changing binary variables to factors

k562_10kb$y <- factor(k562_10kb$y)


##Changing levels of response (y) to yes no

levels(k562_10kb$y) <- c("No", "Yes")

##Removing zero variance predictors

nzv <- nearZeroVar(k562_10kb[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])

nzvar

k562_10kb_f <- k562_10kb[, -which(colnames(k562_10kb) %in% nzvar)]

dim(k562_10kb) 
dim(k562_10kb_f) 

saveRDS(k562_10kb_f, "/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/k562_10kb_f.rds")


setwd("/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/")

cols <- names(k562_10kb_f)[grep("dist", names(k562_10kb_f))]
featdistdf <- k562_10kb_f[,c(1, which(names(k562_10kb_f) %in% cols))]

set.seed(123)
sampids <- sample(which(featdistdf$y=="No"),
                        length(which(featdistdf$y=="Yes")),
                        replace = FALSE)  
  
featdistdf2 <- rbind.data.frame(featdistdf[which(featdistdf$y=="Yes"),],
                           featdistdf[sampids,])
						   
saveRDS(featdistdf2,"/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/featdistdf2.rds")						   
						   
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

saveRDS(distancetab, "/home/stilianoudakisc/TAD_data_analysis/K562/model_filtering/distancetab.rds")



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
#saveRDS(bounds, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/bounds.rds")
#saveRDS(bounds_f, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/bounds_f.rds")


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
#saveRDS(binslist5, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/binslist5.rds")


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

###########################################################################################

binslist5 <- binslist5[seqnames(binslist5)=="chr22"]
y <- countOverlaps(binslist5, bounds_f)
y <- ifelse(y>=1,1,0)
mcols(binslist5)$y <- y


## Creating the data frame with resonse vector for modeling

gm12878_5kb <- data.frame(y = y)


#saveRDS(gm12878_5kb, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/gm12878_5kb.rds")

#creating a granges object from the center of every range in the binned genome granges object
#this will be used to calculate distances from genomic feature to center of genomic bin
binslist5_center <- resize(binslist5, width = 1, fix = "center")

#saveRDS(binslist5_center, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/binslist5_center.rds")


# Creating Feature Space


  
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
	
	w <- unlist(lapply(binslist5, function(x){
							sum(width(
									ranges(
										pintersect(
											findOverlapPairs(x,dat_gr)))))}))

	gm12878_5kb <- cbind.data.frame(gm12878_5kb,w)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- sub(".bed", "", temp[i])
	
	d <- mcols(distanceToNearest(binslist5_center, dat_gr_center))$distance
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,d)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_dist", sep="")

	c <- countOverlaps(binslist5, dat_gr)
	
	gm12878_5kb <- cbind.data.frame(gm12878_5kb,c)
	names(gm12878_5kb)[length(names(gm12878_5kb))] <- paste(sub(".bed", "", temp[i]), "_count", sep="")
	
	}



## Saving the data

dim(gm12878_5kb)

head(gm12878_5kb)

combined <- gm12878_5kb

saveRDS(combined, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/combined.rds")






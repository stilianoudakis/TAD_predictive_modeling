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

gm12878_10kb <- data.frame(y = y)


saveRDS(gm12878_10kb, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/gm12878_10kb.rds")

#creating a granges object from the center of every range in the binned genome granges object
#this will be used to calculate distances from genomic feature to center of genomic bin
binslist5_center <- resize(binslist5, width = 1, fix = "center")

saveRDS(binslist5_center, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/binslist5_center.rds")

# Creating Feature Space

## 3D subcompartments

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/3D_subcompartments/")

subcomp_A <- read.table("GSE63525_GM12878_subcompartments_A.bed",header = FALSE, sep="\t")
subcomp_B <- read.table("GSE63525_GM12878_subcompartments_B.bed",header = FALSE, sep="\t")

A_gr <- GRanges(seqnames=subcomp_A$V1,IRanges(start=subcomp_A$V2, end=subcomp_A$V3))
A_gr_center <- resize(A_gr, width = 1, fix = "center")
B_gr <- GRanges(seqnames=subcomp_B$V1,IRanges(start=subcomp_B$V2, end=subcomp_B$V3))
B_gr_center <- resize(B_gr, width = 1, fix = "center")

  gm12878_10kb$A <- NA
  gm12878_10kb$A[queryHits(findOverlaps(binslist5,A_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,A_gr))))
  gm12878_10kb$A[-queryHits(findOverlaps(binslist5,A_gr))] <- 0
  
  gm12878_10kb$B <- NA
  gm12878_10kb$B[queryHits(findOverlaps(binslist5,B_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,B_gr))))
  gm12878_10kb$B[-queryHits(findOverlaps(binslist5,B_gr))] <- 0
  
  gm12878_10kb$A_dist <- mcols(distanceToNearest(binslist5_center, A_gr_center))$distance
  gm12878_10kb$B_dist <- mcols(distanceToNearest(binslist5_center, B_gr_center))$distance
  
  gm12878_10kb$Gm12878_A_count <- countOverlaps(binslist5, A_gr)
  gm12878_10kb$Gm12878_B_count <- countOverlaps(binslist5, B_gr)
  
##DGV

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DGV/")

temp = list.files()

complex <- read.table(temp[1],header=FALSE,sep="\t")
deletion <- read.table(temp[2],header=FALSE,sep="\t")
duplication <- read.table(temp[3],header=FALSE,sep="\t")
gain_loss <- read.table(temp[4],header=FALSE,sep="\t")
insertion <- read.table(temp[5],header=FALSE,sep="\t")
inversion <- read.table(temp[6],header=FALSE,sep="\t")
mobile_element_insertion <- read.table(temp[7],header=FALSE,sep="\t")
novel_sequence_insertion <- read.table(temp[8],header=FALSE,sep="\t")
sequence_alteration <- read.table(temp[9],header=FALSE,sep="\t")

complex_gr <- GRanges(seqnames=complex$V1,IRanges(start=complex$V2,end=complex$V3))
complex_gr_center <- resize(complex_gr, width = 1, fix = "center")
deletion_gr <- GRanges(seqnames=deletion$V1,IRanges(start=deletion$V2,end=deletion$V3))
deletion_gr_center <- resize(deletion_gr, width = 1, fix = "center")
duplication_gr <- GRanges(seqnames=duplication$V1,IRanges(start=duplication$V2,end=duplication$V3))
duplication_gr_center <- resize(duplication_gr, width = 1, fix = "center")
gain_loss_gr <- GRanges(seqnames=gain_loss$V1,IRanges(start=gain_loss$V2,end=gain_loss$V3))
gain_loss_gr_center <- resize(gain_loss_gr, width = 1, fix = "center")
insertion_gr <- GRanges(seqnames=insertion$V1,IRanges(start=insertion$V2,end=insertion$V3))
insertion_gr_center <- resize(insertion_gr, width = 1, fix = "center")
inversion_gr <- GRanges(seqnames=inversion$V1,IRanges(start=inversion$V2,end=inversion$V3))
inversion_gr_center <- resize(inversion_gr, width = 1, fix = "center")
mobile_element_insertion_gr <- GRanges(seqnames=mobile_element_insertion$V1,IRanges(start=mobile_element_insertion$V2,end=mobile_element_insertion$V3))
mobile_element_insertion_gr_center <- resize(mobile_element_insertion_gr, width = 1, fix = "center")
novel_sequence_insertion_gr <- GRanges(seqnames=novel_sequence_insertion$V1,IRanges(start=novel_sequence_insertion$V2,end=novel_sequence_insertion$V3))
novel_sequence_insertion_gr_center <- resize(novel_sequence_insertion_gr, width = 1, fix = "center")
sequence_alteration_gr <- GRanges(seqnames=sequence_alteration$V1,IRanges(start=sequence_alteration$V2,end=sequence_alteration$V3))
sequence_alteration_gr_center <- resize(sequence_alteration_gr, width = 1, fix = "center")

  gm12878_10kb$complex <- NA
  gm12878_10kb$complex[queryHits(findOverlaps(binslist5,complex_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,complex_gr))))
  gm12878_10kb$complex[-queryHits(findOverlaps(binslist5,complex_gr))] <- 0
  
  gm12878_10kb$deletion <- NA
  gm12878_10kb$deletion[queryHits(findOverlaps(binslist5,deletion_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,deletion_gr))))
  gm12878_10kb$deletion[-queryHits(findOverlaps(binslist5,deletion_gr))] <- 0
  
  gm12878_10kb$duplication <- NA
  gm12878_10kb$duplication[queryHits(findOverlaps(binslist5,duplication_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,duplication_gr))))
  gm12878_10kb$duplication[-queryHits(findOverlaps(binslist5,duplication_gr))] <- 0

  gm12878_10kb$gain_loss <- NA
  gm12878_10kb$gain_loss[queryHits(findOverlaps(binslist5,gain_loss_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,gain_loss_gr))))
  gm12878_10kb$gain_loss[-queryHits(findOverlaps(binslist5,gain_loss_gr))] <- 0

  gm12878_10kb$insertion <- NA
  gm12878_10kb$insertion[queryHits(findOverlaps(binslist5,insertion_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,insertion_gr))))
  gm12878_10kb$insertion[-queryHits(findOverlaps(binslist5,insertion_gr))] <- 0

  gm12878_10kb$inversion <- NA
  gm12878_10kb$inversion[queryHits(findOverlaps(binslist5,inversion_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,inversion_gr))))
  gm12878_10kb$inversion[-queryHits(findOverlaps(binslist5,inversion_gr))] <- 0

  gm12878_10kb$mobile_element_insertion <- NA
  gm12878_10kb$mobile_element_insertion[queryHits(findOverlaps(binslist5,mobile_element_insertion_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,mobile_element_insertion_gr))))
  gm12878_10kb$mobile_element_insertion[-queryHits(findOverlaps(binslist5,mobile_element_insertion_gr))] <- 0

  gm12878_10kb$novel_sequence_insertion <- NA
  gm12878_10kb$novel_sequence_insertion[queryHits(findOverlaps(binslist5,novel_sequence_insertion_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,novel_sequence_insertion_gr))))
  gm12878_10kb$novel_sequence_insertion[-queryHits(findOverlaps(binslist5,novel_sequence_insertion_gr))] <- 0
  
  gm12878_10kb$sequence_alteration <- NA
  gm12878_10kb$sequence_alteration[queryHits(findOverlaps(binslist5,sequence_alteration_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,sequence_alteration_gr))))
  gm12878_10kb$sequence_alteration[-queryHits(findOverlaps(binslist5,sequence_alteration_gr))] <- 0
  
  
  gm12878_10kb$complex_dist <- mcols(distanceToNearest(binslist5_center, complex_gr_center))$distance
  gm12878_10kb$deletion_dist <- mcols(distanceToNearest(binslist5_center, deletion_gr_center))$distance
  gm12878_10kb$duplication_dist <- mcols(distanceToNearest(binslist5_center, duplication_gr_center))$distance
  gm12878_10kb$gain_loss_dist <- mcols(distanceToNearest(binslist5_center, gain_loss_gr_center))$distance
  gm12878_10kb$insertion_dist <- mcols(distanceToNearest(binslist5_center, insertion_gr_center))$distance
  gm12878_10kb$inversion_dist <- mcols(distanceToNearest(binslist5_center, inversion_gr_center))$distance
  gm12878_10kb$mobile_element_insertion_dist <- mcols(distanceToNearest(binslist5_center, mobile_element_insertion_gr_center))$distance
  gm12878_10kb$novel_sequence_insertion_dist <- mcols(distanceToNearest(binslist5_center, novel_sequence_insertion_gr_center))$distance
  gm12878_10kb$sequence_alteration_dist <- mcols(distanceToNearest(binslist5_center, sequence_alteration_gr_center))$distance
  
  gm12878_10kb$complex_count <- countOverlaps(binslist5, complex_gr)
  gm12878_10kb$deletion_count <- countOverlaps(binslist5, deletion_gr)
  gm12878_10kb$duplication_count <- countOverlaps(binslist5, duplication_gr)
  gm12878_10kb$gain_loss_count <- countOverlaps(binslist5, gain_loss_gr)
  gm12878_10kb$insertion_count <- countOverlaps(binslist5, insertion_gr)
  gm12878_10kb$inversion_count <- countOverlaps(binslist5, inversion_gr)
  gm12878_10kb$mobile_element_insertion_count <- countOverlaps(binslist5, mobile_element_insertion_gr)
  gm12878_10kb$novel_sequence_insertion_count <- countOverlaps(binslist5, novel_sequence_insertion_gr)
  gm12878_10kb$sequence_alteration_count <- countOverlaps(binslist5, sequence_alteration_gr)
  
  
## GERP

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/GERP/")

gerp <- read.table("GERP_hg19.bed",header=FALSE,sep="\t")

gerp_gr <- GRanges(seqnames=gerp$V1,IRanges(start=gerp$V2,end=gerp$V3))
mcols(gerp_gr)$score <- gerp$V5
gerp_gr_center <- resize(gerp_gr, width = 1, fix = "center")

  gm12878_10kb$gerp <- NA
  gm12878_10kb$gerp[queryHits(findOverlaps(binslist5,gerp_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,gerp_gr))))
  gm12878_10kb$gerp[-queryHits(findOverlaps(binslist5,gerp_gr))] <- 0
  
  gm12878_10kb$gerp_dist <- mcols(distanceToNearest(binslist5_center, gerp_gr_center))$distance
  
  gm12878_10kb$gerp_count <- countOverlaps(binslist5, gerp_gr)
  
#finding which flanks overlap the gerp file so that we can add a score variable
#all other flanks will have a score of 0
#which(gm12878_10kb$gerp==1)
gm12878_10kb$gerp_score <- 0
gerpoverlap <- findOverlaps(binslist5,gerp_gr)
gerpoverlapdf <- data.frame(queryHits=queryHits(gerpoverlap), score=gerp_gr[subjectHits(gerpoverlap)]$score)
gerpoverlapmean <- aggregate(gerpoverlapdf$score, list(gerpoverlapdf$queryHits), mean)
gm12878_10kb$gerp_score[gerpoverlapmean$Group.1] <- gerpoverlapmean$x


## nestedRepeats

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/nestedRepeats/")

temp = list.files()

DNA <- read.table(temp[1],header=FALSE,sep="\t")
line <- read.table(temp[2],header=FALSE,sep="\t")
low_complexity <- read.table(temp[3],header=FALSE,sep="\t")
LTR  <- read.table(temp[4],header=FALSE,sep="\t")
other <- read.table(temp[5],header=FALSE,sep="\t")
RC <- read.table(temp[6],header=FALSE,sep="\t")
satellite <- read.table(temp[9],header=FALSE,sep="\t")
simple_repeat <- read.table(temp[11],header=FALSE,sep="\t")
SINE <- read.table(temp[12],header=FALSE,sep="\t")

DNA_gr <- GRanges(seqnames=DNA$V1,IRanges(start=DNA$V2,end=DNA$V3))
DNA_gr_center <- resize(DNA_gr, width = 1, fix = "center")
line_gr <- GRanges(seqnames=line$V1,IRanges(start=line$V2,end=line$V3))
line_gr_center <- resize(line_gr, width = 1, fix = "center")
low_complexity_gr <- GRanges(seqnames=low_complexity$V1,IRanges(start=low_complexity$V2,end=low_complexity$V3))
low_complexity_gr_center <- resize(low_complexity_gr, width = 1, fix = "center")
LTR_gr <- GRanges(seqnames=LTR$V1,IRanges(start=LTR$V2,end=LTR$V3))
LTR_gr_center <- resize(LTR_gr, width = 1, fix = "center")
other_gr <- GRanges(seqnames=other$V1,IRanges(start=other$V2,end=other$V3))
other_gr_center <- resize(other_gr, width = 1, fix = "center")
RC_gr <- GRanges(seqnames=RC$V1,IRanges(start=RC$V2,end=RC$V3))
RC_gr_center <- resize(RC_gr, width = 1, fix = "center")
satellite_gr <- GRanges(seqnames=satellite$V1,IRanges(start=satellite$V2,end=satellite$V3))
satellite_gr_center <- resize(satellite_gr, width = 1, fix = "center")
simple_repeat_gr <- GRanges(seqnames=simple_repeat$V1,IRanges(start=simple_repeat$V2,end=simple_repeat$V3))
simple_repeat_gr_center <- resize(simple_repeat_gr, width = 1, fix = "center")
SINE_gr <- GRanges(seqnames=SINE$V1,IRanges(start=SINE$V2,end=SINE$V3))
SINE_gr_center <- resize(SINE_gr, width = 1, fix = "center")

  gm12878_10kb$DNA <- NA
  gm12878_10kb$DNA[queryHits(findOverlaps(binslist5,DNA_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,DNA_gr))))
  gm12878_10kb$DNA[-queryHits(findOverlaps(binslist5,DNA_gr))] <- 0
  
  gm12878_10kb$line <- NA
  gm12878_10kb$line[queryHits(findOverlaps(binslist5,line_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,line_gr))))
  gm12878_10kb$line[-queryHits(findOverlaps(binslist5,line_gr))] <- 0
  
  gm12878_10kb$low_complexity <- NA
  gm12878_10kb$low_complexity[queryHits(findOverlaps(binslist5,low_complexity_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,low_complexity_gr))))
  gm12878_10kb$low_complexity[-queryHits(findOverlaps(binslist5,low_complexity_gr))] <- 0
  
  gm12878_10kb$LTR <- NA
  gm12878_10kb$LTR[queryHits(findOverlaps(binslist5,LTR_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,LTR_gr))))
  gm12878_10kb$LTR[-queryHits(findOverlaps(binslist5,LTR_gr))] <- 0
  
  gm12878_10kb$other <- NA
  gm12878_10kb$other[queryHits(findOverlaps(binslist5,other_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,other_gr))))
  gm12878_10kb$other[-queryHits(findOverlaps(binslist5,other_gr))] <- 0
  
  gm12878_10kb$RC <- NA
  gm12878_10kb$RC[queryHits(findOverlaps(binslist5,RC_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,RC_gr))))
  gm12878_10kb$RC[-queryHits(findOverlaps(binslist5,RC_gr))] <- 0
  
  gm12878_10kb$satellite <- NA
  gm12878_10kb$satellite[queryHits(findOverlaps(binslist5,satellite_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,satellite_gr))))
  gm12878_10kb$satellite[-queryHits(findOverlaps(binslist5,satellite_gr))] <- 0
  
  gm12878_10kb$simple_repeat <- NA
  gm12878_10kb$simple_repeat[queryHits(findOverlaps(binslist5,simple_repeat_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,simple_repeat_gr))))
  gm12878_10kb$simple_repeat[-queryHits(findOverlaps(binslist5,simple_repeat_gr))] <- 0
  
  gm12878_10kb$SINE <- NA
  gm12878_10kb$SINE[queryHits(findOverlaps(binslist5,SINE_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,SINE_gr))))
  gm12878_10kb$SINE[-queryHits(findOverlaps(binslist5,SINE_gr))] <- 0
  
  
  gm12878_10kb$DNA_dist <- mcols(distanceToNearest(binslist5_center, DNA_gr_center))$distance 
  gm12878_10kb$line_dist <- mcols(distanceToNearest(binslist5_center, line_gr_center))$distance 
  gm12878_10kb$low_complexity_dist <- mcols(distanceToNearest(binslist5_center, low_complexity_gr_center))$distance 
  gm12878_10kb$LTR_dist <- mcols(distanceToNearest(binslist5_center, LTR_gr_center))$distance 
  gm12878_10kb$other_dist <- mcols(distanceToNearest(binslist5_center, other_gr_center))$distance 
  gm12878_10kb$RC_dist <- mcols(distanceToNearest(binslist5_center, RC_gr_center))$distance 
  gm12878_10kb$satellite_dist <- mcols(distanceToNearest(binslist5_center, satellite_gr_center))$distance 
  gm12878_10kb$simple_repeat_dist <- mcols(distanceToNearest(binslist5_center, simple_repeat_gr_center))$distance 
  gm12878_10kb$SINE_dist <- mcols(distanceToNearest(binslist5_center, SINE_gr_center))$distance 
  
  gm12878_10kb$DNA_count <- countOverlaps(binslist5, DNA_gr)
  gm12878_10kb$line_count <- countOverlaps(binslist5, line_gr)
  gm12878_10kb$low_complexity_count <- countOverlaps(binslist5, low_complexity_gr)
  gm12878_10kb$LTR_count <- countOverlaps(binslist5, LTR_gr)
  gm12878_10kb$other_count <- countOverlaps(binslist5, other_gr)
  gm12878_10kb$RC_count <- countOverlaps(binslist5, RC_gr)
  gm12878_10kb$satellite_count <- countOverlaps(binslist5, satellite_gr)
  gm12878_10kb$simple_repeat_count <- countOverlaps(binslist5, simple_repeat_gr)
  gm12878_10kb$SINE_count <- countOverlaps(binslist5, SINE_gr)
  

## super_enhancers

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/super_enhancers/")

temp = list.files()

se_GM12878 <- read.table(temp[1],header=FALSE,sep="\t")

se_GM12878_gr <- GRanges(seqnames=se_GM12878$V1,IRanges(start=se_GM12878$V2,end=se_GM12878$V3))
se_GM12878_gr_center <- resize(se_GM12878_gr, width = 1, fix = "center")

  gm12878_10kb$se_GM12878 <- NA
  gm12878_10kb$se_GM12878[queryHits(findOverlaps(binslist5,se_GM12878_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,se_GM12878_gr))))
  gm12878_10kb$se_GM12878[-queryHits(findOverlaps(binslist5,se_GM12878_gr))] <- 0
  
  gm12878_10kb$se_GM12878_dist <- mcols(distanceToNearest(binslist5_center, se_GM12878_gr_center))$distance
  
  gm12878_10kb$se_GM12878_count <- countOverlaps(binslist5, se_GM12878_gr)
  

## VMR

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/VMRs/")

VMR <- read.table("VMR_hg19.bed",header=FALSE,sep="\t")

VMR_gr <- GRanges(seqnames=VMR$V1,IRanges(start=VMR$V2,end=VMR$V3))
VMR_gr_center <- resize(VMR_gr, width = 1, fix = "center")

  gm12878_10kb$VMR <- NA
  gm12878_10kb$VMR[queryHits(findOverlaps(binslist5,VMR_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,VMR_gr))))
  gm12878_10kb$VMR[-queryHits(findOverlaps(binslist5,VMR_gr))] <- 0
  
  gm12878_10kb$VMR_dist <- mcols(distanceToNearest(binslist5_center, VMR_gr_center))$distance
  
  gm12878_10kb$VMR_count <- countOverlaps(binslist5, VMR_gr)

  
## BroadHMM

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/BroadHMM/")

temp <- list.files()

#Gm12878
Gm12878_TxnElongation <- read.table(temp[1],header=FALSE,sep="\t")
Gm12878_WeakTxn <- read.table(temp[2],header=FALSE,sep="\t")
Gm12878_Repressed <- read.table(temp[3],header=FALSE,sep="\t")
Gm12878_Heterochromlo <- read.table(temp[4],header=FALSE,sep="\t")
Gm12878_RepetitiveCNV14 <- read.table(temp[5],header=FALSE,sep="\t")
Gm12878_RepetitiveCNV15 <- read.table(temp[6],header=FALSE,sep="\t")
Gm12878_ActivePromoter <- read.table(temp[7],header=FALSE,sep="\t")
Gm12878_WeakPromoter <- read.table(temp[8],header=FALSE,sep="\t")
Gm12878_PoisedPromoter <- read.table(temp[9],header=FALSE,sep="\t")
Gm12878_StrongEnhancer4 <- read.table(temp[10],header=FALSE,sep="\t")
Gm12878_StrongEnhancer5 <- read.table(temp[11],header=FALSE,sep="\t") 
Gm12878_WeakEnhancer6 <- read.table(temp[12],header=FALSE,sep="\t")
Gm12878_WeakEnhancer7 <- read.table(temp[13],header=FALSE,sep="\t")
Gm12878_Insulator <- read.table(temp[14],header=FALSE,sep="\t")
Gm12878_TxnTransition <- read.table(temp[15],header=FALSE,sep="\t")

Gm12878_TxnElongation_gr <- GRanges(seqnames=Gm12878_TxnElongation$V1,IRanges(start=Gm12878_TxnElongation$V2,end=Gm12878_TxnElongation$V3))
Gm12878_TxnElongation_gr_center <- resize(Gm12878_TxnElongation_gr, width = 1, fix = "center")
Gm12878_WeakTxn_gr <- GRanges(seqnames=Gm12878_WeakTxn$V1,IRanges(start=Gm12878_WeakTxn$V2,end=Gm12878_WeakTxn$V3))
Gm12878_WeakTxn_gr_center <- resize(Gm12878_WeakTxn_gr, width = 1, fix = "center")
Gm12878_Repressed_gr <- GRanges(seqnames=Gm12878_Repressed$V1,IRanges(start=Gm12878_Repressed$V2,end=Gm12878_Repressed$V3))
Gm12878_Repressed_gr_center <- resize(Gm12878_Repressed_gr, width = 1, fix = "center")
Gm12878_Heterochromlo_gr <- GRanges(seqnames=Gm12878_Heterochromlo$V1,IRanges(start=Gm12878_Heterochromlo$V2,end=Gm12878_Heterochromlo$V3)) 
Gm12878_Heterochromlo_gr_center <- resize(Gm12878_Heterochromlo_gr, width = 1, fix = "center")
Gm12878_RepetitiveCNV14_gr <- GRanges(seqnames=Gm12878_RepetitiveCNV14$V1,IRanges(start=Gm12878_RepetitiveCNV14$V2,end=Gm12878_RepetitiveCNV14$V3)) 
Gm12878_RepetitiveCNV14_gr_center <- resize(Gm12878_RepetitiveCNV14_gr, width = 1, fix = "center")
Gm12878_RepetitiveCNV15_gr <- GRanges(seqnames=Gm12878_RepetitiveCNV15$V1,IRanges(start=Gm12878_RepetitiveCNV15$V2,end=Gm12878_RepetitiveCNV15$V3))
Gm12878_RepetitiveCNV15_gr_center <- resize(Gm12878_RepetitiveCNV15_gr, width = 1, fix = "center")
Gm12878_ActivePromoter_gr <- GRanges(seqnames=Gm12878_ActivePromoter$V1,IRanges(start=Gm12878_ActivePromoter$V2,end=Gm12878_ActivePromoter$V3))
Gm12878_ActivePromoter_gr_center <- resize(Gm12878_ActivePromoter_gr, width = 1, fix = "center")
Gm12878_WeakPromoter_gr <- GRanges(seqnames=Gm12878_WeakPromoter$V1,IRanges(start=Gm12878_WeakPromoter$V2,end=Gm12878_WeakPromoter$V3))
Gm12878_WeakPromoter_gr_center <- resize(Gm12878_WeakPromoter_gr, width = 1, fix = "center")
Gm12878_PoisedPromoter_gr <- GRanges(seqnames=Gm12878_PoisedPromoter$V1,IRanges(start=Gm12878_PoisedPromoter$V2,end=Gm12878_PoisedPromoter$V3)) 
Gm12878_PoisedPromoter_gr_center <- resize(Gm12878_PoisedPromoter_gr, width = 1, fix = "center")
Gm12878_StrongEnhancer4_gr <- GRanges(seqnames=Gm12878_StrongEnhancer4$V1,IRanges(start=Gm12878_StrongEnhancer4$V2,end=Gm12878_StrongEnhancer4$V3))
Gm12878_StrongEnhancer4_gr_center <- resize(Gm12878_StrongEnhancer4_gr, width = 1, fix = "center") 
Gm12878_StrongEnhancer5_gr <- GRanges(seqnames=Gm12878_StrongEnhancer5$V1,IRanges(start=Gm12878_StrongEnhancer5$V2,end=Gm12878_StrongEnhancer5$V3))
Gm12878_StrongEnhancer5_gr_center <- resize(Gm12878_StrongEnhancer5_gr, width = 1, fix = "center")
Gm12878_WeakEnhancer6_gr <- GRanges(seqnames=Gm12878_WeakEnhancer6$V1,IRanges(start=Gm12878_WeakEnhancer6$V2,end=Gm12878_WeakEnhancer6$V3)) 
Gm12878_WeakEnhancer6_gr_center <- resize(Gm12878_WeakEnhancer6_gr, width = 1, fix = "center")
Gm12878_WeakEnhancer7_gr <- GRanges(seqnames=Gm12878_WeakEnhancer7$V1,IRanges(start=Gm12878_WeakEnhancer7$V2,end=Gm12878_WeakEnhancer7$V3))
Gm12878_WeakEnhancer7_gr_center <- resize(Gm12878_WeakEnhancer7_gr, width = 1, fix = "center") 
Gm12878_Insulator_gr <- GRanges(seqnames=Gm12878_Insulator$V1,IRanges(start=Gm12878_Insulator$V2,end=Gm12878_Insulator$V3))
Gm12878_Insulator_gr_center <- resize(Gm12878_Insulator_gr, width = 1, fix = "center")
Gm12878_TxnTransition_gr <- GRanges(seqnames=Gm12878_TxnTransition$V1,IRanges(start=Gm12878_TxnTransition$V2,end=Gm12878_TxnTransition$V3)) 
Gm12878_TxnTransition_gr_center <- resize(Gm12878_TxnTransition_gr, width = 1, fix = "center")

  gm12878_10kb$Gm12878_TxnElongation <- NA
  gm12878_10kb$Gm12878_TxnElongation[queryHits(findOverlaps(binslist5,Gm12878_TxnElongation_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_TxnElongation_gr))))
  gm12878_10kb$Gm12878_TxnElongation[-queryHits(findOverlaps(binslist5,Gm12878_TxnElongation_gr))] <- 0
  
  gm12878_10kb$Gm12878_WeakTxn <- NA
  gm12878_10kb$Gm12878_WeakTxn[queryHits(findOverlaps(binslist5,Gm12878_WeakTxn_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_WeakTxn_gr))))
  gm12878_10kb$Gm12878_WeakTxn[-queryHits(findOverlaps(binslist5,Gm12878_WeakTxn_gr))] <- 0
  
  gm12878_10kb$Gm12878_Repressed <- NA
  gm12878_10kb$Gm12878_Repressed[queryHits(findOverlaps(binslist5,Gm12878_Repressed_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_Repressed_gr))))
  gm12878_10kb$Gm12878_Repressed[-queryHits(findOverlaps(binslist5,Gm12878_Repressed_gr))] <- 0
  
  gm12878_10kb$Gm12878_Heterochromlo <- NA
  gm12878_10kb$Gm12878_Heterochromlo[queryHits(findOverlaps(binslist5,Gm12878_Heterochromlo_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_Heterochromlo_gr))))
  gm12878_10kb$Gm12878_Heterochromlo[-queryHits(findOverlaps(binslist5,Gm12878_Heterochromlo_gr))] <- 0
  
  gm12878_10kb$Gm12878_RepetitiveCNV14 <- NA
  gm12878_10kb$Gm12878_RepetitiveCNV14[queryHits(findOverlaps(binslist5,Gm12878_RepetitiveCNV14_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_RepetitiveCNV14_gr))))
  gm12878_10kb$Gm12878_RepetitiveCNV14[-queryHits(findOverlaps(binslist5,Gm12878_RepetitiveCNV14_gr))] <- 0
  
  gm12878_10kb$Gm12878_RepetitiveCNV15 <- NA
  gm12878_10kb$Gm12878_RepetitiveCNV15[queryHits(findOverlaps(binslist5,Gm12878_RepetitiveCNV15_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_RepetitiveCNV15_gr))))
  gm12878_10kb$Gm12878_RepetitiveCNV15[-queryHits(findOverlaps(binslist5,Gm12878_RepetitiveCNV15_gr))] <- 0
  
  gm12878_10kb$Gm12878_ActivePromoter <- NA
  gm12878_10kb$Gm12878_ActivePromoter[queryHits(findOverlaps(binslist5,Gm12878_ActivePromoter_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_ActivePromoter_gr))))
  gm12878_10kb$Gm12878_ActivePromoter[-queryHits(findOverlaps(binslist5,Gm12878_ActivePromoter_gr))] <- 0
  
  gm12878_10kb$Gm12878_WeakPromoter <- NA
  gm12878_10kb$Gm12878_WeakPromoter[queryHits(findOverlaps(binslist5,Gm12878_WeakPromoter_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_WeakPromoter_gr))))
  gm12878_10kb$Gm12878_WeakPromoter[-queryHits(findOverlaps(binslist5,Gm12878_WeakPromoter_gr))] <- 0
  
  gm12878_10kb$Gm12878_PoisedPromoter <- NA
  gm12878_10kb$Gm12878_PoisedPromoter[queryHits(findOverlaps(binslist5,Gm12878_PoisedPromoter_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_PoisedPromoter_gr))))
  gm12878_10kb$Gm12878_PoisedPromoter[-queryHits(findOverlaps(binslist5,Gm12878_PoisedPromoter_gr))] <- 0
  
  gm12878_10kb$Gm12878_StrongEnhancer4 <- NA
  gm12878_10kb$Gm12878_StrongEnhancer4[queryHits(findOverlaps(binslist5,Gm12878_StrongEnhancer4_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_StrongEnhancer4_gr))))
  gm12878_10kb$Gm12878_StrongEnhancer4[-queryHits(findOverlaps(binslist5,Gm12878_StrongEnhancer4_gr))] <- 0
  
  gm12878_10kb$Gm12878_StrongEnhancer5 <- NA
  gm12878_10kb$Gm12878_StrongEnhancer5[queryHits(findOverlaps(binslist5,Gm12878_StrongEnhancer5_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_StrongEnhancer5_gr))))
  gm12878_10kb$Gm12878_StrongEnhancer5[-queryHits(findOverlaps(binslist5,Gm12878_StrongEnhancer5_gr))] <- 0
  
  gm12878_10kb$Gm12878_WeakEnhancer6 <- NA
  gm12878_10kb$Gm12878_WeakEnhancer6[queryHits(findOverlaps(binslist5,Gm12878_WeakEnhancer6_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_WeakEnhancer6_gr))))
  gm12878_10kb$Gm12878_WeakEnhancer6[-queryHits(findOverlaps(binslist5,Gm12878_WeakEnhancer6_gr))] <- 0
  
  gm12878_10kb$Gm12878_WeakEnhancer7 <- NA
  gm12878_10kb$Gm12878_WeakEnhancer7[queryHits(findOverlaps(binslist5,Gm12878_WeakEnhancer7_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_WeakEnhancer7_gr))))
  gm12878_10kb$Gm12878_WeakEnhancer7[-queryHits(findOverlaps(binslist5,Gm12878_WeakEnhancer7_gr))] <- 0
  
  gm12878_10kb$Gm12878_Insulator <- NA
  gm12878_10kb$Gm12878_Insulator[queryHits(findOverlaps(binslist5,Gm12878_Insulator_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_Insulator_gr))))
  gm12878_10kb$Gm12878_Insulator[-queryHits(findOverlaps(binslist5,Gm12878_Insulator_gr))] <- 0
  
  gm12878_10kb$Gm12878_TxnTransition <- NA
  gm12878_10kb$Gm12878_TxnTransition[queryHits(findOverlaps(binslist5,Gm12878_TxnTransition_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_TxnTransition_gr))))
  gm12878_10kb$Gm12878_TxnTransition[-queryHits(findOverlaps(binslist5,Gm12878_TxnTransition_gr))] <- 0
  

  gm12878_10kb$Gm12878_TxnElongation_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_TxnElongation_gr_center))$distance
  gm12878_10kb$Gm12878_WeakTxn_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WeakTxn_gr_center))$distance
  gm12878_10kb$Gm12878_Repressed_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_Repressed_gr_center))$distance
  gm12878_10kb$Gm12878_Heterochromlo_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_Heterochromlo_gr_center))$distance
  gm12878_10kb$Gm12878_RepetitiveCNV14_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_RepetitiveCNV14_gr_center))$distance
  gm12878_10kb$Gm12878_RepetitiveCNV15_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_RepetitiveCNV15_gr_center))$distance
  gm12878_10kb$Gm12878_ActivePromoter_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_ActivePromoter_gr_center))$distance
  gm12878_10kb$Gm12878_WeakPromoter_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WeakPromoter_gr_center))$distance
  gm12878_10kb$Gm12878_PoisedPromoter_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_PoisedPromoter_gr_center))$distance
  gm12878_10kb$Gm12878_StrongEnhancer4_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_StrongEnhancer4_gr_center))$distance
  gm12878_10kb$Gm12878_StrongEnhancer5_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_StrongEnhancer5_gr_center))$distance
  gm12878_10kb$Gm12878_WeakEnhancer6_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WeakEnhancer6_gr_center))$distance
  gm12878_10kb$Gm12878_WeakEnhancer7_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WeakEnhancer7_gr_center))$distance
  gm12878_10kb$Gm12878_Insulator_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_Insulator_gr_center))$distance
  gm12878_10kb$Gm12878_TxnTransition_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_TxnTransition_gr_center))$distance
  
  gm12878_10kb$Gm12878_TxnElongation_count <- countOverlaps(binslist5, Gm12878_TxnElongation_gr)
  gm12878_10kb$Gm12878_WeakTxn_count <- countOverlaps(binslist5, Gm12878_WeakTxn_gr)
  gm12878_10kb$Gm12878_Repressed_count <- countOverlaps(binslist5, Gm12878_Repressed_gr)
  gm12878_10kb$Gm12878_Heterochromlo_count <- countOverlaps(binslist5, Gm12878_Heterochromlo_gr)
  gm12878_10kb$Gm12878_RepetitiveCNV14_count <- countOverlaps(binslist5, Gm12878_RepetitiveCNV14_gr)
  gm12878_10kb$Gm12878_RepetitiveCNV15_count <- countOverlaps(binslist5, Gm12878_RepetitiveCNV15_gr)
  gm12878_10kb$Gm12878_ActivePromoter_count <- countOverlaps(binslist5, Gm12878_ActivePromoter_gr)
  gm12878_10kb$Gm12878_WeakPromoter_count <- countOverlaps(binslist5, Gm12878_WeakPromoter_gr)
  gm12878_10kb$Gm12878_PoisedPromoter_count <- countOverlaps(binslist5, Gm12878_PoisedPromoter_gr)
  gm12878_10kb$Gm12878_StrongEnhancer4_count <- countOverlaps(binslist5, Gm12878_StrongEnhancer4_gr)
  gm12878_10kb$Gm12878_StrongEnhancer5_count <- countOverlaps(binslist5, Gm12878_StrongEnhancer5_gr)
  gm12878_10kb$Gm12878_WeakEnhancer6_count <- countOverlaps(binslist5, Gm12878_WeakEnhancer6_gr)
  gm12878_10kb$Gm12878_WeakEnhancer7_count <- countOverlaps(binslist5, Gm12878_WeakEnhancer7_gr)
  gm12878_10kb$Gm12878_Insulator_count <- countOverlaps(binslist5, Gm12878_Insulator_gr)
  gm12878_10kb$Gm12878_TxnTransition_count <- countOverlaps(binslist5, Gm12878_TxnTransition_gr)

  
## Combined

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/Combined/")

temp <- list.files()

#Gm12878
Gm12878_CTCF <- read.table(temp[1],header=FALSE,sep="\t") 
Gm12878_E <- read.table(temp[2],header=FALSE,sep="\t")
Gm12878_PF <- read.table(temp[3],header=FALSE,sep="\t")
Gm12878_R <- read.table(temp[4],header=FALSE,sep="\t")
Gm12878_T <- read.table(temp[5],header=FALSE,sep="\t")
Gm12878_TSS <- read.table(temp[6],header=FALSE,sep="\t")
Gm12878_WE <- read.table(temp[7],header=FALSE,sep="\t")

Gm12878_CTCF_gr <- GRanges(seqnames=Gm12878_CTCF$V1,IRanges(start=Gm12878_CTCF$V2,end=Gm12878_CTCF$V3))
Gm12878_CTCF_gr_center <- resize(Gm12878_CTCF_gr, width = 1, fix = "center")
Gm12878_E_gr <- GRanges(seqnames=Gm12878_E$V1,IRanges(start=Gm12878_E$V2,end=Gm12878_E$V3))
Gm12878_E_gr_center <- resize(Gm12878_E_gr, width = 1, fix = "center")
Gm12878_PF_gr <- GRanges(seqnames=Gm12878_PF$V1,IRanges(start=Gm12878_PF$V2,end=Gm12878_PF$V3))
Gm12878_PF_gr_center <- resize(Gm12878_PF_gr, width = 1, fix = "center")
Gm12878_R_gr <- GRanges(seqnames=Gm12878_R$V1,IRanges(start=Gm12878_R$V2,end=Gm12878_R$V3))
Gm12878_R_gr_center <- resize(Gm12878_R_gr, width = 1, fix = "center")
Gm12878_T_gr <- GRanges(seqnames=Gm12878_T$V1,IRanges(start=Gm12878_T$V2,end=Gm12878_T$V3))
Gm12878_T_gr_center <- resize(Gm12878_T_gr, width = 1, fix = "center")
Gm12878_TSS_gr <- GRanges(seqnames=Gm12878_TSS$V1,IRanges(start=Gm12878_TSS$V2,end=Gm12878_TSS$V3))
Gm12878_TSS_gr_center <- resize(Gm12878_TSS_gr, width = 1, fix = "center")
Gm12878_WE_gr <- GRanges(seqnames=Gm12878_WE$V1,IRanges(start=Gm12878_WE$V2,end=Gm12878_WE$V3))
Gm12878_WE_gr_center <- resize(Gm12878_WE_gr, width = 1, fix = "center")

  gm12878_10kb$Gm12878_CTCF <- NA
  gm12878_10kb$Gm12878_CTCF[queryHits(findOverlaps(binslist5,Gm12878_CTCF_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_CTCF_gr))))
  gm12878_10kb$Gm12878_CTCF[-queryHits(findOverlaps(binslist5,Gm12878_CTCF_gr))] <- 0
  
  gm12878_10kb$Gm12878_E <- NA
  gm12878_10kb$Gm12878_E[queryHits(findOverlaps(binslist5,Gm12878_E_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_E_gr))))
  gm12878_10kb$Gm12878_E[-queryHits(findOverlaps(binslist5,Gm12878_E_gr))] <- 0
  
  gm12878_10kb$Gm12878_PF <- NA
  gm12878_10kb$Gm12878_PF[queryHits(findOverlaps(binslist5,Gm12878_PF_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_PF_gr))))
  gm12878_10kb$Gm12878_PF[-queryHits(findOverlaps(binslist5,Gm12878_PF_gr))] <- 0
  
  gm12878_10kb$Gm12878_R <- NA
  gm12878_10kb$Gm12878_R[queryHits(findOverlaps(binslist5,Gm12878_R_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_R_gr))))
  gm12878_10kb$Gm12878_R[-queryHits(findOverlaps(binslist5,Gm12878_R_gr))] <- 0
  
  gm12878_10kb$Gm12878_T <- NA
  gm12878_10kb$Gm12878_T[queryHits(findOverlaps(binslist5,Gm12878_T_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_T_gr))))
  gm12878_10kb$Gm12878_T[-queryHits(findOverlaps(binslist5,Gm12878_T_gr))] <- 0
  
  gm12878_10kb$Gm12878_TSS <- NA
  gm12878_10kb$Gm12878_TSS[queryHits(findOverlaps(binslist5,Gm12878_TSS_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_TSS_gr))))
  gm12878_10kb$Gm12878_TSS[-queryHits(findOverlaps(binslist5,Gm12878_TSS_gr))] <- 0
  
  gm12878_10kb$Gm12878_WE <- NA
  gm12878_10kb$Gm12878_WE[queryHits(findOverlaps(binslist5,Gm12878_WE_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_WE_gr))))
  gm12878_10kb$Gm12878_WE[-queryHits(findOverlaps(binslist5,Gm12878_WE_gr))] <- 0
  
  
  gm12878_10kb$Gm12878_CTCF_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_CTCF_gr_center))$distance
  gm12878_10kb$Gm12878_E_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_E_gr_center))$distance
  gm12878_10kb$Gm12878_PF_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_PF_gr_center))$distance
  gm12878_10kb$Gm12878_R_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_R_gr_center))$distance
  gm12878_10kb$Gm12878_T_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_T_gr_center))$distance
  gm12878_10kb$Gm12878_TSS_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_TSS_gr_center))$distance
  gm12878_10kb$Gm12878_WE_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_WE_gr_center))$distance
  
  gm12878_10kb$Gm12878_CTCF_count <- countOverlaps(binslist5, Gm12878_CTCF_gr)
  gm12878_10kb$Gm12878_E_count <- countOverlaps(binslist5, Gm12878_E_gr)
  gm12878_10kb$Gm12878_PF_count <- countOverlaps(binslist5, Gm12878_PF_gr)
  gm12878_10kb$Gm12878_R_count <- countOverlaps(binslist5, Gm12878_R_gr)
  gm12878_10kb$Gm12878_T_count <- countOverlaps(binslist5, Gm12878_T_gr)
  gm12878_10kb$Gm12878_TSS_count <- countOverlaps(binslist5, Gm12878_TSS_gr)
  gm12878_10kb$Gm12878_WE_count <- countOverlaps(binslist5, Gm12878_WE_gr)
  

## DNase I


setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DNaseI/")

temp <- list.files()

Gm12878_DNaseI <- read.table(temp[1],header=FALSE,sep="\t") 

Gm12878_DNaseI_gr <- GRanges(seqnames=Gm12878_DNaseI$V1,IRanges(start=Gm12878_DNaseI$V2,end=Gm12878_DNaseI$V3))
Gm12878_DNaseI_gr_center <- resize(Gm12878_DNaseI_gr, width = 1, fix = "center")

  gm12878_10kb$Gm12878_DNaseI <- NA
  gm12878_10kb$Gm12878_DNaseI[queryHits(findOverlaps(binslist5,Gm12878_DNaseI_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_DNaseI_gr))))
  gm12878_10kb$Gm12878_DNaseI[-queryHits(findOverlaps(binslist5,Gm12878_DNaseI_gr))] <- 0

  gm12878_10kb$Gm12878_DNaseI_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_DNaseI_gr_center))$distance
  
  gm12878_10kb$Gm12878_DNaseI_count <- countOverlaps(binslist5, Gm12878_DNaseI_gr)


## Histone Modifications

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/HistoneModifications/")

temp <- list.files()

Gm12878_H2az <- read.table(temp[1],header=FALSE,sep="\t") 
Gm12878_H3k27ac <- read.table(temp[2],header=FALSE,sep="\t")
Gm12878_H3k27me3 <- read.table(temp[3],header=FALSE,sep="\t")
Gm12878_H3k36me3 <- read.table(temp[4],header=FALSE,sep="\t")
Gm12878_H3k4me1 <- read.table(temp[5],header=FALSE,sep="\t")
Gm12878_H3k4me2 <- read.table(temp[6],header=FALSE,sep="\t")
Gm12878_H3k4me3 <- read.table(temp[7],header=FALSE,sep="\t")
Gm12878_H3k79me2 <- read.table(temp[8],header=FALSE,sep="\t")
Gm12878_H3k9ac <- read.table(temp[9],header=FALSE,sep="\t")
Gm12878_H3k9me3 <- read.table(temp[10],header=FALSE,sep="\t")
Gm12878_H4k20me1 <- read.table(temp[11],header=FALSE,sep="\t")

Gm12878_H2az_gr <- GRanges(seqnames=Gm12878_H2az$V1,IRanges(start=Gm12878_H2az$V2,end=Gm12878_H2az$V3))
Gm12878_H2az_gr_center <- resize(Gm12878_H2az_gr, width = 1, fix = "center")
Gm12878_H3k27ac_gr <- GRanges(seqnames=Gm12878_H3k27ac$V1,IRanges(start=Gm12878_H3k27ac$V2,end=Gm12878_H3k27ac$V3))
Gm12878_H3k27ac_gr_center <- resize(Gm12878_H3k27ac_gr, width = 1, fix = "center")
Gm12878_H3k27me3_gr <- GRanges(seqnames=Gm12878_H3k27me3$V1,IRanges(start=Gm12878_H3k27me3$V2,end=Gm12878_H3k27me3$V3))
Gm12878_H3k27me3_gr_center <- resize(Gm12878_H3k27me3_gr, width = 1, fix = "center")
Gm12878_H3k36me3_gr <- GRanges(seqnames=Gm12878_H3k36me3$V1,IRanges(start=Gm12878_H3k36me3$V2,end=Gm12878_H3k36me3$V3))
Gm12878_H3k36me3_gr_center <- resize(Gm12878_H3k36me3_gr, width = 1, fix = "center")
Gm12878_H3k4me1_gr <- GRanges(seqnames=Gm12878_H3k4me1$V1,IRanges(start=Gm12878_H3k4me1$V2,end=Gm12878_H3k4me1$V3))
Gm12878_H3k4me1_gr_center <- resize(Gm12878_H3k4me1_gr, width = 1, fix = "center")
Gm12878_H3k4me2_gr <- GRanges(seqnames=Gm12878_H3k4me2$V1,IRanges(start=Gm12878_H3k4me2$V2,end=Gm12878_H3k4me2$V3))
Gm12878_H3k4me2_gr_center <- resize(Gm12878_H3k4me2_gr, width = 1, fix = "center")
Gm12878_H3k4me3_gr <- GRanges(seqnames=Gm12878_H3k4me3$V1,IRanges(start=Gm12878_H3k4me3$V2,end=Gm12878_H3k4me3$V3))
Gm12878_H3k4me3_gr_center <- resize(Gm12878_H3k4me3_gr, width = 1, fix = "center")
Gm12878_H3k79me2_gr <- GRanges(seqnames=Gm12878_H3k79me2$V1,IRanges(start=Gm12878_H3k79me2$V2,end=Gm12878_H3k79me2$V3))
Gm12878_H3k79me2_gr_center <- resize(Gm12878_H3k79me2_gr, width = 1, fix = "center")
Gm12878_H3k9ac_gr <- GRanges(seqnames=Gm12878_H3k9ac$V1,IRanges(start=Gm12878_H3k9ac$V2,end=Gm12878_H3k9ac$V3))
Gm12878_H3k9ac_gr_center <- resize(Gm12878_H3k9ac_gr, width = 1, fix = "center")
Gm12878_H3k9me3_gr <- GRanges(seqnames=Gm12878_H3k9me3$V1,IRanges(start=Gm12878_H3k9me3$V2,end=Gm12878_H3k9me3$V3))
Gm12878_H3k9me3_gr_center <- resize(Gm12878_H3k9me3_gr, width = 1, fix = "center")
Gm12878_H4k20me1_gr <- GRanges(seqnames=Gm12878_H4k20me1$V1,IRanges(start=Gm12878_H4k20me1$V2,end=Gm12878_H4k20me1$V3))
Gm12878_H4k20me1_gr_center <- resize(Gm12878_H4k20me1_gr, width = 1, fix = "center")

  gm12878_10kb$Gm12878_H2az <- NA
  gm12878_10kb$Gm12878_H2az[queryHits(findOverlaps(binslist5,Gm12878_H2az_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H2az_gr))))
  gm12878_10kb$Gm12878_H2az[-queryHits(findOverlaps(binslist5,Gm12878_H2az_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k27ac <- NA
  gm12878_10kb$Gm12878_H3k27ac[queryHits(findOverlaps(binslist5,Gm12878_H3k27ac_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k27ac_gr))))
  gm12878_10kb$Gm12878_H3k27ac[-queryHits(findOverlaps(binslist5,Gm12878_H3k27ac_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k27me3 <- NA
  gm12878_10kb$Gm12878_H3k27me3[queryHits(findOverlaps(binslist5,Gm12878_H3k27me3_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k27me3_gr))))
  gm12878_10kb$Gm12878_H3k27me3[-queryHits(findOverlaps(binslist5,Gm12878_H3k27me3_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k36me3 <- NA
  gm12878_10kb$Gm12878_H3k36me3[queryHits(findOverlaps(binslist5,Gm12878_H3k36me3_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k36me3_gr))))
  gm12878_10kb$Gm12878_H3k36me3[-queryHits(findOverlaps(binslist5,Gm12878_H3k36me3_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k4me1 <- NA
  gm12878_10kb$Gm12878_H3k4me1[queryHits(findOverlaps(binslist5,Gm12878_H3k4me1_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k4me1_gr))))
  gm12878_10kb$Gm12878_H3k4me1[-queryHits(findOverlaps(binslist5,Gm12878_H3k4me1_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k4me2 <- NA
  gm12878_10kb$Gm12878_H3k4me2[queryHits(findOverlaps(binslist5,Gm12878_H3k4me2_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k4me2_gr))))
  gm12878_10kb$Gm12878_H3k4me2[-queryHits(findOverlaps(binslist5,Gm12878_H3k4me2_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k4me3 <- NA
  gm12878_10kb$Gm12878_H3k4me3[queryHits(findOverlaps(binslist5,Gm12878_H3k4me3_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k4me3_gr))))
  gm12878_10kb$Gm12878_H3k4me3[-queryHits(findOverlaps(binslist5,Gm12878_H3k4me3_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k79me2 <- NA
  gm12878_10kb$Gm12878_H3k79me2[queryHits(findOverlaps(binslist5,Gm12878_H3k79me2_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k79me2_gr))))
  gm12878_10kb$Gm12878_H3k79me2[-queryHits(findOverlaps(binslist5,Gm12878_H3k79me2_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k9ac <- NA
  gm12878_10kb$Gm12878_H3k9ac[queryHits(findOverlaps(binslist5,Gm12878_H3k9ac_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k9ac_gr))))
  gm12878_10kb$Gm12878_H3k9ac[-queryHits(findOverlaps(binslist5,Gm12878_H3k9ac_gr))] <- 0
  
  gm12878_10kb$Gm12878_H3k9me3 <- NA
  gm12878_10kb$Gm12878_H3k9me3[queryHits(findOverlaps(binslist5,Gm12878_H3k9me3_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H3k9me3_gr))))
  gm12878_10kb$Gm12878_H3k9me3[-queryHits(findOverlaps(binslist5,Gm12878_H3k9me3_gr))] <- 0
  
  gm12878_10kb$Gm12878_H4k20me1 <- NA
  gm12878_10kb$Gm12878_H4k20me1[queryHits(findOverlaps(binslist5,Gm12878_H4k20me1_gr))] <- width(ranges(pintersect(findOverlapPairs(binslist5,Gm12878_H4k20me1_gr))))
  gm12878_10kb$Gm12878_H4k20me1[-queryHits(findOverlaps(binslist5,Gm12878_H4k20me1_gr))] <- 0


  gm12878_10kb$Gm12878_H2az_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H2az_gr_center))$distance
  gm12878_10kb$Gm12878_H3k27ac_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k27ac_gr_center))$distance
  gm12878_10kb$Gm12878_H3k27me3_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k27me3_gr_center))$distance
  gm12878_10kb$Gm12878_H3k36me3_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k36me3_gr_center))$distance
  gm12878_10kb$Gm12878_H3k4me1_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k4me1_gr_center))$distance
  gm12878_10kb$Gm12878_H3k4me2_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k4me2_gr_center))$distance
  gm12878_10kb$Gm12878_H3k4me3_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k4me3_gr_center))$distance
  gm12878_10kb$Gm12878_H3k79me2_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k79me2_gr_center))$distance
  gm12878_10kb$Gm12878_H3k9ac_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k9ac_gr_center))$distance
  gm12878_10kb$Gm12878_H3k9me3_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H3k9me3_gr_center))$distance
  gm12878_10kb$Gm12878_H4k20me1_dist <- mcols(distanceToNearest(binslist5_center, Gm12878_H4k20me1_gr_center))$distance
  
  gm12878_10kb$Gm12878_H2az_count <- countOverlaps(binslist5, Gm12878_H2az_gr)
  gm12878_10kb$Gm12878_H3k27ac_count <- countOverlaps(binslist5, Gm12878_H3k27ac_gr)
  gm12878_10kb$Gm12878_H3k27me3_count <- countOverlaps(binslist5, Gm12878_H3k27me3_gr)
  gm12878_10kb$Gm12878_H3k36me3_count <- countOverlaps(binslist5, Gm12878_H3k36me3_gr)
  gm12878_10kb$Gm12878_H3k4me1_count <- countOverlaps(binslist5, Gm12878_H3k4me1_gr)
  gm12878_10kb$Gm12878_H3k4me2_count <- countOverlaps(binslist5, Gm12878_H3k4me2_gr)
  gm12878_10kb$Gm12878_H3k4me3_count <- countOverlaps(binslist5, Gm12878_H3k4me3_gr)
  gm12878_10kb$Gm12878_H3k79me2_count <- countOverlaps(binslist5, Gm12878_H3k79me2_gr)
  gm12878_10kb$Gm12878_H3k9ac_count <- countOverlaps(binslist5, Gm12878_H3k9ac_gr)
  gm12878_10kb$Gm12878_H3k9me3_count <- countOverlaps(binslist5, Gm12878_H3k9me3_gr)
  gm12878_10kb$Gm12878_H4k20me1_count <- countOverlaps(binslist5, Gm12878_H4k20me1_gr)


## Adding Chromosome information to the data

#gm12878_10kb$CHR <- seqnames(binslist5)

#gm12878_10kb$CHR <- as.character(gm12878_10kb$CHR)


## Saving the data

dim(gm12878_10kb)

head(gm12878_10kb)

saveRDS(gm12878_10kb, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/gm12878_10kb.rds")

## Distance Summaries for each feature from region to TAD boundary

grobjects <- ls()
grobjects <- grobjects[grep("center", grobjects)]

grlist <- rep( list(GRangesList()), length(grobjects) )

for(i in 1:length(grobjects)){
  
  x <- get(grobjects[i])
  
  grlist[[i]] <- x
  
  #grlist[[i]] <- grlist
  
}

grlist <- setNames(grlist, grobjects)

grlist <- grlist[-which(names(grlist)=="binslist5_center")]

length(grlist)

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


saveRDS(grlist, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/grlist.rds")

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

##Taking log2 transform of continous data

#cols <- c(grep("dist",colnames(gm12878_10kb)))
gm12878_10kb[,-1] <- apply(gm12878_10kb[,-1], 2, function(x){log(x + 1, base=2)})


##Changing binary variables to factors

#cols <- c(intersect(grep("score",colnames(gm12878_10kb), invert = TRUE),
#          grep("dist",colnames(gm12878_10kb), invert = TRUE)))
gm12878_10kb$y <- factor(gm12878_10kb$y)


##Changing levels of response (y) to yes no

levels(gm12878_10kb$y) <- c("No", "Yes")

##Removing zero variance predictors

nzv <- nearZeroVar(gm12878_10kb[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])

nzvar

gm12878_10kb_f <- gm12878_10kb[, -which(colnames(gm12878_10kb) %in% nzvar)]

dim(gm12878_10kb) 
dim(gm12878_10kb_f) 

saveRDS(gm12878_10kb_f, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/gm12878_10kb_f.rds")


setwd("/home/stilianoudakisc/TAD_data_analysis/model_filtering/")

cols <- names(gm12878_10kb_f)[grep("dist", names(gm12878_10kb_f))]
featdistdf <- gm12878_10kb_f[,c(1, which(names(gm12878_10kb_f) %in% cols))]

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





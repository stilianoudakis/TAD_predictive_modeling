setwd("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/")

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

domains <- read.table("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/arrowhead_data.txt", header=T)
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

bounds <- GRanges(seqnames=coords$Chromosome, ranges=IRanges(start=coords$coordinate, width=1))

#saving data
saveRDS(bounds, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/bounds.rds")


## Reading in binned genome in the form of contact matrix at 10kb resolution

binslist10 <- read.table("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/binslist10.txt")

dim(binslist10) #268498      2

#ordering the bins according to left endpoint
binslist10 <- binslist10[order(as.numeric(substr(binslist10$V1,4,5)), binslist10$V2, decreasing=FALSE),]

#extracting and renaming first 2 columns
binslist10 <- binslist10[,1:2]
colnames(binslist10) <- c("Chromosome", "Coordinate")


#creating a granges object from binned genome
binslist10 <- GRanges(seqnames = binslist10$Chromosome, ranges = IRanges(start = binslist10$Coordinate,
                                                                         width = 10000))

#saving the binned genome data
saveRDS(binslist10, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/binslist10.rds")


## Creating Response Vector Y (1 if tad boundary is in bin; 0 if not)

#Finding where the TAD boundaries are inside the genomic bins
findOverlaps(binslist10, bounds)

y <- countOverlaps(binslist10, bounds)
length(y) #268498
table(y)
#      0      1      2 
# 252624  14983    891 

#investigating where there are multiple overlaps
head(which(y==2))  #169 736 766 810 857 920

#looking at row 169 of binned genome
findOverlaps(binslist10[169], bounds)
  #     queryHits subjectHits
  #     <integer>   <integer>
  # [1]         1           9
  # [2]         1          10
bounds[9:10]
  #     seqnames             ranges strand
  #        <Rle>          <IRanges>  <Rle>
  # [1]     chr1 [1860000, 1860000]      *
  # [2]     chr1 [1865000, 1865000]      *
#this is a particularly small TAD only 5kb in length
#therefore its starting and ending boundaries are both inside the same genomic bin
#there are 891 of these out of the 15874 TAD boundaries

#recode multiple overlaps as 1
y <- ifelse(y>=1,1,0)
prop.table(table(y)) #0.94087852 0.05912148 
mcols(binslist10)$y <- y


## Creating the data frame with resonse vector for modeling

gm12878_10kb <- data.frame(y = y)


saveRDS(gm12878_10kb, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/gm12878_10kb.rds")

#creating a granges object from the center of every range in the binned genome granges object
#this will be used to calculate distances from genomic feature to center of genomic bin
binslist10_center <- resize(binslist10, width = 1, fix = "center")


# Creating Feature Space

## 3D subcompartments

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/3D_subcompartments/")

subcomp_A <- read.table("GSE63525_GM12878_subcompartments_A.bed",header = FALSE, sep="\t")
subcomp_B <- read.table("GSE63525_GM12878_subcompartments_B.bed",header = FALSE, sep="\t")

A_gr <- GRanges(seqnames=subcomp_A$V1,IRanges(start=subcomp_A$V2, end=subcomp_A$V3))
A_gr_center <- resize(A_gr, width = 1, fix = "center")
B_gr <- GRanges(seqnames=subcomp_B$V1,IRanges(start=subcomp_B$V2, end=subcomp_B$V3))
B_gr_center <- resize(B_gr, width = 1, fix = "center")

  gm12878_10kb$A <- ifelse(countOverlaps(binslist10,A_gr)>=1,1,0)
  gm12878_10kb$B <- ifelse(countOverlaps(binslist10,B_gr)>=1,1,0)
  
  gm12878_10kb$A_dist <- mcols(distanceToNearest(binslist10_center, A_gr_center))$distance
  gm12878_10kb$B_dist <- mcols(distanceToNearest(binslist10_center, B_gr_center))$distance
  
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

  gm12878_10kb$complex <- ifelse(countOverlaps(binslist10,complex_gr)>=1,1,0)
  gm12878_10kb$deletion <- ifelse(countOverlaps(binslist10,deletion_gr)>=1,1,0)
  gm12878_10kb$duplication <- ifelse(countOverlaps(binslist10,duplication_gr)>=1,1,0)
  gm12878_10kb$gain_loss <- ifelse(countOverlaps(binslist10,gain_loss_gr)>=1,1,0)
  gm12878_10kb$insertion <- ifelse(countOverlaps(binslist10,insertion_gr)>=1,1,0)
  gm12878_10kb$inversion <- ifelse(countOverlaps(binslist10,inversion_gr)>=1,1,0)
  gm12878_10kb$mobile_element_insertion <- ifelse(countOverlaps(binslist10,mobile_element_insertion_gr)>=1,1,0)
  gm12878_10kb$novel_sequence_insertion <- ifelse(countOverlaps(binslist10,novel_sequence_insertion_gr)>=1,1,0)
  gm12878_10kb$sequence_alteration <- ifelse(countOverlaps(binslist10,sequence_alteration_gr)>=1,1,0)
  
  gm12878_10kb$complex_dist <- mcols(distanceToNearest(binslist10_center, complex_gr_center))$distance
  gm12878_10kb$deletion_dist <- mcols(distanceToNearest(binslist10_center, deletion_gr_center))$distance
  gm12878_10kb$duplication_dist <- mcols(distanceToNearest(binslist10_center, duplication_gr_center))$distance
  gm12878_10kb$gain_loss_dist <- mcols(distanceToNearest(binslist10_center, gain_loss_gr_center))$distance
  gm12878_10kb$insertion_dist <- mcols(distanceToNearest(binslist10_center, insertion_gr_center))$distance
  gm12878_10kb$inversion_dist <- mcols(distanceToNearest(binslist10_center, inversion_gr_center))$distance
  gm12878_10kb$mobile_element_insertion_dist <- mcols(distanceToNearest(binslist10_center, mobile_element_insertion_gr_center))$distance
  gm12878_10kb$novel_sequence_insertion_dist <- mcols(distanceToNearest(binslist10_center, novel_sequence_insertion_gr_center))$distance
  gm12878_10kb$sequence_alteration_dist <- mcols(distanceToNearest(binslist10_center, sequence_alteration_gr_center))$distance
  
  
## GERP

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/GERP/")

gerp <- read.table("GERP_hg19.bed",header=FALSE,sep="\t")

gerp_gr <- GRanges(seqnames=gerp$V1,IRanges(start=gerp$V2,end=gerp$V3))
mcols(gerp_gr)$score <- gerp$V5
gerp_gr_center <- resize(gerp_gr, width = 1, fix = "center")

  gm12878_10kb$gerp <- ifelse(countOverlaps(binslist10,gerp_gr)>=1,1,0)
  
  gm12878_10kb$gerp_dist <- mcols(distanceToNearest(binslist10_center, gerp_gr_center))$distance
  
#finding which flanks overlap the gerp file so that we can add a score variable
#all other flanks will have a score of 0
#which(gm12878_10kb$gerp==1)
gm12878_10kb$gerp_score <- 0
gerpoverlap <- findOverlaps(binslist10,gerp_gr)
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

  gm12878_10kb$DNA <- ifelse(countOverlaps(binslist10,DNA_gr)>=1,1,0)
  gm12878_10kb$line <- ifelse(countOverlaps(binslist10,line_gr)>=1,1,0)
  gm12878_10kb$low_complexity <- ifelse(countOverlaps(binslist10,low_complexity_gr)>=1,1,0)
  gm12878_10kb$LTR <- ifelse(countOverlaps(binslist10,LTR_gr)>=1,1,0)
  gm12878_10kb$other <- ifelse(countOverlaps(binslist10,other_gr)>=1,1,0)
  gm12878_10kb$RC <- ifelse(countOverlaps(binslist10,RC_gr)>=1,1,0)
  gm12878_10kb$satellite <- ifelse(countOverlaps(binslist10,satellite_gr)>=1,1,0)
  gm12878_10kb$simple_repeat <- ifelse(countOverlaps(binslist10,simple_repeat_gr)>=1,1,0)
  gm12878_10kb$SINE <- ifelse(countOverlaps(binslist10,SINE_gr)>=1,1,0)
  
  
  gm12878_10kb$DNA_dist <- mcols(distanceToNearest(binslist10_center, DNA_gr_center))$distance 
  gm12878_10kb$line_dist <- mcols(distanceToNearest(binslist10_center, line_gr_center))$distance 
  gm12878_10kb$low_complexity_dist <- mcols(distanceToNearest(binslist10_center, low_complexity_gr_center))$distance 
  gm12878_10kb$LTR_dist <- mcols(distanceToNearest(binslist10_center, LTR_gr_center))$distance 
  gm12878_10kb$other_dist <- mcols(distanceToNearest(binslist10_center, other_gr_center))$distance 
  gm12878_10kb$RC_dist <- mcols(distanceToNearest(binslist10_center, RC_gr_center))$distance 
  gm12878_10kb$satellite_dist <- mcols(distanceToNearest(binslist10_center, satellite_gr_center))$distance 
  gm12878_10kb$simple_repeat_dist <- mcols(distanceToNearest(binslist10_center, simple_repeat_gr_center))$distance 
  gm12878_10kb$SINE_dist <- mcols(distanceToNearest(binslist10_center, SINE_gr_center))$distance 
  

## super_enhancers

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/super_enhancers/")

temp = list.files()

se_GM12878 <- read.table(temp[1],header=FALSE,sep="\t")

se_GM12878_gr <- GRanges(seqnames=se_GM12878$V1,IRanges(start=se_GM12878$V2,end=se_GM12878$V3))
se_GM12878_gr_center <- resize(se_GM12878_gr, width = 1, fix = "center")

  gm12878_10kb$se_GM12878 <- ifelse(countOverlaps(binslist10,se_GM12878_gr)>=1,1,0)
  
  gm12878_10kb$se_GM12878_dist <- mcols(distanceToNearest(binslist10_center, se_GM12878_gr_center))$distance
  

## VMR

setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/VMRs/")

VMR <- read.table("VMR_hg19.bed",header=FALSE,sep="\t")

VMR_gr <- GRanges(seqnames=VMR$V1,IRanges(start=VMR$V2,end=VMR$V3))
VMR_gr_center <- resize(VMR_gr, width = 1, fix = "center")

  gm12878_10kb$VMR <- ifelse(countOverlaps(binslist10,VMR_gr)>=1,1,0)
  
  gm12878_10kb$VMR_dist <- mcols(distanceToNearest(binslist10_center, VMR_gr_center))$distance

  
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

  gm12878_10kb$Gm12878_TxnElongation <- ifelse(countOverlaps(binslist10,Gm12878_TxnElongation_gr)>=1,1,0)
  gm12878_10kb$Gm12878_WeakTxn <- ifelse(countOverlaps(binslist10,Gm12878_WeakTxn_gr)>=1,1,0)
  gm12878_10kb$Gm12878_Repressed <- ifelse(countOverlaps(binslist10,Gm12878_Repressed_gr)>=1,1,0)
  gm12878_10kb$Gm12878_Heterochromlo <- ifelse(countOverlaps(binslist10,Gm12878_Heterochromlo_gr)>=1,1,0)
  gm12878_10kb$Gm12878_RepetitiveCNV14 <- ifelse(countOverlaps(binslist10,Gm12878_RepetitiveCNV14_gr)>=1,1,0)
  gm12878_10kb$Gm12878_RepetitiveCNV15 <- ifelse(countOverlaps(binslist10,Gm12878_RepetitiveCNV15_gr)>=1,1,0)
  gm12878_10kb$Gm12878_ActivePromoter <- ifelse(countOverlaps(binslist10,Gm12878_ActivePromoter_gr)>=1,1,0)
  gm12878_10kb$Gm12878_WeakPromoter <- ifelse(countOverlaps(binslist10,Gm12878_WeakPromoter_gr)>=1,1,0)
  gm12878_10kb$Gm12878_PoisedPromoter <- ifelse(countOverlaps(binslist10,Gm12878_PoisedPromoter_gr)>=1,1,0)
  gm12878_10kb$Gm12878_StrongEnhancer4 <- ifelse(countOverlaps(binslist10,Gm12878_StrongEnhancer4_gr)>=1,1,0)
  gm12878_10kb$Gm12878_StrongEnhancer5 <- ifelse(countOverlaps(binslist10,Gm12878_StrongEnhancer5_gr)>=1,1,0)
  gm12878_10kb$Gm12878_WeakEnhancer6 <- ifelse(countOverlaps(binslist10,Gm12878_WeakEnhancer6_gr)>=1,1,0)
  gm12878_10kb$Gm12878_WeakEnhancer7 <- ifelse(countOverlaps(binslist10,Gm12878_WeakEnhancer7_gr)>=1,1,0)
  gm12878_10kb$Gm12878_Insulator <- ifelse(countOverlaps(binslist10,Gm12878_Insulator_gr)>=1,1,0)
  gm12878_10kb$Gm12878_TxnTransition <- ifelse(countOverlaps(binslist10,Gm12878_TxnTransition_gr)>=1,1,0)

  gm12878_10kb$Gm12878_TxnElongation_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_TxnElongation_gr_center))$distance
  gm12878_10kb$Gm12878_WeakTxn_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_WeakTxn_gr_center))$distance
  gm12878_10kb$Gm12878_Repressed_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_Repressed_gr_center))$distance
  gm12878_10kb$Gm12878_Heterochromlo_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_Heterochromlo_gr_center))$distance
  gm12878_10kb$Gm12878_RepetitiveCNV14_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_RepetitiveCNV14_gr_center))$distance
  gm12878_10kb$Gm12878_RepetitiveCNV15_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_RepetitiveCNV15_gr_center))$distance
  gm12878_10kb$Gm12878_ActivePromoter_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_ActivePromoter_gr_center))$distance
  gm12878_10kb$Gm12878_WeakPromoter_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_WeakPromoter_gr_center))$distance
  gm12878_10kb$Gm12878_PoisedPromoter_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_PoisedPromoter_gr_center))$distance
  gm12878_10kb$Gm12878_StrongEnhancer4_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_StrongEnhancer4_gr_center))$distance
  gm12878_10kb$Gm12878_StrongEnhancer5_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_StrongEnhancer5_gr_center))$distance
  gm12878_10kb$Gm12878_WeakEnhancer6_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_WeakEnhancer6_gr_center))$distance
  gm12878_10kb$Gm12878_WeakEnhancer7_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_WeakEnhancer7_gr_center))$distance
  gm12878_10kb$Gm12878_Insulator_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_Insulator_gr_center))$distance
  gm12878_10kb$Gm12878_TxnTransition_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_TxnTransition_gr_center))$distance

  
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

  gm12878_10kb$Gm12878_CTCF <- ifelse(countOverlaps(binslist10,Gm12878_CTCF_gr)>=1,1,0) 
  gm12878_10kb$Gm12878_E <- ifelse(countOverlaps(binslist10,Gm12878_E_gr)>=1,1,0)
  gm12878_10kb$Gm12878_PF <- ifelse(countOverlaps(binslist10,Gm12878_PF_gr)>=1,1,0)
  gm12878_10kb$Gm12878_R <- ifelse(countOverlaps(binslist10,Gm12878_R_gr)>=1,1,0)
  gm12878_10kb$Gm12878_T <- ifelse(countOverlaps(binslist10,Gm12878_T_gr)>=1,1,0)
  gm12878_10kb$Gm12878_TSS <- ifelse(countOverlaps(binslist10,Gm12878_TSS_gr)>=1,1,0)
  gm12878_10kb$Gm12878_WE <- ifelse(countOverlaps(binslist10,Gm12878_WE_gr)>=1,1,0)

  gm12878_10kb$Gm12878_CTCF_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_CTCF_gr_center))$distance
  gm12878_10kb$Gm12878_E_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_E_gr_center))$distance
  gm12878_10kb$Gm12878_PF_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_PF_gr_center))$distance
  gm12878_10kb$Gm12878_R_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_R_gr_center))$distance
  gm12878_10kb$Gm12878_T_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_T_gr_center))$distance
  gm12878_10kb$Gm12878_TSS_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_TSS_gr_center))$distance
  gm12878_10kb$Gm12878_WE_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_WE_gr_center))$distance

## DNase I


setwd("/home/stilianoudakisc/TAD_data_analysis/annotations/DNaseI/")

temp <- list.files()

Gm12878_DNaseI <- read.table(temp[1],header=FALSE,sep="\t") 

Gm12878_DNaseI_gr <- GRanges(seqnames=Gm12878_DNaseI$V1,IRanges(start=Gm12878_DNaseI$V2,end=Gm12878_DNaseI$V3))
Gm12878_DNaseI_gr_center <- resize(Gm12878_DNaseI_gr, width = 1, fix = "center")

gm12878_10kb$Gm12878_DNaseI <- ifelse(countOverlaps(binslist10,Gm12878_DNaseI_gr)>=1,1,0)

gm12878_10kb$Gm12878_DNaseI_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_DNaseI_gr_center))$distance


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

  gm12878_10kb$Gm12878_H2az <- ifelse(countOverlaps(binslist10,Gm12878_H2az_gr)>=1,1,0) 
  gm12878_10kb$Gm12878_H3k27ac <- ifelse(countOverlaps(binslist10,Gm12878_H3k27ac_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H3k27me3 <- ifelse(countOverlaps(binslist10,Gm12878_H3k27me3_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H3k36me3 <- ifelse(countOverlaps(binslist10,Gm12878_H3k36me3_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H3k4me1 <- ifelse(countOverlaps(binslist10,Gm12878_H3k4me1_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H3k4me2 <- ifelse(countOverlaps(binslist10,Gm12878_H3k4me2_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H3k4me3 <- ifelse(countOverlaps(binslist10,Gm12878_H3k4me3_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H3k79me2 <- ifelse(countOverlaps(binslist10,Gm12878_H3k79me2_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H3k9ac <- ifelse(countOverlaps(binslist10,Gm12878_H3k9ac_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H3k9me3 <- ifelse(countOverlaps(binslist10,Gm12878_H3k9me3_gr)>=1,1,0)
  gm12878_10kb$Gm12878_H4k20me1 <- ifelse(countOverlaps(binslist10,Gm12878_H4k20me1_gr)>=1,1,0)

  gm12878_10kb$Gm12878_H2az_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H2az_gr_center))$distance
  gm12878_10kb$Gm12878_H3k27ac_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k27ac_gr_center))$distance
  gm12878_10kb$Gm12878_H3k27me3_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k27me3_gr_center))$distance
  gm12878_10kb$Gm12878_H3k36me3_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k36me3_gr_center))$distance
  gm12878_10kb$Gm12878_H3k4me1_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k4me1_gr_center))$distance
  gm12878_10kb$Gm12878_H3k4me2_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k4me2_gr_center))$distance
  gm12878_10kb$Gm12878_H3k4me3_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k4me3_gr_center))$distance
  gm12878_10kb$Gm12878_H3k79me2_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k79me2_gr_center))$distance
  gm12878_10kb$Gm12878_H3k9ac_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k9ac_gr_center))$distance
  gm12878_10kb$Gm12878_H3k9me3_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H3k9me3_gr_center))$distance
  gm12878_10kb$Gm12878_H4k20me1_dist <- mcols(distanceToNearest(binslist10_center, Gm12878_H4k20me1_gr_center))$distance


## Adding Chromosome information to the data

#gm12878_10kb$CHR <- seqnames(binslist10)

#gm12878_10kb$CHR <- as.character(gm12878_10kb$CHR)


## Saving the data

saveRDS(gm12878_10kb, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/gm12878_10kb.rds")

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

#for(i in 1:length(grlist)){
#  
#  if(is.element("chrY", as.character(seqnames(grlist[[i]])))){
#    grlist[[i]] <- grlist[[i]][-which(as.character(seqnames(grlist[[i]]))=="chrX")]
#  }
#  
#  if(is.element("chrY", as.character(seqnames(grlist[[i]])))){
#    grlist[[i]] <- grlist[[i]][-which(as.character(seqnames(grlist[[i]]))=="chrY")]
#  }
#  
#  mcols(grlist[[i]])$distance <- mcols(distanceToNearest(grlist[[i]], bounds))$distance
#  mcols(grlist[[i]])$logdistance <- log(mcols(distanceToNearest(grlist[[i]], bounds))$distance+1,
#                                          base = 2)
#  
#}


saveRDS(grlist, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/grlist.rds")

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

cols <- c(grep("dist",colnames(gm12878_10kb)))
gm12878_10kb[,cols] <- apply(gm12878_10kb[,cols], 2, function(x){log(x + 1, base=2)})


##Changing binary variables to factors

cols <- c(intersect(grep("score",colnames(gm12878_10kb), invert = TRUE),
          grep("dist",colnames(gm12878_10kb), invert = TRUE)))
gm12878_10kb[,cols] <- lapply(gm12878_10kb[,cols], factor)


##Changing levels of response (y) to yes no

levels(gm12878_10kb$y) <- c("No", "Yes")

##Removing binary features

gm12878_10kb <- gm12878_10kb[,c(1, grep("dist",colnames(gm12878_10kb)))]

##Removing zero variance predictors

nzv <- nearZeroVar(gm12878_10kb[,-1], saveMetrics= TRUE)
nzvar <- rownames(nzv[nzv$nzv,])

nzvar

#gm12878_10kb_f <- gm12878_10kb[, -which(colnames(gm12878_10kb) %in% nzvar)]

gm12878_10kb_f <- gm12878_10kb

saveRDS(gm12878_10kb_f, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/gm12878_10kb_f.rds")


# Step 1 of Pipeline: 

## Evaluating class imbalance

### Testing SMOTE combinations

#### Splitting the data

set.seed(1234)
inTrainingSet <- createDataPartition(gm12878_10kb_f$y,p=.7,list=FALSE)
train <- gm12878_10kb_f[inTrainingSet,]
test <- gm12878_10kb_f[-inTrainingSet,]


#### Establishing tuning/training parameters

fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 3,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)
						   
#### function for roc curves

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

####Testing eight different combinations of perc.over/perc.under using Elastic Nets

#100/200, 200/200, 300/200, 400/200,100/300, 200/300, 300/300, 400/300

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_sm <- list(tpr <- matrix(nrow=length(test$y), 
                                 ncol=8),
                   fpr <- matrix(nrow=length(test$y), 
                                 ncol=8),
                   varimp <- matrix(nrow=dim(gm12878_10kb_f)[2]-1,
                                    ncol=8))
rownames(enetlst_sm[[3]]) <- colnames(gm12878_10kb_f)[-1]

enetperf_sm <- matrix(nrow = 17, ncol=8)
rownames(enetperf_sm) <- c("TN",
                           "FN",
                           "FP",
                           "TP",
                           "Total",
                           "Sensitivity",
                           "Specificity",
                           "Kappa",
                           "Accuracy",
                           "Precision",
                           "FPR",
                           "FNR",
                           "FOR",
                           "NPV",
                           "MCC",
                           "F1",
                           "AUC")


for(i in 1:4){
  set.seed(112)
  #100/200, 200/200, 300/200, 400/200
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 200)
  
  #ENET Model
  enetModel_sm <- train(y ~ ., data=train_smote, 
                        method = "glmnet", 
                        metric="ROC", 
                        trControl = fitControl,
                        family="binomial",
                        tuneLength=5,
                        standardize=FALSE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(enetModel_sm, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_sm[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_sm[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_sm[[3]][,i] <- varImp(enetModel_sm)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_sm,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_sm[1,i] <- TN
  enetperf_sm[2,i] <- FN
  enetperf_sm[3,i] <- FP
  enetperf_sm[4,i] <- TP
  enetperf_sm[5,i] <- sum(confMat$table)
  enetperf_sm[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_sm[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_sm[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_sm[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_sm[10,i] <- TP/(TP+FP)
  enetperf_sm[11,i] <- FP/(FP+TN)
  enetperf_sm[12,i] <- FN/(FN+TN)
  enetperf_sm[13,i] <- FN/(FN+TN)
  enetperf_sm[14,i] <- TN/(TN+FN)
  enetperf_sm[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_sm[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_sm[17,i] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  
  
  #########################################################################################
  
  set.seed(112)
  #100/300, 200/300, 300/300, 400/300
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 300)
  
  #ENET Model
  enetModel_sm <- train(y ~ ., data=train_smote, 
                        method = "glmnet", 
                        metric="ROC", 
                        trControl = fitControl,
                        family="binomial",
                        tuneLength=5,
                        standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.enetModel <- as.vector(predict(enetModel_sm, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_sm[[1]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_sm[[2]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_sm[[3]][,i+4] <- varImp(enetModel_sm)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_sm,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_sm[1,i+4] <- TN
  enetperf_sm[2,i+4] <- FN
  enetperf_sm[3,i+4] <- FP
  enetperf_sm[4,i+4] <- TP
  enetperf_sm[5,i+4] <- sum(confMat$table)
  enetperf_sm[6,i+4] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_sm[7,i+4] <- as.vector(confMat$byClass["Specificity"])
  enetperf_sm[8,i+4] <- as.vector(confMat$overall["Kappa"])
  enetperf_sm[9,i+4] <- as.vector(confMat$overall["Accuracy"])
  enetperf_sm[10,i+4] <- TP/(TP+FP)
  enetperf_sm[11,i+4] <- FP/(FP+TN)
  enetperf_sm[12,i+4] <- FN/(FN+TN)
  enetperf_sm[13,i+4] <- FN/(FN+TN)
  enetperf_sm[14,i+4] <- TN/(TN+FN)
  enetperf_sm[15,i+4] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_sm[16,i+4] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_sm[17,i+4] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  
}
	

saveRDS(enetlst_sm, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_sm_lns.rds")

saveRDS(enetperf_sm, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_sm.rds")

	
#### Evaluating SMOTE

#Plotting Performance
auc.sm <- data.frame(Combination=c("100/200","200/200","300/200","400/200",
                                     "100/300","200/300","300/300","400/300"),
                       AUC=c(enetperf_sm[17,1],enetperf_sm[17,2],enetperf_sm[17,3],enetperf_sm[17,4],
                             enetperf_sm[17,5],enetperf_sm[17,6],enetperf_sm[17,7],enetperf_sm[17,8]))

auc.sm <- auc.sm[order(auc.sm$AUC, decreasing=TRUE),]

auc.sm$Combination <- factor(auc.sm$Combination, levels=auc.sm$Combination)

auc.sm

saveRDS(auc.sm, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/auc.sm.rds")


auc.sm.p<-ggplot(data=auc.sm, aes(x=Combination, y=AUC, fill=Combination)) + 
  xlab("Sampling Combination") + ylab("AUC") +
  geom_bar(stat="identity") + ylim(0,1) +
  scale_fill_manual(values=gray(seq(0,.7,.1)), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Sampling Combinations using SMOTE")

auc.sm.p
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/auc.sm.p.png")

onetwo <- data.frame(fpr=enetlst_sm[[2]][,1],tpr=enetlst_sm[[1]][,1], Combo = "100/200");
twotwo <- data.frame(fpr=enetlst_sm[[2]][,2],tpr=enetlst_sm[[1]][,2], Combo = "200/200");
threetwo <- data.frame(fpr=enetlst_sm[[2]][,3],tpr=enetlst_sm[[1]][,3], Combo = "300/200");
fourtwo <- data.frame(fpr=enetlst_sm[[2]][,4],tpr=enetlst_sm[[1]][,4], Combo = "400/200");
onethree <- data.frame(fpr=enetlst_sm[[2]][,5],tpr=enetlst_sm[[1]][,5], Combo = "100/300");
twothree <- data.frame(fpr=enetlst_sm[[2]][,6],tpr=enetlst_sm[[1]][,6], Combo = "200/300");
threethree <- data.frame(fpr=enetlst_sm[[2]][,7],tpr=enetlst_sm[[1]][,7], Combo = "300/300");
fourthree <- data.frame(fpr=enetlst_sm[[2]][,8],tpr=enetlst_sm[[1]][,8], Combo = "400/300")

allrocdat <- rbind.data.frame(onetwo,
                              twotwo,
                              threetwo,
                              fourtwo,
                              onethree,
                              twothree,
                              threethree,
                              fourthree)

roc.sm <- ggplot(data=allrocdat, aes(x=fpr, y=tpr, color=Combo)) + 
  geom_line(size=1) +
  scale_colour_manual(name="Combination",
    labels=c("100/200", 
             "200/200",
             "300/200",
             "400/200",
             "100/300",
             "200/300",
             "300/300",
             "400/300"),
    values=gray(seq(0,.7,.1))) + 
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curves for Different \n Normalization Techniques")

roc.sm
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.sm.png")


### Elastic Net with balanced classes using bootstrap sampling 

####set number of bootstrap samples

bootsamps = 5

####create a matrix of row ids that represent the zero class

sampids <- matrix(ncol=bootsamps, 
                  nrow=length(gm12878_10kb_f$y[which(gm12878_10kb_f$y=="Yes")]))
				  
####filling in the sample ids matrix

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(gm12878_10kb_f$y=="No"),
                        length(which(gm12878_10kb_f$y=="Yes")),
                        replace = TRUE)
}


#### Elastic-Net With log transform and no standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_b <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   varimp <- matrix(nrow=dim(gm12878_10kb_f)[2]-1,
                                    ncol=bootsamps))
rownames(enetlst_b[[3]]) <- colnames(gm12878_10kb_f)[-1]


enetperf_b <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_b) <- c("TN",
                          "FN",
                          "FP",
                          "TP",
                          "Total",
                          "Sensitivity",
                          "Specificity",
                          "Kappa",
                          "Accuracy",
                          "Precision",
                          "FPR",
                          "FNR",
                          "FOR",
                          "NPV",
                          "MCC",
                          "F1",
                          "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
                           gm12878_10kb_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  enetModel_b <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(enetModel_b, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_b[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_b[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_b[[3]][,i] <- varImp(enetModel_b)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_b,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_b[1,i] <- TN
  enetperf_b[2,i] <- FN
  enetperf_b[3,i] <- FP
  enetperf_b[4,i] <- TP
  enetperf_b[5,i] <- sum(confMat$table)
  enetperf_b[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_b[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_b[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_b[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_b[10,i] <- TP/(TP+FP)
  enetperf_b[11,i] <- FP/(FP+TN)
  enetperf_b[12,i] <- FN/(FN+TN)
  enetperf_b[13,i] <- FN/(FN+TN)
  enetperf_b[14,i] <- TN/(TN+FN)
  enetperf_b[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_b[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_b[17,i] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  
}

saveRDS(enetlst_b,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_b.rds")
saveRDS(enetperf_b,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_b.rds")
 

#### Evaluating Bootstraps

#Mean AUC across n bootstrap samples
mean(enetperf_b[17,])

auc.bs <- round(mean(enetperf_b[17,]),3)
auc.bs

#roc curve
fpr.bs <- rowMeans(enetlst_b[[2]])
tpr.bs <- rowMeans(enetlst_b[[1]])
rocdat.bs <- data.frame(fpr=fpr.bs, tpr=tpr.bs)

roc.bs <- ggplot(rocdat.bs, aes(x=fpr, y=tpr)) + 
  geom_line(size=1, color="black") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve for Balanced Classes \n Using n Bootstrap Samples")

roc.bs
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.bs.png")


### Comparing additional performance metrics between SMOTE and Bootstrap methods

options(scipen = 999)

enetperf_sm <- round(enetperf_sm,2)
enetperf_b <- round(as.matrix(rowMeans(enetperf_b)),2)
enetperf_b[1:5,1] <- round(enetperf_b[1:5,1],0)

perfdat_sm_bs <- cbind.data.frame(rownames(enetperf_b), 
                            enetperf_sm,
                            enetperf_b)
rownames(perfdat_sm_bs) <- NULL
colnames(perfdat_sm_bs) <- c("Metric","100/200", "200/200", "300/200", "400/200", 
                       "100/300", "200/300", "300/300", "400/300",
                       "Bootstraps")

#kable(perfdat_sm_bs)
saveRDS("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/perfdat_sm_bs.rds")

mccf1 <- data.frame(Metric = c(rep("MCC",9), rep("F1",9)),
                    Technique = rep(c("100/200", 
                    "200/200", 
                    "300/200", 
                    "400/200", 
                    "100/300",
                    "200/300", 
                    "300/300", 
                    "400/300",
                    "Bootstraps"), 2),
                    Value = c(as.numeric(perfdat[15,2:10]), as.numeric(perfdat[16,2:10])))

mcc_f1_sm_bs <- ggplot(data=mccf1, aes(x=Technique, y=Value, fill=Metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c('black','lightgray')) +
  xlab("Balancing Technique") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mcc_f1_sm_bs
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mcc_f1_sm_bs.png")

mcc_sm_bs <-ggplot(data=mccf1[1:9,], aes(x=Technique, y=Value, fill=Technique)) + 
  xlab("Balancing Technique") + ylab("MCC") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=gray(rev(seq(0,.8,.1))), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Class Balancing Techniques")
mcc_sm_bs
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mcc_sm_bs.png")


f1_sm_bs <-ggplot(data=mccf1[10:18,], aes(x=Technique, y=Value, fill=Technique)) + 
  xlab("Balancing Technique") + ylab("F1") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=gray(rev(seq(0,.8,.1))), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Class Balancing Techniques")
f1_sm_bs
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/f1_sm_bs.png")



# Comparing 100/200 SMOTE with Bootstrapped model

#roc curves
roc.100200.bs <- ggplot() + 
  geom_line(aes(fpr, tpr, colour=gray(.7)[1]), rocdat.bs) +
  geom_line(aes(fpr, tpr, colour="black"), onetwo) + 
  scale_colour_manual(name="Sampling \n Technique",
    labels=c("Bootstrap","SMOTE: \n 100/200"),
    values=c("black",gray(.7))) +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("100 Bootstrap Samples vs 100/200 SMOTE")
roc.100200.bs
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.100200.bs.png")

  

#variable importance plots
varimp.bs <- as.vector(rowMeans(enetlst_bs[[4]]))
Labels <- rownames(enetlst_bs[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.bs.df <- data.frame(Feature=Labels,
                                 Importance=varimp.bs)
varimp.bs.df <- varimp.bs.df[order(varimp.bs.df$Importance),]
varimp.bs.df <- varimp.bs.df[(dim(varimp.bs.df)[1]-19):dim(varimp.bs.df)[1],]
varimp.bs.df$Feature <- factor(varimp.bs.df$Feature,
                                     levels=varimp.bs.df$Feature)
p.bs <- ggplot(varimp.bs.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="black") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Variable Importance Plot: \n 100 Bootstrap Samples")
p.bs
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/p.bs.png")


varimp.sm <- as.vector(enetlst_sm[[4]][,1])
Labels <- names(enetlst_sm[[4]][,1])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.sm.df <- data.frame(Feature=Labels,
                                 Importance=varimp.sm)
varimp.sm.df <- varimp.sm.df[order(varimp.sm.df$Importance),]
varimp.sm.df <- varimp.sm.df[(dim(varimp.sm.df)[1]-19):dim(varimp.sm.df)[1],]
varimp.sm.df$Feature <- factor(varimp.sm.df$Feature,
                                     levels=varimp.sm.df$Feature)
p.sm <- ggplot(varimp.sm.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill=gray(.7)) +
  coord_flip() +
  theme_minimal() +
  ggtitle("Variable Importance Plot: \n 100/200 SMOTE")
p.sm
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/p.sm.png")


# Step 2 of Pipeline: 

## Evaluating Normalization Techniques

### Setting bootstraps, tuning parameters, etc.

#these are the same as for evaluating class imbalance

####With log transform and standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_ls <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   varimp <- matrix(nrow=dim(gm12878_10kb_f)[2]-1,
                                    ncol=bootsamps))
rownames(enetlst_ls[[3]]) <- colnames(gm12878_10kb_f)[-1]

enetperf_ls <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_ls) <- c("TN",
                           "FN",
                           "FP",
                           "TP",
                           "Total",
                           "Sensitivity",
                           "Specificity",
                           "Kappa",
                           "Accuracy",
                           "Precision",
                           "FPR",
                           "FNR",
                           "FOR",
                           "NPV",
                           "MCC",
                           "F1",
                           "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
                           gm12878_10kb_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=TRUE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_ls[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_ls[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_ls[[3]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_ls[1,i] <- TN
  enetperf_ls[2,i] <- FN
  enetperf_ls[3,i] <- FP
  enetperf_ls[4,i] <- TP
  enetperf_ls[5,i] <- sum(confMat$table)
  enetperf_ls[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_ls[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_ls[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_ls[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_ls[10,i] <- TP/(TP+FP)
  enetperf_ls[11,i] <- FP/(FP+TN)
  enetperf_ls[12,i] <- FN/(FN+TN)
  enetperf_ls[13,i] <- FN/(FN+TN)
  enetperf_ls[14,i] <- TN/(TN+FN)
  enetperf_ls[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_ls[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_ls[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  
}

saveRDS(enetlst_ls,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_ls.rds")
saveRDS(enetperf_ls,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_ls.rds")


####With log transform and no standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_lns <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    varimp <- matrix(nrow=dim(gm12878_10kb_f)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_lns[[3]]) <- colnames(chr1_gm12878_f)[-1]

enetperf_lns <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_lns) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1",
                            "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
                           gm12878_10kb_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_lns[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_lns[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_lns[[3]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_lns[1,i] <- TN
  enetperf_lns[2,i] <- FN
  enetperf_lns[3,i] <- FP
  enetperf_lns[4,i] <- TP
  enetperf_lns[5,i] <- sum(confMat$table)
  enetperf_lns[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_lns[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_lns[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_lns[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_lns[10,i] <- TP/(TP+FP)
  enetperf_lns[11,i] <- FP/(FP+TN)
  enetperf_lns[12,i] <- FN/(FN+TN)
  enetperf_lns[13,i] <- FN/(FN+TN)
  enetperf_lns[14,i] <- TN/(TN+FN)
  enetperf_lns[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_lns[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_lns[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  
}

saveRDS(enetlst_lns,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_lns.rds")
saveRDS(enetperf_lns,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_lns.rds")


####Without log transform and standardization

gm12878_10kb_f_nl <- gm12878_10kb_f

cols <- c(grep("dist",colnames(gm12878_10kb_f_nl)))
gm12878_10kb_f_nl[,cols] <- apply(gm12878_10kb_f_nl[,cols], 2, function(x){2^x})

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_nls <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f_nl$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f_nl$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    varimp <- matrix(nrow=dim(gm12878_10kb_f_nl)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_nls[[3]]) <- colnames(gm12878_10kb_f_nl)[-1]

enetperf_nls <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_nls) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1",
                            "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f_nl[which(gm12878_10kb_f_nl$y=="Yes"),],
                           gm12878_10kb_f_nl[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=TRUE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_nls[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_nls[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_nls[[3]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_nls[1,i] <- TN
  enetperf_nls[2,i] <- FN
  enetperf_nls[3,i] <- FP
  enetperf_nls[4,i] <- TP
  enetperf_nls[5,i] <- sum(confMat$table)
  enetperf_nls[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_nls[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_nls[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_nls[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_nls[10,i] <- TP/(TP+FP)
  enetperf_nls[11,i] <- FP/(FP+TN)
  enetperf_nls[12,i] <- FN/(FN+TN)
  enetperf_nls[13,i] <- FN/(FN+TN)
  enetperf_nls[14,i] <- TN/(TN+FN)
  enetperf_nls[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_nls[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_nls[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  
}

saveRDS(enetlst_nls,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_nls.rds")
saveRDS(enetperf_nls,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_nls.rds")


#### Without log transform and no standardization

enetlst_nlns <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f_nl$y=="Yes"))*2)*.3), 
                                   ncol=bootsamps),
                     fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f_nl$y=="Yes"))*2)*.3), 
                                   ncol=bootsamps),
                     varimp <- matrix(nrow=dim(gm12878_10kb_f_nl)[2]-1,
                                      ncol=bootsamps))
rownames(enetlst_nlns[[3]]) <- colnames(gm12878_10kb_f_nl)[-1]

enetperf_nlns <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_nlns) <- c("TN",
                             "FN",
                             "FP",
                             "TP",
                             "Total",
                             "Sensitivity",
                             "Specificity",
                             "Kappa",
                             "Accuracy",
                             "Precision",
                             "FPR",
                             "FNR",
                             "FOR",
                             "NPV",
                             "MCC",
                             "F1",
                             "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f_nl[which(gm12878_10kb_f_nl$y=="Yes"),],
                           gm12878_10kb_f_nl[sampids[,i],])
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_nlns[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_nlns[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_nlns[[4]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_nlns[1,i] <- TN
  enetperf_nlns[2,i] <- FN
  enetperf_nlns[3,i] <- FP
  enetperf_nlns[4,i] <- TP
  enetperf_nlns[5,i] <- sum(confMat$table)
  enetperf_nlns[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_nlns[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_nlns[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_nlns[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_nlns[10,i] <- TP/(TP+FP)
  enetperf_nlns[11,i] <- FP/(FP+TN)
  enetperf_nlns[12,i] <- FN/(FN+TN)
  enetperf_nlns[13,i] <- FN/(FN+TN)
  enetperf_nlns[14,i] <- TN/(TN+FN)
  enetperf_nlns[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_nlns[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_nlns[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  
}

saveRDS(enetlst_nlns,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_nlns.rds")
saveRDS(enetperf_nlns,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_nlns.rds")

#### Measuring Performance of different Normalization Techniques

# Log tranformed and standardized

auc.ls <- round(mean(enetperf_ls[17,]),3)
auc.ls

#roc curve
fpr.ls <- rowMeans(enetlst_ls[[2]])
tpr.ls <- rowMeans(enetlst_ls[[1]])
rocdat.ls <- data.frame(fpr=fpr.ls, tpr=tpr.ls)
roc.ls <- ggplot(rocdat.ls, aes(x=fpr, y=tpr)) + 
  geom_line(size=1, color="#4D4D4D") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve: Log Transformed & Standardized")

roc.ls 
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.ls.png")
  
varimp.ls <- as.vector(rowMeans(enetlst_ls[[4]]))
Labels <- rownames(enetlst_ls[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.ls.df <- data.frame(Feature=Labels,
                                 Importance=varimp.ls)
varimp.ls.df <- varimp.ls.df[order(varimp.ls.df$Importance),]
varimp.ls.df <- varimp.ls.df[(dim(varimp.ls.df)[1]-9):dim(varimp.ls.df)[1],]
varimp.ls.df$Feature <- factor(varimp.ls.df$Feature,
                                     levels=varimp.ls.df$Feature)
p.ls <- ggplot(varimp.ls.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="#4D4D4D") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Log Transformed \n & Standardized")
p.ls
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/p.ls.png")

#####################################################################

# Log tranformed and un-standardized

auc.lns <- round(mean(enetperf_lns[17,1]),3)
auc.lns

#roc curve
fpr.lns <- rowMeans(enetlst_lns[[2]])
tpr.lns <- rowMeans(enetlst_lns[[1]])
rocdat.lns <- data.frame(fpr=fpr.lns, tpr=tpr.lns)
roc.lns <- ggplot(rocdat.lns, aes(x=fpr.lns, y=tpr.lns)) + 
  geom_line(size=1, color="#000000") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve: Log Transformed & Un-Standardized")
roc.lns
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.lns")

varimp.lns <- as.vector(rowMeans(enetlst_lns[[4]]))
Labels <- rownames(enetlst_lns[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.lns.df <- data.frame(Feature=Labels,
                                 Importance=varimp.lns)
varimp.lns.df <- varimp.lns.df[order(varimp.lns.df$Importance),]
varimp.lns.df <- varimp.lns.df[(dim(varimp.lns.df)[1]-9):dim(varimp.lns.df)[1],]
varimp.lns.df$Feature <- factor(varimp.lns.df$Feature,
                                     levels=varimp.lns.df$Feature)
p.lns <- ggplot(varimp.lns.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="#000000") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Log Transformed \n & Un-Standardized")
p.lns
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/p.lns")

###################################################################

# Not Log tranformed and Standardized

auc.nls <- round(mean(enetperf_nls[17,]),3)
auc.nls

#roc curve
fpr.nls <- rowMeans(enetlst_nls[[2]])
tpr.nls <- rowMeans(enetlst_nls[[1]])
rocdat.nls <- data.frame(fpr=fpr.nls, tpr=tpr.nls)
roc.nls <- ggplot(rocdat.nls, aes(x=fpr.nls, y=tpr.nls)) + 
  geom_line(size=1, color="#E6E6E6") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve: Not Log Transformed & Standardized")
roc.nls
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.nls.png")

varimp.nls <- as.vector(rowMeans(enetlst_nls[[4]]))
Labels <- rownames(enetlst_nls[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.nls.df <- data.frame(Feature=Labels,
                                 Importance=varimp.nls)
varimp.nls.df <- varimp.nls.df[order(varimp.nls.df$Importance),]
varimp.nls.df <- varimp.nls.df[(dim(varimp.nls.df)[1]-9):dim(varimp.nls.df)[1],]
varimp.nls.df$Feature <- factor(varimp.nls.df$Feature,
                                     levels=varimp.nls.df$Feature)
p.nls <- ggplot(varimp.nls.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="#E6E6E6") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Not Log Transformed \n & Standardized")
p.nls
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/p.nls.png")

####################################################################

# Not Log tranformed and Un-Standardized

auc.nlns <- round(mean(enetperf_nlns[17,]),3)
auc.nlns

#roc curve
fpr.nlns <- rowMeans(enetlst_nlns[[2]])
tpr.nlns <- rowMeans(enetlst_nlns[[1]])
rocdat.nlns <- data.frame(fpr=fpr.nlns, tpr=tpr.nlns)
roc.nlns <- ggplot(rocdat.nlns, aes(x=fpr.nlns, y=tpr.nlns)) + 
  geom_line(size=1, color="#999999") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve: Not Log Transformed & Un-Standardized")
roc.nlns
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.nlns.png")

varimp.nlns <- as.vector(rowMeans(enetlst_nlns[[4]]))
Labels <- rownames(enetlst_nlns[[4]])
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.nlns.df <- data.frame(Feature=Labels,
                                 Importance=varimp.nlns)
varimp.nlns.df <- varimp.nlns.df[order(varimp.nlns.df$Importance),]
varimp.nlns.df <- varimp.nlns.df[(dim(varimp.nlns.df)[1]-9):dim(varimp.nlns.df)[1],]
varimp.nlns.df$Feature <- factor(varimp.nlns.df$Feature,
                                     levels=varimp.nlns.df$Feature)
p.nlns <- ggplot(varimp.nlns.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="#999999") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Not Log Transformed \n & Standardized")
p.nlns
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/p.nlns.png")


#### Comparing additional performance metrics across all normalization techniques

options(scipen = 999)

lstab <- round(as.matrix(rowMeans(enetperf_ls)),2)
lnstab <- round(as.matrix(rowMeans(enetperf_lns)),2)
nlstab <- round(as.matrix(rowMeans(enetperf_nls)),2)
nlnstab <- round(as.matrix(rowMeans(enetperf_nlns)),2)

lstab[1:5,1] <- round(lstab[1:5,1],0)
lnstab[1:5,1] <- round(lnstab[1:5,1],0)
nlstab[1:5,1] <- round(nlstab[1:5,1],0)
nlnstab[1:5,1] <- round(nlnstab[1:5,1],0)

perfdat.norm <- cbind.data.frame(rownames(enetperf_ls), 
                            lstab,
                            lnstab,
                            nlstab,
                            nlnstab)
rownames(perfdat.norm) <- NULL
colnames(perfdat.norm) <- c("Metric", "Log/Std", "Log/Un-Std", "No Log/Std", "No Log/Un-Std")

#kable(perfdat.norm)
saveRDS(perfdat.norm,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/perfdat.norm.rds")


mccf1 <- data.frame(Metric = c(rep("MCC",4), rep("F1",4)),
                    Technique = rep(c("Log/Std", 
                                      "Log/Un-Std", 
                                      "No Log/Std", 
                                      "No Log/Un-Std"), 2),
                    Value = c(as.numeric(perfdat[15,2:5]), as.numeric(perfdat[16,2:5])))

mcc_f1_norm <- ggplot(data=mccf1, aes(x=Technique, y=Value, fill=Metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c('black','lightgray')) +
  xlab("Normalization Technique") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mcc_f1_norm
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mcc_f1_norm")

mcc_norm<-ggplot(data=mccf1[1:4,], aes(x=Technique, y=Value, fill=Technique)) + 
  xlab("Normalization Technique") + ylab("MCC") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=gray(rev(c(0,.3,.6,.9))), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Class Normalization Techniques")
mcc_norm
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mcc_norm")

f1_norm<-ggplot(data=mccf1[5:8,], aes(x=Technique, y=Value, fill=Technique)) + 
  xlab("Normalization Technique") + ylab("F1") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=gray(rev(c(0,.3,.6,.9))), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Class Normalization Techniques")
f1_norm
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/f1_norm")


auc.plot <- data.frame("Normalization Technique"=c("Log/Standardaized", 
                                       "Log/Un-Standardaized",
                                       "No Log/Standardaized",
                                       "No Log/Un-Standardaized"),
                       auc=c(auc.ls,
                             auc.lns,
                             auc.nls,
                             auc.nlns))

auc.plot <- auc.plot[order(auc.plot$auc, decreasing=TRUE),]

auc.plot$Normalization.Technique <-factor(auc.plot$Normalization.Technique, 
                                     levels=auc.plot$Normalization.Technique)

auc.norm.p<-ggplot(data=auc.plot, aes(x=Normalization.Technique, y=auc, fill=Normalization.Technique)) + 
  xlab("Normalization Technique") + ylab("AUC") +
  geom_bar(stat="identity") + ylim(0,1) +
  scale_fill_manual(values=grey(c(0,.3,.6,.9)), guide=FALSE) +
  scale_x_discrete(labels= c("Log/ \n Un-Standardaized", 
             "Log/ \n Standardaized",
             "No Log/ \n Un-Standardaized",
              "No Log/ \n Standardaized")) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Normalization Techniques")
auc.norm.p
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/auc.norm.p.png")

auc.plot


rocdat.ls$Technique <- "ls"
rocdat.lns$Technique <- "lns"
rocdat.nls$Technique <- "nls"
rocdat.nlns$Technique <- "nlns"
allrocdat <- rbind.data.frame(rocdat.ls, rocdat.lns, rocdat.nls, rocdat.nlns)

roc.norm <- ggplot(data=allrocdat, aes(x=fpr, y=tpr, color=Technique)) + 
  geom_line(size=1) +
  scale_colour_manual(name="Technique",
    labels=c("Log/ \n Un-Standardaized", 
             "Log/ \n Standardaized",
             "No Log/ \n Un-Standardaized",
              "No Log/ \n Standardaized"),
    values=grey(c(0,.3,.6,.9))) + 
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curves for Different \n Normalization Techniques")
roc.norm
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.norm.png")


all.norm.p <- grid.arrange(p.ls,p.lns,p.nls,p.nlns,ncol=2)
all.norm.p
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/all.norm.p.png")


# Step 3 of Pipeline:

## Evaluating Variable Reduction Techniques

### Setting bootstraps, tuning parameters, etc.

#these are the same as for evaluating class imbalance

###Establish Null and Full models


#### Forward Selection

namesmat <- matrix(NA,nrow=dim(gm12878_10kb_f)[2]-1, ncol=bootsamps)
rownames(namesmat) <- sort(names(gm12878_10kb_f)[-which(names(gm12878_10kb_f)=="y")])

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
                           gm12878_10kb_f[sampids[,i],])
  
  #null model
  glm.null <- glm(y ~ 1, data = data, family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = data, family = binomial)
  
  best.fit.fwd = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="forward",
                      trace=0)
  
  preds <- names(best.fit.fwd$coefficients)[-1]
  preds[grep("_dist",preds,invert = TRUE)] <- unlist(lapply(preds[grep("_dist",preds,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))
  
  namesmat[which(rownames(namesmat) %in% preds),i] <- sort(preds)
  
}

namesmat <- data.frame(namesmat)

namesmat$naperc <- (rowSums(is.na(namesmat))/5)*100

fwd.preds <- rownames(namesmat)[which(namesmat$naperc >= 90)]

gm12878_10kb_fwd <- gm12878_10kb_f[,which((names(gm12878_10kb_f) %in% fwd.preds) | names(gm12878_10kb_f)=="y")]

dim(gm12878_10kb_fwd)

saveRDS(gm12878_10kb_fwd,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/gm12878_10kb_fwd.rds")


#Evaluating performance of reduced data

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_fwd <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_fwd$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_fwd$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    varimp <- matrix(nrow=dim(gm12878_10kb_fwd)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_fwd[[3]]) <- colnames(gm12878_10kb_fwd)[-1]

enetperf_fwd <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_fwd) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1",
                            "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_fwd[which(gm12878_10kb_fwd$y=="Yes"),],
                           gm12878_10kb_fwd[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_fwd[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_fwd[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_fwd[[3]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel_sm,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_fwd[1,i] <- TN
  enetperf_fwd[2,i] <- FN
  enetperf_fwd[3,i] <- FP
  enetperf_fwd[4,i] <- TP
  enetperf_fwd[5,i] <- sum(confMat$table)
  enetperf_fwd[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_fwd[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_fwd[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_fwd[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_fwd[10,i] <- TP/(TP+FP)
  enetperf_fwd[11,i] <- FP/(FP+TN)
  enetperf_fwd[12,i] <- FN/(FN+TN)
  enetperf_fwd[13,i] <- FN/(FN+TN)
  enetperf_fwd[14,i] <- TN/(TN+FN)
  enetperf_fwd[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_fwd[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_fwd[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
}

saveRDS(enetlst_fwd,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_fwd.rds")
saveRDS(enetperf_fwd,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_fwd.rds")


#### Backward Selection

namesmat <- matrix(NA,nrow=dim(gm12878_10kb_f)[2]-1, ncol=bootsamps)
rownames(namesmat) <- sort(names(gm12878_10kb_f)[-which(names(gm12878_10kb_f)=="y")])

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
                           gm12878_10kb_f[sampids[,i],])
  
  #null model
  glm.null <- glm(y ~ 1, data = data, family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = data, family = binomial)
  
  best.fit.fwd = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="backward",
                      trace=0)
  
  preds <- names(best.fit.bwd$coefficients)[-1]
  preds[grep("_dist",preds,invert = TRUE)] <- unlist(lapply(preds[grep("_dist",preds,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))
  
  namesmat[which(rownames(namesmat) %in% preds),i] <- sort(preds)
  
}

namesmat <- data.frame(namesmat)

namesmat$naperc <- (rowSums(is.na(namesmat))/5)*100

bwd.preds <- rownames(namesmat)[which(namesmat$naperc >= 90)]

gm12878_10kb_bwd <- gm12878_10kb_f[,which((names(gm12878_10kb_f) %in% bwd.preds) | names(gm12878_10kb_f)=="y")]

dim(gm12878_10kb_bwd)

saveRDS(gm12878_10kb_bwd,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/gm12878_10kb_bwd.rds")


#Evaluating performance of reduced data

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_bwd <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_bwd$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_bwd$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    varimp <- matrix(nrow=dim(gm12878_10kb_bwd)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_bwd[[3]]) <- colnames(gm12878_10kb_bwd)[-1]

enetperf_bwd <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_bwd) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1",
                            "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_bwd[which(gm12878_10kb_bwd$y=="Yes"),],
                           gm12878_10kb_bwd[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_bwd[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_bwd[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_bwd[[3]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel_sm,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_bwd[1,i] <- TN
  enetperf_bwd[2,i] <- FN
  enetperf_bwd[3,i] <- FP
  enetperf_bwd[4,i] <- TP
  enetperf_bwd[5,i] <- sum(confMat$table)
  enetperf_bwd[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_bwd[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_bwd[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_bwd[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_bwd[10,i] <- TP/(TP+FP)
  enetperf_bwd[11,i] <- FP/(FP+TN)
  enetperf_bwd[12,i] <- FN/(FN+TN)
  enetperf_bwd[13,i] <- FN/(FN+TN)
  enetperf_bwd[14,i] <- TN/(TN+FN)
  enetperf_bwd[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_bwd[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_bwd[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
}

saveRDS(enetlst_bwd,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_bwd.rds")
saveRDS(enetperf_bwd,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_bwd.rds")


#### Stepwise (Both) Selection


namesmat <- matrix(NA,nrow=dim(gm12878_10kb_f)[2]-1, ncol=bootsamps)
rownames(namesmat) <- sort(names(gm12878_10kb_f)[-which(names(gm12878_10kb_f)=="y")])

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
                           gm12878_10kb_f[sampids[,i],])
  
  #null model
  glm.null <- glm(y ~ 1, data = data, family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = data, family = binomial)
  
  best.fit.both = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="both",
                      trace=0)
  
  preds <- names(best.fit.both$coefficients)[-1]
  preds[grep("_dist",preds,invert = TRUE)] <- unlist(lapply(preds[grep("_dist",preds,invert = TRUE)], function(x){substr(x,1,nchar(x)-1)}))
  
  namesmat[which(rownames(namesmat) %in% preds),i] <- sort(preds)
  
}

namesmat <- data.frame(namesmat)

namesmat$naperc <- (rowSums(is.na(namesmat))/5)*100

both.preds <- rownames(namesmat)[which(namesmat$naperc >= 90)]

gm12878_10kb_both <- gm12878_10kb_f[,which((names(gm12878_10kb_f) %in% both.preds) | names(gm12878_10kb_f)=="y")]

dim(gm12878_10kb_both)

saveRDS(gm12878_10kb_both,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/gm12878_10kb_both.rds")


#Evaluating performance of reduced data

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_both <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_both$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_both$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    varimp <- matrix(nrow=dim(gm12878_10kb_both)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_both[[3]]) <- colnames(gm12878_10kb_both)[-1]

enetperf_both <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_both) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1",
                            "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_both[which(gm12878_10kb_both$y=="Yes"),],
                           gm12878_10kb_both[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_both[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_both[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_both[[3]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel_sm,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_both[1,i] <- TN
  enetperf_both[2,i] <- FN
  enetperf_both[3,i] <- FP
  enetperf_both[4,i] <- TP
  enetperf_both[5,i] <- sum(confMat$table)
  enetperf_both[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_both[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_both[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_both[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_both[10,i] <- TP/(TP+FP)
  enetperf_both[11,i] <- FP/(FP+TN)
  enetperf_both[12,i] <- FN/(FN+TN)
  enetperf_both[13,i] <- FN/(FN+TN)
  enetperf_both[14,i] <- TN/(TN+FN)
  enetperf_both[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_both[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_both[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
}

saveRDS(enetlst_both,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_both.rds")
saveRDS(enetperf_both,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_both.rds")


#### Recursive Feature Elimination


#setting rfe parameters
control <- rfeControl(functions=rfFuncs, method="cv", number=10)#, repeats=5)

trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummary)

namesmat <- matrix(NA,nrow=dim(gm12878_10kb_f)[2]-1, ncol=bootsamps)
rownames(namesmat) <- sort(names(gm12878_10kb_f)[-which(names(gm12878_10kb_f)=="y")])


for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
                           gm12878_10kb_f[sampids[,i],])
  rfeModel <- rfe(data[,-1], 
                data[,1], 
                sizes=c(2:50), 
                metric="ROC",
                rfeControl=control,
                trControl = trainctrl)
  
  preds <- predictors(rfeModel)

  namesmat[which(rownames(namesmat) %in% preds),i] <- sort(preds)
  
}

namesmat <- data.frame(namesmat)

namesmat$naperc <- (rowSums(is.na(namesmat))/5)*100

rfe.preds <- rownames(namesmat)[which(namesmat$naperc >= 90)]

gm12878_10kb_rfe <- gm12878_10kb_f[,which((names(gm12878_10kb_f) %in% rfe.preds) | names(gm12878_10kb_f)=="y")]

saveRDS(chr1_gm12878_rfe,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/chr1_gm12878_rfe.rds")


#Evaluating performance of reduced data

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_rfe <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_rfe$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_rfe$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    varimp <- matrix(nrow=dim(gm12878_10kb_rfe)[2]-1,
                                     ncol=bootsamps))
rownames(enetlst_rfe[[3]]) <- colnames(gm12878_10kb_rfe)[-1]

enetperf_rfe <- matrix(nrow = 17, ncol=bootsamps)
rownames(enetperf_rfe) <- c("TN",
                            "FN",
                            "FP",
                            "TP",
                            "Total",
                            "Sensitivity",
                            "Specificity",
                            "Kappa",
                            "Accuracy",
                            "Precision",
                            "FPR",
                            "FNR",
                            "FOR",
                            "NPV",
                            "MCC",
                            "F1",
                            "AUC")

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_rfe[which(gm12878_10kb_rfe$y=="Yes"),],
                           gm12878_10kb_rfe[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_rfe[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_rfe[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_rfe[[3]][,i] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel_sm,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_rfe[1,i] <- TN
  enetperf_rfe[2,i] <- FN
  enetperf_rfe[3,i] <- FP
  enetperf_rfe[4,i] <- TP
  enetperf_rfe[5,i] <- sum(confMat$table)
  enetperf_rfe[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_rfe[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_rfe[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_rfe[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_rfe[10,i] <- TP/(TP+FP)
  enetperf_rfe[11,i] <- FP/(FP+TN)
  enetperf_rfe[12,i] <- FN/(FN+TN)
  enetperf_rfe[13,i] <- FN/(FN+TN)
  enetperf_rfe[14,i] <- TN/(TN+FN)
  enetperf_rfe[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_rfe[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_rfe[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
}

saveRDS(enetlst_rfe,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetlst_rfe.rds")
saveRDS(enetperf_rfe,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/enetperf_rfe.rds")


###Measuring performance of variable reduced data


#forward

auc.fwd <- round(mean(enetperf_fwd[17,]),3)
auc.fwd

rocdat.fwd <- data.frame(sensitivity=rowMeans(enetlst_fwd[[1]]), specificity=rowMeans(enetlst_fwd[[2]]))
rocdat.fwd$Selection <- "fwd"

roc.fwd <- ggplot(rocdat.fwd, aes(x=specificity, y=sensitivity)) + 
  geom_line(size=1, color="red") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve for Forward Selection")
roc.fwd
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.fwd.png")


# Backward Selection

auc.bwd <- round(mean(enetperf_bwd[17,]),3)
auc.bwd

rocdat.bwd <- data.frame(sensitivity=rowMeans(enetlst_bwd[[1]]), specificity=rowMeans(enetlst_bwd[[2]]))
rocdat.bwd$Selection <- "bwd"

roc.bwd <- ggplot(rocdat.bwd, aes(x=specificity, y=sensitivity)) + 
  geom_line(size=1, color="blue") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve for Backward Selection")
roc.bwd
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.bwd.png")



# Both

auc.both <- round(mean(enetperf_both[17,]),3)
auc.both

rocdat.both <- data.frame(sensitivity=rowMeans(enetlst_both[[1]]), specificity=rowMeans(enetlst_both[[2]]))
rocdat.both$Selection <- "both"

roc.both <- ggplot(rocdat.both, aes(x=specificity, y=sensitivity)) + 
  geom_line(size=1, color="green") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve for Both")
roc.both 
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.both.png")



# RFE

auc.rfe <- round(mean(enetperf_rfe[17,]),3)
auc.rfe

rocdat.rfe <- data.frame(sensitivity=rowMeans(enetlst_rfe[[1]]), specificity=rowMeans(enetlst_rfe[[2]]))
rocdat.rfe$Selection <- "both"

roc.rfe <- ggplot(rocdat.rfe, aes(x=specificity, y=sensitivity)) + 
  geom_line(size=1, color="green") +
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curve for Both")
roc.rfe
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.rfe.png")

 

options(scipen = 999)

fwdtab <- round(as.matrix(rowMeans(enetperf_fwd)),2)
bwdtab <- round(as.matrix(rowMeans(enetperf_bwd)),2)
bothtab <- round(as.matrix(rowMeans(enetperf_both)),2)
enetperf_rfe <- round(enetperf_rfe,2)

fwdtab[1:5,1] <- round(fwdtab[1:5,1],0)
bwdtab[1:5,1] <- round(bwdtab[1:5,1],0)
bothtab[1:5,1] <- round(bothtab[1:5,1],0)
enetperf_rfe[1:5,] <- round(enetperf_rfe[1:5,],0)

perfdat.var.sel <- cbind.data.frame(rownames(enetperf_fwd), 
                            fwdtab,
                            bwdtab,
                            bothtab,
                            enetperf_rfe)
rownames(perfdat.var.sel) <- NULL
colnames(perfdat.var.sel) <- c("Metric", "Forward", "Backward", "Both", "RFE")

#kable(perfdat.var.sel)

saveRDS(perfdat.var.sel,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/perfdat.var.sel.rds")


mccf1 <- data.frame(Metric = c(rep("MCC",4), rep("F1",4)),
                    Technique = rep(c("Forward", 
                                      "Backward", 
                                      "Both", 
                                      "RFE"), 2),
                    Value = c(as.numeric(perfdat[15,2:5]), as.numeric(perfdat[16,2:5])))

mcc_f1_var_sel <- ggplot(data=mccf1, aes(x=Technique, y=Value, fill=Metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c('black','lightgray')) +
  xlab("Variable Selection \n Technique") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mcc_f1_var_sel
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mcc_f1_var_sel.png")



mcc_var_sel<-ggplot(data=mccf1[1:4,], aes(x=Technique, y=Value, fill=Technique)) + 
  xlab("Variable Selection \n Technique") + ylab("MCC") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=gray(rev(c(0,.3,.6,.9))), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Class Normalization Techniques")
mcc_var_sel
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mcc_var_sel.png")


f1_var_sel<-ggplot(data=mccf1[5:8,], aes(x=Technique, y=Value, fill=Technique)) + 
  xlab("Variable Selection \n Technique") + ylab("F1") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=gray(rev(c(0,.3,.6,.9))), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Class Normalization Techniques")
f1_var_sel
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/f1_var_sel.png")



auc.plot <- data.frame(Selection=c("Forward", 
                                       "Backward",
                                       "Both",
                                       "RFE"),
                       auc=c(auc.fwd,
                             auc.bwd,
                             auc.both,
                             auc.rfe))

auc.plot <- auc.plot[order(auc.plot$auc, decreasing=TRUE),]

auc.plot$Selection <-factor(auc.plot$Selection, 
                                     levels=auc.plot$Selection)

#datatable(auc.plot)
#kable(auc.plot)
auc.plot


auc.var.sel<-ggplot(data=auc.plot, aes(x=Selection, y=auc, fill=Selection)) + 
  xlab("Variable Selection Technique") + ylab("AUC") +
  geom_bar(stat="identity") + ylim(0,1) +
  scale_fill_manual(values=gray(c(0,.3,.6,.9)), guide=FALSE) +
  annotate("text", x=1, y=1, label= "28", size=6) +
  annotate("text", x=2, y=1, label= "34", size=6) +
  annotate("text", x=3, y=1, label= "26", size=6) +
  annotate("text", x=4, y=1, label= "47", size=6) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Variable Selection Techniques")
auc.var.sel
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/auc.var.sel.png")


allrocdat <- rbind.data.frame(rocdat.fwd, rocdat.bwd, rocdat.both, rocdat.rfe)

roc.var.sel <- ggplot(data=allrocdat, aes(x=specificity, y=sensitivity, color=Selection)) + 
  geom_line(size=1) +
  scale_colour_manual(name="Technique",
    labels=c("Forward", 
             "Backward",
             "Both",
              "RFE"),
    values=c("#999999", "#000000", "#4D4D4D", "#E6E6E6")) + 
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curves for Different \n Variable Selection Techniques")
roc.var.sel
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.var.sel.rds")


intersect(vars.fwd,intersect(vars.bwd,vars.both))
intersect(vars.fwd,intersect(vars.bwd,predictors(rfeModel)))
intersect(intersect(vars.fwd,intersect(vars.bwd,vars.both)),
          intersect(vars.fwd,intersect(vars.bwd,predictors(rfeModel))))

fwd <- (names(chr1_gm12878_f) %in% vars.fwd)
bwd <- (names(chr1_gm12878_f) %in% vars.bwd)
both <- (names(chr1_gm12878_f) %in% vars.both)
rfe <- (names(chr1_gm12878_f) %in% predictors(rfeModel))

#fwd compared to bwd
venndatfb <- cbind(fwd,bwd)
fb <- vennCounts(venndatfb)
vennDiagram(fb, include = "both", 
  names = c("Forward", "Backward"), 
  cex = 1, counts.col = "red")



venndat1 <- cbind(fwd,bwd,both)
venndat2 <- cbind(fwd,bwd,rfe)
venndat3 <- cbind(fwd,bwd,both,rfe)

a <- vennCounts(venndat1)
b <- vennCounts(venndat2)
c <- vennCounts(venndat3)

vennDiagram(a)
vennDiagram(b)

png(filename="/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/fwd.bwd.both.png")
vennDiagram(a, include = "both", 
  names = c("Forward", "Backward", "Both"), 
  cex = 1, counts.col = "red")
dev.off()

png(filename="/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/fwd.bwd.rfe.png")
vennDiagram(b, include = "both", 
  names = c("Forward", "Backward", "RFE"), 
  cex = 1, counts.col = "red")
dev.off()

png(filename="/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/fwd.bwd.both.rfe.png")
vennDiagram(c, include = "both", 
  names = c("Forward", "Backward","Both", "RFE"), 
  cex = 1, counts.col = "red")
dev.off()


# Final Models:

## Random Forest

### set number of bootstrap samples


bootsamps = 10

### set tuning parameters


fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

						
###create a matrix of row ids that represent the zero class

##the number of rows will match the one class
##the number of columns match the number of bootstrap samples

sampids <- matrix(ncol=bootsamps, 
                  nrow=length(gm12878_10kb_fwd$y[which(gm12878_10kb_fwd$y=="Yes")]))


###filling in the sample ids matrix

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(gm12878_10kb_fwd$y=="No"),
                        length(which(gm12878_10kb_fwd$y=="Yes")),
                        replace = TRUE)
}


###Performing model

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_fwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_fwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                varimp <- matrix(nrow=dim(gm12878_10kb_fwd)[2]-1,
                                 ncol=bootsamps))
rownames(rflst[[3]]) <- colnames(gm12878_10kb_fwd)[-1]


# Random Forest

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_fwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_fwd$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                auc <- numeric(bootsamps),
                varimp <- matrix(nrow=dim(gm12878_10kb_fwd)[2]-1,
                                 ncol=bootsamps))
rownames(rflst[[4]]) <- colnames(gm12878_10kb_fwd)[-1]

rfperf <- matrix(nrow = 17, ncol=bootsamps)
rownames(rfperf) <- c("TN",
 							"FN",
 							"FP",
 							"TP",
 							"Total",
							"Sensitivity",
 							"Specificity",
 							"Kappa",
 							"Accuracy",
 							"Precision",
 							"FPR",
 							"FNR",
 							"FOR",
 							"NPV",
 							"MCC",
 							"F1",
							"AUC")

#Performing Random Forest

for(i in 1:bootsamps){

  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_fwd[which(gm12878_10kb_fwd$y=="Yes"),],
                           gm12878_10kb_fwd[sampids[,i],])
  
  #determining the best number of variables randomly sampled as candidates at each split
  set.seed(5430)
  bestmtry <- tuneRF(data[,-1],data$y,
                   improve=.01,trace=0,plot=F) 
  bestmtry <- data.frame(bestmtry)
  bestmtry <- bestmtry[order(bestmtry$OOBError, decreasing = FALSE),]
  #bestmtry$mtry[1]

  #splitting the data
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #determining best number of trees
  tunegrid <- expand.grid(.mtry=bestmtry$mtry[1])
  modellist <- list()
    for (ntree in c(50,200,500,1000)) {
      set.seed(333)
      fit <- train(y~., data=train, 
                   method="rf", 
                   metric="Accuracy",
                   tuneGrid=tunegrid,  
                   ntree=ntree)
      key <- toString(ntree)
      modellist[[key]] <- fit
    }
    # compare results
    results <- resamples(modellist)
    #summary(results)
    #dotplot(results)
    results <- data.frame(summary(results)[3]$statistics$Accuracy)
    results <- results[order(results$Mean, decreasing = TRUE),]

  set.seed(1006)
  rfModel <- train(y~., data=train, 
                    method="rf", 
                    metric="ROC", 
                    tuneGrid=tunegrid, 
                    trControl=fitControl, 
                    ntree=as.numeric(rownames(results)[1]))

  #Prediction vector for ROC and AUC
  pred.rfModel <- as.vector(predict(rfModel, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  rflst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,1]
  rflst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,2]
  rflst[[3]][,i] <- varImp(rfModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.rfModel2 <- predict(rfModel,
                                     newdata=test,
                                     type="raw")
  confMat <- confusionMatrix(data=pred.rfModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  rfperf[1,i] <- TN
  rfperf[2,i] <- FN
  rfperf[3,i] <- FP
  rfperf[4,i] <- TP
  rfperf[5,i] <- sum(confMat$table)
  rfperf[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  rfperf[7,i] <- as.vector(confMat$byClass["Specificity"])
  rfperf[8,i] <- as.vector(confMat$overall["Kappa"])
  rfperf[9,i] <- as.vector(confMat$overall["Accuracy"])
  rfperf[10,i] <- TP/(TP+FP)
  rfperf[11,i] <- FP/(FP+TN)
  rfperf[12,i] <- FN/(FN+TN)
  rfperf[13,i] <- FN/(FN+TN)
  rfperf[14,i] <- TN/(TN+FN)
  rfperf[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  rfperf[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  rfperf[17,i] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))

}


saveRDS(rflst,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/rflst.rds")
saveRDS(rfperf,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/rfperf.rds")


# Mourad Models

##The BorderAnalysisFun function has 4 key parameters:
##1. A list of GRanges objects for each annotation
##2. A character vector listing the annotation variable names
##3. A GRanges object of the TAD data
##4. A Seqinfo object that provides information regarding the genome being used

### Coordinate Data

####1. Creating a list of GRanges for each annotation

##### 3D subcompartments

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/3D_subcompartments")

subcomp_A <- read.table("GSE63525_GM12878_subcompartments_A.BED",header = FALSE, sep="\t",stringsAsFactors = FALSE)
#subcomp_A <- subcomp_A[-which(subcomp_A$V1=="chrX" | subcomp_A$V1=="chrY"),]
subcomp_B <- read.table("GSE63525_GM12878_subcompartments_B.BED",header = FALSE, sep="\t",stringsAsFactors = FALSE)
#subcomp_B <- subcomp_B[-which(subcomp_B$V1=="chrX" | subcomp_B$V1=="chrY"),]

A_gr <- GRanges(seqnames=subcomp_A$V1,IRanges(start=subcomp_A$V2, end=subcomp_A$V3))
B_gr <- GRanges(seqnames=subcomp_B$V1,IRanges(start=subcomp_B$V2, end=subcomp_B$V3))


#####DGV

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/DGV/DGV_files")

temp = list.files()

complex <- read.table(temp[1],header=FALSE,sep="\t",stringsAsFactors = FALSE)
complex <- complex[-which(complex$V1=="chrX" | complex$V1=="chrY"),]
deletion <- read.table(temp[2],header=FALSE,sep="\t",stringsAsFactors = FALSE)
deletion <- deletion[-which(deletion$V1=="chrX" | deletion$V1=="chrY"),]
duplication <- read.table(temp[3],header=FALSE,sep="\t",stringsAsFactors = FALSE)
duplication <- duplication[-which(duplication$V1=="chrX" | duplication$V1=="chrY"),]
gain_loss <- read.table(temp[4],header=FALSE,sep="\t",stringsAsFactors = FALSE)
gain_loss <- gain_loss[-which(gain_loss$V1=="chrX" | gain_loss$V1=="chrY"),]
insertion <- read.table(temp[5],header=FALSE,sep="\t",stringsAsFactors = FALSE)
insertion <- insertion[-which(insertion$V1=="chrX" | insertion$V1=="chrY"),]
inversion <- read.table(temp[6],header=FALSE,sep="\t",stringsAsFactors = FALSE)
inversion <- inversion[-which(inversion$V1=="chrX" | inversion$V1=="chrY"),]
mobile_element_insertion <- read.table(temp[7],header=FALSE,sep="\t",stringsAsFactors = FALSE)
mobile_element_insertion <- mobile_element_insertion[-which(mobile_element_insertion$V1=="chrX" | mobile_element_insertion$V1=="chrY"),]
novel_sequence_insertion <- read.table(temp[8],header=FALSE,sep="\t",stringsAsFactors = FALSE)
novel_sequence_insertion <- novel_sequence_insertion[-which(novel_sequence_insertion$V1=="chrX" | novel_sequence_insertion$V1=="chrY"),]
sequence_alteration <- read.table(temp[9],header=FALSE,sep="\t",stringsAsFactors = FALSE)
sequence_alteration <- sequence_alteration[-which(sequence_alteration$V1=="chrX" | sequence_alteration$V1=="chrY"),]
tandem_duplication <- read.table(temp[10],header=FALSE,sep="\t",stringsAsFactors = FALSE)
tandem_duplication <- tandem_duplication[-which(tandem_duplication$V1=="chrX" | tandem_duplication$V1=="chrY"),]

complex_gr <- GRanges(seqnames=complex$V1,IRanges(start=complex$V2,end=complex$V3))
deletion_gr <- GRanges(seqnames=deletion$V1,IRanges(start=deletion$V2,end=deletion$V3))
duplication_gr <- GRanges(seqnames=duplication$V1,IRanges(start=duplication$V2,end=duplication$V3))
gain_loss_gr <- GRanges(seqnames=gain_loss$V1,IRanges(start=gain_loss$V2,end=gain_loss$V3))
insertion_gr <- GRanges(seqnames=insertion$V1,IRanges(start=insertion$V2,end=insertion$V3))
inversion_gr <- GRanges(seqnames=inversion$V1,IRanges(start=inversion$V2,end=inversion$V3))
mobile_element_insertion_gr <- GRanges(seqnames=mobile_element_insertion$V1,IRanges(start=mobile_element_insertion$V2,end=mobile_element_insertion$V3))
novel_sequence_insertion_gr <- GRanges(seqnames=novel_sequence_insertion$V1,IRanges(start=novel_sequence_insertion$V2,end=novel_sequence_insertion$V3))
sequence_alteration_gr <- GRanges(seqnames=sequence_alteration$V1,IRanges(start=sequence_alteration$V2,end=sequence_alteration$V3))
tandem_duplication_gr <- GRanges(seqnames=tandem_duplication$V1,IRanges(start=tandem_duplication$V2,end=tandem_duplication$V3))

##### GERP

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/GERP/GERP_hg19.bed")

gerp <- read.table("GERP_hg19.BED",header=FALSE,sep="\t",stringsAsFactors = FALSE)
gerp <- gerp[-which(gerp$V1=="chrX" | gerp$V1=="chrY"),]

gerp_gr <- GRanges(seqnames=gerp$V1,IRanges(start=gerp$V2,end=gerp$V3))

##### nestedRepeats

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/nestedRepeats/nestedRepeats_files")

temp = list.files()

DNA <- read.table(temp[1],header=FALSE,sep="\t",stringsAsFactors = FALSE)
DNA <- DNA[-which(DNA$V1=="chrX" | DNA$V1=="chrY"),]
line <- read.table(temp[2],header=FALSE,sep="\t",stringsAsFactors = FALSE)
line <- line[-which(line$V1=="chrX" | line$V1=="chrY"),]
low_complexity <- read.table(temp[3],header=FALSE,sep="\t",stringsAsFactors = FALSE)
low_complexity <- low_complexity[-which(low_complexity$V1=="chrX" | low_complexity$V1=="chrY"),]
LTR <- read.table(temp[4],header=FALSE,sep="\t",stringsAsFactors = FALSE)
LTR <- LTR[-which(LTR$V1=="chrX" | LTR$V1=="chrY"),]
other <- read.table(temp[5],header=FALSE,sep="\t",stringsAsFactors = FALSE)
other <- other[-which(other$V1=="chrX" | other$V1=="chrY"),]
RC <- read.table(temp[6],header=FALSE,sep="\t",stringsAsFactors = FALSE)
RC <- RC[-which(RC$V1=="chrX" | RC$V1=="chrY"),]
RNA <- read.table(temp[7],header=FALSE,sep="\t",stringsAsFactors = FALSE)
RNA <- RNA[-which(RNA$V1=="chrX" | RNA$V1=="chrY"),]
rRNA <- read.table(temp[8],header=FALSE,sep="\t",stringsAsFactors = FALSE)
rRNA <- rRNA[-which(rRNA$V1=="chrX" | rRNA$V1=="chrY"),]
satellite <- read.table(temp[9],header=FALSE,sep="\t",stringsAsFactors = FALSE)
satellite <- satellite[-which(satellite$V1=="chrX" | satellite$V1=="chrY"),]
scRNA <- read.table(temp[10],header=FALSE,sep="\t",stringsAsFactors = FALSE)
scRNA <- scRNA[-which(scRNA$V1=="chrX" | scRNA$V1=="chrY"),]
simple_repeat <- read.table(temp[11],header=FALSE,sep="\t",stringsAsFactors = FALSE)
simple_repeat <- simple_repeat[-which(simple_repeat$V1=="chrX" | simple_repeat$V1=="chrY"),]
SINE <- read.table(temp[12],header=FALSE,sep="\t",stringsAsFactors = FALSE)
SINE <- SINE[-which(SINE$V1=="chrX" | SINE$V1=="chrY"),]
snRNA <- read.table(temp[13],header=FALSE,sep="\t",stringsAsFactors = FALSE)
snRNA <- snRNA[-which(snRNA$V1=="chrX" | snRNA$V1=="chrY"),]
srpRNA <- read.table(temp[14],header=FALSE,sep="\t",stringsAsFactors = FALSE)
srpRNA <- srpRNA[-which(srpRNA$V1=="chrX" | srpRNA$V1=="chrY"),]
tRNA <- read.table(temp[15],header=FALSE,sep="\t",stringsAsFactors = FALSE)
tRNA <- tRNA[-which(tRNA$V1=="chrX" | tRNA$V1=="chrY"),]
unknown <- read.table(temp[16],header=FALSE,sep="\t",stringsAsFactors = FALSE)
unknown <- unknown[-which(unknown$V1=="chrX" | unknown$V1=="chrY"),]

DNA_gr <- GRanges(seqnames=DNA$V1,IRanges(start=DNA$V2,end=DNA$V3))
line_gr <- GRanges(seqnames=line$V1,IRanges(start=line$V2,end=line$V3))
low_complexity_gr <- GRanges(seqnames=low_complexity$V1,IRanges(start=low_complexity$V2,end=low_complexity$V3))
LTR_gr <- GRanges(seqnames=LTR$V1,IRanges(start=LTR$V2,end=LTR$V3))
other_gr <- GRanges(seqnames=other$V1,IRanges(start=other$V2,end=other$V3))
RC_gr <- GRanges(seqnames=RC$V1,IRanges(start=RC$V2,end=RC$V3))
RNA_gr <- GRanges(seqnames=RNA$V1,IRanges(start=RNA$V2,end=RNA$V3))
rRNA_gr <- GRanges(seqnames=rRNA$V1,IRanges(start=rRNA$V2,end=rRNA$V3))
satellite_gr <- GRanges(seqnames=satellite$V1,IRanges(start=satellite$V2,end=satellite$V3))
scRNA_gr <- GRanges(seqnames=scRNA$V1,IRanges(start=scRNA$V2,end=scRNA$V3))
simple_repeat_gr <- GRanges(seqnames=simple_repeat$V1,IRanges(start=simple_repeat$V2,end=simple_repeat$V3))
SINE_gr <- GRanges(seqnames=SINE$V1,IRanges(start=SINE$V2,end=SINE$V3))
snRNA_gr <- GRanges(seqnames=snRNA$V1,IRanges(start=snRNA$V2,end=snRNA$V3))
srpRNA_gr <- GRanges(seqnames=srpRNA$V1,IRanges(start=srpRNA$V2,end=srpRNA$V3))
tRNA_gr <- GRanges(seqnames=tRNA$V1,IRanges(start=tRNA$V2,end=tRNA$V3))
unknown_gr <- GRanges(seqnames=unknown$V1,IRanges(start=unknown$V2,end=unknown$V3))


##### super_enhancers

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/super_enhancers")

temp = list.files()

se_GM12878 <- read.table(temp[1],header=FALSE,sep="\t",stringsAsFactors = FALSE)
se_GM12878 <- se_GM12878[-which(se_GM12878$V1=="chrX" | se_GM12878$V1=="chrY"),]

se_GM12878_gr <- GRanges(seqnames=se_GM12878$V1,IRanges(start=se_GM12878$V2,end=se_GM12878$V3))


##### UCNE

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/UCNEs/UCNE.bed")

UCNE <- read.table("UCNE.BED",header=FALSE,sep="\t",stringsAsFactors = FALSE)
UCNE <- UCNE[-which(UCNE$V1=="chrX" | UCNE$V1=="chrY"),]

UCNE_gr <- GRanges(seqnames=UCNE$V1,IRanges(start=UCNE$V2,end=UCNE$V3))


##### VMR

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/VMRs")

VMR <- read.table("VMR_hg19.BED",header=FALSE,sep="\t",stringsAsFactors = FALSE)
VMR <- VMR[-which(VMR$V1=="chrX" | VMR$V1=="chrY"),]

VMR_gr <- GRanges(seqnames=VMR$V1,IRanges(start=VMR$V2,end=VMR$V3))


##### BroadHMM

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/BroadHmm")

temp <- list.files()

#Gm12878
Gm12878_TxnElongation <- read.table(temp[2],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_TxnElongation <- Gm12878_TxnElongation[-which(Gm12878_TxnElongation$V1=="chrX" | Gm12878_TxnElongation$V1=="chrY"),]
Gm12878_WeakTxn <- read.table(temp[3],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_WeakTxn <- Gm12878_WeakTxn[-which(Gm12878_WeakTxn$V1=="chrX" | Gm12878_WeakTxn$V1=="chrY"),]
Gm12878_Repressed <- read.table(temp[4],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_Repressed <- Gm12878_Repressed[-which(Gm12878_Repressed$V1=="chrX" | Gm12878_Repressed$V1=="chrY"),]
Gm12878_Heterochromlo <- read.table(temp[5],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_Heterochromlo <- Gm12878_Heterochromlo[-which(Gm12878_Heterochromlo$V1=="chrX" | Gm12878_Heterochromlo$V1=="chrY"),]
Gm12878_RepetitiveCNV14 <- read.table(temp[6],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_RepetitiveCNV14 <- Gm12878_RepetitiveCNV14[-which(Gm12878_RepetitiveCNV14$V1=="chrX" | Gm12878_RepetitiveCNV14$V1=="chrY"),]
Gm12878_RepetitiveCNV15 <- read.table(temp[7],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_RepetitiveCNV15 <- Gm12878_RepetitiveCNV15[-which(Gm12878_RepetitiveCNV15$V1=="chrX" | Gm12878_RepetitiveCNV15$V1=="chrY"),]
Gm12878_ActivePromoter <- read.table(temp[8],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_ActivePromoter <- Gm12878_ActivePromoter[-which(Gm12878_ActivePromoter$V1=="chrX" | Gm12878_ActivePromoter$V1=="chrY"),]
Gm12878_WeakPromoter <- read.table(temp[9],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_WeakPromoter <- Gm12878_WeakPromoter[-which(Gm12878_WeakPromoter$V1=="chrX" | Gm12878_WeakPromoter$V1=="chrY"),]
Gm12878_PoisedPromoter <- read.table(temp[10],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_PoisedPromoter <- Gm12878_PoisedPromoter[-which(Gm12878_PoisedPromoter$V1=="chrX" | Gm12878_PoisedPromoter$V1=="chrY"),]
Gm12878_StrongEnhancer4 <- read.table(temp[11],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_StrongEnhancer4 <- Gm12878_StrongEnhancer4[-which(Gm12878_StrongEnhancer4$V1=="chrX" | Gm12878_StrongEnhancer4$V1=="chrY"),]
Gm12878_StrongEnhancer5 <- read.table(temp[12],header=FALSE,sep="\t",stringsAsFactors = FALSE) 
Gm12878_StrongEnhancer5 <- Gm12878_StrongEnhancer5[-which(Gm12878_StrongEnhancer5$V1=="chrX" | Gm12878_StrongEnhancer5$V1=="chrY"),]
Gm12878_WeakEnhancer6 <- read.table(temp[13],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_WeakEnhancer6 <- Gm12878_WeakEnhancer6[-which(Gm12878_WeakEnhancer6$V1=="chrX" | Gm12878_WeakEnhancer6$V1=="chrY"),]
Gm12878_WeakEnhancer7 <- read.table(temp[14],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_WeakEnhancer7 <- Gm12878_WeakEnhancer7[-which(Gm12878_WeakEnhancer7$V1=="chrX" | Gm12878_WeakEnhancer7$V1=="chrY"),]
Gm12878_Insulator <- read.table(temp[15],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_Insulator <- Gm12878_Insulator[-which(Gm12878_Insulator$V1=="chrX" | Gm12878_Insulator$V1=="chrY"),]
Gm12878_TxnTransition <- read.table(temp[16],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_TxnTransition <- Gm12878_TxnTransition[-which(Gm12878_TxnTransition$V1=="chrX" | Gm12878_TxnTransition$V1=="chrY"),]

Gm12878_TxnElongation_gr <- GRanges(seqnames=Gm12878_TxnElongation$V1,IRanges(start=Gm12878_TxnElongation$V2,end=Gm12878_TxnElongation$V3))
Gm12878_WeakTxn_gr <- GRanges(seqnames=Gm12878_WeakTxn$V1,IRanges(start=Gm12878_WeakTxn$V2,end=Gm12878_WeakTxn$V3))
Gm12878_Repressed_gr <- GRanges(seqnames=Gm12878_Repressed$V1,IRanges(start=Gm12878_Repressed$V2,end=Gm12878_Repressed$V3))
Gm12878_Heterochromlo_gr <- GRanges(seqnames=Gm12878_Heterochromlo$V1,IRanges(start=Gm12878_Heterochromlo$V2,end=Gm12878_Heterochromlo$V3)) 
Gm12878_RepetitiveCNV14_gr <- GRanges(seqnames=Gm12878_RepetitiveCNV14$V1,IRanges(start=Gm12878_RepetitiveCNV14$V2,end=Gm12878_RepetitiveCNV14$V3)) 
Gm12878_RepetitiveCNV15_gr <- GRanges(seqnames=Gm12878_RepetitiveCNV15$V1,IRanges(start=Gm12878_RepetitiveCNV15$V2,end=Gm12878_RepetitiveCNV15$V3))
Gm12878_ActivePromoter_gr <- GRanges(seqnames=Gm12878_ActivePromoter$V1,IRanges(start=Gm12878_ActivePromoter$V2,end=Gm12878_ActivePromoter$V3))
Gm12878_WeakPromoter_gr <- GRanges(seqnames=Gm12878_WeakPromoter$V1,IRanges(start=Gm12878_WeakPromoter$V2,end=Gm12878_WeakPromoter$V3))
Gm12878_PoisedPromoter_gr <- GRanges(seqnames=Gm12878_PoisedPromoter$V1,IRanges(start=Gm12878_PoisedPromoter$V2,end=Gm12878_PoisedPromoter$V3)) 
Gm12878_StrongEnhancer4_gr <- GRanges(seqnames=Gm12878_StrongEnhancer4$V1,IRanges(start=Gm12878_StrongEnhancer4$V2,end=Gm12878_StrongEnhancer4$V3)) 
Gm12878_StrongEnhancer5_gr <- GRanges(seqnames=Gm12878_StrongEnhancer5$V1,IRanges(start=Gm12878_StrongEnhancer5$V2,end=Gm12878_StrongEnhancer5$V3))
Gm12878_WeakEnhancer6_gr <- GRanges(seqnames=Gm12878_WeakEnhancer6$V1,IRanges(start=Gm12878_WeakEnhancer6$V2,end=Gm12878_WeakEnhancer6$V3)) 
Gm12878_WeakEnhancer7_gr <- GRanges(seqnames=Gm12878_WeakEnhancer7$V1,IRanges(start=Gm12878_WeakEnhancer7$V2,end=Gm12878_WeakEnhancer7$V3)) 
Gm12878_Insulator_gr <- GRanges(seqnames=Gm12878_Insulator$V1,IRanges(start=Gm12878_Insulator$V2,end=Gm12878_Insulator$V3))
Gm12878_TxnTransition_gr <- GRanges(seqnames=Gm12878_TxnTransition$V1,IRanges(start=Gm12878_TxnTransition$V2,end=Gm12878_TxnTransition$V3)) 


##### Combined

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/Combined")

temp <- list.files()

#Gm12878
Gm12878_CTCF <- read.table(temp[2],header=FALSE,sep="\t",stringsAsFactors = FALSE) 
Gm12878_CTCF <- Gm12878_CTCF[-which(Gm12878_CTCF$V1=="chrX" | Gm12878_CTCF$V1=="chrY"),]
Gm12878_E <- read.table(temp[3],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_E <- Gm12878_E[-which(Gm12878_E$V1=="chrX" | Gm12878_E$V1=="chrY"),]
Gm12878_PF <- read.table(temp[4],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_PF <- Gm12878_PF[-which(Gm12878_PF$V1=="chrX" | Gm12878_PF$V1=="chrY"),]
Gm12878_R <- read.table(temp[5],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_R <- Gm12878_R[-which(Gm12878_R$V1=="chrX" | Gm12878_R$V1=="chrY"),]
Gm12878_T <- read.table(temp[6],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_T <- Gm12878_T[-which(Gm12878_T$V1=="chrX" | Gm12878_T$V1=="chrY"),]
Gm12878_TSS <- read.table(temp[7],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_TSS <- Gm12878_TSS[-which(Gm12878_TSS$V1=="chrX" | Gm12878_TSS$V1=="chrY"),]
Gm12878_WE <- read.table(temp[8],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_WE <- Gm12878_WE[-which(Gm12878_WE$V1=="chrX" | Gm12878_WE$V1=="chrY"),]

Gm12878_CTCF_gr <- GRanges(seqnames=Gm12878_CTCF$V1,IRanges(start=Gm12878_CTCF$V2,end=Gm12878_CTCF$V3))
Gm12878_E_gr <- GRanges(seqnames=Gm12878_E$V1,IRanges(start=Gm12878_E$V2,end=Gm12878_E$V3))
Gm12878_PF_gr <- GRanges(seqnames=Gm12878_PF$V1,IRanges(start=Gm12878_PF$V2,end=Gm12878_PF$V3))
Gm12878_R_gr <- GRanges(seqnames=Gm12878_R$V1,IRanges(start=Gm12878_R$V2,end=Gm12878_R$V3))
Gm12878_T_gr <- GRanges(seqnames=Gm12878_T$V1,IRanges(start=Gm12878_T$V2,end=Gm12878_T$V3))
Gm12878_TSS_gr <- GRanges(seqnames=Gm12878_TSS$V1,IRanges(start=Gm12878_TSS$V2,end=Gm12878_TSS$V3))
Gm12878_WE_gr <- GRanges(seqnames=Gm12878_WE$V1,IRanges(start=Gm12878_WE$V2,end=Gm12878_WE$V3))


##### DNase I

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/DNaseI")

temp <- list.files()

Gm12878_DNaseI <- read.table(temp[1],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_DNaseI <- Gm12878_DNaseI[-which(Gm12878_DNaseI$V1=="chrX" | Gm12878_DNaseI$V1=="chrY"),]

Gm12878_DNaseI_gr <- GRanges(seqnames=Gm12878_DNaseI$V1,IRanges(start=Gm12878_DNaseI$V2,end=Gm12878_DNaseI$V3))

##### Histone Modifications

setwd("C:/Users/Spiro Stilianoudakis/Documents/TAD_data/annotations/HistoneModifications")

temp <- list.files()

Gm12878_H2az <- read.table(temp[1],header=FALSE,sep="\t",stringsAsFactors = FALSE) 
Gm12878_H2az <- Gm12878_H2az[-which(Gm12878_H2az$V1=="chrX" | Gm12878_H2az$V1=="chrY"),]
Gm12878_H3k27ac <- read.table(temp[2],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k27ac <- Gm12878_H3k27ac[-which(Gm12878_H3k27ac$V1=="chrX" | Gm12878_H3k27ac$V1=="chrY"),]
Gm12878_H3k27me3 <- read.table(temp[3],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k27me3 <- Gm12878_H3k27me3[-which(Gm12878_H3k27me3$V1=="chrX" | Gm12878_H3k27me3$V1=="chrY"),]
Gm12878_H3k36me3 <- read.table(temp[4],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k36me3 <- Gm12878_H3k36me3[-which(Gm12878_H3k36me3$V1=="chrX" | Gm12878_H3k36me3$V1=="chrY"),]
Gm12878_H3k4me1 <- read.table(temp[5],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k4me1 <- Gm12878_H3k4me1[-which(Gm12878_H3k4me1$V1=="chrX" | Gm12878_H3k4me1$V1=="chrY"),] 
Gm12878_H3k4me2 <- read.table(temp[6],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k4me2 <- Gm12878_H3k4me2[-which(Gm12878_H3k4me2$V1=="chrX" | Gm12878_H3k4me2$V1=="chrY"),]
Gm12878_H3k4me3 <- read.table(temp[7],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k4me3 <- Gm12878_H3k4me3[-which(Gm12878_H3k4me3$V1=="chrX" | Gm12878_H3k4me3$V1=="chrY"),]
Gm12878_H3k79me2 <- read.table(temp[8],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k79me2 <- Gm12878_H3k79me2[-which(Gm12878_H3k79me2$V1=="chrX" | Gm12878_H3k79me2$V1=="chrY"),]
Gm12878_H3k9ac <- read.table(temp[9],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k9ac <- Gm12878_H3k9ac[-which(Gm12878_H3k9ac$V1=="chrX" | Gm12878_H3k9ac$V1=="chrY"),]
Gm12878_H3k9me3 <- read.table(temp[10],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H3k9me3 <- Gm12878_H3k9me3[-which(Gm12878_H3k9me3$V1=="chrX" | Gm12878_H3k9me3$V1=="chrY"),]
Gm12878_H4k20me1 <- read.table(temp[11],header=FALSE,sep="\t",stringsAsFactors = FALSE)
Gm12878_H4k20me1 <- Gm12878_H4k20me1[-which(Gm12878_H4k20me1$V1=="chrX" | Gm12878_H4k20me1$V1=="chrY"),]

Gm12878_H2az_gr <- GRanges(seqnames=Gm12878_H2az$V1,IRanges(start=Gm12878_H2az$V2,end=Gm12878_H2az$V3))
Gm12878_H3k27ac_gr <- GRanges(seqnames=Gm12878_H3k27ac$V1,IRanges(start=Gm12878_H3k27ac$V2,end=Gm12878_H3k27ac$V3))
Gm12878_H3k27me3_gr <- GRanges(seqnames=Gm12878_H3k27me3$V1,IRanges(start=Gm12878_H3k27me3$V2,end=Gm12878_H3k27me3$V3))
Gm12878_H3k36me3_gr <- GRanges(seqnames=Gm12878_H3k36me3$V1,IRanges(start=Gm12878_H3k36me3$V2,end=Gm12878_H3k36me3$V3))
Gm12878_H3k4me1_gr <- GRanges(seqnames=Gm12878_H3k4me1$V1,IRanges(start=Gm12878_H3k4me1$V2,end=Gm12878_H3k4me1$V3))
Gm12878_H3k4me2_gr <- GRanges(seqnames=Gm12878_H3k4me2$V1,IRanges(start=Gm12878_H3k4me2$V2,end=Gm12878_H3k4me2$V3))
Gm12878_H3k4me3_gr <- GRanges(seqnames=Gm12878_H3k4me3$V1,IRanges(start=Gm12878_H3k4me3$V2,end=Gm12878_H3k4me3$V3))
Gm12878_H3k79me2_gr <- GRanges(seqnames=Gm12878_H3k79me2$V1,IRanges(start=Gm12878_H3k79me2$V2,end=Gm12878_H3k79me2$V3))
Gm12878_H3k9ac_gr <- GRanges(seqnames=Gm12878_H3k9ac$V1,IRanges(start=Gm12878_H3k9ac$V2,end=Gm12878_H3k9ac$V3))
Gm12878_H3k9me3_gr <- GRanges(seqnames=Gm12878_H3k9me3$V1,IRanges(start=Gm12878_H3k9me3$V2,end=Gm12878_H3k9me3$V3))
Gm12878_H4k20me1_gr <- GRanges(seqnames=Gm12878_H4k20me1$V1,IRanges(start=Gm12878_H4k20me1$V2,end=Gm12878_H4k20me1$V3))


##### Concatinating all GRanges into a list

genomicFeatureList.GR <- list(A=A_gr,
B=B_gr, 
complex=complex_gr, 
deletion=deletion_gr, 
duplication=duplication_gr, 
gain_loss=gain_loss_gr, 
insertion=insertion_gr, 
inversion=inversion_gr, 
mobile_element_insertion=mobile_element_insertion_gr, 
novel_sequence_insertion=novel_sequence_insertion_gr, 
sequence_alteration=sequence_alteration_gr, 
tandem_duplication=tandem_duplication_gr, 
gerp=gerp_gr, 
DNA=DNA_gr, 
line=line_gr, 
low_complexity=low_complexity_gr, 
LTR=LTR_gr, 
other=other_gr, 
RC=RC_gr, 
RNA=RNA_gr, 
rRNA=rRNA_gr, 
satellite=satellite_gr, 
scRNA=scRNA_gr, 
simple_repeat=simple_repeat_gr, 
SINE=SINE_gr, 
snRNA=snRNA_gr, 
srpRNA=srpRNA_gr, 
tRNA=tRNA_gr, 
unknown=unknown_gr, 
se=se_GM12878_gr,  
UCNE=UCNE_gr, 
VMR=VMR_gr,
TxnElongation=Gm12878_TxnElongation_gr, 
WeakTxn=Gm12878_WeakTxn_gr, 
Repressed=Gm12878_Repressed_gr, 
Heterochromlo=Gm12878_Heterochromlo_gr,  
RepetitiveCNV14=Gm12878_RepetitiveCNV14_gr, 
RepetitiveCNV15=Gm12878_RepetitiveCNV15_gr, 
ActivePromoter=Gm12878_ActivePromoter_gr, 
WeakPromoter=Gm12878_WeakPromoter_gr, 
PoisedPromoter=Gm12878_PoisedPromoter_gr, 
StrongEnhancer4=Gm12878_StrongEnhancer4_gr,  
StrongEnhancer5=Gm12878_StrongEnhancer5_gr, 
WeakEnhancer6=Gm12878_WeakEnhancer6_gr, 
WeakEnhancer7=Gm12878_WeakEnhancer7_gr, 
Insulator=Gm12878_Insulator_gr, 
TxnTransition=Gm12878_TxnTransition_gr, 
CTCF=Gm12878_CTCF_gr, 
E=Gm12878_E_gr, 
PF=Gm12878_PF_gr, 
R=Gm12878_R_gr, 
T=Gm12878_T_gr, 
TSS=Gm12878_TSS_gr, 
WE=Gm12878_WE_gr, 
DNaseI=Gm12878_DNaseI_gr, 
H2az=Gm12878_H2az_gr, 
H3k27ac=Gm12878_H3k27ac_gr, 
H3k27me3=Gm12878_H3k27me3_gr, 
H3k36me3=Gm12878_H3k36me3_gr, 
H3k4me1=Gm12878_H3k4me1_gr,
H3k4me2=Gm12878_H3k4me2_gr, 
H3k4me3=Gm12878_H3k4me3_gr,
H3k79me2=Gm12878_H3k79me2_gr, 
H3k9ac=Gm12878_H3k9ac_gr, 
H3k9me3=Gm12878_H3k9me3_gr, 
H4k20me1=Gm12878_H4k20me1_gr) 

#removing annotations with less than 22 chromosomes
for(i in 1:length(names(genomicFeatureList.GR))){
   print(length(table(seqnames(genomicFeatureList.GR[[i]])))) 
 }
genomicFeatureList.GR <- genomicFeatureList.GR[-c(20,21,23,26,28,32)]

#removing chromosome X and Y
#chrs <- paste("chr",c(1,10:19,2,20:22,3:9),sep = "")
#for(i in 1:length(names(genomicFeatureList.GR))){
#  if(length(table(seqnames(genomicFeatureList.GR[[i]]))) > 22){
#    genomicFeatureList.GR[[i]] <- genomicFeatureList.GR[[i]][which(as.character(seqnames(genomicFeatureList.GR[[i]])) %in% chrs)]
    #levels(seqnames(genomicFeatureList.GR[[i]])) <- chrs
#  }
#}


####2. A character vector listing the annotation variable names

annotNames <- names(genomicFeatureList.GR)

####3. A GRanges object of the TAD data

#Reading in TAD data for GM12878 cell line
tad <- read.table("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.bed", stringsAsFactors = FALSE)
tad$V1 <- paste("chr",tad$V1,sep = "")
#removing X/Y chromosome
tad <- tad[-which(tad$V1=="chrX" | tad$V1=="chrY"),]
#Formating the data
tad <- tad[order(tad$V1,tad$V2),]

domains.GR <- GRanges(seqnames=tad$V1,IRanges(start=tad$V2,end=tad$V3))


####4. A Seqinfo object that provides information regarding the genome being used

seqn <- unique(tad$V1)

genome <- BSgenome.Hsapiens.UCSC.hg19
seql <- seqlengths(genome)[1:22]
seql <- seql[c(1,10:19,2,20:22,3:9)]
seql <- as.vector(seql)
isCirc <- isCircular(genome)[c(1,10:19,2,20:22,3:9)]
isCirc <- as.vector(isCirc)
seqInfoChr <- Seqinfo(seqn, seqlengths=seql, 
                      isCircular=isCirc, 
                      genome="hg19")

#providing seqlengths to each grange object in the list
for(i in 1:length(annotNames)){
  seqlengths(genomicFeatureList.GR[[i]]) <- seqlengths(seqInfoChr)
}

### Saving objects

saveRDS(genomicFeatureList.GR,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/genomicFeatureList.GR.rds")
saveRDS(annotNames,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/annotNames.rds")
saveRDS(domains.GR,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/domains.GR.rds")
saveRDS(seqInfoChr,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/seqInfoChr.rds")


## Functions associated with Mourad Model

### Function to calculate average bins


averagePerBin <- function(x, binsize, mcolname)
{
 if(!is(x, "GenomicRanges")){stop("'x' must be a GenomicRanges object")}
 if(any(is.na(seqlengths(x)))){stop("'seqlengths(x)' contains NAs")}

 bins <- IRangesList(lapply(seqlengths(x),function(seqlen)IRanges(breakInChunks(seqlen, binsize))))

 cvg <- coverage(x, weight=mcolname)

 views_list <- RleViewsList(lapply(names(cvg),function(seqname)Views(cvg[[seqname]], bins[[seqname]])))

 averageBin=NULL
 for(i in 1:length(views_list)){
  averageBin=c(averageBin,viewMeans(views_list)[[i]])
 }

 return(averageBin)
}

simple_roc <- function(labels, scores){
        labels <- labels[order(scores, decreasing=TRUE)]
        data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
      }
	  
### Function to run MLR (with or without LASSO regularization)


borderAnalysisFun<-function(genomicFeatureList.GR,GFDataType,annotNames,domains.GR,seqInfoChr,analysisMode,binSize=50,borderSize=1000,LRT=FALSE,interactionTerms="",testBorderType=FALSE,verbose=FALSE){


# CHECK INPUT DATA ------------------------------------------------------------ 

if(verbose){print("DATA CHECKING")}

if(class(genomicFeatureList.GR)!="list"){print("genomicFeatureList.GR is not a list object!"); return(0)}
for(i in 1:length(genomicFeatureList.GR)){
 if(class(genomicFeatureList.GR[[i]])!="GRanges"){print("i-th object of genomicFeatureList.GR is not a GenomicRanges object!"); return(0)}
}
if(class(GFDataType)!="character"){print("GFDataType is not a character object!"); return(0)}
if(class(annotNames)!="character"){print("annotNames is not a character object!"); return(0)}
if(class(domains.GR)!="GRanges"){print("domains.GR is not a GenomicRanges object!"); return(0)}
if(class(seqInfoChr)!="Seqinfo"){print("seqInfoChr is not a seqinfo object!");return(0)}
if(class(binSize)!="integer" & class(binSize)!="numeric"){print("binSize is not an integer or numeric object!");return(0)}
if(class(borderSize)!="integer" & class(borderSize)!="numeric"){print("borderSize is not an integer or numeric object!");return(0)}
if(class(LRT)!="logical"){print("LRT is not a logical object!"); return(0)}
if(class(analysisMode)!="character"){print("analysisMode is not a character object!"); return(0)}
if(class(interactionTerms)!="character"){print("interactionTerms is not a character object!"); return(0)}
if(class(testBorderType)!="logical"){print("testBorderType is not a character object!"); return(0)}



# PROCESS DATA ------------------------------------------------------------ 

if(verbose){print("DATA PREPROCESSING")}

chr.V=as.character(seqnames(seqInfoChr))

Borders.GR=NULL
for(i in 1:length(chr.V)){
 if(sum(seqnames(domains.GR)==chr.V[i])){
  domains.GRi=domains.GR[seqnames(domains.GR)==chr.V[i]]

  if(testBorderType){
   Domi1=as.character(domains.GRi$DomainType[1:(length(domains.GRi)-1)])
   Domi2=as.character(domains.GRi$DomainType[2:length(domains.GRi)])
   test=(Domi1>Domi2)
   Domi1flipped=Domi1
   Domi1flipped[test]=Domi2[test]
   Domi2flipped=Domi2
   Domi2flipped[test]=Domi1[test]
   BorderTypei=paste(Domi1flipped,Domi2flipped,sep="_")
  }else{
   BorderTypei=rep("NoType",length(domains.GRi)-1)
  }
  Borders.GRi=GRanges(seqnames=seqnames(domains.GRi[-1]),IRanges(start=start(domains.GRi[-1])-borderSize,end=start(domains.GRi[-1])+borderSize-1),BorderType=BorderTypei,seqinfo=seqInfoChr)
  if(i==1){
   Borders.GR=Borders.GRi
  }else{
   Borders.GR=c(Borders.GR,Borders.GRi)
  }
 }else{print(paste0("No ",chr.V[i]," in Domains.GR"))}
}

# Binned matrix
binMat=NULL
seqLengthChr=seqlengths(seqInfoChr)
for(i in 1:length(chr.V)){
 BordStarti=seq(1,seqLengthChr[i],by=binSize)
 BordEndi=BordStarti+binSize-1
 BordEndi[length(BordEndi)]=seqLengthChr[i]
 binMat=rbind(binMat,cbind(chr.V[i],BordStarti,BordEndi))
 if(verbose){print(chr.V[i])}
}
binMat.GR=GRanges(seqnames=binMat[,1],IRanges(start=as.numeric(binMat[,2]),end=as.numeric(binMat[,3])))
seqinfo(binMat.GR)=seqInfoChr
olBinBorders=findOverlaps(binMat.GR,Borders.GR)
binMat.GR$Border=rep(0,length(binMat.GR))
binMat.GR$Border[queryHits(olBinBorders)]=1
binMat.GR$BorderType=rep("NB",length(binMat.GR))
binMat.GR$BorderType[queryHits(olBinBorders)]=Borders.GR[subjectHits(olBinBorders)]$BorderType


# Annotate borders
binMat.Mat=NULL
for(i in 1:length(genomicFeatureList.GR)){
 if(GFDataType=="bed"){
  genomicFeatureList.GR[[i]]$score=rep(1,length(genomicFeatureList.GR[[i]]))
 }

 binMati=averagePerBin(genomicFeatureList.GR[[i]],binSize,"score")
 binMat.Mat=as(cBind(binMat.Mat,binMati),"Matrix")
 if(verbose){print(paste0(annotNames[i]," annotated"))}
}
borderTypeVec=levels(as.factor(binMat.GR$BorderType))
borderTypeIdx=1:length(borderTypeVec)
binMat.Mat=cBind(cbind(binMat.GR$Border,as.factor(binMat.GR$BorderType)),binMat.Mat)
colnames(binMat.Mat)<-c("Border","BorderType",annotNames)

# Compute correlations among genomic features
corGF=cor(as.matrix(binMat.Mat[,-(1:2)]))


# ENRICHMENT TEST ------------------------------------------------------------ 

if(verbose){print("DATA ANALYSIS")}

#### Analysis for each border type
matCoefMargGLM=NULL
matCoefMultGLM=NULL
matCoefMultLasso=NULL
matCoefInterGLM=NULL
matCoefInterLasso=NULL
list_MultGLM=list()
list_InterGLM=list()
for(bt in borderTypeVec[-which(borderTypeVec=="NB")]){

 # Binned matrix for borderType i
 binMat.Mati=rBind(binMat.Mat[binMat.Mat[,2]==which(borderTypeVec==bt),],binMat.Mat[binMat.Mat[,2]==which(borderTypeVec=="NB"),])
 binMat.mati=as.data.frame(as.matrix(binMat.Mati[,-2]))

 # Enrichment test
 if(sum(analysisMode=="EnrichmentTest")){
  if(verbose){print("Enrichment Test")}
  matCoefMargi=NULL
  pval_LRT=rep(NA,length(annotNames))
  AnalysisMargZero=glm(Border~1,data=binMat.mati,family=binomial())
  for(j in 1:length(annotNames)){
   formGLMMarg=as.formula(paste0("Border~",annotNames[j]))
   AnalysisMarg=glm(formGLMMarg,data=binMat.mati,family=binomial())
   AnalysisMargRes=summary(AnalysisMarg)
   coefj=AnalysisMargRes$coefficients[2,]
   if(LRT){
   #lrtj=lrtest(AnalysisMarg,AnalysisMargZero)
   Dj=(logLik(AnalysisMargZero)[1]-logLik(AnalysisMarg)[1])*-2
   pval_LRT[j]=1-pchisq(Dj,1)
   }
   matCoefMargi=rbind(matCoefMargi,coefj)
   if(verbose){print(annotNames[j])}
  }
  rownames(matCoefMargi)=annotNames
  freqBins=colSums(as.matrix(binMat.mati[,-1]))
  freqPeaks=sapply(genomicFeatureList.GR,length)
  matCoefMargi=cbind(annotNames,as.data.frame(matCoefMargi),pval_LRT,freqBins,freqPeaks)
  matCoefMargGLM=rbind(matCoefMargGLM,cbind(rep(bt,length(annotNames)),matCoefMargi))
 }

 # Multiple Logistic Regression
 if(sum(analysisMode=="MLR")){
  if(verbose){print("Multiple Logistic Regression")}
  
  binMat.mati$Border <- ifelse(binMat.mati$Border==1, "Yes", "No")
  binMat.mati$Border <- factor(binMat.mati$Border)
  
  inTrainingSet <- sample(dim(binMat.mati)[1],floor(dim(binMat.mati)[1]*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- binMat.mati[inTrainingSet,]
  test <- binMat.mati[-inTrainingSet,]
  
  saveRDS(test, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/test.rds")
  
  Analysisi=glm(Border~.,data=train,family=binomial())
  saveRDS(Analysisi, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/Analysisi.rds")
  
  prob <- predict(Analysisi,
					newdata=test,
					type="response")
  saveRDS(prob, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/prob.rds")
  
  pred <- prediction(prob, test$Border)
  saveRDS(pred, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/pred.rds")
  roc <- performance(pred,"tpr","fpr")
  saveRDS(roc, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mourad.roc.rds")
  
  
  list_MultGLM[[bt]]=summary(Analysisi)

  pval_LRT=rep(NA,length(annotNames))
  if(LRT){
   for(j in 1:length(annotNames)){
    formGLMMultj=as.formula(paste0("Border~",paste(annotNames[-j],collapse="+")))
    Analysisij=glm(formGLMMultj,data=binMat.mati,family=binomial())
    #lrtj=lrtest(Analysisij,Analysisi)
    Dj=(logLik(Analysisij)[1]-logLik(Analysisi)[1])*-2
    pval_LRT[j]=1-pchisq(Dj,1)
    if(verbose){print(paste0("LRT multi: ",annotNames[j]))}
   }
  }
  freqBins=colSums(binMat.mati[,-1])
  freqPeaks=sapply(genomicFeatureList.GR,length)
  matCoefMultGLMi=cbind(rownames(summary(Analysisi)$coefficients[-1,]),as.data.frame(summary(Analysisi)$coefficients[-1,]),pval_LRT,freqBins,freqPeaks)
  colnames(matCoefMultGLMi)[c(1,6)]=c("annotNames","pval_LRT")
  matCoefMultGLM=rbind(matCoefMultGLM,cbind(rep(bt,length(annotNames)),matCoefMultGLMi))
 }

 # Multiple Logistic Regression Estimated by Lasso
 if(sum(analysisMode=="MLRLasso")){
  if(verbose){print("Multiple Logistic Regression with Lasso Estimation")}
  
  inTrainingSet <- sample(dim(binMat.Mati)[1],floor(dim(binMat.Mati)[1]*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- binMat.Mati[inTrainingSet,]
  test <- binMat.Mati[-inTrainingSet,]
  
  saveRDS(test, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/test.lasso.rds")
  
  CVLasso=cv.glmnet(train[,-(1:2)],train[,1],family="binomial")
  lambda=CVLasso$lambda.min
  CVLassoError=CVLasso$cvm[which(CVLasso$lambda==lambda)]
  coefLasso=CVLasso$glmnet.fit$beta[,which(CVLasso$lambda==lambda)]
  matCoefMultLassoi=cbind(names(coefLasso),round(coefLasso,5))
  matCoefMultLasso=rbind(matCoefMultLasso,cbind(rep(bt,length(annotNames)),matCoefMultLassoi))
  
  ##########
  saveRDS(CVLasso, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/CVLasso.rds")
  prob <- predict(CVLasso,
					type="response",
					test[,-(1:2)], 
					s = lambda)
  saveRDS(prob, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/prob.lasso.rds")
  pred <- prediction(prob, test[,1])
  saveRDS(pred, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/pred.lasso.rds")
  perf <- performance(pred,"tpr","fpr")
  saveRDS(perf, "/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mourad.roc.lasso.rds")
 }
 
 # Multiple Logistic Regression with Interaction Terms 
 if(sum(analysisMode=="MLRInter")){
  if(verbose){print("Multiple Logistic Regression with Interaction Terms")}
  OneWayTerms=paste(annotNames,collapse='+')
  formInterTesti=as.formula(paste0("Border~",paste(OneWayTerms,collapse="+"),"+",paste(interactionTerms,collapse="+")))
  enrichInterTesti=glm(formInterTesti,data=binMat.mati,family=binomial())
  list_InterGLM[[bt]]=summary(enrichInterTesti)

  pval_LRTInter=rep(NA,length(annotNames)+length(interactionTerms))
  if(LRT){
   for(j in 1:length(interactionTerms)){
    if(length(interactionTerms)>1){
     formGLMInterj=as.formula(paste0("Border~",paste(OneWayTerms,collapse="+"),"+",paste(interactionTerms[-j],collapse="+")))
    }else{
     formGLMInterj=as.formula(paste0("Border~",paste(OneWayTerms,collapse="+")))
    }
    enrichInterTestij=glm(formGLMInterj,data=binMat.mati,family=binomial())
    #lrtInterj=lrtest(enrichInterTestij,enrichInterTesti)
    DInterj=(logLik(enrichInterTestij)[1]-logLik(enrichInterTesti)[1])*-2
    pval_LRTInter[length(annotNames)+j]=1-pchisq(DInterj,1)
    rm(enrichInterTestij)
    if(verbose){print(paste0("LRT multi: ",interactionTerms[j]))}
   }
  }
  matCoefInterGLMi=cbind(rownames(summary(enrichInterTesti)$coefficients[-1,]),as.data.frame(summary(enrichInterTesti)$coefficients[-1,]),pval_LRTInter)
  colnames(matCoefInterGLMi)[c(1,6)]=c("annotNames","pval_LRT")
  matCoefInterGLM=rbind(matCoefInterGLM,cbind(rep(bt,nrow(matCoefInterGLMi)),matCoefInterGLMi))
 }

 # Multiple Logistic Regression with Interaction Terms Estimated by Lasso
 if(sum(analysisMode=="MLRInterLasso")){
  if(verbose){print("Multiple Logistic Regression with Interaction Terms Estimated by Lasso")}
  OneWayTerms=paste(annotNames,collapse='+')
  formInterTesti=as.formula(paste0("Border~",paste(OneWayTerms,collapse="+"),"+",paste(interactionTerms,collapse="+")))
  binMatInter.Mati=sparse.model.matrix(formInterTesti,data=binMat.mati)

  CVLasso=cv.glmnet(binMatInter.Mati[,-1],binMat.Mati[,1],family="binomial")
  lambda=CVLasso$lambda.min
  CVLassoError=CVLasso$cvm[which(CVLasso$lambda==lambda)]
  coefInterLasso=CVLasso$glmnet.fit$beta[,which(CVLasso$lambda==lambda)]
  matCoefInterLassoi=cbind(names(coefInterLasso),round(coefInterLasso,5))
  matCoefInterLasso=rbind(matCoefInterLasso,cbind(rep(bt,nrow(matCoefInterLassoi)),matCoefInterLassoi))
 }


 if(verbose){print(paste0("Border type: ",bt," done!"))}
}

# Rename columns
if(sum(analysisMode=="EnrichmentTest")){
 colnames(matCoefMargGLM)[1:2]=c("BorderType","GenomicFeature")
}
if(sum(analysisMode=="MLR")){
 colnames(matCoefMultGLM)[1:2]=c("BorderType","GenomicFeature")
}
if(sum(analysisMode=="MLRLasso")){
 colnames(matCoefMultLasso)=c("BorderType","GenomicFeature","Estimate")
}
if(sum(analysisMode=="MLRInter")){
 colnames(matCoefInterGLM)[1:2]=c("BorderType","GenomicFeature")
}
if(sum(analysisMode=="MLRInterLasso")){
 colnames(matCoefInterLasso)[1:3]=c("BorderType","GenomicFeature","Estimate")
}



if(verbose){print("All analyses done")}


list2return=list(Enrich=matCoefMargGLM, MLR=matCoefMultGLM, MLRLasso=matCoefMultLasso,MLRInter=matCoefInterGLM,MLRInterLasso=matCoefInterLasso, Mat=binMat.Mat, MLRGLM=list_MultGLM, MLRInterGLM=list_InterGLM,CorGF=corGF)

return(list2return)
}


## Performing Model for regular MLR

GFDataType="bed"
analysisMode="MLR"
binSize=1000;borderSize=1000
LRT=FALSE
interactionTerms=""
testBorderType=FALSE
verbose=FALSE

BA_res=borderAnalysisFun(genomicFeatureList.GR=genomicFeatureList.GR,
                         GFDataType="bed",
                         annotNames=annotNames,
                         domains.GR=domains.GR,
                         seqInfoChr=seqInfoChr,
                         analysisMode="MLR",
                         binSize=1000,borderSize=1000,
                         LRT=FALSE,
                         interactionTerms="",
                         testBorderType=FALSE,
                         verbose=FALSE)

saveRDS(BA_res,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/BA_res.rds")

mourad.auc <- performance(pred,"auc")
mourad.auc <- mourad.auc@y.values[[1]]	
mourad.auc


pred.logit <- ifelse(prob>mean(prob),"Yes","No")
confMat <- confusionMatrix(pred.logit, test$Border, positive="Yes")

mouradperf <- matrix(nrow = 16, ncol=1)
rownames(mouradperf) <- c("TN",
 							"FN",
 							"FP",
 							"TP",
 							"Total",
							"Sensitivity",
 							"Specificity",
 							"Kappa",
 							"Accuracy",
 							"Precision",
 							"FPR",
 							"FNR",
 							"FOR",
 							"NPV",
 							"MCC",
 							"F1",
							"AUC")

  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  mouradperf[1,i] <- TN
  mouradperf[2,i] <- FN
  mouradperf[3,i] <- FP
  mouradperf[4,i] <- TP
  mouradperf[5,i] <- sum(confMat$table)
  mouradperf[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  mouradperf[7,i] <- as.vector(confMat$byClass["Specificity"])
  mouradperf[8,i] <- as.vector(confMat$overall["Kappa"])
  mouradperf[9,i] <- as.vector(confMat$overall["Accuracy"])
  mouradperf[10,i] <- TP/(TP+FP)
  mouradperf[11,i] <- FP/(FP+TN)
  mouradperf[12,i] <- FN/(FN+TN)
  mouradperf[13,i] <- FN/(FN+TN)
  mouradperf[14,i] <- TN/(TN+FN)
  mouradperf[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  mouradperf[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  mouradperf[17,i] <- mourad.auc



saveRDS(mouradperf,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mouradperf.rds")


## Performing Model with LASSO Regularization

GFDataType="bed"
analysisMode="MLRLasso"
binSize=1000;borderSize=1000
LRT=FALSE
interactionTerms=""
testBorderType=FALSE
verbose=FALSE

BA_res_lasso=borderAnalysisFun(genomicFeatureList.GR=genomicFeatureList.GR,
                         GFDataType="bed",
                         annotNames=annotNames,
                         domains.GR=domains.GR,
                         seqInfoChr=seqInfoChr,
                         analysisMode="MLRLasso",
                         binSize=1000,borderSize=1000,
                         LRT=FALSE,
                         interactionTerms="",
                         testBorderType=FALSE,
                         verbose=FALSE)
						 
						 
saveRDS(BA_res_lasso,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/BA_res_lasso.rds")

mourad.lasso.auc <- performance(pred,"auc")
mourad.lasso.auc <- mourad.lasso.auc@y.values[[1]]	
mourad.lasso.auc

head(mourad.roc.lasso@alpha.values[[1]][which(mourad.roc.lasso@y.values[[1]]>=0.7)])
#[1] 0.009212174 0.009212091 0.009212080 0.009212080 0.009212078 0.009212041
head(mourad.roc.lasso@alpha.values[[1]][which(mourad.roc.lasso@x.values[[1]]>=0.3)])
#[1] 0.009914878 0.009914728 0.009914718 0.009914581 0.009914533 0.009914220
pred.logit <- ifelse(prob.lasso>=mean(prob.lasso),"Yes","No")


test.lasso <- as.data.frame(as.matrix(test.lasso))
dim(test.lasso)


test.lasso <- test.lasso[,-2]


test.lasso$Border <- ifelse(test.lasso$Border==1, "Yes", "No")


confMat <- confusionMatrix(pred.logit, test.lasso$Border, positive="Yes")


mourad.lasso.perf <- matrix(nrow = 16, ncol=1)
rownames(mourad.lasso.perf) <- c("TN",
 							"FN",
 							"FP",
 							"TP",
 							"Total",
							"Sensitivity",
 							"Specificity",
 							"Kappa",
 							"Accuracy",
 							"Precision",
 							"FPR",
 							"FNR",
 							"FOR",
 							"NPV",
 							"MCC",
 							"F1",
							"AUC")

  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  mourad.lasso.perf[1,i] <- TN
  mourad.lasso.perf[2,i] <- FN
  mourad.lasso.perf[3,i] <- FP
  mourad.lasso.perf[4,i] <- TP
  mourad.lasso.perf[5,i] <- sum(confMat$table)
  mourad.lasso.perf[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  mourad.lasso.perf[7,i] <- as.vector(confMat$byClass["Specificity"])
  mourad.lasso.perf[8,i] <- as.vector(confMat$overall["Kappa"])
  mourad.lasso.perf[9,i] <- as.vector(confMat$overall["Accuracy"])
  mourad.lasso.perf[10,i] <- TP/(TP+FP)
  mourad.lasso.perf[11,i] <- FP/(FP+TN)
  mourad.lasso.perf[12,i] <- FN/(FN+TN)
  mourad.lasso.perf[13,i] <- FN/(FN+TN)
  mourad.lasso.perf[14,i] <- TN/(TN+FN)
  mourad.lasso.perf[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  mourad.lasso.perf[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  mourad.lasso.perf[17,i] <- mourad.lasso.auc


saveRDS(mourad.lasso.perf,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mourad.lasso.perf.rds")


## Comparing Final Models

m.auc <- performance(mourad.auc,"auc")
m.auc <- m.auc@ y.values[[1]]

m.l.auc <- performance(mourad.lasso.auc,"auc")
m.l.auc <- m.l.auc@ y.values[[1]]

#random forest
rf.auc <- round(mean(rfperf[17,]),3)

#Plotting AUCs
auc.plot <- data.frame(Model=c("Mourad MLR", 
                                       "Mourad MLR w/ LASSO",
                                       "Random Forest"),
                       auc=c(m.auc,
                             m.l.auc,
                             rf.auc))

auc.plot <- auc.plot[order(auc.plot$auc, decreasing=TRUE),]

auc.plot$Model <-factor(auc.plot$Model, 
                                     levels=auc.plot$Model)
									 
auc.plot

auc.final.p<-ggplot(data=auc.plot, aes(x=Model, y=auc, fill=Model)) + 
  xlab("Model") + ylab("AUC") +
  geom_bar(stat="identity") + ylim(0,1) +
  scale_fill_manual(values=grey(c(0,.5,.9)), guide=FALSE) +
  scale_x_discrete(labels= c("Random Forest",
                             "Mourad MLR \n w/ LASSO",
                             "Mourad MLR")) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Model Performance for Different \n Normalization Models")
auc.final.p
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/auc.final.p.png")


# ROC Curves

#mourad model
mourad.roc.fpr <- mourad.roc@ x.values[[1]]
mourad.roc.tpr <- mourad.roc@ y.values[[1]]
mourad.roc.df <- cbind.data.frame(fpr=mourad.roc.fpr, 
                                  tpr=mourad.roc.tpr,
                                  Model = rep("M", length(mourad.roc.tpr)))

#mourad model w/ lasso
mourad.lasso.roc.fpr <- mourad.roc@ x.values[[1]]
mourad.lasso.roc.tpr <- mourad.roc@ y.values[[1]]
mourad.lasso.roc.df <- cbind.data.frame(fpr=mourad.lasso.roc.fpr, 
                                  tpr=mourad.lasso.roc.tpr,
                                  Model = rep("MwL", length(mourad.lasso.roc.fpr)))

#random forest
rf.fpr <- rowMeans(rflst[[2]])
rf.tpr <- rowMeans(rflst[[1]])
rf.roc.df <- cbind.data.frame(fpr=rf.fpr, 
                                  tpr=rf.tpr,
                                  Model = rep("RF", length(rf.fpr)))

#concatenating data frames
allrocdat <- rbind.data.frame(mourad.roc.df, mourad.lasso.roc.df, rf.roc.df)

roc.final <- ggplot(data=allrocdat, aes(x=fpr, y=tpr, color=Model)) + 
  geom_line(size=1) +
  scale_colour_manual(name="Model",
    labels=c("Mourad MLR", 
             "Mourad MLR \n w/ LASSO",
             "Random Forest"),
    values=grey(c(.9,.5,0))) + 
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curves for Different \n Models")
roc.final
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.final.png")


#Variable Importance Plot

#RF
varimp.rf <- as.vector(rowMeans(rflst[[4]]))

rownames(rflst[[4]])[grep("Gm12878_", rownames(rflst[[4]]))] <- gsub("Gm12878_", "", rownames(rflst[[4]])[grep("Gm12878_", rownames(rflst[[4]]))])

#rownames(rflst[[4]])[grep("_dist", rownames(rflst[[4]]))] <- gsub("_dist", "", rownames(rflst[[4]])[grep("_dist", rownames(rflst[[4]]))])

varimp.rf.df <- data.frame(Feature=rownames(rflst[[4]]),
                           Importance=varimp.rf)
varimp.rf.df <- varimp.rf.df[order(varimp.rf.df$Importance),]
numvarrf <- dim(varimp.rf.df)[1]
varimp.rf.df <- varimp.rf.df[(numvarrf-19):numvarrf,]
varimp.rf.df$Feature <- factor(varimp.rf.df$Feature,levels=varimp.rf.df$Feature)

rfp <- ggplot(varimp.rf.df, aes(x=Feature, 
                                y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill="black") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Variable Importance Plot \n for Random Forest")
rfp
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/rfp.png")


# Estimates from Mourad Models

#mourad model
dim(mourad.summary)

sig.vars <- mourad.summary[mourad.summary$`Pr(>|z|)` < 0.05,]
dim(sig.vars)

sig.vars <- sig.vars[order(abs(sig.vars$Estimate), decreasing = TRUE),]
rownames(sig.vars) <- NULL
sig.vars <- sig.vars[1:20, which(colnames(sig.vars) %in% c("GenomicFeature","Estimate", "Pr(>|z|)"))]

#kable(sig.vars)
saveRDS(sig.vars,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/sig.vars.rds")


#mourad model w/ lasso
dim(mourad.lasso.summary)

mourad.lasso.summary <- as.data.frame(mourad.lasso.summary)
mourad.lasso.summary$Estimate <- as.numeric(as.character(mourad.lasso.summary$Estimate))

sig.vars.lasso <- mourad.lasso.summary[order(abs(mourad.lasso.summary$Estimate), decreasing = TRUE),]
rownames(sig.vars.lasso) <- NULL
sig.vars.lasso <- sig.vars.lasso[1:20, which(colnames(sig.vars.lasso) %in% c("GenomicFeature","Estimate"))]

#kable(sig.vars.lasso)
saveRDS(sig.vars.lasso,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/sig.vars.lasso.rds")


options(scipen = 999)

mourad.lasso.perf <- round(mourad.lasso.perf,2)
mourad.lasso.perf[1:5,1] <- round(mourad.lasso.perf[1:5,1],0)

#Our pipeline
rfperf <- round(as.matrix(rowMeans(rfperf)),2)
rfperf[1:5,1] <- round(rfperf[1:5,1],0)


perfdat.final <- cbind.data.frame(rownames(mouradperf), 
                            mouradperf,
                            mourad.lasso.perf,
                            rfperf)
rownames(perfdat.final) <- NULL
colnames(perfdat.final) <- c("Metric", "MLR", "MLR w/ LASSO Regularization", "Our Pipeline")

#kable(perfdat.final)
saveRDS(perfdat.final,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/perfdat.final.rds")



mccf1 <- data.frame(Metric = c(rep("MCC",3), rep("F1",3)),
                    Model = rep(c("MLR", 
                                      "MLR w/ LASSO", 
                                      "Our Pipeline"), 2),
                    Value = c(as.numeric(perfdat[15,2:4]), as.numeric(perfdat[16,2:4])))

mcc_f1_final <- ggplot(data=mccf1, aes(x=Model, y=Value, fill=Metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c('black','lightgray')) +
  xlab("Model") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mcc_f1_final 
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mcc_f1_final.png")


mcc_final <-ggplot(data=mccf1[1:3,], aes(x=Model, y=Value, fill=Model)) + 
  xlab("Model") + ylab("MCC") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=gray(rev(c(0,.4,.9))), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Measuring MCC \n Our Pipeline vs MLR (with and without LASSO)")
mcc_final 
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/mcc_final.png")


f1_final<-ggplot(data=mccf1[4:6,], aes(x=Model, y=Value, fill=Model)) + 
  xlab("Model") + ylab("F1") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=gray(rev(c(0,.4,.9))), guide=FALSE) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Measuring F1 \n Our Pipeline vs MLR (with and without LASSO)")
f1_final
ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/f1_final.png")



#finding common features between the models

#remove "_dist" from feature list of random forest
rffeat <- varimp.rf.df$Feature[order(varimp.rf.df$Importance, decreasing = TRUE)]
rffeat <- gsub("_dist", "", rffeat)
rffeat <- factor(rffeat)
rfrank <- 1:20

mrank <- match(rffeat, sig.vars$GenomicFeature)

mwlrank <- match(rffeat, sig.vars.lasso$GenomicFeature)

rankdf <- cbind.data.frame(Feature=rffeat, 
                           "Random Forest" = rfrank <- 1:20, 
                           "Mourad" = mrank,
                           "Mourad w/ LASSO" = mwlrank)

#kable(rankdf)
saveRDS(rankdf,"/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/rankdf.rds")


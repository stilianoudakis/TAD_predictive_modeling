# Measuring Performances and Providing Summaries

# Loading Libraries

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
library(limma)

gm12878_10kb <- readRDS("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/gm12878_10kb.rds")
gm12878_10kb_f <- readRDS("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/gm12878_10kb_f.rds")
gm12878_10kb_reduced <- readRDS("/home/stilianoudakisc/TAD_data_analysis/final_models/gm12878_10kb_reduced.rds")


distancetab <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/distancetab.rds")
featdistdf2 <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/featdistdf2.rds")

bounds <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/bounds.rds")
grlist <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/grlist.rds")

binslist10_center <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/binslist10_center.rds")

#############################################################################

#data exploration

domains <- read.table("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/arrowhead_data.txt", header=T)
domains <- domains[,1:3]

#removing chrX from boundary data
domains <- domains[-which(domains$Chromosome=="chrX"),]

#adding distance variable
domains$distance <- domains$End - domains$Start
#log2 transform of distance
domains$logdist <- log(domains$distance, base = 2)

#distance of TADs

tadlength <- ggplot() + geom_density(data=domains, aes(x=distance, group=Chromosome, color=Chromosome)) + 
	xlab("Length of TADs") + 
	xlim(0,1000000) + 
	ylab("Density")+ 
	scale_colour_manual(name="Chromosome",
		labels=as.character(
				unique(
					sort(
						as.numeric(
							substr(as.character(domains$Chromosome),4,nchar(as.character(domains$Chromosome))))))),
		values=rainbow(22)) +
	theme_minimal() 
tadlength
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/tadlength.png")

	
tadloglength <- ggplot() + geom_density(data=domains, aes(x=logdist, group=Chromosome, color=Chromosome)) + 
	xlab("Log Length of TADs") + 
	ylab("Density") + 
	scale_colour_manual(name="Chromosome",
		labels=as.character(
				unique(
					sort(
						as.numeric(
							substr(as.character(domains$Chromosome),4,nchar(as.character(domains$Chromosome))))))),
		values=rainbow(22)) +
	theme_minimal()
tadloglength
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/tadloglength.png")



x <- grlist

x <- x[-which(names(x)=="binslist10_center")]

for(i in 1:length(x)){
  
  if(is.element("chrX", as.character(seqnames(x[[i]])))){
    x[[i]] <- x[[i]][-which(as.character(seqnames(x[[i]]))=="chrX")]
  }
  
  if(is.element("chrY", as.character(seqnames(x[[i]])))){
    x[[i]] <- x[[i]][-which(as.character(seqnames(x[[i]]))=="chrY")]
  }
  
  mcols(x[[i]])$distance <- mcols(distanceToNearest(x[[i]], bounds))$distance
  mcols(x[[i]])$logdistance <- log(mcols(distanceToNearest(x[[i]], bounds))$distance+1,
                                          base = 2)
  
}


#Plots

# (Log) Distance from Bin with TAD in it to genomic feature
datalist = list()
tadbins <- binslist10_center[which(mcols(binslist10_center)$y==0)]
for(i in 1:length(x)){

data <- data.frame(Feature=rep(names(x)[i], length(tadbins)),
					Distance = mcols(distanceToNearest(tadbins, x[[i]]))$distance,
					LogDist = log(mcols(distanceToNearest(tadbins, x[[i]]))$distance+1, base=2))

datalist[[i]] <- data					
					
}

featdata <- dplyr::bind_rows(datalist)

dist2tad <- ggplot() + geom_density(data=featdata, aes(x=Distance, group=Feature, color=Feature)) + 
	xlab("Distances from TAD Boundary to Region") +
	xlim(0, 80000) + 
	#ylim(0,0.0004) +
	ylab("Density") + 
	theme_minimal() +
	theme(legend.position="none")
dist2tad
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/dist2tad.png")
	

	
logdist2tad <- ggplot() + geom_density(data=featdata, aes(x=LogDist, group=Feature, color=Feature)) + 
	xlab("Log Distance from TAD Boundary to Nearest Region") +
	#xlim(0, 1000000) + 
	ylab("Density") + 
	theme_minimal() +
	theme(legend.position="none")
logdist2tad 
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/logdist2tad.png")
	


nearTAD <- aggregate(featdata[, 3], list(featdata$Feature), mean)
colnames(nearTAD) <- c("Feature", "MeanLogDist")
nearTAD$Feature <- gsub("_gr_center","_dist",nearTAD$Feature)

nearTAD2 <- aggregate(featdata[, 3], list(featdata$Feature), sd)
colnames(nearTAD2) <- c("Feature", "SDLogDist")
nearTAD2$Feature <- gsub("_gr_center","_dist",nearTAD2$Feature)

nearTADdist <- cbind.data.frame(nearTAD, nearTAD2$SDLogDist)
colnames(nearTADdist)[3] <- "SDLogDist"


# (Log) Distance from Bin with no TAD in it to genomic feature
datalist = list()
set.seed(123)
notadbins <- binslist10_center[sample(which(mcols(binslist10_center)$y==0),15874)]
for(i in 1:length(x)){

data <- data.frame(Feature=rep(names(x)[i], length(notadbins)),
					Distance = mcols(distanceToNearest(notadbins, x[[i]]))$distance,
					LogDist = log(mcols(distanceToNearest(notadbins, x[[i]]))$distance+1, base=2))

datalist[[i]] <- data					
					
}

featdatant <- dplyr::bind_rows(datalist)

dist2notad <- ggplot() + geom_density(data=featdatant, aes(x=Distance, group=Feature, color=Feature)) + 
	xlab("Distances from Bin with No TAD Boundary to Region") +
	xlim(0, 80000) + 
	#ylim(0,0.0004) +
	ylab("Density") + 
	theme_minimal() +
	theme(legend.position="none")
dist2notad
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/dist2notad.png")
	

	
logdist2notad <- ggplot() + geom_density(data=featdatant, aes(x=LogDist, group=Feature, color=Feature)) + 
	xlab("Log Distance from Bin with No TAD to Nearest Region") +
	#xlim(0, 1000000) + 
	ylab("Density") + 
	theme_minimal() +
	theme(legend.position="none")
logdist2notad 
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/logdist2notad.png")
	

noTAD <- aggregate(featdatant[, 3], list(featdatant$Feature), mean)
colnames(noTAD) <- c("Feature", "MeanLogDist")
noTAD$Feature <- gsub("_gr_center","_dist",noTAD$Feature)

noTAD2 <- aggregate(featdatant[, 3], list(featdatant$Feature), sd)
colnames(noTAD2) <- c("Feature", "SDLogDist")
noTAD2$Feature <- gsub("_gr_center","_dist",noTAD2$Feature)

noTADdist <- cbind.data.frame(noTAD, noTAD2$SDLogDist)
colnames(noTADdist)[3] <- "SDLogDist"



# (Log) Distance from every Bin to genomic feature
datalist = list()
for(i in 1:length(x)){

data <- data.frame(Feature=rep(names(x)[i], length(binslist10_center)),
					Distance = mcols(distanceToNearest(binslist10_center, x[[i]]))$distance,
					LogDist = log(mcols(distanceToNearest(binslist10_center, x[[i]]))$distance+1, base=2))

datalist[[i]] <- data					
					
}

featdataallbins <- dplyr::bind_rows(datalist)

dist2alltad <- ggplot() + geom_density(data=featdataallbins, aes(x=Distance, group=Feature, color=Feature)) + 
	xlab("Distances from Bin to Region") +
	xlim(0, 40000) + 
	#ylim(0,0.0004) +
	ylab("Density") + 
	theme_minimal() +
	theme(legend.position="none")
dist2alltad
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/dist2alltad.png")
	

	
logdist2alltad <- ggplot() + geom_density(data=featdataallbins, aes(x=LogDist, group=Feature, color=Feature)) + 
	xlab("Log Distance from Bin to Nearest Region") +
	#xlim(0, 1000000) + 
	ylab("Density") + 
	theme_minimal() +
	theme(legend.position="none")
logdist2alltad 
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/logdist2alltad.png")
	

#propdata <- data.frame(Chromosome = sort(rep(1:22,2)), TAD = rep(c(0,1),22), Percent = c(20995,  1529,
#																						22525 , 1305,
#																						18402  ,1080,
#																						17978   ,793,
#																						16846,   904,
#																						15762 ,  984,
#																						14730  , 818,
#																						13553,   737,
#																						11346 ,  655,
#																						12334  , 798,
#																						12253   ,864,
#																						12170,   885,
#																						9113  ,451,
#																						8311,  519,
#																						7602 , 562,
#																						7412  ,476,
#																						7085 , 700,
#																						7071  ,399,
#																						5073  ,511,
#																						5518  ,439,
#																						3350  ,167,
#																						3195  ,298))
#																						
#propdata$TAD <- factor(propdata$TAD)
#
#tadprops <- ggplot(data=propdata, aes(x=Chromosome, y=Percent, fill=TAD)) + 
#	geom_bar(position = "fill", stat="identity") + 
#	theme_minimal() +
#	theme(legend.position="none")
#tadprops
#ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/tadprops.png")


#############################################################################

# Listing Features Before and After Filtering 

dim(gm12878_10kb)
dim(gm12878_10kb_f)

total_feats <- names(gm12878_10kb)[-1]
saveRDS(total_feats,"/home/stilianoudakisc/TAD_data_analysis/measuring_performances/total_feats.rds")

filtered_feats <- setdiff(names(gm12878_10kb)[-1], names(gm12878_10kb_f)[-1])
saveRDS(filtered_feats,"/home/stilianoudakisc/TAD_data_analysis/measuring_performances/filtered_feats.rds")

remaining_feats <- names(gm12878_10kb_f)[-1]
saveRDS(remaining_feats,"/home/stilianoudakisc/TAD_data_analysis/measuring_performances/remaining_feats.rds")

#############################################################################

# Normalization Performance

setwd("/home/stilianoudakisc/TAD_data_analysis/comparing_normalization/")

enetlst_ls <- readRDS("enetlst_ls.rds")
enetperf_ls <- readRDS("enetperf_ls.rds")
enetlst_lns <- readRDS("enetlst_lns.rds")
enetperf_lns <- readRDS("enetperf_lns.rds")
enetlst_nls <- readRDS("enetlst_nls.rds")
enetperf_nls <- readRDS("enetperf_nls.rds")
enetlst_nlns <- readRDS("enetlst_nlns.rds")
enetperf_nlns <- readRDS("enetperf_nlns.rds")

options(scipen = 999)

lstab <- round(enetperf_ls,3)
lnstab <- round(enetperf_lns,3)
nlstab <- round(enetperf_nls,3)
nlnstab <- round(enetperf_nlns,3)

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

saveRDS(perfdat.norm,"/home/stilianoudakisc/TAD_data_analysis/measuring_performances/perfdat.norm.rds")

#datatable(perfdat.norm)


mccf1 <- data.frame(Metric = c(rep("MCC",4), rep("F1",4)),
                    Technique = rep(c("Log/Std", 
                                      "Log/Un-Std", 
                                      "No Log/Std", 
                                      "No Log/Un-Std"), 2),
                    Value = c(as.numeric(perfdat.norm[15,2:5]), as.numeric(perfdat.norm[16,2:5])))

mcc_f1_norm <- ggplot(data=mccf1, aes(x=Technique, y=Value, fill=Metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #scale_fill_manual(values=c('black','lightgray')) +
  xlab("Normalization Technique") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mcc_f1_norm
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/mcc_f1_norm.png")


# Log tranformed and standardized

varimp.ls <- as.vector(enetlst_ls[[3]])
Labels <- colnames(gm12878_10kb_f)[-1]
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
           fill=rainbow(4)[1]) +
  coord_flip() +
  theme_minimal() +
  ggtitle("Log Transformed \n & Standardized")
p.ls
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/p.ls.png")

# Log tranformed and un-standardized

varimp.lns <- as.vector(enetlst_lns[[3]])
Labels <- colnames(gm12878_10kb_f)[-1]
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
           fill=rainbow(4)[2]) +
  coord_flip() +
  theme_minimal() +
  ggtitle("Log Transformed \n & Un-Standardized")
p.lns
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/p.lns.png")

# Not Log tranformed and Standardized

varimp.nls <- as.vector(enetlst_nls[[3]])
Labels <- colnames(gm12878_10kb_f)[-1]
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
           fill=rainbow(4)[3]) +
  coord_flip() +
  theme_minimal() +
  ggtitle("Not Log Transformed \n & Standardized")
p.nls
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/p.nls.png")

# Not Log tranformed and Un-Standardized

varimp.nlns <- as.vector(enetlst_nlns[[3]])
Labels <- colnames(gm12878_10kb_f)[-1]
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
           fill=rainbow(4)[4]) +
  coord_flip() +
  theme_minimal() +
  ggtitle("Not Log Transformed \n & Standardized")
p.nlns
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/p.nlns.png")

all.norm.p <- grid.arrange(p.ls,p.lns,p.nls,p.nlns,ncol=2)
all.norm.p
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/all.norm.p.png")


#############################################################################

# Class Imbalance Performance

setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/")

enetlst_sm <- readRDS("enetlst_sm_lns.rds")
enetlst_bs<-readRDS("enetlst_bs.rds")
perfdat_sm_bs <- readRDS("perfdat_sm_bs.rds")


onetwo <- data.frame(fpr=enetlst_sm[[2]][,1],tpr=enetlst_sm[[1]][,1], Combo = "100/200");
twotwo <- data.frame(fpr=enetlst_sm[[2]][,2],tpr=enetlst_sm[[1]][,2], Combo = "200/200");
threetwo <- data.frame(fpr=enetlst_sm[[2]][,3],tpr=enetlst_sm[[1]][,3], Combo = "300/200");
fourtwo <- data.frame(fpr=enetlst_sm[[2]][,4],tpr=enetlst_sm[[1]][,4], Combo = "400/200");
bootsamps <- data.frame(fpr=rowMeans(enetlst_bs[[2]]),tpr=rowMeans(enetlst_bs[[1]]), Combo = "Under-Sampling");

allrocdat <- rbind.data.frame(onetwo,
                              twotwo,
                              threetwo,
                              fourtwo,
							  bootsamps)

roc.class_unbalance <- ggplot(data=allrocdat, aes(x=fpr, y=tpr, color=Combo)) + 
  geom_line(size=1) +
  scale_colour_manual(name="Sampling \n Technique",
    labels=c("100/200", 
             "200/200",
             "300/200",
             "400/200",
             "Under-Sampling"),
    values=rainbow(5)) + 
  xlab("1-Specificity") + 
  ylab("Sensitivity") + 
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept=0, slope=1) +
  theme_minimal() +
  ggtitle("ROC Curves for Different \n Normalization Techniques")
  
  
mccf1 <- data.frame(Metric = c(rep("MCC",5), rep("F1",5)),
                    Technique = rep(c("100/200", 
                    "200/200", 
                    "300/200", 
                    "400/200", 
                    "Under-Sampling"), 2),
                    Value = c(as.numeric(perfdat_sm_bs[15,2:6]), 
						as.numeric(perfdat_sm_bs[16,2:6])))

mcc_f1_class_imb <- ggplot(data=mccf1, aes(x=Technique, y=Value, fill=Metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #scale_fill_manual(values=c("black" ,"gray")) +
  xlab("Balancing Technique") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mcc_f1_class_imb
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/mcc_f1_class_imb.png")



#############################################################################

# Measuring Variable Selection Performance

setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/")


gm12878_10kb_f<-readRDS("gm12878_10kb_f.rds")
gm12878_10kb_fwd<-readRDS("gm12878_10kb_fwd.rds")
gm12878_10kb_bwd<-readRDS("gm12878_10kb_bwd.rds")
gm12878_10kb_both<-readRDS("gm12878_10kb_both.rds")

vars.fwd <- names(gm12878_10kb_fwd)[-1]
vars.bwd <- names(gm12878_10kb_bwd)[-1]
vars.both <- names(gm12878_10kb_both)[-1]

vars.fwd
vars.bwd
vars.both

enetperf_fwd <-readRDS("enetperf_fwd.rds")
enetperf_bwd <-readRDS("enetperf_bwd.rds")
enetperf_both <-readRDS("enetperf_both.rds")

enetlst_fwd <-readRDS("enetlst_fwd.rds")
enetlst_bwd <-readRDS("enetlst_bwd.rds")
enetlst_both <-readRDS("enetlst_both.rds")


options(scipen = 999)

fwdtab <- round(as.matrix(rowMeans(enetperf_fwd)),3)
bwdtab <- round(as.matrix(rowMeans(enetperf_bwd)),3)
bothtab <- round(as.matrix(rowMeans(enetperf_both)),3)

fwdtab[1:5,1] <- round(fwdtab[1:5,1],0)
bwdtab[1:5,1] <- round(bwdtab[1:5,1],0)
bothtab[1:5,1] <- round(bothtab[1:5,1],0)

perfdat.var.sel <- cbind.data.frame(rownames(enetperf_fwd), 
                            fwdtab,
                            bwdtab,
                            bothtab)
							
rownames(perfdat.var.sel) <- NULL
colnames(perfdat.var.sel) <- c("Metric", "Forward", "Backward", "Both")
perfdat.var.sel$Metric <- as.character(perfdat.var.sel$Metric)

nvc <- c("NumVarsChosen",length(vars.fwd),length(vars.bwd),length(vars.both))

perfdat.var.sel <- rbind.data.frame(perfdat.var.sel, nvc)
perfdat.var.sel


mccf1 <- data.frame(Metric = c(rep("MCC",3), rep("F1",3)),
                    Technique = rep(c("Forward", 
                                      "Backward", 
                                      "Both"), 2),
                    Value = c(as.numeric(perfdat.var.sel[15,2:4]), as.numeric(perfdat.var.sel[16,2:4])))

mcc_f1_var_sel <- ggplot(data=mccf1, aes(x=Technique, y=Value, fill=Metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylim(0,1) +
  #scale_fill_manual(values=c('black','lightgray')) +
  xlab("Variable Selection \n Technique") +
  annotate("text", x=1, y=.85, label= "26", size=6) + 
  annotate("text", x=2, y=.85, label= "25", size=6) + 
  annotate("text", x=3, y=.85, label= "27", size=6) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
mcc_f1_var_sel
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/mcc_f1_var_sel.png")


intersect(vars.fwd,intersect(vars.bwd,vars.both))
#intersect(vars.fwd,intersect(vars.bwd,predictors(rfeModel)))
#intersect(intersect(vars.fwd,intersect(vars.bwd,vars.both)),
#          intersect(vars.fwd,intersect(vars.bwd,predictors(rfeModel))))

fwd <- (names(gm12878_10kb_f) %in% vars.fwd)
bwd <- (names(gm12878_10kb_f) %in% vars.bwd)
both <- (names(gm12878_10kb_f) %in% vars.both)
#rfe <- (names(chr1_gm12878_f) %in% predictors(rfeModel))

#fwd compared to bwd
#venndatfb <- cbind(fwd,bwd)
#fb <- vennCounts(venndatfb)
#vennDiagram(fb, include = "both", 
#  names = c("Forward", "Backward"), 
#  cex = 1, counts.col = "red")



venndat1 <- cbind(fwd,bwd,both)
#venndat2 <- cbind(fwd,bwd,rfe)
#venndat3 <- cbind(fwd,bwd,both,rfe)

a <- vennCounts(venndat1)
#b <- vennCounts(venndat2)
#c <- vennCounts(venndat3)

vennd <- vennDiagram(a, include = "both", 
  names = c("Forward", "Backward", "Both"), 
  cex = 1, counts.col = "red")
  
png(filename="/home/stilianoudakisc/TAD_data_analysis/measuring_performances/vennd.png")
vennd
dev.off()

#################################################################################

# Final Models Performance

setwd("/home/stilianoudakisc/TAD_data_analysis/mourad_model/")

#Mourad

mouradperf <- readRDS("mouradperf.rds")
BA_res <- readRDS("BA_res.rds")

mourad.lasso.perf <- readRDS("mourad.lasso.perf.rds")
BA_res_lasso <- readRDS("BA_res_lasso.rds")

#Our Pipeline

setwd("/home/stilianoudakisc/TAD_data_analysis/final_models/")

#enet unbalanced

enetlst_unb <- readRDS("enetlst_unb.rds")
enetperf_unb <- readRDS("enetperf_unb.rds")

#enet

enetlst <- readRDS("enetlst.rds")
enetperf <- readRDS("enetperf.rds")


#rf unbalanced

#rflst_unb <- readRDS("rflst_unb.rds")
#rfperf_unb <- readRDS("rfperf_unb.rds")

#rf 

#rflst <- readRDS("rflst.rds")
#rfperf <- readRDS("rfperf.rds")

#full data

setwd("/home/stilianoudakisc/TAD_data_analysis/full_data_model/")

enetlst_full <- readRDS("enetlst_full.rds")
enetperf_full <- readRDS("enetperf_full.rds")

#rflst_full <- readRDS("rflst_full.rds")
#rfperf_full <- readRDS("rfperf_full.rds")



# comparing models

perfdat.final <- data.frame(Metric = rownames(mouradperf),
							#Full_Enet = enetperf_full[,1],
							#Full_RF = rfperf_full[,1],
							MLR = mouradperf[,1],
							MLR_Lasso = mourad.lasso.perf[,1],
							#Enet_unb = enetperf_unb[,1],
							Enet = rowMeans(enetperf) #,
							#RF_unb = rfperf_unb[,1],
							#RF = rowMeans(rfperf)
							)
							
#perfdat.final <- data.frame(Metric = rownames(mouradperf),
#							MLR = mouradperf[,1],
#							MLR_Lasso = mourad.lasso.perf[,1],
#							Enet_unb = enetperf_unb[,1],
#							Enet = rowMeans(enetperf),
#							RF = rowMeans(rfperf))
							
perfdat.final[,-1] <- round(perfdat.final[,-1], 3)
			
							
mccf1 <- data.frame(Metric = c(rep("MCC",3), rep("F1",3)),
                    Technique = rep(c("MLR",
										"MLR w/ LASSO",
										"E-Net"), 2),
                    Performance = c(as.numeric(perfdat.final[15,2:4]), as.numeric(perfdat.final[16,2:4])))

mcc_f1_final <- ggplot(data=mccf1, aes(x=Technique, y=Value, fill=Metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #scale_fill_manual(values=c('black','lightgray')) +
  xlab("Model") +
  theme_minimal() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
		axis.text.x = element_text(angle = 90, hjust = 1),
		legend.title=element_text(size=20), 
		legend.text=element_text(size=15),
		legend.position = c(.7,.9))
mcc_f1_final
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/mcc_f1_final.png")						
							
# variable importance

#mlr

mlr <- data.frame(Feature = BA_res$MLR$GenomicFeature,
				Estimate = BA_res$MLR$Estimate,
				sigma = apply(BA_res$Mat[,-c(1,2)], 2, sd))
mlr$standpred <- mlr$Estimate*mlr$sigma	
mlr$absstandpred <- abs(mlr$standpred)
mlr <- mlr[order(mlr$absstandpred, decreasing=TRUE),]
mlr <- mlr[-which(mlr$Feature=="A" | mlr$Feature=="B"),]

varimp.mlr <- mlr$absstandpred
Labels <- mlr$Feature
#Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.mlr.df <- data.frame(Feature=Labels,
                                 Importance=varimp.mlr)
varimp.mlr.df <- varimp.mlr.df[order(varimp.mlr.df$Importance),]
varimp.mlr.df <- varimp.mlr.df[(dim(varimp.mlr.df)[1]-9):dim(varimp.mlr.df)[1],]
varimp.mlr.df$Feature <- factor(varimp.mlr.df$Feature,
                                     levels=varimp.mlr.df$Feature)
p.mlr <- ggplot(varimp.mlr.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill=rainbow(3)[1]) +
  coord_flip() +
  theme_minimal() +
  ggtitle("MLR")
p.mlr
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/p.mlr.png")


#mlr w/ lasso

mlrlasso <- data.frame(Feature = BA_res_lasso$MLRLasso[,2],
						Estimate = as.numeric(BA_res_lasso$MLRLasso[,3]),
						sigma = apply(BA_res_lasso$Mat[,-c(1,2)], 2, sd))
mlrlasso$standpred <- mlrlasso$Estimate*mlrlasso$sigma
mlrlasso$absstandpred <- abs(mlrlasso$standpred)
mlrlasso <- mlrlasso[order(mlrlasso$absstandpred, decreasing=TRUE),]
mlrlasso <- mlrlasso[-which(mlrlasso$Feature=="A" | mlrlasso$Feature=="B"),]
			
varimp.mlrlasso <- mlrlasso$absstandpred
Labels <- mlrlasso$Feature
#Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.mlrlasso.df <- data.frame(Feature=Labels,
                                 Importance=varimp.mlrlasso)
varimp.mlrlasso.df <- varimp.mlrlasso.df[order(varimp.mlrlasso.df$Importance),]
varimp.mlrlasso.df <- varimp.mlrlasso.df[(dim(varimp.mlrlasso.df)[1]-9):dim(varimp.mlrlasso.df)[1],]
varimp.mlrlasso.df$Feature <- factor(varimp.mlrlasso.df$Feature,
                                     levels=varimp.mlrlasso.df$Feature)
p.mlrlasso <- ggplot(varimp.mlrlasso.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill=rainbow(3)[2]) +
  coord_flip() +
  theme_minimal() +
  ggtitle("MLR w/ LASSO")
p.mlrlasso
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/p.mlrlasso.png")

			
# enet

enet <- data.frame(Feature = rownames(enetlst[[4]]),
					Estimate = rowMeans(enetlst[[4]]),
					sigma = apply(gm12878_10kb_reduced[,-1], 2, sd))
enet$standpred <- enet$Estimate*enet$sigma
enet$absstandpred <- abs(enet$standpred)
enet <- enet[order(enet$absstandpred, decreasing=TRUE),]
enet <- enet[-which(enet$Feature=="A" | enet$Feature=="B"),]
rownames(enet) <- NULL

varimp.enet <- enet$absstandpred
Labels <- as.character(enet$Feature)
Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
varimp.enet.df <- data.frame(Feature=Labels,
                                 Importance=varimp.enet)
varimp.enet.df <- varimp.enet.df[order(varimp.enet.df$Importance),]
varimp.enet.df <- varimp.enet.df[(dim(varimp.enet.df)[1]-9):dim(varimp.enet.df)[1],]
varimp.enet.df$Feature <- factor(varimp.enet.df$Feature,
                                     levels=varimp.enet.df$Feature)
p.enet <- ggplot(varimp.enet.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Importance Plot for Gradient Boosting Machine") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill=rainbow(3)[3]) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
		#axis.text.x = element_text(angle = 90, hjust = 1),
		#legend.title=element_text(size=20), 
		#legend.text=element_text(size=15),
		#legend.position = c(.7,.9)
  ggtitle("Elastic-Net")
p.enet	
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/p.enet.png")


varimp.final <- grid.arrange(p.mlr, p.mlrlasso,p.enet, ncol=1) 
varimp.final
ggsave("/home/stilianoudakisc/TAD_data_analysis/measuring_performances/varimp.final.png")


#featdata2 <- aggregate(featdata[, 3], list(featdata$Feature), mean)
#colnames(featdata2) <- c("Feature", "MeanLogDist")
#featdata2$Feature <- gsub("_gr_center","_dist",featdata2$Feature)

#commoncoefs <- na.omit(match(enet$Feature, featdata2$Feature))
#commoncoefs2 <- na.omit(match(featdata2$Feature, enet$Feature))

#featdata3 <- featdata2[commoncoefs,]
#enet2 <- enet[commoncoefs2,]

#enetcoefs <- cbind.data.frame(featdata3, Coefficient=abs(enet2$Estimate))

#plot(enetcoefs$MeanDist, enetcoefs$Coefficient)

#featdata2 <- featdata2[order(featdata2$MeanDist, decreasing=FALSE),]




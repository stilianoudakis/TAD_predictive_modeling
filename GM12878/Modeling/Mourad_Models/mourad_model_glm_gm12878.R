setwd("/home/stilianoudakisc/TAD_data_analysis/mourad_model/")

library(caret)
library(GenomicRanges)
library(Matrix)
library(glmnet)
library(pROC)
library(ROCR)
#library(zoo)
#library(pracma)

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
  
  saveRDS(test, "test.rds")
  
  Analysisi=glm(Border~.,data=train,family=binomial())
  saveRDS(Analysisi, "Analysisi.rds")
  
  prob <- predict(Analysisi,
					newdata=test,
					type="response")
  saveRDS(prob, "prob.rds")
  
  pred <- prediction(prob, test$Border)
  saveRDS(pred, "pred.rds")
  roc <- performance(pred,"tpr","fpr")
  saveRDS(roc, "mourad.roc.rds")
  
  
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
  CVLasso=cv.glmnet(binMat.Mati[,-(1:2)],binMat.Mati[,1],family="binomial")
  lambda=CVLasso$lambda.min
  CVLassoError=CVLasso$cvm[which(CVLasso$lambda==lambda)]
  coefLasso=CVLasso$glmnet.fit$beta[,which(CVLasso$lambda==lambda)]
  matCoefMultLassoi=cbind(names(coefLasso),round(coefLasso,5))
  matCoefMultLasso=rbind(matCoefMultLasso,cbind(rep(bt,length(annotNames)),matCoefMultLassoi))
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


seqInfoChr <- readRDS("seqInfoChr.rds")
annotNames <- readRDS("annotNames.rds")
genomicFeatureList.GR <- readRDS("genomicFeatureList.GR.rds")
domains.GR <- readRDS("domains.GR.rds")
GFDataType="bed"
analysisMode="MLR"
binSize=10000;borderSize=1
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
                         binSize=10000,borderSize=1,
                         LRT=FALSE,
                         interactionTerms="",
                         testBorderType=FALSE,
                         verbose=FALSE)

saveRDS(BA_res, "BA_res.rds")	



pred <- readRDS("pred.rds")
mourad.auc <- performance(pred,"auc")
mourad.auc <- mourad.auc@y.values[[1]]	
mourad.auc

prob <- readRDS("prob.rds")
test <- readRDS("test.rds") 
	 
pred.logit <- ifelse(prob>mean(prob),"Yes","No")
confMat <- confusionMatrix(pred.logit, test$Border, positive="Yes")

mouradperf <- matrix(nrow = 17, ncol=1)
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
  mouradperf[1,1] <- confMat$table[1,1]
  mouradperf[2,1] <- confMat$table[1,2]
  mouradperf[3,1] <- confMat$table[2,1]
  mouradperf[4,1] <- confMat$table[2,2]
  mouradperf[5,1] <- sum(confMat$table)
  mouradperf[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  mouradperf[7,1] <- as.vector(confMat$byClass["Specificity"])
  mouradperf[8,1] <- as.vector(confMat$overall["Kappa"])
  mouradperf[9,1] <- as.vector(confMat$overall["Accuracy"])
  mouradperf[10,1] <- TP/(TP+FP)
  mouradperf[11,1] <- FP/(FP+TN)
  mouradperf[12,1] <- FN/(FN+TN)
  mouradperf[13,1] <- FN/(FN+TN)
  mouradperf[14,1] <- TN/(TN+FN)
  mouradperf[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  mouradperf[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  mouradperf[17,1] <- mourad.auc

saveRDS(mouradperf, "mouradperf.rds")
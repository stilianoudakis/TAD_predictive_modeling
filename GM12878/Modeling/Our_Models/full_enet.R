#Elastic Net on Full Data

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


setwd("/home/stilianoudakisc/TAD_data_analysis/full_data_model/")

gm12878_10kb <- readRDS("gm12878_10kb.rds")

gm12878_10kb$y <- factor(gm12878_10kb$y)

levels(gm12878_10kb$y) <- c("No","Yes")

#set tuning parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 3,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


#splitting data
inTrainingSet <- sample(length(gm12878_10kb$y),floor(length(gm12878_10kb$y)*.7))
#inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
train <- gm12878_10kb[inTrainingSet,]
test <- gm12878_10kb[-inTrainingSet,]
  
  
#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_full <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   varimp <- matrix(nrow=dim(gm12878_10kb)[2]-1,
                                    ncol=1))
rownames(enetlst_full[[3]]) <- colnames(gm12878_10kb)[-1]


enetperf_full <- matrix(nrow = 17, ncol=1)
rownames(enetperf_full) <- c("TN",
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
  enetlst_full[[1]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_full[[2]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_full[[3]] <- varImp(eNetModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_full[1,1] <- TN
  enetperf_full[2,1] <- FN
  enetperf_full[3,1] <- FP
  enetperf_full[4,1] <- TP
  enetperf_full[5,1] <- sum(confMat$table)
  enetperf_full[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_full[7,1] <- as.vector(confMat$byClass["Specificity"])
  enetperf_full[8,1] <- as.vector(confMat$overall["Kappa"])
  enetperf_full[9,1] <- as.vector(confMat$overall["Accuracy"])
  enetperf_full[10,1] <- TP/(TP+FP)
  enetperf_full[11,1] <- FP/(FP+TN)
  enetperf_full[12,1] <- FN/(FN+TN)
  enetperf_full[13,1] <- FN/(FN+TN)
  enetperf_full[14,1] <- TN/(TN+FN)
  enetperf_full[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_full[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_full[17,1] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  



saveRDS(enetlst_full, "enetlst_full.rds")
saveRDS(enetperf_full, "enetperf_full.rds")
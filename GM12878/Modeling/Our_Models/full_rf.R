#RF on Full Data

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


#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst_full <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   varimp <- matrix(nrow=dim(gm12878_10kb)[2]-1,
                                    ncol=1))
rownames(rflst_full[[3]]) <- colnames(gm12878_10kb)[-1]


rfperf_full <- matrix(nrow = 17, ncol=1)
rownames(rfperf_full) <- c("TN",
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

data <- gm12878_10kb
  
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
  rflst_full[[1]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,1]
  rflst_full[[2]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,2]
  rflst_full[[3]] <- varImp(rfModel)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.rfModel2 <- predict(rfModel,
                                     newdata=test,
                                     type="raw")
  confMat <- confusionMatrix(data=pred.rfModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  rfperf_full[1,1] <- TN
  rfperf_full[2,1] <- FN
  rfperf_full[3,1] <- FP
  rfperf_full[4,1] <- TP
  rfperf_full[5,1] <- sum(confMat$table)
  rfperf_full[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  rfperf_full[7,1] <- as.vector(confMat$byClass["Specificity"])
  rfperf_full[8,1] <- as.vector(confMat$overall["Kappa"])
  rfperf_full[9,1] <- as.vector(confMat$overall["Accuracy"])
  rfperf_full[10,1] <- TP/(TP+FP)
  rfperf_full[11,1] <- FP/(FP+TN)
  rfperf_full[12,1] <- FN/(FN+TN)
  rfperf_full[13,1] <- FN/(FN+TN)
  rfperf_full[14,1] <- TN/(TN+FN)
  rfperf_full[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  rfperf_full[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  rfperf_full[17,1] <- pROC::auc(pROC::roc(test$y, pred.rfModel))




saveRDS(rflst_full,"/home/stilianoudakisc/TAD_data_analysis/full_data_model/rflst_full.rds")
saveRDS(rfperf_full,"/home/stilianoudakisc/TAD_data_analysis/full_data_model/rfperf_full.rds")


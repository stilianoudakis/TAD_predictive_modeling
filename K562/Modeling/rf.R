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



# Reading in data

setwd("/home/stilianoudakisc/TAD_data_analysis/K562/final_models/")



k562_10kb_reduced <- readRDS("k562_10kb_reduced.rds")



# Final Models:

## Random Forest

### set number of bootstrap samples


bootsamps = 5

### set tuning parameters


fitControl <- trainControl(method = "cv",
                           number = 10,
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
						
###create a matrix of row ids that represent the zero class

##the number of rows will match the one class
##the number of columns match the number of bootstrap samples

sampids <- matrix(ncol=bootsamps, 
                  nrow=length(k562_10kb_reduced$y[which(k562_10kb_reduced$y=="Yes")]))


###filling in the sample ids matrix

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(k562_10kb_reduced$y=="No"),
                        length(which(k562_10kb_reduced$y=="Yes")),
                        replace = TRUE)
}


###Performing model

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst <- list(tpr <- matrix(nrow=ceiling((length(which(k562_10kb_reduced$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                fpr <- matrix(nrow=ceiling((length(which(k562_10kb_reduced$y=="Yes"))*2)*.3), 
                              ncol=bootsamps),
                varimp <- matrix(nrow=dim(k562_10kb_reduced)[2]-1,
                                 ncol=bootsamps))
rownames(rflst[[3]]) <- colnames(k562_10kb_reduced)[-1]

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
  data <- rbind.data.frame(k562_10kb_reduced[which(k562_10kb_reduced$y=="Yes"),],
                           k562_10kb_reduced[sampids[,i],])
  
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
  rfperf[17,i] <- pROC::auc(pROC::roc(test$y, pred.rfModel))

}


saveRDS(rflst,"/home/stilianoudakisc/TAD_data_analysis/final_models/rflst.rds")
saveRDS(rfperf,"/home/stilianoudakisc/TAD_data_analysis/final_models/rfperf.rds")




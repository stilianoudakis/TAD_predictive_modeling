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

setwd("/home/stilianoudakisc/TAD_data_analysis/final_models/")



gm12878_10kb_reduced <- readRDS("gm12878_10kb_reduced.rds")



# Final Models:

## Elastic Net

### set tuning parameters


fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
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
						

###Performing model

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_reduced$y=="Yes"))*2)*.3), 
                              ncol=1),
                fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_reduced$y=="Yes"))*2)*.3), 
                              ncol=1),
                varimp <- matrix(nrow=dim(gm12878_10kb_reduced)[2]-1,
                                 ncol=1),
				coefs <- matrix(nrow=dim(gm12878_10kb_reduced)[2]-1,
									ncol=1))
rownames(enetlst[[3]]) <- colnames(gm12878_10kb_reduced)[-1]
rownames(enetlst[[4]]) <- colnames(gm12878_10kb_reduced)[-1]				 
								 
enetperf <- matrix(nrow = 17, ncol=1)
rownames(enetperf) <- c("TN",
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

#Performing Elastic Net


  data <- gm12878_10kb_reduced
  
  #splitting the data
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
  enetlst[[1]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst[[2]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst[[3]] <- varImp(eNetModel)$importance[,1]
  enetlst[[4]] <- coef(eNetModel$finalModel, eNetModel$bestTune$lambda)[-1,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                                     newdata=test,
                                     type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf[1,1] <- TN
  enetperf[2,1] <- FN
  enetperf[3,1] <- FP
  enetperf[4,1] <- TP
  enetperf[5,1] <- sum(confMat$table)
  enetperf[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf[7,1] <- as.vector(confMat$byClass["Specificity"])
  enetperf[8,1] <- as.vector(confMat$overall["Kappa"])
  enetperf[9,1] <- as.vector(confMat$overall["Accuracy"])
  enetperf[10,1] <- TP/(TP+FP)
  enetperf[11,1] <- FP/(FP+TN)
  enetperf[12,1] <- FN/(FN+TN)
  enetperf[13,1] <- FN/(FN+TN)
  enetperf[14,1] <- TN/(TN+FN)
  enetperf[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf[17,1] <- pROC::auc(pROC::roc(test$y, pred.enetModel))




saveRDS(enetlst,"/home/stilianoudakisc/TAD_data_analysis/final_models/enetlst_unb.rds")
saveRDS(enetperf,"/home/stilianoudakisc/TAD_data_analysis/final_models/enetperf_unb.rds")


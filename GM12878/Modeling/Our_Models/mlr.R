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


gm12878_5kb <- readRDS("gm12878_5kb.rds")

#set tuning parameters
#fitControl <- trainControl(method = "cv",
#                           number = 10,
#                           ## Estimate class probabilities
#                           classProbs = TRUE,
#                           ## Evaluate performance using 
#                           ## the following function
#                           summaryFunction = twoClassSummary)
						   
#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


						   
#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
glmlst <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_5kb$y=="Yes"))*2)*.3), 
                              ncol=1),
                fpr <- matrix(nrow=ceiling((length(which(gm12878_5kb$y=="Yes"))*2)*.3), 
                              ncol=1),
                varimp <- matrix(nrow=dim(gm12878_5kb)[2]-1,
                                 ncol=1),
				coefs <- matrix(nrow=dim(gm12878_5kb)[2]-1,
									ncol=1))
rownames(glmlst[[3]]) <- colnames(gm12878_5kb)[-1]
rownames(glmlst[[4]]) <- colnames(gm12878_5kb)[-1]				 
								 
glmperf <- matrix(nrow = 17, ncol=1)
rownames(glmperf) <- c("TN",
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
							

#Performing MLR

data <- gm12878_5kb
  
  set.seed(123)
  #splitting the data
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #GLM Model
  #glmModel <- train(y ~ ., data=train, 
  #                   method = "glm", 
  #                   metric="ROC", 
  #                   trControl = fitControl, 
  #                   family="binomial")
					 
  glmModel <- glm(y ~ ., data=train, family=binomial(link="logit"))
					 
  #Prediction vector for ROC and AUC
  pred.glmModel <- predict(glmModel, newdata = test, type = "response")
  pred.logit <- ifelse(pred.glmModel>.50,1,0)
  #pred.glmModel <- predict(glmModel, 
  #                                  newdata=test, 
  #                                  type="response")
  glmlst[[1]] <- simple_roc(test$y,pred.logit)[,1]
  glmlst[[2]] <- simple_roc(test$y,pred.logit)[,2]
  glmlst[[3]] <- varImp(glmModel)$importance[,1]
  glmlst[[4]] <- coef(glmModel$finalModel)[-1,1]
  
  #Prediction vector for other performance metrics
  confMat <- confusionMatrix(pred.logit, test$y, positive="1")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  glmperf[1,1] <- TN
  glmperf[2,1] <- FN
  glmperf[3,1] <- FP
  glmperf[4,1] <- TP
  glmperf[5,1] <- sum(confMat$table)
  glmperf[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  glmperf[7,1] <- as.vector(confMat$byClass["Specificity"])
  glmperf[8,1] <- as.vector(confMat$overall["Kappa"])
  glmperf[9,1] <- as.vector(confMat$overall["Accuracy"])
  glmperf[10,1] <- TP/(TP+FP)
  glmperf[11,1] <- FP/(FP+TN)
  glmperf[12,1] <- FN/(FN+TN)
  glmperf[13,1] <- FN/(FN+TN)
  glmperf[14,1] <- TN/(TN+FN)
  glmperf[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  glmperf[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  glmperf[17,1] <- pROC::auc(pROC::roc(test$y, pred.logit))




saveRDS(glmlst,"/home/stilianoudakisc/TAD_data_analysis/final_models/glmlst.rds")
saveRDS(glmperf,"/home/stilianoudakisc/TAD_data_analysis/final_models/glmperf.rds")




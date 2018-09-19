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
fitControl <- trainControl(method = "cv",
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)
						   
#set up grid of alpha and lambda values
lambda.grid = 10^seq(2,-2,length=100)
alpha.grid = 1
srchGrid <- expand.grid(.alpha=alpha.grid, .lambda=lambda.grid)					   

#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
lassolst <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_5kb$y=="Yes"))*2)*.3), 
                              ncol=1),
                fpr <- matrix(nrow=ceiling((length(which(gm12878_5kb$y=="Yes"))*2)*.3), 
                              ncol=1),
                varimp <- matrix(nrow=dim(gm12878_5kb)[2]-1,
                                 ncol=1),
				coefs <- matrix(nrow=dim(gm12878_5kb)[2]-1,
									ncol=1))
rownames(lassolst[[3]]) <- colnames(gm12878_5kb)[-1]
rownames(lassolst[[4]]) <- colnames(gm12878_5kb)[-1]				 
								 
lassoperf <- matrix(nrow = 17, ncol=1)
rownames(lassoperf) <- c("TN",
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
							

#Performing MLR w/ LASSO

data <- gm12878_5kb
data$y <- factor(ifelse(data$y==1, "Yes", "No"))
  
  set.seed(123)
  #splitting the data
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #LASSO Model

  lassoModel <- train(y ~ ., data=train,
					method = "glmnet", 
					trControl = fitControl,
					metric = "ROC",
					tuneGrid = srchGrid,
					standardize = TRUE)
		
  #Prediction vector for ROC and AUC
  pred.lassoModel <- as.vector(predict(lassoModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  lassolst[[1]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.lassoModel)[,1]
  lassolst[[2]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.lassoModel)[,2]
  lassolst[[3]] <- varImp(lassoModel)$importance[,1]
  lassolst[[4]] <- coef(lassoModel$finalModel)[-1,1]
  
  #Prediction vector for other performance metrics
  pred.lassoModel2 <- predict(lassoModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.lassoModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  lassoperf[1,1] <- TN
  lassoperf[2,1] <- FN
  lassoperf[3,1] <- FP
  lassoperf[4,1] <- TP
  lassoperf[5,1] <- sum(confMat$table)
  lassoperf[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  lassoperf[7,1] <- as.vector(confMat$byClass["Specificity"])
  lassoperf[8,1] <- as.vector(confMat$overall["Kappa"])
  lassoperf[9,1] <- as.vector(confMat$overall["Accuracy"])
  lassoperf[10,1] <- TP/(TP+FP)
  lassoperf[11,1] <- FP/(FP+TN)
  lassoperf[12,1] <- FN/(FN+TN)
  lassoperf[13,1] <- FN/(FN+TN)
  lassoperf[14,1] <- TN/(TN+FN)
  lassoperf[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  lassoperf[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  lassoperf[17,1] <- pROC::auc(pROC::roc(test$y, pred.lassoModel))




saveRDS(lassolst,"/home/stilianoudakisc/TAD_data_analysis/final_models/lassolst.rds")
saveRDS(lassoperf,"/home/stilianoudakisc/TAD_data_analysis/final_models/lassoperf.rds")


 		

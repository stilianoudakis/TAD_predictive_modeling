library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)

setwd("/home/stilianoudakisc/TAD_data_analysis/comparing_normalization/")

gm12878_5kb_f <- readRDS("gm12878_5kb_f.rds")

						   
#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


########################################################################################


#With log transform and standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_lns <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_5kb_f$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   fpr <- matrix(nrow=ceiling((length(which(gm12878_5kb_f$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   varimp <- matrix(nrow=dim(gm12878_5kb_f)[2]-1,
                                    ncol=1))
rownames(enetlst_lns[[3]]) <- colnames(gm12878_5kb_f)[-1]

enetperf_ls <- matrix(nrow = 17, ncol=1)
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
						   
data <- gm12878_5kb_f

#split data
inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
train <- data[inTrainingSet,]
test <- data[-inTrainingSet,]

#set tuning parameters
fitControl <- trainControl(method = "cv",
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)
						   
#set up grid of alpha and lambda values
lambda.grid = 10^seq(2,-2,length=100)
alpha.grid = seq(0,1,length=10)
srchGrid <- expand.grid(.alpha=alpha.grid, .lambda=lambda.grid)

#ENET Model
  eNetModel <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneGrid = srchGrid,
                     standardize=TRUE)
					 
#Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_lns[[1]]<- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_lns[[2]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_lns[[3]] <- varImp(eNetModel)$importance[,1]
  
#Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_ls[1,1] <- TN
  enetperf_ls[2,1] <- FN
  enetperf_ls[3,1] <- FP
  enetperf_ls[4,1] <- TP
  enetperf_ls[5,1] <- sum(confMat$table)
  enetperf_ls[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_ls[7,1] <- as.vector(confMat$byClass["Specificity"])
  enetperf_ls[8,1] <- as.vector(confMat$overall["Kappa"])
  enetperf_ls[9,1] <- as.vector(confMat$overall["Accuracy"])
  enetperf_ls[10,1] <- TP/(TP+FP)
  enetperf_ls[11,1] <- FP/(FP+TN)
  enetperf_ls[12,1] <- FN/(FN+TN)
  enetperf_ls[13,1] <- FN/(FN+TN)
  enetperf_ls[14,1] <- TN/(TN+FN)
  enetperf_ls[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_ls[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_ls[17,1] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  


saveRDS(enetlst_lns, "enetlst_lns.rds")
saveRDS(enetperf_ls, "enetperf_ls.rds")


###################################################################################



					 



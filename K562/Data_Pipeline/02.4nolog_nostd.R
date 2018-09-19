library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)

setwd("/home/stilianoudakisc/TAD_data_analysis/K562/comparing_normalization/")

k562_10kb_f <- readRDS("k562_10kb_f.rds")

k562_10kb_f_nl <- k562_10kb_f

cols <- c(grep("dist",colnames(k562_10kb_f_nl)))
k562_10kb_f_nl[,cols] <- apply(k562_10kb_f_nl[,cols], 2, function(x){2^x})


#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

########################################################################################


#With No log transform and No standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_nlns <- list(tpr <- matrix(nrow=ceiling((length(which(k562_10kb_f_nl$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   fpr <- matrix(nrow=ceiling((length(which(k562_10kb_f_nl$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   varimp <- matrix(nrow=dim(k562_10kb_f_nl)[2]-1,
                                    ncol=1))
rownames(enetlst_nlns[[3]]) <- colnames(k562_10kb_f_nl)[-1]

enetperf_nlns <- matrix(nrow = 17, ncol=1)
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
						   
data <- k562_10kb_f_nl

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
                     standardize=FALSE)
					 
#Prediction vector for ROC and AUC
  pred.eNetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_nlns[[1]]<- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,1]
  enetlst_nlns[[2]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.eNetModel)[,2]
  enetlst_nlns[[3]] <- varImp(eNetModel)$importance[,1]
  
#Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_nlns[1,1] <- TN
  enetperf_nlns[2,1] <- FN
  enetperf_nlns[3,1] <- FP
  enetperf_nlns[4,1] <- TP
  enetperf_nlns[5,1] <- sum(confMat$table)
  enetperf_nlns[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_nlns[7,1] <- as.vector(confMat$byClass["Specificity"])
  enetperf_nlns[8,1] <- as.vector(confMat$overall["Kappa"])
  enetperf_nlns[9,1] <- as.vector(confMat$overall["Accuracy"])
  enetperf_nlns[10,1] <- TP/(TP+FP)
  enetperf_nlns[11,1] <- FP/(FP+TN)
  enetperf_nlns[12,1] <- FN/(FN+TN)
  enetperf_nlns[13,1] <- FN/(FN+TN)
  enetperf_nlns[14,1] <- TN/(TN+FN)
  enetperf_nlns[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_nlns[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_nlns[17,1] <- pROC::auc(pROC::roc(test$y, pred.eNetModel))
  


saveRDS(enetlst_nlns, "enetlst_nlns.rds")
saveRDS(enetperf_nlns, "enetperf_nlns.rds")


###################################################################################


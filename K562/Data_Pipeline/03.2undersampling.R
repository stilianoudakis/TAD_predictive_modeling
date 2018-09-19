#Elastic Net with balanced classes using bootstrap sampling 

library(caret)
#library(data.table)
library(gbm)
library(pROC)
library(plyr)
library(dplyr)
#library(DMwR)
library(gridExtra)
library(ggplot2)

setwd("/home/stilianoudakisc/TAD_data_analysis/K562/evaluating_class_imbalance/")

k562_10kb_f <- readRDS("k562_10kb_f.rds")

#### function for roc curves

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


#Randomly undersample from majority class to get balanced classes
set.seed(7215)
sampids <- sample(which(k562_10kb_f$y=="No"),
                        length(which(k562_10kb_f$y=="Yes")),
                        replace = TRUE)
						
#combining the two classes to create balanced data
  data <- rbind.data.frame(k562_10kb_f[which(k562_10kb_f$y=="Yes"),],
                           k562_10kb_f[sampids,])
						   
#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_b <- list(tpr <- matrix(nrow=ceiling((length(which(data$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   fpr <- matrix(nrow=ceiling((length(which(data$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   varimp <- matrix(nrow=dim(data)[2]-1,
                                    ncol=1))
rownames(enetlst_b[[3]]) <- colnames(data)[-1]


enetperf_b <- matrix(nrow = 17, ncol=1)
rownames(enetperf_b) <- c("TN",
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

#splitting data	
set.seed(7215)					  
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
  enetModel_b <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneGrid=srchGrid,
                     standardize=TRUE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(enetModel_b, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_b[[1]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_b[[2]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_b[[3]] <- varImp(enetModel_b)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_b,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_b[1,1] <- TN
  enetperf_b[2,1] <- FN
  enetperf_b[3,1] <- FP
  enetperf_b[4,1] <- TP
  enetperf_b[5,1] <- sum(confMat$table)
  enetperf_b[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_b[7,1] <- as.vector(confMat$byClass["Specificity"])
  enetperf_b[8,1] <- as.vector(confMat$overall["Kappa"])
  enetperf_b[9,1] <- as.vector(confMat$overall["Accuracy"])
  enetperf_b[10,1] <- TP/(TP+FP)
  enetperf_b[11,1] <- FP/(FP+TN)
  enetperf_b[12,1] <- FN/(FN+TN)
  enetperf_b[13,1] <- FN/(FN+TN)
  enetperf_b[14,1] <- TN/(TN+FN)
  enetperf_b[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_b[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_b[17,1] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  
saveRDS(enetlst_b, "enetlst_bs.rds")

saveRDS(enetperf_b, "enetperf_bs.rds")

  
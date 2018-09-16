#Testing SMOTE combinations with Elastic Net

library(caret)
library(data.table)
library(gbm)
library(glmnet)
library(pROC)
library(plyr)
library(dplyr)
library(DMwR)
library(gridExtra)


setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/")

gm12878_5kb_f <- readRDS("gm12878_5kb_f.rds")

# Splitting the data
set.seed(5228)
inTrainingSet <- createDataPartition(gm12878_5kb_f$y,p=.7,list=FALSE)
train <- gm12878_5kb_f[inTrainingSet,]
test <- gm12878_5kb_f[-inTrainingSet,]

#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

# Establishing tuning/training parameters
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


#########################################################################################


#SMOTE

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_sm <- list(tpr <- matrix(nrow=length(test$y), 
                                 ncol=4),
                   fpr <- matrix(nrow=length(test$y), 
                                 ncol=4),
                   varimp <- matrix(nrow=dim(gm12878_5kb_f)[2]-1,
                                    ncol=4))
rownames(enetlst_sm[[3]]) <- colnames(gm12878_5kb_f)[-1]

enetperf_sm <- matrix(nrow = 17, ncol=4)
rownames(enetperf_sm) <- c("TN",
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
						   

for(i in 1:4){
  set.seed(112)
  #100/200, 200/200, 300/200, 400/200
  train_smote <- SMOTE(y ~ ., 
                       data=train, 
                       perc.over = i*100, 
                       perc.under = 200)
  
  #ENET Model
  enetModel_sm <- train(y ~ ., data=train_smote, 
                        method = "glmnet", 
                        metric="ROC", 
                        trControl = fitControl,
                        family="binomial",
                        tuneGrid = srchGrid,
                        standardize=TRUE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(enetModel_sm, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_sm[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_sm[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_sm[[3]][,i] <- varImp(enetModel_sm)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_sm,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_sm[1,i] <- TN
  enetperf_sm[2,i] <- FN
  enetperf_sm[3,i] <- FP
  enetperf_sm[4,i] <- TP
  enetperf_sm[5,i] <- sum(confMat$table)
  enetperf_sm[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_sm[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_sm[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_sm[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_sm[10,i] <- TP/(TP+FP)
  enetperf_sm[11,i] <- FP/(FP+TN)
  enetperf_sm[12,i] <- FN/(FN+TN)
  enetperf_sm[13,i] <- FN/(FN+TN)
  enetperf_sm[14,i] <- TN/(TN+FN)
  enetperf_sm[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_sm[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_sm[17,i] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  
 }
 
saveRDS(enetlst_sm, "enetlst_sm_lns.rds")

saveRDS(enetperf_sm, "enetperf_sm.rds")



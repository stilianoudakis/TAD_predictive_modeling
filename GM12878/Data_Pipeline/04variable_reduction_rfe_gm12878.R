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

setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/")

#### Forward Selection

gm12878_10kb_f <- readRDS("gm12878_10kb_f.rds")

####set number of bootstrap samples

#bootsamps = 5

#### function for roc curves

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


#### Establishing tuning/training parameters

fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

####create a matrix of row ids that represent the zero class

#sampids <- matrix(ncol=bootsamps, 
#                  nrow=length(gm12878_10kb_f$y[which(gm12878_10kb_f$y=="Yes")]))
				  
####filling in the sample ids matrix

#set.seed(123)
#for(j in 1:bootsamps){
#  sampids[,j] <- sample(which(gm12878_10kb_f$y=="No"),
#                        length(which(gm12878_10kb_f$y=="Yes")),
#                        replace = TRUE)
#}


#### Recursive Feature Elimination


#setting rfe parameters
control <- rfeControl(functions=rfFuncs, method="cv", number=10, repeats=5)

trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummary)

#namesmat <- matrix(NA,nrow=dim(gm12878_10kb_f)[2]-1, ncol=bootsamps)
#rownames(namesmat) <- sort(names(gm12878_10kb_f)[-which(names(gm12878_10kb_f)=="y")])


#for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  #data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
  #                         gm12878_10kb_f[sampids[,i],])
  
  data <- gm12878_10kb_f
  
  rfeModel <- rfe(data[,-1], 
                data[,1], 
                sizes=c(2:50), 
                metric="ROC",
                rfeControl=control,
                trControl = trainctrl)
  
  preds <- predictors(rfeModel)

  #namesmat[which(rownames(namesmat) %in% preds),i] <- sort(preds)
  
#}

#namesmat <- data.frame(namesmat)

#namesmat$naperc <- (rowSums(!is.na(namesmat))/5)*100

#rfe.preds <- rownames(namesmat)[which(namesmat$naperc >= 90)]

gm12878_10kb_rfe <- gm12878_10kb_f[,which((names(gm12878_10kb_f) %in% preds) | names(gm12878_10kb_f)=="y")]

saveRDS(gm12878_10kb_rfe,"gm12878_10kb_rfe.rds")


#Evaluating performance of reduced data

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
#enetlst_rfe <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_rfe$y=="Yes"))*2)*.3), 
#                                  ncol=bootsamps),
#                    fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_rfe$y=="Yes"))*2)*.3), 
#                                  ncol=bootsamps),
#                    varimp <- matrix(nrow=dim(gm12878_10kb_rfe)[2]-1,
#                                     ncol=bootsamps))
#rownames(enetlst_rfe[[3]]) <- colnames(gm12878_10kb_rfe)[-1]
#
#enetperf_rfe <- matrix(nrow = 17, ncol=bootsamps)
#rownames(enetperf_rfe) <- c("TN",
#                            "FN",
#                            "FP",
#                            "TP",
#                            "Total",
#                            "Sensitivity",
#                            "Specificity",
#                            "Kappa",
#                            "Accuracy",
#                            "Precision",
#                            "FPR",
#                            "FNR",
#                            "FOR",
#                            "NPV",
#                            "MCC",
#                            "F1",
#                            "AUC")
#
#for(i in 1:bootsamps){
#  set.seed(7215)
#  #combining the two classes to create balanced data
#  data <- rbind.data.frame(gm12878_10kb_rfe[which(gm12878_10kb_rfe$y=="Yes"),],
#                           gm12878_10kb_rfe[sampids[,i],])
#  
#  
#  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
#  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
#  train <- data[inTrainingSet,]
#  test <- data[-inTrainingSet,]
#  
#  #ENET Model
#  eNetModel <- train(y ~ ., data=train, 
#                     method = "glmnet", 
#                     metric="ROC", 
#                     trControl = fitControl, 
#                     family="binomial", 
#                     tuneLength=5,
#                    standardize=FALSE)
#  
#  #Prediction vector for ROC and AUC					  
#  pred.enetModel <- as.vector(predict(eNetModel, 
#                                      newdata=test, 
#                                      type="prob")[,"Yes"])
#  enetlst_rfe[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
#  enetlst_rfe[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
#  enetlst_rfe[[3]][,i] <- varImp(eNetModel)$importance[,1]
#  
#  #Prediction vector for other performance metrics
#  pred.enetModel2 <- predict(eNetModel,
#                             newdata=test,
#                             type="raw")
#  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
#  TN = as.numeric(confMat$table[1,1])
#  FN = as.numeric(confMat$table[1,2])
#  FP = as.numeric(confMat$table[2,1])
#  TP = as.numeric(confMat$table[2,2])
#  TP = as.numeric(confMat$table[2,2])
#  enetperf_rfe[1,i] <- TN
#  enetperf_rfe[2,i] <- FN
#  enetperf_rfe[3,i] <- FP
#  enetperf_rfe[4,i] <- TP
#  enetperf_rfe[5,i] <- sum(confMat$table)
#  enetperf_rfe[6,i] <- as.vector(confMat$byClass["Sensitivity"])
#  enetperf_rfe[7,i] <- as.vector(confMat$byClass["Specificity"])
#  enetperf_rfe[8,i] <- as.vector(confMat$overall["Kappa"])
#  enetperf_rfe[9,i] <- as.vector(confMat$overall["Accuracy"])
#  enetperf_rfe[10,i] <- TP/(TP+FP)
#  enetperf_rfe[11,i] <- FP/(FP+TN)
#  enetperf_rfe[12,i] <- FN/(FN+TN)
#  enetperf_rfe[13,i] <- FN/(FN+TN)
#  enetperf_rfe[14,i] <- TN/(TN+FN)
#  enetperf_rfe[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
#  enetperf_rfe[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + #((TP/(TP+FP))))
#  enetperf_rfe[17,i] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
#}
#
#saveRDS(enetlst_rfe,"enetlst_rfe.rds")
#saveRDS(enetperf_rfe,"enetperf_rfe.rds")


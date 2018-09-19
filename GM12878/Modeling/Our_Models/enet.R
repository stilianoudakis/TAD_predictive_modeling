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


gm12878_5kb_reduced <- readRDS("gm12878_5kb_reduced.rds")

### set number of bootstrap samples


bootsamps = 5

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
					   
#### function for roc curves

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}
						
###create a matrix of row ids that represent the zero class

##the number of rows will match the one class
##the number of columns match the number of bootstrap samples

sampids <- matrix(ncol=bootsamps, 
                  nrow=length(gm12878_5kb_reduced$y[which(gm12878_5kb_reduced$y=="Yes")]))


###filling in the sample ids matrix

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(gm12878_5kb_reduced$y=="No"),
                        length(which(gm12878_5kb_reduced$y=="Yes")),
                        replace = TRUE)
}



# Final Models:

# Elastic Net

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_5kb_reduced$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(gm12878_5kb_reduced$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    varimp <- matrix(nrow=dim(gm12878_5kb_reduced)[2]-1,
                                     ncol=bootsamps),
					coefs <- matrix(nrow=dim(gm12878_5kb_reduced)[2]-1,
									ncol=bootsamps))
rownames(enetlst[[3]]) <- colnames(gm12878_5kb_reduced)[-1]
rownames(enetlst[[4]]) <- colnames(gm12878_5kb_reduced)[-1]

enetperf <- matrix(nrow = 17, ncol=bootsamps)
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

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_5kb_reduced[which(gm12878_5kb_reduced$y=="Yes"),],
                           gm12878_5kb_reduced[sampids[,i],])
  
  
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
                     tuneGrid=srchGrid,
                     standardize=TRUE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(eNetModel, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst[[3]][,i] <- varImp(eNetModel)$importance[,1]
  enetlst[[4]][,i] <- coef(eNetModel$finalModel, eNetModel$bestTune$lambda)[-1,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(eNetModel,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf[1,i] <- TN
  enetperf[2,i] <- FN
  enetperf[3,i] <- FP
  enetperf[4,i] <- TP
  enetperf[5,i] <- sum(confMat$table)
  enetperf[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf[10,i] <- TP/(TP+FP)
  enetperf[11,i] <- FP/(FP+TN)
  enetperf[12,i] <- FN/(FN+TN)
  enetperf[13,i] <- FN/(FN+TN)
  enetperf[14,i] <- TN/(TN+FN)
  enetperf[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf[17,i] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
}

saveRDS(enetlst,"enetlst.rds")
saveRDS(enetperf,"enetperf.rds")

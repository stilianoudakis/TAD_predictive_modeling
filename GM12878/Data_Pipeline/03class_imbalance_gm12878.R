#Testing SMOTE combinations for Elastic Net


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

gm12878_10kb_f <- readRDS("gm12878_10kb_f.rds")

# Splitting the data
set.seed(5228)
inTrainingSet <- createDataPartition(gm12878_10kb_f$y,p=.7,list=FALSE)
train <- gm12878_10kb_f[inTrainingSet,]
test <- gm12878_10kb_f[-inTrainingSet,]

# Establishing tuning/training parameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#Testing eight different combinations of perc.over/perc.under
#100/200, 200/200, 300/200, 400/200, 
#100/300, 200/300, 300/300, 400/300

#########################################################################################


#SMOTE

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_sm <- list(tpr <- matrix(nrow=length(test$y), 
                                 ncol=4),
                   fpr <- matrix(nrow=length(test$y), 
                                 ncol=4),
                   varimp <- matrix(nrow=dim(gm12878_10kb_f)[2]-1,
                                    ncol=4))
rownames(enetlst_sm[[3]]) <- colnames(gm12878_10kb_f)[-1]

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
                        tuneLength=5,
                        standardize=FALSE)
  
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
  
  
  #########################################################################################
  
  #set.seed(112)
  #100/300, 200/300, 300/300, 400/300
  #train_smote <- SMOTE(y ~ ., 
  #                     data=train, 
  #                     perc.over = i*100, 
  #                     perc.under = 300)
 # 
  #ENET Model
  #enetModel_sm <- train(y ~ ., data=train_smote, 
  #                      method = "glmnet", 
  #                      metric="ROC", 
  #                      trControl = fitControl,
  #                      family="binomial",
  #                      tuneLength=5,
  #                      standardize=FALSE)
  
  #Prediction vector for ROC and AUC
  #pred.enetModel <- as.vector(predict(enetModel_sm, 
  #                                    newdata=test, 
  #                                    type="prob")[,"Yes"])
  #enetlst_sm[[1]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  #enetlst_sm[[2]][,i+4] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  #enetlst_sm[[3]][,i+4] <- varImp(enetModel_sm)$importance[,1]
  
  #Prediction vector for other performance metrics
  #pred.enetModel2 <- predict(enetModel_sm,
  #                           newdata=test,
  #                           type="raw")
  #confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  #TN = as.numeric(confMat$table[1,1])
  #FN = as.numeric(confMat$table[1,2])
  #FP = as.numeric(confMat$table[2,1])
  #TP = as.numeric(confMat$table[2,2])
  #enetperf_sm[1,i+4] <- TN
  #enetperf_sm[2,i+4] <- FN
  #enetperf_sm[3,i+4] <- FP
  #enetperf_sm[4,i+4] <- TP
  #enetperf_sm[5,i+4] <- sum(confMat$table)
  #enetperf_sm[6,i+4] <- as.vector(confMat$byClass["Sensitivity"])
  #enetperf_sm[7,i+4] <- as.vector(confMat$byClass["Specificity"])
  #enetperf_sm[8,i+4] <- as.vector(confMat$overall["Kappa"])
  #enetperf_sm[9,i+4] <- as.vector(confMat$overall["Accuracy"])
  #enetperf_sm[10,i+4] <- TP/(TP+FP)
  #enetperf_sm[11,i+4] <- FP/(FP+TN)
  #enetperf_sm[12,i+4] <- FN/(FN+TN)
  #enetperf_sm[13,i+4] <- FN/(FN+TN)
  #enetperf_sm[14,i+4] <- TN/(TN+FN)
  #enetperf_sm[15,i+4] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  #enetperf_sm[16,i+4] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  #enetperf_sm[17,i+4] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  
}

saveRDS(enetlst_sm, "/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/enetlst_sm_lns.rds")

saveRDS(enetperf_sm, "/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/enetperf_sm.rds")



#### Evaluating SMOTE

#Plotting Performance
#auc.sm <- data.frame(Combination=c("100/200","200/200","300/200","400/200",
#                                     "100/300","200/300","300/300","400/300"),
#                       AUC=c(enetperf_sm[17,1],enetperf_sm[17,2],enetperf_sm[17,3],enetperf_sm[17,4],
#                             enetperf_sm[17,5],enetperf_sm[17,6],enetperf_sm[17,7],enetperf_sm[17,8]))

#auc.sm <- auc.sm[order(auc.sm$AUC, decreasing=TRUE),]

#auc.sm$Combination <- factor(auc.sm$Combination, levels=auc.sm$Combination)

#auc.sm

#saveRDS(auc.sm, "/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/auc.sm.rds")


#auc.sm.p<-ggplot(data=auc.sm, aes(x=Combination, y=AUC, fill=Combination)) + 
#  xlab("Sampling Combination") + ylab("AUC") +
#  geom_bar(stat="identity") + ylim(0,1) +
#  scale_fill_manual(values=gray(seq(0,.7,.1)), guide=FALSE) +
#  theme_minimal() + 
#  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#  ggtitle("Model Performance for Different \n Sampling Combinations using SMOTE")

#auc.sm.p
#ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/auc.sm.p.png")

#onetwo <- data.frame(fpr=enetlst_sm[[2]][,1],tpr=enetlst_sm[[1]][,1], Combo = "100/200");
#twotwo <- data.frame(fpr=enetlst_sm[[2]][,2],tpr=enetlst_sm[[1]][,2], Combo = "200/200");
#threetwo <- data.frame(fpr=enetlst_sm[[2]][,3],tpr=enetlst_sm[[1]][,3], Combo = "300/200");
#fourtwo <- data.frame(fpr=enetlst_sm[[2]][,4],tpr=enetlst_sm[[1]][,4], Combo = "400/200");
#onethree <- data.frame(fpr=enetlst_sm[[2]][,5],tpr=enetlst_sm[[1]][,5], Combo = "100/300");
#twothree <- data.frame(fpr=enetlst_sm[[2]][,6],tpr=enetlst_sm[[1]][,6], Combo = "200/300");
#threethree <- data.frame(fpr=enetlst_sm[[2]][,7],tpr=enetlst_sm[[1]][,7], Combo = "300/300");
#fourthree <- data.frame(fpr=enetlst_sm[[2]][,8],tpr=enetlst_sm[[1]][,8], Combo = "400/300")

#allrocdat <- rbind.data.frame(onetwo,
#                              twotwo,
#                              threetwo,
#                              fourtwo,
#                              onethree,
#                              twothree,
#                              threethree,
#                              fourthree)

#roc.sm <- ggplot(data=allrocdat, aes(x=fpr, y=tpr, color=Combo)) + 
#  geom_line(size=1) +
#  scale_colour_manual(name="Combination",
#    labels=c("100/200", 
#             "200/200",
#            "300/200",
#             "400/200",
#             "100/300",
#             "200/300",
#             "300/300",
#             "400/300"),
#   values=gray(seq(0,.7,.1))) + 
#  xlab("1-Specificity") + 
#  ylab("Sensitivity") + 
#  xlim(0, 1) +
#  ylim(0, 1) +
#  geom_abline(intercept=0, slope=1) +
#  theme_minimal() +
#  ggtitle("ROC Curves for Different \n Normalization Techniques")

#roc.sm
#ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.sm.png")

#####################################################################################33

#Bootstrap samples

#set number of bootstrap samples
bootsamps = 5

#set tuning parameters
fitControl <- trainControl(method = "repeatedcv",
                           number =5,
                           repeats = 3,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary)

#create a matrix of row ids that represent the zero class
#the number of rows will match the one class
#the number of columns match the number of bootstrap samples
sampids <- matrix(ncol=bootsamps, 
                  nrow=length(gm12878_10kb_f$y[which(gm12878_10kb_f$y=="Yes")]))


#filling in the sample ids matrix
set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(gm12878_10kb_f$y=="No"),
                        length(which(gm12878_10kb_f$y=="Yes")),
                        replace = TRUE)
}



#function for roc curves
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}



#With log transform and standardization

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_b <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                 ncol=bootsamps),
                   varimp <- matrix(nrow=dim(gm12878_10kb_f)[2]-1,
                                    ncol=bootsamps))
rownames(enetlst_b[[3]]) <- colnames(gm12878_10kb_f)[-1]


enetperf_b <- matrix(nrow = 17, ncol=bootsamps)
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

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_10kb_f[which(gm12878_10kb_f$y=="Yes"),],
                           gm12878_10kb_f[sampids[,i],])
  
  
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #ENET Model
  enetModel_b <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=TRUE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(enetModel_b, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_b[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_b[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_b[[3]][,i] <- varImp(enetModel_b)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_b,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_b[1,i] <- TN
  enetperf_b[2,i] <- FN
  enetperf_b[3,i] <- FP
  enetperf_b[4,i] <- TP
  enetperf_b[5,i] <- sum(confMat$table)
  enetperf_b[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_b[7,i] <- as.vector(confMat$byClass["Specificity"])
  enetperf_b[8,i] <- as.vector(confMat$overall["Kappa"])
  enetperf_b[9,i] <- as.vector(confMat$overall["Accuracy"])
  enetperf_b[10,i] <- TP/(TP+FP)
  enetperf_b[11,i] <- FP/(FP+TN)
  enetperf_b[12,i] <- FN/(FN+TN)
  enetperf_b[13,i] <- FN/(FN+TN)
  enetperf_b[14,i] <- TN/(TN+FN)
  enetperf_b[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_b[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_b[17,i] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  
}


saveRDS(enetlst_b, "/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/enetlst_bs.rds")

saveRDS(enetperf_b, "/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/enetperf_bs.rds")



#Mean AUC across n bootstrap samples
mean(enetperf_b[17,])

auc.bs <- round(mean(enetperf_b[17,]),3)
auc.bs

#roc curve
#fpr.bs <- rowMeans(enetlst_b[[2]])
#tpr.bs <- rowMeans(enetlst_b[[1]])
#rocdat.bs <- data.frame(fpr=fpr.bs, tpr=tpr.bs)

#roc.bs <- ggplot(rocdat.bs, aes(x=fpr, y=tpr)) + 
#  geom_line(size=1, color="black") +
#  xlab("1-Specificity") + 
#  ylab("Sensitivity") + 
#  xlim(0, 1) +
#  ylim(0, 1) +
#  geom_abline(intercept=0, slope=1) +
#  theme_minimal() +
#  ggtitle("ROC Curve for Balanced Classes \n Using n Bootstrap Samples")

#roc.bs
#ggsave("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb_bins/roc.bs.png")


############################################################################################


#Unbalanced data


#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
enetlst_unb <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   fpr <- matrix(nrow=ceiling((length(which(gm12878_10kb_f$y=="Yes"))*2)*.3), 
                                 ncol=1),
                   varimp <- matrix(nrow=dim(gm12878_10kb_f)[2]-1,
                                    ncol=1))
rownames(enetlst_unb[[3]]) <- colnames(gm12878_10kb_f)[-1]


enetperf_unb <- matrix(nrow = 17, ncol=1)
rownames(enetperf_unb) <- c("TN",
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


 
  
  #ENET Model
  enetModel_unb <- train(y ~ ., data=train, 
                     method = "glmnet", 
                     metric="ROC", 
                     trControl = fitControl, 
                     family="binomial", 
                     tuneLength=5,
                     standardize=FALSE)
  
  #Prediction vector for ROC and AUC					  
  pred.enetModel <- as.vector(predict(enetModel_unb, 
                                      newdata=test, 
                                      type="prob")[,"Yes"])
  enetlst_unb[[1]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,1]
  enetlst_unb[[2]] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.enetModel)[,2]
  enetlst_unb[[3]] <- varImp(enetModel_unb)$importance[,1]
  
  #Prediction vector for other performance metrics
  pred.enetModel2 <- predict(enetModel_unb,
                             newdata=test,
                             type="raw")
  confMat <- confusionMatrix(data=pred.enetModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  enetperf_unb[1,1] <- TN
  enetperf_unb[2,1] <- FN
  enetperf_unb[3,1] <- FP
  enetperf_unb[4,1] <- TP
  enetperf_unb[5,1] <- sum(confMat$table)
  enetperf_unb[6,1] <- as.vector(confMat$byClass["Sensitivity"])
  enetperf_unb[7,1] <- as.vector(confMat$byClass["Specificity"])
  enetperf_unb[8,1] <- as.vector(confMat$overall["Kappa"])
  enetperf_unb[9,1] <- as.vector(confMat$overall["Accuracy"])
  enetperf_unb[10,1] <- TP/(TP+FP)
  enetperf_unb[11,1] <- FP/(FP+TN)
  enetperf_unb[12,1] <- FN/(FN+TN)
  enetperf_unb[13,1] <- FN/(FN+TN)
  enetperf_unb[14,1] <- TN/(TN+FN)
  enetperf_unb[15,1] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  enetperf_unb[16,1] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  enetperf_unb[17,1] <- pROC::auc(pROC::roc(test$y, pred.enetModel))
  



saveRDS(enetlst_unb, "/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/enetlst_unb.rds")

saveRDS(enetperf_unb, "/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/enetperf_unb.rds")



#Mean AUC across n bootstrap samples
enetperf_unb[17,1]

auc.bs <- round(enetperf_unb[17,1],3)
auc.bs

#roc curve
#fpr.unb <- enetlst_unb[[2]]
#tpr.unb <- enetlst_unb[[1]]
#rocdat.unb <- data.frame(fpr=fpr.unb, tpr=tpr.unb)

#roc.unb <- ggplot(rocdat.unb, aes(x=fpr, y=tpr)) + 
#  geom_line(size=1, color="black") +
#  xlab("1-Specificity") + 
#  ylab("Sensitivity") + 
#  xlim(0, 1) +
#  ylim(0, 1) +
#  geom_abline(intercept=0, slope=1) +
#  theme_minimal() +
#  ggtitle("ROC Curve for Balanced Classes \n Using n Bootstrap Samples")

#roc.bs
#ggsave("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/roc.unb.png")


###########################################################################################

### Comparing additional performance metrics between SMOTE and Bootstrap methods

options(scipen = 999)

enetperf_sm <- round(enetperf_sm,3)
enetperf_b <- round(as.matrix(rowMeans(enetperf_b)),3)
enetperf_b[1:5,1] <- round(enetperf_b[1:5,1],0)
enetperf_unb <- round(enetperf_unb, 3)

perfdat_sm_bs <- cbind.data.frame(rownames(enetperf_b), 
                            enetperf_sm,
                            enetperf_b,
							enetperf_unb)
rownames(perfdat_sm_bs) <- NULL
colnames(perfdat_sm_bs) <- c("Metric","100/200", "200/200", "300/200", "400/200",
                       "Bootstraps","Unbalanced")
					   
perfdat_sm_bs

#kable(perfdat_sm_bs)
saveRDS(perfdat_sm_bs, "/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/perfdat_sm_bs.rds")

#mccf1 <- data.frame(Metric = c(rep("MCC",10), rep("F1",10)),
#                    Technique = rep(c("100/200", 
#                    "200/200", 
#                    "300/200", 
#                    "400/200", 
#                    "100/300",
#                    "200/300", 
#                    "300/300", 
#                    "400/300",
#                    "Bootstraps",
#					"Unbalanced"), 2),
#                    Value = c(as.numeric(perfdat[15,2:11]), as.numeric(perfdat[16,2:11])))

#mcc_f1_sm_bs <- ggplot(data=mccf1, aes(x=Technique, y=Value, fill=Metric)) +
#  geom_bar(stat="identity", position=position_dodge()) +
#  scale_fill_manual(values=c('black','lightgray')) +
#  xlab("Balancing Technique") + 
#  theme_minimal() + 
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#mcc_f1_sm_bs
#ggsave("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/mcc_f1_sm_bs.png")

#mcc_sm_bs <-ggplot(data=mccf1[1:10,], aes(x=Technique, y=Value, fill=Technique)) + 
#  xlab("Balancing Technique") + ylab("MCC") +
#  geom_bar(stat="identity") +
  #scale_fill_manual(values=gray(rev(seq(0,.8,.1))), guide=FALSE) +
#  theme_minimal() + 
#  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#  ggtitle("Model Performance for Different \n Class Balancing Techniques")
#mcc_sm_bs
#ggsave("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/mcc_sm_bs.png")


#f1_sm_bs <-ggplot(data=mccf1[10:18,], aes(x=Technique, y=Value, fill=Technique)) + 
#  xlab("Balancing Technique") + ylab("F1") +
#  geom_bar(stat="identity") +
#  scale_fill_manual(values=gray(rev(seq(0,.8,.1))), guide=FALSE) +
#  theme_minimal() + 
#  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#  ggtitle("Model Performance for Different \n Class Balancing Techniques")
#f1_sm_bs
#ggsave("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/f1_sm_bs.png")



# Comparing 100/200 SMOTE with Bootstrapped model

#roc curves
#roc.100200.bs <- ggplot() + 
#  geom_line(aes(fpr, tpr, colour=gray(.3)[1]), rocdat.bs) +
#  geom_line(aes(fpr, tpr, colour=gray(.7)[1]), onetwo) + 
#  geom_line(aes(fpr, tpr, colour="black"), roc.unb)
#  scale_colour_manual(name="Sampling \n Technique",
#    labels=c("Bootstrap","SMOTE: \n 100/200"),
#    values=c("black",gray(.7))) +
#  xlab("1-Specificity") + 
#  ylab("Sensitivity") + 
#  xlim(0, 1) +
#  ylim(0, 1) +
#  geom_abline(intercept=0, slope=1) +
#  theme_minimal() +
#  ggtitle("100 Bootstrap Samples vs 100/200 SMOTE")
#roc.100200.bs
#ggsave("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/roc.100200.bs.png")

  

#variable importance plots
#varimp.bs <- as.vector(rowMeans(enetlst_bs[[4]]))
#Labels <- rownames(enetlst_bs[[4]])
#Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
#varimp.bs.df <- data.frame(Feature=Labels,
#                                 Importance=varimp.bs)
#varimp.bs.df <- varimp.bs.df[order(varimp.bs.df$Importance),]
#varimp.bs.df <- varimp.bs.df[(dim(varimp.bs.df)[1]-19):dim(varimp.bs.df)[1],]
#varimp.bs.df$Feature <- factor(varimp.bs.df$Feature,
#                                     levels=varimp.bs.df$Feature)
#p.bs <- ggplot(varimp.bs.df, aes(x=Feature, y=Importance)) +
#  xlab("Predictors") +
#  ylab("Importance") +
#  #ggtitle("Importance Plot for Gradient Boosting Machine") +
#  geom_bar(stat="identity", 
#           width=.5, 
#           position="dodge",
#           fill="black") +
#  coord_flip() +
#  theme_minimal() +
#  ggtitle("Variable Importance Plot: \n 100 Bootstrap Samples")
#p.bs
#ggsave("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/p.bs.png")


#varimp.sm <- as.vector(enetlst_sm[[4]][,1])
#Labels <- names(enetlst_sm[[4]][,1])
#Labels[grep("Gm12878_", Labels)] <- gsub("Gm12878_","",Labels[grep("Gm12878_", Labels)])
#varimp.sm.df <- data.frame(Feature=Labels,
#                                 Importance=varimp.sm)
#varimp.sm.df <- varimp.sm.df[order(varimp.sm.df$Importance),]
#varimp.sm.df <- varimp.sm.df[(dim(varimp.sm.df)[1]-19):dim(varimp.sm.df)[1],]
#varimp.sm.df$Feature <- factor(varimp.sm.df$Feature,
#                                     levels=varimp.sm.df$Feature)
#p.sm <- ggplot(varimp.sm.df, aes(x=Feature, y=Importance)) +
#  xlab("Predictors") +
#  ylab("Importance") +
#  #ggtitle("Importance Plot for Gradient Boosting Machine") +
#  geom_bar(stat="identity", 
#           width=.5, 
#           position="dodge",
#           fill=gray(.7)) +
#  coord_flip() +
#  theme_minimal() +
#  ggtitle("Variable Importance Plot: \n 100/200 SMOTE")
#p.sm
#ggsave("/home/stilianoudakisc/TAD_data_analysis/evaluating_SMOTE/p.sm.png")


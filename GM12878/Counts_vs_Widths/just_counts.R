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

setwd("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/")


gm12878_5kb <- readRDS("gm12878_5kb.rds")

gm12878_5kb <- gm12878_5kb[,c(1,grep("count", names(gm12878_5kb)))]

dim(gm12878_5kb)

##Changing levels of response (y) to yes no
gm12878_5kb$y <- factor(gm12878_5kb$y)

levels(gm12878_5kb$y) <- c("No", "Yes")

### set number of bootstrap samples

bootsamps = 10

### set tuning parameters

fitControl <- trainControl(method = "cv",
                           number = 10,
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
						
###create a matrix of row ids that represent the zero class

##the number of rows will match the one class
##the number of columns match the number of bootstrap samples

sampids <- matrix(ncol=bootsamps, 
                  nrow=length(gm12878_5kb$y[which(gm12878_5kb$y=="Yes")]))


###filling in the sample ids matrix

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(gm12878_5kb$y=="No"),
                        length(which(gm12878_5kb$y=="Yes")),
                        replace = TRUE)
}

		
# Final Models:

# Random Forest

#set length of list objects that will be filled in with specificities
#and sensitivities and aucs and variable importance
rflst <- list(tpr <- matrix(nrow=ceiling((length(which(gm12878_5kb$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    fpr <- matrix(nrow=ceiling((length(which(gm12878_5kb$y=="Yes"))*2)*.3), 
                                  ncol=bootsamps),
                    varimp <- matrix(nrow=dim(gm12878_5kb)[2]-1,
                                     ncol=bootsamps) #,
					#coefs <- matrix(nrow=dim(gm12878_5kb)[2]-1,
					#				ncol=bootsamps)
					)
rownames(rflst[[3]]) <- colnames(gm12878_5kb)[-1]
#rownames(rflst[[4]]) <- colnames(gm12878_5kb)[-1]

rfperf <- matrix(nrow = 17, ncol=bootsamps)
rownames(rfperf) <- c("TN",
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
  data <- rbind.data.frame(gm12878_5kb[which(gm12878_5kb$y=="Yes"),],
                           gm12878_5kb[sampids[,i],])
  
  
  #determining the best number of variables randomly sampled as candidates at each split
  set.seed(5430)
  bestmtry <- tuneRF(data[,-1],data$y,
                   improve=.01,trace=0,plot=F) 
  bestmtry <- data.frame(bestmtry)
  bestmtry <- bestmtry[order(bestmtry$OOBError, decreasing = FALSE),]
  #bestmtry$mtry[1]

  #splitting the data
  inTrainingSet <- sample(length(data$y),floor(length(data$y)*.7))
  #inTrainingSet <- createDataPartition(data$y,p=.7,list=FALSE)
  train <- data[inTrainingSet,]
  test <- data[-inTrainingSet,]
  
  #determining best number of trees
  tunegrid <- expand.grid(.mtry=bestmtry$mtry[1])
  modellist <- list()
    for (ntree in c(50,200,500,1000)) {
      set.seed(333)
      fit <- train(y~., data=train, 
                   method="rf", 
                   metric="Accuracy",
                   tuneGrid=tunegrid,  
                   ntree=ntree)
      key <- toString(ntree)
      modellist[[key]] <- fit
    }
    # compare results
    results <- resamples(modellist)
    #summary(results)
    #dotplot(results)
    results <- data.frame(summary(results)[3]$statistics$Accuracy)
    results <- results[order(results$Mean, decreasing = TRUE),]

  set.seed(1006)
  rfModel <- train(y~., data=train, 
                    method="rf", 
                    metric="ROC", 
                    tuneGrid=tunegrid, 
                    trControl=fitControl, 
                    ntree=as.numeric(rownames(results)[1]))

  
  #Prediction vector for ROC and AUC
  pred.rfModel <- as.vector(predict(rfModel, 
                                    newdata=test, 
                                    type="prob")[,"Yes"])
  rflst[[1]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,1]
  rflst[[2]][,i] <- simple_roc(ifelse(test$y=="Yes",1,0),pred.rfModel)[,2]
  rflst[[3]][,i] <- varImp(rfModel)$importance[,1]
  #rflst[[4]][,i] <- coef(eNetModel$finalModel, eNetModel$bestTune$lambda)[-1,1]
  
  #Prediction vector for other performance metrics
  pred.rfModel2 <- predict(rfModel,
                                     newdata=test,
                                     type="raw")
  confMat <- confusionMatrix(data=pred.rfModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  rfperf[1,i] <- TN
  rfperf[2,i] <- FN
  rfperf[3,i] <- FP
  rfperf[4,i] <- TP
  rfperf[5,i] <- sum(confMat$table)
  rfperf[6,i] <- as.vector(confMat$byClass["Sensitivity"])
  rfperf[7,i] <- as.vector(confMat$byClass["Specificity"])
  rfperf[8,i] <- as.vector(confMat$overall["Kappa"])
  rfperf[9,i] <- as.vector(confMat$overall["Accuracy"])
  rfperf[10,i] <- TP/(TP+FP)
  rfperf[11,i] <- FP/(FP+TN)
  rfperf[12,i] <- FN/(FN+TN)
  rfperf[13,i] <- FN/(FN+TN)
  rfperf[14,i] <- TN/(TN+FN)
  rfperf[15,i] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  rfperf[16,i] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  rfperf[17,i] <- pROC::auc(pROC::roc(test$y, pred.rfModel))

  
}

saveRDS(rflst,"rflst_c.rds")
saveRDS(rfperf,"rfperf_c.rds")	


####################################################################################

varimp.c <- rowMeans(rflst_c[[3]])
varimp.c.df <- data.frame(Feature = rownames(rflst_c[[3]]), Importance =as.vector( varimp.c))
varimp.c.df <- varimp.c.df[order(varimp.c.df$Importance, decreasing=FALSE),]
varimp.c.df <- varimp.c.df[(dim(varimp.c.df)[1]-9):dim(varimp.c.df)[1],]
varimp.c.df$Feature <- factor(varimp.c.df$Feature,
                                     levels=varimp.c.df$Feature)
p.c <- ggplot(varimp.c.df, aes(x=Feature, y=Importance)) +
  xlab("Predictors") +
  ylab("Importance") +
  #ggtitle("Using Both Counts and Widths of Overlaps") +
  geom_bar(stat="identity", 
           width=.5, 
           position="dodge",
           fill=rainbow(3)[3]) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
		#axis.text.x = element_text(angle = 90, hjust = 1),
		#legend.title=element_text(size=20), 
		#legend.text=element_text(size=15),
		#legend.position = c(.7,.9)
  ggtitle("Elastic-Net")
p.c	
ggsave("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/p.c.png")

	

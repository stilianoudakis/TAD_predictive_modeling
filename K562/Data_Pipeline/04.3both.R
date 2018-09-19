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

setwd("/home/stilianoudakisc/TAD_data_analysis/K562/evaluating_variable_reduction/")

k562_10kb_f <- readRDS("k562_10kb_f.rds")


####set number of bootstrap samples

bootsamps = 10

#### function for roc curves

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


#### Establishing tuning/training parameters

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


####create a matrix of row ids that represent the zero class

sampids <- matrix(ncol=bootsamps, 
                  nrow=length(k562_10kb_f$y[which(k562_10kb_f$y=="Yes")]))
				  
####filling in the sample ids matrix

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(k562_10kb_f$y=="No"),
                        length(which(k562_10kb_f$y=="Yes")),
                        replace = TRUE)
}

## Stepwise (Both) Selection

namesmat <- matrix(NA,nrow=dim(k562_10kb_f)[2]-1, ncol=bootsamps)
rownames(namesmat) <- sort(names(k562_10kb_f)[-which(names(k562_10kb_f)=="y")])

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(k562_10kb_f[which(k562_10kb_f$y=="Yes"),],
                           k562_10kb_f[sampids[,i],])
  
  #null model
  glm.null <- glm(y ~ 1, data = data, family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = data, family = binomial)
  
  best.fit.both = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="both",
                      trace=0)
  
  preds <- names(best.fit.both$coefficients)[-1]
  preds <- gsub("`","",preds)
  
  namesmat[which(rownames(namesmat) %in% preds),i] <- sort(preds)
  
}

namesmat <- data.frame(namesmat)

namesmat$naperc <- (rowSums(!is.na(namesmat))/bootsamps)*100

both.preds <- rownames(namesmat)[which(namesmat$naperc >= 95)]

both.preds

k562_10kb_both <- k562_10kb_f[,which((names(k562_10kb_f) %in% both.preds) | names(k562_10kb_f)=="y")]

dim(k562_10kb_both)

saveRDS(k562_10kb_both,"k562_10kb_both.rds")

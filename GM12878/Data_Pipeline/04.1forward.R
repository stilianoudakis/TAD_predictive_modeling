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

gm12878_5kb_f <- readRDS("gm12878_5kb_f.rds")


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
                  nrow=length(gm12878_5kb_f$y[which(gm12878_5kb_f$y=="Yes")]))
				  
####filling in the sample ids matrix

set.seed(123)
for(j in 1:bootsamps){
  sampids[,j] <- sample(which(gm12878_5kb_f$y=="No"),
                        length(which(gm12878_5kb_f$y=="Yes")),
                        replace = FALSE)
}

## Forward

namesmat <- matrix(NA,nrow=dim(gm12878_5kb_f)[2]-1, ncol=bootsamps)
rownames(namesmat) <- sort(names(gm12878_5kb_f)[-which(names(gm12878_5kb_f)=="y")])

for(i in 1:bootsamps){
  set.seed(7215)
  #combining the two classes to create balanced data
  data <- rbind.data.frame(gm12878_5kb_f[which(gm12878_5kb_f$y=="Yes"),],
                           gm12878_5kb_f[sampids[,i],])
  
  #null model
  glm.null <- glm(y ~ 1, data = data, family = binomial)
  #full model
  glm.full <- glm(y ~ ., data = data, family = binomial)
  
  best.fit.fwd = step(glm.null,
                      scope=list(lower=formula(glm.null),
                                 upper=formula(glm.full)), 
                      direction="forward",
                      trace=0)
  
  preds <- names(best.fit.fwd$coefficients)[-1]
  preds <- gsub("`","",preds)
  
  namesmat[which(rownames(namesmat) %in% preds),i] <- sort(preds)
  
}

namesmat <- data.frame(namesmat)

namesmat$naperc <- (rowSums(!is.na(namesmat))/bootsamps)*100

fwd.preds <- rownames(namesmat)[which(namesmat$naperc >= 95)]

fwd.preds

gm12878_5kb_fwd <- gm12878_5kb_f[,which((names(gm12878_5kb_f) %in% fwd.preds) | names(gm12878_5kb_f)=="y")]

dim(gm12878_5kb_fwd)

saveRDS(gm12878_5kb_fwd,"gm12878_5kb_fwd.rds")

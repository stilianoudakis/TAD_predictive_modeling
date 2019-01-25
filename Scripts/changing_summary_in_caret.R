library(caret)
library(mlbench)
library(Hmisc)
library(randomForest)
library(MLmetrics)


#n <- 100
#p <- 40
#sigma <- 1
#set.seed(1)
#sim <- mlbench.friedman1(n, sd = sigma)
#colnames(sim$x) <- c(paste("real", 1:5, sep = ""),
#                     paste("bogus", 1:5, sep = ""))
#bogus <- matrix(rnorm(n * p), nrow = n)
#colnames(bogus) <- paste("bogus", 5+(1:ncol(bogus)), sep = "")
#x <- cbind(sim$x, bogus)
#y <- sim$y
#normalization <- preProcess(x)
#x <- predict(normalization, x)
#x <- as.data.frame(x)
#subsets <- c(1:5, 10, 15, 20, 25)

#set.seed(10)

#ctrl <- rfeControl(functions = lmFuncs,
#                   method = "repeatedcv",
#                   repeats = 5,
#                   verbose = FALSE)

#lmProfile <- rfe(x, y,
#                 sizes = subsets,
#                 rfeControl = ctrl)

#lmProfile

set.seed(346)
dat <- twoClassSim(200, noiseVars = 5)

f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
  c(F1 = f1_val)
}

rfRFE <-  list(summary = twoClassSummary,
               #summary = f1,
               fit = function(x, y, first, last, ...){
                 library(randomForest)
                 randomForest(x, y, importance = first, ntree=1000, mtry=sqrt(ncol(x)))
               },
               pred = function(object, x)  predict(object, x),
               rank = function(object, x, y) {
                 vimp <- varImp(object)
                 vimp <- vimp[order(vimp$Overall,decreasing = TRUE),,drop = FALSE]
                 vimp$var <- rownames(vimp)                  
                 vimp
               },
               selectSize = pickSizeBest,
               selectVar = pickVars)

#rfRFE$summary <- f1

ctrl <- rfeControl(functions = rfRFE,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)


#ctrl$functions <- rfRFE
ctrl$returnResamp <- "all"

subsets <- c(1:5, 10, 15, 20)

set.seed(10)

rfProfile <- rfe(dat[,-21], dat$Class, sizes = subsets, rfeControl = ctrl)
#rfProfile <- rfe(x, y, sizes = subsets, rfeControl = ctrl)
rfProfile

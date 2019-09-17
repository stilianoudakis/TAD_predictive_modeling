#function to predict cl, resolution, and chromosome specific tad boundaries using optimal parameters
predict_tad_bounds_func <- function(resolution, bounds.GR, datamatrix, chromosome, sampling="smote", metric="MCC", seed=123, crossvalidation=TRUE, number=10, featureSelection=FALSE){
  resolution=as.integer(resolution)
  
  bounds.GR <- flank(bounds.GR, width=resolution, both=TRUE)
  
  y <- countOverlaps(GRanges(seqnames = tolower(chromosome),
                             IRanges(start = as.numeric(rownames(datamatrix)),
                                     #end = as.numeric(rownames(datamatrix)),
                                     width = resolution)), bounds.GR)
  y <- ifelse(y>=1,1,0)
  
  full_data <- cbind.data.frame(y, datamatrix)
  ##Changing levels of response (y) to yes no
  full_data$y <- factor(full_data$y)
  levels(full_data$y) <- c("No", "Yes")
  
  #splitting into training and testing
  set.seed(seed)
  inTrainingSet <- sample(length(full_data$y),floor(length(full_data$y)*.7))
  train <- full_data[inTrainingSet,]
  test <- full_data[-inTrainingSet,]
  
  if(sampling=="ros"){
    #assign sample indeces
    set.seed(seed)
    sampids.train <- sample(x = which(train$y=="Yes"),
                            size = length(which(train$y=="No")),
                            replace = TRUE)
    
    train <- rbind.data.frame(train[which(train$y=="No"),],
                              train[sampids.train,])
    
    #Randomly shuffle the data
    set.seed(seed)
    train <- train[sample(nrow(train)),]
    
  }else if(sampling=="rus"){
    set.seed(123)
    sampids.train <- sample(x = which(train$y=="No"),
                            size = length(which(train$y=="Yes")),
                            replace = FALSE)
    
    train <- rbind.data.frame(train[which(train$y=="Yes"),],
                              train[sampids.train,])
    
    #Randomly shuffle the data
    set.seed(321)
    train <- train[sample(nrow(train)),]
    
  }else if(sampling=="smote"){
    set.seed(seed)
    train <- SMOTE(y ~ ., 
                   data=train, 
                   perc.over = 100, 
                   perc.under = 200)
    
    #Randomly shuffle the data
    set.seed(seed)
    train <- train[sample(nrow(train)),]
  }else{train=train}
  
  
  ##defining one summary function for roc,auprc,f1,and mcc
  allSummary <- function (data, lev = NULL, model = NULL) {
    lvls <- levels(data$obs)
    
    #mcc
    mcc <- ModelMetrics::mcc(ifelse(data$obs == lev[2], 0, 1), data[, lvls[1]], cutoff = .5)
    
    #roc
    b1 <- twoClassSummary(data, lev, model)
    
    #auprc & f1
    c1 <- prSummary(data, lev, model)
    
    out <- c(mcc, b1, c1)
    names(out)[1] <- c("MCC")
    out
  }
  
  set.seed(seed)
  seeds <- vector(mode = "list", length = (number+1))
  for(i in 1:number) seeds[[i]]<- sample.int(n=1000, 10)
  #for the last model
  seeds[[11]]<-sample.int(1000, 1)
  
  if(featureSelection==FALSE){
    ##setting contols for model
    fitControl <- trainControl(seeds = seeds,
                               method = if(crossvalidation==TRUE){"cv"}else{"none"},
                               number = if(crossvalidation==TRUE){number}else{NULL},
                               ## Estimate class probabilities
                               classProbs = TRUE,
                               ## Evaluate performance using 
                               ## the following function
                               summaryFunction = allSummary)
    
    #performing model
    tadModel <- train(y~., data=train, 
                      method="rf", 
                      metric=metric, 
                      tuneLength=10,
                      trControl=fitControl)
  }else{
    rfectrl <- rfFuncs
    rfectrl$summary <- allSummary
    control <- rfeControl(functions=rfectrl, 
                          method=if(crossvalidation==TRUE){"cv"}else{"none"}, 
                          number=if(crossvalidation==TRUE){number}else{NULL},  
                          verbose=TRUE)
    control$returnResamp <- "final"
    
    n=dim(train)[2]-1
    z <- numeric()
    x=0
    i=1
    while(x < n){
      x = 2^(i);i=i+1;z <- c(z,x)
    }
    z[length(z)] <- n
    
    set.seed(seed)
    tadModel <- rfe(train[,-1], 
                    train[,1],
                    metric = metric,
                    sizes=z,
                    rfeControl=control)
  }
  
  return(tadModel)
}

#function to predict cl, resolution, and chromosome specific tad boundaries using optimal parameters
predict_tad_bounds_func <- function(bounds.GR, datamatrix, chromosome, model, sampling, metric, seed, crossvalidation, number){
  
  y <- countOverlaps(GRanges(seqnames = tolower(chromosome),
                             IRanges(start = as.numeric(rownames(datamatrix)),
                                     end = as.numeric(rownames(datamatrix)))), bounds.GR)
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
  
  #reduce predictor space
  if(predictor=="distance"){reduce_predictor = "_dist"}else if(predictor=="oc"){reduce_predictor = "_count"}else{reduce_predictor = "_perc"}
  
  train <- train[,c(1,grep(reduce_predictor, names(train)))]
  test <- test[,c(1,grep(reduce_predictor, names(test)))]
  
  # set number of iterations
  samps = 50
  
  if(sampling=="ros"){
    #assign sample indeces
    sampids <- matrix(ncol=samps, 
                      nrow=length(train$y[which(train$y=="No")]))
    
    #filling in the sample ids matrix
    set.seed(123)
    for(j in 1:samps){
      sampids[,j] <- sample(x = which(train$y=="Yes"),
                            size = length(which(train$y=="No")),
                            replace = TRUE)
    }
    train <- rbind.data.frame(train[which(train$y=="No"),],
                              train[sampids[,1],])
    #Randomly shuffle the data
    set.seed(321)
    train <- train[sample(nrow(train)),]
  }else if(sampling=="rus"){
    #assign sample indeces
    sampids <- matrix(ncol=samps, 
                      nrow=length(train$y[which(train$y=="Yes")]))
    
    #filling in the sample ids matrix
    set.seed(123)
    for(j in 1:samps){
      sampids[,j] <- sample(x = which(train$y=="No"),
                            size = length(which(train$y=="Yes")),
                            replace = FALSE)
    }
    train <- rbind.data.frame(train[which(train$y=="Yes"),],
                              train[sampids[,1],])
    
    #Randomly shuffle the data
    set.seed(321)
    train <- train[sample(nrow(train)),]
  }else if(sampling=="smote"){
    set.seed(123)
    train <- SMOTE(y ~ ., 
                   data=train, 
                   perc.over = 100, 
                   perc.under = 200)
    
    #Randomly shuffle the data
    set.seed(123)
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
  
  ##setting contols for elastic net
  fitControl <- trainControl(#seeds = seeds,
    method = crossvalidation,
    number = number,
    ## Estimate class probabilities
    classProbs = TRUE,
    ## Evaluate performance using 
    ## the following function
    summaryFunction = allSummary)
  
  #performing model
  tadModel <- train(y~., data=train, 
                   method=model, 
                   metric=metric, 
                   tuneLength=10,
                   trControl=fitControl)
  return(tadModel)
}
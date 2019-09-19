#function to predict cl, resolution, and chromosome specific tad boundaries using optimal parameters
predict_tad_bounds_func_old <- function(resolution, bounds.GR, annotationListGR, chromosome, predictortype=predictor, sampling="smote", metric="MCC", seed=123, crossvalidation=TRUE, number=10){
  resolution=as.integer(resolution)
  
  seqData <- c(seqData[1]-1, seqData)
  
  start=0
  end=seqData[length(seqData)] - (seqData[length(seqData)] %% resolution)
  rows = seqData[seqData %in% seq(start, end, resolution)]
  
  data_mat <- matrix(nrow=length(rows),
                     ncol=length(annotationListGR))
  rownames(data_mat) <- rows
  
  dat_mat_gr <- GRanges(seqnames = tolower(chromosome),
                        IRanges(start = as.numeric(rownames(data_mat)),
                                #end = as.numeric(rownames(data_mat)),
                                width=resolution))
  
  dat_mat_gr_dist <- GRanges(seqnames = tolower(chromosome),
                             IRanges(start = as.numeric(rownames(data_mat)),
                                     #end = as.numeric(rownames(data_mat)),
                                     width=1))
  
  if(predictortype=="distance"){
    for(i in 1:length(annotationListGR)){
      d <- distance_func(dat_mat_gr_dist, annotationListGR[[i]])
      data_mat[,i] <- d
    }
    data_mat <- apply(data_mat,2,function(x){log(x + 1, base=2)})
    colnames(data_mat) <- names(annotationListGR)
  }else if(predictortype=="binary"){
    for(i in 1:length(annotationListGR)){
      cb <- binary_func(dat_mat_gr, annotationListGR[[i]])
      data_mat[,i] <- cb
    }
    colnames(data_mat) <- names(annotationListGR)
  }else if(predictortype=="oc"){
    for(i in 1:length(annotationListGR)){
      co <- count_func(dat_mat_gr, annotationListGR[[i]])
      data_mat[,i] <- co
    }
    colnames(data_mat) <- names(annotationListGR)
  }else{
    for(i in 1:length(annotationListGR)){
      cp <- percent_func(dat_mat_gr, annotationListGR[[i]])
      data_mat[,i] <- cp
    }
    colnames(data_mat) <- names(annotationListGR)
  }
  
  bounds.GR <- flank(bounds.GR, width=resolution, both=TRUE)
  
  y <- countOverlaps(GRanges(seqnames = tolower(chromosome),
                             IRanges(start = as.numeric(rownames(data_mat)),
                                     #end = as.numeric(rownames(data_mat)),
                                     width = resolution)), bounds.GR)
  y <- ifelse(y>=1,1,0)
  
  full_data <- cbind.data.frame(y, data_mat)
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
  
  rfperf <- data.frame(Metric = c("TN",
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
                                  "AUC",
                                  "Youden",
                                  "AUPRC"), 
                       Performance=NA)
  
  pred.tadModel <- as.vector(predict(tadModel, 
                                     newdata=test[,-1], 
                                     type="prob")[,"Yes"])
  
  fg <- pred.tadModel[test$y == "Yes"]
  bg <- pred.tadModel[test$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  pred.tadModel2 <- predict(tadModel,
                            newdata=test,
                            type="raw")
  
  
  confMat <- confusionMatrix(data=pred.tadModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  rfperf[1,2] <- TN
  rfperf[2,2] <- FN
  rfperf[3,2] <- FP
  rfperf[4,2] <- TP
  rfperf[5,2] <- sum(confMat$table)
  rfperf[6,2] <- as.vector(confMat$byClass["Sensitivity"])
  rfperf[7,2] <- as.vector(confMat$byClass["Specificity"])
  rfperf[8,2] <- as.vector(confMat$overall["Kappa"])
  rfperf[9,2] <- as.vector(confMat$overall["Accuracy"])
  rfperf[10,2] <- TP/(TP+FP)
  rfperf[11,2] <- FP/(FP+TN)
  rfperf[12,2] <- FN/(FN+TN)
  rfperf[13,2] <- FN/(FN+TN)
  rfperf[14,2] <- TN/(TN+FN)
  rfperf[15,2] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  rfperf[16,2] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  rfperf[17,2] <- pROC::auc(pROC::roc(test$y, pred.tadModel))
  rfperf[18,2] <- (TP/(TP + FN)) + (TN/(TN + FP)) - 1
  rfperf[19,2] <- pr$auc.integral
  
  
  tadModelResults <- list(tadModel, rfperf)
  names(tadModelResults) <- c("RFModel", "RFPerformances")
  
  return(tadModelResults)
}

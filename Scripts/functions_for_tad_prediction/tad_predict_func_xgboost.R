
tad_predict_func_noenet_xgboost <- function(CL, resolution, sampling, predictor, chromosome){
  
  ##gbm
  ###model results
  save_data_directory_gbmmodel <- paste0("/home/stilianoudakisc/TAD_data_analysis/",
                                         CL,
                                         "/",
                                         resolution,
                                         "/results_by_chr/",
                                         chromosome,
                                         "/",
                                         sampling,
                                         "/",
                                         predictor,
                                         "/enet/gbmModel_noenet.rds")
  ###variable importances
  save_data_directory_gbmimpvars <- paste0("/home/stilianoudakisc/TAD_data_analysis/",
                                           CL,
                                           "/",
                                           resolution,
                                           "/results_by_chr/",
                                           chromosome,
                                           "/",
                                           sampling,
                                           "/",
                                           predictor,
                                           "/enet/gbmimpvars_noenet.rds")
  ###performances							  
  save_data_directory_gbmperf <- paste0("/home/stilianoudakisc/TAD_data_analysis/",
                                        CL,
                                        "/",
                                        resolution,
                                        "/results_by_chr/",
                                        chromosome,
                                        "/",
                                        sampling,
                                        "/",
                                        predictor,
                                        "/enet/gbmperf_noenet.rds")								
  
  #reading in whole genome data
  ##setting directory specific to cell line and resolution 
  if(CL=="GM12878"){dataset = "gm12878"}else{dataset = "K562"}
  wgd_directory <- paste0("/home/stilianoudakisc/TAD_data_analysis/",
                          CL,
                          "/",
                          resolution,
                          "/training_testing/",
                          dataset,
                          "_",
                          resolution,
                          ".rds")
  whole_genome_dat <- readRDS(wgd_directory)
  
  #reading in bin data on the whole genome that contains chromosome information
  ##setting directory specific to cell line and resolution 
  bd_directory <- paste0("/home/stilianoudakisc/TAD_data_analysis/",
                         CL,
                         "/",
                         resolution,
                         "/training_testing/binslist",
                         gsub("kb","",resolution),
                         "_center.rds")
  binslist_dat <- readRDS(bd_directory)
  
  #removing duplicate annotations
  if(CL=="GM12878"){
    varstorm = c(names(whole_genome_dat)[grep("Ctcf-", names(whole_genome_dat))][-c(1:4)],
                 names(whole_genome_dat)[grep("Ebf1sc137065-", names(whole_genome_dat))][-c(1:4)],
                 names(whole_genome_dat)[grep("P300-", names(whole_genome_dat))][-c(1:4)],
                 names(whole_genome_dat)[grep("Pol2-", names(whole_genome_dat))][-c(1:4)],
                 names(whole_genome_dat)[grep("Rad21-", names(whole_genome_dat))][-c(1:4)])} else {
                   varstorm = c(names(whole_genome_dat)[grep("Atf3-", names(whole_genome_dat))][-c(1:4)],
                                names(whole_genome_dat)[grep("Cmyc-", names(whole_genome_dat))][-c(1:4)],
                                names(whole_genome_dat)[grep("Ctcf-", names(whole_genome_dat))][-c(1:4)],
                                names(whole_genome_dat)[grep("E2f6-", names(whole_genome_dat))][-c(1:4)],
                                names(whole_genome_dat)[grep("Max-", names(whole_genome_dat))][-c(1:4)],
                                names(whole_genome_dat)[grep("P300-", names(whole_genome_dat))][-c(1:4)],
                                names(whole_genome_dat)[grep("Pol2-", names(whole_genome_dat))][-c(1:4)],
                                names(whole_genome_dat)[grep("Rad21-", names(whole_genome_dat))][-c(1:4)],
                                names(whole_genome_dat)[grep("Yy1-", names(whole_genome_dat))][-c(1:4)])
                 }
  
  whole_genome_dat <- whole_genome_dat[,-which(names(whole_genome_dat) %in% varstorm)]
  
  #reducing the whole genome data to only specific chromosome
  reduce_chr <- tolower(chromosome)
  chr_specific_dat <- whole_genome_dat[which(as.character(seqnames(binslist_dat))==reduce_chr),]
  
  #splitting into training and testing
  set.seed(123)
  inTrainingSet <- sample(length(chr_specific_dat$y),floor(length(chr_specific_dat$y)*.7))
  train <- chr_specific_dat[inTrainingSet,]
  test <- chr_specific_dat[-inTrainingSet,]
  
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
  
  #defining one summary function for roc,auprc,f1,and mcc
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
  
  #setting contols xgboost
  fitControl <- trainControl(#seeds = seeds,
    method = "cv",
    number = 5,
    ## Estimate class probabilities
    classProbs = TRUE,
    ## Evaluate performance using 
    ## the following function
    summaryFunction = allSummary)
  
  #Setting tunning parameters
  
  ##Step 1: Number of Iterations and the Learning Rate
  tune_grid <- expand.grid(
	nrounds = c(50, 100, 500, 1000), #Number of trees, default: 100
	eta = c(0.025, 0.05, 0.1, 0.3), #Learning rate, default: 0.3
	max_depth = c(2, 4, 6, 8), #Maximum tree depth, default: 6
	gamma = 0, #Used for tuning of Regularization, default: 0
	colsample_bytree = 1, #Column sampling, default: 1
	min_child_weight = 1, #Minimum leaf weight, default: 1
	subsample = 1 #Row sampling, default: 1
  )
  
  ###set seeds for each cross folds for Step 1
  set.seed(123)
  seeds <- vector(mode = "list", length = 6) #6 refers (n_repeats*nresampling)+1
  for(i in 1:5) seeds[[i]]<- sample.int(n=1000, 64) #5 refers to (n_repeats*nresampling); 54 refers to the number of tuning parameter combinations
  #for the last model
  seeds[[6]]<-sample.int(1000, 1) #is used to set the seed for the last (optimum) model fit
  
 
  ###performing xgboost
  gbmModel_tune_1 <- train(y~., data=train, 
                    method="xgbTree", 
                    metric="MCC", 
                    tuneGrid=tune_grid, 
                    trControl=fitControl)
					
  ##Step 2: Fixing learning rate, maximum depth, and nrounds; tuning min_child_weight
  tune_grid2 <- expand.grid(
	nrounds = gbmModel_tune_1$bestTune$nrounds, #Number of trees, default: 100
	eta = gbmModel_tune_1$bestTune$eta, #Learning rate, default: 0.3
	max_depth = gbmModel_tune_1$bestTune$max_depth, #Maximum tree depth, default: 6
	gamma = 0, #Used for tuning of Regularization, default: 0
	colsample_bytree = 1, #Column sampling, default: 1
	min_child_weight = c(1, 2, 3), #Minimum leaf weight, default: 1
	subsample = 1 #Row sampling, default: 1
  )
  
  ###set seeds for each cross folds for Step 1
  set.seed(123)
  seeds <- vector(mode = "list", length = 6) #6 refers (n_repeats*nresampling)+1
  for(i in 1:5) seeds[[i]]<- sample.int(n=1000, 3) #5 refers to (n_repeats*nresampling); 54 refers to the number of tuning parameter combinations
  #for the last model
  seeds[[6]]<-sample.int(1000, 1) #is used to set the seed for the last (optimum) model fit
  
  ###performing xgboost
  gbmModel_tune_2 <- train(y~., data=train, 
                    method="xgbTree", 
                    metric="MCC", 
                    tuneGrid=tune_grid2, 
                    trControl=fitControl)
  
  ##Step 3: Tuning Column and Row Sampling
  tune_grid3 <- expand.grid(
	nrounds = gbmModel_tune_2$bestTune$nrounds, #Number of trees, default: 100
	eta = gbmModel_tune_2$bestTune$eta, #Learning rate, default: 0.3
	max_depth = gbmModel_tune_2$bestTune$max_depth, #Maximum tree depth, default: 6
	gamma = 0, #Used for tuning of Regularization, default: 0
	colsample_bytree = c(0.4, 0.6, 0.8, 1.0), #Column sampling, default: 1
	min_child_weight = gbmModel_tune_2$bestTune$min_child_weight, #Minimum leaf weight, default: 1
	subsample = c(0.5, 0.75, 1.0) #Row sampling, default: 1
  )
  
  ###set seeds for each cross folds for Step 1
  set.seed(123)
  seeds <- vector(mode = "list", length = 6) #6 refers (n_repeats*nresampling)+1
  for(i in 1:5) seeds[[i]]<- sample.int(n=1000, 12) #5 refers to (n_repeats*nresampling); 54 refers to the number of tuning parameter combinations
  #for the last model
  seeds[[6]]<-sample.int(1000, 1) #is used to set the seed for the last (optimum) model fit
  
  ###performing xgboost
  gbmModel_tune_3 <- train(y~., data=train, 
                    method="xgbTree", 
                    metric="MCC", 
                    tuneGrid=tune_grid3, 
                    trControl=fitControl)
  
  ##Step 4: Tunning gamma
  #tune_grid4 <- expand.grid(
  #	nrounds = gbmModel_tune_3$bestTune$nrounds, #Number of trees, default: 100
  #	eta = gbmModel_tune_3$bestTune$eta, #Learning rate, default: 0.3
  #	max_depth = gbmModel_tune_3$bestTune$max_depth, #Maximum tree depth, default: 6
  # gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0), #Used for tuning of Regularization, default: 0
  #	colsample_bytree = gbmModel_tune_3$bestTune$colsample_bytree, #Column sampling, default: 1
  #	min_child_weight = gbmModel_tune_3$bestTune$min_child_weight, #Minimum leaf weight, default: 1
  #	subsample = gbmModel_tune_3$bestTune$subsample #Row sampling, default: 1
  #)
  
  ###set seeds for each cross folds for Step 1
  #set.seed(123)
  #seeds <- vector(mode = "list", length = 6) #6 refers (n_repeats*nresampling)+1
  #for(i in 1:5) seeds[[i]]<- sample.int(n=1000, 7) #5 refers to (n_repeats*nresampling); 54 refers to the number of tuning parameter combinations
  #for the last model
  #seeds[[6]]<-sample.int(1000, 1) #is used to set the seed for the last (optimum) model fit
  
  ###performing xgboost
  #gbmModel_tune_4 <- train(y~., data=train, 
  #                  method="xgbTree", 
  #                  metric="MCC", 
  #                  tuneGrid=tune_grid4, 
  #                  trControl=fitControl)
  
  ##Fitting the final Model
  final_grid <- expand.grid(
	nrounds = gbmModel_tune_3$bestTune$nrounds,
	eta = gbmModel_tune_3$bestTune$eta,
	max_depth = gbmModel_tune_3$bestTune$max_depth,
	gamma = 0,
	colsample_bytree = gbmModel_tune_3$bestTune$colsample_bytree,
	min_child_weight = gbmModel_tune_3$bestTune$min_child_weight,
	subsample = gbmModel_tune_3$bestTune$subsample
  )
  
  ###performing xgboost
  gbmModel_tune_final <- train(y~., data=train, 
                    method="xgbTree", 
                    metric="MCC", 
                    tuneGrid=final_grid, 
                    trControl=fitControl)
  
  
  saveRDS(gbmModel_tune_final, save_data_directory_gbmmodel)
  
  #calculating variable importance
  gbmimpvars <- data.frame(Feature=rownames(varImp(gbmModel_tune_final)$importance), Importance=varImp(gbmModel_tune_final)$importance[,1])
  gbmimpvars <- gbmimpvars[order(gbmimpvars$Importance, decreasing=TRUE),]
  
  saveRDS(gbmimpvars, save_data_directory_gbmimpvars)
  
  #create matrix of performance metrics
  gbmperf <- data.frame(Metric = c("TN",
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
  
  pred.gbmModel <- as.vector(predict(gbmModel_tune_final, 
                                     newdata=test, 
                                     type="prob")[,"Yes"])
  
  fg <- pred.gbmModel[test$y == "Yes"]
  bg <- pred.gbmModel[test$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  pred.gbmModel2 <- predict(gbmModel_tune_final,
                            newdata=test,
                            type="raw")
  
  
  confMat <- confusionMatrix(data=pred.gbmModel2, test$y, positive="Yes")
  TN = as.numeric(confMat$table[1,1])
  FN = as.numeric(confMat$table[1,2])
  FP = as.numeric(confMat$table[2,1])
  TP = as.numeric(confMat$table[2,2])
  gbmperf[1,2] <- TN
  gbmperf[2,2] <- FN
  gbmperf[3,2] <- FP
  gbmperf[4,2] <- TP
  gbmperf[5,2] <- sum(confMat$table)
  gbmperf[6,2] <- as.vector(confMat$byClass["Sensitivity"])
  gbmperf[7,2] <- as.vector(confMat$byClass["Specificity"])
  gbmperf[8,2] <- as.vector(confMat$overall["Kappa"])
  gbmperf[9,2] <- as.vector(confMat$overall["Accuracy"])
  gbmperf[10,2] <- TP/(TP+FP)
  gbmperf[11,2] <- FP/(FP+TN)
  gbmperf[12,2] <- FN/(FN+TN)
  gbmperf[13,2] <- FN/(FN+TN)
  gbmperf[14,2] <- TN/(TN+FN)
  gbmperf[15,2] <- (TP*TN - FP*FN)/( sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ) )
  gbmperf[16,2] <- (2*(TP/(TP+FN))*(TP/(TP+FP)))/((TP/(TP+FN)) + ((TP/(TP+FP))))
  gbmperf[17,2] <- pROC::auc(pROC::roc(test$y, pred.gbmModel))			
  gbmperf[18,2] <- (TP/(TP + FN)) + (TN/(TN + FP)) - 1
  gbmperf[19,2] <- pr$auc.integral
  
  saveRDS(gbmperf, save_data_directory_gbmperf)
  
}



annot_specific_predict_func <- function(CL, resolution, sampling="smote", predictor="distance", annotation, chromosome){
	##annotation-specific random forest
	###model results
	save_data_directory_rfmodel <- paste0("/home/stilianoudakisc/TAD_data_analysis/",
								  CL,
								  "/",
								  resolution,
								  "/results_by_chr/",
								  chromosome,
								  "/",
								  sampling,
								  "/",
								  predictor,
								  "/enet/rfModel_",
								  annotation,
								  ".rds")
	###variable importances
	save_data_directory_rfimpvars <- paste0("/home/stilianoudakisc/TAD_data_analysis/",
								  CL,
								  "/",
								  resolution,
								  "/results_by_chr/",
								  chromosome,
								  "/",
								  sampling,
								  "/",
								  predictor,
								  "/enet/rfimpvars_",
								  annotation,
								  ".rds")
	###performances							  
	save_data_directory_rfperf <- paste0("/home/stilianoudakisc/TAD_data_analysis/",
								  CL,
								  "/",
								  resolution,
								  "/results_by_chr/",
								  chromosome,
								  "/",
								  sampling,
								  "/",
								  predictor,
								  "/enet/rfperf_",
								  annotation,
								  ".rds")	

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
	if(annotation=="tfbs"){
		vars <- list.files(paste0("/home/stilianoudakisc/TAD_data_analysis/annotations/tfbs/cell_type_specific/",
						   tolower(CL)))
		vars <- gsub(".bed", "", vars)
		vars <- paste0(vars,"_dist")
		#reduce training features
		train <- train[,c(1,which(names(train) %in% vars))]
		}else if(annotation=="chromstates"){
			vars <- list.files("/home/stilianoudakisc/TAD_data_analysis/annotations/BroadHMM")
				if(CL=="GM12878"){
					vars <- vars[grep("Gm12878", vars)]
					vars <- gsub(".bed", "", vars)
					vars <- paste0(vars,"_dist")
					#reduce training features
					train <- train[,c(1,which(names(train) %in% vars))]
					}else{
						vars <- vars[grep("K562", vars)]
						vars <- gsub(".bed", "", vars)
						vars <- paste0(vars,"_dist")
						#reduce training features
						train <- train[,c(1,which(names(train) %in% vars))]
						}
			}else{
				vars <- list.files("/home/stilianoudakisc/TAD_data_analysis/annotations/HistoneModifications")
				if(CL=="GM12878"){
					vars <- vars[grep("Gm12878", vars)]
					vars <- gsub(".bed", "", vars)
					vars <- paste0(vars,"_dist")
					#reduce training features
					train <- train[,c(1,which(names(train) %in% vars))]
					}else{
						vars <- vars[grep("K562", vars)]
						vars <- gsub(".bed", "", vars)
						vars <- paste0(vars,"_dist")
						#reduce training features
						train <- train[,c(1,which(names(train) %in% vars))]
						}
				}
	
	#resampling
	set.seed(123)
	train <- SMOTE(y ~ ., 
					data=train, 
                    perc.over = 100, 
                    perc.under = 200)
				
	#Randomly shuffle the data
	set.seed(123)
	train <- train[sample(nrow(train)),]
	
	##set tuning parameters
	##set seeds for each cross folds
	set.seed(123)
	#length is = (n_repeats*nresampling)+1
	seeds <- vector(mode = "list", length = 11)
	#(1 is the number of tuning parameter lambda
	for(i in 1:10) seeds[[i]]<- sample.int(n=1000, 100)
	#for the last model
	seeds[[11]]<-sample.int(1000, 1)

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
                           method = "cv",
                           number = 10,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = allSummary)
	
	#set number of predictors to consider at each node
	tunegrid <- expand.grid(mtry=ceiling(sqrt(dim(train)[2] - 1)))
	
	#performing random forest
	rfModel <- train(y~., data=train, 
                 method="rf", 
                 metric="MCC", 
                 #tuneGrid=tunegrid, 
				 tuneLength=10,
                 trControl=fitControl, 
                 ntree=500) #as.numeric(rownames(results)[1]))
	
	saveRDS(rfModel, save_data_directory_rfmodel)
	
	#calculating variable importance
	rfimpvars <- data.frame(Feature=rownames(varImp(rfModel)$importance), Importance=varImp(rfModel)$importance[,1])
	rfimpvars <- rfimpvars[order(rfimpvars$Importance, decreasing=TRUE),]
	
	saveRDS(rfimpvars, save_data_directory_rfimpvars)
	
	#create matrix of performance metrics
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
	
	pred.rfModel <- as.vector(predict(rfModel, 
                                  newdata=test, 
                                  type="prob")[,"Yes"])

	fg <- pred.rfModel[test$y == "Yes"]
	bg <- pred.rfModel[test$y == "No"]
	pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

	pred.rfModel2 <- predict(rfModel,
                         newdata=test,
                         type="raw")


	confMat <- confusionMatrix(data=pred.rfModel2, test$y, positive="Yes")
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
	rfperf[17,2] <- pROC::auc(pROC::roc(test$y, pred.rfModel))			
	rfperf[18,2] <- (TP/(TP + FN)) + (TN/(TN + FP)) - 1
	rfperf[19,2] <- pr$auc.integral
	
	saveRDS(rfperf, save_data_directory_rfperf)
}

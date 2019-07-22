library(GenomicRanges)
library(ggplot2)
library(caret)
library(GenomeInfoDbData)
library(GenomeInfoDb)

#function to create grangeslist object from functional genomic annotation .bed files
annots_to_granges_func <- function(filepath){
  annots_gr_list <- GRangesList() 
  #if(grepl(paste0(names, collapse = "|"), list.files(filepath)))
  for(i in 1:length(list.files(filepath))){
    dat <- read.table(paste0(filepath,"/",list.files(filepath)[i]),
                      header = FALSE, 
                      sep="\t")
    
    annots_gr_list[[i]] <- GRanges(seqnames=dat$V1,
                                   IRanges(start=dat$V2, 
                                           end=dat$V3))
  }
  names(annots_gr_list) <- gsub(".bed", "", list.files(filepath))
  return(annots_gr_list)
}

###############################################################################3333

#Functions for calculating different types of predictors/overlaps
####calculating binary overlaps
binary_func <- function(binned_data_gr, annot_data_gr){
  #Finding the total number of overlaps between genomic bins and the specific genomic annotation
  count_binary <- countOverlaps(binned_data_gr, annot_data_gr)
  
  #Binarizing the overlap (1 or 0)
  count_binary <- ifelse(count_binary>0,1,0)
  
  return(count_binary)
}

####calculating count overlaps
count_func <- function(binned_data_gr, annot_data_gr){
  #Finding the total number of overlaps between genomic bins and the specific genomic annotation
  count_total <- countOverlaps(binned_data_gr, annot_data_gr)
  
  return(count_total)
}

####calculating percent overlaps
percent_func <- function(binned_data_gr, annot_data_gr){
  count_percent <- numeric(length(binned_data_gr))
  
  #Finding the total number of overlaps between genomic bins and the specific genomic annotation
  c <- countOverlaps(binned_data_gr, annot_data_gr)
  
  #places where c=0 denotes no overlap 
  #places where c>0 denotes some type of overlap, could be partial or within
  
  #for c=0 assign percent overlap as 0
  count_percent[which(c==0)] <- 0
  
  #for c=1:
  count_percent[which(c==1)] <- width(pintersect(findOverlapPairs(annot_data_gr,binned_data_gr[which(c==1)])))/(width(binned_data_gr[1]))
  
  #for c>1:
  #iterate through all bins with multiple overlaps with the annotation of interest
  #find how many and which annotations overlap within each iterate
  #calculate the width of each overlap
  #sum the total widths and divide by bin width
  #bins with some type of overlap
  mo <- which(c>1)
  count_percent[mo] <- unlist(
    lapply(
      lapply(
        lapply(
          mo, function(x){width(pintersect(findOverlapPairs(annot_data_gr,binned_data_gr[x])))}),
        sum),
      function(x){x/width(binned_data_gr[1])}))
  
  return(count_percent)
}

####calculating distance
distance_func <- function(binned_data_center_gr, annot_data_center_gr){
  #distance from center of genomic bin to nearest genomic region of interest
  dist <- mcols(distanceToNearest(binned_data_center_gr, annot_data_center_gr))$distance
  
  return(dist)
}

########################################################################################

#function to create data matrix for training predictive model
annots_to_datamatrix_func <- function(resolution, predictortype="distance", annotationListGR, chromosome){
  #determining dimensions of chromosome specific data matrix
  genome <- getBSgenome("hg19")
  seqlength <- genome@seqinfo@seqlengths[which(genome@seqinfo@seqnames==tolower(chromosome))]
  ##start = first chromosome specific coordinate
  ##end = last chromosome specific coordinate that is a factor of resolution
  start = 0
  end = seqlength - (seqlength %% resolution)
  data_mat <- matrix(nrow=length(seq.int(start, 
                                  end,
                                  by=resolution)), 
                     ncol=length(annotationListGR))
  rownames(data_mat) <- seq.int(start, 
                                end,
                                by=resolution)
  
  dat_mat_gr <- GRanges(seqnames = tolower(chromosome),
                        IRanges(start = as.numeric(rownames(data_mat)),
                                end = as.numeric(rownames(data_mat))))
  
  if(predictortype=="distance"){
    for(i in 1:length(annotationListGR)){
      d <- distance_func(dat_mat_gr, annotationListGR[[i]])
      data_mat[,i] <- d
    }
    data_mat <- apply(data_mat,2,function(x){log(x + 1, base=2)})
    colnames(data_mat) <- names(annotationListGR)
  }
  else if(predictortype=="binary"){
    for(i in 1:length(annotationListGR)){
      cb <- binary_func(dat_mat_gr, annotationListGR[[i]])
      data_mat[,i] <- cb
    }
    colnames(data_mat) <- names(annotationListGR)
  }
  else if(predictortype=="count"){
    for(i in 1:length(annotationListGR)){
      co <- count_func(dat_mat_gr, annotationListGR[[i]])
      data_mat[,i] <- co
    }
    colnames(data_mat) <- names(annotationListGR)
  }
  else{
    for(i in 1:length(annotationListGR)){
      cp <- percent_func(dat_mat_gr, annotationListGR[[i]])
      data_mat[,i] <- cp
    }
    colnames(data_mat) <- names(annotationListGR)
  }
  
  return(data_mat)
}

################################################################################################################

#function to extract unique chromosome specific TAD boundaries
extract_boundaries_func <- function(domains.mat, preprocess, chromosome, resolution){
  domains.mat <- domains.mat[,1:3]
  domains.mat[,1] <- paste0("chr",domains.mat[,1])
  colnames(domains.mat) <- c("Chromosome", "Start", "End")
  if(preprocess==TRUE){
    ##restricting domain data to TADs > 2*resolution and < 2,000,000
    domains.mat <- domains.mat[which((domains.mat$End - domains.mat$Start)>(2*resolution) & (domains.mat$End - domains.mat$Start)<2000000),]
    }
  
  ##concatenating boundary coordinates into one long vector
  coords <- domains.mat
  colnames(coords)[2:3] <- c("coordinate", "coordinate")
  coords <- rbind.data.frame(coords[,c(1,2)],coords[,c(1,3)])
  coords$Chromosome <- as.character(coords$Chromosome)
  
  ##sorting by chromosome and coordinate
  coords <- coords[order(as.numeric(substr(coords$Chromosome,4,5)), coords$coordinate, decreasing = FALSE),]
  
  ##remove duplicates for coordinates that are conjoined
  coords <- coords[!duplicated(coords),]
  
  bounds <- GRanges(seqnames=coords$Chromosome, ranges=IRanges(start=coords$coordinate, width=1))
  
  #chromsome specific boundary data
  bounds_chr_specific <- bounds[which(seqnames(bounds)==tolower(chromosome))]
  
  return(bounds_chr_specific)
}


##################################################################################################

#function to predict cl, resolution, and chromosome specific tad boundaries using optimal parameters
predict_tad_bounds_func <- function(bounds.GR, datamatrix, chromosome, sampling="smote", metric="MCC", seed=123, crossvalidation=TRUE, number=10){
  
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
  
  # set number of iterations
  samps = 50
  
  if(sampling=="ros"){
    #assign sample indeces
    sampids <- matrix(ncol=samps, 
                      nrow=length(train$y[which(train$y=="No")]))
    
    #filling in the sample ids matrix
    set.seed(seed)
    for(j in 1:samps){
      sampids[,j] <- sample(x = which(train$y=="Yes"),
                            size = length(which(train$y=="No")),
                            replace = TRUE)
    }
    train <- rbind.data.frame(train[which(train$y=="No"),],
                              train[sampids[,1],])
    #Randomly shuffle the data
    set.seed(seed)
    train <- train[sample(nrow(train)),]
  }else if(sampling=="rus"){
    #assign sample indeces
    sampids <- matrix(ncol=samps, 
                      nrow=length(train$y[which(train$y=="Yes")]))
    
    #filling in the sample ids matrix
    set.seed(seed)
    for(j in 1:samps){
      sampids[,j] <- sample(x = which(train$y=="No"),
                            size = length(which(train$y=="Yes")),
                            replace = FALSE)
    }
    train <- rbind.data.frame(train[which(train$y=="Yes"),],
                              train[sampids[,1],])
    
    #Randomly shuffle the data
    set.seed(seed)
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
    method = if(crossvalidation==TRUE){print("cv")}else{print("none")},
    number = if(crossvalidation==TRUE){number}else{print(NULL)},
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
  return(tadModel)
}

############################################################################################

#function to predict TAD boundaries at bp resolution
predict_at_bp_resolution_func <- function(bounds.GR, resolution, chromosome, annotationListGR, tadModel, predwindow=10, region=TRUE){
  #creating chromosome specific test data at bp resolution
  
  #chromosome-specific test data is at basepair resolution
  genome <- getBSgenome("hg19")
  seqlength <- genome@seqinfo@seqlengths[which(genome@seqinfo@seqnames==tolower(chromosome))]
  ##start = first chromosome specific coordinate
  ##end = last chromosome specific coordinate that is a factor of resolution
  start = 0
  end = seqlength - (seqlength %% resolution)
  test_data <- data.frame(baseNum = seq.int(from=start, to=end, by=1))
  
  #annotate where the observed TAD boundary is according to TAD caller
  test_data$y <- ifelse(test_data$baseNum %in% start(bounds.GR), 1, 0)
  
  #create separate granges object to use to calculate distance to center of annotation
  test_bins_center <- GRanges(seqnames=tolower(chromosome),
                              IRanges(start=test_data$baseNum,
                                      end=test_data$baseNum))
  
  for(i in 1:length(annotationListGR)){
    #print(temp[i])
    dat_gr <- annotationListGR[[i]]
    
    dat_gr_center <- resize(dat_gr, width = 1, fix = "center")
    
    ####calculating distance
    d <- distance_func(test_bins_center, dat_gr_center)
    test_data <- cbind.data.frame(test_data,d)
  }
  
  colnames(test_data)[-c(1,2)] <- colnames(tadModel$trainingData)[-1]
  
  #take log of features
  test_data[,-c(1,2)] <- apply(test_data[,-c(1,2)],2,function(x){log(x + 1, base=2)})
  
  #factor the observed TAD boundary indicator
  test_data$y <- factor(test_data$y)
  levels(test_data$y) <- c("No", "Yes")
  
  #predict location of TAD boundary at basepair resolution using optimal model
  ##only include the data matrix part of the test data
  pred.tadModel <- predict(tadModel,newdata=test_data[,-c(1,2)],type="prob")[,"Yes"]
  
  #create data frame with basepair coordinates, predicted value of TAD boundary, observed TAD boundary
  predprob_df <- data.frame(baseNum = test_data$baseNum,
                            PredProb = pred.tadModel,
                            trueBound = test_data$y)
  
  #indentifying where the model gives 100% probability of TAD boundary location
  ##taking the center of areas where multiple adjacent basepairs correspond to 100% probability
  predprob_df_1 <- predprob_df[which(predprob_df$PredProb==1),]
  lagpad <- function(x, k) {
    if (k>0) {
      return (c(rep(NA, k), x)[1 : length(x)] );
    }
    else {
      return (c(x[(-k+1) : length(x)], rep(NA, -k)));
    }
  }
  predprob_df_1$lagbaseNum <- lagpad(predprob_df_1$baseNum, 1)
  predprob_df_1$predprob1length <- predprob_df_1$baseNum - predprob_df_1$lagbaseNum  
  
  predprob_df_1$retain <-  NA
  predprob_df_1$retain[1] <- 1
  x <- 1
  #consider where groups of probabilities==1 are not adjacent
  ##set the window as the average number of basepairs for nonadjacent probabilities==1
  ###all basepairs within that window are considered a predicted TAD boundary region (PTBR)
  for(i in 2:(dim(predprob_df_1)[1])){
    if(predprob_df_1$predprob1length[i]<=predwindow){predprob_df_1$retain[i] = x}else{x = x+1; predprob_df_1$retain[i] = x}
  }
  
  #consider the middle of a PTR as a predicted TAD boundary point (PTBP)
  prob1intervals <- unique(predprob_df_1$retain)
  mid <- numeric()
  for(i in prob1intervals){
    dat <- predprob_df_1[which(predprob_df_1$retain==prob1intervals[i]),]
    n=dim(dat)[1]
    mid[i] <- ceiling((dat$baseNum[1]+dat$baseNum[n])/2)
  }
  
  predprob_df$predBound <- ifelse(predprob_df$baseNum %in% mid, "Yes", "No")
  
  #creating a list with two granges objects
  ##trueBound_gr: chromosome-specific observed TAD boundaries 
  ##predBound_gr: chromosome-specific predicted TAD boundaries 
  trueBound_gr <- GRanges(seqnames=tolower(chromosome),
                          IRanges(start=predprob_df$baseNum[which(predprob_df$trueBound=="Yes")],
                                  end=predprob_df$baseNum[which(predprob_df$trueBound=="Yes")]))
  predBound_gr <- GRanges(seqnames=tolower(chromosome),
                          IRanges(start=predprob_df$baseNum[which(predprob_df$predBound=="Yes")],
                                  end=predprob_df$baseNum[which(predprob_df$predBound=="Yes")]))
  
  #creating a list with two granges objects
  ##trueBound_gr_region: chromosome-specific observed TAD boundaries flanked by resolution
  ##predBound_gr_region: chromosome-specific predicted TAD boundary regions given by predwindow
  trueBound_gr_region <- flank(trueBound_gr, width=10000, both=TRUE)
  predBound_gr_region <- GRangesList()
  for(i in 1:(length(table(predprob_df_1$retain)))){
    predBound_gr_region[[i]] <- GRanges(seqnames=tolower(chromosome), IRanges(start=predprob_df_1$baseNum[which(predprob_df_1$retain==i)][1],
                                                                              end=predprob_df_1$baseNum[which(predprob_df_1$retain==i)][length(predprob_df_1$baseNum[which(predprob_df_1$retain==i)])]))
  }
  predBound_gr_region <- do.call("c", predBound_gr_region)
  
  if(region==TRUE){
    finer_res_list <- list(trueBound_gr_region,predBound_gr_region)
  }else{finer_res_list <- list(trueBound_gr,predBound_gr)}
  
  return(finer_res_list)
  
}


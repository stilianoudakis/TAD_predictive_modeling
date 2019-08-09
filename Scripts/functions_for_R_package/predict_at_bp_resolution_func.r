#function to predict TAD boundaries at bp resolution
predict_at_bp_resolution_func <- function(bounds.GR, resolution, chromosome, seqList, annotationListGR, tadModel){
  #creating chromosome specific test data at bp resolution
  
  #chromosome-specific test data is at basepair resolution
  #genome <- getBSgenome("hg19")
  #seqlength <- genome@seqinfo@seqlengths[which(genome@seqinfo@seqnames==tolower(chromosome))]
  ##start = first chromosome specific coordinate
  ##end = last chromosome specific coordinate that is a factor of resolution
  #start = 0
  #end = seqlength - (seqlength %% resolution)
  #test_data <- data.frame(baseNum = seq.int(from=start, to=end, by=1))
  
  
  test_data <- data.frame(baseNum = seqList[[as.numeric(gsub("CHR", "", chromosome))]])
  
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
  lagpad <- function(x, k){
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
  ##set the window as 1 base pair 
  ###all basepairs within that window are considered a predicted TAD boundary region (PTBR)
  for(i in 2:(dim(predprob_df_1)[1])){
    if(predprob_df_1$predprob1length[i]<=1){predprob_df_1$retain[i] = x}else{x = x+1; predprob_df_1$retain[i] = x}
  }
  
  #consider the middle of a PTR as a predicted TAD boundary point (PTBP)
  prob1intervals <- unique(predprob_df_1$retain)
  mid <- numeric()
  for(i in prob1intervals){
    dat <- predprob_df_1[which(predprob_df_1$retain==prob1intervals[i]),]
    n=dim(dat)[1]
    mid[i] <- ceiling((dat$baseNum[1]+dat$baseNum[n])/2)
  }
  
  #cluster the differences between midpoints
  x <- dist(mid, method = "euclidean")
  hc1 <- hclust(x, method = "complete")
  z <- sapply(2:(length(mid)-1), function(i) { 
    mean(silhouette(cutree(hc1, i), dist=x)[,"sil_width"]) })
  k=which.max(z)+1
  p <- pam(x, diss=TRUE, k=k, cluster.only=TRUE, trace.lev = 0)
  clustdf <- data.frame(mid=mid, clust=cutree(hc1, k=which.max(z)))
  
  prob1intervals <- unique(clustdf$clust)
  mid2 <- numeric()
  for(i in prob1intervals){
    dat <- clustdf[which(clustdf$clust==prob1intervals[i]),]
    n=dim(dat)[1]
    mid2[i] <- ceiling((dat$mid[1]+dat$mid[n])/2)
  }
  
  predprob_df$predBound <- ifelse(predprob_df$baseNum %in% mid2, "Yes", "No")
  
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
  #trueBound_gr_region <- flank(trueBound_gr, width=10000, both=TRUE)
  #predBound_gr_region <- GRangesList()
  #for(i in 1:(length(table(predprob_df_1$retain)))){
  #  predBound_gr_region[[i]] <- GRanges(seqnames=tolower(chromosome), IRanges(start=predprob_df_1$baseNum[which(predprob_df_1$retain==i)][1],
  #                                                                            end=predprob_df_1$baseNum[which(predprob_df_1$retain==i)][length(predprob_df_1$baseNum[which(predprob_df_1$retain==i)])]))
  #}
  #predBound_gr_region <- do.call("c", predBound_gr_region)
  
  #if(region==TRUE){
  #  finer_res_list <- list(predprob_df_1, predprob_df, trueBound_gr_region,predBound_gr_region)
  #}else{finer_res_list <- list(predprob_df_1, predprob_df, trueBound_gr,predBound_gr)}
  #names(finer_res_list)[3:4] <- c("Called", "Predicted")
  
  trueBound_gr_region <- flank(trueBound_gr, width=resolution, both=TRUE)
  
  finer_res_list <- list(predprob_df_1, predprob_df, trueBound_gr_region, predBound_gr)
  return(finer_res_list)
  
}

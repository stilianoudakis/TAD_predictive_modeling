#function to predict TAD boundaries at bp resolution
predict_at_bp_resolution_func <- function(bounds.GR, resolution, chromosome, seqData, annotationListGR, tadModel){
  resolution=as.integer(resolution)
  
  test_data <- matrix(nrow=length(seqData),
                      ncol=length(annotationListGR))
  
  seqData <- c(seqData[1]-1, seqData)
  
  d <- lapply(annotationListGR, function(x){distance_func(GRanges(seqnames=tolower(chromosome),
                                                                  IRanges(start=seqData,
                                                                          end=seqData)), 
                                                          x)})
  for(i in 1:length(annotationListGR)){
    test_data[,i] <- d[[i]]
  }
  
  colnames(test_data) <- names(annotationListGR)
  test_data <- apply(test_data,2,function(x){log(x + 1, base=2)})
  
  pred.tadModel <- predict(tadModel,newdata=test_data,type="prob")[,"Yes"]
  
  prob1basenumdiff <- diff(seqData[which(pred.tadModel==1)])
  retain <- numeric()
  x <- 1
  for(i in 1:length(prob1basenumdiff)){
    if(prob1basenumdiff[i]<=1){retain[i] = x}else{x = x+1; retain[i] = x}
  }
  retain = c(1,retain)
  
  mid <- numeric()
  for(i in unique(retain)){
    #print(i)
    bpdat <- seqData[which(pred.tadModel==1)][which(retain==i)]
    n=length(bpdat)
    mid[i] <- ceiling((bpdat[1]+bpdat[n])/2)
  }
  
  x <- dist(mid, method = "euclidean")
  hc1 <- hclust(x, method = "complete")
  z <- sapply(2:(length(mid)-1), function(i) { 
    mean(silhouette(cutree(hc1, i), dist=x)[,"sil_width"]) })
  k=which.max(z)+1
  c <- clara(mid, 
             k=k, 
             samples=100,
             metric = "euclidean",
             stand = FALSE, 
             trace = 2, 
             medoids.x = TRUE)
  
  predBound_gr <- GRanges(seqnames=tolower(chromosome),
                          IRanges(start=seqData[which(seqData %in% c$medoids)],
                                  end=seqData[which(seqData %in% c$medoids)]))
  trueBound_gr <- GRanges(seqnames=tolower(chromosome),
                          IRanges(start=seqData[which(seqData %in% start(bounds.GR))],
                                  end=seqData[which(seqData %in% start(bounds.GR))]))
  
  trueBound_gr_region <- flank(trueBound_gr, width=resolution, both=TRUE)
  
  finer_res_list <- list(trueBound_gr_region, predBound_gr)
  return(finer_res_list)
}


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
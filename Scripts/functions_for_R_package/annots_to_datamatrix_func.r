#function to create data matrix for training predictive model
annots_to_datamatrix_func <- function(seqData, resolution, predictortype="distance", annotationListGR, chromosome){
  resolution=as.integer(resolution)
  
  seqData <- c(seqData[1]-1, seqData)
  
  start=seqData[1] 
  start=0
  end=seqData[length(seqData)] - (seqData[length(seqData)] %% resolution)
  rows = seqData[seqData %in% seq(start, end, resolution)]
  
  data_mat <- matrix(nrow=length(rows),
                     ncol=length(annotationListGR))
  rownames(data_mat) <- rows
  
  dat_mat_gr <- flank(GRanges(seqnames = tolower(chromosome),
                              IRanges(start = as.numeric(rownames(data_mat)),
                                      #end = as.numeric(rownames(data_mat)),
                                      width=1)), (resolution/2), both=TRUE)
  
  if(predictortype=="distance"){
    for(i in 1:length(annotationListGR)){
      d <- distance_func(dat_mat_gr, annotationListGR[[i]])
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
  }else if(predictortype=="count"){
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
  
  return(data_mat)
}

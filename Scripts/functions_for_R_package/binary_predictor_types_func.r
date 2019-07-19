####calculating binary overlaps
binary_func <- function(binned_data_gr, annot_data_gr){
  #Finding the total number of overlaps between genomic bins and the specific genomic annotation
  count_binary <- countOverlaps(binned_data_gr, annot_data_gr)
  
  #Binarizing the overlap (1 or 0)
  count_binary <- ifelse(count_binary>0,1,0)
  
  return(count_binary)
}

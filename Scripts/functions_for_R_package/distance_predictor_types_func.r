####calculating distance
distance_func <- function(binned_data_center_gr, annot_data_center_gr){
  #distance from center of genomic bin to nearest genomic region of interest
  dist <- mcols(distanceToNearest(binned_data_center_gr, annot_data_center_gr))$distance
  
  return(dist)
}
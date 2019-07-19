####calculating count overlaps
count_func <- function(binned_data_gr, annot_data_gr){
  #Finding the total number of overlaps between genomic bins and the specific genomic annotation
  count_total <- countOverlaps(binned_data_gr, annot_data_gr)
  
  return(count_total)
}

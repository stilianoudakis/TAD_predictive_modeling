
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
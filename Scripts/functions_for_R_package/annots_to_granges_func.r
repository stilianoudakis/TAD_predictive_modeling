
#function to create grangeslist object from functional genomic annotation .bed files
annots_to_granges_func <- function(filepath){
  annots_gr_list <- GRangesList() 
  #if(grepl(paste0(names, collapse = "|"), list.files(filepath)))
  for(i in 1:length(list.files(filepath, pattern = "*.bed"))){
    dat <- read.table(paste0(filepath,"/",list.files(filepath, pattern = "*.bed")[i]),
                      header = FALSE)
    
    annots_gr_list[[i]] <- GRanges(seqnames=dat$V1,
                                   IRanges(start=dat$V2, 
                                           end=dat$V3))
    mcols(annots_gr_list[[i]])$coverage <- dat[,4]
  }
  names(annots_gr_list) <- gsub(".bed", "", list.files(filepath, pattern = "*.bed"))
  return(annots_gr_list)
}

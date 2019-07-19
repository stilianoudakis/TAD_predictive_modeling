extract_boundaries_func <- function(domains.mat, preprocess, chromosome){
  if(preprocess==TRUE){
    ##restricting domain data to TADs > 2*resolution and < 2,000,000
    domains <- domains[which((domains$End - domains$Start)>(2*resolution) & (domains$End - domains$Start)<2000000),]
    }
  
  ##concatenating boundary coordinates into one long vector
  coords <- domains
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
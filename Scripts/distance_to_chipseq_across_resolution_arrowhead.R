source("C:/Users/stili/Documents/preciseTAD_functions/functions/annots_to_granges_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/binary_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/count_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/distance_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/extract_boundaries_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/percent_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/TADRF.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/preciseTAD.r")


file_path="Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS"
annots.GR <- annots_to_granges_func(filepath=file_path, signal=4)
chromosomes <- paste0("CHR", c(1:8,10:22))

#l <- list()
#res <- c("5000","10000","25000","50000","100000")
#for(i in 1:5){
#  domains <- read.table(paste0("Z:/TAD_data_analysis/GSE63525_data/arrowhead_output/GM12878/GM12878_domain_data_",
#                               res[i],
#                               ".b.txt"))
#  
#  bounds <- extract_boundaries_func(domains.mat=domains, 
#                                    preprocess=FALSE, 
#                                    CHR=chromosomes, 
#                                    resolution=res[i])
#  
#  d <- matrix(nrow=length(bounds), ncol=4)
#  colnames(d) <- c("CTCF","RAD21","SMC3","ZNF143") 
#  for(j in 1:4){
#    #print(c(i,j))
#    d[,j] <- log2(mcols(distanceToNearest(bounds, annots.GR[[j]]))$distance+1)
#  }
#  
#  l[[i]] <- d
#  
#}
#names(l) <- c("5 kb","10 kb","25 kb","50 kb","100 kb")


boundslist <- list()
for(i in 1:5){
  domains <- read.table(paste0("Z:/TAD_data_analysis/GSE63525_data/arrowhead_output/GM12878/GM12878_domain_data_",
                               res[i],
                               ".b.txt"))
  
  bounds <- extract_boundaries_func(domains.mat=domains, 
                                    preprocess=FALSE, 
                                    CHR=chromosomes, 
                                    resolution=res[i])
  boundslist[[i]] <- bounds
}  

l <- list()
for(i in 1:5){
  for(j in 1:4){
    s <- length(log2(mcols(distanceToNearest(annots.GR[[j]],boundslist[[i]]))$distance+1))
    d <- matrix(nrow=s, ncol=5)
    colnames(d) <- c("5 kb","10 kb","25 kb","50 kb","100 kb")
    for(h in 1:5){
      d[,h] <- log2(mcols(distanceToNearest(annots.GR[[j]],boundslist[[h]]))$distance+1)
    }
    l[[j]] <- d
  }
}
names(l) <- c("CTCF","RAD21","SMC3","ZNF143") 


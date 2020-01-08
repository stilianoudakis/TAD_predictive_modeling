source("C:/Users/stili/Documents/preciseTAD_functions/functions/annots_to_granges_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/binary_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/count_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/distance_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/extract_boundaries_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/percent_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/TADRF.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/preciseTAD.r")

# GM12878

bounds.list <- list()
jaccard.mat <- matrix(nrow=5,ncol=5)
rownames(jaccard.mat) <- colnames(jaccard.mat) <- c("5 kb", "10 kb", "25 kb", "50 kb", "100 kb")

for(i in c("5000","10000","25000","50000","100000")){
  domains <- read.table(paste0("Z:/TAD_data_analysis/GSE63525_data/arrowhead_output/GM12878/GM12878_domain_data_",
                               i,
                               ".b.txt"))
  
  bounds <- extract_boundaries_func(domains.mat=domains, 
                                    preprocess=FALSE, 
                                    CHR=CHR, 
                                    resolution=i)
  
  bounds.list[[i]] <- flank(bounds, 5000, both=TRUE)
  
}

for(i in 1:5){
  #	for(j in 1:5){
  
  jaccard.mat[,i] <- unlist(
    lapply(
      bounds.list,
      function(x){eval(R_bedtools_jaccard(x,bounds.list[[i]]))@listData$jaccard}))
  
  #	}
}


col3 <- colorRampPalette(c("blue", "white", "red"))
#col3 <- colorRampPalette(c("white", "red"))		
corrplot(jaccard.mat, 
         method = "color",
         type="upper", 
         #order="hclust", 
         tl.col="black", 
         tl.srt=90,
         #tl.cex = .4,
         col=col3(200),
         #na.label.col = "gray",
         na.label = "1",
         diag=FALSE,
         outline = TRUE,
         addCoef.col = "black", 
         number.digits = 2,
         cl.pos="n",
         mar = c(0, 0, 0, 0))


jaccard.mat.heat <- jaccard.mat

diag(jaccard.mat.heat) <- NA

distance.col = as.dist(jaccard.mat.heat)
cluster.col = hclust(distance.col, method = "average")
my_palette <- colorRampPalette(c("white", "pink", "red"))(ncol(jaccard.mat.heat))
heatmap.2(t(jaccard.mat.heat),
          labRow=rownames(jaccard.mat.heat),
          labCol=colnames(jaccard.mat.heat),
          #cellnote = ,  # same data set for cell labels
          #main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins=c(12,8),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          na.color = "gray",
          #breaks=col_breaks,    # enable color transition at specified limits
          key=FALSE,
          #keysize = 1.5,
          #key.title = "Elastic-Net Coefficient",
          dendrogram="none",     # only draw a row dendrogram,
          Rowv = FALSE,
          Colv = FALSE,
          cexRow=1,
          cexCol=1,
          rowsep=0:nrow(t(jaccard.mat.heat)),
          colsep=0:ncol(t(jaccard.mat.heat)),
          sepwidth=c(0.00001,0.00001),
          sepcolor="black")

################################################################################################

# K562

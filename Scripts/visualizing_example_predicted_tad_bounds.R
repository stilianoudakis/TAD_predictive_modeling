library(BSgenome)
library(randomForest)
library(DMwR)
library(cluster)
library(pROC)
library(PRROC)
library(caret)
library(GenomicRanges)

# Loading necessary functions

source("C:/Users/stili/Documents/preciseTAD_functions/functions/annots_to_granges_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/binary_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/count_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/distance_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/signal_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/extract_boundaries_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/percent_func.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/TADRF.r")
source("C:/Users/stili/Documents/preciseTAD_functions/functions/preciseTAD.r")

# Extracting TAD boundaries using the extract_boundaries_func function for ARROWHEAD defined TADs at 10kb resolution

domains <- read.table("Z:/TAD_data_analysis/GM12878/10kb/GM12878_domain_data_10000.b.txt")
bounds.GR <- extract_boundaries_func(domains.mat=domains, 
                                     preprocess=FALSE, 
                                     CHR="CHR22", 
                                     resolution=10000)

# Obtaining list of GRanges objects of ChIP-seq data used for modelling

genomicElements.GR <- annots_to_granges_func(filepath = "Z:/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS",
                                             pattern="*.bed",
                                             signal=4)

tadModel <- TADRF(bounds.GR=bounds.GR,
                  resolution=10000,
                  genomicElements.GR=genomicElements.GR,
                  featureType="distance",
                  resampling="smote",
                  trainCHR="CHR22",
                  predictCHR="CHR22",
                  number=5,
                  ntrees=500,
                  metric="MCC",
                  verbose=TRUE,
                  seed=123,
                  model=TRUE,
                  importances=TRUE,
                  performances=TRUE)

bounds.GR=bounds.GR;
resolution=10000;
genomicElements.GR=genomicElements.GR;
featureType="distance";
CHR="CHR22";
chromCoords=list(17389000,18011000);
tadModel=tadModel[[1]];
threshold="roc";
samples=100;
flank=NULL;
verbose=TRUE;
seed=123

resolution=as.integer(resolution)

getChrLength <- function(genome = "BSgenome.Hsapiens.UCSC.hg19"){
  g <- getBSgenome(genome, masked=FALSE)
  data.frame(chrom=1:24, length=seqlengths(g)[1:24])
}
.chrAsNum <- function(tbl){
  tbl$chrom <- gsub("chr", "", tbl$chrom)
  tbl$chrom[tbl$chrom=="X"] <- 23
  tbl$chrom[tbl$chrom=="Y"] <- 24
  tbl$chrom <- as.numeric(tbl$chrom)
  tbl[order(tbl$chrom),]
}
getCentromeres <- function( genome="hg19" ){
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case of failure, try another mirror
  if(inherits(mySession, "try-error"))
    mySession <- browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  obj <- ucscTableQuery(mySession, table="gap")
  tbl <- getTable(obj)
  tbl <- tbl[tbl$type=="centromere", c("chrom", "chromStart", "chromEnd")]
  colnames(tbl)[2:3] <- c("centromerStart", "centromerEnd")
  .chrAsNum(tbl)
}
makeHg19 <- function(){
  tbl <- merge(getChrLength(), getCentromeres(), by="chrom")
  cumlen <- c(0, cumsum(as.numeric(tbl$length))[-nrow(tbl)])
  cbind.data.frame(tbl, cumlen=cumlen)    
}
hg19 <- makeHg19()
hg19$chrom <- paste0("CHR",hg19$chrom)

##################################
#CREATING BP RESOLUTION TEST DATA#
##################################

if(class(chromCoords)=="list"){
  seqDataTest <- c(chromCoords[[1]]:chromCoords[[2]])
  
  if("TRUE" %in% table(seqDataTest %in% c(hg19$centromerStart[hg19$chrom==CHR]:hg19$centromerEnd[hg19$chrom==CHR]))){
    centromereTestStart <- hg19$centromerStart[hg19$chrom==CHR]
    centromereTestEnd <- hg19$centromerEnd[hg19$chrom==CHR]
    seqDataTest <- seqDataTest[-which(seqDataTest %in% c(centromereTestStart:centromereTestEnd))]
  }
  
  
  
}else{
  seqLengthTest <- hg19$length[hg19$chrom==CHR]
  
  seqDataTest <- c(0:seqLengthTest)
  centromereTestStart <- hg19$centromerStart[hg19$chrom==CHR]
  centromereTestEnd <- hg19$centromerEnd[hg19$chrom==CHR]
  seqDataTest <- seqDataTest[-which(seqDataTest %in% c(centromereTestStart:centromereTestEnd))]
  
}

test_data <- matrix(nrow=length(seqDataTest),
                    ncol=length(genomicElements.GR))

if(featureType=="distance"){
  p <- lapply(genomicElements.GR, function(x){distance_func(GRanges(seqnames=tolower(CHR),
                                                                    IRanges(start=seqDataTest,
                                                                            end=seqDataTest)), 
                                                            x)})
}else if(featureType=="binary"){
  p <- lapply(genomicElements.GR, function(x){binary_func(GRanges(seqnames=tolower(CHR),
                                                                  IRanges(start=seqDataTest,
                                                                          end=seqDataTest)), 
                                                          x)})
}else if(featureType=="oc"){
  p <- lapply(genomicElements.GR, function(x){count_func(GRanges(seqnames=tolower(CHR),
                                                                 IRanges(start=seqDataTest,
                                                                         end=seqDataTest)), 
                                                         x)})
}else if(featureType=="op"){
  p <- lapply(genomicElements.GR, function(x){percent_func(GRanges(seqnames=tolower(CHR),
                                                                   IRanges(start=seqDataTest,
                                                                           end=seqDataTest)), 
                                                           x)})
}else{
  p <- lapply(genomicElements.GR, function(x){signal_func(GRanges(seqnames=tolower(CHR),
                                                                  IRanges(start=seqDataTest,
                                                                          end=seqDataTest)), 
                                                          x)})
}

for(i in 1:length(genomicElements.GR)){
  test_data[,i] <- p[[i]]
}
colnames(test_data) <- names(genomicElements.GR)

if(featureType=="distance"){test_data <- apply(test_data,2,function(x){log(x + 1, base=2)})}

##################################
#PREDICTING AT BP RESOLUTION     #
##################################

chunk <- 20000000
n <- nrow(test_data)
r <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
d <- split.data.frame(test_data,r)

counter=0
if(verbose==TRUE){
  predictions <- lapply(d, function(x){
    counter <<- counter + 1;
    print(paste0("Obtaining predicted probabilities for chunck #", names(d)[counter], " of test data out of ", length(d), " chunks, each with 20 million rows"));
    p <- predict(tadModel,newdata=x,type="prob")[,"Yes"];
    return(p)})
}else{
  predictions <- lapply(d, function(x){
    #counter <<- counter + 1;
    #print(paste0(names(d)[counter]));
    p <- predict(tadModel,newdata=x,type="prob")[,"Yes"];
    return(p)})
}

predictions <- do.call("c", predictions)

test_data_Y <- ifelse(seqDataTest %in% start(bounds.GR), 1, 0)

if(threshold=="roc"){
  ss <- pROC::roc(test_data_Y, predictions, quiet=TRUE)
  t <- ss$thresholds[which.max(ss$specificities+ss$sensitivities)]
}else{
  t <- threshold 
}

rm("test_data")

prob1basenumdiff <- diff(seqDataTest[which(predictions>=t)])
retain <- numeric()
x <- 1
for(i in 1:length(prob1basenumdiff)){
  if(prob1basenumdiff[i]<=1){retain[i] = x}else{x = x+1; retain[i] = x}
}
retain = c(1,retain)

mid <- numeric()
for(i in unique(retain)){
  #print(i)
  bpdat <- seqDataTest[which(predictions>=t)][which(retain==i)]
  n=length(bpdat)
  mid[i] <- ceiling((bpdat[1]+bpdat[n])/2)
}

rm("predictions")

x <- dist(mid, method = "euclidean")
hc1 <- hclust(x, method = "complete")
z <- sapply(2:(length(mid)-1), function(i) { 
  mean(silhouette(cutree(hc1, i), dist=x)[,"sil_width"]) })
k=which.max(z)+1

if(verbose==TRUE){
  print(paste0("Initializing CLARA with ", k, " clusters"));
  c <- clara(mid, 
             k=k, 
             samples=samples,
             metric = "euclidean",
             stand = FALSE, 
             trace = 2, 
             medoids.x = TRUE)
}else{
  c <- clara(mid, 
             k=k, 
             samples=samples,
             metric = "euclidean",
             stand = FALSE, 
             trace = 0, 
             medoids.x = TRUE)
}


predprobdata <- data.frame(basenum=c(17389000:18011000),
                           prob=predictions)
predprobdata$predicted <- ifelse(predprobdata$basenum %in% c$medoids, "Yes", "No")
predprobdata$called <- ifelse(predprobdata$basenum %in% start(trueBound_gr), "Yes", "No")

calleddata <- predprobdata[which(predprobdata$called=="Yes"),]

predicteddata <- predprobdata[which(predprobdata$prob>=t),]

predicteddata2 <- predprobdata[which(predprobdata$predicted=="Yes"),]

ggplot() + 
  geom_vline(data=predicteddata, aes(xintercept=basenum, 
                                     color="red"),
             size=.5,
             show.legend = TRUE)+
  geom_line(data=predprobdata, aes(x=basenum, y=prob), 
            color='black',
            size=1)+
  geom_vline(data=calleddata, aes(xintercept=basenum, 
                                  color="blue"),
             size=.75,
             show.legend = TRUE)+
  geom_vline(data=predicteddata2, aes(xintercept=basenum, 
                                      color="green"),
             size=.75,
             show.legend = TRUE)+
  xlab("Base Pair Coordinate") +
  ylab("Probability") +
  scale_color_manual(name = '', 
                     values =c("blue", "green", "red"),
                     labels = c("Called Boundary",
                                "PTBP",
                                "PTBA"))+
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        #strip.text.x = element_text(size = 15),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


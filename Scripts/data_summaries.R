
# GM12878

CHR = paste0("CHR", c(1:8,10:22))
numtads <- numeric()
numbounds <- numeric()
numbins <- numeric()
imbalance <- numeric()
r = c("5000","10000","25000","50000","100000")

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

seqLength <- hg19$length[hg19$chrom %in% CHR]
seqDataList <- list()
for(k in 1:length(seqLength)){
  seqData <- c(0:seqLength[k])
  centromereStart <- hg19$centromerStart[hg19$chrom==CHR[k]]
  centromereEnd <- hg19$centromerEnd[hg19$chrom==CHR[k]]
  seqDataList[[k]] <- seqData[-which(seqData %in% c(centromereStart:centromereEnd))]
}
names(seqDataList) <- CHR

for(i in 1:length(r)){
  resolution = as.integer(r[i])
  
  domains <- read.table(paste0("Z:/TAD_data_analysis/GSE63525_data/arrowhead_output/GM12878/GM12878_domain_data_",
                               r[i],
                               ".b.txt"))
  
  domains <- domains[which(domains$V1 %in% as.numeric(gsub("CHR","",CHR))),]
  
  numtads[i] <- nrow(domains)
  
  bounds <- extract_boundaries_func(domains.mat=domains, 
                                    preprocess=FALSE, 
                                    CHR=CHR, 
                                    resolution=resolution)
  
  numbounds[i] <- length(bounds)
  
  bounds.GR.flank <- flank(bounds, width=(resolution/2), both=TRUE)
  
  data_mat_list <- list()
  outcome_list <- list()
  #full_list <- list()
  #train_list <- list()
  #test_list <- list()
  
  for(j in 1:length(CHR)){
    #start = 0
    start = (resolution/2)
    end = seqDataList[[j]][length(seqDataList[[j]])] - (seqDataList[[j]][length(seqDataList[[j]])] %% resolution) + (resolution/2)
    #end = seqDataList[[j]][length(seqDataList[[j]])] - (seqDataList[[j]][length(seqDataList[[j]])] %% (resolution/2))
    rows = seqDataList[[j]][seqDataList[[j]] %in% seq(start, end, resolution)]
    
    #data_mat_list[[j]] <- matrix(nrow=length(rows),
    #                             ncol=length(genomicElements.GR))
    #rownames(data_mat_list[[j]]) <- rows
    
    outcome_list[[j]] <- countOverlaps(GRanges(seqnames = tolower(CHR[j]),
                                               IRanges(#start = as.numeric(rownames(data_mat_list[[j]])),
                                                       start = rows,
                                                       width = resolution)), bounds.GR.flank)
    outcome_list[[j]] <- ifelse(outcome_list[[j]],1,0)
    outcome_list[[j]] <- factor(outcome_list[[j]])
    levels(outcome_list[[j]]) <- c("No", "Yes")
    
    #full_list[[j]] <- cbind.data.frame(outcome_list[[j]], data_mat_list[[j]])
    #names(full_list[[j]])[1] <- "y"
  }
  ci <- unlist(outcome_list)
  numbins[i] <- length(ci)
  imbalance[i] <- as.numeric(round(table(ci)[2]/sum(table(ci)),2))
}

data.frame(resolution = c("5 kb",
                          "10 kb",
                          "25 kb",
                          "50 kb",
                          "100 kb"),
           numtads=numtads,
           numbounds=numbounds,
           numbins=numbins,
           imbalance=imbalance)


#########################################################################################

# K562 

CHR = paste0("CHR", c(1:8,10:22))
numtads <- numeric()
numbounds <- numeric()
numbins <- numeric()
imbalance <- numeric()
r = c("5000","10000","25000","50000","100000")

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

seqLength <- hg19$length[hg19$chrom %in% CHR]
seqDataList <- list()
for(k in 1:length(seqLength)){
  seqData <- c(0:seqLength[k])
  centromereStart <- hg19$centromerStart[hg19$chrom==CHR[k]]
  centromereEnd <- hg19$centromerEnd[hg19$chrom==CHR[k]]
  seqDataList[[k]] <- seqData[-which(seqData %in% c(centromereStart:centromereEnd))]
}
names(seqDataList) <- CHR

for(i in 1:length(r)){
  resolution = as.integer(r[i])
  
  domains <- read.table(paste0("Z:/TAD_data_analysis/GSE63525_data/arrowhead_output/K562/K562_domain_data_",
                               r[i],
                               ".b.txt"))
  
  domains <- domains[which(domains$V1 %in% as.numeric(gsub("CHR","",CHR))),]
  
  numtads[i] <- nrow(domains)
  
  bounds <- extract_boundaries_func(domains.mat=domains, 
                                    preprocess=FALSE, 
                                    CHR=CHR, 
                                    resolution=resolution)
  
  numbounds[i] <- length(bounds)
  
  bounds.GR.flank <- flank(bounds, width=(resolution/2), both=TRUE)
  
  data_mat_list <- list()
  outcome_list <- list()
  #full_list <- list()
  #train_list <- list()
  #test_list <- list()
  
  for(j in 1:length(CHR)){
    #start = 0
    start = (resolution/2)
    end = seqDataList[[j]][length(seqDataList[[j]])] - (seqDataList[[j]][length(seqDataList[[j]])] %% resolution) + (resolution/2)
    #end = seqDataList[[j]][length(seqDataList[[j]])] - (seqDataList[[j]][length(seqDataList[[j]])] %% (resolution/2))
    rows = seqDataList[[j]][seqDataList[[j]] %in% seq(start, end, resolution)]
    
    #data_mat_list[[j]] <- matrix(nrow=length(rows),
    #                             ncol=length(genomicElements.GR))
    #rownames(data_mat_list[[j]]) <- rows
    
    outcome_list[[j]] <- countOverlaps(GRanges(seqnames = tolower(CHR[j]),
                                               IRanges(#start = as.numeric(rownames(data_mat_list[[j]])),
                                                 start = rows,
                                                 width = resolution)), bounds.GR.flank)
    outcome_list[[j]] <- ifelse(outcome_list[[j]],1,0)
    outcome_list[[j]] <- factor(outcome_list[[j]])
    levels(outcome_list[[j]]) <- c("No", "Yes")
    
    #full_list[[j]] <- cbind.data.frame(outcome_list[[j]], data_mat_list[[j]])
    #names(full_list[[j]])[1] <- "y"
  }
  ci <- unlist(outcome_list)
  numbins[i] <- length(ci)
  imbalance[i] <- as.numeric(round(table(ci)[2]/sum(table(ci)),2))
}

data.frame(resolution = c("5 kb",
                          "10 kb",
                          "25 kb",
                          "50 kb",
                          "100 kb"),
           numtads=numtads,
           numbounds=numbounds,
           numbins=numbins,
           imbalance=imbalance)

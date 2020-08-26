source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/preciseTAD.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/createTADdata.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/TADrandomForest.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/TADrfe.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/distance_func.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/percent_func.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/count_func.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/binary_func.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/signal_func.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/annots_to_granges_func.R")
source("/home/stilianoudakisc/TAD_data_analysis/functions_for_R_package/extract_boundaries_func.R")

#arrowhead

domains <- read.table("/home/stilianoudakisc/TAD_data_analysis/GSE63525_data/arrowhead_output/GM12878/GM12878_domain_data_5000.b.txt", header=F)
bounds.GR_gm12878 <- extract_boundaries_func(domains.mat=domains, 
                                     preprocess=FALSE, 
                                     CHR=paste0("CHR", c(1:8,10:22)), 
                                     resolution=5000)
genomicElements.GR_gm12878 <- annots_to_granges_func(filepath = "/home/stilianoudakisc/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS_nonames",
                                             pattern="*.bed",
                                             signal=4)

domains <- read.table("/home/stilianoudakisc/TAD_data_analysis/GSE63525_data/arrowhead_output/K562/K562_domain_data_5000.b.txt", header=F)
bounds.GR_k562 <- extract_boundaries_func(domains.mat=domains, 
                                             preprocess=FALSE, 
                                             CHR=paste0("CHR", c(1:8,10:22)), 
                                             resolution=5000)
genomicElements.GR_k562 <- annots_to_granges_func(filepath = "/home/stilianoudakisc/TAD_data_analysis/annotations/all_common_annotations/k562/topTFBS_nonames",
                                                     pattern="*.bed",
                                                     signal=4)

##################

chrs <- paste0("CHR", c(1:8,10:22))
aggPR <- list()
aggAUPRC <- list()
for(i in 1:length(chrs)){
  print(i)
  m <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                      chrs[i],
                      "/rus/distance/TADRF_holdout_gm12878_on_gm12878.rds"))
  
  testData_gm12878  <- createTADdata(bounds.GR=bounds.GR_gm12878,
                                    resolution=5000,
                                    genomicElements.GR=genomicElements.GR_gm12878,
                                    featureType="distance",
                                    resampling="none",
                                    trainCHR=chrs[i],
                                    predictCHR=NULL,
                                    seed=123)
  
  pred <- as.vector(predict(m[[1]],newdata=testData_gm12878[[1]][,-1],type="prob")[,"Yes"])
  fg <- pred[testData_gm12878[[1]]$y == "Yes"]
  bg <- pred[testData_gm12878[[1]]$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  aggPR[[i]] <- pr$curve[,-1]
  
  aggAUPRC[[i]] <- m[[3]][20,2]
}

saveRDS(aggPR, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/aggPR.rds")
saveRDS(aggAUPRC, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/aggAUPRC.rds")

################

chrs <- paste0("CHR", c(1:8,10:22))
agkPR <- list()
agkAUPRC <- list()
for(i in 1:length(chrs)){
  print(i)
  m <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb/results_by_chr/",
                      chrs[i],
                      "/rus/distance/TADRF_holdout_gm12878_on_k562.rds"))
  
  testData_k562  <- createTADdata(bounds.GR=bounds.GR_k562,
                                     resolution=5000,
                                     genomicElements.GR=genomicElements.GR_k562,
                                     featureType="distance",
                                     resampling="none",
                                     trainCHR=chrs[i],
                                     predictCHR=NULL,
                                     seed=123)
  
  pred <- as.vector(predict(m[[1]],newdata=testData_k562[[1]][,-1],type="prob")[,"Yes"])
  fg <- pred[testData_k562[[1]]$y == "Yes"]
  bg <- pred[testData_k562[[1]]$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  agkPR[[i]] <- pr$curve[,-1]
  
  agkAUPRC[[i]] <- m[[3]][20,2]
}

saveRDS(agkPR, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/agkPR.rds")
saveRDS(agkAUPRC, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/agkAUPRC.rds")

##################

chrs <- paste0("CHR", c(1:8,10:22))
akkPR <- list()
akkAUPRC <- list()
for(i in 1:length(chrs)){
  print(i)
  m <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/5kb/results_by_chr/",
                      chrs[i],
                      "/rus/distance/TADRF_holdout_k562_on_k562.rds"))
  
  testData_k562  <- createTADdata(bounds.GR=bounds.GR_k562,
                                     resolution=5000,
                                     genomicElements.GR=genomicElements.GR_k562,
                                     featureType="distance",
                                     resampling="none",
                                     trainCHR=chrs[i],
                                     predictCHR=NULL,
                                     seed=123)
  
  pred <- as.vector(predict(m[[1]],newdata=testData_k562[[1]][,-1],type="prob")[,"Yes"])
  fg <- pred[testData_k562[[1]]$y == "Yes"]
  bg <- pred[testData_k562[[1]]$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  akkPR[[i]] <- pr$curve[,-1]
  
  akkAUPRC[[i]] <- m[[3]][20,2]
}

saveRDS(akkPR, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/akkPR.rds")
saveRDS(akkAUPRC, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/akkAUPRC.rds")

################

chrs <- paste0("CHR", c(1:8,10:22))
akgPR <- list()
akgAUPRC <- list()
for(i in 1:length(chrs)){
  print(i)
  m <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/5kb/results_by_chr/",
                      chrs[i],
                      "/rus/distance/TADRF_holdout_k562_on_gm12878.rds"))
  
  testData_gm12878  <- createTADdata(bounds.GR=bounds.GR_gm12878,
                                  resolution=5000,
                                  genomicElements.GR=genomicElements.GR_gm12878,
                                  featureType="distance",
                                  resampling="none",
                                  trainCHR=chrs[i],
                                  predictCHR=NULL,
                                  seed=123)
  
  pred <- as.vector(predict(m[[1]],newdata=testData_gm12878[[1]][,-1],type="prob")[,"Yes"])
  fg <- pred[testData_gm12878[[1]]$y == "Yes"]
  bg <- pred[testData_gm12878[[1]]$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  akgPR[[i]] <- pr$curve[,-1]
  
  akgAUPRC[[i]] <- m[[3]][20,2]
}

saveRDS(akgPR, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/akgPR.rds")
saveRDS(akgAUPRC, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/akgAUPRC.rds")

#################################################################################################

#peakachu

domains <- readRDS("/home/stilianoudakisc/TAD_data_analysis/peakachu/peakachu_GM12878.rds")
bounds.GR_gm12878 <- extract_boundaries_func(domains.mat=domains, 
                                             preprocess=FALSE, 
                                             CHR=paste0("CHR", c(1:8,10:22)), 
                                             resolution=5000)
genomicElements.GR_gm12878 <- annots_to_granges_func(filepath = "/home/stilianoudakisc/TAD_data_analysis/annotations/all_common_annotations/gm12878/topTFBS_nonames",
                                                     pattern="*.bed",
                                                     signal=4)

domains <- readRDS("/home/stilianoudakisc/TAD_data_analysis/peakachu/peakachu_K562.rds")
bounds.GR_k562 <- extract_boundaries_func(domains.mat=domains, 
                                          preprocess=FALSE, 
                                          CHR=paste0("CHR", c(1:8,10:22)), 
                                          resolution=5000)
genomicElements.GR_k562 <- annots_to_granges_func(filepath = "/home/stilianoudakisc/TAD_data_analysis/annotations/all_common_annotations/k562/topTFBS_nonames",
                                                  pattern="*.bed",
                                                  signal=4)

################

chrs <- paste0("CHR", c(1:8,10:22))
pggPR <- list()
pggAUPRC <- list()
for(i in 1:length(chrs)){
  print(i)
  m <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                      chrs[i],
                      "/rus/distance/TADRF_holdout_gm12878_on_gm12878_peakachu.rds"))
  
  testData_gm12878  <- createTADdata(bounds.GR=bounds.GR_gm12878,
                                     resolution=5000,
                                     genomicElements.GR=genomicElements.GR_gm12878,
                                     featureType="distance",
                                     resampling="none",
                                     trainCHR=chrs[i],
                                     predictCHR=NULL,
                                     seed=123)
  
  pred <- as.vector(predict(m[[1]],newdata=testData_gm12878[[1]][,-1],type="prob")[,"Yes"])
  fg <- pred[testData_gm12878[[1]]$y == "Yes"]
  bg <- pred[testData_gm12878[[1]]$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  pggPR[[i]] <- pr$curve[,-1]
  
  pggAUPRC[[i]] <- m[[3]][20,2]
}

saveRDS(pggPR, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/pggPR.rds")
saveRDS(pggAUPRC, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/pggAUPRC.rds")

##################

chrs <- paste0("CHR", c(1:8,10:22))
pgkPR <- list()
pgkAUPRC <- list()
for(i in 1:length(chrs)){
  print(i)
  m <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/10kb/results_by_chr/",
                      chrs[i],
                      "/rus/distance/TADRF_holdout_gm12878_on_k562_peakachu.rds"))
  
  testData_k562  <- createTADdata(bounds.GR=bounds.GR_k562,
                                  resolution=5000,
                                  genomicElements.GR=genomicElements.GR_k562,
                                  featureType="distance",
                                  resampling="none",
                                  trainCHR=chrs[i],
                                  predictCHR=NULL,
                                  seed=123)
  
  pred <- as.vector(predict(m[[1]],newdata=testData_k562[[1]][,-1],type="prob")[,"Yes"])
  fg <- pred[testData_k562[[1]]$y == "Yes"]
  bg <- pred[testData_k562[[1]]$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  pgkPR[[i]] <- pr$curve[,-1]
  
  pgkAUPRC[[i]] <- m[[3]][20,2]
}

saveRDS(pgkPR, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/pgkPR.rds")
saveRDS(pgkAUPRC, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/pgkAUPRC.rds")

##################

chrs <- paste0("CHR", c(1:8,10:22))
pkkPR <- list()
pkkAUPRC <- list()
for(i in 1:length(chrs)){
  print(i)
  m <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/10kb/results_by_chr/",
                      chrs[i],
                      "/rus/distance/TADRF_holdout_k562_on_k562_peakachu.rds"))
  
  testData_k562  <- createTADdata(bounds.GR=bounds.GR_k562,
                                  resolution=5000,
                                  genomicElements.GR=genomicElements.GR_k562,
                                  featureType="distance",
                                  resampling="none",
                                  trainCHR=chrs[i],
                                  predictCHR=NULL,
                                  seed=123)
  
  pred <- as.vector(predict(m[[1]],newdata=testData_k562[[1]][,-1],type="prob")[,"Yes"])
  fg <- pred[testData_k562[[1]]$y == "Yes"]
  bg <- pred[testData_k562[[1]]$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  pkkPR[[i]] <- pr$curve[,-1]
  
  pkkAUPRC[[i]] <- m[[3]][20,2]
}

saveRDS(pkkPR, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/pkkPR.rds")
saveRDS(pkkAUPRC, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/pkkAUPRC.rds")

##################

chrs <- paste0("CHR", c(1:8,10:22))
pkgPR <- list()
pkgAUPRC <- list()
for(i in 1:length(chrs)){
  print(i)
  m <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/10kb/results_by_chr/",
                      chrs[i],
                      "/rus/distance/TADRF_holdout_k562_on_gm12878_peakachu.rds"))
  
  testData_gm12878  <- createTADdata(bounds.GR=bounds.GR_gm12878,
                                     resolution=5000,
                                     genomicElements.GR=genomicElements.GR_gm12878,
                                     featureType="distance",
                                     resampling="none",
                                     trainCHR=chrs[i],
                                     predictCHR=NULL,
                                     seed=123)
  
  pred <- as.vector(predict(m[[1]],newdata=testData_gm12878[[1]][,-1],type="prob")[,"Yes"])
  fg <- pred[testData_gm12878[[1]]$y == "Yes"]
  bg <- pred[testData_gm12878[[1]]$y == "No"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  pkgPR[[i]] <- pr$curve[,-1]
  
  pkgAUPRC[[i]] <- m[[3]][20,2]
}

saveRDS(pkgPR, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/pkgPR.rds")
saveRDS(pkgAUPRC, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/pkgAUPRC.rds")

###########################################################################################

#plotting

## gm12878

aggPR <- readRDS("Z:/TAD_data_analysis/miscellaneous/aggPR.rds")
aggAUPRC <- readRDS("Z:/TAD_data_analysis/miscellaneous/aggAUPRC.rds")

agkPR <- readRDS("Z:/TAD_data_analysis/miscellaneous/agkPR.rds")
agkAUPRC <- readRDS("Z:/TAD_data_analysis/miscellaneous/agkAUPRC.rds")

pggPR <- readRDS("Z:/TAD_data_analysis/miscellaneous/pggPR.rds")
pggAUPRC <- readRDS("Z:/TAD_data_analysis/miscellaneous/pggAUPRC.rds")

pgkPR <- readRDS("Z:/TAD_data_analysis/miscellaneous/pgkPR.rds")
pgkAUPRC <- readRDS("Z:/TAD_data_analysis/miscellaneous/pgkAUPRC.rds")


d <- rbind.data.frame(aggPR[[1]],agkPR[[1]],pggPR[[1]],pgkPR[[1]])
d$train_test <- c(rep("GM12878 on GM12878", nrow(aggPR[[1]])),
                  rep("GM12878 on K562", nrow(agkPR[[1]])),
                  rep("GM12878 on GM12878", nrow(pggPR[[1]])),
                  rep("GM12878 on K562", nrow(pgkPR[[1]])))
d$tool <- c(rep("Arrowhead", nrow(aggPR[[1]]) + nrow(agkPR[[1]])),
            rep("Peakachu", nrow(pggPR[[1]]) + nrow(pgkPR[[1]])))

ggplot(d, aes(x=1-V2, y=V1, linetype=train_test, color=tool)) +
  facet_grid(. ~ tool, scales = "free_y") +
  geom_line(size=1.5) +
  geom_text(
    label="AUPRC: ", 
    x=.67,
    y=.30,
    size = 4.5,
    color="black"
  ) +
  geom_rect(
    xmin = .50,
    xmax=.995,
    ymin=.27,
    ymax=.33,
    alpha=0,
    #color=Tool,
    fill="white"
  ) +
  geom_text(
    label="AUPRC: ", 
    x=.67,
    y=.20,
    size = 4.5,
    color="black"
  ) +
  geom_rect(
    xmin = .50,
    xmax=.995,
    ymin=.17,
    ymax=.23,
    alpha=0,
    #color=Tool,
    fill="white",
    linetype="dashed"
  ) +
  xlab("Recall") + 
  ylab("Precision") + 
  scale_linetype_discrete(name="Training/Testing") +
  scale_color_manual(values = c("blue","red")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  guides(color=FALSE, linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")

################################################################################

#plotting

## k562

akkPR <- readRDS("Z:/TAD_data_analysis/miscellaneous/akkPR.rds")
akkAUPRC <- readRDS("Z:/TAD_data_analysis/miscellaneous/akkAUPRC.rds")

akgPR <- readRDS("Z:/TAD_data_analysis/miscellaneous/akgPR.rds")
agkAUPRC <- readRDS("Z:/TAD_data_analysis/miscellaneous/akgAUPRC.rds")

pkkPR <- readRDS("Z:/TAD_data_analysis/miscellaneous/pkkPR.rds")
pkkAUPRC <- readRDS("Z:/TAD_data_analysis/miscellaneous/pkkAUPRC.rds")

pkgPR <- readRDS("Z:/TAD_data_analysis/miscellaneous/pkgPR.rds")
pkgAUPRC <- readRDS("Z:/TAD_data_analysis/miscellaneous/pkgAUPRC.rds")


d <- rbind.data.frame(akkPR[[1]],akgPR[[1]],pkkPR[[1]],pkgPR[[1]])
d$train_test <- c(rep("K562 on K562", nrow(akkPR[[1]])),
                  rep("K562 on GM12878", nrow(akgPR[[1]])),
                  rep("K562 on K562", nrow(pkkPR[[1]])),
                  rep("K562 on GM12878", nrow(pkgPR[[1]])))
d$train_test <- factor(d$train_test, levels=c("K562 on K562","K562 on GM12878"))
d$tool <- c(rep("Arrowhead", nrow(akkPR[[1]]) + nrow(akgPR[[1]])),
            rep("Peakachu", nrow(pkkPR[[1]]) + nrow(pkgPR[[1]])))

ggplot(d, aes(x=1-V2, y=V1, linetype=train_test, color=tool)) +
  facet_grid(. ~ tool, scales = "free_y") +
  geom_line(size=1.5) +
  geom_text(
    label="AUPRC: ", 
    x=.67,
    y=.25,
    size = 4.5,
    color="black"
  ) +
  geom_rect(
    xmin = .50,
    xmax=.995,
    ymin=.22,
    ymax=.28,
    alpha=0,
    #color=Tool,
    fill="white"
  ) +
  geom_text(
    label="AUPRC: ", 
    x=.67,
    y=.15,
    size = 4.5,
    color="black"
  ) +
  geom_rect(
    xmin = .50,
    xmax=.995,
    ymin=.12,
    ymax=.18,
    alpha=0,
    #color=Tool,
    fill="white",
    linetype="dashed"
  ) +
  xlab("Recall") + 
  ylab("Precision") + 
  scale_linetype_discrete(name="Training/Testing") +
  scale_color_manual(values = c("blue","red")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  guides(color=FALSE, linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme_minimal() +
  theme_bw()+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        #panel.spacing = unit(2, "lines"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.title = element_text(size=20),
        legend.position = "bottom")


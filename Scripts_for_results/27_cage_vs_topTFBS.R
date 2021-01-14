#27_cage_vs_topTFBS

library(dplyr)
library(ggplot2)

#arrowhead

##gm12878

#badata_gm12878_5kb <- numeric()
#for(i in paste0("CHR", c(1:8,10:22))){
# for(j in c("ctcf",
#             "rad21",
#             "smc3",
#             "znf143",
#             "ctcf_rad21",
#             "ctcf_smc3",
#             "ctcf_znf143",
#             "rad21_smc3",
#             "rad21_znf143",
#             "smc3_znf143",
#             "ctcf_rad21_smc3",
#             "ctcf_rad21_znf143",
#             "rad21_smc3_znf143",
#             "smc3_znf143_ctcf",
#             "topTFBS")){
#    print(c(i,j))
#    
#    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb/results_by_chr/",
#                               i, 
#                               "/",
#                               "rus", 
#                               "/", 
#                               "distance", 
#                               "/TADRF_holdout_",
#                               j,
#                              ".rds"))
#    perfdata <- perfdata[[3]]						   
#    
#    badata_gm12878_5kb <- rbind(badata_gm12878_5kb,
#                                perfdata[10,2])
#  }
#}


#saveRDS(badata_gm12878_5kb, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_gm12878_topTFBS_combinations.rds")

badata_gm12878_5kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_gm12878_topTFBS_combinations.rds")

all_rfperfs <- cbind.data.frame(V1=badata_gm12878_5kb,
                                Combination=rep(c("ctcf",
                                                  "rad21",
                                                  "smc3",
                                                  "znf143",
                                                  "ctcf_rad21",
                                                  "ctcf_smc3",
                                                  "ctcf_znf143",
                                                  "rad21_smc3",
                                                  "rad21_znf143",
                                                  "smc3_znf143",
                                                  "ctcf_rad21_smc3",
                                                  "ctcf_rad21_znf143",
                                                  "rad21_smc3_znf143",
                                                  "smc3_znf143_ctcf",
                                                  "topTFBS"),21),
                                Chromosome=c(rep("CHR1", 15),
                                             rep("CHR2", 15),
                                             rep("CHR3", 15),
                                             rep("CHR4", 15),
                                             rep("CHR5", 15),
                                             rep("CHR6", 15),
                                             rep("CHR7", 15),
                                             rep("CHR8", 15),
                                             rep("CHR10", 15),
                                             rep("CHR11", 15),
                                             rep("CHR12", 15),
                                             rep("CHR13", 15),
                                             rep("CHR14", 15),
                                             rep("CHR15", 15),
                                             rep("CHR16", 15),
                                             rep("CHR17", 15),
                                             rep("CHR18", 15),
                                             rep("CHR19", 15),
                                             rep("CHR20", 15),
                                             rep("CHR21", 15),
                                             rep("CHR22", 15)))

all_rfperfs$Chromosome <- factor(all_rfperfs$Chromosome, levels = paste0("CHR", c(1:8,10:22)))
all_rfperfs$Combination <- factor(all_rfperfs$Combination, levels = c("ctcf",
                                                                      "rad21",
                                                                      "smc3",
                                                                      "znf143",
                                                                      "ctcf_rad21",
                                                                      "ctcf_smc3",
                                                                      "ctcf_znf143",
                                                                      "rad21_smc3",
                                                                      "rad21_znf143",
                                                                      "smc3_znf143",
                                                                      "ctcf_rad21_smc3",
                                                                      "ctcf_rad21_znf143",
                                                                      "rad21_smc3_znf143",
                                                                      "smc3_znf143_ctcf",
                                                                      "topTFBS"))

all_rfperfs_ba <- all_rfperfs %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))

##############################################################################################

#badata_gm12878_5kb <- numeric()
#for(i in paste0("CHR", c(1:8,10:22))){
#  for(j in c("cage",
#             "cage_ctcf",
#             "cage_rad21",
#             "cage_smc3",
#             "cage_znf143",
#             "cage_ctcf_rad21",
#             "cage_ctcf_smc3",
#             "cage_ctcf_znf143",
#             "cage_rad21_smc3",
#             "cage_rad21_znf143",
#             "cage_smc3_znf143",
#             "cage_ctcf_rad21_smc3",
#             "cage_rad21_smc3_znf143",
#             "cage_smc3_znf143_ctcf",
#             "cage_znf143_ctcf_rad21",
#             "cage_ctcf_rad21_smc3_znf143")){
#    print(c(i,j))
#    
#    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/5kb/results_by_chr/",
#                              i, 
#                               "/",
#                               "rus", 
#                               "/", 
#                               "distance", 
#                               "/TADRF_holdout_",
#                               j,
#                               ".rds"))
#    perfdata <- perfdata[[3]]						   
#    
#    badata_gm12878_5kb <- rbind(badata_gm12878_5kb,
#                                perfdata[10,2])
#  }
#}

#saveRDS(badata_gm12878_5kb, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_gm12878_topTFBS_combinations_cage.rds")

badata_gm12878_5kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_gm12878_topTFBS_combinations_cage.rds")

all_rfperfs_cage <- cbind.data.frame(V1=badata_gm12878_5kb,
                                     Combination=rep(c("cage",
                                                       "cage_ctcf",
                                                       "cage_rad21",
                                                       "cage_smc3",
                                                       "cage_znf143",
                                                       "cage_ctcf_rad21",
                                                       "cage_ctcf_smc3",
                                                       "cage_ctcf_znf143",
                                                       "cage_rad21_smc3",
                                                       "cage_rad21_znf143",
                                                       "cage_smc3_znf143",
                                                       "cage_ctcf_rad21_smc3",
                                                       "cage_rad21_smc3_znf143",
                                                       "cage_smc3_znf143_ctcf",
                                                       "cage_znf143_ctcf_rad21",
                                                       "cage_ctcf_rad21_smc3_znf143"),21),
                                     Chromosome=c(rep("CHR1", 16),
                                                  rep("CHR2", 16),
                                                  rep("CHR3", 16),
                                                  rep("CHR4", 16),
                                                  rep("CHR5", 16),
                                                  rep("CHR6", 16),
                                                  rep("CHR7", 16),
                                                  rep("CHR8", 16),
                                                  rep("CHR10", 16),
                                                  rep("CHR11", 16),
                                                  rep("CHR12", 16),
                                                  rep("CHR13", 16),
                                                  rep("CHR14", 16),
                                                  rep("CHR15", 16),
                                                  rep("CHR16", 16),
                                                  rep("CHR17", 16),
                                                  rep("CHR18", 16),
                                                  rep("CHR19", 16),
                                                  rep("CHR20", 16),
                                                  rep("CHR21", 16),
                                                  rep("CHR22", 16)))

all_rfperfs_cage$Chromosome <- factor(all_rfperfs_cage$Chromosome, levels = paste0("CHR", c(1:8,10:22)))
all_rfperfs_cage$Combination <- factor(all_rfperfs_cage$Combination, levels = c("cage",
                                                                                "cage_ctcf",
                                                                                "cage_rad21",
                                                                                "cage_smc3",
                                                                                "cage_znf143",
                                                                                "cage_ctcf_rad21",
                                                                                "cage_ctcf_smc3",
                                                                                "cage_ctcf_znf143",
                                                                                "cage_rad21_smc3",
                                                                                "cage_rad21_znf143",
                                                                                "cage_smc3_znf143",
                                                                                "cage_ctcf_rad21_smc3",
                                                                                "cage_rad21_smc3_znf143",
                                                                                "cage_smc3_znf143_ctcf",
                                                                                "cage_znf143_ctcf_rad21",
                                                                                "cage_ctcf_rad21_smc3_znf143"))

all_rfperfs_cage_ba <- all_rfperfs_cage %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))


all_rfperfs_point <- rbind.data.frame(all_rfperfs_ba,all_rfperfs_cage_ba)
all_rfperfs_box <- rbind.data.frame(all_rfperfs,all_rfperfs_cage)


a2 <- ggplot(all_rfperfs_box, aes(x=reorder(Combination, V1, mean), y=V1, fill=Combination)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(.50,.85) +
  guides(fill=FALSE) +
  theme_minimal() +
  theme_bw()+
  xlab("Predictor Combination") + 
  ylab("Balanced Accuracy") + 
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  strip.text.x = element_text(size = 15),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.text=element_text(size=15),
  legend.title=element_text(size=15),
  legend.position = "bottom",
  plot.title = element_text(size=15))


################################################################################
################################################################################

#k562

#badata_k562_5kb <- numeric()
#for(i in paste0("CHR", c(1:8,10:22))){
#  for(j in c("ctcf",
#             "rad21",
#             "smc3",
#             "znf143",
#             "ctcf_rad21",
#             "ctcf_smc3",
#             "ctcf_znf143",
#             "rad21_smc3",
#             "rad21_znf143",
#             "smc3_znf143",
#             "ctcf_rad21_smc3",
#             "ctcf_rad21_znf143",
#             "rad21_smc3_znf143",
#             "smc3_znf143_ctcf",
#             "topTFBS")){
#    print(c(i,j))
#    
#    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/5kb/results_by_chr/",
#                               i, 
#                               "/",
#                               "rus", 
#                               "/", 
#                               "distance", 
#                               "/TADRF_holdout_",
#                               j,
#                               ".rds"))
#    perfdata <- perfdata[[3]]						   
#   
#    badata_k562_5kb <- rbind(badata_k562_5kb,
#                             perfdata[10,2])
#  }
#}


#saveRDS(badata_k562_5kb, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_k562_topTFBS_combinations.rds")

badata_k562_5kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_k562_topTFBS_combinations.rds")

all_rfperfs <- cbind.data.frame(V1=badata_k562_5kb,
                                Combination=rep(c("ctcf",
                                                  "rad21",
                                                  "smc3",
                                                  "znf143",
                                                  "ctcf_rad21",
                                                  "ctcf_smc3",
                                                  "ctcf_znf143",
                                                  "rad21_smc3",
                                                  "rad21_znf143",
                                                  "smc3_znf143",
                                                  "ctcf_rad21_smc3",
                                                  "ctcf_rad21_znf143",
                                                  "rad21_smc3_znf143",
                                                  "smc3_znf143_ctcf",
                                                  "topTFBS"),21),
                                Chromosome=c(rep("CHR1", 15),
                                             rep("CHR2", 15),
                                             rep("CHR3", 15),
                                             rep("CHR4", 15),
                                             rep("CHR5", 15),
                                             rep("CHR6", 15),
                                             rep("CHR7", 15),
                                             rep("CHR8", 15),
                                             rep("CHR10", 15),
                                             rep("CHR11", 15),
                                             rep("CHR12", 15),
                                             rep("CHR13", 15),
                                             rep("CHR14", 15),
                                             rep("CHR15", 15),
                                             rep("CHR16", 15),
                                             rep("CHR17", 15),
                                             rep("CHR18", 15),
                                             rep("CHR19", 15),
                                             rep("CHR20", 15),
                                             rep("CHR21", 15),
                                             rep("CHR22", 15)))

all_rfperfs$Chromosome <- factor(all_rfperfs$Chromosome, levels = paste0("CHR", c(1:8,10:22)))
all_rfperfs$Combination <- factor(all_rfperfs$Combination, levels = c("ctcf",
                                                                      "rad21",
                                                                      "smc3",
                                                                      "znf143",
                                                                      "ctcf_rad21",
                                                                      "ctcf_smc3",
                                                                      "ctcf_znf143",
                                                                      "rad21_smc3",
                                                                      "rad21_znf143",
                                                                      "smc3_znf143",
                                                                      "ctcf_rad21_smc3",
                                                                      "ctcf_rad21_znf143",
                                                                      "rad21_smc3_znf143",
                                                                      "smc3_znf143_ctcf",
                                                                      "topTFBS"))

all_rfperfs_ba <- all_rfperfs %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))

##############################################################################################

#badata_k562_5kb <- numeric()
#for(i in paste0("CHR", c(1:8,10:22))){
#  for(j in c("cage",
#             "cage_ctcf",
#             "cage_rad21",
#             "cage_smc3",
#             "cage_znf143",
#             "cage_ctcf_rad21",
#             "cage_ctcf_smc3",
#            "cage_ctcf_znf143",
#            "cage_rad21_smc3",
#             "cage_rad21_znf143",
#             "cage_smc3_znf143",
#             "cage_ctcf_rad21_smc3",
#             "cage_rad21_smc3_znf143",
#             "cage_smc3_znf143_ctcf",
#             "cage_znf143_ctcf_rad21",
#             "cage_ctcf_rad21_smc3_znf143")){
#    print(c(i,j))
#    
#    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/5kb/results_by_chr/",
#                               i, 
#                               "/",
#                               "rus", 
#                               "/", 
#                               "distance", 
#                               "/TADRF_holdout_",
#                               j,
#                               ".rds"))
#    perfdata <- perfdata[[3]]						   
#    
#    badata_k562_5kb <- rbind(badata_k562_5kb,
#                             perfdata[10,2])
#  }
#}

#saveRDS(badata_k562_5kb, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_k562_topTFBS_combinations_cage.rds")

badata_k562_5kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_k562_topTFBS_combinations_cage.rds")

all_rfperfs_cage <- cbind.data.frame(V1=badata_k562_5kb,
                                     Combination=rep(c("cage",
                                                       "cage_ctcf",
                                                       "cage_rad21",
                                                       "cage_smc3",
                                                       "cage_znf143",
                                                       "cage_ctcf_rad21",
                                                       "cage_ctcf_smc3",
                                                       "cage_ctcf_znf143",
                                                       "cage_rad21_smc3",
                                                       "cage_rad21_znf143",
                                                       "cage_smc3_znf143",
                                                       "cage_ctcf_rad21_smc3",
                                                       "cage_rad21_smc3_znf143",
                                                       "cage_smc3_znf143_ctcf",
                                                       "cage_znf143_ctcf_rad21",
                                                       "cage_ctcf_rad21_smc3_znf143"),21),
                                     Chromosome=c(rep("CHR1", 16),
                                                  rep("CHR2", 16),
                                                  rep("CHR3", 16),
                                                  rep("CHR4", 16),
                                                  rep("CHR5", 16),
                                                  rep("CHR6", 16),
                                                  rep("CHR7", 16),
                                                  rep("CHR8", 16),
                                                  rep("CHR10", 16),
                                                  rep("CHR11", 16),
                                                  rep("CHR12", 16),
                                                  rep("CHR13", 16),
                                                  rep("CHR14", 16),
                                                  rep("CHR15", 16),
                                                  rep("CHR16", 16),
                                                  rep("CHR17", 16),
                                                  rep("CHR18", 16),
                                                  rep("CHR19", 16),
                                                  rep("CHR20", 16),
                                                  rep("CHR21", 16),
                                                  rep("CHR22", 16)))

all_rfperfs_cage$Chromosome <- factor(all_rfperfs_cage$Chromosome, levels = paste0("CHR", c(1:8,10:22)))
all_rfperfs_cage$Combination <- factor(all_rfperfs_cage$Combination, levels = c("cage",
                                                                                "cage_ctcf",
                                                                                "cage_rad21",
                                                                                "cage_smc3",
                                                                                "cage_znf143",
                                                                                "cage_ctcf_rad21",
                                                                                "cage_ctcf_smc3",
                                                                                "cage_ctcf_znf143",
                                                                                "cage_rad21_smc3",
                                                                                "cage_rad21_znf143",
                                                                                "cage_smc3_znf143",
                                                                                "cage_ctcf_rad21_smc3",
                                                                                "cage_rad21_smc3_znf143",
                                                                                "cage_smc3_znf143_ctcf",
                                                                                "cage_znf143_ctcf_rad21",
                                                                                "cage_ctcf_rad21_smc3_znf143"))

all_rfperfs_cage_ba <- all_rfperfs_cage %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))


all_rfperfs_point <- rbind.data.frame(all_rfperfs_ba,all_rfperfs_cage_ba)
all_rfperfs_box <- rbind.data.frame(all_rfperfs,all_rfperfs_cage)


b2 <- ggplot(all_rfperfs_box, aes(x=reorder(Combination, V1, mean), y=V1, fill=Combination)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(.50,.85) +
  guides(fill=FALSE) +
  theme_minimal() +
  theme_bw()+
  xlab("Predictor Combination") + 
  ylab("Balanced Accuracy") + 
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  strip.text.x = element_text(size = 15),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.text=element_text(size=15),
  legend.title=element_text(size=15),
  legend.position = "bottom",
  plot.title = element_text(size=15))


################################################################################
################################################################################
################################################################################
################################################################################

#peakachu

##gm12878

#badata_gm12878_10kb <- numeric()
#for(i in paste0("CHR", c(1:8,10:22))){
#  for(j in c("ctcf",
#             "rad21",
#             "smc3",
#             "znf143",
#             "ctcf_rad21",
#             "ctcf_smc3",
#             "ctcf_znf143",
#             "rad21_smc3",
#             "rad21_znf143",
#             "smc3_znf143",
#             "ctcf_rad21_smc3",
#             "ctcf_rad21_znf143",
#             "rad21_smc3_znf143",
#             "smc3_znf143_ctcf",
#             "topTFBS")){
#    print(c(i,j))
#    
#    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/10kb/results_by_chr/",
#                               i, 
#                               "/",
#                               "rus", 
#                               "/", 
#                               "distance", 
#                               "/TADRF_holdout_peakachu_",
#                               j,
#                               ".rds"))
#    perfdata <- perfdata[[3]]						   
#    
#    badata_gm12878_10kb <- rbind(badata_gm12878_10kb,
#                                 perfdata[10,2])
#  }
#}


#saveRDS(badata_gm12878_10kb, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_gm12878_topTFBS_combinations_peakachu.rds")

badata_gm12878_10kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_gm12878_topTFBS_combinations_peakachu.rds")

all_rfperfs <- cbind.data.frame(V1=badata_gm12878_10kb,
                                Combination=rep(c("ctcf",
                                                  "rad21",
                                                  "smc3",
                                                  "znf143",
                                                  "ctcf_rad21",
                                                  "ctcf_smc3",
                                                  "ctcf_znf143",
                                                  "rad21_smc3",
                                                  "rad21_znf143",
                                                  "smc3_znf143",
                                                  "ctcf_rad21_smc3",
                                                  "ctcf_rad21_znf143",
                                                  "rad21_smc3_znf143",
                                                  "smc3_znf143_ctcf",
                                                  "topTFBS"),21),
                                Chromosome=c(rep("CHR1", 15),
                                             rep("CHR2", 15),
                                             rep("CHR3", 15),
                                             rep("CHR4", 15),
                                             rep("CHR5", 15),
                                             rep("CHR6", 15),
                                             rep("CHR7", 15),
                                             rep("CHR8", 15),
                                             rep("CHR10", 15),
                                             rep("CHR11", 15),
                                             rep("CHR12", 15),
                                             rep("CHR13", 15),
                                             rep("CHR14", 15),
                                             rep("CHR15", 15),
                                             rep("CHR16", 15),
                                             rep("CHR17", 15),
                                             rep("CHR18", 15),
                                             rep("CHR19", 15),
                                             rep("CHR20", 15),
                                             rep("CHR21", 15),
                                             rep("CHR22", 15)))

all_rfperfs$Chromosome <- factor(all_rfperfs$Chromosome, levels = paste0("CHR", c(1:8,10:22)))
all_rfperfs$Combination <- factor(all_rfperfs$Combination, levels = c("ctcf",
                                                                      "rad21",
                                                                      "smc3",
                                                                      "znf143",
                                                                      "ctcf_rad21",
                                                                      "ctcf_smc3",
                                                                      "ctcf_znf143",
                                                                      "rad21_smc3",
                                                                      "rad21_znf143",
                                                                      "smc3_znf143",
                                                                      "ctcf_rad21_smc3",
                                                                      "ctcf_rad21_znf143",
                                                                      "rad21_smc3_znf143",
                                                                      "smc3_znf143_ctcf",
                                                                      "topTFBS"))

all_rfperfs_ba <- all_rfperfs %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))

##############################################################################################

#badata_gm12878_10kb <- numeric()
#for(i in paste0("CHR", c(1:8,10:22))){
#  for(j in c("cage",
#             "cage_ctcf",
#             "cage_rad21",
#             "cage_smc3",
#             "cage_znf143",
#            "cage_ctcf_rad21",
#            "cage_ctcf_smc3",
#             "cage_ctcf_znf143",
#             "cage_rad21_smc3",
#             "cage_rad21_znf143",
#             "cage_smc3_znf143",
#             "cage_ctcf_rad21_smc3",
#             "cage_rad21_smc3_znf143",
#             "cage_smc3_znf143_ctcf",
#             "cage_znf143_ctcf_rad21",
#             "cage_ctcf_rad21_smc3_znf143")){
#   print(c(i,j))
#   
#   perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/GM12878/10kb/results_by_chr/",
#                               i, 
#                               "/",
#                               "rus", 
#                               "/", 
#                               "distance", 
#                              "/TADRF_holdout_peakachu_",
#                               j,
#                               ".rds"))
#    perfdata <- perfdata[[3]]						   
#    
#    badata_gm12878_10kb <- rbind(badata_gm12878_10kb,
#                                 perfdata[10,2])
#  }
#}

#saveRDS(badata_gm12878_10kb, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_gm12878_topTFBS_combinations_peakachu_cage.rds")

badata_gm12878_10kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_gm12878_topTFBS_combinations_peakachu_cage.rds")

all_rfperfs_cage <- cbind.data.frame(V1=badata_gm12878_10kb,
                                     Combination=rep(c("cage",
                                                       "cage_ctcf",
                                                       "cage_rad21",
                                                       "cage_smc3",
                                                       "cage_znf143",
                                                       "cage_ctcf_rad21",
                                                       "cage_ctcf_smc3",
                                                       "cage_ctcf_znf143",
                                                       "cage_rad21_smc3",
                                                       "cage_rad21_znf143",
                                                       "cage_smc3_znf143",
                                                       "cage_ctcf_rad21_smc3",
                                                       "cage_rad21_smc3_znf143",
                                                       "cage_smc3_znf143_ctcf",
                                                       "cage_znf143_ctcf_rad21",
                                                       "cage_ctcf_rad21_smc3_znf143"),21),
                                     Chromosome=c(rep("CHR1", 16),
                                                  rep("CHR2", 16),
                                                  rep("CHR3", 16),
                                                  rep("CHR4", 16),
                                                  rep("CHR5", 16),
                                                  rep("CHR6", 16),
                                                  rep("CHR7", 16),
                                                  rep("CHR8", 16),
                                                  rep("CHR10", 16),
                                                  rep("CHR11", 16),
                                                  rep("CHR12", 16),
                                                  rep("CHR13", 16),
                                                  rep("CHR14", 16),
                                                  rep("CHR15", 16),
                                                  rep("CHR16", 16),
                                                  rep("CHR17", 16),
                                                  rep("CHR18", 16),
                                                  rep("CHR19", 16),
                                                  rep("CHR20", 16),
                                                  rep("CHR21", 16),
                                                  rep("CHR22", 16)))

all_rfperfs_cage$Chromosome <- factor(all_rfperfs_cage$Chromosome, levels = paste0("CHR", c(1:8,10:22)))
all_rfperfs_cage$Combination <- factor(all_rfperfs_cage$Combination, levels = c("cage",
                                                                                "cage_ctcf",
                                                                                "cage_rad21",
                                                                                "cage_smc3",
                                                                                "cage_znf143",
                                                                                "cage_ctcf_rad21",
                                                                                "cage_ctcf_smc3",
                                                                                "cage_ctcf_znf143",
                                                                                "cage_rad21_smc3",
                                                                                "cage_rad21_znf143",
                                                                                "cage_smc3_znf143",
                                                                                "cage_ctcf_rad21_smc3",
                                                                                "cage_rad21_smc3_znf143",
                                                                                "cage_smc3_znf143_ctcf",
                                                                                "cage_znf143_ctcf_rad21",
                                                                                "cage_ctcf_rad21_smc3_znf143"))

all_rfperfs_cage_ba <- all_rfperfs_cage %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))


all_rfperfs_point <- rbind.data.frame(all_rfperfs_ba,all_rfperfs_cage_ba)
all_rfperfs_box <- rbind.data.frame(all_rfperfs,all_rfperfs_cage)



a4 <- ggplot(all_rfperfs_box, aes(x=reorder(Combination, V1, mean), y=V1, fill=Combination)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(.50,.85) +
  guides(fill=FALSE) +
  theme_minimal() +
  theme_bw()+
  xlab("Predictor Combination") + 
  ylab("Balanced Accuracy") + 
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  strip.text.x = element_text(size = 15),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.text=element_text(size=15),
  legend.title=element_text(size=15),
  legend.position = "bottom",
  plot.title = element_text(size=15))


################################################################################
################################################################################

#k562

#badata_k562_10kb <- numeric()
#for(i in paste0("CHR", c(1:8,10:22))){
#  for(j in c("ctcf",
#             "rad21",
#             "smc3",
#             "znf143",
#             "ctcf_rad21",
#             "ctcf_smc3",
#             "ctcf_znf143",
#             "rad21_smc3",
#             "rad21_znf143",
#             "smc3_znf143",
#             "ctcf_rad21_smc3",
#             "ctcf_rad21_znf143",
#             "rad21_smc3_znf143",
#             "smc3_znf143_ctcf",
#             "topTFBS")){
#    print(c(i,j))
#   
#    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/10kb/results_by_chr/",
#                               i, 
#                               "/",
#                               "rus", 
#                               "/", 
#                               "distance", 
#                               "/TADRF_holdout_peakachu_",
#                               j,
#                               ".rds"))
#    perfdata <- perfdata[[3]]						   
#    
#    badata_k562_10kb <- rbind(badata_k562_10kb,
#                              perfdata[10,2])
#  }
#}


#saveRDS(badata_k562_10kb, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_k562_topTFBS_combinations_peakachu.rds")

badata_k562_10kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_k562_topTFBS_combinations_peakachu.rds")

all_rfperfs <- cbind.data.frame(V1=badata_k562_10kb,
                                Combination=rep(c("ctcf",
                                                  "rad21",
                                                  "smc3",
                                                  "znf143",
                                                  "ctcf_rad21",
                                                  "ctcf_smc3",
                                                  "ctcf_znf143",
                                                  "rad21_smc3",
                                                  "rad21_znf143",
                                                  "smc3_znf143",
                                                  "ctcf_rad21_smc3",
                                                  "ctcf_rad21_znf143",
                                                  "rad21_smc3_znf143",
                                                  "smc3_znf143_ctcf",
                                                  "topTFBS"),21),
                                Chromosome=c(rep("CHR1", 15),
                                             rep("CHR2", 15),
                                             rep("CHR3", 15),
                                             rep("CHR4", 15),
                                             rep("CHR5", 15),
                                             rep("CHR6", 15),
                                             rep("CHR7", 15),
                                             rep("CHR8", 15),
                                             rep("CHR10", 15),
                                             rep("CHR11", 15),
                                             rep("CHR12", 15),
                                             rep("CHR13", 15),
                                             rep("CHR14", 15),
                                             rep("CHR15", 15),
                                             rep("CHR16", 15),
                                             rep("CHR17", 15),
                                             rep("CHR18", 15),
                                             rep("CHR19", 15),
                                             rep("CHR20", 15),
                                             rep("CHR21", 15),
                                             rep("CHR22", 15)))

all_rfperfs$Chromosome <- factor(all_rfperfs$Chromosome, levels = paste0("CHR", c(1:8,10:22)))
all_rfperfs$Combination <- factor(all_rfperfs$Combination, levels = c("ctcf",
                                                                      "rad21",
                                                                      "smc3",
                                                                      "znf143",
                                                                      "ctcf_rad21",
                                                                      "ctcf_smc3",
                                                                      "ctcf_znf143",
                                                                      "rad21_smc3",
                                                                      "rad21_znf143",
                                                                      "smc3_znf143",
                                                                      "ctcf_rad21_smc3",
                                                                      "ctcf_rad21_znf143",
                                                                      "rad21_smc3_znf143",
                                                                      "smc3_znf143_ctcf",
                                                                      "topTFBS"))

all_rfperfs_ba <- all_rfperfs %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))

##############################################################################################

#badata_k562_10kb <- numeric()
#for(i in paste0("CHR", c(1:8,10:22))){
#  for(j in c("cage",
#             "cage_ctcf",
#             "cage_rad21",
#             "cage_smc3",
#             "cage_znf143",
#             "cage_ctcf_rad21",
#             "cage_ctcf_smc3",
#             "cage_ctcf_znf143",
#             "cage_rad21_smc3",
#             "cage_rad21_znf143",
#             "cage_smc3_znf143",
#             "cage_ctcf_rad21_smc3",
#             "cage_rad21_smc3_znf143",
#             "cage_smc3_znf143_ctcf",
#             "cage_znf143_ctcf_rad21",
#             "cage_ctcf_rad21_smc3_znf143")){
#   print(c(i,j))
#    
#    perfdata <- readRDS(paste0("/home/stilianoudakisc/TAD_data_analysis/K562/10kb/results_by_chr/",
#                               i, 
#                               "/",
#                               "rus", 
#                               "/", 
#                               "distance", 
#                               "/TADRF_holdout_peakachu_",
#                               j,
#                               ".rds"))
#    perfdata <- perfdata[[3]]						   
#    
#    badata_k562_10kb <- rbind(badata_k562_10kb,
#                              perfdata[10,2])
#  }
#}

#saveRDS(badata_k562_10kb, "/home/stilianoudakisc/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_k562_topTFBS_combinations_peakachu_cage.rds")

badata_k562_10kb <- readRDS("Z:/TAD_data_analysis/miscellaneous/saved_rds_files/TADRF_holdout_k562_topTFBS_combinations_peakachu_cage.rds")

all_rfperfs_cage <- cbind.data.frame(V1=badata_k562_10kb,
                                     Combination=rep(c("cage",
                                                       "cage_ctcf",
                                                       "cage_rad21",
                                                       "cage_smc3",
                                                       "cage_znf143",
                                                       "cage_ctcf_rad21",
                                                       "cage_ctcf_smc3",
                                                       "cage_ctcf_znf143",
                                                       "cage_rad21_smc3",
                                                       "cage_rad21_znf143",
                                                       "cage_smc3_znf143",
                                                       "cage_ctcf_rad21_smc3",
                                                       "cage_rad21_smc3_znf143",
                                                       "cage_smc3_znf143_ctcf",
                                                       "cage_znf143_ctcf_rad21",
                                                       "cage_ctcf_rad21_smc3_znf143"),21),
                                     Chromosome=c(rep("CHR1", 16),
                                                  rep("CHR2", 16),
                                                  rep("CHR3", 16),
                                                  rep("CHR4", 16),
                                                  rep("CHR5", 16),
                                                  rep("CHR6", 16),
                                                  rep("CHR7", 16),
                                                  rep("CHR8", 16),
                                                  rep("CHR10", 16),
                                                  rep("CHR11", 16),
                                                  rep("CHR12", 16),
                                                  rep("CHR13", 16),
                                                  rep("CHR14", 16),
                                                  rep("CHR15", 16),
                                                  rep("CHR16", 16),
                                                  rep("CHR17", 16),
                                                  rep("CHR18", 16),
                                                  rep("CHR19", 16),
                                                  rep("CHR20", 16),
                                                  rep("CHR21", 16),
                                                  rep("CHR22", 16)))

all_rfperfs_cage$Chromosome <- factor(all_rfperfs_cage$Chromosome, levels = paste0("CHR", c(1:8,10:22)))
all_rfperfs_cage$Combination <- factor(all_rfperfs_cage$Combination, levels = c("cage",
                                                                                "cage_ctcf",
                                                                                "cage_rad21",
                                                                                "cage_smc3",
                                                                                "cage_znf143",
                                                                                "cage_ctcf_rad21",
                                                                                "cage_ctcf_smc3",
                                                                                "cage_ctcf_znf143",
                                                                                "cage_rad21_smc3",
                                                                                "cage_rad21_znf143",
                                                                                "cage_smc3_znf143",
                                                                                "cage_ctcf_rad21_smc3",
                                                                                "cage_rad21_smc3_znf143",
                                                                                "cage_smc3_znf143_ctcf",
                                                                                "cage_znf143_ctcf_rad21",
                                                                                "cage_ctcf_rad21_smc3_znf143"))

all_rfperfs_cage_ba <- all_rfperfs_cage %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(Performance = mean(V1, na.rm = TRUE),
                   PerformanceSD = sd(V1, na.rm = TRUE))


all_rfperfs_point <- rbind.data.frame(all_rfperfs_ba,all_rfperfs_cage_ba)
all_rfperfs_box <- rbind.data.frame(all_rfperfs,all_rfperfs_cage)


b4 <- ggplot(all_rfperfs_box, aes(x=reorder(Combination, V1, mean), y=V1, fill=Combination)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(.50,.85) +
  guides(fill=FALSE) +
  theme_minimal() +
  theme_bw()+
  xlab("Predictor Combination") + 
  ylab("Balanced Accuracy") + 
  theme(axis.text.x = element_text(size=15,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  strip.text.x = element_text(size = 15),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.text=element_text(size=15),
  legend.title=element_text(size=15),
  legend.position = "bottom",
  plot.title = element_text(size=15))

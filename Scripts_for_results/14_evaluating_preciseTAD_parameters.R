#evaluating_preciseTAD_parameters

library(S4Vectors)
library(GenomicRanges)
library(ggplot2)
library(caret)
library(randomForest)
library(glmnet)
library(gbm)
library(DMwR)
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(cluster)
library(pROC)
library(SuperLearner)
library(data.table)
library(PRROC)
library(bigmemory)
library(ranger)
library(dplyr)
library(parallel)
library(doSNOW)
library(foreach)
library(pbapply)
library(ggpubr)
library(scales)
library(GGally)
library(VennDiagram)
library(network)
library(sna)
library(gtable)
library(grid)
library(magrittr) 
library(DT)
library(tidyverse)
library(reshape)
library(ROSE)
library(lattice)
library(corrplot)
library(cluster)
library(RColorBrewer)
library(GenometriCorr)
library(ggsignif)
library(Vennerable)
library(VennDiagram)
library(ChIPpeakAnno)
library(EnrichedHeatmap)

source("Z:/TAD_data_analysis/functions_for_R_package/preciseTAD.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADRF.R")
source("Z:/TAD_data_analysis/functions_for_R_package/TADrfe.R")
source("Z:/TAD_data_analysis/functions_for_R_package/distance_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/percent_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/count_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/binary_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/signal_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/annots_to_granges_func.R")
source("Z:/TAD_data_analysis/functions_for_R_package/extract_boundaries_func.R")

# arrowhead

pt.975.500_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.500.rds")
pt.975.1000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.1000.rds")
pt.975.2500_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.2500.rds")
pt.975.5000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.5000.rds")
pt.975.10000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.10000.rds")
pt.975.15000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.15000.rds")
pt.975.20000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.20000.rds")
pt.975.25000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.25000.rds")
#pt.975.30000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.975.30000.rds")

pt.99.500_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.500.rds")
pt.99.1000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.1000.rds")
pt.99.2500_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.2500.rds")
pt.99.5000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.5000.rds")
pt.99.10000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.10000.rds")
pt.99.15000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.15000.rds")
pt.99.20000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.20000.rds")
pt.99.25000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.25000.rds")
#pt.99.30000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt.99.30000.rds")

pt1.0.500_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.500.rds")
pt1.0.1000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.1000.rds")
pt1.0.2500_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.2500.rds")
pt1.0.5000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.5000.rds")
pt1.0.10000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.10000.rds")
pt1.0.15000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.15000.rds")
pt1.0.20000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.20000.rds")
pt1.0.25000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.25000.rds")
#pt1.0.30000_a <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/arrowhead_gt/pt1.0.30000.rds")

# peakachu

pt.975.500_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.500.rds")
pt.975.1000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.1000.rds")
pt.975.2500_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.2500.rds")
pt.975.5000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.5000.rds")
pt.975.10000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.10000.rds")
pt.975.15000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.15000.rds")
pt.975.20000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.20000.rds")
pt.975.25000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.25000.rds")
#pt.975.30000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.975.30000.rds")

pt.99.500_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.500.rds")
pt.99.1000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.1000.rds")
pt.99.2500_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.2500.rds")
pt.99.5000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.5000.rds")
pt.99.10000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.10000.rds")
pt.99.15000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.15000.rds")
pt.99.20000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.20000.rds")
pt.99.25000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.25000.rds")
#pt.99.30000 <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt.99.30000.rds")

pt1.0.500_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.500.rds")
pt1.0.1000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.1000.rds")
pt1.0.2500_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.2500.rds")
pt1.0.5000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.5000.rds")
pt1.0.10000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.10000.rds")
pt1.0.15000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.15000.rds")
pt1.0.20000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.20000.rds")
pt1.0.25000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.25000.rds")
#pt1.0.30000_p <- readRDS("Z:/TAD_data_analysis/miscellaneous/evaluating_preciseTAD_parameters/chr22/peakachu_gt/pt1.0.30000.rds")


#######################################################################################################################################

# plotting

## arrowhead

df <- data.frame(NE=c(pt.975.500_a$NormilizedEnrichment,
                      pt.975.1000_a$NormilizedEnrichment,
                      pt.975.2500_a$NormilizedEnrichment,
                      pt.975.5000_a$NormilizedEnrichment,
                      pt.975.10000_a$NormilizedEnrichment,
                      pt.975.15000_a$NormilizedEnrichment,
                      pt.975.20000_a$NormilizedEnrichment,
                      pt.975.25000_a$NormilizedEnrichment,
                      pt.99.500_a$NormilizedEnrichment,
                      pt.99.1000_a$NormilizedEnrichment,
                      pt.99.2500_a$NormilizedEnrichment,
                      pt.99.5000_a$NormilizedEnrichment,
                      pt.99.10000_a$NormilizedEnrichment,
                      pt.99.15000_a$NormilizedEnrichment,
                      pt.99.20000_a$NormilizedEnrichment,
                      pt.99.25000_a$NormilizedEnrichment,
                      pt1.0.500_a$NormilizedEnrichment,
                      pt1.0.1000_a$NormilizedEnrichment,
                      pt1.0.2500_a$NormilizedEnrichment,
                      pt1.0.5000_a$NormilizedEnrichment,
                      pt1.0.10000_a$NormilizedEnrichment,
                      pt1.0.15000_a$NormilizedEnrichment,
                      pt1.0.20000_a$NormilizedEnrichment,
                      pt1.0.25000_a$NormilizedEnrichment),
                 t=c(rep('t=0.975', 4*8),
                     rep('t=0.99', 4*8),
                     rep('t=1.0', 4*8)),
                 eps=rep(c(rep(500, 4),
                       rep(1000, 4),
                       rep(2500, 4),
                       rep(5000, 4),
                       rep(10000, 4),
                       rep(15000, 4),
                       rep(20000, 4),
                       rep(25000, 4)),3),
                 a=rep(c("CTCF","RAD21","SMC3","ZNF143"),8*3))

p1 <- ggplot(df, aes(x=eps, y=NE, group=a, color=a)) +
  facet_grid(. ~ t) +
  geom_line(size=2, position=position_dodge(0.5)) +
  geom_point(size=5, position=position_dodge(0.5)) +
  #geom_vline(xintercept = 20000, linetype="dashed") +
  theme_minimal() +
  theme_bw()+
  xlab("Epsilon") +
  ylab("Normalized Enrichment") +
  scale_color_discrete(name="Annotation") +
  theme(axis.text.x = element_text(size=20,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 20),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  strip.text.x = element_text(size = 20),
  legend.text=element_text(size=20),
  legend.title=element_text(size=20),
  plot.title = element_text(size=20),
  legend.position="bottom")

## peakachu

df <- data.frame(NE=c(pt.975.500_p$NormilizedEnrichment,
                      pt.975.1000_p$NormilizedEnrichment,
                      pt.975.2500_p$NormilizedEnrichment,
                      pt.975.5000_p$NormilizedEnrichment,
                      pt.975.10000_p$NormilizedEnrichment,
                      pt.975.15000_p$NormilizedEnrichment,
                      pt.975.20000_p$NormilizedEnrichment,
                      pt.975.25000_p$NormilizedEnrichment,
                      pt.99.500_p$NormilizedEnrichment,
                      pt.99.1000_p$NormilizedEnrichment,
                      pt.99.2500_p$NormilizedEnrichment,
                      pt.99.5000_p$NormilizedEnrichment,
                      pt.99.10000_p$NormilizedEnrichment,
                      pt.99.15000_p$NormilizedEnrichment,
                      pt.99.20000_p$NormilizedEnrichment,
                      pt.99.25000_p$NormilizedEnrichment,
                      pt1.0.500_p$NormilizedEnrichment,
                      pt1.0.1000_p$NormilizedEnrichment,
                      pt1.0.2500_p$NormilizedEnrichment,
                      pt1.0.5000_p$NormilizedEnrichment,
                      pt1.0.10000_p$NormilizedEnrichment,
                      pt1.0.15000_p$NormilizedEnrichment,
                      pt1.0.20000_p$NormilizedEnrichment,
                      pt1.0.25000_p$NormilizedEnrichment),
                 t=c(rep('t=0.975', 4*8),
                     rep('t=0.99', 4*8),
                     rep('t=1.0', 4*8)),
                 eps=rep(c(rep(500, 4),
                           rep(1000, 4),
                           rep(2500, 4),
                           rep(5000, 4),
                           rep(10000, 4),
                           rep(15000, 4),
                           rep(20000, 4),
                           rep(25000, 4)),3),
                 a=rep(c("CTCF","RAD21","SMC3","ZNF143"),8*3))


p2 <- ggplot(df, aes(x=eps, y=NE, group=a, color=a)) +
  facet_grid(. ~ t) +
  geom_line(size=2, position=position_dodge(0.5)) +
  geom_point(size=5, position=position_dodge(0.5)) +
  #geom_vline(xintercept = 20000, linetype="dashed") +
  theme_minimal() +
  theme_bw()+
  xlab("Epsilon") +
  ylab("Normalized Enrichment") +
  scale_color_discrete(name="Annotation") +
  theme(axis.text.x = element_text(size=20,
                                   angle = 45, 
                                   #margin = ggplot2::margin(t = 35),
                                   hjust = 1
  ),
  axis.text.y = element_text(size = 20),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  strip.text.x = element_text(size = 20),
  legend.text=element_text(size=20),
  legend.title=element_text(size=20),
  plot.title = element_text(size=20),
  legend.position="bottom")

ggarrange(p1,p2,ncol = 2, labels = "AUTO", common.legend = TRUE, legend = "bottom")



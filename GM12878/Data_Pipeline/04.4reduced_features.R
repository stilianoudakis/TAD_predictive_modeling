library(GenomicRanges)
library(caret)
library(data.table)
library(gbm)
library(randomForest)
library(glmnet)
library(pROC)
library(plyr)
library(dplyr)
library(ggplot2)
library(DMwR)
library(gridExtra)
library(pROC)
library(ROCR)
library(leaps)

setwd("/home/stilianoudakisc/TAD_data_analysis/evaluating_variable_reduction/")

gm12878_5kb_fwd <- readRDS("gm12878_5kb_fwd.rds")

gm12878_5kb_bwd <- readRDS("gm12878_5kb_bwd.rds")

gm12878_5kb_both <- readRDS("gm12878_5kb_both.rds")

fwdfeats <- names(gm12878_5kb_fwd)[-1]
bwdfeats <- names(gm12878_5kb_bwd)[-1]
bothfeats <- names(gm12878_5kb_both)[-1]
unionfeats <- union(union(fwdfeats, bwdfeats), bothfeats)

unionfeats <- c("y", unionfeats)

gm12878_5kb_f <- readRDS("gm12878_5kb_f.rds")

gm12878_5kb_reduced <- gm12878_5kb_f[, names(gm12878_5kb_f) %in% unionfeats]

dim(gm12878_5kb_reduced)

saveRDS(gm12878_5kb_reduced, "gm12878_5kb_reduced.rds")
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

setwd("/home/stilianoudakisc/TAD_data_analysis/K562/evaluating_variable_reduction/")

k562_10kb_fwd <- readRDS("k562_10kb_fwd.rds")

k562_10kb_bwd <- readRDS("k562_10kb_bwd.rds")

k562_10kb_both <- readRDS("k562_10kb_both.rds")

fwdfeats <- names(k562_10kb_fwd)[-1]
bwdfeats <- names(k562_10kb_bwd)[-1]
bothfeats <- names(k562_10kb_both)[-1]
unionfeats <- union(union(fwdfeats, bwdfeats), bothfeats)

unionfeats <- c("y", unionfeats)

k562_10kb_f <- readRDS("k562_10kb_f.rds")

k562_10kb_reduced <- k562_10kb_f[, names(k562_10kb_f) %in% unionfeats]

dim(k562_10kb_reduced)

saveRDS(k562_10kb_reduced, "k562_10kb_reduced.rds")
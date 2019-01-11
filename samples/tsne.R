#!/usr/bin/env Rscript

library(tidyverse)
library(glmnet)
library(Rtsne)

dir_work <- '/mnt/HA/groups/rosenGrp/embed_ag_samples/out'

seed <- 3245
args <- commandArgs(trailingOnly=TRUE)

model <- args[1]
nc <- as.integer(args[2])

path_out <- gsub('lasso','tsne',model)
cat(sprintf('Output will be saved to %s\n.',path_out))

lasso <- readRDS(model)

tsne_sample <- Rtsne(lasso$dat %>% select(-body_site,-SampleID) %>% as.matrix(),dims=2,
                     perplexity=250,check_duplicates=FALSE,theta=.5,max_iter=1000,
                     pca=FALSE,initial_dims=50,pca_center=TRUE,pca_scale=TRUE,  
                     verbose=TRUE)$Y
colnames(tsne_sample) <- c('D1','D2')

lasso$dat <- lasso$dat %>%
    left_join(as_tibble(tsne_sample) %>% mutate(SampleID=lasso$dat$SampleID),by='SampleID')

saveRDS(lasso,path_out)


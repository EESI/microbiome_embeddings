#!/usr/bin/env Rscript

.libPaths('/mnt/HA/groups/rosenGrp/Rpkgs')

library(tidyverse)
library(glmnet)
library(doParallel)

dir_work <- '/mnt/HA/groups/rosenGrp/embed_ag_samples/out'

seed <- 3245
args <- commandArgs(trailingOnly=TRUE)

nc <- as.integer(args[1])
kmer <- args[2]

path_out <- sprintf(file.path(dir_work,'seqs_ag_k%s_lasso_new.rds'),kmer)
cat(sprintf('Output will be saved to %s\n.',path_out))

train_test <- readRDS('/mnt/HA/groups/rosenGrp/embed_ag_samples/tables/train_test_ids_new.rds')

train_ids <- read_csv('/mnt/HA/groups/rosenGrp/embed_ag_samples/out/lasso_train_ids.csv.gz')
test_ids <- read_csv('/mnt/HA/groups/rosenGrp/embed_ag_samples/out/lasso_test_ids.csv.gz')

X_train <- read.csv(sprintf('/mnt/HA/groups/rosenGrp/embed_ag_samples/kmertbls/ag_%s_mer_train_X.csv',kmer),header=FALSE)
X_train <- data.frame(SampleID=train_ids$SampleID,X_train,stringsAsFactors = FALSE)
rownames(X_train) <- X_train$SampleID

X_test <- read.csv(sprintf('/mnt/HA/groups/rosenGrp/embed_ag_samples/kmertbls/ag_%s_mer_test_X.csv',kmer),header=FALSE)
X_test <- data.frame(SampleID=test_ids$SampleID,X_test,stringsAsFactors = FALSE)
rownames(X_test) <- X_test$SampleID

kmer_tbl <- rbind(X_train,X_test)

dat <- read_rds(file.path(dir_work,'ag_metadata.rds')) %>%
  select(PRIMARY_ID,body_site) %>%
  left_join(read_delim(file.path(dir_work,'ag_PRJEB11419.txt'),delim='\t') %>%
              select(PRIMARY_ID=secondary_sample_accession,SampleID=run_accession),
            by='PRIMARY_ID') %>%
  left_join(read_csv(file.path(dir_work,'ag_total_kmers.csv.gz')),by='SampleID') %>%
  filter(nreads >= 10000) %>%
  select(-PRIMARY_ID,-nreads) %>%
  filter(!is.na(body_site),
         body_site %in% c('UBERON:feces','UBERON:skin of hand','UBERON:skin of head','UBERON:tongue')) %>%
  mutate(body_site=as.character.factor(body_site),
         body_site=ifelse(body_site %in% c('UBERON:skin of hand','UBERON:skin of head'),'UBERON:skin',body_site)) %>%
  inner_join(kmer_tbl)

train <- dat %>% filter(SampleID %in% train_test$train)
test <- dat %>% filter(SampleID %in% train_test$test)

cat(sprintf('Creating cluster with %s cores.\n',nc))
registerDoParallel(nc)
cat('Performing lasso cross validation.\n')
set.seed(seed)

cv <- cv.glmnet(train %>% select(-body_site,-SampleID) %>% as.matrix(),train$body_site,
                family='multinomial',type.measure='class',parallel=TRUE,
                standardize=FALSE,nfolds=10)

y <- test$body_site
yhat <- as.vector(predict(cv,newx=test %>% select(-body_site,-SampleID) %>% as.matrix(),
                          type='class',s=cv$lambda.min))

out <- list(lasso=cv,train=train,test=test,results=data.frame(y=y,yhat=yhat))

cat('Saving results.\n')
saveRDS(out,path_out)


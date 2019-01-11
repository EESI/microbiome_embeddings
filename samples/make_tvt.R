#!/usr/bin/env Rscript

library(tidyverse)
library(caret)

dir_work <- '/mnt/HA/groups/rosenGrp/embed_ag_samples/out'

seed <- 3245

otu <- readRDS(file.path(dir_work,'ag_otu_table.rds'))

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
  filter(SampleID %in% otu$SampleID)

set.seed(seed)
idx_train <- createDataPartition(dat$body_site,times=1,p=0.9,list=TRUE)$Resample1

train <- dat$SampleID[idx_train]
test <- dat$SampleID[-idx_train]

saveRDS(list(train=train,test=test),file.path(dir_work,'train_test_ids.rds'))

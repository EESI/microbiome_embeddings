#!/usr/bin/env Rscript

library(tidyverse)
library(glmnet)
library(doMC)

dir_work <- '/mnt/HA/groups/rosenGrp/embed_ag_samples/out'

seed <- 3245
args <- commandArgs(trailingOnly=TRUE)

model <- args[1]
nc <- as.integer(args[2])

path_out <- gsub('remb','lasso',model)
path_out <- gsub('csv.gz','rds',path_out)
cat(sprintf('Output will be saved to %s\n.',path_out))

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
  inner_join(read_csv(model))

cat(sprintf('Creating cluster with %s cores.\n',nc))
registerDoMC(nc)
cat('Performing lasso cross validation.\n')
set.seed(seed)
cv <- cv.glmnet(dat %>% select(-body_site,-SampleID) %>% as.matrix(),dat$body_site,
                family='multinomial',type.measure='class',parallel=TRUE)

cat('Saving results.\n')
saveRDS(list(dat=dat,cv=cv),path_out)


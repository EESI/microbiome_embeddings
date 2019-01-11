#!/usr/bin/env Rscript

library(tidyverse)
library(Biostrings)
library(DECIPHER)

dir_work <- '/home/sw424/embed_samples'
dir_grp <- '/mnt/HA/groups/rosenGrp'

args <- commandArgs(trailingOnly=TRUE)
nc <- as.integer(args[1])
bs <- args[2]

path <- file.path(dir_grp,'embed_ag_samples',sprintf('alignment_forproteus_UBERON.%s.rds',bs))
seqs <- readRDS(path)

Nss <- 25000

set.seed(3453)
ssidx <- sample(seq_along(seqs),Nss)
seqs <- seqs[ssidx]

set.seed(432)
cat(sprintf('Aligning %s sequences from %s\n\n',length(seqs),bs))
aligned <- AlignSeqs(seqs,processors=nc,verbose=TRUE,
                     useStructures=FALSE,gapOpening=-25,gapExtension=-10)

saveRDS(aligned,
        file.path(dir_grp,'embed_ag_samples',sprintf('alignment_forstaph_%s.rds',bs)))

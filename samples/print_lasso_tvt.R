#!/usr/bin/env Rscript

.libPaths('/mnt/HA/groups/rosenGrp/Rpkgs')

library(caret)

get_class_acc <- function(x,classes){
    out <- sapply(classes,function(cl){
        z <- x[as.character(x[,1]) == cl,]
        mean(as.character(z[,1]) == as.character(z[,2]))
    })
    names(out) <- classes
    out
}

dir_work <- '/mnt/HA/groups/rosenGrp/embed_ag_samples/out'

fn_cluster <- file.path(dir_work,'seqs_ag_cluster_lasso_new.rds')
fns_otu <- list.files(dir_work,pattern='seqs_ag_[a-z]+_lasso_new.rds',full.names=TRUE)
fns_kmer <- list.files(dir_work,pattern='seqs_ag_k[0-9]+_lasso_new.rds',full.names=TRUE)
fns_emb <- list.files(dir_work,pattern='seqs_ag_sample_.*_lasso_tvt_new.*.rds',full.names=TRUE)

params <- gsub('^.*seqs_ag_sample_(.*).rds','\\1',fns_emb)

results <- vector(mode='numeric',length=length(fns_emb) + length(fns_otu) + length(fns_kmer))
results_class <- vector(mode='list',length=length(fns_emb) + length(fns_otu) + length(fns_kmer))

tax <- gsub('^.*seqs_ag_([a-z]+)_lasso_new.rds$','\\1',fns_otu)
kmer <- gsub('^.*seqs_ag_(k[0-9]+)_lasso_new.rds$','\\1',fns_kmer)

names(results) <- c(params,tax,kmer)
names(results_class) <- c(params,tax,kmer)


r <- readRDS(fns_otu[1])$results
classes <- unique(as.character(r[,1]))

for (i in seq_along(fns_emb)){
    fn <- fns_emb[i]
    r <- readRDS(fn)$results
    results[i] <- mean(as.character(r[,1]) == as.character(r[,2]))
    results_class[[i]] <- get_class_acc(r,classes)
    cat(sprintf('%s\n',names(results_class)[i]))
    print(confusionMatrix(r$yhat,r$y,mode='prec_recall'))
    cat('\n\n')
}

for (i in seq_along(fns_otu)){
    fn <- fns_otu[i]
    r <- readRDS(fn)$results
    results[i+length(fns_emb)] <- mean(as.character(r[,1]) == as.character(r[,2]))
    results_class[[i+length(fns_emb)]] <- get_class_acc(r,classes)
    cat(sprintf('%s\n',names(results_class)[i]))
    print(confusionMatrix(r$yhat,r$y,mode='prec_recall'))
    cat('\n\n')
}

for (i in seq_along(fns_kmer)){
    fn <- fns_kmer[i]
    r <- readRDS(fn)$results
    results[i+length(fns_emb)] <- mean(as.character(r[,1]) == as.character(r[,2]))
    results_class[[i+length(fns_emb)+length(fns_otu)]] <- get_class_acc(r,classes)
    cat(sprintf('%s\n',names(results_class)[i]))
    print(confusionMatrix(r$yhat,r$y,mode='prec_recall'))
    cat('\n\n')
}

cat('Overall Performance\n')
cat('Model\tAccuracy\n')
cat(sprintf('%s:\t%.4f\n',names(results),results))
cat('\nClass Performance\n')
cat('Model\tClass\tAccuracy\n')
for (i in seq_along(results_class)){
    for (j in seq_along(results_class[[i]])){
        cat(sprintf('%s:\t%s\t%.4f\n',
                    names(results_class)[i],
                    names(results_class[[i]])[j],
                    results_class[[i]][j]))
    }
    cat(sprintf('Mean performance: %s\n',mean(results_class[[i]])))
    cat('\n')
}

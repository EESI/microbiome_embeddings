#!/usr/bin/env Rscript

get_class_acc <- function(x,classes){
    out <- sapply(classes,function(cl){
        z <- x[as.character(x[,1]) == cl,]
        mean(as.character(z[,1]) == as.character(z[,2]))
    })
    names(out) <- classes
    out
}

dir_work <- '/mnt/HA/groups/rosenGrp/embed_ag_samples/out'

fn_otu <- file.path(dir_work,'seqs_ag_otu_lasso.rds')
fn_pca <- file.path(dir_work,'seqs_ag_pca_lasso.rds')
fn_cluster <- file.path(dir_work,'seqs_ag_cluster_lasso.rds')
fns_emb <- list.files(dir_work,pattern='seqs_ag_sample_.*_lasso_tvt.*.rds',full.names=TRUE)
params <- gsub('^.*seqs_ag_sample_(.*).rds','\\1',fns_emb)

results <- vector(mode='numeric',length=length(fns_emb) + 3)
results_class <- vector(mode='list',length=length(fns_emb) + 3)

names(results) <- c('otu_table','pca_table','cluster_table',params)
names(results_class) <- c('otu_table','pca_table','cluster_table',params)


r <- readRDS(fn_otu)$results
classes <- unique(as.character(r[,1]))
results[1] <- mean(as.character(r[,1]) == as.character(r[,2]))
results_class[[1]] <- get_class_acc(r,classes)

r <- readRDS(fn_pca)$results
results[2] <- mean(as.character(r[,1]) == as.character(r[,2]))
results_class[[2]] <- get_class_acc(r,classes)

r <- readRDS(fn_cluster)$results
results[3] <- mean(as.character(r[,1]) == as.character(r[,2]))
results_class[[3]] <- get_class_acc(r,classes)

for (i in seq_along(fns_emb)){
    fn <- fns_emb[i]
    r <- readRDS(fn)$results
    results[i+3] <- mean(as.character(r[,1]) == as.character(r[,2]))
    results_class[[i+3]] <- get_class_acc(r,classes)
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

#!/bin/bash

GRP=/mnt/HA/groups/rosenGrp/embed_cons/
OUT=$GRP/out/kegg/
params=($(ls $GRP/out/kegg/embeddings_split))
a="1e-05"

for p in ${params[@]}
do
    embeds=($(ls $GRP/out/kegg/embeddings_split/$p/*))
    out=$OUT/seqs_cluster_${p}_${a}_remb_raw.csv
    zcat ${embeds[0]} | head -1 >> $out
    for e in ${embeds[@]}
    do
        zcat $e | tail -n +2 >> $out
    done
    gzip $out
done

exit 0

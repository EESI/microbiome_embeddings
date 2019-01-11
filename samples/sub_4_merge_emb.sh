#!/bin/bash

GRP=/mnt/HA/groups/rosenGrp/embed_ag_samples
EMB=$GRP/out/embeddings_split

params=($(ls $EMB))

py=/mnt/HA/opt/python/3.6.1/bin/python3
script=$GRP/gen_emb.py
a="1e-05"

for p in ${params[@]}
do
    echo "Merging ${p}."
    embeddings=($(ls $EMB/$p/*))
    out="$GRP/out/seqs_ag_sample_${p}_${a}_remb_raw.csv"
    rm $out
    i=0
    n=1
    echo -ne "[$n]\t"
    for e in ${embeddings[@]}
    do
        zcat $e >> $out
        echo -ne "."
        ((i++))
        [ $i -eq 100 ] && { ((n+=100)); echo -ne "\n[$n]\t" ; i=0; }
    done
    echo "Compressing."
    gzip -f $out
done

exit 0

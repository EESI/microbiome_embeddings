#!/bin/bash

PATH=/mnt/HA/groups/rosenGrp/embed_samples/data/ag/embeddings_split/10_256_5_50_10_1e-06_100
out=/mnt/HA/groups/rosenGrp/embed_ag_samples/ag_sample_file_split.txt

i=1
for p in $PATH/*_raw.csv.gz
do
    if ((i == 1)); then
        echo -ne $p > $out
    elif ! ((i % 20)); then
        echo -ne "\n${p}" >> $out
    else
        echo -ne ",${p}" >> $out
    fi
    ((i++))
done

exit 0

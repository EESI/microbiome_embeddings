#!/bin/bash

GRP=/mnt/HA/groups/rosenGrp/embed_cons
CLUSTERS=$GRP/out/kegg/clusters

clusters=($(ls $CLUSTERS/cl_*))
out=$GRP/out/kegg/clusters.fasta

for cl in ${clusters[@]}
do
    echo ">${cl##*/}" >> $out
    awk '!/^>/ { printf "%s", $0; n = "\n" }
    END { printf "%s", n }
    ' $cl >> $out
done

exit 0

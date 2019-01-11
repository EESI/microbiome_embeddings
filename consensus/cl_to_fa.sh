#!/bin/bash

GRP=/mnt/HA/groups/rosenGrp/embed_cons
CLUSTERS=$GRP/out/kegg/clusters

clusters=($(ls $CLUSTERS/cl_*))
OUT=$GRP/out/kegg/clusters_fasta

mkdir -p $OUT

for cl in ${clusters[@]}
do
    fn=${cl##*/}
    awk '!/^>/ { printf "%s", $0; n = "\n" }
    /^>/ { print n $0; n = "" }
    END { printf "%s", n }
    ' $cl >> $OUT/${fn}.fasta
done

exit 0

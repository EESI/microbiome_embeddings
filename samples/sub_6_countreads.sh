#!/bin/bash

GRP=/mnt/HA/groups/rosenGrp
FAS=$GRP/embed_samples/data/ag/fasta

out=$GRP/embed_ag_samples/out/ag_total_kmers.csv
fas=($(ls $FAS/*))

i=0
n=1
echo "Counting reads per sample"
echo -ne "[$n]\t"
echo "SampleID,nreads" > $out 
for fa in ${fas[@]}
do
    basename=${fa##*/}
    sampleid=${basename%.*}
    nlines=$(cat $fa | wc -l)
    nreads=$(($nlines / 2))

    echo "$sampleid,$nreads" >> $out
    echo -ne "."
    ((i++))
    [ $i -eq 100 ] && { ((n+=100)); echo -ne "\n[$n]\t" ; i=0; }
done

echo -ne "\nCompressing"
gzip -f $out

exit 0

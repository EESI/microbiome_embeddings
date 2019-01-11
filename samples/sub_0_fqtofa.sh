#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=12:00:00
#$ -P rosenPrj
#$ -pe shm 1
#$ -l mem_free=30G
#$ -l h_vmem=32G
#$ -q all.q

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load python/3.6-current

GRP=/mnt/HA/groups/rosenGrp
FQS=$GRP/embed_reads/ag/data
FAS=$GRP/embed_samples/data/ag/fasta

mkdir $FAS

fqs=($(ls $FQS/*))

for fq in ${fqs[@]};
do
    fn=$(basename "$fq")
    id=${fn%.fastq.gz}
    echo "Converting ${fn}."
    zcat $fq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > ${FAS}/${id}.fasta
done

exit 0

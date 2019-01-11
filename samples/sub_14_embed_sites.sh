#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=01:00:00
#$ -P rosenPrj
#$ -pe shm 1
#$ -l mem_free=7G
#$ -l h_vmem=8G
#$ -q all.q
#$ -t 1-749

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load python/3.6-current

GRP=/mnt/HA/groups/rosenGrp/embed_ag_samples
OUT=$GRP/out/site_embedding_split

mkdir $OUT
rm $OUT/seqs_10_256_5_50_10_1e-06_100_embedding_np_*.csv

py=/mnt/HA/opt/python/3.6.1/bin/python3
script=create_site_embeddings_split.py

$py $script $SGE_TASK_ID

exit 0

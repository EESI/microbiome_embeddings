#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=48:00:00
#$ -P rosenPrj
#$ -pe shm 1
#$ -l mem_free=60G
#$ -l h_vmem=64G
#$ -q all.q
#$ -t 1-2

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load python/3.6-current

GRP=/mnt/HA/groups/rosenGrp/embed_cons
MODELS=/home/sw424/embed_samples/data/models

py=/mnt/HA/opt/python/3.6.1/bin/python3
script=$GRP/embed_reads_clust.py

fas=$GRP/out/kegg/clusters.fasta
models=(gg_6_256_5_50_10_1e-06_100_model.pkl gg_10_256_5_50_10_1e-06_100_model.pkl)

$py $script $fas $MODELS/${models[$SGE_TASK_ID-1]}

exit 0

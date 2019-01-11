#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=02:00:00
#$ -P rosenPrj
#$ -pe shm 1
#$ -l mem_free=14G
#$ -l h_vmem=16G
#$ -q all.q

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load python/3.6-current

GRP=/mnt/HA/groups/rosenGrp/embed_cons
MODELS=/home/sw424/embed_samples/data/models

py=/mnt/HA/opt/python/3.6.1/bin/python3
script=$GRP/gen_clust_emb.py

$py $script

exit 0

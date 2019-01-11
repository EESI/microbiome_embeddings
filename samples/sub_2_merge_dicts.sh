#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=05:00:00
#$ -P rosenPrj
#$ -pe shm 1
#$ -l mem_free=30G
#$ -l h_vmem=32G
#$ -q all.q
#$ -t 1-4

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load python/3.6-current

GRP=/mnt/HA/groups/rosenGrp
PARAMS=$GRP/embed_samples/data/ag/total_kmers_split

py=/mnt/HA/opt/python/3.6.1/bin/python3
script=$GRP/embed_samples/merge_dicts.py

params=($(ls $PARAMS))

$py $script ${params[$SGE_TASK_ID-1]}

exit 0

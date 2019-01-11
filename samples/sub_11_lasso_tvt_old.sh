#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=02:00:00
#$ -P rosenPrj
#$ -l ua=haswell
#$ -pe shm 10
#$ -l mem_free=5G
#$ -l h_vmem=6G
#$ -t 1-4

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load samtools/1.3.1

WRK=/mnt/HA/groups/rosenGrp/embed_ag_samples

script=$WRK/lasso_tvt.R
models=($(ls $WRK/out/seqs_ag_sample*remb*csv.gz))

$script ${models[$SGE_TASK_ID-1]} $NSLOTS

exit 0

#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=24:00:00
#$ -P rosenPrj
#$ -pe shm 4
#$ -l mem_free=180G
#$ -l h_vmem=200G
#$ -t 4-10:2

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load samtools/1.3.1

WRK=/mnt/HA/groups/rosenGrp/embed_ag_samples

script=$WRK/lasso_kmer.R
tables=($(ls $WRK/tables/ag_*_table_new.rds))

$script $NSLOTS $SGE_TASK_ID

exit 0

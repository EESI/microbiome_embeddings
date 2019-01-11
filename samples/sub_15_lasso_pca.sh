#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=24:00:00
#$ -P rosenPrj
#$ -pe shm 10
#$ -l mem_free=15G
#$ -l h_vmem=16G

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load samtools/1.3.1

WRK=/mnt/HA/groups/rosenGrp/embed_ag_samples

script=$WRK/lasso_otu_pca.R

$script $NSLOTS

exit 0

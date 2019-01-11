#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=48:00:00
#$ -P rosenPrj
#$ -pe shm 48
#$ -l mem_free=239G
#$ -l h_vmem=240G

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load samtools/1.3.1

WRK=/mnt/HA/groups/rosenGrp/embed_ag_samples

script=$WRK/align_site.R

BS="tongue"

$script $NSLOTS $BS

exit 0

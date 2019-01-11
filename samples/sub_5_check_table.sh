#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M zz374@drexel.edu
#$ -l h_rt=01:00:00
#$ -P rosenPrj
#$ -pe shm 1
#$ -l mem_free=15G
#$ -l h_vmem=16G
#$ -q all.q

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load python/3.6-current

GRP=/mnt/HA/groups/rosenGrp/embed_ag_samples

py=/mnt/HA/opt/python/3.6.1/bin/python3
script=$GRP/check_table.py

$py $script

exit 0

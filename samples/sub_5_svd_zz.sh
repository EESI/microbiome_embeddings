#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M zz374@drexel.edu
#$ -l h_rt=50:00:00
#$ -P rosenPrj
#$ -pe shm 16
#$ -l ua=sandybridge
#$ -l mem_free=3G
#$ -l h_vmem=2G
#$ -q long.q

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load python/intelpython/3

GRP=/mnt/HA/groups/rosenGrp/embed_ag_samples

script=$GRP/svd.py

python $script

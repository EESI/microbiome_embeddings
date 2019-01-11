#!/bin/bash
#$ -S bin/bash
#$ -j y
#$ -cwd
#$ -M sw424@drexel.edu
#$ -l h_rt=01:00:00
#$ -P rosenPrj
#$ -pe shm 1
#$ -l mem_free=30G
#$ -l h_vmem=32G
#$ -q all.q

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
module load vsearch/gcc/2.7.1

USERNAME=sw424

SCRATCH=/scratch/$USERNAME/vsearch
WORK=/home/$USERNAME/embed_samples
OUT=/mnt/HA/groups/rosenGrp/embed_cons/out/kegg

name=kegg
reads=seqs.fasta
results=seqs_cluster_results.txt
vcons=seqs_cluster_vcons.txt

mkdir -p $SCRATCH
mkdir -p $OUT/clusters

cp $WORK/data/$name/$reads $SCRATCH/

cd $SCRATCH

vsearch --cluster_smallmem $reads --consout $SCRATCH/$vcons \
    --clusters $SCRATCH/cl_ --clusterout_id --id 0.8 \
    --usersort --uc $SCRATCH/$results

mv $results $vcons $OUT/
mv cl_* $OUT/clusters/

exit 0

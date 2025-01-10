#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=SEACRMakeBedgraph
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log

echo "[$0] $SLURM_JOB_NAME $@" # log the command line

module purge
source activate base
conda activate rnaPseudo

filename=AgoSorted.txt # $1

linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
    #filename=$(basename -- "$FILE")
      EXT=${line::-9}
      echo $EXT
      samtools view -f 0x2 -b ${EXT}.Read.bam > ${EXT}.paired.bam
      bedtools bamtobed -bedpe -i ${EXT}.paired.bam > ${EXT}.bed
      awk '$1==$4 && $6-$2 < 1000 {print $0}' ${EXT}.bed > ${EXT}.clean.bed
      cut -f 1,2,6 ${EXT}.clean.bed | sort -k1,1 -k2,2n -k3,3n -T TempSort > ${EXT}.fragments.bed
      bedtools genomecov -bg -i ${EXT}.fragments.bed -g genome.fa.fai > ${EXT}.fragments.bedgraph
      rm ${EXT}.paired.bam
      rm ${EXT}.bed
      rm ${EXT}.clean.bed
      rm ${EXT}.fragments.bed
    fi
    linenum=$((linenum + 1))
done < $filename

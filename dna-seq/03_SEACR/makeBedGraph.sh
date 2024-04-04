#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=sortBam
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

module purge
module load anaconda
conda activate rnaPseudo


for FILE in *.sorted.bam
do
  #filename=$(basename -- "$FILE")

  EXT=${FILE::-11}
  echo $EXT
  samtools view -f 0x2 -b ${EXT}.sorted.bam > ${EXT}.paired.bam
  bedtools bamtobed -bedpe -i ${EXT}.paired.bam > ${EXT}.bed
  awk '$1==$4 && $6-$2 < 1000 {print $0}' ${EXT}.bed > ${EXT}.clean.bed
  cut -f 1,2,6 ${EXT}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${EXT}.fragments.bed
  bedtools genomecov -bg -i ${EXT}.fragments.bed -g genome.fa.fai > ${EXT}.fragments.bedgraph
done

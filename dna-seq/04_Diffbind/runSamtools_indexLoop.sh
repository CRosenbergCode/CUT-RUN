#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --job-name=samtools_index
#SBATCH --partition=day-long-cpu
#SBATCH --output=%x.%j.log # gives jobname.ID.log
# Available partitions 
# day-long-cpu 
# day-long-highmem
# exp-gpu
# short-cpu
# short-gpu 
# short-highmem
# week-long-cpu
# week-long-gpufile
# week-long-highmem

#migrate to Hisat folder with bam files

# to run this
#sbatch runSamtools_indexLoop.sh
set -e #creates log at first error

module purge
source activate base
conda activate samtools

mkdir -p $HOME/tmp
export TMP=$HOME/tmp
export TMPDIR=$TMP

linenum=0
bams=$(ls *.bam)

for file in $bams
do
  TWO=${file/.bam/.bai}
  echo ${TWO}
  samtools sort -o sorted_${file} $file
  samtools index sorted_${file} > $TWO
done

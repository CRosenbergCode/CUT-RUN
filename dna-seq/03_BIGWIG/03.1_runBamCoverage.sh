#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=BamCoverage
#SBATCH --output=%x.%j.out
#SBATCH --time=3:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

module purge
module load anaconda
conda activate deeptools_kernel

for FILE in *.bam
do
  pref=${FILE::-4}
  #pref=${pref:21}
  echo ${pref}
  bamCoverage --bam ${FILE} -o bin10/${pref}.bw \
      --binSize 10 \
      --normalizeUsing RPGC \
      --effectiveGenomeSize 1220000000 \ #Current size is for Aedes Aegypti
      --extendReads 
done

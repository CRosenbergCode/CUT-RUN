#!/bin/bash
#SBATCH --partition=amilan
#SBATCH --job-name=BamCoverage
#SBATCH --output=%x.%j.out
#SBATCH --time=3:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

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
      --effectiveGenomeSize 12200000000 \
      --extendReads
done
#Adjust --binSize to adjust how granular the result bigwig file is, higher is lower resolution

#RPGC (reads per genome coverage) is the normalization metric chosen. Normalization methods are generally equivalent other than scaling factor.

#effectiveGenomeSize of 12200000000 is based upon the Aedes aegypti genome and should be changed for other organisms

#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=hisatDNALoop
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
module load anaconda
conda activate rnaPseudo


for FILE in ../01_fastp/trimmed/trimmed.*R1_001.fastq.gz
do
  TWO=${FILE/R1/R2} #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
  OUT=${TWO:8:(${#TWO}-23)} #Remove trimmed, and end. 23 is 8 (trimmed.) + 15 (R1_001.fastq.gz)
  #echo ${TWO}
  hisat2 --phred33 -p 8 -x aedesGenome -1 ${FILE} -2 ${TWO} -S HisatOutput/${OUT}.sam
  #Increase -p to increase the number of threads
  #--phred33 is for illumina encoding, but other encoding types for fastq scores may require differnt parameters
done

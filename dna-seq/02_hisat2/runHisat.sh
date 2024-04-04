#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=hisatDNAPairedOnly
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

module purge
module load anaconda
conda activate rnaPseudo


for FILE in trimmed.*R1_001.fastq.gz
do
  TWO=$(echo $FILE | rev | sed s/./2/14 | rev) #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
  OUT=${TWO:8:(${#TWO}-23)} #Remove trimmed, and end. 23 is 8 (trimmed.) + 15 (R1_001.fastq.gz)
  #echo ${TWO}
  hisat2 --phred33 -p 16 -x genome -1 ${FILE} -2 ${TWO} -S hisatPairStrict/${OUT}.sam --no-mixed --no-discordant
done


#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=sortBam
#SBATCH --output=%x.%j.out
#SBATCH --time=0:40:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

module purge
module load anaconda
conda activate rnaPseudo


for FILE in ../02_hisat2/*.bam
do
  EXT=${FILE:12:-4}
  echo $EXT
  samtools sort -N -@ 8 $FILE -o $EXT.Read.bam
done

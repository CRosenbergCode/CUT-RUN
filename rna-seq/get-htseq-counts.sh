#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=htseqCount
#SBATCH --output=%x.%j.out
#SBATCH --time=6:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pegah.eizadkhah@colostate.edu

module purge
module load anaconda
conda activate cutruntools2.1

for FILE in hisatAligned/*.sam
do
  fileName=$(basename $FILE .sam)
  htseq-count -n 5 --stranded=reverse $FILE VectorBase-66_AaegyptiLVP_AGWG.gtf > htseqCounts/${fileName}_HSCounts.txt
done

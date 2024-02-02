#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=srainstallAedesFaire
#SBATCH --output=%x.%j.out
#SBATCH --time=05:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colorado.edu

module purge
module load anaconda
conda activate cutruntools2.1
module load sra-toolkit/3.0.0

fasterq-dump  -p  -O /scratch/alpine/$USER/sra SRR2530419
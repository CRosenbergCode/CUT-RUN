#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=hisatAedes68Build
#SBATCH --output=%x.%j.out
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

module purge
source activate base
conda activate rnaPseudo

hisat2-build -p 8 VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta hisatIndex/Aedesae68Index

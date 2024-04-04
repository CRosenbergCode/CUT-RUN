#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=hisatEColiBuild
#SBATCH --output=%x.%j.out
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

module purge
module load anaconda
conda activate rnaPseudo


hisat2-build -p 8 EColiIndex/mergedGenome.fa genomeEColi

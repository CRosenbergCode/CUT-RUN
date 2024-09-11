#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=hisatYeastBuild
#SBATCH --output=%x.%j.out
#SBATCH --time=1:30:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
module load anaconda
conda activate rnaPseudo


hisat2-build -p 8 genomeCombFull.fa aedesGenome
#.fa file is the fasta file for the genome
#Second parameter (aedesGenome) is the name for the indexes

#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=star-indices
#SBATCH --output=%x.%j.out
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pegah.eizadkhah@colostate.edu

module load anaconda
conda activate cutruntools2.1

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /scratch/alpine/c836951585@colostate.edu/DAS/STAR_indices --genomeFastaFiles /projects/c836951585@colostate.edu/CutnRunRNAPipelineTesting/VectorBase-66_AaegyptiLVP_AGWG_Genome.fasta --sjdbGTFfile /projects/c836951585@colostate.edu/CutnRunRNAPipelineTesting/VectorBase-66_AaegyptiLVP_AGWG_Orf50.gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 149
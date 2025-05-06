#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --job-name=STAR_indices
#SBATCH --partition=day-long-cpu
#SBATCH --output=%x.%j.log

# Instructions for how to use:
if [[ -z "$@" ]]; then
	echo -e "\nTo use this script, run the following command:"
	echo -e "sbatch generate_STAR_indices.sh path/to/your_genome.fasta path/to/your.gff\n"
	echo -e "Example:\nsbatch generate_STAR_indices.sh VectorBase-66_AaegyptiLVP_AGWG_Genome.fasta VectorBase-66_AaegyptiLVP_AGWG_Orf50.gff\n"
	exit 1
fi

source $HOME/miniconda3/bin/activate base
conda init
conda activate rMATS

mkdir STAR_indices

STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir ./STAR_indices \
	--genomeFastaFiles $1 \
	--sjdbGTFfile $2 \
	--sjdbGTFtagExonParentTranscript Parent \
	--sjdbOverhang 149


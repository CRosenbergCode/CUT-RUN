#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --job-name=STAR
#SBATCH --partition=day-long-cpu
#SBATCH --output=%x.%j.log

source $HOME/miniconda3/bin/activate base
conda init
conda activate rMATS

curdir=$(pwd)
fastqdir=$1
cd $fastqdir
for fastq in $(ls *.fastq.gz | grep R1_[0-9]*.fastq.gz); do
	fast2=${fastq//R1/R2}
	fastpref=${fastq/_R1_[0-9]*.fastq.gz/}
	fastpref=${fastpref/R1_[0-9]*.fastq.gz/}
	STAR --runThreadN 8 \
		--genomeDir ${curdir}/STAR_indices \
		--readFilesIn $fastq $fast2 \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--readFilesCommand zcat \
		--outFileNamePrefix $fastpref
done


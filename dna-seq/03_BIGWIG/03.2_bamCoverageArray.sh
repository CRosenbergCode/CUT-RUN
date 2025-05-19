#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --job-name=bamcoverage_array
#SBATCH --partition=short-cpu
#SBATCH --output=%x.%j.log # gives jobname.ID.log
# Available partitions
# day-long-cpu
# day-long-gpu
# day-long-highmem
# exp-gpu
# short-cpu*
# short-gpu
# short-highmem
# week-long-cpu
# week-long-gpu
# week-long-highmem

# ================ 03.2_bamCoverageArray.sh ==================
#
#  DESCRIPTION: A Slurm array script to generate bigwigs from
#            sorted bam files listed in the input file.
#
#  USAGE: sbatch --array=0-(number of bams minus 1) 03.2_bamCoverageArray.sh $1
#
#  INPUTS: $1 - A plaintext file containing the paths to the
#               (sorted) bams. Each line should be the path to
#               a single bam.
#
#  OUTPUT: Bigwig (.bw) files in the same directory as the bams
#
# ============================================================

source activate base
conda init
conda activate deeptools_kernel

export TMPDIR="$2"

# Run bamCoverage for the assigned sorted bam
# 4 processors may or may not be optimal, but the original scripts used nthreads=4.
linenum=0
while read -r line
do
	if [ $SLURM_ARRAY_TASK_ID -eq $linenum ] ; then
		filename=$(echo $line | rev | cut -d/ -f1 | rev)
		filedir=$(dirname -- "$line")
		cd filedir
		FILE="sorted_${filename}"
		pref=${filename::-4}
		echo ${pref}
		bamCoverage --bam ${FILE} -o ${pref}.bw \
			--binSize 10 \
			--normalizeUsing RPGC \
			--effectiveGenomeSize 12200000000 \
			--extendReads \
			--numberOfProcessors 4
	fi
	linenum=$((linenum + 1))
done < $1


#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --job-name=sort_index_bams
#SBATCH --partition=day-long-cpu
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


# This is a helper script and not intended to be run on its own.
# It sorts/indexes all the bams in the current directory.

# There's no reason you CAN'T run it on its own, though.
# To do so, sbatch this script from whatever directory
# you have your unsorted bams in, like so:

# sbatch path/to/sort_index_helper.sh


source $HOME/miniconda3/bin/activate base
conda init
conda activate samtools

echo "pwhuhughuhfhfghGUHRUHUH: $(pwd)"

# Sort and index bams
# Currently I have this configured to use 16 threads.
# If there is a more optimal number, replace the #SBATCH -n 16
# with it.
for bam in *.bam; do
	samtools sort -@$SLURM_NTASKS $bam -o sorted_$bam
	samtools index -@$SLURM_NTASKS sorted_$bam
done
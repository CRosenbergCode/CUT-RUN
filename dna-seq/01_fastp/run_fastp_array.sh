# I am writing this in a bit of a rush, so here are the basics:
#
# 1. This script should be in the 01_fastp directory
#
# 2. Your fastq files should be in the "raw" folder.
#
# 3. Run this script like so:
#
#	bash run_fastp_array.sh
#
# 4. That's it! Check that it's running properly with sacct or squeue.

ls raw/*_R1_*.fastq.gz > fastparray.tmp

num_arr=$(cat fastparray.tmp | wc -l)
sbatch --array=0-$((num_arr - 1)) fastp_array.slurm fastparray.tmp


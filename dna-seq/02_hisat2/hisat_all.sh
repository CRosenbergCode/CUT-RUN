# ===================== hisat_all.sh =======================
#
#  DESCRIPTION: A short non-slurm script to successively
#            run the three hisat scripts with no needed
#            input from the user.
#
#  USAGE (DON'T RUN WITH SLURM): bash hisat_all.sh
#
#  OUTPUT: .sam and .bam files in the HisatAligned folder.
#
# =========================================================

# If there are indices already, do not re-generate them.
# If there aren't, then run 02.1_runBuildAedesae68Genome. 
if ls hisatIndex | grep -q "\.ht2" ; then
	# First, detect how many array jobs to dispatch using the same logic as the 02.2 script.
	num_arr=$(ls ../01_fastp/trimmed/ -p | grep -v / | grep ".fastq.gz" | grep "R1" | grep "trimmed." | wc -l)
	# Then, run Hisat Align Array, get the job ID for making 02.3 a dependent job.
	array_message=$(sbatch --array=0-$((num_arr - 1)) 02.2_runAlignArrayHisat.sh)
	array_id=$(echo $array_message | grep -oh "[0-9]*$")
else
	# Run build genome, get the job ID for making 02.2 a dependent job.
	build_message=$(sbatch 02.1_runBuildAedesae68Genome.sh)
	build_id=$(echo $build_message | grep -oh "[0-9]*$")

	# First, detect how many array jobs to dispatch using the same logic as the 02.2 script.
	num_arr=$(ls ../01_fastp/trimmed/ -p | grep -v / | grep ".fastq.gz" | grep "R1" | grep "trimmed." | wc -l)
	# Then, run Hisat Align Array, get the job ID for making 02.3 a dependent job.
	array_message=$(sbatch --depend=afterok:${build_id} --array=0-$((num_arr - 1)) 02.2_runAlignArrayHisat.sh)
	array_id=$(echo $array_message | grep -oh "[0-9]*$")
fi

# Run 02.3_runBamLoop.sh
cd HisatAligned
sbatch --depend=afterok:${array_id} ../02.3_runBamLoop.sh
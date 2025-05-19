# =================== 01_fastqToBamBai.sh =====================
#
#  DESCRIPTION: A short NON-SLURM script to successively
#            run fastp and the three hisat scripts with no
#            input from the user, and sort/index resulting bams.
#
#  USAGE (DON'T RUN WITH SLURM): bash 01_fastqToBamBai.sh
#
#  BEFORE RUNNING: This script should be in the 01_fastp folder.
#            Your raw .fastq.gz files should be in the
#            01_fastp/raw folder. Ensure that the 01_fastp/trimmed
#            and the 02_hisat2/HisatAligned folders are empty
#            to avoid mixing your runs. You can always move old
#            files to another folder if you want to save them.
#
#  OUTPUT: sorted .bam and .bai files in the 02_hisat2/HisatAligned
#          directory, and other products of the fastq and hisat2
#          scripts. (trimmed fastq.gz files, .ht2, .bam files)
#
# =============================================================

module load slurm

# Run fastp, get the slurm job ID setting dependencies
fastp_message=$(sbatch runFastpHTML.sh)
fastp_id=$(echo $fastp_message | grep -oh "[0-9]*$")

# When fastp finishes, run build E. coli & get job ID
cd ../02_hisat2
# Detect how many array jobs to dispatch using the same logic as the 02.2 script.
num_arr=$(ls ../01_fastp/raw/ -p | grep -v / | grep ".fastq.gz" | grep "R1" | wc -l)

# Check if we need to run the 02.1 build genome script.
if ls hisatIndex | grep -q "\.ht2" ; then
	# Run Hisat Align Array & get job ID
	array_message=$(sbatch --depend=afterok:${fastp_id} --array=0-$((num_arr - 1)) 02.2_runAlignArrayHisat.sh)
	array_id=$(echo $array_message | grep -oh "[0-9]*$")
else
	# Run 02.1 Build Genome and get the job ID.
	build_message=$(sbatch --depend=afterok:${fastp_id} 02.1_runBuildAedesae68Genome.sh)
	build_id=$(echo $build_message | grep -oh "[0-9]*$")
	
	# When build genome finishes, run Hisat Align Array & get job ID
	array_message=$(sbatch --depend=afterok:${build_id} --array=0-$((num_arr - 1)) 02.2_runAlignArrayHisat.sh)
	array_id=$(echo $array_message | grep -oh "[0-9]*$")
fi

# When all array jobs finish, run the bam loop
cd HisatAligned
bamloop_message=$(sbatch --depend=afterok:${array_id} ../02.3_runBamLoop.sh)
bamloop_id=$(echo $bamloop_message | grep -oh "[0-9]*$")

sbatch --depend=afterok:${bamloop_id} ../../01_fastp/sort_index_helper.sh

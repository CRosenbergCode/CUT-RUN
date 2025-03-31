# =================== 01_fastqToBam.sh =====================
#
#  DESCRIPTION: A short non-slurm script to successively
#            run fastp and the three hisat scripts with no
#            input from the user.
#
#  USAGE (DON'T RUN WITH SLURM): ./01_fastqToBam.sh
#
#  OUTPUT: .bam files in the 02_hisat2/HisatAligned directory,
#          and other products of the fastq and hisat2 scripts.
#          (trimmed fastq.gz files, .ht2, .sam files)
#
# ==========================================================

module load slurm

# Run fastp, get the slurm job ID setting dependencies
fastp_message=$(sbatch runFastpHTML.sh)
fastp_id=$(echo $fastp_message | grep -oh "[0-9]*$")

# When fastp finishes, run build E. coli & get job ID
cd ../02_hisat2
build_message=$(sbatch --depend=afterok:${fastp_id} 02.1_runBuildEColi.sh)
build_id=$(echo $build_message | grep -oh "[0-9]*$")

# Detect how many array jobs to dispatch using the same logic as the 02.2 script.
num_arr=$(ls ../01_fastp/raw/ -p | grep -v / | grep ".fastq.gz" | grep "R1" | wc -l)

# When build E. coli finishes, run Hisat Align Array & get job ID
array_message=$(sbatch --depend=afterok:${build_id} --array=0-$((num_arr - 1)) 02.2_runAlignArrayHisat.sh)
array_id=$(echo $array_message | grep -oh "[0-9]*$")

# When all array jobs finish, run the bam loop
cd HisatAligned
sbatch --depend=afterok:${array_id} ../02.3_runBamLoop.sh

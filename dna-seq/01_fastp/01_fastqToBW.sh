# ==================== 01_fastqToBW.sh ======================
#
#  DESCRIPTION: A short non-slurm script to successively
#            run fastp, three hisat scripts, and bigwig with
#            no input from the user.
#
#  USAGE (DON'T RUN WITH SLURM): bash 01_fastqToBW.sh $1
#
#     ARGUMENTS:
#                $1 - Replace this with the name of your comparisons file.
#                     This is an optional argument. If you have
#                     a list of runs you would like to compare,
#                     include them in this file using the following
#                     syntax:
#                     
#                     1. Each line has a set of comma-separated file
#                        names, with the last file being the "input",
#                        or what you're comparing against:
#                        
#              Example:  Treatment.fastq.gz,Input.fastq.gz
#                     
#                     2. You can compare multiple treatments to the
#                        same input by adding them to the beginning
#                        of the line. The input should be at the end.
#                        
#              Example:  Treatment_Ac.fastq.gz,Treatment_Me.fastq.gz,Treatment_Ago2.fastq.gz,Input.fastq.gz
#
#                     3. Example line:
#
#              BF-d3-merged_me.fastq.gz,BF-d3-merged_ac.fastq.gz,BF-d3-merged_input.fastq.gz
#
#  BEFORE RUNNING: This script should be in the 01_fastp folder.
#            Your raw .fastq.gz files should be in the
#            01_fastp/raw folder. Ensure that the 01_fastp/trimmed
#            and the 02_hisat2/HisatAligned folders are empty
#            to avoid mixing your runs. You can always move old
#            files to another folder if you want to save them.
#
#  OUTPUT: .bw files in the 02_hisat2/HisatAligned directory,
#          (THIS IS NOT A TYPO, NOT IN THE BIGWIG DIRECTORY)
#          and other products of the fastq and hisat2 scripts.
#          (trimmed fastq.gz files, .ht2, .bam, bai files)
#
# ===========================================================

# If a comparisons file is given, make sure it's formatted properly.
if [[ -n "$1" ]]; then
  # Ensure file is not DOS-formatted, because MobaRTE is actually awful
  if cat -v $1 | grep "\^M" ; then
    echo -e "DOS Line endings found in your comparisons file: ${1}\nIf you're using MobaXTerm, open your comparisons file with the\nMobaXTerm text editor and click the little penguin button at the top, then save!"
    echo "Exiting. Please fix your comparisons file."
    exit
  fi

  # Get global path for later.
  COMP=$(readlink -f $1)

  # Make sure it ends in a newline; if not then add one. (POSIX formatting)
  if [ -n "$(tail -c1 $1)" ]; then
    echo "" >> $1;
  fi
  
  sed -e "s:_R[1 2]_[0-9][0-9][0-9].fastq.gz:.bam:g" $1 > convBam_${1}
fi

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

cd ../../03_BIGWIG
if [[ -n $1 ]]; then
	sbatch --depend=afterok:${bamloop_id} SICC_array.sh ../02_hisat2/HisatAligned 1 $COMP
else
	sbatch --depend=afterok:${bamloop_id} SICC_array.sh ../02_hisat2/HisatAligned 1
fi

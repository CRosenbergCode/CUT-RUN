#!/bin/bash

#SBATCH --partition=short-cpu
#SBATCH --job-name=SICC_scheduled
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Although this has a slurm header, you don't need to run this
# as a slurm job. It only has a slurm header so the 01_fastqToBW
# script can schedule it using dependencies.

# ==================== SICC_array.sh ======================
#
#  DESCRIPTION: A short script to successively run array-style
#            bigwig scripts.
#
#  USAGE (don't use slurm): bash SICC_array.sh $1 $2 $3
#
#  INPUT: $1 - Path to the folder with bams to be
#              converted into bigwigs.
#         $2 - Whether to sort and index bams. Set this to 0
#              if your bams are already sorted/indexed, or 1
#              if they still need to be indexed.
#         $3 - Path to the file containing bigwig
#              comparisons to be made. This is optional, and
#              only if you want to do a comparison step.
#              You may simply not pass a second argument.
#
#              Comparison data is written using the following rules:
#                    
#                      1. Each line is comma-separated, with the input at the end:
#
#                          Treatment.bam,Input.bam
#
#                      2. You can compare multiple treatments to the same control
#                          like so. They will each produce a different bigwig:
#
#                          Treatment_Ac.bam,Treatment_Me.bam,Treatment_Ago2.bam,Input.bam
#
#                      3. Example line:
#
#                          BF-d3-merged_me.bam,BF-d3-merged_ac.bam,BF-d3-merged_input.bam
#
#  OUTPUT: Sorted .bam files, .bai files, .bw files and
#          (optionally) comparison .bw files.
#
# =========================================================

curDir=$(pwd)

# If a comparisons file is given, make sure it's formatted properly.
if [[ -n "$3" ]]; then
  # Ensure file is not DOS-formatted, because MobaRTE is actually awful
  if cat -v $3 | grep "\^M" ; then
    echo -e "DOS Line endings found in your comparisons file: ${2}\nIf you're using MobaXTerm, click the little penguin button\nwhen editing your comparisons file in the text editor and save your file!"
    echo "Exiting. Please fix your comparisons file."
    exit
  fi

  # Get global path for later.
  COMP=$(readlink -f $3)

  # Make sure it ends in a newline; if not then add one. (POSIX formatting)
  if [ -n "$(tail -c1 $3)" ]; then
    echo "" >> $3;
  fi
fi

# Change into the relevant directory
cd $1

# Setting a custom temp directory. DeepTools will otherwise default to Riviera's /tmp
# directory, which is teeny (~2GB) and I think also shared between everyone.
mkdir dt_tmp
pwdstring=$(pwd)
TMPDIR="${pwdstring}/dt_tmp"

# Detect how big our sort/index job array should be
num_si_arr=$(ls *.bam | wc -l)
# Make the file with bams to operate on, for 03.1
ls *.bam > SICC.tmp
sed -e "s:^:sorted_:g" SICC.tmp > SICC_sorted.tmp
if [[ $2 -eq 0 ]] ; then
	# Detect how big our coverage job array should be
	num_coverage_arr=$num_si_arr
	# Schedule the coverage array job to run after the sort/indexing is done.
	coverage_message=$(sbatch --array=0-$((num_coverage_arr - 1)) ${curDir}/03.2_bamCoverageArray.sh SICC_sorted.tmp $TMPDIR)
	coverage_id=$(echo $coverage_message | grep -oh "[0-9]*$")
else
	# Run the sort/index array job
	sortindex_message=$(sbatch --array=0-$(($num_si_arr - 1)) ${curDir}/03.1_sortIndexArray.sh SICC.tmp)
	sortindex_id=$(echo $sortindex_message | grep -oh "[0-9]*$")

	# Detect how big our coverage job array should be
	num_coverage_arr=$num_si_arr
	# Schedule the coverage array job to run after the sort/indexing is done.
	coverage_message=$(sbatch --depend=afterok:${sortindex_id} --array=0-$((num_coverage_arr - 1)) ${curDir}/03.2_bamCoverageArray.sh SICC_sorted.tmp $TMPDIR)
	coverage_id=$(echo $coverage_message | grep -oh "[0-9]*$")
fi

# If there is no comparisons file input, quit.
if [[ -z "$3" ]]; then
  exit
fi

# Detect how big our comparison job array should be (this is a little overcomplicated but it's fine)
num_compare_arr=0
while IFS=',' read -r -a line; do
    for ((i=0; i < ${#line[@]}-1; i++)); do
	    num_compare_arr=$((num_compare_arr + 1))
    done
done < $COMP

mkdir bw_compare
sbatch --depend=afterok:${coverage_id} --array=0-$((num_compare_arr - 1)) ${curDir}/03.3_bwCompareArray.sh . $COMP

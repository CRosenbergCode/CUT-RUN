#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --job-name=bigwig_SICC
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


# =============================== 03_SICC.sh =================================
#
#  DESCRIPTION: A (totally SICC) script that runs samtools sort and index for
#          each .bam in the specified directory, then runs bamCoverage and
#          bigwigCompare for each. (SICC = Sort, Index, Coverage, Compare)
#          This process takes several hours, so it should be run as a slurm job
#          on Riviera with sbatch.
#
#  USAGE: sbatch 03_SICC.sh $1 $2
#
#      Argument 1 - Path to the directory containing your .bam files.
#      Argument 2 - Path to the text file containing comparison data. If not
#                    specified, it will not run the comparison. Comparison
#                    data is written using the following rules:
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
#
#  OUTPUT: Sorted .bam files, bam index (.bai) files, bigwig (.bw) files for
#          each .bam, and (optionally) the compared bigwig files.
#
# ============================================================================

source activate base
conda init
conda activate deeptools_kernel

# Setting a custom temp directory. DeepTools will otherwise default to Riviera's /tmp
# directory, which is teeny (~2GB) and I think also shared between everyone.
mkdir dt_tmp
pwdstring=$(pwd)
export TMPDIR="${pwdstring}/dt_tmp"

# If a comparisons file is given, make sure it's formatted properly.
if [[ -n "$2" ]]; then
  # Ensure file is not DOS-formatted, because MobaRTE is actually awful
  if cat -v $2 | grep "\^M" ; then
    echo -e "DOS Line endings found in your comparisons file: ${2}\nIf you're using MobaXTerm, click the little penguin button\nwhen editing your comparisons file in the text editor and save your file!"
    echo "Exiting. Please fix your comparisons file."
    exit
  fi
  
  # Get global path for later.
  COMP=$(readlink -f $2)
  
  # Make sure it ends in a newline; if not then add one. (POSIX formatting)
  if [ -n "$(tail -c1 $2)" ]; then
    echo "" >> $2;
  fi
fi

# Change into relevant directory.
cd $1

# Sort and index bams
for bam in *.bam; do
	samtools sort $bam -o sorted_$bam
	samtools index sorted_$bam
done

# Run bamCoverage for all sorted bams.
# 4 processors may or may not be optimal, but the original scripts used nthreads=4.
for FILE in sorted_*.bam
do
  prefnobam=${FILE::-4}
  pref=${prefnobam:7}
  echo ${pref}
  bamCoverage --bam ${FILE} -o ${pref}.bw \
      --binSize 10 \
      --normalizeUsing RPGC \
      --effectiveGenomeSize 12200000000 \
      --extendReads \
      --numberOfProcessors 4
done

# If there is no second argument, quit.
if [[ -z "$2" ]]; then
  exit
fi

# Run bigwigCompare for all .bw files, put them in the new bw_compare directory.
# This only works for the naming convention for the merged bams I was given. When we
# establish a naming convention in the future, this will likely need to be updated.

mkdir bw_compare

# Loop over each line, run bigwigCompare for each treatment/input pair per line.
while IFS=',' read -r -a line; do
    for ((i=0; i < ${#line[@]}-1; i++)); do
        treat=${line[i]}
        input=${line[-1]}
        treat=${treat/.bam/.bw}
        input=${input/.bam/.bw}
        bigwigCompare -b1 $treat -b2 $input -o ./bw_compare/compare_$treat --binSize 10 --numberOfProcessors 4
    done
done < $COMP
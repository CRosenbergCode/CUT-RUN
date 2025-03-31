#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=BamCoverage
#SBATCH --output=%x.%j.out
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4

# ========================= 03.1_runBamCoverage.sh ===========================
#
#  DESCRIPTION: Creates bigwigs of all bams in the current directory
#
#  USAGE: sbatch 03.1_runBamCoverage.sh
#
#  OUTPUT: Bigwig (.bw) files in the current directory. They will
#          be used in 03.2_runBWCompare.sh.
#
# ============================================================================

source $HOME/miniconda3/bin/activate base
conda init
conda activate deeptools_kernel

# Setting a custom temp directory. This is needed so deepTools doesn't use
# the default /tmp directory on Riviera, which is tiny. (~2GB)
mkdir dt_tmp
pwdstring=$(pwd)
export TMPDIR="${pwdstring}/dt_tmp"

for FILE in *.bam
do
  pref=${FILE::-4}
  #pref=${pref:21}
  echo ${pref}
  bamCoverage --bam ${FILE} -o ${pref}.bw \
      --binSize 10 \
      --normalizeUsing RPGC \
      --effectiveGenomeSize 12200000000 \
      --extendReads \
      --numberOfProcessors 4
done
#Adjust --binSize to adjust how granular the result bigwig file is, higher is lower resolution
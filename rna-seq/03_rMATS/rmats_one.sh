#!/bin/bash -l

#SBATCH --job-name=rmats
#SBATCH --output=%x.%j.log # gives jobname.ID.log
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --partition=short-cpu

#=========================== rmats_one.sh ==============================
#
# DESCRIPTION: runs rMATS once, comparing two sets of bam files.
# 
# USAGE: rmats_one.sh $1 $2 $3
#      Argument 1 - Path to text file listing paths to first set of bam files. (comma-delimited)
#      Argument 2 - Path to text file listing paths to second set of bam files. (comma-delimited)
#      Argument 3 - Name of this current run. Will be used for output directory names.
#
#
# Notes:
#      This script only works in my directory on Riviera. I was running into some issues
#      with rMATS being picky about absolute paths, but I'm less bad at this now, so I
#      should be able to fix this and make it more generalizable. For now, be sure to
#      rewrite paths for your own conda installations, GTF file, etc.
#
#=======================================================================

source activate base
conda init
conda activate rMATS

rmats.py --b1 $1 --b2 $2 --gtf ./VectorBase-68_AaegyptiLVP_AGWG.gtf --od ./rmats_${3}_out/ -t paired --readLength 150 --tmp ./rmats_${3}_tmp/ --nthread 16

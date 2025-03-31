#!/bin/bash -l

#SBATCH --job-name=makematrix
#SBATCH --output=%x.%j.log
#SBATCH -N 1
#SBATCH -n 6 # The number of processors used in the example
#SBATCH --partition=short-cpu

# INSTRUCTIONS:
#
# Create a text file with the paths to each of the bigwig files
# you want to include in the matrix. These are the groups you want
# to graph together in heatmaps and profiles.
# THE PATHS SHOULD BE SEPARATED BY SPACES, NOT NEWLINES.
#
# Be sure to name the file something descriptive, as it will be used
# in naming the resulting region bedfile and matrix .gz file.
#
# WHAT TO CHANGE:
#
#   Change the --referencePoint option to focus on different parts of the genes. (options are TSS, TES, and center)
#   Change the --samplesLabel option to properly label your different groups. This can also be done in the make_profile and heatmap scripts.
#   Change the -a and -b options to change the number of bases after (a) and before (b) the reference point.
#   Change the -p option if you want a different number of processors. (be sure to update #SBATCH -n at the top of the script too)
#
# When you are ready to run it, run the following command:
#
#    sbatch make_matrices.sh path/to/bigwiglist.txt
#
# Your outputs will be in the directory you ran this in.

source activate base
conda init
conda activate deeptools_kernel

# Grabbing everything after the last forward slash...
filename=$(echo $1 | sed 's:.*/::')
# ...and before the last period.
filename=${filename%.*}

computeMatrix reference-point --referencePoint TSS \ 
  -b 2000 -a 2000 \
  -R ./Aedes68.Bed \ 
  -S $(cat $1) \ 
  --skipZeros \ 
  -o ./${filename}.gz \
  -p 6 \ 
  --outFileSortedRegions ./region_${filename}.bed
#!/bin/bash -l

#SBATCH --job-name=makeprofile
#SBATCH --output=%x.%j.log
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition=short-cpu

# Instructions:
#
# You will need to edit these options a lot according to your needs.
# Read the documentation for different options: https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#
# This will create a profile of your matrix's bigwigs, showing the input
# read density around the transcription start site. Currently, this is
# set up for BF vs RVFV. Here is what you need to change to create different plots:
# 
# Change the --colors option to include the desired colors for all your groups. If you have three groups, you will need to write 3 colors.
# Change the --refPointLabel option if you are focusing on something other than the TSS. (see make_matrices.sh for how to change the focus)
# Change the --samplesLabel option to properly label your different groups. This can also be done in make_matrices.sh.
# 
# The -T option sets the title, and is set to the second argument, since you may
# run the script multiple times with the same settings, but for different matrices.
#
# How to run (slurm recommended since it takes a minute otherwise):
#
#    sbatch make_profile.sh path/to/matrix.gz "Your title, in quotes"
#
# Your outputs will be in a new folder called figures.

source activate base
conda init
conda activate deeptools_kernel

# Grabbing everything after the last forward slash...
filename=$(echo $1 | sed 's:.*/::')
# ...and before the last period.
filename=${filename%.*}

title=$2

mkdir figures

plotProfile -m $1 -out ./figures/${filename}_profile.png \ 
    --perGroup \ 
    --colors blue red \ 
    --refPointLabel "TSS" \ 
    -T "$title" \ 
    --samplesLabel "Bloodfed" "RVFV" \ 
    -y "Reads per Genomic Content" \ 
    --regionsLabel " "
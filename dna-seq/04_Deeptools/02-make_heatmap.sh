#!/bin/bash -l

#SBATCH --job-name=makeheatmap
#SBATCH --output=%x.%j.log
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition=short-cpu

# Instructions:
#
# You will need to edit these options a lot according to your needs.
# Read the documentation for different options: https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html
#
# Change the --colorMap option if you want different coloring. (check the docs for different options)
# Change the --whatToShow option to change what elements are displayed. (this can include the plot, heatmap, and/or color legend. Check docs)
# Change the --refPointLabel option if you are focusing on something other than the TSS. (see make_matrices.sh for how to change the focus)
# Change the --samplesLabel option to properly label your different groups. This can also be done in make_matrices.sh.
# Change the --zMin and --zMax values to specify the min and max cutoffs displayed by the colors in the heatmap.
# 
# The -T option sets the title, and is set to the second argument, since you may
# run the script multiple times with the same settings, but for different matrices.
#
# How to run (slurm recommended since it takes a minute otherwise):
#
#    sbatch make_heatmap.sh path/to/matrix.gz "Your title, in quotes"
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

plotHeatmap -m $1 -out ./figures/${filename}_heatmap.png \
	--colorMap RdBu \
	--whatToShow 'heatmap and colorbar' \
	--refPointLabel "TSS" \
	--samplesLabel "Bloodfed" "RVFV" \
	-T "$title" \
	--zMin "-1" --zMax "3"
#	If you wanted automatic z-scaling:
#	--zMin auto --zMax auto
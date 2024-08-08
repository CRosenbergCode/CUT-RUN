#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=MACS2PeakCall
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
source activate base
conda activate rnaPseudo

#-t flag refers
#The --broad flag should generally be used for H3K9me3, while H3K27ac and TFs do not need a flag
#-n refers to the name of the directory which will contain the output files
#--slocal refers to the
#'-g' refers to the effective genome size, approximately 1.3e+9 in Aedes aegypti.
#

#'norm' should generally be used for experimetns
#'non' can be used if data is already adjusted via spike-in dna

#Prior to this, bedgraphs must have been created using the '' script
#The first parameter is the bedgraph to call peaks for

#The second parameter should generally be a bedgraph of the corresponding input for the first sample
#If not using input normalization, a number corresponding to the top percentage of peaks can be used (i.e. 0.01 to call the top 1% of peaks)
#This is generally not advised unless a robust blacklisting method has already been used

#'relaxed' allows for more peaks, calls the

SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output

SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output

SEACR_1.3.sh target.bedgraph 0.01 non stringent output

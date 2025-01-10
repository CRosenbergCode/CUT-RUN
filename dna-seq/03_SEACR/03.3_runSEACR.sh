#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=SEACRAgoD1
#SBATCH --output=%x.%j.out
#SBATCH --time=4:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=8
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

bash SEACR_1.3.sh trimmed_BF-d1-rep1AGo2_R2_001.fragments.bedgraph trimmed_BF-d1-rep1input_R2_001.fragments.bedgraph norm relaxed BF_Ago2_D1_Rep1.relaxed.bed

bash SEACR_1.3.sh trimmed_BF-d1-rep2AGo2_R2_001.fragments.bedgraph trimmed_BF-d1-rep2-input_R2_001.fragments.bedgraph norm relaxed BF_Ago2_D1_Rep2.relaxed.bed

bash SEACR_1.3.sh trimmed_SF-d1-rep1AGo2_R2_001.fragments.bedgraph trimmed_SF-d1-rep1-input_R2_001.fragments.bedgraph norm relaxed SF_Ago2_D1_Rep1.relaxed.bed

bash SEACR_1.3.sh trimmed_SF-d1-rep2AGo2_R2_001.fragments.bedgraph trimmed_SF-d1-rep2-input_R2_001.fragments.bedgraph norm relaxed SF_Ago2_D1_Rep2.relaxed.bed

bash SEACR_1.3.sh trimmed_SF-d1-rep3AGo2_R2_001.fragments.bedgraph trimmed_SF-d1-rep3-input_R2_001.fragments.bedgraph norm relaxed SF_Ago2_D1_Rep3.relaxed.bed

bash SEACR_1.3.sh trimmed_BF-d1-rep4-Ago2-60_R2_001.fragments.bedgraph trimmed_BF-d1-rep3-input_R2_001.fragments.bedgraph norm relaxed BF_Ago2_D1_Rep4.relaxed.bed

bash SEACR_1.3.sh trimmed_SF-d1-rep4-Ago2-60_R2_001.fragments.bedgraph trimmed_SF-d1-rep3-input_R2_001.fragments.bedgraph norm relaxed SF_Ago2_D1_Rep4.relaxed.bed





bash SEACR_1.3.sh trimmed_BF-d1-rep1AGo2_R2_001.fragments.bedgraph trimmed_BF-d1-rep1input_R2_001.fragments.bedgraph norm stringent BF_Ago2_D1_Rep1.stringent.bed

bash SEACR_1.3.sh trimmed_BF-d1-rep2AGo2_R2_001.fragments.bedgraph trimmed_BF-d1-rep2-input_R2_001.fragments.bedgraph norm stringent BF_Ago2_D1_Rep2.stringent.bed

bash SEACR_1.3.sh trimmed_SF-d1-rep1AGo2_R2_001.fragments.bedgraph trimmed_SF-d1-rep1-input_R2_001.fragments.bedgraph norm stringent SF_Ago2_D1_Rep1.stringent.bed

bash SEACR_1.3.sh trimmed_SF-d1-rep2AGo2_R2_001.fragments.bedgraph trimmed_SF-d1-rep2-input_R2_001.fragments.bedgraph norm stringent SF_Ago2_D1_Rep2.stringent.bed

bash SEACR_1.3.sh trimmed_SF-d1-rep3AGo2_R2_001.fragments.bedgraph trimmed_SF-d1-rep3-input_R2_001.fragments.bedgraph norm stringent SF_Ago2_D1_Rep3.stringent.bed

bash SEACR_1.3.sh trimmed_BF-d1-rep4-Ago2-60_R2_001.fragments.bedgraph trimmed_BF-d1-rep3-input_R2_001.fragments.bedgraph norm stringent BF_Ago2_D1_Rep4.stringent.bed

bash SEACR_1.3.sh trimmed_SF-d1-rep4-Ago2-60_R2_001.fragments.bedgraph trimmed_SF-d1-rep3-input_R2_001.fragments.bedgraph norm stringent SF_Ago2_D1_Rep4.stringent.bed

#SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output

#SEACR_1.3.sh target.bedgraph 0.01 non stringent output

#trimmed_BF-d1-rep1AGo2_R2_001.paired.bam
#03.3_runSEACR.sh             trimmed_BF-d1-rep1AGo2_R2_001.fragments.bedgraph
#genome.fa                    trimmed_BF-d1-rep1input_R2_001.fragments.bedgraph
#SEACR_1.3.R                  trimmed_BF-d1-rep2AGo2_R2_001.fragments.bedgraph
#SEACR_1.3.sh                 trimmed_SF-d1-rep1AGo2_R2_001.fragments.bedgraph
#SEACR_Directions.md          trimmed_SF-d1-rep1-input_R2_001.fragments.bedgraph
#SEACRMakeBedgraph.31608.out  trimmed_SF-d1-rep2AGo2_R2_001.fragments.bedgraph
#SEACRMakeBedgraph.31609.out  trimmed_SF-d1-rep3AGo2_R2_001.fragments.bedgraph

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

#'norm' should generally be used for experimetns
#'non' can be used if data is already adjusted via spike-in dna

#Prior to this, bedgraphs must have been created using the '03.2_makeBedGraph' script
#The first parameter is the bedgraph to call peaks for, generally the sample with the antibody of interest (Ex: BF_Ago2_D1_Rep1)

#The second parameter should generally be a bedgraph of the corresponding input for the first sample (Ex: BF_Input_D1_Rep1)
#If not using input normalization, a number corresponding to the top percentage of peaks can be used (i.e. 0.01 to call the top 1% of peaks)
#This is generally not advised unless a robust blacklisting method has already been used

#'relaxed' allows for more peaks
#'stringent' provides a smaller set of peaks but with higher confidence

#The final parameter is the prefix of the output file
#.*.bed will be added as a suffix, where * is the stringency of the call

#Examples Provided below

#General Format

#bash SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output

#bash SEACR_1.3.sh target.bedgraph 0.01 non stringent output

#Call input adjusted peaks with the relaxed thresholds
bash SEACR_1.3.sh trimmed_BF-d1-rep1AGo2_R2_001.fragments.bedgraph trimmed_BF-d1-rep1input_R2_001.fragments.bedgraph norm relaxed BF_Ago2_D1_Rep1

#Call input adjusted peaks with the stringent threshold
bash SEACR_1.3.sh trimmed_BF-d1-rep1AGo2_R2_001.fragments.bedgraph trimmed_BF-d1-rep1input_R2_001.fragments.bedgraph norm stringent BF_Ago2_D1_Rep1

#Call non-adjusted peaks and return results which represent the top 0.01 percent of peaks found
bash SEACR_1.3.sh trimmed_BF-d1-rep2AGo2_R2_001.fragments.bedgraph 0.01 non stringent BF_Ago2_D1_Rep2

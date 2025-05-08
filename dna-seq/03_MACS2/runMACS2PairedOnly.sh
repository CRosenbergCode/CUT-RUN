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

macs2 callpeak -t hisatOnlyPaired/BF_Ac_Rep1_S10_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatOnlyPaired/Input_BF_rep1_S2_L006_.bam \
  -n BF_Ac_1 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/BF_Ac_1

macs2 callpeak -t hisatOnlyPaired/BF_Ac_Rep2_S16_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatOnlyPaired/Input_BF_rep2_S8_L006_.bam \
  -n BF_Ac_2 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/BF_Ac_2

macs2 callpeak -t hisatOnlyPaired/BF_Me_Rep1_S11_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatOnlyPaired/Input_BF_rep1_S2_L006_.bam \
  -n BF_Me_1 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/BF_Me_1

macs2 callpeak -t hisatOnlyPaired/BF_Me_Rep2_S13_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatOnlyPaired/Input_BF_rep2_S8_L006_.bam \
  -n BF_Me_2 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/BF_Me_2

macs2 callpeak -t hisatOnlyPaired/BF_neg_Rep1_S12_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatOnlyPaired/Input_BF_rep1_S2_L006_.bam \
  -n BF_Neg_1 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/BF_Neg_1

macs2 callpeak -t hisatOnlyPaired/BF_neg_Rep2_S6_L006_.bam  \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatOnlyPaired/Input_BF_rep2_S8_L006_.bam \
  -n BF_Neg_2 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/BF_Neg_2


macs2 callpeak -t hisatOnlyPaired/RVFV_Ac_Rep1-1_S9_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatOnlyPaired/Input_RVFV_rep1_S1_L006_.bam \
  -n RVFV_Ac_1 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/RVFV_Ac_1_1

macs2 callpeak -t hisatOnlyPaired/RVFV_Ac_Rep1-2_S17_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatOnlyPaired/Input_RVFV_rep1_S1_L006_.bam \
  -n RVFV_Ac_1 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/RVFV_Ac_1_2

macs2 callpeak -t hisatOnlyPaired/RVFV_Ac_Rep2_S14_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatOnlyPaired/Input_RVFV_rep2_S7_L006_.bam \
  -n RVFV_Ac_2 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/RVFV_Ac_2

macs2 callpeak -t hisatOnlyPaired/RVFV_Me_Rep1_S3_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatOnlyPaired/Input_RVFV_rep1_S1_L006_.bam \
  -n RVFV_Me_1 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/RVFV_Me_1

macs2 callpeak -t hisatOnlyPaired/RVFV_Me_Rep2_S5_L006_.bamm \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatOnlyPaired/Input_RVFV_rep2_S7_L006_.bam \
  -n RVFV_Me_2 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/RVFV_Me_2


macs2 callpeak -t hisatOnlyPaired/RVFV_neg_Rep1_S4_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatOnlyPaired/Input_RVFV_rep1_S1_L006_.bam \
  -n RVFV_Neg_1 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/RVFV_Neg_1

macs2 callpeak -t hisatOnlyPaired/RVFV_neg_Rep2_S15_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatOnlyPaired/Input_RVFV_rep2_S7_L006_.bam \
  -n RVFV_Neg_2 \
  --keep-dup all \
  --slocal 2000 \
  --outdir MACS2Peaks/RVFV_Neg_2


#Broad is H3K9me3
#Narrow is H3K27ac

#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=BWCompare
#SBATCH --output=%x.%j.out
#SBATCH --time=3:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
module load anaconda
conda activate deeptools_kernel

bigwigCompare -b1 BF_Ac_Rep1_S10_L006_.bw -b2 Input_BF_rep1_S2_L006_.bw -o adBF_Ac_1.bw
bigwigCompare -b1 BF_Me_Rep1_S11_L006_.bw -b2 Input_BF_rep1_S2_L006_.bw -o adBF_Me_1.bw
bigwigCompare -b1 BF_neg_Rep1_S12_L006_.bw -b2 Input_BF_rep1_S2_L006_.bw -o adBF_Neg_1.bw

bigwigCompare -b1 BF_Ac_Rep2_S16_L006_.bw -b2 Input_BF_rep2_S8_L006_.bw -o adBF_Ac_2.bw
bigwigCompare -b1 BF_Me_Rep2_S13_L006_.bw -b2 Input_BF_rep2_S8_L006_.bw -o adBF_Me_2.bw
bigwigCompare -b1 BF_neg_Rep2_S6_L006_.bw  -b2 Input_BF_rep2_S8_L006_.bw -o adBF_Neg_2.bw

bigwigCompare -b1 RVFV_Ac_Rep1-1_S9_L006_.bw -b2 Input_RVFV_rep1_S1_L006_.bw -o adRVFV_Ac_1_1.bw
bigwigCompare -b1 RVFV_Ac_Rep1-2_S17_L006_.bw -b2 Input_RVFV_rep1_S1_L006_.bw -o adRVFV_Ac_1_2.bw
bigwigCompare -b1 RVFV_Me_Rep1_S3_L006_.bw -b2 Input_RVFV_rep1_S1_L006_.bw -o adRVFV_Me_1.bw
bigwigCompare -b1 RVFV_neg_Rep1_S4_L006_.bw -b2 Input_RVFV_rep1_S1_L006_.bw -o adRVFV_Neg_1.bw

bigwigCompare -b1 RVFV_Ac_Rep2_S14_L006_.bw -b2 Input_RVFV_rep2_S7_L006_.bw-o adRVFV_Ac_2.bw
bigwigCompare -b1 RVFV_Me_Rep2_S5_L006_.bw -b2 Input_RVFV_rep2_S7_L006_.bw -o adRVFV_Me_2.bw
bigwigCompare -b1 RVFV_neg_Rep2_S15_L006_.bw -b2 Input_RVFV_rep2_S7_L006_.bw -o adRVFV_Neg_2.bw

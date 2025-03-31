#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=BWCompare
#SBATCH --output=%x.%j.out
#SBATCH --time=3:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4

# ========================== 03.2_runBWCompare.sh ============================
#
#  DESCRIPTION: Creates compared bigwigs of each specified treatment/control
#                pair. You will have to edit this script yourself. -b1 will
#                be your treatment bigwig, and -b2 will be the input. (control)
#
#  USAGE: sbatch 03.2_runBWCompare.sh
#
#  OUTPUT: Bigwig (.bw) files in a directory called "bw_compare". They will
#            be comparisons (by default, Log2 ratios) of the signals described
#            by the 03.1 bigwigs.
#
# ===========================================================================

source $HOME/miniconda3/bin/activate base
conda init
conda activate deeptools_kernel

# Setting a custom temp directory. This is needed so deepTools doesn't use
# the default /tmp directory on Riviera, which is tiny. (~2GB)
mkdir dt_tmp
pwdstring=$(pwd)
export TMPDIR="${pwdstring}/dt_tmp"

mkdir bw_compare

bigwigCompare -b1 BF_Ac_Rep1_S10_L006_.bw -b2 Input_BF_rep1_S2_L006_.bw -o adBF_Ac_1.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 BF_Me_Rep1_S11_L006_.bw -b2 Input_BF_rep1_S2_L006_.bw -o adBF_Me_1.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 BF_neg_Rep1_S12_L006_.bw -b2 Input_BF_rep1_S2_L006_.bw -o adBF_Neg_1.bw --binSize 10 --numberOfProcessors 4

bigwigCompare -b1 BF_Ac_Rep2_S16_L006_.bw -b2 Input_BF_rep2_S8_L006_.bw -o adBF_Ac_2.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 BF_Me_Rep2_S13_L006_.bw -b2 Input_BF_rep2_S8_L006_.bw -o adBF_Me_2.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 BF_neg_Rep2_S6_L006_.bw  -b2 Input_BF_rep2_S8_L006_.bw -o adBF_Neg_2.bw --binSize 10 --numberOfProcessors 4

bigwigCompare -b1 RVFV_Ac_Rep1-1_S9_L006_.bw -b2 Input_RVFV_rep1_S1_L006_.bw -o adRVFV_Ac_1_1.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 RVFV_Ac_Rep1-2_S17_L006_.bw -b2 Input_RVFV_rep1_S1_L006_.bw -o adRVFV_Ac_1_2.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 RVFV_Me_Rep1_S3_L006_.bw -b2 Input_RVFV_rep1_S1_L006_.bw -o adRVFV_Me_1.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 RVFV_neg_Rep1_S4_L006_.bw -b2 Input_RVFV_rep1_S1_L006_.bw -o adRVFV_Neg_1.bw --binSize 10 --numberOfProcessors 4

bigwigCompare -b1 RVFV_Ac_Rep2_S14_L006_.bw -b2 Input_RVFV_rep2_S7_L006_.bw-o adRVFV_Ac_2.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 RVFV_Me_Rep2_S5_L006_.bw -b2 Input_RVFV_rep2_S7_L006_.bw -o adRVFV_Me_2.bw --binSize 10 --numberOfProcessors 4
bigwigCompare -b1 RVFV_neg_Rep2_S15_L006_.bw -b2 Input_RVFV_rep2_S7_L006_.bw -o adRVFV_Neg_2.bw --binSize 10 --numberOfProcessors 4
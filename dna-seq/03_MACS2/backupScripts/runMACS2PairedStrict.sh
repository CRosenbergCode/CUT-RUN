#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=MACS2PairedStrict
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hunter.ogg@colostate.edu

module purge
module load anaconda
conda activate rnaPseudo

for FILE in hisatPairStrict/*.sam
do
  EXT=${FILE::-4}
  echo $EXT
  samtools sort $FILE -o $EXT.bam -@ 6
done

macs2 callpeak -t hisatPairStrict/BF_Ac_Rep1_S10_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatPairStrict/Input_BF_rep1_S2_L006_.bam \
  -n BF_Ac_1 \
  --slocal 2000 \
  --outdir macs2PairedStrict/BF_Ac_1

macs2 callpeak -t hisatPairStrict/BF_Ac_Rep2_S16_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatPairStrict/Input_BF_rep2_S8_L006_.bam \
  -n BF_Ac_2 \
  --slocal 2000 \
  --outdir macs2PairedStrict/BF_Ac_2

macs2 callpeak -t hisatPairStrict/BF_Me_Rep1_S11_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatPairStrict/Input_BF_rep1_S2_L006_.bam \
  -n BF_Me_1 \
  --slocal 2000 \
  --outdir macs2PairedStrict/BF_Me_1

macs2 callpeak -t hisatPairStrict/BF_Me_Rep2_S13_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatPairStrict/Input_BF_rep2_S8_L006_.bam \
  -n BF_Me_2 \
  --slocal 2000 \
  --outdir macs2PairedStrict/BF_Me_2

macs2 callpeak -t hisatPairStrict/BF_neg_Rep1_S12_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatPairStrict/Input_BF_rep1_S2_L006_.bam \
  -n BF_Neg_1 \
  --slocal 2000 \
  --outdir macs2PairedStrict/BF_Neg_1

macs2 callpeak -t hisatPairStrict/BF_neg_Rep2_S6_L006_.bam  \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatPairStrict/Input_BF_rep2_S8_L006_.bam \
  -n BF_Neg_2 \
  --slocal 2000 \
  --outdir macs2PairedStrict/BF_Neg_2


macs2 callpeak -t hisatPairStrict/RVFV_Ac_Rep1-1_S9_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatPairStrict/Input_RVFV_rep1_S1_L006_.bam \
  -n RVFV_Ac_1 \
  --slocal 2000 \
  --outdir macs2PairedStrict/RVFV_Ac_1_1

macs2 callpeak -t hisatPairStrict/RVFV_Ac_Rep1-2_S17_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatPairStrict/Input_RVFV_rep1_S1_L006_.bam \
  -n RVFV_Ac_1 \
  --slocal 2000 \
  --outdir macs2PairedStrict/RVFV_Ac_1_2

macs2 callpeak -t hisatPairStrict/RVFV_Ac_Rep2_S14_L006_.bam \
  -f BAMPE -g 1.3e+9 \
  -c hisatPairStrict/Input_RVFV_rep2_S7_L006_.bam \
  -n RVFV_Ac_2 \
  --slocal 2000 \
  --outdir macs2PairedStrict/RVFV_Ac_2

macs2 callpeak -t hisatPairStrict/RVFV_Me_Rep1_S3_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatPairStrict/Input_RVFV_rep1_S1_L006_.bam \
  -n RVFV_Me_1 \
  --slocal 2000 \
  --outdir macs2PairedStrict/RVFV_Me_1

macs2 callpeak -t hisatPairStrict/RVFV_Me_Rep2_S5_L006_.bamm \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatPairStrict/Input_RVFV_rep2_S7_L006_.bam \
  -n RVFV_Me_2 \
  --slocal 2000 \
  --outdir macs2PairedStrict/RVFV_Me_2



macs2 callpeak -t hisatPairStrict/RVFV_neg_Rep1_S4_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatPairStrict/Input_RVFV_rep1_S1_L006_.bam \
  -n RVFV_Neg_1 \
  --slocal 2000 \
  --outdir macs2PairedStrict/RVFV_Neg_1

macs2 callpeak -t hisatPairStrict/RVFV_neg_Rep2_S15_L006_.bam \
  -f BAMPE -g 1.3e+9 --broad \
  -c hisatPairStrict/Input_RVFV_rep2_S7_L006_.bam \
  -n RVFV_Neg_2 \
  --slocal 2000 \
  --outdir macs2PairedStrict/RVFV_Neg_2


#Broad is H3K9me3
#Narrow is H3K27ac

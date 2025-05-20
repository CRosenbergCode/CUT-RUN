#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=ParamComparesMACSAgo2
#SBATCH --output=%x.%j.out
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log

module purge
source activate base
conda activate rnaPseudo

#-t flag refers
#The --broad flag should generally be used for H3K9me3, while H3K27ac and TFs do not need a flag
#-n refers to the name of the directory which will contain the output files
#--slocal refers to the
#'-g' refers to the effective genome size, approximately 1.3e+9 in Aedes aegypti.
#

echo "[$0] $SLURM_JOB_NAME $@" # log the command line

#FE
#Q
#slocal
#llocal

declare -a dmin=(20 30 40 50)
declare -a qvals=(0.01 0.005 0.00001)
declare -a slocals=(750 1000 2000 5000)
declare -a llocals=(5000 7500 10000 20000 50000)


## now loop through the above array
linenum=0
for d in "${dmin[@]}"
do
   for q in "${qvals[@]}"
   do
      for s in "${slocals[@]}"
      do
         for l in "${llocals[@]}"
         do
           if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
           then
             macs2 callpeak -t BF_Ago2_Merged_D7_Aae.bam \
               -f BAMPE -g 1.3e+9 \
               -n BF_Ago2_D7_d${d}_l${l}_s${s}_q${q}_lamba \
               -c BF_Input_Merged_D7_Aae.bam \
               --d-min $d \
               --llocal $l \
               --slocal $s \
               -q $q \
               --outdir ParamTests

             macs2 callpeak -t SF_Ago2_Merged_D7_Aae.bam \
               -f BAMPE -g 1.3e+9 \
               -n SF_Ago2_D7_d${d}_l${l}_s${s}_q${q}_lamba \
               -c SF_Input_Merged_D7_Aae.bam \
               --d-min $d \
               --llocal $l \
               --slocal $s \
               -q $q \
               --outdir ParamTests
           fi
           linenum=$((linenum + 1))
         done
      done
   done
done

#--cutoff-analysis

#Broad is H3K9me3
#Narrow is H3K27ac

#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=samtoolsMergeAutomated
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


#echo "[$0] $SLURM_JOB_NAME $@" # log the command line


#declare -a conditions=("BF" "SF" "RVFV")

#find /folder1/folder2/ -name file.txt0*


#If identical input does not exist, will use one with different rep number

nthreads=$(($SLURM_NTASKS - 1))

outdir=""

merged=""

ls -l | head -1

for FILE in *.bam; do
  if [[ $merged =~ $FILE ]]; then
    continue
  fi
  mergedName="${FILE%%Rep*}Merged_${FILE##*Rep*[0-9]_}"
  echo -e "Merged Name: $mergedName"
  matches=$(ls ${FILE%%Rep*}*${FILE##*Rep*[0-9]_})
  echo -e "Matches: $matches \n"
  merged="${merged} ${matches}"
  samtools merge ${outdir}$mergedName $matches -@$nthreads
  #echo -e "Called peaks for $FILE using $input as input control and saved results with prefix ${nameOnly}_Control\n"
done

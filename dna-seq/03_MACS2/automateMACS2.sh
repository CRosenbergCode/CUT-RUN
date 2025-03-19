#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=MACS2Automated
#SBATCH --output=%x.%j.out
#SBATCH --time=8:00:00
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


echo "[$0] $SLURM_JOB_NAME $@" # log the command line


declare -a antis=("Ago2" "H3K9Me3" "H3K27Ac" "Neg")
#declare -a conditions=("BF" "SF" "RVFV")

#find /folder1/folder2/ -name file.txt0*


#If identical input does not exist, will use one with different rep number

outdir="MACS2Peaks"

match_input() {
    matched=$(echo "$1" | sed "s/-[0-9]\+[A-Z]*//")
    for anti in "${antis[@]}"
    do
      matched="${matched/$anti/"Input"}"
    done
    if [ ! -f $matched ]; then
        echo "Cannot find $matched for input control"
        matched=$(ls ${matched%%Rep*}*${matched##*Rep*[0-9]_}* | head -1)
        echo -e "Used file $matched instead\n"
    fi
    echo $matched
}

ls -l | head -1

for FILE in *.bam; do
  nameOnly=$(basename $FILE .bam)
  if [[ $FILE =~ "Input" ]]; then
     echo -e "Skipped $FILE\n"
     continue
  fi
  input=$(match_input "$FILE")
  if [[ $FILE =~ "H3K9Me3" ]]; then
    macs2 callpeak -t $FILE \
      -f BAMPE -g 1.3e+9 \
      -c $input \
      -n ${nameOnly}_Control \
      --keep-dup all \
      --broad \
      --outdir $outdir
      echo "Using broad parameter"
  else
    macs2 callpeak -t $FILE \
      -f BAMPE -g 1.3e+9 \
      -c $input \
      -n ${nameOnly}_Control \
      --keep-dup all \
      --outdir $outdir
  fi

  echo -e "Called peaks for $FILE using $input as input control and saved results with prefix ${nameOnly}_Control\n"
done

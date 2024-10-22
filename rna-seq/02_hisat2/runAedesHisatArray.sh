#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=100
#SBATCH --time=96:00:00
#SBATCH --qos=normal
#SBATCH --partition=week-long-cpu
#SBATCH --job-name=AedesArrayHisat
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log

echo "[$0] $SLURM_JOB_NAME $@" # log the command line

module purge
source activate base
conda activate rnaPseudo

date # timestamp

filename=AedesTrimmed.txt # $1

#Run this first to generate a list of files.
#ls | grep -v / | grep ".fastq.gz" | grep "R1" | grep "trimmed" >> AedesTrimmed.txt

# to run this script use
#sbatch --array=0-<N> <script.sh>
#Where N is the number of lines in your AedesTrimmed.txt file minus 1. You can get that easily using the following command
#wc -l AedesTrimmed.txt
#sbatch --array=0-25 runAedesHisatArray.sh

linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then

      TWO=${line/_R1_/_R2_}
      prefix="../01_Fastp/DeDupTrimmed/DeDupTrimmed."
      OUT=${TWO#$prefix}
      echo ${TWO}
      echo $OUT
      hisat2 --phred33 --rna-strandness RF -p 100 -x HisatIndex/AedesTranscriptome68 -1 ${line} -2 ${TWO} -S HisatAligned/$OUT.sam
    fi
    linenum=$((linenum + 1))
done < $filename

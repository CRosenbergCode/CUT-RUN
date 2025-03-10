#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=128:00:00
#SBATCH --qos=normal
#SBATCH --partition=week-long-cpu
#SBATCH --job-name=GoArray
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log

echo "[$0] $SLURM_JOB_NAME $@" # log the command line

module purge
source activate base
conda activate rnaPseudo

date # timestamp

filename=Segments.txt # $1

linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
      Rscript GoArray.R $line
    fi
    linenum=$((linenum + 1))
done < $filename

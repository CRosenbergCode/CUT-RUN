#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=04:00:00
#SBATCH --qos=normal
#SBATCH --partition=day-long-cpu
#SBATCH --job-name=concordanthscountArrayAedes
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log

echo "[$0] $SLURM_JOB_NAME $@" # log the command line

module purge
source activate base
conda activate rnaPseudo

date # timestamp

filename=AedesAligned.txt # $1

linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
      pref=${line/..\/02_Hisat2\/HisatAligned\//hsCountsAedes\/}
      suff=${pref/.bam/.HSCounts.txt}
      echo ${pref}
      echo ${suff}
      htseq-count --stranded=reverse $line 68genes.gtf > $suff
    fi
    linenum=$((linenum + 1))
done < $filename

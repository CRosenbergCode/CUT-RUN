#!/usr/bin/env bash
# template array script- parallelize over lines in a file 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:01:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --job-name=slurm-job
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log
set -e # quit on error
echo "[$0] $SLURM_JOB_NAME $@" # log the command line

# set temporary directories to the node's tmp space
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR

filename=example-metadata.txt # $1

date # timestamp
module load anaconda
conda activate cutruntools2.1
linenum=0
while read -a fields
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
        cmd="command ${fields[1]} ${fields[0]} -n $SLURM_NTASKS" # example command that uses the requested ntasks
        echo "$cmd"  
        #time eval $cmd
        break
    fi
    linenum=$((linenum + 1))

done < $filename

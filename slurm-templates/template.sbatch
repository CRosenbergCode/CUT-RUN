#!/usr/bin/env bash
# template script with Alpine settings
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --job-name=slurm-job
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%j.log # gives slurm.ID.log
set -e # quit on failure
echo "[$0] $SLURM_JOB_NAME $@" # log the command line

# set temporary directories to the node's tmp space
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR

date # timestamp
module load anaconda
conda activate cutruntools2.1

cmd="command -n $SLURM_NTASKS" # example command that uses the requested ntasks
echo "$cmd" # 
time eval $cmd

#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=20:00:00
#SBATCH --qos=normal
#SBATCH --partition=day-long-cpu
#SBATCH --job-name=BamConvert
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%A.log # gives slurm.ID.log

module purge
source activate base
conda activate rnaPseudo

date # timestamp
nthreads=$(($SLURM_NTASKS - 1))


linenum=0
for file in *.sam
do
  TWO=${file/.sam/.bam}
  echo ${TWO}
  samtools sort -@$nthreads <(samtools view -bS -@$nthreads $file)  > $TWO
done

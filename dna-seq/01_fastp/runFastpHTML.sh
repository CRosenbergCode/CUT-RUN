#!/bin/bash
#SBATCH --partition=day-long-cpu
#SBATCH --job-name=fastpLoop
#SBATCH --output=%x.%j.out
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=17
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
#source activate base
#conda activate rnaPseudo
source $HOME/miniconda3/bin/activate cutrun

mkdir -vp trimmed htmls

for FILE in raw/*R1_001.fastq.gz
do
  echo $FILE
  #TWO=$(echo $FILE | rev | sed s/./2/14 | rev) #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
  TWO=${FILE/R1/R2}
  TRIM1="trimmed/trimmed."$(basename $FILE)
  TRIM2="trimmed/trimmed."$(basename $TWO)
  # fastp uses 16 threads max (version 0.23.4)
  cmd="fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2} -h htmls/$(basename $TWO).html -j htmls/$(basename $TWO).json -w $((SLURM_NTASKS-1)) --dedup"
  echo $cmd
  time eval $cmd
done

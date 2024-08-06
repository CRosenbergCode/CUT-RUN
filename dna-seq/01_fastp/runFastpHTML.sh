#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=fastpLoop
#SBATCH --output=%x.%j.out
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

module purge
source activate base
conda activate rnaPseudo

for FILE in rawFastq/*R1_001.fastq.gz
do
  echo $FILE
  TWO=$(FILE/R1/R2) #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
  TRIM1=$(FILE/rawFastq\//trimmedFastq\/.trimmed)
  TRIM2=$(TWO/rawFastq\//trimmedFastq\/.trimmed)
  #TRIM1="trimmedFastq/trimmed."${FILE}
  #TRIM2="trimmedFastq/trimmed."${TWO}
  echo ${TWO}
  echo ${TRIM1}
  fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2} -h HTMLOutput/${TWO}.html -w 5
done

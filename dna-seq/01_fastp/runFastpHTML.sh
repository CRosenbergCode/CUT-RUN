#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=fastpLoop
#SBATCH --output=%x.%j.out
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

module purge
module load anaconda
conda activate rnaPseudo

for FILE in *R1_001.fastq.gz
do
  echo $FILE
  TWO=$(echo $FILE | rev | sed s/./2/14 | rev) #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
  TRIM1="trimmed/trimmed."${FILE}
  TRIM2="trimmed/trimmed."${TWO}
  #echo ${TWO}
  #echo ${TRIM1}
  fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2} -h ${TWO}.html -w 5
done

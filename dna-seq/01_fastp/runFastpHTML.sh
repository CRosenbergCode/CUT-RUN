#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=fastpLoop
#SBATCH --output=%x.%j.out
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
source activate base
conda activate rnaPseudo

for FILE in raw/*R1_001.fastq.gz
do
  echo $FILE
  #TWO=$(echo $FILE | rev | sed s/./2/14 | rev) #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
  TWO=${FILE/R1/R2}
  TRIM1="trimmed/trimmed."${FILE}
  TRIM2="trimmed/trimmed."${TWO}
  fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2} -h htmls/${TWO}.html -w 5
  #increase -w to increase number of threads, but be sure to increase the number of threas as well
done

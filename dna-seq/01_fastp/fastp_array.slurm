#!/bin/bash
#SBATCH --partition=day-long-cpu
#SBATCH --job-name=fastpArray
#SBATCH --output=%x.%j.out
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=17
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module purge
source $HOME/miniconda3/bin/activate
#conda activate fastp
conda activate rnaPseudo

mkdir -vp trimmed htmls

fqnum=0
while read -r line
do
	if [ $SLURM_ARRAY_TASK_ID -eq $fqnum ]
	then
		echo $line
		#TWO=$(echo $FILE | rev | sed s/./2/14 | rev) #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
		TWO=${line/_R1_/_R2_}
		TRIM1="trimmed/trimmed."$(basename $line)
		TRIM2="trimmed/trimmed."$(basename $TWO)
		# fastp uses 16 threads max (version 0.23.4)
		cmd="fastp -i ${line} -I ${TWO} -o ${TRIM1} -O ${TRIM2} -h htmls/$(basename $TWO).html -j htmls/$(basename $TWO).json -w $((SLURM_NTASKS-1))"
		echo $cmd
		time eval $cmd
	fi
	fqnum=$((fqnum + 1))
done < $1


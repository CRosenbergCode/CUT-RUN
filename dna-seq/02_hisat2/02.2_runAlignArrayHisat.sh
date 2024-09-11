#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=128:00:00
#SBATCH --qos=normal
#SBATCH --partition=week-long-cpu
#SBATCH --job-name=HisatArrayCxt
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log

echo "[$0] $SLURM_JOB_NAME $@" # log the command line

module purge
source activate base
conda activate rnaPseudo

date # timestamp

filename=FastqTrimmed.txt # $1
ls ../fastp_01/trimmed/ -p | grep -v / | grep ".fastq.gz" | grep "R1" | grep "trimmed." >> ${filename}
sed -e ‘s/^\/..\/fastp_01\/trimmed\//’ -i ${filename}

linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
      TWO=${line/_R1_/_R2_}
      OUT=${TWO:16:(${#TWO}-25)} #Remove trimmed, and end. 23 is 8 (trimmed.) + 15 (R1_001$
      echo $OUT
      #STRANDNESS MAY NEED TO CHANGE FOR DIFFERENT LIBRARY PREP KITS
      #Does not matter for DNA :D
      hisat2 --phred33 -p 32 -x hisatIndex/genomeAedesEColi -1 ${line} -2 ${TWO} -S HisatAligned/$OUT.sam
    fi
    linenum=$((linenum + 1))
done < $filename

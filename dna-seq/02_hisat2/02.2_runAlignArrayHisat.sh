#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=128:00:00
#SBATCH --qos=normal
#SBATCH --partition=week-long-cpu
#SBATCH --job-name=HisatArrayCxt
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log
#SBATCH --error=%x.%A-%a.err

echo "[$0] $SLURM_JOB_NAME $@" # log the command line

module purge
source $HOME/miniconda3/bin/activate base
conda activate rnaPseudo

date # timestamp

filename=FastqTrimmed.txt # $1
fastqfolder=../01_fastp/trimmed/
ls ${fastqfolder} -p | grep -v / | grep ".fastq.gz" | grep "R1" | grep "trimmed." > ${filename}
# This isn't necessary; if you ls in a directory and don't specify parts of the filename, you won't get
# the directory prefix. Also it seemed to be suffering from some sort of Word formatting problems with "smart" quotes, so I changed that.
#sed -e 's/^\/..\/01_fastp\/trimmed\///' -i ${filename}

linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
      TWO=${line/_R1_/_R2_}
      OUT=${TWO:8:-16} # Cut off "trimmed." (8 characters) and "_R#_001.fastq.gz" (16 characters)
      echo $OUT
      #STRANDNESS MAY NEED TO CHANGE FOR DIFFERENT LIBRARY PREP KITS
      #Does not matter for DNA :D
      hisat2 --phred33 -p 32 -x hisatIndex/genomeAedesEColi -1 ${fastqfolder}${line} -2 ${fastqfolder}${TWO} -S HisatAligned/$OUT.sam
    fi
    linenum=$((linenum + 1))
done < $filename

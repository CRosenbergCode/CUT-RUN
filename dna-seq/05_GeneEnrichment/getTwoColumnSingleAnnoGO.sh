#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=week-long-cpu
#SBATCH --job-name=TwoColumnGOAnnos
#SBATCH --mail-user=$USER
#SBATCH --mail-type=all
#SBATCH --output=%x.%A-%a.log # gives slurm.ID.log

echo "[$0] $SLURM_JOB_NAME $@" # log the command line
date # timestamp


module purge
source activate base
conda activate rnaPseudo

module load slurm


#--wait

echo "Extracting annotations from $1"

grep "GO:[0-9]\+[,.*GO:[0-9]\+]*" $1 > IntermediateEggNog.tsv

grep -o "GO:[0-9]\+,.*GO:[0-9]\+" IntermediateEggNog.tsv > GoTermsOnlyEggNog.txt

grep -o "AAEL[0-9]\+" IntermediateEggNog.tsv > GeneBaseOnlyEggNog.txt

paste -d "\t" GeneBaseOnlyEggNog.txt GoTermsOnlyEggNog.txt > MiniEggNog.tsv

split -l 1000 MiniEggNog.tsv "GoTermsSegments" --additional-suffix=".tsv" -d

ls GoTermsSegments*.tsv > Segments.txt

sed -i '1i \Gene.ID\tGo.Terms'  GoTermsSegments*.tsv


num_arr=$(wc -l < Segments.txt)

# Create a slurm job using the GoArray.sh script for each segment
sbatch -W --array=0-$((num_arr - 1)) GoArray.sh

#Do not continue script until slurm array jobs have terminated
wait

cat GoTermsSegments*.tsv.sorted > AllGoSorted.tsv

#rm GoTermsSegments*.tsv*
#rm IntermediateEggNog.tsv
#rm GoTermsOnlyEggNog.txt
#rm MiniEggNog.tsv


#Add header to file
#sed -i '1i \Gene.ID\tGo.Term' AllGoSorted.tsv

#You did it!

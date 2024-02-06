#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --job-name=jobName
#SBATCH --output=%x.%j.out
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colorado.edu

module purge
module load anaconda
conda activate cutruntools2.1

<code>

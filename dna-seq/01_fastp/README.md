# fastp
**Purpose -** Filters and trims fastq files to remove adapters, duplicates, and low-quality reads. [GitHub](https://github.com/OpenGene/fastp). 

**Input -** fastq.gz file from sequencing.

**Output -** fastq.gz file after trimming, deduplication, filtering, etc. Will be used as input for Hisat2, the alignment tool.

## Usage
These scripts are designed to be used on Riviera with Slurm.

1. Put all fastq files to be filtered/trimmed in the `01_fastp/raw` folder. They should end with `R1_001.fastq.gz` or `R2_001.fastq.gz`. Ensure that each R1 file has an R2 equivalent by name.
2. Run the `runFastpHTML.sh` script with Slurm like so: `sbatch runFastpHTML.sh`
3. Wait a few seconds, and then make sure fastp didn't fail by running `sacct`, which shows your recent jobs and their state. (completed, running, or failed)
4. Once the job has completed, move on to the `02_hisat2` folder for running Hisat2. Your trimmed fastq files will be in `01_fastp/trimmed`

## Alternative: Going from .fastq.gz to .bam or .bw, `01_fastqToX.sh`
There are three scripts that schedule a pipeline to convert .fastq.gz files to bams, sorted bams/bais, or bigwigs (in this corresponding order): `01_fastqToBam.sh`, `01_fastqToBamBai.sh`, and `01_fastqToBW.sh`

### 01_fastqToBam.sh and 01_fastqToBamBai.sh
1. Ensure that all of your fastq.gz files are in the `01_fastp/raw` folder.
2. Ensure that the `02_hisat/HisatAligned` folder is clear.
3. Run the relevant script like so:
``` bash
bash 01_fastqToBam.sh
```
or
``` bash
bash 01_fastqToBamBai.sh
```
4. Wait about 10 seconds, then check that everything is running correctly with the following command. You should see the fastp job running, as well as other jobs currently waiting on their dependencies to be met:
``` bash
squeue | grep YourRivieraUsername
```

### 01_fastqToBW.sh
1. Ensure that all of your fastq.gz files are in the `01_fastp/raw` folder.
2. Ensure that the `02_hisat/HisatAligned` folder is clear.
3. *(Optional)* Create your comparisons file, if you are running input-adjustment. Comparisons files are created with the following rules:
	- Each line should have a set of comma-separated names of **fastq.gz** files, with the final fastq file being the one to run comparisons against. (AKA, the input)  Which "R" number (R1 or R2) you use is irrelevant 
	- Here is an example line, comparing both `d1_ac_R1_001.fastq.gz` and `d1_me_R1_001.fastq.gz` to `d1_input_R1_001.fastq.gz`. This will generate two output bigwig files, one for `d1_ac` and one for `d1_me`:
		- `d1_ac.bam,d1_me.bam,d1_input.bam`

**There is also an example comparisons file available in the 01_fastp folder.**

4. Run 01_fastqToBW like so:
``` bash
bash 01_fastqToBW.sh
```
or, if using input-adjustment,
``` bash
bash 01_fastqToBW.sh YourComparisons.txt
```
5. Wait about 10 seconds, then check that everything is running correctly with the following command. You should see the fastp job running, as well as other jobs currently waiting on their dependencies to be met:
``` bash
squeue | grep YourRivieraUsername
```

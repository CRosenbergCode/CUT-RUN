# Hisat2
**Purpose -** Aligns trimmed fastq reads to a reference genome or transcriptome. Generates files to be used as inputs for MACS2 or SEACR.

**Input -** Trimmed fastq.gz files from fastp.

**Output -** SAM file, which will then be converted to BAM format for use with a peak caller like MACS2 or SEACR.

## Usage
These scripts are designed to be used on Riviera with Slurm.

1. If you ran fastp already, your trimmed fastq.gz files should already be in `01_fastp/trimmed`  If you did not run fastp but have trimmed fastq.gz files, make sure the `01_fastp/trimmed` folder is empty then add them to it.
2. If Slurm is not loaded, load it with the command `module load slurm`
3. Download a merged genome file from the NAS to run the 02.1 script with. The current merged genome as of May 16, 2025 is `VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta` and it is located in `NAS_Server_Shared/InformaticsResources/`. Put the genome file in the `02_hisat2` directory.
4. Run `02.1_runBuildAedesae68Genome.sh`. You only need to do this once, and then you will have .ht2 files generated in the `hisatIndices` folder. **You do not need to run it again in subsequent runs unless you change genome files. If you do that, then clear out the hisatIndices folder, then re-run the 02.1 script with the new genome.**
``` bash
sbatch 02.1_runBuildAedes68Genome.sh
```
5. Wait a few seconds, and then make sure 02.1 didn't fail by running the `sacct` command.
6. Once 02.1 is finished, run the `02.2_runAlignArrayHisat.sh` script.
``` bash
sbatch 02.2_runAlignArrayHisat.sh
```
7. Wait a few seconds, and then make sure 02.2 didn't fail by running the `sacct` command.
8. Change directories into the `HisatAligned` folder and run the `02.3_runBamLoop.sh` script. *(Note: this is currently untested but it should work. I'm a bit unsure about the 02.3 script, I think it may only work for the current directory)*
``` bash
sbatch ../02.3_runBamLoop.sh
```
9. Check that this script didn't fail either with the `sacct` command. Once it has finished, move on to `03_MACS2` or `03_SEACR`

## Alternative: `hisat_all.sh` 
1. If you ran fastp already, your trimmed fastq.gz files should already be in `01_fastp/trimmed`  If you did not run fastp but have trimmed fastq.gz files, make sure the `01_fastp/trimmed` folder is empty then add them to it.
2. Download a merged genome (fasta) file from the NAS to run the 02.1 script with. Put the genome file in the `02_hisat2` directory. You only need to do this step once.
	- The current merged genome as of May 16, 2025 is `VectorBase-68_AaegyptiLVP_AGWG_Genome.fasta` and it is located in `NAS_Server_Shared/InformaticsResources/`.
	- **If the genome file has changed, you will need to change the name of the fasta file in the** `02.1_runBuildAedesae68Genome.sh` **script.**
	- **If you want to use different hisat indices, move/delete the old index files in the hisatIndex folder.**
3. Run the `hisat_all.sh` script with the following command:
``` bash
bash hisat_all.sh
```
4. Wait about 10 seconds, then check that everything is running correctly with the following command. You should see the buildEColi job running, as well as other jobs currently waiting on their dependencies to be met:
``` bash
squeue | grep YourRivieraUsername
```

# MACS2
**Purpose -** Identifies peaks in alignment maps by comparing to a background. Can also identify transcription factor binding site locations. Outputs are used as inputs to Diffbind.
**Input -** BAM files from Hisat2 outputs.
**Output -** Bed and .xls files that describe the alignment map peaks, used as inputs for Diffbind.

## Usage
These scripts are designed to be used on Riviera with Slurm. There is currently no pipeline script for MACS2; follow these steps to create your own.

 1. Open the `runMACS2PairedOnly.sh` script. If you are using MobaRTE, make sure the script is using UNIX line endings, as MobaRTE sometimes defaults to DOS. (Click the penguin icon in the bar above the text editing area)
 2. Replace the `macs2 callpeak` commands (starting on line 24) with your own, using the following steps. If you need to check the documentation for the callpeak function, check [the MACS3 docs](https://macs3-project.github.io/MACS/docs/callpeak.html).
	1. For better legibility, spread the command across multiple lines using backslashes, as is shown in the original script.
	2. The `-t` argument should be the path to the bam file from your Hisat2 run.
	3. The `-f` argument should be `BAMPE`, since we are using paired-end reads and .bam files as inputs.
	4. The `-g` argument is the genome length. Set it to `1.3e+9` for Aedes aegypti.
	5. Use the `--broad` argument if your histone modification produces broad peaks, like H3K9me3. Do not use it for narrow peaks, like H3K27ac.
	6. **(NOTE BEFORE PUSHING: The `-c` argument is supposed to be the "control," which I can only assume is the "input" .bam file for the same run. But, I'm not sure about this, so I want to double-check first.)**
	7. The `-n` argument is the name and will be used as names for the output .xls and .bed files. Set it to something that makes sense (like the group, histone modification, and replicate. `BF_Me_2` for example) and avoid duplicates.
	8. The `--slocal` argument should be `2000` for Aedes aegypti.
	9. The `--outdir` argument is the directory that the MACS2 outputs will go to. This should be in the `MACS2Peaks` folder.
	10. Example of all of this put together for the second replicate of methylation in the bloodfed group:
		`macs2 callpeak -t ../02_hisat2/HisatAligned/BF_Me_Rep2.bam \`
		`-f BAMPE -g 1.3e+9 \`
		`-c ../02_hisat2/HisatAligned/BF_Me_Rep2_Control.bam`
		`-n BF_Me_2`
		`--slocal 2000`
		`--outdir MACS2Peaks/BF_Me_2`
3. Repeat the substeps of step 2 for each pair of bam files in the Hisat2 outputs.
4. When you are finished, use slurm to run the new `runMACS2PairedOnly.sh` script using `sbatch runMACS2PairedOnly.sh` and check that it's working by running `sacct`
5. After the script has finished, you can move on to the Diffbind step.

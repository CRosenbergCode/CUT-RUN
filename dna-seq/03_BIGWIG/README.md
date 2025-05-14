# Bigwig (DeepTools)
**Purpose -** Creates bigwig (.bw) files for viewing in IGV. For our purposes, these files will show up as "tracks" in IGV that show changes in immunoprecipitation at different regions in the genome.
**Input -** .bam files from HISAT2.
**Output -** .bw files describing the changes in immunoprecipitation.

## Usage
These scripts are designed to be used on Riviera with Slurm. Note that this step also requires the deeptools_kernel environment, but the scripts should automatically handle the changing of environments for you.

1. Move the .bam files from the HISAT folder into the current (03_BIGWIG) directory. This can be achieved with the following command (you must be in the 03_BIGWIG directory for this to work):
``` bash
mv ../02_hisat2/HisatAligned/*.bam .
```
2. Run the following command:
``` bash
sbatch 03.1_runBamCoverage.sh
```
3. While the slurm job is running (you can check its progress with the `sacct` command), you can rewrite the `03.2_runBWCompare.sh` script to work with your files. Follow these steps to do so:
	1. Remove the `bigwigCompare` lines in the script. You will be re-making these for your own run.
	2. Find the two bigwigs (.bw files) you want to compare. The current runBamCoverage script will create .bw files with the exact same name as your .bam files, but inside a directory called `bin10`. One bigwig should be a treatment (acetylation, methylation, etc.) and the other should be the control, or "input" that does not use immunoprecipitation.
	3. Those two bigwigs will be your `-b1` and `-b2` arguments. Start your first `bigwigCompare` call by writing the following, replacing the placeholder file names with your own bigwigs:
		`bigwigCompare -b1 FIRST_BIGWIG.bw -b2 SECOND_BIGWIG.bw`
	4. Select an output file name (should be named after the treatment) and specify it as the -o directory. The most recent version of the script creates a directory called `bw_compare`, so you may want to put it in that. Your command should now be:
		`bigwigCompare -b1 FIRST_BIGWIG.bw -b2 SECOND_BIGWIG.bw -o bw_compare/OUTPUT_BIGWIG.bw`
	5. Finally, specify a bin size and processor count. We currently use 10 and 4 respectively. Your final command should look like this:
		`bigwigCompare -b1 FIRST_BIGWIG.bw -b2 SECOND_BIGWIG.bw -o ./bw_compare/OUTPUT_BIGWIG.bw --binSize 10 --numberOfProcessors 4`
	6. Repeat these steps for each treatment/input pair.
4. Once the `03.1_runBamCoverage.sh` slurm job has finished, run your new `03.2_runBWCompare.sh` script. When it finishes, your outputs will be in the `bw_compare` folder.
``` bash
sbatch 03.2_runBWCompare.sh
```

## Alternative: Faster Array-Based Approach, `SICC_array.sh`
1. *(Optional)* If you will be running comparisons, create your comparisons file in the directory your bams are in:
	- Each line should have a set of comma-separated paths to bam files, with the final bam file being the one to run comparisons against. The paths should be relative to the comparisons file, so you should only need the names of the bams, not the folders they are in.
	- Here is an example line, comparing both `d1_ac.bam` and `d1_me.bam` to `d1_input.bam`. This will generate two output bigwig files, one for `d1_ac` and one for `d1_me`:
		- `d1_ac.bam,d1_me.bam,d1_input.bam`

2. Run SICC_array.sh with the following slurm command: `bash SICC_array.sh path/to/bams/folder (0 or 1) path/to/bams/folder/comparisons.txt` (assuming your comparisons file is in your bams folder)
	- The first argument is the path to the folder containing all the bams you would like to convert to bigwigs.
	- The second argument is whether the bams should be sorted and indexed (If they are already sorted, and you already have the `.bai` index files in the folder, you can enter `0`. Otherwise, enter `1`. *Do not include parentheses*)
	- The third argument is optional. If you have a comparisons file in your bams folder, using the instructions written above in step 1, 

3. Wait about 10 seconds, then check that all the jobs were properly scheduled using the following command:
``` bash
squeue | grep YourRivieraUsername
```

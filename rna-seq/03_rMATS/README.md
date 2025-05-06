# DAS with rMATS-turbo
**Purpose -** Identifies differential alternative splicing (DAS) events between two different RNA-seq samples. (such as bloodfed vs RVFV)

**Input -** .fastq.gz files from fastp outputs.

**Output -** .txt files containing information on the different alternative splicing events, as well as summary and metadata information. These can be used with subsequent scripts to make .csv files with DAS information.

## Usage
1. Move your trimmed .fastq files from the fastq directory to a folder of your choice in the rMATS directory. (optional, but makes writing the comparisons files easier) You can do this by changing into the rMATS directory and running the following command: `mkdir fastq && mv ../01_fastp/trimmed/*.fastq.gz ./fastq/`
2. Run STAR, the alignment tool:
	1. First, you must generate the indices. Run the `generate_STAR_indices.sh` script like so, giving it the path to your genome fasta and the organism GFF/GTF:
		- `sbatch generate_STAR_indices.sh path/to/your/genome.fasta path/to/your/gfffile.gff`
	2. Then, run STAR, providing it with the path to the directory with your .fastq.gz files:
		- `sbatch run_STAR.sh path/to/fastqs_directory` (if you used the command in step 1, the path to your fastqs directory will just be `fastq`
	3. Ensure it ran properly using `sacct`. Your bam files will be in the same directory as your fastq files.
3. Create your "b" files. Write the paths/names of the bams for each replicate in a group, separated by commas, into a text file.
	- For example, if we had samples `BF_d1_r1.bam`, `BF_d1_r2.bam`, and `BF_d1_r3.bam`, all in the `STAR_BF_out` folder, then we would write `STAR_BF_out/BF_d1_r1.bam,STAR_BF_out/BF_d1_r2.bam,STAR_BF_out/BF_d1_r3.bam` into a text file, which we might name `BF_d1.txt` or something similar.
4. For each pair of "b" files that you want to compare, run rMATS. This is done with the `rmats_one.sh` slurm script:
	- Run it like so: `sbatch rmats_one.sh X Y Z` where X is your first "b" file, Y is your second "b" file, and Z is the name of the run. The name of the run will be used as directory prefixes, so it is advisable to follow directory naming conventions, such as avoiding spaces.
5. Alternatively, if you have multiple comparisons to run, you can use the `rmats_batch.sh` script. Here's how to name your "b" files for rmats to auto-detect and name the runs:
	- The first part of the "b" file name should be the run name.
	- The second part should either be `_b1.txt` or `_b2.txt` depending on which "b" file it is in the pair.
	- For example, if we had a sugarfed vs bloodfed run we wanted to name "SFvsBF_d1_2025", we would put our sugarfed bam file names in `SFvsBF_d1_2025_b1.txt` and our bloodfed bams in `SFvsBF_d1_2025_b2.txt`
	- When you're done naming your bams, run `bash rmats_batch.sh` from the same directory as your "b" files.

## Making significant CSVs from rMATS output
Here we go over how to use the `make_deduped_csv.sh` and `make_all_deduped_csvs.sh` scripts with a genome annotations file (CSV or TSV) to create CSV files with all significant rMATS outputs. A genome annotations file is available on the NAS, currently under the file path `NAS_Server_Shared/InformaticsResources/Aedes_aegypti/Aedes_Ae_gene_annotations_TRANSCRIPT-LEVEL_GO2_14_2025.tsv`
1. If the current version of the annotations file is not tab-delimited (a TSV), then you must make a tab-delimited version. This is very easy to do. Open it in Excel or your preferred spreadsheet editor and re-save or export as a text (TSV) file. Be sure that auto-convert dates is off, so your genes do not turn into dates. You can also use an online CSV to TSV converter.
2. To convert one rMATS .txt output into a CSV, run `make_deduped_csv.sh` like so:
	- `bash make_deduped_csv.sh path/to/rmats_output_file.txt path/to/annotations.tsv`
	- The rMATS outputs you should convert all end in `MATS.JC.txt` or `MATS.JCEC.txt`. Files ending in JC only use exon junction reads to determine alternative splicing. JCEC files use both exon junction reads as well as exon body reads.
3. To convert all JC and JCEC .txt files into CSVs, run `make_all_deduped_csvs.sh` like so:
	- `bash make_all_deduped_csvs.sh path/to/rmats_out path/to/annotations.tsv`, where `path/to/rmats_out` is simply the output directory of your rMATS run.
4. Your output CSVs will be in the directory you ran the scripts in.

#!/bin/bash

#SBATCH --partition=short-cpu
#SBATCH --job-name=htmlSummaryLoop
#SBATCH --output=%x.%j.out
#SBATCH --time=0:30:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
##SBATCH --mail-user=haogg@colostate.edu

# Last modified 2/3/25 by Zoey Mikol, from Hunter Ogg's original script
# Usage: Put it in the same directory as your html and log files and run
# 	with command: ./summarizeChipMetadata.sh

# Function to convert a string like "1.25 G" to 1250000000 or "37.548 K" to 37548.
# Still works if there are other terms after the number. For example, an input
# of "12.47 M (97.8434%)" just yields the output "12470000", no percentage.
convert_units() {
    echo "$1" | awk '{
        if ($2 == "G") {
            $1 = $1 * 1000000000;
        } else if ($2 == "M") {
            $1 = $1 * 1000000;
        } else if ($2 == "K") {
            $1 = $1 * 1000;
        }
	print int(strtonum($1));
    }'
}

# Function to pull the percentage out of a line and convert it to a proportion.
convert_pct() {
	# The alignment rate is sometimes unspecified and written as a hash sign,
	# so we don't want to change it in that case.
	if test "$1" = "#" ; then
		printf "#"
		exit
	fi

	# Grab everything after the first open parenthesis...
	pct=${1#*(}
	# ... and before the first percent after that.
	pct=${pct%\%*}
	# Multiply the number by 0.01. 8 sig figs.
	awk -vx=$pct 'BEGIN{printf "%.8f" ,x * 0.01}'
}

# Create the output tsv and add the column headers.
outfile="ChipMetadata.tsv"
echo -e "Sample Name\tDuplication Rate\tReads Passed\tReads Passed (Proportion)\tTotal\tQ20 Bases\tQ20 Bases (Proportion)\tQ30 Bases\tQ30 Bases (Proportion)\tInsert Size Peak\tMean Insert Size (after filtering)\tAlignment Rate" > ${outfile}

# Not sure what this comment is for but I'm leaving it in just in case it's relevant.
#echo -e "Peak Sample\tTotal Peaks\tPercentage Within 1000\tReads Passed\tTotal\tQ20 Bases\tQ30 Bases" > ${outfile}

for FILE in *.html
do
	base="${FILE/_R1_001.fastq.gz.html/}"
	base2="${base/_R2_001.fastq.gz.html/}"
	echo ${base2}

	# Find the line with the duplication rate and extract the number from it.
	# Run it through the percentage -> proportion function.
	dup=$(grep "duplication rate:" $FILE | grep -o "[0-9]\+\.[0-9]\+")
	dup=$(convert_pct "$dup")

	# Extract the value from the reads passed line, get the proportion and number of reads.
	reads_pass=$(grep "reads passed filters:" $FILE | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")
	reads_pass_pct=$(convert_pct "$reads_pass")
	reads_pass=$(convert_units "$reads_pass")

	# Extract the value of the total reads line, convert it to a number.
	tot_reads=$(grep -A 2 "id='before_filtering_summary'" $FILE | grep "total reads:" | grep -o "[0-9]\+\.[0-9]\+\s[GMK]")
	tot_reads=$(convert_units "$tot_reads")
	# Another old comment I will be leaving here because I'm not sure what it means.
	#Need to determine before or after filtering

	# Get the number/proportion of Q20's and Q30's
	q_20=$(grep "Q20 bases" $FILE -m 1 | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")
	q_20_pct=$(convert_pct "$q_20")
	q_20=$(convert_units "$q_20")

	q_30=$(grep "Q30 bases" $FILE -m 1 | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")
	q_30_pct=$(convert_pct "$q_30")
	q_30=$(convert_units "$q_30")

	# Get the insert size peak and mean insert size after filtering. (These
	# don't need to be run through the conversion functions)
	insert_peak=$(grep "Insert size peak:" $FILE -m 1 | grep -o ">[0-9]\+<" | grep -o "[0-9]\+")
	mean_insert=$(grep "mean length after filtering" $FILE -m 1 | grep -o "[0-9]\+bp," | grep -o "[0-9]\+")

	# Get the overall alignment rate...
	align_rate=""#$(grep "overall alignment rate" $FILE | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")

	# Loop through the hisat log files in the current directory until a match is found.
	for LOG in *.log
	do
		if grep -q $base2 "$LOG"; then
			# Get the error file, and pull the overall alignment rate from that.
			ERR=${LOG/.log/.err}
			align_rate=$(grep "overall alignment rate" $ERR | grep -o "[0-9]\+\.[0-9]\+%")
		fi
	done
	echo "${base2}, $align_rate"
	# Finally, convert it to a proportion.
	align_rate=$(convert_pct "$align_rate")

	# Write all the pulled values to the console and the output tsv.
	write_line="${base2}\t${dup}\t${reads_pass}\t${reads_pass_pct}\t${tot_reads}\t${q_20}\t${q_20_pct}\t${q_30}\t${q_30_pct}\t${insert_peak}\t${mean_insert}\t${align_rate}"

	echo -e ${write_line}

	echo -e ${write_line} >> ${outfile}
done

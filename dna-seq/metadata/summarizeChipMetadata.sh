#!/bin/bash

#SBATCH --partition=day-long-cpu
#SBATCH --job-name=htmlSummaryLoop
#SBATCH --output=%x.%j.out
#SBATCH --time=0:30:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haogg@colostate.edu

outfile="ChipMetadata.tsv"
#echo "Duplication Rate  Reads Passed  Total Reads Q20 Bases Q30 Bases" > ${outfile}
echo -e "Sample Name\tDuplication Rate\tReads Passed\tTotal\tQ20 Bases\tQ30 Bases\tInsert Size Peak\tMean Insert Size (after filtering)\tAlignment Rate" > ${outfile}

#echo -e "Peak Sample\tTotal Peaks\tPercentage Within 1000\tReads Passed\tTotal\tQ20 Bases\tQ30 Bases" > ${outfile}

write_line="${FILE}\t${dup}\t${reads_pass}\t${tot_reads}\t${q_20}\t${q_30}"

for FILE in *.html
do
	base="${FILE/_R1_001.fastq.gz.html/}"
	base2="${base/_R2_001.fastq.gz.html/}"
	#base="${base/R1/R2}"
	echo ${base2}
	dup=$(grep "duplication rate:" $FILE | grep -o "[0-9]\+\.[0-9]\+%")

	#Don't want the percentage here....

	reads_pass=$(grep "reads passed filters:" $FILE | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")
	#Split above into 2? Just get percentage?

	tot_reads=$(grep -A 2 "id='before_filtering_summary'" $FILE | grep "total reads:" | grep -o "[0-9]\+\.[0-9]\+\s[GMK]")
	#tot_reads=$(grep "total reads:" $FILE | grep -o "[0-9]\+\.[0-9]\+\s[GMK]")
	#Need to determine before or after filtering

	q_20=$(grep "Q20 bases" $FILE -m 1 | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")

	q_30=$(grep "Q30 bases" $FILE -m 1 | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")

	#write_line="${FILE} ${dup} ${reads_pass}  ${tot_reads}  ${q_20} ${q_30}"

	insert_peak=$(grep "Insert size peak:" $FILE -m 1 | grep -o ">[0-9]\+<" | grep -o "[0-9]\+")
	#echo $insert_peak
	mean_insert=$(grep "mean length after filtering" $FILE -m 1 | grep -o "[0-9]\+bp," | grep -o "[0-9]\+")

	#mean length after filtering

	align_rate=""#$(grep "overall alignment rate" $FILE | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")

	#extension="${FILE##*.}"
	for LOG in *.log
	do
		if grep -q $base2 "$LOG"; then
			#echo "Found File!"
			ERR=${LOG/.log/.err}
			align_rate=$(grep "overall alignment rate" $ERR | grep -o "[0-9]\+\.[0-9]\+%")
		fi
		#align_rate=$(grep "overall alignment rate" $FILE | grep -o "[0-9]\+\.[0-9]\+\s[GMK].*%)")
	done

	write_line="${base2}\t${dup}\t${reads_pass}\t${tot_reads}\t${q_20}\t${q_30}\t${insert_peak}\t${mean_insert}\t${align_rate}"

	echo -e ${write_line}

	echo -e ${write_line} >> ${outfile}
done

#if grep -q SomeString "$File"; then
#  Some Actions # SomeString was found
#fi

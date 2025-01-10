first_line="\tTotal Peaks\tAverage Peak Width"
for PROM in *.gtf
do
	first_line+="\t"
	first_line+=$PROM
done

write_out=$first_line

for PEAK in *.bed
do
	av_len=$(awk '{ total += ($3-$2) } END { print total/NR }' $PEAK)
	tot=$(wc -l $PEAK | cut -d " " -f 1)

	next_line=$PEAK
	next_line+="\t"
	next_line+=$tot
	next_line+="\t"
	next_line+=$av_len

	for PROM in *.gtf
	do
		next_line+="\t"
		within=$(bedtools intersect -a $PROM -b $PEAK | wc -l | cut -d " " -f 1)
		perc=$(echo "scale=4 ; $within / $tot" | bc)

		#ad_perc=$(echo "scale=8 ; $perc / $av_len" | bc)
		#mult_perc=$(echo "scale=4 ; $ad_perc * 10000" | bc)

		perc=$(echo "scale=4 ; $within / $tot * 10000 / $av_len" | bc)
		next_line+=$perc
	done
	write_out+="\n"
	write_out+=$next_line
done


for PEAK in *.narrowPeak
do
	av_len=$(awk '{ total += ($3-$2) } END { print total/NR }' $PEAK)
	tot=$(wc -l $PEAK | cut -d " " -f 1)

	next_line=$PEAK
	next_line+="\t"
	next_line+=$tot
	next_line+="\t"
	next_line+=$av_len

	for PROM in *.gtf
	do
		next_line+="\t"
		within=$(bedtools intersect -a $PROM -b $PEAK | wc -l | cut -d " " -f 1)
		perc=$(echo "scale=4 ; $within / $tot" | bc)

		#ad_perc=$(echo "scale=8 ; $perc / $av_len" | bc)
		#mult_perc=$(echo "scale=4 ; $ad_perc * 10000" | bc)

		perc=$(echo "scale=4 ; $within / $tot * 10000 / $av_len" | bc)
		next_line+=$perc
	done
	write_out+="\n"
	write_out+=$next_line
done
echo -e $write_out > "PromoterSummaryLengthAdjusted.tsv"


#echo "scale=4 ; 0.0123 / 9145.61" | bc

#s9145.61

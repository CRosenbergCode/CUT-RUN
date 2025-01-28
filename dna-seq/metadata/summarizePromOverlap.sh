outfile="ParamSummaryLengthAdjusted.tsv"

first_line="\tTotalPeaks\tMeanPeakWidth\tMedianPeakWidth\tQValue\tD\tSlocal\tLlocal\tLambda"
for PROM in *.gtf
do
	first_line+="\t"
	first_line+=$PROM
done
#first_line+="\n"
echo -e $first_line > $outfile
#write_out=$first_line

for PEAK in *.narrowPeak
do
	av_len=$(awk '{ total += ($3-$2) } END { print total/NR }' $PEAK)
	tot=$(wc -l $PEAK | cut -d " " -f 1)
  med_len=$(awk '{ print $3-$2; }' $PEAK | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
  #echo $PEAK
	q_val=$(echo $PEAK | grep -o "q[0-9]\+\.[0-9]\+" | grep -o "[0-9]\+\.[0-9]\+")
	l_val=$(echo $PEAK | grep -o "l[0-9]\+" | grep -o "[0-9]\+")
	s_val=$(echo $PEAK | grep -o "s[0-9]\+" | grep -o "[0-9]\+")
	d_val=$(echo $PEAK | grep -o "d[0-9]\+" | grep -o "[0-9]\+")

  lamd=0
	if [[ $PEAK  == *"lamba"* ]]; then
	   lamd=1
	fi

	next_line=$PEAK
	next_line+="\t"
	next_line+=$tot
	next_line+="\t"
	next_line+=$av_len
  next_line+="\t"
	next_line+=$med_len

	next_line+="\t"
	next_line+=$q_val
	next_line+="\t"
	next_line+=$d_val
  next_line+="\t"
	next_line+=$s_val
	next_line+="\t"
	next_line+=$l_val
	next_line+="\t"
	next_line+=$lamd



	for PROM in *.gtf
	do
		next_line+="\t"
		within=$(bedtools intersect -b $PROM -a $PEAK | wc -l | cut -d " " -f 1) #b then a, otherwise peaks overlapping multiple promoters can be incorrectly called
		perc=$(echo "scale=4 ; $within / $tot" | bc)

		perc=$(echo "scale=4 ; $within / $tot * 10000 / $av_len" | bc)
		next_line+=$perc
	done
	#next_line+="\n"
	#write_out+="\n"
	#write_out+=$next_line
	echo -e $next_line >> $outfile
done
#echo -e $write_out > "ParamSummaryLengthAdjusted.tsv"

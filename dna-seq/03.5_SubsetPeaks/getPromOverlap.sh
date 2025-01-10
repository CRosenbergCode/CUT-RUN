for PEAK in *.relaxed.bed
do
	next_line+="\t"
	bedtools intersect -wa -b AedesGenesProm10000.gtf -a $PEAK > test.prom
  cp test.prom $PEAK
  rm test.prom
done

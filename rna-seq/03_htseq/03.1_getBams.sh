ls ../02_hisat2/ConcordantOutput/ -p | grep -v / | grep ".bam" > ConcordantBams.txt
sed -e "s/^/..\/02_hisat2\/ConcordantOutput\//" -i ConcordantBams.txt



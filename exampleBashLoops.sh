for FILE in *R1_001.fastq.gz
do
  TWO=$(echo $FILE | rev | sed s/./2/14 | rev) #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
  TRIM1="trimmed/trimmed."${FILE}
  TRIM2="trimmed/trimmed."${TWO}
  #echo ${TWO}
  #echo ${TRIM1}
  fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2}
done


#b=${a:12:5}

for FILE in trimmed.*R1_001.fastq.gz
do
  TWO=$(echo $FILE | rev | sed s/./2/14 | rev) #Read in R1 file, reverse it, replace the 12th character with a 2, and reverse for R2
  OUT=${TWO:8:$(${#TWO}-23)} #Remove trimmed, and end. 23 is 8 (trimmed.) + 15 (R1_001.fastq.gz)
  #echo ${TWO}
  hisat2 --phred33 -x genome -1 ${FILE} -2 ${TWO} -S ${OUT}.sam
done

#Example file name for this sequencing run RNA-BF_Rep1_S21_L006_R1_001.fastq.gz


for FILE in *R1_001.fastq.gz
do
  TWO=${FILE/_R1_/_R2_} # make read pair name
  TRIM1="trimmed/trimmed.${FILE}"
  TRIM2="trimmed/trimmed.${TWO}"
  #echo ${TWO}
  #echo ${TRIM1}
  fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2}
done


#b=${a:12:5}

for FILE in trimmed.*R1_001.fastq.gz
do
  TWO=${FILE/_R1_/_R2_} # make read pair name
  OUT=${TWO/trimmed./} # remove trimmed. prefix
  OUT=${OUT/.fastq.gz/} # remove .fastq.gz
  #echo ${TWO}
  hisat2 --phred33 -x $1 -1 ${FILE} -2 ${TWO} -S ${OUT}.sam
done

#Example file name for this sequencing run RNA-BF_Rep1_S21_L006_R1_001.fastq.gz


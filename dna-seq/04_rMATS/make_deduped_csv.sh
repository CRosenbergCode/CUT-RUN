# ========================== make_deduped_csv.sh ==========================
# 
#  DESCRIPTION: Makes a CSV of all significant alternative splicing events
#          in an rMATS output file, duplicates excluded. Echoes CSV name to console.
# 
#  USAGE: make_deduped_csv.sh $1 $2
#        Argument 1 - Path to the rMATS output file, like A3SS.MATS.JC.txt.
#        Argument 2 - Path to your genome annotations file. (an A. aegypti
#                      genome annotations file is in the NAS folder)
#
#  (Ensure that annotations are tab-delimited. This is an option in Excel)
#
# =========================================================================

annotations=$2

# Get the filename without any file extension or previous URL.
filename=$(printf $1 | rev | cut -d/ -f1 | cut -d. -f2- | rev)

# Get tab-delimited index of p-value.
pfield=`expr $(sed -n '1p' $1 | sed 's/FDR.*$//' | tr -cd '\t' | wc -c) + 1`

# A3SS, A5SS, MXE, RI, or SE, from beginning of file name.
splicetype=$(printf $filename | cut -d. -f1)

# Select the relevant columns to enter into CSV, store in $indices
indices=$(printf "1\n2\n4\n5\n${pfield}")
columnNames=$(printf "ID\nGeneID\nChromosome\nStrand\nFDR")
case $splicetype in

  'A3SS' | 'A5SS')
    indices=$(echo -e $indices"\n6\n7\n8\n9")
    columnNames=$columnNames$(printf "\nLongES\nLongEE\nShortES\nShortEE")
    ;;
  'MXE')
    indices=$(echo -e $indices"\n6\n7\n8\n9")
    columnNames=$columnNames$(printf "\nFirstES\nFirstEE\nSecondES\nSecondEE")
    ;;
  'RI' | 'SE')
    indices=$(echo -e $indices"\n6\n7")
    columnNames=$columnNames$(printf "\nExonStart\nExonEnd")
    ;;
  *)
    echo -e "Error! Unknown rMATS splicing type \"${splicetype}\" detected! Should be A3SS, A5SS, MXE, SE, or RI!"
    exit
    ;;
esac

# Create first line of CSV with $columnNames
touch significant_${filename}.csv
firstLine=$(echo $columnNames | tr ' ' ',')",Annotation"
echo $firstLine > significant_${filename}.csv

# Get indices of all lines with p-value less than 0.05
line_nums=$(awk -v pindex="$pfield" -F'\t' '$pindex < 0.05 {print NR}' $1)

# Used later for comparing duplicates
compareIndex=${indices[@]: -1}

# Add in data, pulling from the fields indicated by $indices
for line in $line_nums; do
  csv_line=""
  match=$(sed -n "${line}p" $1 | cut -f"6-$compareIndex" | tr '\t' ',')
  if grep -q -- "$match" significant_${filename}.csv; then
    continue
  fi
  
  geneid=$(sed -n "${line}p" $1 | cut -f2)
  geneid=${geneid//\"/}
  annot=$(grep $geneid $annotations | head -1 | cut -f7)
  
  if [[ "$annot" == "unspecified product" ]]; then
    annot=$(grep $geneid $annotations | head -1 | cut -f10)
  fi
  
	for index in $indices; do
		csv_line=${csv_line}$(sed -n "${line}p" $1 | cut -f"$index")","
	done
  echo -e "$csv_line$annot" >> significant_${filename}.csv
done

# Echoes file name that it put the data into, so you can pipe it to another script if necessary.
echo "significant_${filename}.csv"

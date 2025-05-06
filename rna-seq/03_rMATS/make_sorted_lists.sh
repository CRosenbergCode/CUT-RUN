#!/bin/bash -l

#SBATCH --job-name=makelist
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition=short-cpu

# ======================== make_sorted_lists.sh =========================
#
#  DESCRIPTION: Makes two sorted CSVs (one by gene ID, the other by
#          description) containing all splicing events detailed in the
#          CSVs from rMATS outputs. (Generated from make_all_deduped_csvs.sh)
#          Can be run with sbatch as well, since it's a little slow.
#
#  USAGE: make_sorted_lists.sh $1 $2 $3
#        Argument 1 - Path to rMATS CSVs directory.
#        Argument 2 - Prefix for sorted CSV files.
#        Argument 3 - Path to tab-delimited annotation file.
#
#  OUTPUT: Sorted CSV files named $2_sort_gene.csv and $2_sort_desc.csv.
#
#  TODO: Include splice sites. This should be easy, but will require a
#      little work in get_unique_pairs to ensure output is clean. (Reverse
#      line, remove first 2 or 4 comma-delimited values depending on the
#      splice type. Very easy)
#
# =======================================================================

# Set up inputs
csvs=$(ls $1/*.JC.csv)
new_prefix=$2
annotfile=$3

# Make CSV files, ensure they're empty.
sortedbygene="${new_prefix}_sort_gene.csv"
sortedbydesc="${new_prefix}_sort_desc.csv"
tmpfile="${new_prefix}_tmp.csv"
printf "" > $tmpfile
printf "" > $sortedbygene
printf "" > $sortedbydesc

# Loop over each CSV in the rMATS directory...
for csv in $csvs; do
  isfirst=1
  # Skipping the first line in the CSV, (because it's the header)
  # Construct a new line with the desired information (see sed commands at bottom)
	while read line; do
		if [[ $isfirst -ne 1 ]] ; then
      gene=$(echo "$line" | cut -d',' -f2 | tr -d '\"')
      annotline=$(grep $gene $annotfile | head -1)
      
      # Get description. Prioritize Product.Description over PFam.Description)
      desc=$(echo "$annotline" | cut -f7)
      if [[ "$desc" == "unspecified product" ]]; then
        desc=$(echo "$annotline" | cut -f10)
      fi
      
      # Construct the next line starting with the GeneID then the gene name...
      nextline="${gene},$(echo "$annotline" | cut -f8)"
      # Then the gene's function... (DIV, CYT/STR, MIT, etc.)
      nextline="${nextline},$(echo "$annotline" | cut -f4)"
      # Then the splice type... (A3SS, A5SS, MXE, etc.)
			nextline="${nextline},$(echo ${csv#*/significant_} | cut -d'.' -f1)"
      # Then the description...
      nextline="${nextline},$desc"
      # Then the gene's start and end locations in the genome.
      nextline="${nextline},$(echo "$annotline" | cut -f6)"
      # Then add the line to a temporary file.
      echo $nextline >> $tmpfile
		else
			isfirst=0
		fi
	done < $csv
done

# Sort the temporary CSV by gene, put into the gene-sorted CSV.
sort -t',' -k1,4 $tmpfile > $sortedbygene
# Do the same but sorting by function and description.
sort -t',' -k3,3 -k5,5 -k1,1 -k4,4 $tmpfile > $sortedbydesc

# Insert header lines at the top.
sed -i '1iGeneID,GeneName,Function,SpliceType,Description,GenomicLocation' $sortedbygene
sed -i '1iGeneID,GeneName,Function,SpliceType,Description,GenomicLocation' $sortedbydesc

# Delete the temporary CSV.
rm $tmpfile
# ========================== count_types.sh ==========================
#
#  DESCRIPTION: Counts the hits of each splicing type (Alt. 3' site,
#          alt 5' site, mutually exclusive, retained intron, skipped exon)
#          in a CSV file in the summary directory.
#
#  USAGE: count_types.sh $1
#        Argument 1 - path to the CSV file to count hits in.
#
#  OUTPUT: Console message detailing hits for each type.
#
# ====================================================================

# List of each of the 5 splice types.
types="A3SS A5SS MXE RI SE"

# Header, shows the file you're counting hits in.
echo -e "\nType counts for ${1}:"

# Loops over each type, counts hits, adds hits to running total, prints hits.
sum=0
for t in $types; do
  count=$(grep "^[^,]*,[^,]*,[^,]*,${t}," $1 -c)
  sum=$(($sum + $count))
  echo -e "${t}: $count"
done

# Prints running total.
echo -e "Total: ${sum}\n"

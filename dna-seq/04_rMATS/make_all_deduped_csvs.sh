# ======================= make_all_deduped_csvs.sh =======================
# 
#  DESCRIPTION: Runs make_deduped_csv.sh for every relevant file in the
#          specified directory. (first argument)
#
#  USAGE: make_all_deduped_csvs.sh $1 $2
#        Argument 1 - Directory with rMATS outputs.
#        Argument 2 - Path to your genome annotations file. (Must be tab-
#                      delimited. an A. aegypti genome annotations file
#                      is in the NAS folder
#
#  make_deduped_csv.sh must be in the current working directory and NOT renamed
#
# ========================================================================

rmatsfiles=$(ls $1/*.txt | grep -v 'summary.txt$' | grep -v 'fromGTF' | grep -v 'raw')
annotationfile=$2

#echo $rmatsfiles

for file in $rmatsfiles; do
  ./make_deduped_csv.sh $file $annotationfile
done
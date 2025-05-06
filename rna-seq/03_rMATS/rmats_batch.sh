#===============================================================
#
#  rmats_batch.sh
#
#  Description: runs rMATS for all "b" files in the CWD.
#
#  Usage: bash rmats_batch.sh
#      (no arguments needed)
#
#  Notes:
#      You must have files ending in "_b1.txt" and "_b2.txt"
#      in your current working directory in order to run this
#      script properly. This is another script I made a while
#      ago and thus needs to be rewritten with better practices.
#
#===============================================================

bfiles=$(ls ./*_b1.txt)

for bfile in $bfiles; do
  run_name=$(echo $bfile | sed 's/.*\///g' | sed 's/_b[1-2]\.txt//g')
  b1=$bfile
  b2=${bfile//b1/b2}
  
  #echo "Name: ${run_name}, B1: ${b1}, B2: ${b2}"
  
  sbatch ./rmats_one.sh $b1 $b2 ${run_name}
done

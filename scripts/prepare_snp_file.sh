#!/bin/bash

# Halt on error.
set -euo pipefail

# Check for the correct number of command-line arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 ped_file map_file output_file [duplicated]"
  exit 1
fi

ped="$1"
map="$2"
output_file="$3"
# duplicated="$4"
# 
# # If duplicated is provided
# if [ -z "duplicated" ]; then
#   # Remove duplicate individuals
#   awk '!seen[$2]++' $ped > ped.temp
#   mv ped.temp $ped
# fi

# Convert genotypes to 0125
awk '{print $2}' $map > outfile.temp
sed 's/$/ 2/' outfile.temp > recode_alleles.txt

plink --ped $ped --map $map --cow --recode A --recode-allele recode_alleles.txt --out plink_0125

rm outfile.temp

# Format snp file to BLUPF90
awk '(NR>1)' plink_0125.raw > plink.temp
awk '{$1=$3=$4=$5=$6="";gsub(FS "+",FS)}1' plink.temp > plink.temp2
sed 's/^[ \t]*//' plink.temp2 > plink.temp3
sed 's/ /,/' plink.temp3 > plink.temp4 
sed 's/ //g' plink.temp4 > plink.temp5
column -s',' -t plink.temp5 > $output_file

rm *temp*

printf "Routine complete."

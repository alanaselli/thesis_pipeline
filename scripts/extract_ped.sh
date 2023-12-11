#!/bin/bash

# Check for the correct number of command-line arguments
if [ "$#" -lt 3 ]; then
  echo "Usage: $0 input_file output_file pedigree_file [num_generations]"
  exit 1
fi

input_file="$1"
output_file="$2"
pedigree_file="$3"
filter_value="$4"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file not found: $input_file"
  exit 1
fi

# If filter_value is provided
if [ -z "$filter_value" ]; then
  cp $input_file out_temp
else
  #awk -v col="8" -v value="$filter_value" '{if ($col > value) print}' "$input_file" > "out_temp"
  awk -v col="8" '{if ($col >= (max - value + 1)) print}' max=$(awk -v col="8" 'NR > 1 {print $col}' "$input_file" | sort -n | tail -n 1) value="$filter_value" "$input_file"> "out_temp"

fi

# Delete header and save to output_file
sed '1d' out_temp > $output_file

# Use awk to extract the first 3 columns and remove the header, then save to pedigree file
awk 'NR > 1 {print $1, $2, $3}' out_temp > "$pedigree_file"

rm out_temp

echo "Data extracted and saved to $output_file and $pedigree_file"


#!/bin/bash

# Check for the correct number of command-line arguments
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 input_file output_file"
  exit 1
fi

input_file="$1"
output_file="$2"

# Use awk create first 6 fields
awk '{split($1, arr, "_"); print "1", $1, arr[1], arr[2], 0, 0}' "$input_file" > temp.file

# Prepare genotypes
awk '!($1="")' "$input_file" > temp.file7
sed -i '' -e 's/1/12/g' temp.file7
sed -i '' -e 's/0/11/g' temp.file7
sed -i '' -e 's/3/21/g' temp.file7
sed -i '' -e 's/4/22/g' temp.file7
sed -i '' -e 's/./& /g' temp.file7

# Merge files
paste -d' ' temp.file temp.file7 > temp.file.final
cat temp.file.final | tr -s ' ' > "$output_file"

rm *temp.file*

echo "$output_file was created sucessfully."
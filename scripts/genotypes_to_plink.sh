#!/bin/bash

# Check for the correct number of command-line arguments
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 input_file output_file"
  exit 1
fi

input_file="$1"
output_file="$2"

# Use awk to process each line
#awk '{split($1, arr, "_"); print "1", $1, arr[1], arr[2], 0, 0, $2}' "$input_file" > temp.file
awk '{split($1, arr, "_"); print "1", $1, arr[1], arr[2], 0, 0}' "$input_file" > temp.file

# Insert spaces in the last column
#awk '{print $7}' temp.file > temp.file7
awk '{print $2}' "$input_file" > temp.file7
sed 's/./& /g' temp.file7 > temp.file_spaces
paste -d' ' temp.file temp.file_spaces > "$output_file"

rm *temp.file*

echo "$output_file was created sucessfully."
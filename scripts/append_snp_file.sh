#!/bin/bash

# Halt on error.
set -euo pipefail

# Check for the correct number of command-line arguments
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 original_file file_to_append"
  exit 1
fi

original_file="$1"
file_to_append="$2"

cat $file_to_append >> $original_file
column -s' ' -t $original_file > temp.file
mv temp.file $original_file

echo "$file_to_append was successfully appended to $original_file."
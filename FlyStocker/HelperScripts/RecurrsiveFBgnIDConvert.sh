#!/usr/bin/env bash

if [ -z "$1" ]; then
  echo "Usage: $0 <BASE_DIR>"
  exit 1
fi

BASE_DIR="$1"

# Loop through each subdirectory in BASE_DIR
for folder in "$BASE_DIR"/*; do
  # If it's a directory (and not just a file)
  if [ -d "$folder" ]; then
    echo "Running GetFBgnIDs.py on: $folder"
    python3 GetFBgnIDs.py "$folder" ext_gene
  fi
done

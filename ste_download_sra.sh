# Script to downlad SRA files with the split them if they are Paired end.
# Input is one or more txt-file(s) with SRA numbers, in the ROOT_DIR folder.

#!/usr/bin/env bash
set -euo pipefail

# Directory where .txt files with accession lists live
ROOT_DIR="/data/stefan/stress_granules/paired/Curdy_dataset"

# Make sure ROOT_DIR exists
if [[ ! -d $ROOT_DIR ]]; then
  echo "Error: ROOT_DIR '$ROOT_DIR' does not exist." >&2
  exit 1
fi

# Loop over every .txt file in ROOT_DIR
for list_path in "$ROOT_DIR"/*.txt; do
  # If no .txt files match, echo msg
  if [[ ! -e $list_path ]]; then
    echo "No .txt files found in $ROOT_DIR"
    exit 0
  fi

  echo "Processing accession list: $list_path"
  while IFS= read -r accession || [[ -n $accession ]]; do
    # skip empty lines and comments
    [[ -z $accession || $accession =~ ^# ]] && continue
    echo "  Downloading $accession into $ROOT_DIR …"
    fastq-dump --split-files --outdir "$ROOT_DIR" "$accession"
  done < "$list_path"
done


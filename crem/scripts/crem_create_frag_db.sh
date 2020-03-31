#!/usr/bin/env bash

# $1 - input SMILES file
# $2 - output dir which will contain all intermediate files as well as fragment.db file
# $3 - number of CPUs to use

if [[ $# -le 2 || $# -gt 3 ]]; then
  echo "Incorrect number of arguments:
  the first argument should be input SMILES file;
  the second argument should be output directory where all intermediate files as well as output fragment.db file will be stored;
  the third argument is a number of CPUs to use (optional, default 1)."
  exit 1
fi

if [[ $# -eq 3 ]]; then
  CPU=$3
else
  CPU=1
fi

INPUT_FILE=$1
OUTPUT_DIR=$2

mkdir -p $2

fragmentation -i $INPUT_FILE -o $OUTPUT_DIR/frags.txt -c $CPU -v
# sort -o $OUTPUT_DIR/frags.txt $OUTPUT_DIR/frags.txt

for i in 1 2 3; do
  frag_to_env -i $OUTPUT_DIR/frags.txt -o $OUTPUT_DIR/r$i.txt -r $i -c $CPU -v
#  sort -o $OUTPUT_DIR/r$i.txt $OUTPUT_DIR/r$i.txt
  sort $OUTPUT_DIR/r$i.txt | uniq -c > $OUTPUT_DIR/r${i}_c.txt
  env_to_db -i $OUTPUT_DIR/r${i}_c.txt -o $OUTPUT_DIR/fragments.db -r $i -c -v
done


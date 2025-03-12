#!/bin/bash

if [ $# -eq 0 ]; then
  echo "No GFA file spencified"
  echo "Usage: $0 GFA_FILE"
  exit
fi

if [ $# -gt 1 ]; then
  echo "Too many argument"
  echo "Usage: $0 GFA_FILE"
  exit
fi

inp_file="$1"
new_file="${inp_file%.*}-ap.gfa"
cp -v "$inp_file" "$new_file"
awk '/^S/ && $5 ~ /SN:Z:chr/ && $6 ~ /SO:i:0/ {sub("SN:Z:", "", $5); sub("LN:i:", "", $4); print "P\t" $5 "\t" $2 "+\t" $4 "M"}' "$inp_file" >> "$new_file"

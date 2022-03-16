#!/bin/bash

# read -p "Enter Directory to read the fasta files: " input_dir

# echo "INPUT DIR: $input_dir"

# read -p "Enter Directory to save the output: " output_dir

# echo "OUTPUT DIR: $output_dir"
input_dir="input_fasta"
output_dir="test_output"
# command="python3 main.py $filename $output_dir"

for file in input_fasta/*; do python3 main.py "$file" "$output_dir"; done
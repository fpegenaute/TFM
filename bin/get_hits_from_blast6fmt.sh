#!/usr/bin/bash

# Headers for the format 5 in BLAST:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

# Print the lines whose evalue=0 (perfect identity)
awk -F'\t' '$11=="0.0" {print $0}' test_blast.out > exact_matches.out


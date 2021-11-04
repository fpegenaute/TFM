#!/usr/bin/bash

# Get all the templates for a query from the reports of alphafold
grep -w -R "template match" *.out reports | awk 'NF>1{print $NF}'| sed 's/.$//' | awk '!seen[$0]++' | sed 's/_/\t/g'

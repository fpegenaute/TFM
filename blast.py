from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
import os

# Set custom db path
# export BLASTDB="/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa"

os.environ["BLASTDB"]="/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa"

db="pdbaa"
fasta_dir="./templates/FASTA/"
fasta_dir="./input_fasta/"
query="CREB"
extension=".fa"

with open(fasta_dir+query+extension, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        print(record.id)
        print(record.seq)

blastp_cline = NcbiblastpCommandline(query=fasta_dir+query+extension, db=db, matrix="BLOSUM80",outfmt=6, out=query+"_blast.out")

stdout, stderr = blastp_cline()

# Catch hits with evalue 0 and equal length as the query
# Headers for the format 5 in BLAST:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore


from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
import os
from Bio.Blast import NCBIXML

def run_blast_local(blastdb, db, fasta_dir, query, extension=".fa"):
    """
    Run BLAST Locally
    """

    # Set custom db path
    os.environ["BLASTDB"]=blastdb
    
    # Call BLAST
    blastp_cline = NcbiblastpCommandline(query=fasta_dir+query+extension, db=db, matrix="BLOSUM80",outfmt=5, out=query+"_blast.out")
    stdout, stderr = blastp_cline()

    outblast = "The results are stored in "+query+"_blast.out"
    return outblast





def exact_match_retriever(filename):
    """
    Given a XML (format 5) formatted BLAST output, 
    return a list of exact matches, meaning e-value == 0.0
    AND length of the query == length of the match.
    """

    with open(filename) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    # Get exact matches Eval+length
    matches = []
    E_VALUE_THRESH = 0.0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect == E_VALUE_THRESH:
                if alignment.length == len(hsp.query):
                    ID = alignment.title[4:8]
                    matches.append(ID)
    return(matches)


if __name__ == "__main__":
    # Set vars for BLAST
    blastdb = "/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa"
    db="pdbaa"
    fasta_dir="./templates/FASTA/"
    fasta_dir="./input_fasta/"
    query="test"
    extension=".fa"

    # Run BLAST
    outblast =run_blast_local(blastdb, db, fasta_dir, query)
    print(outblast)

    matches = exact_match_retriever("test_blast.out")
    print(matches)
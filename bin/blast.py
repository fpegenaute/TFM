from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
import os
from Bio.Blast import NCBIXML
from pathlib import Path

def run_blast_local(fasta, blastdb, db="pdbaa"):
    """
    Run BLAST Locally

    - blastdb = Path were the BLAST databases are stored (echo $BLASTDB)
    - fasta_dir = Path where your target fasta sequence is
    - query = name of the query (has to correspond to the name of the fasta file)
    - db = Which database to search against (default pdbaa)

    Returns a string announcing where the results are
    """

    # Extract query name from filename
    query = Path(fasta).stem

    # Set environment variable to the DB path
    os.environ["BLASTDB"]=blastdb
    
    # Call BLAST
    blastp_cline = NcbiblastpCommandline(query=fasta, 
    db=db, matrix="BLOSUM80",outfmt=5, out=query+"_blast.out")
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
    
    fasta= "test.fa"
  
    # Run BLAST
    outblast =run_blast_local(fasta, blastdb)
    print(outblast)

    # matches = exact_match_retriever("test_blast.out")
    # print(matches)
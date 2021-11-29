from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
import os
from Bio.Blast import NCBIXML
from pathlib import Path

from Bio.PDB.Dice import extract

def run_blast_local(fasta, blastdb, outdir, db="pdbaa"):
    """
    Run BLAST Locally

    - blastdb = Path were the BLAST databases are stored (echo $BLASTDB)
    - fasta_dir = Path where your target fasta sequence is
    - query = name of the query (has to correspond to the name of the fasta file)
    - db = Which database to search against (default pdbaa)

    Returns a string announcing where the results are
    """

    # Extract query name from filename
    query = Path(fasta).stem.split('.')[0]

    # Set environment variable to the DB path
    os.environ["BLASTDB"]=blastdb
    
    # Call BLAST
    blastp_cline = NcbiblastpCommandline(query=fasta, 
    db=db, matrix="BLOSUM80",outfmt=5, out=os.path.join(outdir, f"{query}_blast.out"), evalue=0.00000005)
    stdout, stderr = blastp_cline()

    outblast = os.path.join(outdir, f"{query}_blast.out")
    return outblast



def exact_match_retriever(filename):
    """
    Given a XML (format 5) formatted BLAST output, 
    return a dictionary of exact matches, meaning e-value == 0.0
    AND length of the query == length of the match.
    
    Each dictionary entry has the format "PDB ID : Chain"
    """

    with open(filename) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    # Get exact matches Eval+length
    matches = {}
    E_VALUE_THRESH = 0.00000005
    i = 0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if i < 4:
                if alignment.length == len(hsp.query):
                    ID = alignment.title[4:8]
                    Chain = alignment.title[9]
                    i+=1

                    matches.update({ID : Chain})
    return(matches)


if __name__ == "__main__":
    # Set vars for BLAST
    blastdb = "/home/gallegolab/Desktop/TFM/databases/BLAST/pdbaa"
    
    fasta= "input_fasta/test.fa"
  
    # Run BLAST
    outblast =run_blast_local(fasta, blastdb)
    print(outblast)
    matches = exact_match_retriever(outblast)
    print(matches)

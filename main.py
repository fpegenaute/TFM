###############################################################################
# Developed by Ferran Pegenaute as part of the Master Thesis (MSc in Bioinfor-
################## matics for Health Sciences, UPF) ###########################
###################### ferran.pegenaute@upf.edu ###############################
#################### ferran.pegenaute01@gmail.com #############################

import Bio.SeqIO as IO
import argparse
from bin.blast import *
from bin.PDB_retriever import *
import  bin.config as cfg
import logging as l
from pathlib import Path, PurePath, PurePosixPath
import os
import shutil
from bin.extract_flexible_residues import extract_residue_list
from bin.process_predicted_model import *
from bin.graphical_summary import plot_coverage, plot_dfi_summary
from bin.dfi.DFI_plotter import run_dfi , extract_flexible_residues, plot_dfi, plot_peaks

parser = argparse.ArgumentParser(description="""This program retrieves
                        Structural information from a sequence in a fasta file
                        """, usage="main.py input_fasta outdir AF_preset [options]")

parser.add_argument("FASTA", 
                    help="Input sequence in FASTA format")
parser.add_argument("outdir", 
                    help="Output directory to store the retrieved PDBs", 
                    default=".")
# parser.add_argument("-m", "--model_preset", 
#                     help="model preset for AlphaFold2", 
#                     default="monomer")
parser.add_argument("-v", "--verbose", 
                    help="Increase output verbosity", 
                    action="store_true")
# parser.add_argument("-c", "--custom", 
#                     help="Use custom slurm batch script for AlphaFold2")
# parser.add_argument("-n", "--noslurm",
#                     help="""run AlphaFold2 locally, without using SLURM.
#                     Only recommended if this python script is submitted on a 
#                     batch script on its own""", 
#                     default=False, 
#                     action="store_true")




args = parser.parse_args()
fasta = args.FASTA
query_name = Path(fasta).stem.split('.')[0]

### Initializing the LOG system ###

logdir = os.path.join(args.outdir, query_name, "LOG", "")
Path(logdir).mkdir(parents=True, exist_ok=True)

l.basicConfig(format = "%(levelname)s:%(message)s", 
                        filename = os.path.join(logdir, f"{query_name}.log"), 
                        level = l.DEBUG)

l.debug("...STARTING...\n")		

# If verbose is set, the LOG file is also printed in STDOUT
if args.verbose:		
	l.getLogger().addHandler(l.StreamHandler())		


## Info about the query

record_dict = IO.to_dict(IO.parse(fasta, "fasta"))
if len(record_dict.keys()) > 1:
    raise NotImplemented("This program only works with one sequence at a time")
if len(record_dict.keys()) == 1:
    if list(record_dict.keys())[0] != query_name:
        raise NameError(f"""Please, make sure your filename and fasta identifier 
        coincide. filename: {query_name} / ID name: {record_dict.keys()}""")
    
    query_length = len(record_dict[query_name].seq)
    l.info(f"Query length: {query_length}")


# For now, it only forks with one squence at a time
# for key in record_dict.items():
#     print(key[0],"\n ",len(key[1].seq))

# Create folders. last path element is empty to add a slash
pdb_dir = os.path.join(args.outdir, query_name, "PDB", "" )
fasta_dir = os.path.join(args.outdir, query_name, "FASTA", "" )

Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(fasta_dir).mkdir(parents=True, exist_ok=True)

## 1. Check if the input sequence is already in the PDB  

# Locate the Database
blastdb = cfg.blastconfig["blastdb"]
l.info(f"BLAST database is located at: {blastdb}")

# Create folders. last path element is empty to add a slash
blast_dir = os.path.join(args.outdir, query_name,  "BLAST", "" )
Path(blast_dir).mkdir(parents=True, exist_ok=True)
l.info(f"The BLAST output will be stored in:{blast_dir}")

# Run BLAST
l.info(f"Starting BLAST. Query: {fasta}, Database: {Path(blastdb).stem}")
outblast = run_blast_local(fasta, blastdb, blast_dir)
l.info(f"BLAST results stored in : {outblast}")

# Catch exact matches
exact_matches = exact_match_retriever(outblast)
l.info(f""" The target sequence has close homologs in the PDB with 
    code/s: {exact_matches.keys()}""")

structures_for_query = []
# Retrieve exact matches from the PDB

# Retrieve
if exact_matches:
    retrieve_pdb_info(exact_matches, pdb_dir, fasta_dir)
    # Check lengths of the actual PDB Chains and store them accordingly
    for file in os.listdir(pdb_dir):
        current = os.path.join(pdb_dir, file)
        if os.path.isfile(current):
            l.info(f"File being processed: {file}")
            identifier = file.split(".")[0].upper()
            # Check which chain of the hit should it choose
            
            # Make the directory for the chains
            chain_dir = os.path.join(pdb_dir, "CHAINS", "" )
            l.info(f"Making directory for the chains at {chain_dir}")
            Path(chain_dir).mkdir(parents=True, exist_ok=True)
            
            # Extract the desired chain
            l.info(f"Extracting the chain")
            splitter = ChainSplitter(mmcif=True, out_dir=os.path.join(pdb_dir, "CHAINS"))
            chain_path = splitter.make_pdb(os.path.join(pdb_dir, file), exact_matches[identifier] )
            structures_for_query.append(chain_path)
            l.debug(f"CHAIN PATH: {chain_path}")

            
            # Store partial matches (<95% of the query length)
            pdb_len = check_PDB_len(chain_path, exact_matches[identifier])
            l.info(f"Length of the template {PurePosixPath(chain_path).name}: {pdb_len}")
            
            l.info(f"PDB_LEN: {pdb_len} . QUERY_LEN: {0.95*query_length}")
            if pdb_len > 10 and pdb_len < (query_length):
                l.info(f"""{PurePosixPath(chain_path).name} has length {pdb_len}, it will be stored 
                    as a partial match""")
                try:
                    shutil.move(os.path.join(pdb_dir, file), os.path.join(pdb_dir,"partial", file))
                except Exception:
                    directory = os.path.join(pdb_dir,"partial")
                    l.info(f"\"{directory}\" does not exist, it will be created")
                    os.mkdir(os.path.join(pdb_dir,"partial"))
                    shutil.move(os.path.join(pdb_dir, file), os.path.join(pdb_dir,"partial", file))
            if pdb_len > 10 and pdb_len > (0.95*query_length):
                l.info(f"""{PurePosixPath(chain_path).name} has length {pdb_len}, it will be stored 
                    as a full-length match""")
                try:
                    shutil.move(os.path.join(pdb_dir, file), os.path.join(pdb_dir,"total", file))
                    structures_for_query.append(os.path.join(pdb_dir,"total", file))
                except Exception:
                    directory = os.path.join(pdb_dir,"total")
                    l.info(f"\"{directory}\" does not exist, it will be created")
                    os.mkdir(os.path.join(pdb_dir,"total"))
                    shutil.move(os.path.join(pdb_dir, file), os.path.join(pdb_dir,"total", file))
                    structures_for_query.append(os.path.join(pdb_dir,"total", file))




### Submit a sob in Slurm with the AlphaFold run

l.debug("Alphafold will not run, this is a test")


af_dir = os.path.join(args.outdir,  query_name, "ALPHAFOLD", "" )
Path(af_dir).mkdir(parents=True, exist_ok=True)
# Make folder for the AF2 output
l.info(f"Creating folder for AF2 output in:{af_dir}")

# if args.custom:
#     slurm_script = args.custom
# else:
#     slurm_dir = f"{af_dir}"
#     Path(slurm_dir).mkdir(parents=True, exist_ok=True)





### Extract confident regions


# Setting up the parameters for the PHENIX library
master_phil = iotbx.phil.parse(master_phil_str)
params = master_phil.extract()
master_phil.format(python_object=params).show(out=sys.stdout)
p = params.process_predicted_model

p.b_value_field_is = 'lddt'
p.domain_size = 15
p.remove_low_confidence_residues = True
p.maximum_rmsd = 1.5
p.split_model_by_compact_regions = True

from iotbx.data_manager import DataManager
dm = DataManager()
dm.set_overwrite(True)

## PROVISIONAL
PAE_dir = os.path.join(af_dir,  "PAE", "" )
Path(PAE_dir).mkdir(parents=True, exist_ok=True)
PAE_json = os.path.join(PAE_dir, "sec3_PAE.json")


l.info("Extracting high confidence domains")
domains_dir = os.path.join(af_dir, "DOMAINS", "")
Path(domains_dir).mkdir(parents=True, exist_ok=True)
l.info(f"Domains will be stored in:{domains_dir}")
af_conficent_regions = []

for filename in os.listdir(af_dir):
    if os.path.isfile(os.path.join(af_dir, filename)):
        l.info(f"Processing file: {filename}")
        print("\nProcessing and splitting model into domains")
        
        m = dm.get_model(os.path.join(af_dir, filename))
        pae_matrix = pae_matrix = parse_json_PAE(PAE_json )
        model_info = process_predicted_model(m,  params, pae_matrix)

        chainid_list = model_info.chainid_list
        print("Segments found: %s" %(" ".join(chainid_list)))

        mmm = model_info.model.as_map_model_manager()
        mmm.write_model(os.path.join(domains_dir, f"{PurePosixPath(filename).stem}_domains.pdb"))
        
        structures_for_query.append(os.path.join(domains_dir, f"{PurePosixPath(filename).stem}_domains.pdb"))

        conf_domains = extract_residue_list(os.path.join(domains_dir, f"{PurePosixPath(filename).stem}_domains.pdb"), domains_dir)
        l.info(f"Residue list of confident domains: {conf_domains}")




  # START DFI   
   
from bin.graphical_summary import plot_dfi_summary
    
plot_dfi_summary(structures_for_query)



## Launch graphical summary 
print(f"CONFIDENT FILES: {structures_for_query}")
plot_coverage(fasta, structures_for_query)




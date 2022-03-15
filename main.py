###############################################################################
# Developed by Ferran Pegenaute as part of the Master Thesis (MSc in Bioinfor-
################## matics for Health Sciences, UPF) ###########################
###################### ferran.pegenaute@upf.edu ###############################
#################### ferran.pegenaute01@gmail.com #############################

from doctest import REPORT_CDIFF
from re import T
import matplotlib.pyplot as plt

import Bio.SeqIO as IO
import argparse
from bin.blast import *
from bin.PDB_retriever import *
import  bin.config as cfg
import logging as l
from pathlib import Path, PurePosixPath
import os
import shutil
from bin.extract_flexible_residues import extract_residue_list
from bin.process_predicted_model import *
from bin.graphical_summary import  plot_dfi_hinge_summary
from matplotlib import pyplot as plt
import fnmatch
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError





parser = argparse.ArgumentParser(description="""This program retrieves
                        Structural information from a sequence in a fasta file
                        """, usage="main.py input_fasta outdir [options]")

parser.add_argument("FASTA", 
                    help="Input sequence in FASTA format")
parser.add_argument("outdir", 
                    help="Output directory to store the retrieved PDBs", 
                    default=".")
parser.add_argument("-af", "--alphamodel",  nargs='?',
                    help="AlphaFold2 model in PDB format", 
                    default=None)
parser.add_argument("-rf", "--rosettamodel",  nargs='?',
                    help="RoseTTaFold model in PDB format", 
                    default=None)
parser.add_argument("-j", "--PAE_json",  nargs='?',
                    help="AlphaFold2 PAE JSON file from the AF-EBI server",
                    default=None)
parser.add_argument("-c", "--custom_templates",  nargs='?',
                    help="A custom experimentally solved PDB provided by the user",
                    default=None)
parser.add_argument("-ra", "--run_alphafold", 
                    help="Send an batch script using SLURM (you need to be in a cluster with slurm and AF2 installed)", 
                    action="store_true")
parser.add_argument("-rr", "--run_rosettafold", 
                    help="Send an batch script using SLURM (you need to be in a cluster with slurm and RoseTTafold installed)", 
                    action="store_true")
parser.add_argument("-v", "--verbose", 
                    help="Increase output verbosity", 
                    action="store_true")


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


# Create folders. last path element is empty to add a slash
pdb_dir = os.path.join(args.outdir, query_name, "PDB", "" )
fasta_dir = os.path.join(args.outdir, query_name, "FASTA", "" )
report_dir = os.path.join(args.outdir, query_name, "REPORT", "" )
hinges_dir = os.path.join(args.outdir, query_name, "HINGES", "" )
IMP_dir = os.path.join(args.outdir, query_name, "IMP", "" )

Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(fasta_dir).mkdir(parents=True, exist_ok=True)
Path(report_dir).mkdir(parents=True, exist_ok=True)
Path(hinges_dir).mkdir(parents=True, exist_ok=True)
Path(IMP_dir).mkdir(parents=True, exist_ok=True)


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
            splitter = ChainSplitter(mmcif=True, out_dir=chain_dir)
            chain_path = splitter.make_pdb(os.path.join(pdb_dir, file), 
                                    exact_matches[identifier], overwrite=True )
            # structures_for_query.append(chain_path)
            # print(f"STRUCTURES FOR QUERY: {structures_for_query}")
            
            # Store partial matches (<95% of the query length)
            pdb_len = check_PDB_len(chain_path, exact_matches[identifier])
            l.info(f"Length of the template {PurePosixPath(chain_path).name}: {pdb_len}")
            
            l.info(f"PDB_LEN: {pdb_len} . QUERY_LEN: {query_length}")
            if pdb_len > 10 and pdb_len < (0.95*query_length):
                l.info(f"""{PurePosixPath(chain_path).name} has length {pdb_len}, it will be stored 
                    as a partial match""")
                newpath = os.path.join(pdb_dir,"partial", f"{PurePosixPath(chain_path).name}")
                try:
                    l.info(f"MOVING {chain_path} TO {newpath}")
                    shutil.move(chain_path, newpath)
                    structures_for_query.append(newpath)
                except Exception:
                    directory = os.path.join(pdb_dir,"partial", "")
                    l.info(f"\"{directory}\" does not exist, it will be created")
                    Path(directory).mkdir(parents=True, exist_ok=True)
                    shutil.move(chain_path, newpath)
                    structures_for_query.append(newpath)
            if pdb_len > 10 and pdb_len > (0.95*query_length):
                l.info(f"""{PurePosixPath(chain_path).name} has length {pdb_len}, it will be stored 
                    as a full-length match""")
                newpath = os.path.join(pdb_dir,"total", f"{PurePosixPath(chain_path).name}")
                try:
                    shutil.move(chain_path, newpath)
                    structures_for_query.append(newpath)
                except Exception:
                    directory = os.path.join(pdb_dir, "total", "")
                    l.info(f"\"{directory}\" does not exist, it will be created")
                    Path(directory).mkdir(parents=True, exist_ok=True)          
                    shutil.move(chain_path, newpath)
                    structures_for_query.append(newpath)

# Don't forget the user's templates!
if args.custom_templates:
        structures_for_query.append(args.custom_templates)
        


### ALPHAFOLD & PAE

# Make folder for the AF2 output
af_dir = os.path.join(args.outdir,  query_name, "ALPHAFOLD", "" )
l.info(f"Creating folder for AF2 output in:{af_dir}")
Path(af_dir).mkdir(parents=True, exist_ok=True)

# Make the directory for PAE
PAE_dir = os.path.join(af_dir,  "PAE", "" )
Path(PAE_dir).mkdir(parents=True, exist_ok=True)

   
# If you want to use your AF model and PAE file:
if args.alphamodel:
    # Store the AF model and the PAE file in the correct folders
    AF_server_model = args.alphamodel
    shutil.copy(AF_server_model, os.path.join(af_dir, PurePosixPath(AF_server_model).name))
if args.PAE_json:
    PAE_json = args.PAE_json

from bin.utilities import submit_AF_to_SLURM, submit_RF_to_SLURM

if args.run_alphafold:
    submit_AF_to_SLURM(fasta, af_dir, workload_manager="sbatch", dummy_dir=".", max_jobs_in_queue=None )


### ROSETTAFOLD

# Make folder for the RFoutput
rf_dir = os.path.join(args.outdir,  query_name, "ROSETTAFOLD", "" )
l.info(f"Creating folder for RoseTTaFold output in:{rf_dir}")
Path(rf_dir).mkdir(parents=True, exist_ok=True)


   
# If you want to use your RF model and PAE file:
if args.rosettamodel:
    # Store the AF model and the PAE file in the correct folders
    RF_custom_model = args.rosettamodel
    shutil.copy(RF_custom_model, os.path.join(rf_dir, PurePosixPath(RF_custom_model).name))
# If you want to send a RF job to a HPC Cluster
if args.run_rosettafold:
    submit_RF_to_SLURM(fasta, rf_dir, workload_manager="sbatch", dummy_dir=".", max_jobs_in_queue=None)

### Extract confident regions


# Setting up the parameters for the PHENIX library
master_phil = iotbx.phil.parse(master_phil_str)
params = master_phil.extract()
master_phil.format(python_object=params).show(out=sys.stdout)
p = params.process_predicted_model


p.domain_size = 15
p.remove_low_confidence_residues = True
p.maximum_rmsd = 1.5
p.split_model_by_compact_regions = True

from iotbx.data_manager import DataManager
dm = DataManager()
dm.set_overwrite(True)


l.info("Extracting high confidence domains")
domains_dir = os.path.join(af_dir, "DOMAINS", "")
Path(domains_dir).mkdir(parents=True, exist_ok=True)
l.info(f"Domains will be stored in:{domains_dir}")
af_conficent_regions = []

if (args.alphamodel and args.PAE_json) or (args.run_alphafold):
    p.b_value_field_is = 'lddt'
    for filename in os.listdir(af_dir):
        if os.path.isfile(os.path.join(af_dir, filename)):
            l.info(f"Processing file: {filename}")
            print("\nProcessing and splitting model into domains")
            
            m = dm.get_model(os.path.join(af_dir, filename))
            pae_matrix = pae_matrix = parse_json_PAE(PAE_json)
            model_info = process_predicted_model(m,  params, pae_matrix)

            chainid_list = model_info.chainid_list
            print("Segments found: %s" %(" ".join(chainid_list)))

            mmm = model_info.model.as_map_model_manager()
            
            # Write all the domains in one file
            mmm.write_model(os.path.join(domains_dir, 
                            f"{PurePosixPath(filename).stem}_domains.pdb"))
            
            # Write different domains in different files
            for chainid in chainid_list:
                selection_string = "chain %s" %(chainid)
                ph = model_info.model.get_hierarchy()
                asc1 = ph.atom_selection_cache()
                sel = asc1.selection(selection_string)
                m1 = model_info.model.select(sel)
                # dm.write_model_file(m1, '%s_%s.pdb' %(output_file_name[:-4],chainid))
                filepath = os.path.join(domains_dir, 
                        f"{PurePosixPath(filename).stem}_{chainid}_AF.pdb")
                dm.write_model_file(m1, filepath)
                structures_for_query.append(filepath)
                

            
            # structures_for_query.append(os.path.join(domains_dir, 
            #                 f"{PurePosixPath(filename).stem}_domains.pdb"))

            conf_domains = extract_residue_list(os.path.join(domains_dir,  
                            f"{PurePosixPath(filename).stem}_domains.pdb"), 
                            domains_dir)
            l.info(f"Residue list of confident domains: {conf_domains}")


## Launch graphical summary 
l.info(f"CONFIDENT FILES: {structures_for_query}")
nrow = len(structures_for_query)
l.info(f"NROW: {nrow}")
# plot_coverage(fasta, structures_for_query, nrow)

plt.show()

### HINGE DETECTION ###
from bin.graphical_summary import StructuReport


for structure in structures_for_query:
    reporter = StructuReport(structure)
    # Get coverage of the structure
    coverage_df = reporter.get_coverage(fasta, save_csv=True, outdir=report_dir)
    # Get hinges and save the .hng files
    hinges = reporter.get_hinges(alpha_range=None, 
                                            save_hng=True, outdir=hinges_dir)
    # Get DFI
    dfi_df = reporter.get_dfi_coverage(reference_fasta=fasta, save_csv=True, outdir=report_dir)
                                    
  
# plot_dfi_hinge_summary(structures_for_query, fasta, hinges_dir)


## Write Topology file
from bin.custom_top import RigidBody, write_custom_topology
from bin.utilities import get_chain_names, get_residue_range



# Make rhe RigidBody objects from the structures
i = 1
rigid_bodies = []
for structure in structures_for_query:
    
    filename, extension = get_filename_ext(structure)
    
    # Extract chain name
    chain_IDs = get_chain_names(structure)
    if len(chain_IDs) > 1 or fnmatch.fnmatch(structure, "*AF.pdb"):
        print(f"""{structure}: Assuming a AlphaFold model""")

        for chain in chain_IDs:
            # Extract residue range
            res_range = get_residue_range(structure, chain=chain)    
            # Create the RigidBody instance
            rigid_body = RigidBody(resolution="all",
            molecule_name= f"{filename}_{chain}", 
            color="orange" , 
            fasta_fn=fasta, 
            # fasta_id=fasta, 
            pdb_fn=structure, 
            chain=chain,
            residue_range=res_range , 
            rigid_body=i, 
            super_rigid_body="", 
            chain_of_super_rigid_bodies="", 
            bead_size=20,
            em_residues_per_gaussian=0, 
            type="AF_model")
            # Add the rigid body to a list
            rigid_bodies.append(rigid_body)
            i +=1
    
    if  fnmatch.fnmatch(structure, "*RF.pdb"):
        print(f"""{structure}: Assuming a RoseTTaFold model""")

        for chain in chain_IDs:
            # Extract residue range
            res_range = get_residue_range(structure, chain=chain)    
            # Create the RigidBody instance
            rigid_body = RigidBody(resolution="all",
            molecule_name= f"{filename}_{chain}", 
            color="orange" , 
            fasta_fn=fasta, 
            # fasta_id=fasta, 
            pdb_fn=structure, 
            chain=chain,
            residue_range=res_range , 
            rigid_body=i, 
            super_rigid_body="", 
            chain_of_super_rigid_bodies="", 
            bead_size=20,
            em_residues_per_gaussian=0, 
            type="RF_model")
            # Add the rigid body to a list
            rigid_bodies.append(rigid_body)
            i +=1
    if len(chain_IDs) > 1:
        print(f"Skipping {structure}, since it has more than one chain")



    # Extract residue range
    res_range = get_residue_range(structure)    
    # Create the RigidBody instance
    rigid_body = RigidBody(resolution="all",
    molecule_name= filename, 
    color="blue" , 
    fasta_fn=fasta, 
    # fasta_id=fasta, 
    pdb_fn=structure, 
    chain=chain_IDs[0],
    residue_range=res_range , 
    rigid_body=i, 
    super_rigid_body="", 
    chain_of_super_rigid_bodies="", 
    bead_size=10,
    em_residues_per_gaussian=0, 
    type="experimental")

    # Add the rigid body to a list
    rigid_bodies.append(rigid_body)
    i +=1

## MAKE THE COMPOSITE
from bin.custom_top import make_composite
import pandas as pd

composite_rb = make_composite(rigid_bodies)


## TO DO: export the coverage composite csv 
i=0
for rb in composite_rb:
    coverage = rb.get_coverage()
    if i == 0:
        composite_coverage = coverage
        i+=1
    else:
        if i == 1:
            merged_left = pd.merge(left=composite_coverage, right=coverage, 
                        how="left", left_on="ResID", right_on="ResID")
        else:
            merged_left = pd.merge(left=merged_left, right=coverage, 
                        how="left", left_on="ResID", right_on="ResID")

out_path = os.path.join(report_dir, "COVERAGE", f"{query_name}_composite_coverage.csv")
merged_left.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
    
        
    



# Convert to list and sort by the ones who start earlier in the sequence
rigid_bodies.sort(key=lambda x: x.residue_range[0])

# Write the topology file
write_custom_topology(os.path.join(IMP_dir, f"{query_name}.topology"), rigid_bodies)

# Finally, write the automatic report
report_template = os.path.join(args.outdir, query_name, "report_template.ipynb")
shutil.copy("src/report_template.ipynb", report_template)



notebook_dir = os.path.join(args.outdir, query_name)
final_report = os.path.join(notebook_dir, f"{query_name}_report.ipynb")

with open(report_template) as f:
    nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=600)
    
    # Execute/run the notebook
    try:
        out = ep.preprocess(nb, {'metadata': {'path': notebook_dir}})
    except CellExecutionError:
        out = None
        msg = f"Error executing the notebook"
        msg += f"See notebook  for the traceback.'"
        print(msg)
        raise
    finally:
        with open(final_report, mode='w', encoding='utf-8') as f:
            nbformat.write(nb, f)
        os.remove(report_template)
        

import subprocess
subprocess.run(["jupyter", "nbconvert", "--to", "html", "--no-input", "--no-prompt ", f"{final_report}",])





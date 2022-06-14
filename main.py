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
from bin.graphical_summary import StructuReport
from bin.extract_flexible_residues import extract_residue_list
from bin.process_predicted_model import *
from matplotlib import pyplot as plt
import fnmatch
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
from bin.custom_top import make_rb_list, make_composite, write_custom_topology
import pandas as pd




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



record_dict = IO.to_dict(IO.parse(fasta, "fasta"))
if len(record_dict.keys()) > 1:
    raise NotImplemented("This program only works with one sequence at a time")
if len(record_dict.keys()) == 1:
    if list(record_dict.keys())[0] != query_name:
        raise NameError(f"""Please, make sure your filename and fasta identifier 
        coincide. filename: {query_name} / ID name: {record_dict.keys()}""")
    
    query_length = len(record_dict[query_name].seq)
    l.info(f"Query length: {query_length}")


l.info(f"Create folders:")
blast_dir = os.path.join(args.outdir, query_name,  "BLAST", "" )
pdb_dir = os.path.join(args.outdir, query_name, "PDB", "" )
fasta_dir = os.path.join(args.outdir, query_name, "FASTA", "" )
report_dir = os.path.join(args.outdir, query_name, "REPORT", "" )
hinges_dir = os.path.join(args.outdir, query_name, "HINGES", "" )
IMP_dir = os.path.join(args.outdir, query_name, "IMP", "" )
l.info(f"{pdb_dir}, {fasta_dir}, {report_dir}, {hinges_dir}, {IMP_dir} .")

Path(blast_dir).mkdir(parents=True, exist_ok=True)
Path(pdb_dir).mkdir(parents=True, exist_ok=True)
Path(fasta_dir).mkdir(parents=True, exist_ok=True)
Path(report_dir).mkdir(parents=True, exist_ok=True)
Path(hinges_dir).mkdir(parents=True, exist_ok=True)
Path(IMP_dir).mkdir(parents=True, exist_ok=True)




## 1. Check if the input sequence is already in the PDB  

l.info("### BLAST ###")

l.info(f"The BLAST output will be stored in:{blast_dir}")

# Run BLAST
outblast = run_blast_local(fasta, blast_dir)

# Catch exact matches
exact_matches = exact_match_retriever(outblast)
l.info(f""" The target sequence has close homologs in the PDB with 
    code/s: {exact_matches.keys()}""")

structures_for_query = []
l.info(f"Structures for query: {structures_for_query}")

# Retrieve from the PDB
l.info("Retrieving structures from the pdb")
if exact_matches:
    retrieve_pdb_info(exact_matches, pdb_dir, fasta_dir)
    # Check lengths of the actual PDB Chains and store them accordingly
    for file in os.listdir(pdb_dir):
        current = os.path.join(pdb_dir, file)
        if os.path.isfile(current):
            l.info(f"File being processed: {file}")
            identifier = file.split(".")[0].upper()
            
            # Make the directory for the chains
            chain_dir = os.path.join(pdb_dir, "CHAINS", "" )
            l.info(f"Making directory for the chains at {chain_dir}")
            Path(chain_dir).mkdir(parents=True, exist_ok=True)
            
            # Extract the desired chain
            l.info(f"Extracting the chain")
            splitter = ChainSplitter(mmcif=True, out_dir=chain_dir)
            chain_path = splitter.make_pdb(os.path.join(pdb_dir, file), 
                                    exact_matches[identifier], overwrite=True )
            
            # Store partial matches (<95% of the query length)
            pdb_len = check_PDB_len(chain_path, exact_matches[identifier])
            l.info(f"""Length of the template {PurePosixPath(chain_path).name}: 
                                {pdb_len}""")
            
            l.info(f"PDB_LEN: {pdb_len} . QUERY_LEN: {query_length}")
            if pdb_len > 10 and pdb_len < (0.95*query_length):
                l.info(f"""{PurePosixPath(chain_path).name} has length {pdb_len}, 
                it will be stored as a partial match""")
                newpath = os.path.join(pdb_dir,"partial", 
                                f"{PurePosixPath(chain_path).name}")
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
                l.info(f"""{PurePosixPath(chain_path).name} has length {pdb_len}, 
                it will be stored as a full-length match""")
                newpath = os.path.join(pdb_dir,"total", 
                                f"{PurePosixPath(chain_path).name}")
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
l.info("Checking user's templates")
if args.custom_templates:
        structures_for_query.append(args.custom_templates)
        

### ALPHAFOLD & PAE
  
# If you want to use your AF model and PAE file:
if args.alphamodel:
    l.info(f"custom AF model detected: {args.alphamodel}")
    af_dir = os.path.join(args.outdir,  query_name, "ALPHAFOLD", "" )
    Path(af_dir).mkdir(parents=True, exist_ok=True)
    AF_server_model = args.alphamodel
    shutil.copy(AF_server_model, os.path.join(af_dir, 
                    PurePosixPath(AF_server_model).name))
if args.PAE_json:
    l.info(f"custom PAE matrix detected: {args.PAE_json}")
    PAE_dir = os.path.join(af_dir,  "PAE", "" )
    Path(PAE_dir).mkdir(parents=True, exist_ok=True)
    PAE_json = args.PAE_json

from bin.utilities import submit_AF_to_SLURM, submit_RF_to_SLURM

if args.run_alphafold:
    l.info(f"Submitting AF2 SLURM batch script")
    af_dir = os.path.join(args.outdir,  query_name, "ALPHAFOLD", "" )
    Path(af_dir).mkdir(parents=True, exist_ok=True)
    submit_AF_to_SLURM(fasta, af_dir, workload_manager="sbatch", 
                    dummy_dir=".", max_jobs_in_queue=None )


### ROSETTAFOLD

   
# If you want to use your RF model:
if args.rosettamodel:
    l.info(f"custom RF model detected: {args.rosettamodel}")
    rf_dir = os.path.join(args.outdir,  query_name, "ROSETTAFOLD", "" )
    custom_rf_dir = os.path.join(rf_dir, "CUSTOM")
    Path(rf_dir).mkdir(parents=True, exist_ok=True)
    Path(custom_rf_dir).mkdir(parents=True, exist_ok=True)
    RF_custom_model = args.rosettamodel
    shutil.copy(RF_custom_model, os.path.join(custom_rf_dir, PurePosixPath(RF_custom_model).name))
if args.run_rosettafold:
    l.info(f"Submitting RF SLURM batch script")
    rf_dir = os.path.join(args.outdir,  query_name, "ROSETTAFOLD", "" )
    Path(rf_dir).mkdir(parents=True, exist_ok=True)
    submit_RF_to_SLURM(fasta, rf_dir, workload_manager="sbatch", dummy_dir=".", max_jobs_in_queue=None)

### Extract confident regions

if (args.alphamodel and args.PAE_json) or (args.run_alphafold):
    # Setting up the parameters for the PHENIX library
    master_phil = iotbx.phil.parse(master_phil_str)
    params = master_phil.extract()
    master_phil.format(python_object=params).show(out=sys.stdout)
    p = params.process_predicted_model
    p.domain_size = cfg.CCTBXconfig["AF2_domain_size"]
    p.maximum_rmsd = cfg.CCTBXconfig["AF2_maximum_rmsd"]
    p.b_value_field_is = 'lddt'

    from iotbx.data_manager import DataManager
    
    dm = DataManager()
    dm.set_overwrite(True)


    l.info("Extracting AF2 high confidence domains")
    domains_dir = os.path.join(af_dir, "DOMAINS", "")
    Path(domains_dir).mkdir(parents=True, exist_ok=True)
    l.info(f"Domains will be stored in:{domains_dir}")
    
    af_conficent_regions = []
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
                filepath = os.path.join(domains_dir, 
                        f"{PurePosixPath(filename).stem}_{chainid}_AF.pdb")
                dm.write_model_file(m1, filepath)
                structures_for_query.append(filepath)
                

            conf_domains = extract_residue_list(os.path.join(domains_dir,  
                            f"{PurePosixPath(filename).stem}_domains.pdb"), 
                            domains_dir)
            l.info(f"Residue list of confident domains: {conf_domains}")

# For Rosettafold models
if os.path.exists(os.path.join(args.outdir,  query_name, "ROSETTAFOLD", "" )):
    rf_dir = os.path.join(args.outdir,  query_name, "ROSETTAFOLD", "" )
    rf_models_dir = os.path.join(rf_dir, "model", "")
    # Setting up the parameters for the PHENIX library
    master_phil = iotbx.phil.parse(master_phil_str)
    params = master_phil.extract()
    master_phil.format(python_object=params).show(out=sys.stdout)
    p = params.process_predicted_model
    p.domain_size = cfg.CCTBXconfig["RF_domain_size"]
    p.maximum_rmsd = cfg.CCTBXconfig["RF_maximum_rmsd"]
    p.b_value_field_is = 'rmsd'

    from iotbx.data_manager import DataManager
    
    dm = DataManager()
    dm.set_overwrite(True)


    l.info("Extracting RoseTTaFold high confidence domains")
    domains_dir = os.path.join(rf_dir, "DOMAINS", "")
    Path(domains_dir).mkdir(parents=True, exist_ok=True)
    l.info(f"Domains will be stored in:{domains_dir}")
    abs_rf_dir = os.path.abspath(rf_models_dir)
    abs_custom_rf_dir = os.path.abspath(rf_models_dir)
    rf_conficent_regions = []
    for filename in os.listdir(abs_rf_dir):
        if os.path.isfile(os.path.join(abs_rf_dir,filename)) and \
            (fnmatch.fnmatch(filename, "model_*.crderr.pdb") or \
                fnmatch.fnmatch(filename, "model_*-crderr.pdb") ): 
            newname = filename.split(".")
            l.info(f"LIST NEWNAME: {newname}")
            noext = newname[0:-1]
            noext = "-".join(noext)
            ext = newname[-1]
            newname = noext+"."+ext
            l.info(f"NEWNAME: {newname}")
            filename = os.path.join(abs_rf_dir, filename)
            newname = os.path.join(abs_rf_dir, newname)
            os.rename(filename, newname)

            l.info(f"Processing file: {newname}")
            print("\nProcessing and splitting model into domains")
            
            m = dm.get_model(newname)
            model_info = process_predicted_model(m,  params)

            chainid_list = model_info.chainid_list
            print("Segments found: %s" %(" ".join(chainid_list)))

            mmm = model_info.model.as_map_model_manager()
            
            # Write all the domains in one file
            mmm.write_model(os.path.join(domains_dir, 
                            f"{PurePosixPath(newname).stem}_domains.pdb"))
            
            # Write different domains in different files
            for chainid in chainid_list:
                selection_string = "chain %s" %(chainid)
                ph = model_info.model.get_hierarchy()
                asc1 = ph.atom_selection_cache()
                sel = asc1.selection(selection_string)
                m1 = model_info.model.select(sel)
                filepath = os.path.join(domains_dir, 
                        f"{PurePosixPath(newname).stem}_{chainid}_RF.pdb")
                dm.write_model_file(m1, filepath)
                structures_for_query.append(filepath)
                   

            conf_domains = extract_residue_list(os.path.join(domains_dir,  
                            f"{PurePosixPath(newname).stem}_domains.pdb"), 
                            domains_dir)           

            conf_domains = extract_residue_list(os.path.join(domains_dir,  
                            f"{PurePosixPath(newname).stem}_domains.pdb"), 
                            domains_dir)
            l.info(f"Residue list of confident domains: {conf_domains}")

l.info(f"CONFIDENT FILES: {structures_for_query}")
nrow = len(structures_for_query)



l.info("### HINGE DETECTION and DFI ###")
for structure in structures_for_query:
    reporter = StructuReport(structure)
    # Get coverage of the structure
    coverage_df = reporter.get_coverage(fasta, save_csv=True, outdir=report_dir)
    # Get hinges and save the .hng files
    hinges = reporter.get_hinges(alpha_range=cfg.PACKMANconfig["alpha_range"], 
                                            save_hng=True, outdir=hinges_dir)
    # Get DFI
    dfi_df = reporter.get_dfi_coverage(reference_fasta=fasta, save_csv=True, outdir=report_dir)
                                    
  


## Write Topology file ##

rigid_bodies = make_rb_list(structures_for_query, fasta)


composite_rb = make_composite(rigid_bodies)


# Exprt the composite coverage in .csv
out_path = os.path.join(report_dir, "COVERAGE", f"{query_name}_composite_coverage.csv")
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
            i+=1
        else:
            merged_left = pd.merge(left=merged_left, right=coverage, 
                        how="left", left_on="ResID", right_on="ResID")
            i+=1

if i == 1:
    composite_coverage.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
elif i > 1:
    merged_left.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')






# Convert to list and sort by the ones who start earlier in the sequence
composite_rb.sort(key=lambda x: x.residue_range[0])

# Write the topology file
write_custom_topology(os.path.join(IMP_dir, f"{query_name}.topology"), composite_rb)
os.rmdir("obsolete")
os.remove("DCI_pymol_output.txt")

from make_dashboard import *

app.run_server()
exit(0)



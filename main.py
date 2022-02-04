###############################################################################
# Developed by Ferran Pegenaute as part of the Master Thesis (MSc in Bioinfor-
################## matics for Health Sciences, UPF) ###########################
###################### ferran.pegenaute@upf.edu ###############################
#################### ferran.pegenaute01@gmail.com #############################

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
from bin.graphical_summary import plot_coverage, plot_dfi_summary
from bin.graphical_summary import plot_dfi_summary
from bin.packman import predict_hinge, write_hng_file
from packman.anm import hdANM
from packman import molecule
from packman import ANM





parser = argparse.ArgumentParser(description="""This program retrieves
                        Structural information from a sequence in a fasta file
                        """, usage="main.py input_fasta outdir [options]")

parser.add_argument("FASTA", 
                    help="Input sequence in FASTA format")
parser.add_argument("outdir", 
                    help="Output directory to store the retrieved PDBs", 
                    default=".")
parser.add_argument("-a", "--alphamodel",  nargs='?',
                    help="AlphaFold2 model in PDB format", 
                    default=None)
parser.add_argument("-j", "--PAE_json",  nargs='?',
                    help="AlphaFold2 PAE JSON file from the AF-EBI server",
                    default=None)
parser.add_argument("-r", "--run_alphafold", 
                    help="Send an batch script using SLURM (you need to be in a cluster with slurm and AF2 installed)", 
                    action="store_true")
parser.add_argument("-v", "--verbose", 
                    help="Increase output verbosity", 
                    action="store_true")


# parser.add_argument("-m", "--model_preset", 
#                     help="model preset for AlphaFold2", 
#                     default="monomer")
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
            splitter = ChainSplitter(mmcif=True, out_dir=chain_dir)
            chain_path = splitter.make_pdb(os.path.join(pdb_dir, file), exact_matches[identifier], overwrite=True )
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

from bin.utilities import submit_AF_to_SLURM

if args.run_alphafold:
    submit_AF_to_SLURM(fasta, af_dir, workload_manager="sbatch", dummy_dir=".", max_jobs_in_queue=None )


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





l.info("Extracting high confidence domains")
domains_dir = os.path.join(af_dir, "DOMAINS", "")
Path(domains_dir).mkdir(parents=True, exist_ok=True)
l.info(f"Domains will be stored in:{domains_dir}")
af_conficent_regions = []

if (args.alphamodel and args.PAE_json) or (args.run_alphafold):
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
            mmm.write_model(os.path.join(domains_dir, f"{PurePosixPath(filename).stem}_domains.pdb"))
            
            structures_for_query.append(os.path.join(domains_dir, f"{PurePosixPath(filename).stem}_domains.pdb"))

            conf_domains = extract_residue_list(os.path.join(domains_dir, f"{PurePosixPath(filename).stem}_domains.pdb"), domains_dir)
            l.info(f"Residue list of confident domains: {conf_domains}")




## Launch graphical summary 
l.info(f"CONFIDENT FILES: {structures_for_query}")
nrow = len(structures_for_query)
l.info(f"NROW: {nrow}")
plot_coverage(fasta, structures_for_query, nrow)
plot_dfi_summary(structures_for_query, fasta)

plt.show()

### HINGE DETECTION ###



for structure in structures_for_query:
    Protein = molecule.load_structure(structure)
    filename, ext = get_filename_ext(structure)
    try:
        Protein[0]
    except Exception:
        print("Make sure your filename is  of the form: XXXXX.pdb/XXXX.cif")

    chains = [chain for chain in Protein[0].get_chains()]
    backbone = [j for i in Protein[0][chains[0].get_id()].get_backbone() for j in i if j is not None]

        
    
    ##### For running iteratively several values of alpha:
    # alpha_start, alpha_stop, step_size = 2.5 , 4.5 , 0.5 # Previously from 1 to 10
    # for i in np.arange(alpha_start, alpha_stop, step_size):
    #     i = np.around(i, decimals=1)
    #     try:
    #         predict_hinge(backbone, Alpha=i, outputfile=open(str(i)+'.txt', 'w'))
    #         # predict_hinge(backbone, outfile, Alpha=4,method='alpha_shape',filename='Output.pdb',MinimumHingeLength=5,nclusters=2)

    #     except:
    #         continue    

    predict_hinge(backbone, Alpha=3.65, outputfile=open(str(f"{filename}_packman_output")+'.txt', 'w'))

    hinges = []
    hinges_nosig = []
    l.info("Significant hinges")
    for hinge in backbone[0].get_parent().get_parent().get_hinges():
        resids = [x.get_id() for x in hinge.get_elements()]
        if hinge.get_pvalue() < 0.05: 
            hinges.append(hinge)
            print(f"""HINGE {hinge.get_id()}
            \t p-value: {hinge.get_pvalue()}
            \t alpha-value: {hinge.get_alpha_value()}
            \t Location: {resids[0]} - {resids[-1]}  
            \t Length: {resids[-1] - resids[-0]}
            """)
        else:
            hinges_nosig.append(hinge)

    write_hng_file(structure, hinges, f"{filename}_hinges.hng")  

    l.info("### MOTION MOVIE ###")
    calpha=[i for i in Protein[0][chains[0].get_id()].get_calpha() if i is not None]
    Model=hdANM(calpha,dr=15,power=0,hng_file=f"{filename}_hinges.hng")
    Model.calculate_hessian(mass_type='residue')
    Model.calculate_decomposition()
    Model.get_eigenvalues()
    Model.get_eigenvectors()
    Model.calculate_movie(6,scale=2,n=40, ftype="pdb")


    print("### STRUCTURAL COMPLIANCE ###")


    resids = [j.get_id() for  j in Protein[0][chains[0].get_id()].get_residues() if j is not None]

    #Step 2.3
    ANM_MODEL = ANM( calpha, pf=True, dr=float('Inf'), power=3 )

    #Step 3
    ANM_MODEL.calculate_hessian()
    ANM_MODEL.calculate_decomposition()
    ANM_MODEL.calculate_stiffness_compliance()

    stiffness_map  = ANM_MODEL.get_stiffness_map()
    compliance_map = ANM_MODEL.get_compliance_map()

    b_factors          = [i.get_bfactor() for i in calpha]
    fluctuations       = ANM_MODEL.get_fluctuations()
    stiffness_profile  = ANM_MODEL.get_stiffness_profile()
    compliance_profile = ANM_MODEL.get_compliance_profile()


    # PLOTTING

    from matplotlib import pyplot as plt
    import numpy

    plt.plot(resids, b_factors/numpy.linalg.norm(b_factors), color= 'blue', alpha=0.4)
    plt.plot(resids, compliance_profile/numpy.linalg.norm(compliance_profile),color=  'black')
    # plt.plot(resids, stiffness_profile/numpy.linalg.norm(stiffness_profile),color='green')
    for hinge in hinges:
        resid = [x.get_id() for x in hinge.get_elements()]
        plt.axvspan(resid[0], resid[-1], color='green', alpha=0.4)
    for hinge in hinges_nosig:
        resid = [x.get_id() for x in hinge.get_elements()]
        plt.axvspan(resid[0], resid[-1], color='red', alpha=0.4)

    plt.show()



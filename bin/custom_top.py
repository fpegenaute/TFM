"""
Python file with:
    - RigidBody class: info for IMP topology files

    - generate_molec_dict function

	- Molecule dictionary (to modify if needed)

	- Python process to write custom topology file using data in molecule dict. 
"""

from pathlib import PurePosixPath
from tkinter.messagebox import NO
from bin.utilities import get_filename_ext
from Bio.PDB import MMCIFParser, PDBParser
from Bio import SeqIO
import itertools
import pandas as pd
import os
from pathlib import Path
import logging as l
import fnmatch
import copy

class RigidBody():
    """
    Object containing info for the IMP topology file. 
    
    Arguments: Bio PDB Structure Object. Note that this class is NOT a Child of
    Bio.PDB.Structure.
    
    Attributes:
     - resolution (int)
     - Chain_ID (str)
     - color (str)
     - fasta_file (str/path)
     - fasta_ID (str)
     - pdb_file (str/path)
     - chain_ID
     - residue_range(tuple)
     - pdb_offset (int)
     - bead_size (int)
     - em_residues_per_gaussian (int). def. 0 (leave if no EM fitting will be done)
     - rigid_body (int)
     - super_rigid_body (int)
     - chain_of_super_rigid_bodies (int) 
     - overlap (list)
     - type (str)
     - include (Bool)

    """
    def __init__(self, resolution, molecule_name, color, fasta_fn, 
        pdb_fn, chain, residue_range, 
        rigid_body, super_rigid_body, chain_of_super_rigid_bodies, type,  
        bead_size=10, em_residues_per_gaussian=0):

        self.resolution = resolution
        self.molecule_name = molecule_name
        self.color = color
        self.fasta_fn = fasta_fn
        self.fasta_id = PurePosixPath(fasta_fn).stem
        self.pdb_fn = pdb_fn
        self.chain = chain
        self.residue_range = residue_range
        
        if self.residue_range[0] == 1:
            self.pdb_offset = 0
        else:
            self.pdb_offset = - (self.residue_range[0] - 1)  
             
        self.bead_size = bead_size
        self.em_residues_per_gaussian = em_residues_per_gaussian
        self.rigid_body = rigid_body
        self.super_rigid_body = super_rigid_body
        self.chain_of_super_rigid_bodies = chain_of_super_rigid_bodies
        self.overlap = []
        self.type = type
        self.include = True
        self.attributes = [self.resolution, self.molecule_name, self.color, 
            self.fasta_fn, self.fasta_id, self.pdb_fn, self.chain,
            self.residue_range, self.pdb_offset, self.bead_size, 
            self.em_residues_per_gaussian, self.rigid_body,
            self.super_rigid_body, self.chain_of_super_rigid_bodies, self.overlap, 
            self.type, self.include]
    def get_full_PDB(self):
        """
        Get the name of the full pdb belonging to the chain of the RB
        and the parent folder. e.g. PDB/5yfp.pdb
        """
        filename = str(Path(*Path(self.pdb_fn).parts[-1:]))
        full_pdbs = []
       
        if "PDB" in str(self.pdb_fn):
            filename = filename[0:-6]
            path = os.path.abspath(Path(*Path(self.pdb_fn).parts[:-2]))
            for child in Path(path).iterdir():
                    if Path(child).is_file() and filename in str(child):
                        full_pdbs.append(Path(*Path(child).parts[-2:]))
            
        if "-crderr" in filename:
            filename = filename[0:-10]+filename[-4:]
            filename = filename.replace("-", ".")
            path = os.path.abspath(Path(*Path(self.pdb_fn).parts[:-2]))
            for child in Path(path).iterdir():
                    if Path(child).is_dir() and "DOMAINS" in str(child):
                        for rosetta_child in Path(child).iterdir():
                            if Path(rosetta_child).is_file() and filename in str(rosetta_child):
                                full_pdbs.append(Path(*Path(rosetta_child).parts[-3:]))
                                


        if "ALPHAFOLD" in str(self.pdb_fn):
            filename = filename[0:-10]+filename[-4:]
            path = os.path.abspath(Path(*Path(self.pdb_fn).parts[:-2]))
            for child in Path(path).iterdir():
                if Path(child).is_dir() and "DOMAINS" in str(child):
                    for alpha_child in Path(child).iterdir():
                        if "domains" in str(alpha_child):
                            full_pdbs.append(Path(*Path(alpha_child).parts[-3:]))
        #     
        return full_pdbs


    def get_structure(self):
        """
        Parse the structure file (self.pdb_fn) and return a structure object
        (Bio.PDB)
        """
        filename, ext = get_filename_ext(self.pdb_fn)

        if ext == "pdb" or ext == "ent" :
            parser = PDBParser(QUIET=True)
        elif ext == "cif":
            parser = MMCIFParser(QUIET=True)
        else:
            raise NameError("""Your file must have \"pdb\", \"ent\" or \"cif\" as 
            an extension""")
        structure = parser.get_structure(filename, self.pdb_fn) 
        return structure
    
    def get_resIDs(self):
        """
        Returns a list with the residue numbers of the rigid body
        """
        structure = self.get_structure()
        chainID = self.chain
        residues = []

        chain = structure[0][chainID]
        for residue in chain.get_residues():
            if residue.get_full_id()[3][0] == " ": # Exclude hetatm and h20
                residues.append(residue.id[1])
        
    
        return (residues)
    
    def update_overlap(self, rigid_body):
        """
        Given another RigidBody instance, check their structural overlap, and 
        update their overlap atttributes accordingly. Update also the residue range
        """
        overlap = list(set(self.get_resIDs()) & set(rigid_body.get_resIDs()))
        overlap.sort()
        
        
        if len(overlap) >= 1:
            # Prioritize the experimental structures
            if fnmatch.fnmatch(self.pdb_fn, "*[AR]F.pdb") and \
            ("AF.pdb" not in str(rigid_body.pdb_fn) and "RF.pdb" not in str(rigid_body.pdb_fn)):
                self.overlap = overlap
                self.residue_range = (self.overlap[0], self.overlap[-1])
                l.info(f"""Overlap attribute updated for {self.pdb_fn} 
                (predicted structure)""")
                return
            elif fnmatch.fnmatch(rigid_body.pdb_fn, "*[AR]F.pdb") and \
                 ("AF.pdb" not in str(self.pdb_fn) and "RF.pdb" not in str(self.pdb_fn)):
                rigid_body.overlap = overlap
                rigid_body.residue_range = \
                    (rigid_body.overlap[0], rigid_body.overlap[-1])
                l.info(f"""overlap attribute updated for {self.pdb_fn}  
                (predicted structure)""")
                return
            # If both are models, choose the biggest one
            elif fnmatch.fnmatch(self.pdb_fn, "*[AR]F.pdb") and \
                fnmatch.fnmatch(rigid_body.pdb_fn, "*[AR]F.pdb"):
                    if len(self.get_resIDs()) < len(rigid_body.get_resIDs()):
                        self.overlap = overlap
                        self.residue_range = \
                            (self.overlap[0], self.overlap[-1])
                        l.info(f"overlap attribute updated for {self.pdb_fn}")
                        return
                    elif len(self.get_resIDs()) > len(rigid_body.get_resIDs()):
                        rigid_body.overlap = overlap
                        rigid_body.residue_range = \
                            (rigid_body.overlap[0], rigid_body.overlap[-1])
                    l.info(f"""overlap attribute updated for 
                    {rigid_body.pdb_fn}""")
                    return
            # If both are experimental, choose the biggest one
            elif len(self.get_resIDs()) < len(rigid_body.get_resIDs()):
                self.overlap = overlap
                self.residue_range = (self.overlap[0], self.overlap[-1])
                l.info(f"overlap attribute updated for {self.pdb_fn}")
                return
            elif len(self.get_resIDs()) > len(rigid_body.get_resIDs()):
                rigid_body.overlap = overlap
                rigid_body.residue_range = \
                    (rigid_body.overlap[0], rigid_body.overlap[-1])
                l.info(f"overlap attribute updated for {rigid_body.pdb_fn}")
                return
            elif len(self.get_resIDs()) == len(rigid_body.get_resIDs()):
                if self.extract_avg_BFactor() < rigid_body.extract_avg_BFactor():
                    self.overlap = overlap
                    self.residue_range = (self.overlap[0], self.overlap[-1])
                    l.info(f"""Full overlap between {self.pdb_fn} and 
                    {rigid_body.pdb_fn} Assigning the full overlap to {self.pdb_fn} 
                    (It has lower average b factor). This means it will be 
                    discarded later, as overlap==length""")
                    return

                elif self.extract_avg_BFactor() > rigid_body.extract_avg_BFactor():
                    rigid_body.overlap = overlap
                    rigid_body.residue_range = (rigid_body.overlap[0], 
                                    rigid_body.overlap[-1])
                    l.info(f"""Full overlap between {self.pdb_fn} and 
                    {rigid_body.pdb_fn} Assigning the full overlap to {rigid_body.pdb_fn}  
                    (It has lower average b factor). This means it will be 
                    discarded later, as overlap==length""")
                    return
                
                else:
                    self.overlap = overlap
                    self.residue_range = (self.overlap[0], self.overlap[-1])
                    l.info(f"""Full overlap between {self.pdb_fn} and 
                    {rigid_body.pdb_fn} and equal average BFactors. Assigning arbitrarely
                    the full overlap to {self.pdb_fn}. This means it will be 
                    discarded later, as overlap==length""")
                    return 
        else:
            l.info(f"No overlap between {self.pdb_fn} and {rigid_body.pdb_fn}")
            return 
               

    def get_length(self):
        """
        Get the total length of the rigid body (no discotinuities taken into
        account)
        """
        return self.residue_range[1] - self.residue_range[0]
    
    def extract_avg_BFactor(self) :
        """
        Given a structure (Bio.PDB) return the averaged B factor
        """
        total_res_BFactor = float()

        for model in self.get_structure():
            for chain in model:
                j = 0
                for residue in chain:
                    if residue.get_full_id()[3][0] == " " :
                        sum_BFactors = 0
                        i = 0
                        for atom in residue:
                            sum_BFactors += atom.get_bfactor()
                            i += 1
                        # Add the residue Bfactor to the totel
                        total_res_BFactor += sum_BFactors/i
                    j += 1
        return total_res_BFactor/j

    def get_coverage(self, save_csv=False, outdir=None):
        """
        Return a pandas dataframe with the 
        coverage of the structure w.r.t it. 

        The df will have two columns 'ResID' (int) and 'filename' (0/1)
        save_csv: save the df as a csv
        outfile: Name of the file path to save the csv 
        """
        ref_ids, covered_ids = extract_coincident_positions(self.fasta_fn, 
                                                                self.pdb_fn)
        coverage_df = pd.DataFrame({"ResID":ref_ids,
                                    f"{self.pdb_fn}" :covered_ids})

        if save_csv:
            outdir = os.path.join(outdir, "COVERAGE")
            out_path = os.path.join(outdir, f"{self.structure_ID}_coverage.csv")
            try:
                coverage_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')
            except Exception:
                Path(outdir).mkdir(parents=True, exist_ok=True)
                coverage_df.to_csv(out_path, encoding='utf-8', 
                                            index=False, float_format='%.3f')

        return coverage_df

    def split_rb_hinges(self, hinges_list):
        """
        Given a RigidBody instance and a list of hinges, split that RigidBody.

        - self: rigidBody instance
        - hinges_list: list of hinges in tuple format, e.g. [ (2,5), (124,456) ]

        """
        rb_list = []

        # Sort them by appearence on the protein seq
        hinges_list.sort(key=lambda tup: tup[0])

        

        
        if len(hinges_list) == 0:
            rb_list.append(self)
            return rb_list    
        elif len(hinges_list) == 1:
            if hinges_list[0][0] < self.residue_range[0] and \
                    hinges_list[0][1] < self.residue_range[0]:
                    rb_list.append(self)
                    return rb_list 
            elif hinges_list[0][0] > self.residue_range[1] and \
                hinges_list[0][1] > self.residue_range[1]:
                    rb_list.append(self)
                    return rb_list 
            else:
                rb1 = copy.deepcopy(self)
                # hinge covering the beginning of the structure
                if hinges_list[0][0] < self.residue_range[0] and\
                    hinges_list[0][1] > self.residue_range[0] and \
                     hinges_list[0][1] < self.residue_range[1]:
                        rb1.residue_range = (hinges_list[0][1], self.residue_range[1])
                        rb_list = rb_list + [rb1]
                # hinge starting at the beginning of the structure
                if hinges_list[0][0] == self.residue_range[0] and\
                    hinges_list[0][1] < self.residue_range[1]:
                        rb1.residue_range = (hinges_list[0][1], self.residue_range[1])
                        rb_list = rb_list + [rb1]
                # hinge in the middle of the structure
                if hinges_list[0][0] > self.residue_range[0] and \
                     hinges_list[0][1] < self.residue_range[1]:
                        rb1.residue_range = (self.residue_range[0], hinges_list[0][0], )
                        rb2 = copy.deepcopy(self)
                        rb2.residue_range = (hinges_list[0][1], self.residue_range[1])
                        rb_list = rb_list + [rb1, rb2]
                # hinge ending at the end of the structure
                if hinges_list[0][1] == self.residue_range[1] and\
                    hinges_list[0][0] < self.residue_range[1] and \
                        hinges_list[0][0] > self.residue_range[0]:
                        rb1.residue_range = (self.residue_range[0], hinges_list[0][0])
                        rb_list = rb_list + [rb1]
                # hinge covering the end of the structure
                if hinges_list[0][0] > self.residue_range[0] and\
                    hinges_list[0][0] < self.residue_range[1] and\
                     hinges_list[0][1] > self.residue_range[1]:
                        rb1.residue_range = (self.residue_range[0],hinges_list[0][0])
                        rb_list = rb_list + [rb1]
                return rb_list
            
        elif len(hinges_list) > 1:
            for index in range(len(hinges_list)):                
                if index == 0:
                    print("index0")
                    if hinges_list[index][0] < self.residue_range[0] and \
                        hinges_list[index][1] < self.residue_range[0]:
                            rb_list.append(self)
                            continue
                    elif hinges_list[index][0] > self.residue_range[1] and \
                        hinges_list[index][1] > self.residue_range[1]:
                            rb_list.append(self)
                            continue
                    else:
                        rb1 = copy.deepcopy(self)
                        # hinge covering the beginning of the structure
                        if hinges_list[index][0] < self.residue_range[0] and\
                            hinges_list[index][1] > self.residue_range[0] and \
                            hinges_list[index][1] < self.residue_range[1]:
                                rb1.residue_range = (hinges_list[index][1], self.residue_range[1])
                                rb_list = rb_list + [rb1]
                        # hinge starting at the beginning of the structure
                        if hinges_list[index][0] == self.residue_range[0] and\
                            hinges_list[index][1] < self.residue_range[1]:
                                rb1.residue_range = (hinges_list[index][1], self.residue_range[1])
                                rb_list = rb_list + [rb1]
                        # hinge in the middle of the structure
                        if hinges_list[index][0] > self.residue_range[0] and \
                            hinges_list[index][1] < self.residue_range[1]:
                                rb1.residue_range = (self.residue_range[0], hinges_list[index][0], )
                                rb2 = copy.deepcopy(self)
                                rb2.residue_range = (hinges_list[index][1], self.residue_range[1])
                                rb_list = rb_list + [rb1, rb2]
                        # hinge ending at the end of the structure
                        if hinges_list[index][1] == self.residue_range[1] and\
                            hinges_list[index][0] < self.residue_range[1] and \
                                hinges_list[index][0] > self.residue_range[0]:
                                rb1.residue_range = (self.residue_range[0], hinges_list[index][0])
                                rb_list = rb_list + [rb1]
                        # hinge covering the end of the structure
                        if hinges_list[index][0] > self.residue_range[0] and\
                            hinges_list[index][0] < self.residue_range[1] and\
                            hinges_list[index][1] > self.residue_range[1]:
                                rb1.residue_range = (self.residue_range[0],hinges_list[index][0])
                                rb_list = rb_list + [rb1]
                    
                else:
                    rb1 = rb_list[-1]
                    if hinges_list[index][0] < rb1.residue_range[0] and \
                        hinges_list[index][1] < rb1.residue_range[0]:
                        rb_list.append(rb1)
                        continue
                    elif hinges_list[index][0] > rb1.residue_range[1] and \
                        hinges_list[index][1] > rb1.residue_range[1]:
                        rb_list.append(rb1)
                        continue
                    # hinge covering the beginning of the structure
                    if hinges_list[index][0] < rb1.residue_range[0] and\
                        hinges_list[index][1] > rb1.residue_range[0] and \
                        hinges_list[index][1] < rb1.residue_range[1]:
                            rb1.residue_range = (hinges_list[index][1], rb1.residue_range[1])
                            rb_list = rb_list + [rb1]
                            continue
                    # hinge starting at the beginning of the structure
                    if hinges_list[index][0] == rb1.residue_range[0] and\
                        hinges_list[index][1] < rb1.residue_range[1]:
                            rb1.residue_range = (hinges_list[index][1], rb1.residue_range[1])
                            rb_list = rb_list + [rb1]
                            continue
                    # hinge in the middle of the structure
                    # TO DO: apply the same correction of the rb_list[-1].residue_range to the rest of cases
                    if hinges_list[index][0] > rb1.residue_range[0] and \
                        hinges_list[index][1] < rb1.residue_range[1]:
                            rb2 = copy.deepcopy(rb1)
                            rb2.residue_range = (hinges_list[index][1], rb1.residue_range[1])
                            rb_list[-1].residue_range = (rb1.residue_range[0], hinges_list[index][0])
                            rb_list = rb_list + [rb2]
                            continue
                    # hinge ending at the end of the structure
                    if hinges_list[index][1] == rb1.residue_range[1] and\
                        hinges_list[index][0] < rb1.residue_range[1] and \
                            hinges_list[index][0] > rb1.residue_range[0]:
                            rb1.residue_range = (rb1.residue_range[0], hinges_list[index][0])
                            rb_list = rb_list + [rb1]
                            continue
                    # hinge covering the end of the structure
                    if hinges_list[index][0] > rb1.residue_range[0] and\
                        hinges_list[index][0] < rb1.residue_range[1] and\
                        hinges_list[index][1] > rb1.residue_range[1]:
                            rb1.residue_range = (rb1.residue_range[0],hinges_list[index][0])
                            rb_list = rb_list + [rb1]
                            continue



            return rb_list    
        


### FUNCTIONS
from bin.graphical_summary import PDB_get_resid_set, extract_coincident_positions 

def make_composite(rb_list, reference_fasta=None, save_csv=False, outdir=None):
    """
    Given a RigidBody list, retunr a RigidBody list containing the same RigidBodies, 
    but with adjusted residue range, so they achieve the highest coverage possible
    of the reference query fasta file avoiding overlaps.

    Args:
        - rb_list: List of RigidBody objects
        - reference_fasta: Reference sequence
        - save_csv: Save a csv with the columns: "Reference", "RigidBody1","RigidBodyN"  
             - Reference: Reference resIDs
             - RigidBodyN: ResID present in the RigidBodies (0 absent/1 present)
        - outdir: Directory to store the composite .csv file

    """
    print(f"ORIGINAL RB LIST LEN: {len(rb_list)}")

    # calculate the overlaps, return the fragments with no overlaps


    for pair in itertools.combinations(rb_list, 2):
        rb1, rb2 = pair     
        rb1.update_overlap(rb2)
        print(f"OVERLAP {rb1.pdb_fn}: {len(rb1.overlap)}, OL {rb2.pdb_fn}: {len(rb2.overlap)}")
        
    clean_rb_list = [rb for rb in rb_list if len(rb.overlap) == 0]

    print(f"FINAL RB LIST LEN: {len(clean_rb_list)}")  

    return clean_rb_list

from bin.utilities import get_chain_names, get_residue_range

def make_rb_list(structures_list, fasta):
    """
    Given a list of paths of PDB or mmcif files, return a list of all of them as 
    RigidBody Objects. Classify them as AlphaFold models, RoseTTaFold models, or 
    experimental structures.

    Include as an argument the reference fasta file that all the structures 
    correspond to
    """
    i = 1
    rigid_bodies = []
    for structure in structures_list:    
        filename, extension = get_filename_ext(structure)
        
        # Extract chain name
        chain_IDs = get_chain_names(structure)
        if len(chain_IDs) > 1 or fnmatch.fnmatch(structure, "*AF.pdb"):
            l.info(f"""{structure}: Assuming a AlphaFold model""")

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
                bead_size=4,
                em_residues_per_gaussian=0, 
                type="AF_model")
                # Add the rigid body to a list
                rigid_bodies.append(rigid_body)
                i +=1
        elif fnmatch.fnmatch(structure, "*RF.pdb"):
            l.info(f"""{structure}: Assuming a RoseTTaFold model""")

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
                bead_size=4,
                em_residues_per_gaussian=0, 
                type="RF_model")
                # Add the rigid body to a list
                rigid_bodies.append(rigid_body)
                i +=1
        elif len(chain_IDs) > 1:
            l.info(f"Assuming {structure},experimental with multiple chains")
            for chain in chain_IDs:
                # Extract residue range
                res_range = get_residue_range(structure, chain=chain)    
                # Create the RigidBody instance
                rigid_body = RigidBody(resolution="all",
                molecule_name= f"{filename}_{chain}", 
                color="blue" , 
                fasta_fn=fasta, 
                # fasta_id=fasta, 
                pdb_fn=structure, 
                chain=chain,
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
        else:
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

    return rigid_bodies

def write_custom_topology(path_to_file, rigid_body_list):
    """
    Method to write a custom topology.txt file for IMP.pmi modeling
    :param path_to_file: path to write custom topology
    :rigid_body_list: list of RigidBody objects
    """

    top_file = open(path_to_file, "w")
    # Write header of topology
    header = ["molecule_name", "color", "fasta_fn", "fasta_id", "pdb_fn", 
    "chain", "residue_range", "pdb_offset", "bead_size", 
    "em_residues_per_gaussian", "rigid_body", "super_rigid_body",
              "chain_of_super_rigid_bodies"]

    
    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                   "{:28}|{:12}|{:19}|{:27}|\n".format(header[0], header[1], 
                    header[2], header[3], header[4], header[5], header[6], 
                    header[7], header[8], header[9], header[10], header[11], 
                    header[12]))
    # Write molecule lines from dictionary
    # Notice that you can modify the following line acording to the desired output
    rigid_body_counter = 1
    first = True
    for rb in  rigid_body_list:
        resolution = rb.resolution
        mol_name = rb.fasta_id
        color = rb.color
        fasta_fn = PurePosixPath(rb.fasta_fn).name
        fasta_id = rb.fasta_id #+":"+rb.chain
        pdb_fn = str(rb.get_full_PDB()[0])
        chain = rb.chain
        start_residue = rb.residue_range[0]
        last_residue = rb.residue_range[1]
        # offset = rb.pdb_offset
        offset = "" # The offset is not necessary if the res range is correct
        bead_size = rb.bead_size
        em_gaussian = rb.em_residues_per_gaussian
        print(f"RB COUNT {rigid_body_counter}")
        if rigid_body_counter == 1:
            if start_residue == 1:
                top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                            "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, 
                            fasta_fn, fasta_id, str(pdb_fn), chain,"{},{}".format(start_residue, 
                                last_residue), offset, bead_size, em_gaussian, 
                                rigid_body_counter, "", ""))
                top_file.write("\n") 
                rigid_body_counter += 1
            elif start_residue > 1:
                top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                            "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, "brown", 
                            fasta_fn, fasta_id, "BEADS", chain,"{},{}".format(1, 
                                start_residue-1), offset, bead_size, em_gaussian, 
                                rigid_body_counter, "", ""))
                top_file.write("\n") 
                top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                            "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, 
                            fasta_fn, fasta_id, str(pdb_fn), chain,"{},{}".format(start_residue, 
                                last_residue), offset, bead_size, em_gaussian, 
                                rigid_body_counter, "", ""))
                top_file.write("\n") 
                rigid_body_counter += 1
            continue
        if rigid_body_counter > 1 :
            print(f"STARR{start_residue} RB  end + 2{rigid_body_list[rigid_body_counter-2].residue_range[1]+1}")
            if rigid_body_counter < len(rigid_body_list):
                print("A")
                if start_residue == rigid_body_list[rigid_body_counter-2].residue_range[1]+1:
                    print("A1")
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, 
                                fasta_fn, fasta_id, str(pdb_fn), chain,"{},{}".format(start_residue, 
                                    last_residue), offset, bead_size, em_gaussian, 
                                    rigid_body_counter, "", ""))
                    top_file.write("\n") 
                    rigid_body_counter += 1
                elif start_residue > rigid_body_list[rigid_body_counter-2].residue_range[1]+1:
                    print(f"A2")
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, "brown", 
                                fasta_fn, fasta_id, "BEADS", chain,"{},{}".format(rigid_body_list[rigid_body_counter-2].residue_range[1]+1, 
                                    start_residue-1), offset, bead_size, em_gaussian, 
                                    rigid_body_counter, "", ""))
                    top_file.write("\n") 
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, 
                                fasta_fn, fasta_id, str(pdb_fn), chain,"{},{}".format(start_residue, 
                                    last_residue), offset, bead_size, em_gaussian, 
                                    rigid_body_counter, "", ""))
                    top_file.write("\n") 
                    rigid_body_counter += 1
                continue
            if rigid_body_counter == len(rigid_body_list):
                print("B")
                records = list(SeqIO.parse(rb.fasta_fn, "fasta"))
                sequence = str(records[0].seq)
               
                if start_residue == rigid_body_list[rigid_body_counter-2].residue_range[1]+1:
                    print("B1")
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, 
                                fasta_fn, fasta_id, str(pdb_fn), chain,"{},{}".format(start_residue, 
                                    last_residue), offset, bead_size, em_gaussian, 
                                    rigid_body_counter, "", ""))
                    top_file.write("\n") 
                    rigid_body_counter += 1
                elif start_residue > rigid_body_list[rigid_body_counter-2].residue_range[1]+1:
                    print("B2")
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, "brown", 
                                fasta_fn, fasta_id, "BEADS", chain,"{},{}".format(rigid_body_list[rigid_body_counter-2].residue_range[1]+1, 
                                    start_residue-1), offset, bead_size, em_gaussian, 
                                    rigid_body_counter, "", ""))
                    top_file.write("\n") 
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, 
                                fasta_fn, fasta_id, str(pdb_fn), chain,"{},{}".format(start_residue, 
                                    last_residue), offset, bead_size, em_gaussian, 
                                    rigid_body_counter, "", ""))
                    top_file.write("\n") 
                    rigid_body_counter += 1
                if last_residue == len(sequence):
                    print("B3")
                    continue
                if last_residue < len(sequence):
                    print("B4")
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, "brown", 
                                fasta_fn, fasta_id, "BEADS", chain,"{},{}".format(last_residue+1, 
                                    len(sequence)), offset, bead_size, em_gaussian, 
                                    rigid_body_counter, "", ""))
                    top_file.write("\n") 
                    rigid_body_counter += 1
            continue

            top_file.write("\n") 
      
    top_file.close()
      




if __name__ == '__main__':

    

    rb1 = RigidBody("all", "DNA_A", "red", "complex.fasta", "DNA1,DNA", 
        "complex.pdb", "a", (1, 250), "0", "1", "0")
    rb2 = RigidBody("all", "DNA_B", "blue", "complex.fasta", "DNA2,DNA", 
        "complex.pdb", "b", (251, 500), "0", "1", "0")

    rb_list = [rb1, rb2]

    write_custom_topology("test_topology", rb_list)
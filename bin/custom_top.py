"""
Python file with:
    - RigidBody class: info for IMP topology files

    - generate_molec_dict function

	- Molecule dictionary (to modify if needed)

	- Python process to write custom topology file using data in molecule dict. 
"""

from pathlib import PurePosixPath
from bin.utilities import get_filename_ext
from Bio.PDB import MMCIFParser, PDBParser
import itertools
import pandas as pd
import os
from pathlib import Path

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
        self.attributes = [self.resolution, self.molecule_name, self.color, 
            self.fasta_fn, self.fasta_id, self.pdb_fn, self.chain,
            self.residue_range, self.pdb_offset, self.bead_size, 
            self.em_residues_per_gaussian, self.rigid_body,
            self.super_rigid_body, self.chain_of_super_rigid_bodies, self.overlap]
    
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
            if len(self.get_resIDs()) < len(rigid_body.get_resIDs()):
                self.overlap = overlap
                self.residue_range = (self.overlap[0], self.overlap[-1])
                print(f"overlap attribute updated for {self.pdb_fn}")
            elif len(self.get_resIDs()) > len(rigid_body.get_resIDs()):
                rigid_body.overlap = overlap
                rigid_body.residue_range = (rigid_body.overlap[0], rigid_body.overlap[-1])
                print(f"overlap attribute updated for {rigid_body.pdb_fn}")
            else:
                if self.extract_avg_BFactor() < rigid_body.extract_avg_BFactor():
                    self.overlap = overlap
                    self.residue_range = (self.overlap[0], self.overlap[-1])
                    print(f"""Full overlap between {self.pdb_fn} and 
                {rigid_body.pdb_fn} Assigning the full overlap to {self.pdb_fn} 
                (It has lower average b factor). This means it will be 
                discarded later, as overlap==length""")

                elif self.extract_avg_BFactor() > rigid_body.extract_avg_BFactor():
                    rigid_body.overlap = overlap
                    rigid_body.residue_range = (rigid_body.overlap[0], 
                                    rigid_body.overlap[-1])
                    print(f"""Full overlap between {self.pdb_fn} and 
                {rigid_body.pdb_fn} Assigning the full overlap to {rigid_body.pdb_fn}  
                (It has lower average b factor). This means it will be 
                discarded later, as overlap==length""")
                
                else:
                    self.overlap = overlap
                    self.residue_range = (self.overlap[0], self.overlap[-1])
                    print(f"""Full overlap between {self.pdb_fn} and 
                {rigid_body.pdb_fn} and equal average BFactors. Assigning arbitrarely
                the full overlap to {self.pdb_fn}. This means it will be 
                discarded later, as overlap==length""")
        else:
            print(f"No overlap between {self.pdb_fn} and {rigid_body.pdb_fn}")
               

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

        The df will have two columns 'ResID' (int) and 'Structure' (0/1)
        save_csv: save the df as a csv
        outfile: Name of the file path to save the csv 
        """
        ref_ids, covered_ids = extract_coincident_positions(self.fasta_fn, 
                                                                self.pdb_fn)
        coverage_df = pd.DataFrame({"ResID":ref_ids,
                                    f"Structure{self.pdb_fn}" :covered_ids})

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

        

        


### FUNCTIONS
from bin.graphical_summary import extract_coincident_positions

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
    
    for pair in itertools.combinations(rb_list, 2):
        rb1, rb2 = pair
        rb1.update_overlap(rb2)
    
    # Discard RigidBodies fully overlapped
    for rb in rb_list:
        if len(rb.get_resIDs()) == len(rb.overlap):
            rb_list.remove(rb)
    
    # Discard RBs whose overlap is fully covered by another structure
    for pair in itertools.combinations(rb_list, 2):
        rb1, rb2 = pair
        if len(set(rb2.overlap) & set(rb1.get_resIDs())) == len(set(rb2.overlap)):
            rb_list.remove(rb2)
            continue
        if len(set(rb1.overlap) & set(rb2.get_resIDs())) == len(set(rb1.overlap)):
            rb_list.remove(rb1)
            



    return rb_list

def write_custom_topology(path_to_file, rigid_body_list):
    """
    Method to write a custom topology.txt file for IMP.pmi modeling
    :param path_to_file: path to write custom topology
    :rigid_body_list: list of RigidBody objects
    """
    top_file = open(path_to_file, "w")
    # Write header of topology
    header = ["molecule_name", "color", "fasta_fn", "fasta_id", "pdb_fn", "chain", "residue_range",
              "pdb_offset", "bead_size", "em_residues_per_gaussian", "rigid_body", "super_rigid_body",
              "chain_of_super_rigid_bodies"]
    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                   "{:28}|{:12}|{:19}|{:27}|\n".format(header[0], header[1], header[2], header[3], header[4],
                                                       header[5], header[6], header[7], header[8], header[9],
                                                       header[10], header[11], header[12]))
    # Write molecule lines from dictionary
    # Notice that you can modify the following line acording to the desired output
    rigid_body_counter = 1
    for rb in rigid_body_list:
        resolution = rb.resolution
        mol_name = rb.molecule_name
        color = rb.color
        fasta_fn = rb.fasta_fn
        fasta_id = rb.fasta_id
        pdb_fn = rb.pdb_fn
        chain = rb.chain
        start_residue = rb.residue_range[0]
        last_residue = rb.residue_range[1]
        offset = rb.pdb_offset
        bead_size = rb.bead_size
        em_gaussian = rb.em_residues_per_gaussian
        if resolution == "all":
            top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                           "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, fasta_fn, fasta_id, pdb_fn,
                                                                chain, "all", offset, bead_size, em_gaussian,
                                                                rigid_body_counter, "", ""))
            rigid_body_counter += 1
        else:
            # If you want to add different options for X molecules (DNA or Proteins) this
            # is a way
            for n in range(start_residue, last_residue + 1, resolution):
                if start_residue + resolution <= last_residue:
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                   "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, fasta_fn, fasta_id, pdb_fn,
                                                                        chain, "{},{}".format(start_residue,
                                                                                              start_residue + resolution),
                                                                        offset, bead_size, em_gaussian,
                                                                        rigid_body_counter, "", ""))
                    start_residue += resolution + 1
                    rigid_body_counter += 1
                else:
                    top_file.write("|{:15}|{:10}|{:20}|{:15}|{:16}|{:7}|{:15}|{:11}|{:11}|"
                                "{:28}|{:<12}|{:19}|{:27}|\n".format(mol_name, color, fasta_fn, fasta_id, pdb_fn,
                                                                    chain, "{},{}".format(start_residue,
                                                                                            last_residue),
                                                                    offset, bead_size, em_gaussian,
                                                                    rigid_body_counter, "", ""))
                    rigid_body_counter += 1
                    break
        top_file.write("\n")  # write a blank line between different molecules
    top_file.close()

if __name__ == '__main__':

    

    rb1 = RigidBody("all", "DNA_A", "red", "complex.fasta", "DNA1,DNA", "complex.pdb", "a", (1, 250), "0", "1", "0")
    rb2 = RigidBody("all", "DNA_B", "blue", "complex.fasta", "DNA2,DNA", "complex.pdb", "b", (251, 500), "0", "1", "0")

    rb_list = [rb1, rb2]

    write_custom_topology("test_topology", rb_list)
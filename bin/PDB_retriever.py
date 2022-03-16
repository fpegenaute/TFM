from posixpath import join
from Bio.PDB import MMCIFParser, PDBParser, PDBList, PDBIO, Select
import requests
import os
import logging as l
from bin.utilities import choose_parser
  
def retrieve_pdb_info(hit_dict, pdb_dir, fasta_dir):
    """
    Given a list of PDB codes, it downloads the PDB/MMcif files, and gets the
    data from the pdb entry, from where it extracts the sequence and writes it 
    into a fasta file.

    Inputs:

    - hit_dict: List of PDB codes
    - pdb_dir: Directory to store the PDB files
    - fasta_dir: Directory to store the Fasta files

    ** note tha the fasta file contains all the sequence, not only the solved
    portion of the protein.

    It returns the request
    """

    # Initialize object PDBList
    pdbl = PDBList() 
    pdb_list = hit_dict.keys()

    #input(f"PDB List = {len(pdb_list)} Press Enter to continue...")
    for template in pdb_list:
        # Dowload the PDB file from the web
        pdbl.retrieve_pdb_file(template, pdir=pdb_dir)     
        req = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{template}').json()[template.lower()]
        
        fasta_template_name = os.path.join(fasta_dir, f"{template}.fa")
        with open(fasta_template_name, "w") as f:
            # Print chain (atm A)
            header = [template,req[0]["in_chains"][0]]
            headerstr = "_".join(header)
            sequence = req[0]['sequence']
            f.write(">"+headerstr+"\n")
            f.write(sequence)

    return




def check_PDB_len(pdbfile, chain):
    """
    Given a PDB/mmCIF file, check the actual length of the desired chain, 
    return the value
    """
    parser, identifier = choose_parser(pdbfile)
    # Get the structure
    structure = parser.get_structure(identifier, pdbfile) 
    model = structure[0]
    
    res_number = 0
    non_resi = 0

    for r in model[chain].get_residues():
        if r.id[0] == ' ':
            res_number +=1
        else:
            non_resi +=1

    return res_number




class ChainSplitter:
    def __init__(self, mmcif=True, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        if mmcif ==True:
            self.parser = MMCIFParser()
            self.writer = PDBIO()

        if mmcif ==False:
            self.parser = PDBParser()
            self.writer = PDBIO()
        
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, chain_letters, overwrite=False, struct=None):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """
        chain_letters = [chain.upper() for chain in chain_letters]

        # Input/output files
        (pdb_dir, pdb_filename) = os.path.split(pdb_path)
        pdb_id = os.path.basename(pdb_filename).split('.')[0]
        out_name = "%s_%s.pdb" % (pdb_id, "".join(chain_letters))
        out_path = os.path.join(self.out_dir,  out_name)
        l.info(f"OUT PATH:{out_path}")
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            l.info("Chain%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(chain_letters), pdb_id, out_name))
            return out_path

        l.info("Extracting chain%s %s from %s..." % (plural,
                ", ".join(chain_letters), pdb_filename))

        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters))

        return out_path


class SelectChains(Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)




if __name__ == "__main__":
    
    # Example
    pdb_list = ["1pa2","5tnt"]
    pdb_dir = "./templates/PDB/"
    fasta_dir="./templates/FASTA/"
    req = retrieve_pdb_info(pdb_list, pdb_dir, fasta_dir)



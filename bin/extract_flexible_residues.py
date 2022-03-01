from Bio.PDB.Structure import Structure
from bin.process_predicted_model import *

from Bio.PDB import MMCIFParser, PDBParser

def get_resid_list(structure_file, outdir):
    """
    Given a PDB or MMCif File, return a list of tuples of (ResID, Chain), and
    write a file of tab separated Resid Chain columns.
    """
    identifier, extension = get_filename_ext(structure_file)

    if extension == "pdb":
        parser = PDBParser(QUIET=True)
    elif extension == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise NameError("Your file has to have \"pdb\" or \"cif\" as an extension")

    structure = parser.get_structure(identifier, structure_file)    

    model = structure[0]
    resid_list = []
    with open(f"{outdir}/resid_chain.tsv") as resid_file:

        for chain in model:
            for i in chain.get_residues():
                # Print to file
                # The chain is the first element of the tuple in the 3rd position of the element
                resid_file.write(f"{i.get_full_id()[2]}\t{i.get_full_id()[3][1]}")
                # Add to a list of tuples
                resid_tuple = (i.get_full_id()[2],i.get_full_id()[3][1])
                resid_list.append(resid_tuple)

    return resid_list            


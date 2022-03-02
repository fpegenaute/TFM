"""
Python file with:
    - RigidBody class: info for IMP topology files

    - generate_molec_dict function

	- Molecule dictionary (to modify if needed)

	- Python process to write custom topology file using data in molecule dict. 
"""

class RigidBody():
    """
    Object containing info for the IMP topology file. Attributes:
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

    """
    def __init__(self, resolution, molecule_name, color, fasta_fn, fasta_id, 
        pdb_fn, chain, residue_range, 
        rigid_body, super_rigid_body, chain_of_super_rigid_bodies,  pdb_offset=0, 
        bead_size=10, em_residues_per_gaussian=0):

        self.attributes = [resolution, molecule_name, color, fasta_fn, fasta_id, pdb_fn, chain, 
        residue_range, pdb_offset, bead_size, em_residues_per_gaussian, rigid_body,
        super_rigid_body, chain_of_super_rigid_bodies]
        self.resolution = resolution
        self.molecule_name = molecule_name
        self.color = color
        self.fasta_fn = fasta_fn
        self.fasta_id = fasta_id
        self.pdb_fn = pdb_fn
        self.chain = chain
        self.residue_range = residue_range
        if self.residue_range[0] == 1:
            self.pdb_offset = pdb_offset
        else:
            self.pdb_offset = - (self.residue_range[0] - 1)   
        self.bead_size = bead_size
        self.em_residues_per_gaussian = em_residues_per_gaussian
        self.rigid_body = rigid_body
        self.super_rigid_body = super_rigid_body
        self.chain_of_super_rigid_bodies = chain_of_super_rigid_bodies



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
                    break
        top_file.write("\n")  # write a blank line between different molecules
    top_file.close()

if __name__ == '__main__':

    

    rb1 = RigidBody("all", "DNA_A", "red", "complex.fasta", "DNA1,DNA", "complex.pdb", "a", (1, 250), "0", "1", "0")
    rb2 = RigidBody("all", "DNA_B", "blue", "complex.fasta", "DNA2,DNA", "complex.pdb", "b", (251, 500), "0", "1", "0")

    rb_list = [rb1, rb2]

    write_custom_topology("test_topology", rb_list)
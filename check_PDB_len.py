from Bio import PDB


# Check length of the PDBs
parser = PDB.PDBParser()

# for pdb in pdb_dir:

pdb1 ='./1bfg.pdb' 
structure = parser.get_structure("1bfg", "pdb") 
model = structure[0]
res_no = 0
non_resi = 0

for model in structure:
    for chain in model:
        for r in chain.get_residues():
            if r.id[0] == ' ':
                res_no +=1
            else:
                non_resi +=1

print ("Residues:  %i" % (res_no))
print ("Other:     %i" % (non_resi))
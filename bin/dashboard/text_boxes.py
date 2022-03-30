# Text boxes
intro_md = '''
# Report

This is a summary of the results 

## Coverage

In this grapph,you can see whic parts of the reference FASTA sequence are 
covered by structures in the PDB
'''

flex_md = """
## Hinges and flexibility

Flexibility is an important feature of proteins, since they need to move to 
perform their function and interact with their substrates. In the following 
section, we provide you with two types of flexibility prediction: the Dynamic 
Flexibility Index and Hinge Prediction.

*Dynamic Flexibility Index*  
This is per-residue index indicating the contribution of each residue to the 
overall flexibility of the protein. It uses a method based in an Elastic Network 
Model, which is a more lightweight (but less precise, obviously) alternative to 
Molecular Dynamics. for ore info, 
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673471/) is the original 
paper.

*Hinge Prediction*  
Hinges are the regions of the protein that allow it to move and change 
conformations. Using 
[this tool](https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbac007/6525212?login=true) 
the predicted hinge regions are showed on top of the DFI plot, with the 
significative ones colored in green, and  the non-significative ones in red.

"""

composite_md = """
## Composite

From all the structures retrieved by the program and provided by the user, 
make a composite with all of them that covers as much of the reference sequence 
as possible (avoiding overlaps).

This comopsite will be used to automatically build an 
[IMP topology file](https://integrativemodeling.org/2.5.0/doc/ref/classIMP_1_1pmi_1_1topology_1_1TopologyReader.html)
"""
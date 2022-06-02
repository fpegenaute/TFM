# Text boxes
intro_md = '''
# Report

This is a summary of the results. It works for all the output directories inside
the output/ directory.

Select which of your proteins of interest  want to explore in the dropdown 
below

## Coverage

In this graph,you can see whic parts of the reference FASTA sequence are 
covered by structure. This structures come from either the [Protein Data Bank]
(https://www.rcsb.org/), [AlphaFold](https://www.deepmind.com/blog/alphafold-a-solution-to-a-50-year-old-grand-challenge-in-biology)
 models or [RoseTTaFold](https://www.ipd.uw.edu/2021/07/rosettafold-accurate-protein-structure-prediction-accessible-to-all/)
 models. 

'''

flex_md = """
## Hinges and flexibility

Flexibility is an important feature of proteins, since they need to move to 
perform their function and interact with their substrates. In the following 
section, we provide you with two types of flexibility prediction: the Dynamic 
Flexibility Index and Hinge Prediction. The overlap of these two measures
might be helpful for you, in case you wanted to modify the final topology file
with some of these hinges.

*Dynamic Flexibility Index*  
This is per-residue index indicating the contribution of each residue to the 
overall flexibility of the protein. It uses a method based in an Elastic Network 
Model (ENM), which is a more lightweight (but less precise, obviously) 
alternative to Molecular Dynamics. for more info, 
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673471/) is the original 
paper.

*Hinge Prediction*  
Hinges are the regions of the protein that allow it to move and change 
conformations. Using 
[this tool](https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbac007/6525212?login=true) 
 we provide you with some suggested hinge regions. Note that this information
 is only available for experimental structures. This is due to the use of ENM, 
 it is not designed to work with predicted models that might contain important
 artifacts, and, in this case, that are split into the highly confidently
 predicted regions.

"""

composite_md = """
## Composite and Topology File

From all the structures retrieved by the program and provided by the user, the
program generates this composite, trying to  cover as much of the reference 
sequence as possible, avoiding overlaps.

This comopsite can be used to automatically build an 
[IMP topology file](https://integrativemodeling.org/2.5.0/doc/ref/classIMP_1_1pmi_1_1topology_1_1TopologyReader.html)

If you want to generate another topology file, with your custom choice of 
fragments, simply select the fragments woy want to include and click the button "Create 
Topology File". it will save yhe output inj the folder IMP of your selected 
output folder.
"""

hinges_md = """
### Custom hinges
In this section you can introduce hinge regions. Hinge regions are those regions
of the protein that bend, allowing the movement of the more rigid domains, which 
is essential for the interaciton of the proteins with other biomolecules.

**How are hinges encoded?**  
Let's imagine you have a proterin of 200 amino acids. The DFI and PACKMAN hinge 
prediction are indicating a putative hinge region between positions 50 and 100 
and another  one between 120 and 130.


In the box, you will need to introduce the hinges with the following format:
50:100,120:130

The program will split the structures, in the topology file, according to the 
hinges introduced.
"""

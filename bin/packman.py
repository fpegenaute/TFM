
from packman.apps import predict_hinge
from packman import molecule
import numpy as np
from pathlib import Path


def write_hng_file(pdbfile, hinges, outfile):
    """
    Generate a hinge file
    """
    filename = Path(pdbfile).stem
    ALL_RESIDUES = {}
    Protein = molecule.load_structure(pdbfile)
    for i in Protein[0].get_chains():
        try:
            ALL_RESIDUES[i.get_id()] = sorted([i.get_id() for i in Protein[0][i.get_id()].get_residues() if i!=None])
        except:
            None

    select_count = 0
    last_hinge_end = 0
    fh = open(outfile, 'w')
    for numi, i in enumerate(hinges):
        current_hinge = hinges[numi]
        ChainOfHinge = current_hinge.get_elements()[0].get_parent().get_id()

        if select_count==0:
            hinge_res_ids = sorted([j.get_id() for j in current_hinge.get_elements()])
            
            select_count += 1
            if(ALL_RESIDUES[ChainOfHinge][0]!=hinge_res_ids[0]):
                fh.write(filename+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(ALL_RESIDUES[ChainOfHinge][0])+':'+str(hinge_res_ids[0]-1)+'\n' )
                fh.write(filename+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
            else:
                fh.write(filename+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
            last_hinge_end = hinge_res_ids[-1]
        else:
            hinge_res_ids = sorted([j.get_id() for j in current_hinge.get_elements()])
            select_count += 1
            fh.write(filename+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end+1)+':'+str(hinge_res_ids[0]-1)+'\n' )
            fh.write(filename+'_'+ChainOfHinge +'\t'+ 'H'+str(select_count) +'\t'+ str(hinge_res_ids[0])+':'+str(hinge_res_ids[-1])+'\n' )
            last_hinge_end = hinge_res_ids[-1]
        try:
            if(ChainOfHinge != hinges[numi+1].get_elements()[0].get_parent().get_id() ):
                select_count += 1
                fh.write(filename+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end+1)+':'+str(ALL_RESIDUES[ChainOfHinge][-1])+'\n' )
                last_hinge_end = 0
        except:
            None

        if(last_hinge_end != ALL_RESIDUES[ChainOfHinge][-1]):
            select_count += 1
            fh.write(filename+'_'+ChainOfHinge +'\t'+ 'D'+str(select_count) +'\t'+ str(last_hinge_end+1)+':'+str(ALL_RESIDUES[ChainOfHinge][-1])+'\n' )
    fh.flush()
    fh.close()

# def softmax(x):
#     """Compute softmax values for a list of values in x."""
#     e_x = np.exp(x - np.max(x))
#     return e_x / e_x.sum()

if __name__ == "__main__":

    pdbfile = input("Input Protein: ")
    Protein = molecule.load_structure(pdbfile)

    try:
        Protein[0]
    except Exception:
        print("Make sure your filename is  of the form: XXXXX.pdb/XXXX.cif")

    # In the tutorial they also use this, it gives similar results, try using backbone or Catoms
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

    predict_hinge(backbone, Alpha=3.65, outputfile=open(str("hinge_output")+'.txt', 'w'))
    
    hinges = []
    hinges_nosig = []
    print("Significant hinges")
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
    
    write_hng_file(pdbfile, hinges, "hinges.hng")  

    print("### MOTION MOVIE ###")
    from packman.anm import hdANM
    calpha=[i for i in Protein[0][chains[0].get_id()].get_calpha() if i is not None]
    Model=hdANM(calpha,dr=15,power=0,hng_file="hinges.hng")
    Model.calculate_hessian(mass_type='residue')
    Model.calculate_decomposition()
    Model.get_eigenvalues()
    Model.get_eigenvectors()
    Model.calculate_movie(6,scale=2,n=40, ftype="pdb")


    print("### STRUCTURAL COMPLIANCE ###")
    #Step 1.1
    from packman import molecule

    #Step 1.2
    #(Default is CIF format; to change it to PDB)
    #molecule.download_structure('1LF7',ftype = 'pdb')

    #Step 1.3
    # mol=molecule.load_structure(input("Load a PDB/mmCif file: "))

    mol=Protein
    mol=molecule.load_structure("4hu2.pdb")


    #Step 2.1
    from packman import anm

    #Step 2.2
    # c_alpha = mol[0].get_calpha()
    resids = [j.get_id() for  j in Protein[0][chains[0].get_id()].get_residues() if j is not None]

    

    #Step 2.3
    ANM_MODEL = anm.ANM( calpha, pf=True, dr=float('Inf'), power=3 )

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

    # If you want to check the correlation btw Bfactors and compliance
    # from scipy.stats import pearsonr
    # pearsonr(b_factors,compliance_profile)
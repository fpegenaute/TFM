
import logging as l

try:
    from modeller import *
    from modeller.scripts import complete_pdb
    from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions
    modeller_ok = True
except Exception:
    modeller_ok = False
    l.warn("Something is wrong with Modeller. The program will not have optimiaztion available. Check that it is installed and in your PYHONPATH.  If the problem is persistent, go to the Modeller website for more info")

from operator import itemgetter

def Optimizemodel(pdb_file):
	"""
	It creates a PDB file with the optimized model from the input pdb, 
	with its energies and restraint contributions. Also it will create pdbs on 
	every step of the Molecular Dynamics optimization. The energy is returned as
	the total value	of Modeller's objective function, molpdf. 
	It also shows the topt 10 contributors to the molpdf before and after optimization
	"""
	# Setting up
	env = environ()
	env.io.atom_files_directory = ['../atom_files']
	env.edat.dynamic_sphere = True

	env.libs.topology.read(file='$(LIB)/top_heav.lib')
	env.libs.parameters.read(file='$(LIB)/par.lib')

	code, ext = pdb_file.split('.')

	# Complete the pdb and make a model
	mdl = complete_pdb(env, pdb_file)
	mdl.write(file=code+'.ini')

	# Select all atoms from the model, make restraints and save them in a file
	atmsel = selection(mdl)
	mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
	mdl.restraints.write(file=code+'.rsr')

	# Check the energy before optimization and save it in a var
	initial_mpdf = atmsel.energy()


	# Create optimizer objects
	cg = conjugate_gradients(output='REPORT')
	md = molecular_dynamics(output='REPORT')

	# Open a file to get basic stats on each optimization
	stats_file = open(code+'_opt.stats', 'w')

	
	# Run MD. Write out a PDB structure every 10 steps during the run. 
	# Write stats every 10 steps
	md.optimize(atmsel, temperature=300, max_iterations=50,
            actions=[actions.write_structure(10, code+'.MD%04d.pdb'),
                     actions.trace(10, stats_file)])

	# Run CG, and write stats every 5 steps
	cg.optimize(atmsel, max_iterations=50,
            actions=[actions.trace(5, stats_file)])

	# Final energy
	final_mpdf = atmsel.energy()

	# Assess DOPE
	atmsel.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='TvLDH.profile',
              normalize_profile=True, smoothing_window=15)

	
	# Print the energies and the contributions
	initial_cont_all = dict(initial_mpdf[1])
	top_init_conts = dict(sorted(initial_cont_all.items(), key = itemgetter(1), reverse = True)[:5])
	
	l.info("\n\nThe initial energy of " + code + " is " + str(initial_mpdf[0]))
	print("\n\nThe initial energy of " + code + " is " + str(initial_mpdf[0]))
	print("The top 10 initial contributions the restraints are:\n")
	for keys, values in top_init_conts.items():
		print(keys, ":", values)

	final_cont_all = dict(final_mpdf[1])
	top_final_conts = dict(sorted(final_cont_all.items(), key = itemgetter(1), reverse = True)[:5])
	
	l.info("\n\nThe final energy of " + code + " is " + str(final_mpdf[0]))
	print("\n\nThe final energy of " + code + " is " + str(final_mpdf[0]))
	print("Final contributions the restraints are:\n")
	for keys, values in top_final_conts.items():
		print(keys, ":", values)

	
	mdl.write(file=code+'_optimized.pdb')
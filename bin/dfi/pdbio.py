#!/usr/bin/env python
# For reading and writing pdb files
from __future__ import print_function
from collections import namedtuple


def pdb_reader(filename, CAonly=False, noalc=True, chainA=False,
               chain_name='A', Verbose=False):
    """
    Reads in the ATOM entry of a pdb file. In the case of an NMR structure,
    the function reads in the first model.

    Input
    -----
    filename: file
       Filename of pdb file
    CAonly: bool
       Flag to only read the alpha-carbons.
    noalc: bool
       Flag to not read an alc
    chainA: bool
       Read only a single chain
    chain_name: str
       Name of chain to select (e.g., 'A')
    Verbose: str
       Flag for debugging

    Output
    ------
    ATOMS: ls
       ls of ATOM objects that make up the pdb.
    """

    ATOM = namedtuple('ATOM', ['atom_index', 'atom_name', 'alc',
                               'res_name', 'chainID', 'res_index',
                               'insert_code', 'x', 'y', 'z',
                               'occupancy', 'temp_factor', 'atom_type'])

    ATOMS = []
    readatoms = 0
    with open(filename) as pdb:
        for line in pdb:
            if line.startswith('ENDMDL'):
                if(Verbose):
                    print("MULTIPLE MODELS...USING MODEL1")
                return ATOMS

            if line.startswith('ATOM'):
                atom_index = int(line[7:11])
                atom_name = line[12:16]
                if CAonly and not(atom_name.strip() == 'CA'):
                    continue
                alc = line[16]  # alternate location
                if noalc and not((alc == ' ' or alc == 'A')):
                    continue
                res_name = line[17:20]
                chainID = line[21]
                if chainA and not(chainID == chain_name):
                    continue
                res_index = line[22:27].strip(' ')
                insert_code = line[26]
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                if(len(line) > 54):
                    occupancy = line[55:60]
                    temp_factor = float(line[61:66])
                    atom_type = line[77]
                else:
                    occupancy = 1.
                    temp_factor = 1.
                    atom_type = ' '
                ATOMS.append(ATOM(atom_index, atom_name, alc,
                                  res_name, chainID, res_index, insert_code,
                                  x, y, z, occupancy, temp_factor, atom_type))
                readatoms += 1

    print("Read %d atoms from the %s" % (readatoms, filename))
    return ATOMS


def pdb_writer(ATOMS, msg="HEADER  FROM PDBIO\n", filename="out.pdb",
               modelnum=1, atomoffset=0, residueoffset=0, mode="w"):

    with open(filename, mode) as pdb:
        pdb.write(msg)
        pdb.write("MODEL %d\n" % modelnum)
        pdb.write("PARENT N/A\n")
        for atom in ATOMS:
            record = 'ATOM  '
            atom_index = atom.atom_index + atomoffset
            atom_name = atom.atom_name
            alc = atom.alc
            res_name = atom.res_name
            chainID = atom.chainID
            res_index = atom.res_index
            iCode = atom.insert_code
            x = atom.x
            y = atom.y
            z = atom.z
            occupancy = 1.00
            temp_factor = atom.temp_factor
            atom_type = atom.atom_type
            atom1 = "{}{:5d} {:>4s}{:<1s}{:3s} {:<1s}{:4s}{:<1s}".format(
                record,
                atom_index,
                atom_name,
                alc,
                res_name,
                chainID,
                res_index,
                iCode)
            atom2 = "   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}  \n".format(
                x,
                y,
                z,
                occupancy,
                temp_factor,
                atom_type)
            pdb.write(atom1 + atom2)
        pdb.write("TER\n")
        pdb.write("END\n")
    print("Wrote out to file, %s" % filename)

#!/usr/bin/env python
"""
DFI (Dynamic Flexibility Index)
===============================

Description
------------
DFI Calculates the dynamics flexibility index
in order to study the conformational dynamics
of a protein.


Usage
-----
dfi_calc.py --pdb PDBFILE [--covar COVARFILE --chain CHAINID --fdfi RESNUMS]

Input
-----
PDBFILE:     PDBFILE
COVARFILE:    Covariance (Inverse Hessian) Matrix in a [NxN] ascii format
RESNUMS:     Chain + Residues number in the pdb, e.g. A15 B21

Output
------
* Structure used for DFI: -dficolor.pdb
* Master DFI: -dfianalysis.csv

Example
-------
```
./dfi_calc.py --pdb 1l2y.pdb [--covar covariance.dat --chain] A --fdfi A10]
```
"""
from __future__ import print_function, division
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import linalg as LA
from scipy import stats
from six.moves import range
from .pdbio import *
from .colordfi import colorbydfi


if __name__ == "__main__" and len(sys.argv) < 2:
    print(__doc__)
    exit()


def getcoords(ATOMS, Verbose=False):
    """
    Returns x,y and z numpy arrays of coordinates from
    ATOM object

    Input
    -----
    ATOMS: ATOM object
       Object for holding ATOM entries of pdb file
    Verbose: bool:
       Flag for Verbose Output

    Output
    ------
    (x,y,z): numpy
       numpy arrays of x,y,z

    """
    x = np.array(
        [atom.x for atom in ATOMS if atom.atom_name.strip() == 'CA'],
        dtype=float)
    y = np.array(
        [atom.y for atom in ATOMS if atom.atom_name.strip() == 'CA'],
        dtype=float)
    z = np.array(
        [atom.z for atom in ATOMS if atom.atom_name.strip() == 'CA'],
        dtype=float)

    return x, y, z


def calchessian(resnum, x, y, z, gamma, cutoff=None, Verbose=False):
    """
    Calculates the hessian and retuns the result

    Input
    ------
    resnum: int
       Number of residues
    x,y,z: numpy arrays
       Numpy array of coordinates
    gamma: int
       Value of spring constant (default set to 100)
    cutoff: float
       value of cutoff when using a distance based Hessian (default None)
    Verbose: bool
       Verbose Output for debug mode (default False).

    Output
    ------
    hess: numpy
       numpy array of the Hessian 3Nx3N shape
    """
    numresthree = 3 * resnum
    hess = np.zeros((numresthree, numresthree))

    if(Verbose):
        print("i,j,x1,y1,z1,x2,y2,z2,x_ij,y_ij,z_ij,r,k,g,cut")
    for i in range(resnum):
        for j in range(resnum):
            if i == j:
                continue
            x_i = x[i]
            y_i = y[i]
            z_i = z[i]
            x_j = x[j]
            y_j = y[j]
            z_j = z[j]
            x_ij = x_i - x_j
            y_ij = y_i - y_j
            z_ij = z_i - z_j
            r = x_ij * x_ij + y_ij * y_ij + z_ij * z_ij
            sprngcnst = (gamma * gamma * gamma) / (r * r * r)
            if(cutoff):
                if np.sqrt(r) > cutoff:
                    sprngcnst = 0.
            if(Verbose):
                print(','.join(np.array([i, j, x_i, y_i, z_i, x_j, y_j,
                                         z_j, x_ij, y_ij, z_ij, r,
                                         sprngcnst, gamma, cutoff],
                                        dtype=str)))

            # creation of Hii
            hess[3 * i, 3 * i] += sprngcnst * (x_ij * x_ij / r)
            hess[3 * i + 1, 3 * i + 1] += sprngcnst * (y_ij * y_ij / r)
            hess[3 * i + 2, 3 * i + 2] += sprngcnst * (z_ij * z_ij / r)

            hess[3 * i, 3 * i + 1] += sprngcnst * (x_ij * y_ij / r)
            hess[3 * i, 3 * i + 2] += sprngcnst * (x_ij * z_ij / r)
            hess[3 * i + 1, 3 * i] += sprngcnst * (y_ij * x_ij / r)

            hess[3 * i + 1, 3 * i + 2] += sprngcnst * (y_ij * z_ij / r)
            hess[3 * i + 2, 3 * i] += sprngcnst * (x_ij * z_ij / r)
            hess[3 * i + 2, 3 * i + 1] += sprngcnst * (y_ij * z_ij / r)

            # creation of Hij
            hess[3 * i, 3 * j] -= sprngcnst * (x_ij * x_ij / r)
            hess[3 * i + 1, 3 * j + 1] -= sprngcnst * (y_ij * y_ij / r)
            hess[3 * i + 2, 3 * j + 2] -= sprngcnst * (z_ij * z_ij / r)

            hess[3 * i, 3 * j + 1] -= sprngcnst * (x_ij * y_ij / r)
            hess[3 * i, 3 * j + 2] -= sprngcnst * (x_ij * z_ij / r)
            hess[3 * i + 1, 3 * j] -= sprngcnst * (y_ij * x_ij / r)

            hess[3 * i + 1, 3 * j + 2] -= sprngcnst * (y_ij * z_ij / r)
            hess[3 * i + 2, 3 * j] -= sprngcnst * (x_ij * z_ij / r)
            hess[3 * i + 2, 3 * j + 1] -= sprngcnst * (y_ij * z_ij / r)

    return hess


def flatandwrite(matrix, outfile):
    """Flattens out a matrix to a Nx1 column and write out to a file. """
    np.savetext(outfile, matrix.flatten())


def dfianal(fname, Array=False):
    """Calculate various dfi quantities and then output"""
    if not(Array):
        with open(fname, 'r') as infile:
            dfi = np.array([x.strip('\n') for x in infile], dtype=float)
    else:
        dfi = fname

    dfirel = dfi / np.mean(dfi)
    dfizscore = stats.zscore(dfi)
    dfiperc = pctrank(dfi)
    return dfi, dfirel, dfiperc, dfizscore


def pctrank(dfi, inverse=False):
    """
    Calculate %rank of DFI values

    Input
    -----
    dfi: numpy
       Array of dfi values
    inverse: bool
       Invert the pct ranking (default False)
    """
    if type(dfi).__module__ != 'numpy':
        raise ValueError('Input needs to be a numpy array')

    dfiperc = []
    lendfi = float(len(dfi))

    for m in dfi:
        if inverse:
            amt = np.sum(dfi >= m)
        else:
            amt = np.sum(dfi <= m)
        dfiperc.append(amt / lendfi)
    return np.array(dfiperc, dtype=float)


def calcperturbMat(invHrs, direct, resnum, Normalize=True):
    """
    Caclulates perturbation matrix used for dfi calculation.

    Input
    -----
    invHRS: numpy matrix
       covariance matrix (3N,3N), where N is the number of residues
    direct: numpy matrix
       matrix of peturbation directions
    resnum: int
       number of residues in protein
    Normalize: bool
       Normalize peturbation matrix

    Output
    ------
    peturbMat: numpy matrix
       NxN peturbation matrix where N is the number of residues

    """
    perturbMat = np.zeros((resnum, resnum))
    for k in range(len(direct)):
        peturbDir = direct[k, :]
        for j in range(int(resnum)):
            delforce = np.zeros(3 * resnum)
            delforce[3 * j:3 * j + 3] = peturbDir
            delXperbVex = np.dot(invHrs, delforce)
            delXperbMat = delXperbVex.reshape((resnum, 3))
            delRperbVec = np.sqrt(np.sum(delXperbMat * delXperbMat, axis=1))
            perturbMat[:, j] += delRperbVec[:]
    perturbMat /= 7

    if(Normalize):
        nrmlperturbMat = perturbMat / np.sum(perturbMat)
    else:
        print("WARNING: The perturbation matrix is not NORMALIZED")
        nrmlperturbMat = perturbMat

    return nrmlperturbMat


def chainresmap(ATOMS, Verbose=False):
    """
    Returns a dict object with the chainResNum as the key and the index
    of the atom
    """
    table = {}
    for i in range(len(ATOMS)):
        if ATOMS[i].chainID == ' ':
            entry = ATOMS[i].res_index
        else:
            # str(ATOMS[i].res_index)
            entry = ATOMS[i].chainID + ATOMS[i].res_index
        table[entry] = i
    if(Verbose):
        print(table)
    return table


def fdfiresf(ls_chain, table):
    """Returns numpy array of f-dfi res"""
    ls_ind = []

    for res in ls_chain:
        if res in table:
            ls_ind.append(table[res])
        else:
            print("WARNING: Can't find %s" % res)
            continue
    return np.array(ls_ind, dtype=int)


def fdfires_cords(fdfires, x, y, z):
    """
    Pull out the fdfi coordinates from (x,y,z)

    Input
    -----
    fdfires: ls
       indices of fdfires
    x,y,z: numpy
       numpy arrays of the coordinates

    Output
    ------
    nx3 matrix of f-DFI coordinates
    """

    return np.column_stack((x[fdfires], y[fdfires], z[fdfires]))


def rdist(r, fr):
    """
    Calculate the distance or r from fr

    Input
    ------
    r: numpy
       array of coordinates
    fr numpy
       nx3 matrix of f-DFI coordinates

    Output
    ------
    return rdist: array of distances from f-DFI sites
    """
    r_ij = fr - r
    rr = r_ij * r_ij
    return np.sqrt(rr.sum(axis=1))


def outputToDF(ATOMS, dfi, pctdfi, fdfi=None, pctfdfi=None, ls_ravg=None,
               ls_rmin=None, outfile=None, Verbose=True, writetofile=False):
    """
    Outputs the results of the DFI calculation to a DataFrame and csv file

    Input
    -----
    ATOMS: ATOM object
       Object to hold ATOM entries of PDB files
    dfi: numpy
       numpy array of dfi values
    pctdfi: numpy
       numpy array of pctdfi values
    fdfi: numpy
       numpy array of fdfi values
    pctdfi: numpy
       numpy array of pctdfi values
    ls_ravg: ls
       list of the average distance of a residue
       to all f-dfi residues
    ls_rmin: ls
       list of the min distance of a residue
       to all f-dfi residues
    outfile: str
       Name of file to write out the DataFrame in csv format
    Verbose: bool
       Output for debugging
    writetofile: bool
       If True will write out to file, otherwise just return the
       df.

    Output
    ------
    df_dfi: DataFrame
       DataFrame object containing all the inputs

    """
    mapres = {'ALA': 'A',
              'CYS': 'C',
              'ASP': 'D',
              'GLU': 'E',
              'PHE': 'F',
              'GLY': 'G',
              'HIS': 'H',
              'ILE': 'I',
              'LYS': 'K',
              'LEU': 'L',
              'MET': 'M',
              'PRO': 'P',
              'ARG': 'R',
              'GLN': 'Q',
              'ASN': 'N',
              'SER': 'S',
              'THR': 'T',
              'TRP': 'W',
              'TYR': 'Y',
              'VAL': 'V',
              'MSE': 'M'}
    dfx = pd.DataFrame()
    dfx['ResID'] = [ATOMS[i].res_index.strip(' ') for i in range(len(ATOMS))]
    dfx['dfi'] = dfi
    dfx['pctdfi'] = pctdfi
    dfx['ChainID'] = [ATOMS[i].chainID for i in range(len(ATOMS))]
    dfx['Res'] = [ATOMS[i].res_name for i in range(len(ATOMS))]
    dfx['R'] = dfx['Res'].map(mapres)

    if type(fdfi).__module__ == 'numpy':
        dfx['fdfi'] = fdfi
        dfx['pctfdfi'] = pctfdfi
        dfx['ravg'] = ls_ravg
        dfx['rmin'] = ls_rmin
        mask = (dfx['rmin'] > 8.0) & (dfx['pctfdfi'] > 0.75)
        dfx['A'] = mask.map(lambda x: 'A' if x else 'NotA')

    if(writetofile):
        dfx.to_csv(outfile, index=False)
        print("Wrote out to %s" % (outfile))
    return dfx


def top_quartile_pos(pctfdfi):
    """
    returns a list of indices of positions in the top quartile of pctdfi

    Input
    -----
    pctdfi: numpy
       numpy array of pctdfi or pctfdfi values

    Output
    ------
    top_quartile: ls
       ls of indices thare in the the top quartile
    """
    return [i for i, val in enumerate(pctfdfi) if val > 0.75]


def check_args(args=None):
    """
    Parse command lines input

    Input
    -----
    ls of command line input

    Output
    ------
    pdbfile: file
       name of pdb file to run dfi calculation
    pdbid: 'str'
       4 Letter PDB code
    mdhess: file
        name of file that contains hessian matrix from MD
    ls_reschain: list
       list of f-dfi Residues (e.g., ['A17','A19'])
    chain_name: str
       list of chain to select (Depracated)

    """
    parser = argparse.ArgumentParser(
        description='DFI CLI')
    parser.add_argument('--pdb',
                        help='PDB File to run DFI',
                        required=True)
    parser.add_argument('--covar',
                        help='3Nx3N covariance matrix')
    parser.add_argument('--chain',
                        help='port of the web server')
    parser.add_argument('--fdfi',
                        help='f-DFI residues',
                        nargs='+')

    results = parser.parse_args(args)
    pdbid = results.pdb.split('.')[0]
    return results.pdb, pdbid, results.covar, results.fdfi, results.chain


def _writeout_eigevalues(e_vals, eigenfile):
    """
    Write out eigenvalues from the Hessian Matrix

    Input:
    e_vals: numpy
       array of eigenfiles
    eigenfiel: str
       eigenfile name
    """
    with open(eigenfile, 'w') as outfile:
        for i, val in enumerate(np.sort(e_vals)):
            outfile.write("%d\t%f\n" % (i, np.real(val)))


def calc_covariance(numres, x, y, z, invhessfile=None, Verbose=False,
                    eigenfile=None):
    """
    Calculates the covariance matrix by first
    calculating the hessian from coordinates and then
    inverting it.

    Input
    -----
    numres: int
       number of residues
    x,y,z: numpy
       numpy array of coordinates
    Verbose: bool
       flag for debugging and writing out contents
       of covariance or the inverse Hessian

    Output
    ------
    invHRS: numpy
       (3*numres,3*numres) matrix

    """
    gamma = 100
    hess = calchessian(numres, x, y, z, gamma, Verbose)
    if(Verbose):
        print("Hessian")
        print(hess)
        flatandwrite(hess, 'hesspy.debug')
        e_vals, e_vecs = LA.eig(hess)
        _writeout_eigevalues(e_vals, eigenfile)

    U, w, Vt = LA.svd(hess, full_matrices=False)
    S = LA.diagsvd(w, len(w), len(w))
    assert np.allclose(hess, np.dot(U, np.dot(S, Vt))), "SVD didn't go well"
    if(Verbose):
        flatandwrite(U, 'Upy-test.debug')
        flatandwrite(w, 'wpy-test.debug')
        flatandwrite(Vt, 'Vtpy-test.debug')

    # the near zero eigenvalues blowup the inversion so
    # we will truncate them and add a small amount of bias
    tol = 1e-6
    singular = w < tol
    invw = 1 / w
    invw[singular] = 0.
    invHrs = np.dot(np.dot(U, np.diag(invw)), Vt)
    if(Verbose):
        flatandwrite(invHrs, invhessfile)
    assert np.sum(
       singular) == 6., "# of near-singular eigenvals: %f" % np.sum(singular)
    return invHrs


def calc_dfi(pdbfile, pdbid=None, covar=None, ls_reschain=[], chain_name=None,
             Verbose=False, writetofile=False, colorpdb=False,
             dfianalfile=None):
    """Main function for calculating DFI

    Inputs
    ------
    pdbfile: file
       PDB File for dfi calculation
    pdbid: str
       4 character PDBID from PDB
    covar: file
       hessian file obtained from MD
    ls_reschain: ls
       list of f-dfi residues by chain and index (e.g., ['A19','A20']
    chain_name: str
       chain name (e.g., A) to pull out specific chain of the PDB
    Verbose: bool
       switch for debugging
    writefofile: bool
       If True will writeou to a csv file
    colorpdb: bool
       If True will output a colorpdb
    dfianalfile: str
       Name of custom output file. This is useful for when you may
       want to number outputs using different covariance matrices
       that correspond to different time windows.

    Output
    ------
    df_dfi: DataFrame
       DataFrame object for DFI values
    """
    if(not(pdbid)):
        pdbid = pdbfile.split('.')[0]
    eigenfile = pdbid + '-eigenvalues.txt'
    invhessfile = pdbid + '-pinv_svd.debug'
    if(not(dfianalfile)):
        dfianalfile = pdbid + '-dfianalysis.csv'

    # read in the pdb file
    ATOMS = pdb_reader(pdbfile, CAonly=True, noalc=True, chainA=False,
                             chain_name=chain_name, Verbose=False)
    x, y, z = getcoords(ATOMS)
    numres = len(ATOMS)

    # create covariance matrix or read it in if provided
    if not(covar):
        invHrs = calc_covariance(numres, x, y, z, Verbose=False,
                                 eigenfile=eigenfile,
                                 invhessfile=invhessfile)
    else:  # this is where we load the Hessian if provided
        invHrs = np.loadtxt(covar)

    # RUN DFI
    directions = np.vstack(([1, 0, 0], [0, 1, 0], [0, 0, 1], [
                           1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]))
    normL = np.linalg.norm(directions, axis=1)
    direct = directions / normL[:, None]
    nrmlperturbMat = calcperturbMat(invHrs, direct, numres)
    dfi = np.sum(nrmlperturbMat, axis=1)
    dfi, reldfi, pctdfi, zscoredfi = dfianal(dfi, Array=True)

    # f-dfi
    if ls_reschain:
        # find the f-dfi residues
        fdfiset = set(ls_reschain)
        ls_reschain = list(fdfiset)
        ls_reschain.sort()
        fdfires = np.sort(fdfiresf(ls_reschain, chainresmap(ATOMS)))
        # calculate f-dfi
        fdfitop = np.sum(nrmlperturbMat[:, fdfires], axis=1) / len(fdfires)
        fdfibot = np.sum(nrmlperturbMat, axis=1) / len(nrmlperturbMat)
        fdfi, relfdfi, pctfdfi, zscorefdfi = dfianal(
            fdfitop / fdfibot, Array=True)
        rlist = np.column_stack((x, y, z))  # dump into a list
        # get coordinates of the f-dfi residues
        fr = fdfires_cords(fdfires, x, y, z)
        ls_ravg = np.array([rdist(r, fr).mean() for r in rlist])
        ls_rmin = np.array([rdist(r, fr).min() for r in rlist])

    # output to dataframe
    if len(ls_reschain) > 0:
        df_dfi = outputToDF(ATOMS, dfi, pctdfi, fdfi=fdfi, pctfdfi=pctfdfi,
                            ls_ravg=ls_ravg, ls_rmin=ls_rmin,
                            outfile=dfianalfile,
                            writetofile=writetofile)
    else:
        df_dfi = outputToDF(ATOMS, dfi, pctdfi,
                            outfile=dfianalfile, writetofile=writetofile)

    # output to ColoredDFI Files
    if(colorpdb):
        colorbydfi(
            dfianalfile, pdbfile, colorbyparam='pctdfi',
            outfile=pdbid + '-dficolor.pdb')

        if len(ls_reschain) > 0:
            colorbydfi(
                dfianalfile, pdbfile, colorbyparam='pctfdfi',
                outfile=pdbid + '-fdficolor.pdb')

    if not(writetofile):
        return df_dfi


if __name__ == "__main__":
    pdbfile, pdbid, covar, ls_reschain, chain_name = check_args(
        sys.argv[1:])
    print("Processing %s" % pdbfile)
    df_dfi = calc_dfi(pdbfile, pdbid, covar=covar, ls_reschain=ls_reschain,
                      chain_name=chain_name,
                      writetofile=True, colorpdb=True)

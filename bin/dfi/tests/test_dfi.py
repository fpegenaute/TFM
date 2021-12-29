import numpy as np
import dfi
from dfi.datafiles import example_pdb
from dfi.dfi_calc import check_args


def test_pctrank_ascending():
    a = np.array([1, 2, 3, 4, 5])
    assert np.all(dfi.dfi_calc.pctrank(a) == np.array(
        [0.2, 0.4, 0.6, 0.8, 1.0], dtype=float))


def test_pctrank_descending():
    a = np.array([1, 2, 3, 4, 5])
    assert np.all(dfi.dfi_calc.pctrank(a, inverse=True) ==
                  np.array([1., 0.8, 0.6, 0.4, 0.2], dtype=float))


def test_comlineargs():
    comline = ['--pdb', example_pdb, '--fdfi',
               'A19', 'A10', '--covar', 'mwcovar.dat']

    pdbfile, pdbid, covar, ls_reschain, chain_name = check_args(comline)

    assert pdbfile == example_pdb
    assert covar == 'mwcovar.dat'
    assert np.all(ls_reschain == np.array(['A19', 'A10']))

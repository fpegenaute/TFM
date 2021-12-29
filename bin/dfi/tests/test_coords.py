import dfi
from dfi.datafiles import xyz_npy, example_pdb

import numpy as np


def test_coords():
    xyz = np.load(xyz_npy)
    ATOMS = dfi.pdbio.pdb_reader(example_pdb, CAonly=True)
    x, y, z = dfi.dfi_calc.getcoords(ATOMS)
    xyz_test = np.vstack((x, y, z))
    assert np.all(xyz == xyz_test)

#!/usr/bin/env python
import dfi.pdbio
import numpy as np
from dfi.datafiles import example_pdb, test_pdb


class TestPDBIO():
    ATOMS = dfi.pdbio.pdb_reader(example_pdb,
                                 CAonly=True,
                                 Verbose=True)

    def test_coordinates(self):
        x = -8.608
        y = 3.135
        z = -1.618
        assert np.isclose(self.ATOMS[0].x, x)
        assert np.isclose(self.ATOMS[0].y, y)
        assert np.isclose(self.ATOMS[0].z, z)
        assert 'CA' == self.ATOMS[0].atom_name.strip()
        assert 'ASN' == self.ATOMS[0].res_name
        assert 'A' == self.ATOMS[0].chainID

    def test_writer(self):
        ATOMS = dfi.pdbio.pdb_reader(test_pdb, CAonly=True)
        x = -8.608
        y = 3.135
        z = -1.618
        print(ATOMS[0].x, x)
        assert np.isclose(ATOMS[0].x, x)
        assert np.isclose(ATOMS[0].y, y)
        assert np.isclose(ATOMS[0].z, z)
        assert 'CA' == ATOMS[0].atom_name.strip()
        assert 'ASN' == ATOMS[0].res_name
        assert 'A' == ATOMS[0].chainID

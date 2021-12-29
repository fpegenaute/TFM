from dfi.download_pdb import fetch_pdb


class TestPDBDownloader():
    """
    Class for downloading pdb
    """
    pdbid = '1l2y'

    def test_pdb_stringIO(self):
        pdbid = TestPDBDownloader.pdbid
        pdb = fetch_pdb(pdbid, writetofile=False, Verbose=True)

    def test_pdb_file(self):
        pdbid = TestPDBDownloader.pdbid
        fetch_pdb(pdbid, writetofile=True, Verbose=True)

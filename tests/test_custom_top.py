from dis import get_instructions
import unittest
import pandas as pd
##  run this as python -m tests.test_custom_top on the root dir

import sys
sys.path.append(".")

from bin.custom_top import RigidBody


class TestRigidbody(unittest.TestCase):

    def setUp(self):
        self.rb1 = RigidBody(resolution="all",
            molecule_name= "5lg4_B", 
            color="blue" , 
            fasta_fn="input_fasta/SEC3.fasta", 
            # fasta_id=fasta, 
            pdb_fn="RESULTS/output/SEC3/PDB/partial/5lg4_B.pdb", 
            chain="B",
            residue_range= (74,249), 
            rigid_body=1, 
            super_rigid_body="", 
            chain_of_super_rigid_bodies="", 
            bead_size=10,
            em_residues_per_gaussian=0, 
            type="Experimental")
        
        self.rb2 = RigidBody(resolution="all",
            molecule_name= "SEC3_AF_A2_AF_2", 
            color="orange" , 
            fasta_fn="input_fasta/SEC3.fasta", 
            # fasta_id=fasta, 
            pdb_fn="RESULTS/output/SEC3/ALPHAFOLD/DOMAINS/SEC3_AF_A2_AF.pdb", 
            chain="2",
            residue_range=(325,450) , 
            rigid_body=2, 
            super_rigid_body="", 
            chain_of_super_rigid_bodies="", 
            bead_size=20,
            em_residues_per_gaussian=0, 
            type="AF_model")
        
        self.rb3 = RigidBody(resolution="all",
            molecule_name= "SEC3_RF_A1_RF", 
            color="orange" , 
            fasta_fn="input_fasta/SEC3.fasta", 
            # fasta_id=fasta, 
            pdb_fn="RESULTS/output/SEC3/ROSETTAFOLD/DOMAINS/model_4-crderr_A1_RF.pdb", 
            chain="1",
            residue_range=(135,230) , 
            rigid_body=3, 
            super_rigid_body="", 
            chain_of_super_rigid_bodies="", 
            bead_size=20,
            em_residues_per_gaussian=0, 
            type="RF_model")
    

    # Testing some of the automatically generated attributes
    def test_pdb_offset(self):
        self.assertEqual(-73, self.rb1.pdb_offset)
        self.assertEqual(-324, self.rb2.pdb_offset)
        self.assertEqual(-134, self.rb3.pdb_offset)

    def test_fasta_fn(self):
        expected_fasta = "input_fasta/SEC3.fasta"
        self.assertEqual(expected_fasta, self.rb1.fasta_fn)

    def test_fasta_id(self):
        expected_fasta = "SEC3"
        self.assertEqual(expected_fasta, self.rb1.fasta_id)

    # Tests for the methods
    def test_get_ResIDs(self):
        resids = self.rb1.get_resIDs()
        self.assertGreater(len(resids), 0)
        self.assertEqual(len(self.rb3.get_resIDs()), 62)
    
    def test_update_overlap(self):
        self.rb1.update_overlap(self.rb3)
        self.assertEqual(len(self.rb1.overlap), 0)
        self.assertEqual(len(self.rb2.overlap), 0)
        self.assertEqual(len(self.rb3.overlap), 62)
    
    def test_get_length(self):
        self.assertEqual(self.rb1.get_length(), 175)
    
    def test_extract_avg_BFactor(self):
        self.assertAlmostEqual(
            self.rb2.extract_avg_BFactor(), 24.993650793650797, places=6
            )
    
    def test_get_coverage(self):
        coverage_df = self.rb1.get_coverage()
        self.assertIsInstance(coverage_df, pd.DataFrame)
        self.assertEqual(coverage_df.columns[0], "ResID")
        self.assertCountEqual(
            coverage_df[coverage_df.columns[1]].unique(), [0,1]
            )
        self.assertEqual(
            len(coverage_df[coverage_df.columns[0]]), 
            len(coverage_df[coverage_df.columns[1]])
            )
    def test_split_rb_hinges(self):
        #Use case: 0 hinges, same RB
        rb_list = self.rb1.split_rb_hinges([])
        res_ranges = [rb.residue_range for rb in rb_list]
        self.assertCountEqual(res_ranges, [(74, 249)])

        # Use case: 1 hinge, 2 rbs
        rb_list = self.rb1.split_rb_hinges([(100,130)])
        res_ranges = [rb.residue_range for rb in rb_list]
        self.assertCountEqual(res_ranges, [(0, 100), (130, 249)])

        # Use case: >1 hinges
        rb_list = self.rb1.split_rb_hinges([(76,86), (100,130), (140,220)])
        res_ranges = [rb.residue_range for rb in rb_list] 
        res_ranges.sort(key=lambda tup: tup[0])
        self.assertCountEqual(res_ranges, [(0, 76), (86, 100), (130, 140), (220, 249)])


    
        
if __name__ == '__main__':
    print(__name__)
    unittest.main()
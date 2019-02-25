"""
Unit and regression test for the kinomodel package.
"""

# Import package, test suite, and other packages as needed
import unittest
import numpy as np

class KinomodelTestCase(unittest.TestCase):

    #def test_kinomodel_installed():
    # make sure kinomodel has been correctly installed and importable
    #self.assertTrue("kinomodel" in sys.modules)

    def test_command(self):
        # absolute import (with kinomodel installed)        
        #from kinomodel.features import featurize
        #JG (temperary)
        from features import featurize
        with self.assertRaises(ValueError):
            featurize(chain=0, coord='pdb', feature='interact', pdb='3PP0')

    def test_features(self):
        # absolute import (with kinomodel installed)
        #from kinomodel.features import featurize
        #JG (temperary)
        from features import featurize

        # example 1: a kinase with no gap(s) in the binding pocket residues
        mean_dist = featurize(chain='A', coord='pdb', feature='interact', pdb='3PP0')
        self.assertEqual(
            round(
                np.asscalar(mean_dist[0]), 7),
            1.3685706)  # mean ligand-pocket distance

        # example 2: a kinase with gap(s) in the binding pocket residues
        mean_dist = featurize(chain='A', coord='pdb', feature='interact', pdb='3RCD')
        self.assertEqual(
            round(
                np.asscalar(mean_dist[0]), 7),
            1.4069413)  # mean ligand-pocket distance

        # example 3: a kinase with multiple occupancy
        mean_dist = featurize(chain='A', coord='pdb', feature='interact', pdb='1M17')
        self.assertEqual(
            round(
                np.asscalar(mean_dist[0]), 7),
            1.6548363)  # mean ligand-pocket distance

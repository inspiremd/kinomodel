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
        # make sure the input command line is expected
        with self.assertRaises(ValueError):
            featurize(chain=0, coord='pdb', feature='conf', pdb='3PP0')

    def test_features(self):
        # absolute import (with kinomodel installed)
        #from kinomodel.features import featurize

        #JG (temperary)
        from features import featurize

        # example 1: a kinase with no gap(s) in the binding pocket residues
        (key_res, dihedrals, distances) = featurize(chain='A', coord='pdb', feature='conf', pdb='3PP0')
        self.assertEqual(key_res,
            [767, 775, 836, 838, 753, 770, 774, 864, 862, 863, 873, 814])
        self.assertEqual(
            round(
                np.asscalar(dihedrals[0][0]), 7),
            -2.0228872)  # the first dihedral value
        self.assertEqual(
            round(
                np.asscalar(distances[0][0]), 7),
            0.7770488)  # the first distance value

        # example 2: a kinase with gap(s) in the binding pocket residues
        (key_res, dihedrals, distances) = featurize(chain='A', coord='pdb', feature='conf', pdb='3RCD')
        self.assertEqual(key_res,
            [767, 775, 836, 838, 753, 770, 774, 864, 862, 863, 873, 814])
        self.assertEqual(
            round(
                np.asscalar(dihedrals[0][0]), 7),
            -2.1615183)  # the first dihedral value
        self.assertEqual(
            round(
                np.asscalar(distances[0][0]), 7),
            1.0546854)  # the first distance value

        # example 3: a kinase with multiple occupancy
        (key_res, dihedrals, distances) = featurize(chain='A', coord='pdb', feature='conf', pdb='1M17')
        self.assertEqual(key_res,
            [735, 743, 804, 806, 721, 738, 742, 832, 830, 831, 841, 782])
        self.assertEqual(
            round(
                np.asscalar(dihedrals[0][0]), 7),
            -2.4006705)  # the first dihedral value
        self.assertEqual(
            round(
                np.asscalar(distances[0][0]), 7),
            0.3558538)  # the first distance value

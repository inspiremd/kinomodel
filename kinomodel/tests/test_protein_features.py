'''
Unit and regression test for the protein_features module.
'''
# Import package, test suite, and other packages as needed
import sys
import numpy as np
import protein_features as pf
import unittest
from unittest.mock import MagicMock, patch

# clean traceback msgs
sys.tracebacklimit = 0


class ProteinTestCase(unittest.TestCase):
    def test_type(self):
        # make sure the input fails to be read in when it is not a string
        pf.input = MagicMock(return_value=7)
        with self.assertRaises(ValueError):
            pf.basics()

    def test_comma(self):
        # make sure the input fails to be read in when it's not separated by comma
        pf.input = MagicMock(return_value='3pp0A')
        with self.assertRaises(ValueError):
            pf.basics()

    def test_basics(self):
        ## make sure the basic info of a kinase is retrieved based on its pdb code and chain index
        # example 1: a kinase with no gap(s) in the binding pocket residues
        pf.input = MagicMock(return_value='3pp0,A')
        self.assertEqual(pf.basics()[1], 407)
        self.assertEqual(pf.basics()[2], 'ErbB2')
        self.assertEqual(pf.basics()[3], 4820)
        self.assertEqual(
            pf.basics()[4],
            'KVLGSGAFGTVYKVAIKVLEILDEAYVMAGVGPYVSRLLGIQLVTQLMPYGCLLDHVREYLEDVRLVHRDLAARNVLVITDFGLA'
        )
        self.assertEqual(pf.basics()[5], [
            724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736,
            750, 751, 752, 753, 754, 755, 766, 767, 768, 769, 770, 771, 772,
            773, 774, 775, 776, 777, 778, 780, 781, 782, 783, 784, 785, 786,
            787, 788, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805,
            806, 807, 808, 809, 810, 811, 812, 835, 836, 837, 838, 839, 840,
            841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,
            861, 862, 863, 864, 865, 866, 867
        ])
        self.assertEqual(
            pf.basics()[6],
            [767, 775, 836, 838, 753, 770, 774, 864, 862, 863, 873, 814])

        # example 2: a kinase with gap(s) in the binding pocket residues
        pf.input = MagicMock(return_value='3rcd,A')
        self.assertEqual(pf.basics()[1], 407)
        self.assertEqual(pf.basics()[2], 'ErbB2')
        self.assertEqual(pf.basics()[3], 9325)
        self.assertEqual(
            pf.basics()[4],
            'KVLGSGAFGTVYKVAIKVLEILDEAYVMAGVGPYVSRLLGIQLVTQLMPYGCLLDHVREYLEDVRLVHRDLAARNVLVITDFGL_'
        )
        self.assertEqual(pf.basics()[5], [
            724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736,
            750, 751, 752, 753, 754, 755, 766, 767, 768, 769, 770, 771, 772,
            773, 774, 775, 776, 777, 778, 780, 781, 782, 783, 784, 785, 786,
            787, 788, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805,
            806, 807, 808, 809, 810, 811, 812, 835, 836, 837, 838, 839, 840,
            841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,
            861, 862, 863, 864, 865, 866, 0
        ])
        self.assertEqual(
            pf.basics()[6],
            [767, 775, 836, 838, 753, 770, 774, 864, 862, 863, 873, 814])

    def test_features(self):
        ## make sure the features (dihedrals, distances) are correctly calculated for pdb structures
        # example 1: a kinase with no gap(s) in the binding pocket residues
        pf.input = MagicMock(return_value='pdb')
        pdb_chainid = ('3pp0', 'A')
        numbering = [
            724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736,
            750, 751, 752, 753, 754, 755, 766, 767, 768, 769, 770, 771, 772,
            773, 774, 775, 776, 777, 778, 780, 781, 782, 783, 784, 785, 786,
            787, 788, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805,
            806, 807, 808, 809, 810, 811, 812, 835, 836, 837, 838, 839, 840,
            841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,
            861, 862, 863, 864, 865, 866, 867
        ]
        self.assertEqual(
            round(
                np.asscalar(pf.features(pdb_chainid, numbering)[0][0][0]), 7),
            -2.0228872)  # the first dihedral value
        self.assertEqual(
            round(
                np.asscalar(pf.features(pdb_chainid, numbering)[1][0][0]), 7),
            0.7770488)  # the first distance value

        # example 2: a kinase with gap(s) in the binding pocket residues
        pf.input = MagicMock(return_value='pdb')
        pdb_chainid = ('3rcd', 'A')
        numbering = [
            724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736,
            750, 751, 752, 753, 754, 755, 766, 767, 768, 769, 770, 771, 772,
            773, 774, 775, 776, 777, 778, 780, 781, 782, 783, 784, 785, 786,
            787, 788, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805,
            806, 807, 808, 809, 810, 811, 812, 835, 836, 837, 838, 839, 840,
            841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,
            861, 862, 863, 864, 865, 866, 0
        ]
        self.assertEqual(
            round(
                np.asscalar(pf.features(pdb_chainid, numbering)[0][0][0]), 7),
            -2.1615183)  # the first dihedral value
        self.assertEqual(
            round(
                np.asscalar(pf.features(pdb_chainid, numbering)[1][0][0]), 7),
            1.0546854)  # the first distance value
        self.assertEqual(len(pf.features(pdb_chainid, numbering)[1][0]), 4)


if __name__ == '__main__':
    unittest.main()

'''
Unit and regression test for the interact_features module.
'''
# Import package, test suite, and other packages as needed
import sys
import numpy as np
import interact_features as inf
import unittest
from unittest.mock import MagicMock, patch

# clean traceback msgs
sys.tracebacklimit = 0


class InteractTestCase(unittest.TestCase):
    def test_type(self):
        # make sure the input fails to be read in when it is not a string
        inf.input = MagicMock(return_value=7)
        with self.assertRaises(ValueError):
            inf.basics()

    def test_comma(self):
        # make sure the input fails to be read in when it's not separated by comma
        inf.input = MagicMock(return_value='3pp0A')
        with self.assertRaises(ValueError):
            inf.basics()

    def test_basics(self):
        ## make sure the basic info of a kinase is retrieved based on its pdb code and chain index
        # example 1: a kinase with no gap(s) in the binding pocket residues
        inf.input = MagicMock(return_value='3pp0,A')
        self.assertEqual(inf.basics()[1], 407)
        self.assertEqual(inf.basics()[2], 'ErbB2')
        self.assertEqual(inf.basics()[3], 4820)
        self.assertEqual(inf.basics()[4], '03Q')
        self.assertEqual(
            inf.basics()[5],
            'KVLGSGAFGTVYKVAIKVLEILDEAYVMAGVGPYVSRLLGIQLVTQLMPYGCLLDHVREYLEDVRLVHRDLAARNVLVITDFGLA'
        )
        self.assertEqual(inf.basics()[6], [
            724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736,
            750, 751, 752, 753, 754, 755, 766, 767, 768, 769, 770, 771, 772,
            773, 774, 775, 776, 777, 778, 780, 781, 782, 783, 784, 785, 786,
            787, 788, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805,
            806, 807, 808, 809, 810, 811, 812, 835, 836, 837, 838, 839, 840,
            841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,
            861, 862, 863, 864, 865, 866, 867
        ])

        # example 2: a kinase with gap(s) in the binding pocket residues
        inf.input = MagicMock(return_value='3rcd,A')
        self.assertEqual(inf.basics()[1], 407)
        self.assertEqual(inf.basics()[2], 'ErbB2')
        self.assertEqual(inf.basics()[3], 9325)
        self.assertEqual(inf.basics()[4], '03P')
        self.assertEqual(
            inf.basics()[5],
            'KVLGSGAFGTVYKVAIKVLEILDEAYVMAGVGPYVSRLLGIQLVTQLMPYGCLLDHVREYLEDVRLVHRDLAARNVLVITDFGL_'
        )
        self.assertEqual(inf.basics()[6], [
            724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736,
            750, 751, 752, 753, 754, 755, 766, 767, 768, 769, 770, 771, 772,
            773, 774, 775, 776, 777, 778, 780, 781, 782, 783, 784, 785, 786,
            787, 788, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805,
            806, 807, 808, 809, 810, 811, 812, 835, 836, 837, 838, 839, 840,
            841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,
            861, 862, 863, 864, 865, 866, 0
        ])

    def test_features(self):
        ## make sure the features (mean distance between ligand and pocket residues) are correctly calculated for pdb structures
        # example 1: a kinase with no gap(s) in the binding pocket residues
        inf.input = MagicMock(return_value='pdb')
        pdb_chainid = ('3pp0', 'A')
        ligand = '03Q'
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
                np.asscalar(inf.features(pdb_chainid, ligand, numbering)), 7),
            2.5292468)  # mean ligand-pocket distance

        # example 2: a kinase with gap(s) in the binding pocket residues
        inf.input = MagicMock(return_value='pdb')
        pdb_chainid = ('3rcd', 'A')
        ligand = '03P'
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
                np.asscalar(inf.features(pdb_chainid, ligand, numbering)), 7),
            1.5608727)  # mean ligand-pocket distance


if __name__ == '__main__':
    unittest.main()

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

        # example 3: a kinase with multiple occupancy
        inf.input = MagicMock(return_value='1m17,A')
        self.assertEqual(inf.basics()[1], 406)
        self.assertEqual(inf.basics()[2], 'EGFR')
        self.assertEqual(inf.basics()[3], 873)
        self.assertEqual(inf.basics()[4], 'AQ4')
        self.assertEqual(
            inf.basics()[5],
            'KVLGSGAFGTVYKVAIKELEILDEAYVMASVDPHVCRLLGIQLITQLMPFGCLLDYVREYLEDRRLVHRDLAARNVLVITDFGLA'
        )
        self.assertEqual(inf.basics()[6], [
            692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704,
            718, 719, 720, 721, 722, 723, 734, 735, 736, 737, 738, 739, 740,
            741, 742, 743, 744, 745, 746, 748, 749, 750, 751, 752, 753, 754,
            755, 756, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773,
            774, 775, 776, 777, 778, 779, 780, 803, 804, 805, 806, 807, 808,
            809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821,
            829, 830, 831, 832, 833, 834, 835
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
            1.3685706)  # mean ligand-pocket distance

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
            1.4236906)  # mean ligand-pocket distance

        # example 3: a kinase with multiple occupancy
        inf.input = MagicMock(return_value='pdb')
        pdb_chainid = ('1m17', 'A')
        ligand = 'AQ4'
        numbering = [
            692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704,
            718, 719, 720, 721, 722, 723, 734, 735, 736, 737, 738, 739, 740,
            741, 742, 743, 744, 745, 746, 748, 749, 750, 751, 752, 753, 754,
            755, 756, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773,
            774, 775, 776, 777, 778, 779, 780, 803, 804, 805, 806, 807, 808,
            809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821,
            829, 830, 831, 832, 833, 834, 835
        ]
        self.assertEqual(
            round(
                np.asscalar(inf.features(pdb_chainid, ligand, numbering)), 7),
            1.6548363)  # mean ligand-pocket distance


if __name__ == '__main__':
    unittest.main()

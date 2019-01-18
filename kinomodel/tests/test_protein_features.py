"""
Unit and regression test for the kinomodel package.
"""

# Import package, test suite, and other packages as needed
import kinomodel as ki
import unittest
import sys
import argparse as ap
import numpy as np

# clean traceback msgs
sys.tracebacklimit = 0


class KinomodelTestCase(unittest.TestCase):

    #def test_kinomodel_installed():
    # make sure kinomodel has been correctly installed and importable
    #self.assertTrue("kinomodel" in sys.modules)

    def test_command(self):
        # make sure the input command line is expected
        with self.assertRaises(ValueError):
            ki.main(
                ap.Namespace(chain=0, coord='pdb', feature='conf', pdb='3PP0'))

    def test_basics(self):
        # example 1: a kinase with no gap(s) in the binding pocket residues
        self.kinase = ki.main(
            ap.Namespace(chain='A', coord='pdb', feature='conf', pdb='3PP0'))
        self.assertEqual(self.kinase.kinase_id, 407)
        self.assertEqual(self.kinase.name, 'ErbB2')
        self.assertEqual(self.kinase.struct_id, 4820)
        self.assertEqual(self.kinase.ligand, 'N/A')
        self.assertEqual(
            self.kinase.pocket_seq,
            'KVLGSGAFGTVYKVAIKVLEILDEAYVMAGVGPYVSRLLGIQLVTQLMPYGCLLDHVREYLEDVRLVHRDLAARNVLVITDFGLA'
        )
        self.assertEqual(self.kinase.numbering, [
            724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736,
            750, 751, 752, 753, 754, 755, 766, 767, 768, 769, 770, 771, 772,
            773, 774, 775, 776, 777, 778, 780, 781, 782, 783, 784, 785, 786,
            787, 788, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805,
            806, 807, 808, 809, 810, 811, 812, 835, 836, 837, 838, 839, 840,
            841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,
            861, 862, 863, 864, 865, 866, 867
        ])
        self.assertEqual(
            self.kinase.key_res,
            [767, 775, 836, 838, 753, 770, 774, 864, 862, 863, 873, 814])

        # example 2: a kinase with gap(s) in the binding pocket residues
        self.kinase = ki.main(
            ap.Namespace(chain='A', coord='pdb', feature='conf', pdb='3RCD'))
        self.assertEqual(self.kinase.kinase_id, 407)
        self.assertEqual(self.kinase.name, 'ErbB2')
        self.assertEqual(self.kinase.struct_id, 9325)
        self.assertEqual(self.kinase.ligand, 'N/A')
        self.assertEqual(
            self.kinase.pocket_seq,
            'KVLGSGAFGTVYKVAIKVLEILDEAYVMAGVGPYVSRLLGIQLVTQLMPYGCLLDHVREYLEDVRLVHRDLAARNVLVITDFGL_'
        )
        self.assertEqual(self.kinase.numbering, [
            724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736,
            750, 751, 752, 753, 754, 755, 766, 767, 768, 769, 770, 771, 772,
            773, 774, 775, 776, 777, 778, 780, 781, 782, 783, 784, 785, 786,
            787, 788, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 805,
            806, 807, 808, 809, 810, 811, 812, 835, 836, 837, 838, 839, 840,
            841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853,
            861, 862, 863, 864, 865, 866, 0
        ])
        self.assertEqual(
            self.kinase.key_res,
            [767, 775, 836, 838, 753, 770, 774, 864, 862, 863, 873, 814])

        # example 3: a kinase with multiple occupancy
        self.kinase = ki.main(
            ap.Namespace(chain='A', coord='pdb', feature='conf', pdb='1M17'))
        self.assertEqual(self.kinase.kinase_id, 406)
        self.assertEqual(self.kinase.name, 'EGFR')
        self.assertEqual(self.kinase.struct_id, 873)
        self.assertEqual(self.kinase.ligand, 'N/A')
        self.assertEqual(
            self.kinase.pocket_seq,
            'KVLGSGAFGTVYKVAIKELEILDEAYVMASVDPHVCRLLGIQLITQLMPFGCLLDYVREYLEDRRLVHRDLAARNVLVITDFGLA'
        )
        self.assertEqual(self.kinase.numbering, [
            692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704,
            718, 719, 720, 721, 722, 723, 734, 735, 736, 737, 738, 739, 740,
            741, 742, 743, 744, 745, 746, 748, 749, 750, 751, 752, 753, 754,
            755, 756, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773,
            774, 775, 776, 777, 778, 779, 780, 803, 804, 805, 806, 807, 808,
            809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821,
            829, 830, 831, 832, 833, 834, 835
        ])
        self.assertEqual(
            self.kinase.key_res,
            [735, 743, 804, 806, 721, 738, 742, 832, 830, 831, 841, 782])

    def test_features(self):
        # example 1: a kinase with no gap(s) in the binding pocket residues
        self.kinase = ki.main(
            ap.Namespace(chain='A', coord='pdb', feature='conf', pdb='3PP0'))
        self.assertEqual(
            round(
                np.asscalar(self.kinase.dihedrals[0][0]), 7),
            -2.0228872)  # the first dihedral value
        self.assertEqual(
            round(
                np.asscalar(self.kinase.distances[0][0]), 7),
            0.7770488)  # the first distance value
        
        # example 2: a kinase with gap(s) in the binding pocket residues
        self.kinase = ki.main(
            ap.Namespace(chain='A', coord='pdb', feature='conf', pdb='3RCD'))
        self.assertEqual(
            round(
                np.asscalar(self.kinase.dihedrals[0][0]), 7),
            -2.1615183)  # the first dihedral value
        self.assertEqual(
            round(
                np.asscalar(self.kinase.distances[0][0]), 7),
            1.0546854)  # the first distance value

        # example 3: a kinase with multiple occupancy
        self.kinase = ki.main(
            ap.Namespace(chain='A', coord='pdb', feature='conf', pdb='1M17'))
        self.assertEqual(
            round(
                np.asscalar(self.kinase.dihedrals[0][0]), 7),
            -2.4006705)  # the first dihedral value
        self.assertEqual(
            round(
                np.asscalar(self.kinase.distances[0][0]), 7),
            0.3558538)  # the first distance value

if __name__ == '__main__':
    unittest.main()

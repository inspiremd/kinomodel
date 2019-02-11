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
        # make sure the input command line is expected
        from kinomodel.features import featurize
        with self.assertRaises(ValueError):
            featurize(chain=0, coord='pdb', feature='interact', pdb='3PP0')

    def test_basics(self):
        # example 1: a kinase with no gap(s) in the binding pocket residues
        from kinomodel.features import featurize
        self.kinase = featurize(chain='A', coord='pdb', feature='interact', pdb='3PP0')
        self.assertEqual(self.kinase.kinase_id, 407)
        self.assertEqual(self.kinase.name, 'ErbB2')
        self.assertEqual(self.kinase.struct_id, 4820)
        self.assertEqual(self.kinase.ligand, '03Q')
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
            [])
        # example 2: a kinase with gap(s) in the binding pocket residues
        self.kinase = featurize(chain='A', coord='pdb', feature='interact', pdb='3RCD')
        self.assertEqual(self.kinase.kinase_id, 407)
        self.assertEqual(self.kinase.name, 'ErbB2')
        self.assertEqual(self.kinase.struct_id, 9325)
        self.assertEqual(self.kinase.ligand, '03P')
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
            [])

        # example 3: a kinase with multiple occupancy
        self.kinase = featurize(chain='A', coord='pdb', feature='interact', pdb='1M17')
        self.assertEqual(self.kinase.kinase_id, 406)
        self.assertEqual(self.kinase.name, 'EGFR')
        self.assertEqual(self.kinase.struct_id, 873)
        self.assertEqual(self.kinase.ligand, 'AQ4')
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
            [])

    def test_features(self):
        # example 1: a kinase with no gap(s) in the binding pocket residues
        self.kinase = featurize(chain='A', coord='pdb', feature='interact', pdb='3PP0')
        self.assertEqual(self.kinase.dihedrals,
            [])  # the first dihedral value
        self.assertEqual(self.kinase.distances,
            [])  # the first distance value
        self.assertEqual(
            round(
                np.asscalar(self.kinase.mean_dist), 7),
            1.3685706)  # mean ligand-pocket distance

        # example 2: a kinase with gap(s) in the binding pocket residues
        self.kinase = featurize(chain='A', coord='pdb', feature='interact', pdb='3RCD')
        self.assertEqual(self.kinase.dihedrals,
            [])  # the first dihedral value
        self.assertEqual(self.kinase.distances,
            [])  # the first distance value
        self.assertEqual(
            round(
                np.asscalar(self.kinase.mean_dist), 7),
            1.4236906)  # mean ligand-pocket distance

        # example 3: a kinase with multiple occupancy
        self.kinase = featurize(chain='A', coord='pdb', feature='interact', pdb='1M17')
        self.assertEqual(self.kinase.dihedrals,
            [])  # the first dihedral value
        self.assertEqual(self.kinase.distances,
            [])  # the first distance value
        self.assertEqual(
            round(
                np.asscalar(self.kinase.mean_dist), 7),
            1.6548363)  # mean ligand-pocket distance

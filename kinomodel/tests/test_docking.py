"""
Test docking
"""

# Import package, test suite, and other packages as needed
import unittest
import tempfile
from pkgutil import get_data
import os, sys

def get_data_filename(filename):
    """
    Return the path to specified data file

    .. todo :: Replace this with utils.get_data_filename

    """
    # this only works when kinomodel is installed (not conda-installable yet)
    #path = os.path.join(os.path.dirname(sys.modules['kinomodel'].__file__), 'data', filename)
    #JG (temperary)
    path = os.path.join('/home/guoj1/projects/INSPIRE/kinomodel/kinomodel/', 'data/', filename)
    return path

class DockingTestCase(unittest.TestCase):

    def test_hybrid_docking(self):
        "Test hybrid docking to Abl:nilotinib (PDBID:3CS9)"
        # absolute import (with kinomodel installed)
        # from kinomodel.docking import hybrid_docking
        #JG (temprary)
        from docking import hybrid
        receptor_path = get_data_filename('docking/3cs9.pdb')
        molecules_path = get_data_filename('fda-approved.smi')

        with tempfile.TemporaryDirectory() as docked_molecules_directory:
            docked_molecules_path = os.path.join(docked_molecules_directory, 'docked.mol2')
            hybrid.hybrid_docking(receptor_path, molecules_path, docked_molecules_path)

        # we also have the option to keep the docked files
        #docked_molecules_path = os.path.join(get_data_filename('docked.mol2'))
        #hybrid.hybrid_docking(receptor_path, molecules_path, docked_molecules_path)

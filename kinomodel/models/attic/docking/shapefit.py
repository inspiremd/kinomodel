from docking import create_receptor, load_receptor, dock_molecule, pose_molecule
import mdtraj as md
from openeye import oedocking, oechem
import openmoltools as moltools
import os.path
import os
import glob

pdb_directory = '../hauser-abl-benchmark/receptors'
ligand_directory = '../hauser-abl-benchmark/ligands'

pdbs = glob.glob(os.path.join(pdb_directory, '*'))
pdbs_list = []
for pdb in pdbs:
    pdbs_list.append(os.path.basename(pdb))

molecules_to_dock = ['C5=C(C1=CN=CC=C1)N=C(NC2=C(C=CC(=C2)NC(C3=CC=C(C=C3)CN4CCN(CC4)C)=O)C)N=C5 Imatinib',
                    'C1=C(SC(=N1)NC2=CC(=NC(=N2)C)N3CCN(CC3)CCO)C(=O)NC4=C(C=CC=C4Cl)C Dasatinib',
                    'Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)C(=O)Nc4cc(cc(c4)n5cc(nc5)C)C(F)(F)F Nilotinib',
                    'C1=C(C(=CC2=C1C(=C(C=N2)C#N)NC3=C(C=C(C(=C3)OC)Cl)Cl)OCCCN4CCN(CC4)C)OC Bosutinib',
                    'Cc1ccc(cc1C#Cc2cnc3n2nccc3)C(=O)Nc4ccc(c(c4)C(F)(F)F)CN5CCN(CC5)C Ponatinib',
                    'COc1cc2c(cc1OCCCN3CCOCC3)c(ncn2)Nc4ccc(c(c4)Cl)F Gefitinib',
                    'COCCOc1cc2c(cc1OCCOC)ncnc2Nc3cccc(c3)C#C Erlotinib',
                    'CNC(=O)c1ccccc1Sc2ccc3c(c2)[nH]nc3/C=C/c4ccccn4 Axitinib']

for pdb in pdbs_list:
    cocrystal = pdb.split('-')[1].split('_')[0]
    if cocrystal == 'gefitinib' or cocrystal ==  'erlotinib':
        cocrystal = pdb.split('-')[1].split('.')[0]
    cocrystal_path = os.path.join(ligand_directory, cocrystal + '.mol2')
    receptor_path = os.path.join(pdb_directory, pdb)

    # Load in protein
    input_mol_stream = oechem.oemolistream(receptor_path)
    protein_oemol = oechem.OEGraphMol()
    oechem.OEReadMolecule(input_mol_stream, protein_oemol)

    # Load in ligand
    input_mol_stream = oechem.oemolistream(cocrystal_path)
    ligand_oemol = oechem.OEGraphMol()
    oechem.OEReadMolecule(input_mol_stream, ligand_oemol)

    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein_oemol, ligand_oemol)

    for mol in molecules_to_dock:
        smiles, name = mol.split()
        posed_molecule_file_path = '%s-posed-%s.mol2' % (pdb.split('.')[0], name)
        posed_oemol = pose_molecule(receptor, smiles, n_poses=2)
        moltools.openeye.molecule_to_mol2(posed_oemol, posed_molecule_file_path, residue_name='MOL')

import simtk.openmm as mm
import simtk.unit as unit
import numpy as np


class Kinase(object):
    def __init__(self, pdb_chainid, kinase_id, name, struct_id, pocket_seq,
                 numbering, key_res, dihedrals, distances):
        '''
	Any kinase can be represented as a "kinase" object with the following parameters:
	----Parameters----
	pdb_chainid: [list of tuples] User-specified kinase info: (pdb code, chain index).
	kinase_id: [int] The standard ID of a kinase enforced by the KLIFS database.  This may be found in a kinaseId-lookUp table associated to this package.
	name: [str] The standard name of the kinase used by the KLIFS database.
        struct_id: [int] The ID associated with a specific chain in the pdb structure of a kinase.
	pocket_seq: [str] The 85 discontinuous residues (from multisequence alignment) that define the binding pocket of a kinase.
        numbering: [list of int] The residue indices of the 85 pocket residues specific to the structure.
	key_res: [list of int] A list of residue indices that are relevant to the collective variables.
        dihedrals: [list of floats] A list (one frame) or lists (multiple frames) of dihedrals relevant to kinase conformation.
        distances: [list of floats] A list (one frame) or lists (multiple frames) of intramolecular distances relevant to kinase conformation.
        '''
        self.pdb_chainid = pdb_chainid
        self.kinase_id = kinase_id
        self.name = name
        self.struct_id = struct_id
        self.pocket_seq = pocket_seq
        self.numbering = numbering
        self.key_res = key_res
        self.dihedrals = dihedrals
        self.distances = distances


class Interact(object):
    def __init__(self, pdb_chainid, kinase_id, name, struct_id, ligand,
                 pocket_seq, numbering, mean_dist):
        '''
        Any kinase can be represented as a "kinase" object with the following parameters:
        ----Parameters----
        pdb_chainid: [list of tuples] User-specified kinase info: (pdb code, chain index).
        kinase_id: [int] The standard ID of a kinase enforced by the KLIFS database.  This may be found in a kinaseId-lookUp table associated to this package.
        name: [str] The standard name of the kinase used by the KLIFS database.
        struct_id: [int] The ID associated with a specific chain in the pdb structure of a kinase.
        ligand: [str] The ligand name as it appears in the pdb file.
        pocket_seq: [str] The 85 discontinuous residues (from multisequence alignment) that define the binding pocket of a kinase.
        numbering: [list of int] The residue indices of the 85 pocket residues specific to the structure.
        mean_dist: [float] A float (one frame) or a list of floats (multiple frames), which is the mean pairwise distance between ligand heavy atoms and the CAs of the 85 pocket residues.  
        '''
        self.pdb_chainid = pdb_chainid
        self.kinase_id = kinase_id
        self.name = name
        self.struct_id = struct_id
        self.ligand = ligand
        self.pocket_seq = pocket_seq
        self.numbering = numbering
        self.mean_dist = mean_dist

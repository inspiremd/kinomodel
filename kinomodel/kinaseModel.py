import simtk.openmm as mm
import simtk.unit as unit
import numpy as np

class kinase(object):
    def __init__(self, pdbChain, kinaseID, name, structID, pocketSeq, numbering, keyRes):
        '''
	Any kinase can be represented as a "kinase" object with the following parameters:
	----Parameters----
	pdbChain: [list of tuples] User-specified kinase info: (pdb code, chain index).
	kinaseID: [int] The standard ID of a kinase enforced by the KLIFS database.  This may be found in a kinaseId-lookUp table associated to this package.
	name: [str] The standard name of the kinase used by the KLIFS database.
        structID: [int] The ID associated with a specific chain in the pdb structure of a kinase.
	pocketSeq: [str] The 85 discontinuous residues (from multisequence alignment) that define the binding pocket of a kinase.
        numbering: [list of int] The residue indices of the 85 pocket residues specific to the structure.
	keyRes: [list of int] A list of residue indices that are relevant to the collective variables.
        '''
        self.pdbChain = pdbChain
        self.kinaseID = kinaseID
        self.name = name
        self.structID = structID
        self.pocketSeq = pocketSeq
        self.numbering = numbering
        self.keyRes = keyRes


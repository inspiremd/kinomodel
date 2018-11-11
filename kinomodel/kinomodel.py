"""
kinomodel.py
This is a tool  for kinase conformational change/ligand binding pose modeling within the entire Kinome.

Handles the primary functions
"""
import sys
from kinase_model import *
import requests
import ast
import urllib.request

def main(args=None):
    """
    When given a PDB code plus a chain index, this script collect information including (1) the kinase ID, (2) the standard name, (3) the KLIFS-defined structure ID, (4) the KLIFS-defined sequence of 85 binding pocket residues (http://klifs.vu-compmedchem.nl/index.php), (5) the indices of the 85 residues in the corresponding structure, and (6) residue indices involved in collective variables for kinase conformational changes and ligand interactions.

    """
    # get the user input (a PDB code and a chain index) 
    input_info = input('Please input kinase info (PDB code, chain index) e.g. (3pp0, A): ').replace(' ','').split(',')
    pdb_chainid = tuple(input_info)

    # get information of the query kinase from the KLIFS database and gives values of kinase_id, name and pocket_seq (numbering)
    url="http://klifs.vu-compmedchem.nl/api/structures_pdb_list?pdb-codes="+str(pdb_chainid[0]) # form the query command
    clean = requests.get(url).text.replace('true','True').replace('false','False') # clean up the info from KLIFS
    for structure in ast.literal_eval(clean): # each pdb code corresponds to multiple structures
        if structure['chain'] == str(pdb_chainid[1]): # find the specific chain
            kinase_id = int(structure['kinase_ID'])
            name = str(structure['kinase'])
            pocket_seq = str(structure['pocket'])
            struct_id = int(structure['structure_ID'])

    # Get the numbering of the 85 pocket residues
    cmd = "http://klifs.vu-compmedchem.nl/details.php?structure_id="+str(struct_id)
    info = urllib.request.urlopen(cmd)
    for line_number, line in enumerate(info):
        line = line.decode()
        if 'pocketResidues=[' in line:
            numbering = ast.literal_eval((line[line.find('=')+1:line.find(';')]))
   
    # check if there is gaps/missing residues among the pocket residues. If so, enforce their indices as 0 and avoid using them to compute collective variables.
    for i in range(len(numbering)):
        if numbering[i] == -1:
            print ("Warning: There is a gap/missing residue at position: "+str(i+1)+". Its index will be enforced as 0 and it will not be used to compute collective variables.")
            numbering[i] = 0

    # define indices of the residues relevant to a list of 12 collective variables relevant to kinase conformational changes. These variables include: angle between aC and aE helices, the key K-E salt bridge, DFG-Phe conformation (two distances), X-DFG-Phi, X-DFG-Psi, DFG-Asp-Phi, DFG-Asp-Psi, DFG-Phe-Phi, DFG-Phe-Psi, DFG-Phe-Chi1, and the FRET L-S distance. All features are under the current numbering of the structure provided.
    key_res = []
    # angle between aC and aE helices
    key_res.append(numbering[20]) # residue 21 (res1 in aC)
    key_res.append(numbering[28]) # res29 (res2 in aC)
    key_res.append(numbering[60]) # res61 (res1 in aE)
    key_res.append(numbering[62]) # res63 (res2 in aE)

    # key salt bridge
    key_res.append(numbering[16]) # res17 (K in beta3) 
    key_res.append(numbering[23]) # res24 (E in aC)

    # DFG conformation and Phe conformation
    key_res.append(numbering[27]) # res28 (ExxxX)
    key_res.append(numbering[81]) # res82 (DFG-Phe)

    # X-DFG Phi/Psi
    key_res.append(numbering[79]) # res80 (X-DFG)

    # DFG-Asp conformation
    key_res.append(numbering[80]) # res81 (DFG-Asp)

    # FRET distance
    key_res.append(numbering[84]+6 if numbering[84] else 0) # not in the list of 85 (equivalent to Aura"S284"), only infer if the reference is non-zero

    key_res.append(numbering[58]+2 if numbering[58] else 0) # not in the list of 85 (equivalent to Aura"L225"), only infer if the reference is non-zero

    # print out kinase information
    print ("---------------------Results----------------------")
    print ("Kinase ID: "+str(kinase_id))
    print ("Kinase name: "+str(name))
    print ("Pocket residues: "+str(pocket_seq))
    print ("Structure ID: "+str(struct_id))
    print ("Numbering: "+str(numbering))
    print ("Residues involved in collective variables: "+str(key_res))

    # Pupolate a kinase object
    my_kinase = Kinase(pdb_chainid, kinase_id, name, struct_id, pocket_seq, numbering, key_res)

    return my_kinase


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    main()

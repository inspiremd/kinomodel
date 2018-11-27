'''
kinomodel.py
This is a tool for kinase conformation/ligand binding pose modeling through the entire Kinome.
Handles the primary functions

'''
import sys
from kinase_model import * 
from protein_features import *
from interact_features import *

def main():

# get parameters from the protein or interact featurization modules to populate a Kinase or Interact object based on user intention
    input_info = input('Please specify whether you would like to compute collective variables related to protein conformation or protein-ligand interaction (by typing "conf" or "interact"): ')

    if input_info == "conf":
        (pdb_chainid, kinase_id, name, struct_id, pocket_seq, numbering, key_res, dihedrals, distances) = protein() 
        my_kinase = Kinase(pdb_chainid, kinase_id, name, struct_id, pocket_seq,
                       numbering, key_res, dihedrals, distances)
    elif input_info == "interact":
        (pdb_chainid, kinase_id, name, struct_id, ligand, pocket_seq, numbering, mean_dist) = interact() 
        my_interact = Interact(pdb_chainid, kinase_id, name, struct_id, ligand, pocket_seq, numbering, mean_dist)
        
    return

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    main()

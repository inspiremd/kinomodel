'''
kinomodel.py
This is a tool for kinase conformation/ligand binding pose modeling through the entire Kinome.
Handles the primary functions

'''
import sys
from kinase_model import *
import protein_features as pf
import interact_features as inf
# clean traceback msgs
sys.tracebacklimit = 0


def main():

    # get parameters from the protein or interact featurization modules to populate a Kinase or Interact object based on user intention
    input_info = input(
        'Please specify whether you would like to compute collective variables related to protein conformation or protein-ligand interaction (by typing "conf" or "interact"): '
    )

    # make sure the input format is expect
    if str(input_info) != 'conf' and str(input_info) != 'interact':
        raise ValueError(
            'You must choose from kinase conformation ("conf") and ligand interaction ("interact").'
        )

    if input_info == "conf":
        (pdb_chainid, kinase_id, name, struct_id, pocket_seq, numbering,
         key_res) = pf.basics()
        (dihedrals, distances) = pf.features(pdb_chainid, numbering)
        my_kinase = Kinase(pdb_chainid, kinase_id, name, struct_id, pocket_seq,
                           numbering, key_res, dihedrals, distances)
    elif input_info == "interact":
        (pdb_chainid, kinase_id, name, struct_id, ligand, pocket_seq,
         numbering) = inf.basics()
        mean_dist = inf.features(pdb_chainid, ligand, numbering)
        my_interact = Interact(pdb_chainid, kinase_id, name, struct_id, ligand,
                               pocket_seq, numbering, mean_dist)

    return


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    main()

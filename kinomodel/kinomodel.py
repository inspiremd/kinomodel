'''
kinomodel.py
This is a tool for kinase conformational/ligand binding featurization through the entire Kinome.
Handles the primary functions

'''
import sys
import argparse
import kinase_model

# main function
def main(args):
    """
    This is the main function for kinomodel. It takes the PDB code, chain id and certain coordinates of 
    a kinase from a command line and returns its basic information and structural and/or interation features.

    Parameters
    ----------
    args: a Namespace object from argparse
        Information from parsing the command line
        e.g. Namespace(chain='A', coord='pdb', feature='conf', pdb='3PP0')
    Returns
    -------
    my_kinase: a Kinase object from kinase_model
        Contains basic information of a kinase including PDB code, chain id, kinase id, name, structure id, 
        ligand name, the 85 pocket residues and their numbering, key dihedrals and distances that define 
        the conformational state of a kinase, the mean distance between ligand heavy atoms and the pocket residues etc.

    """
    if args.feature == "conf":
        import protein_features as pf
        (kinase_id, name, struct_id, pocket_seq, numbering,
         key_res) = pf.basics(args.pdb, args.chain)
        ligand = 'N/A'
        mean_dist = 0
        (dihedrals, distances) = pf.features(args.pdb, args.chain, args.coord,
                                             numbering)
        my_kinase = kinase_model.Kinase(
            args.pdb, args.chain, kinase_id, name, struct_id, ligand,
            pocket_seq, numbering, key_res, dihedrals, distances, mean_dist)
    elif args.feature == "interact":
        import interact_features as inf
        (kinase_id, name, struct_id, ligand, pocket_seq, numbering) = inf.basics(args.pdb, args.chain)
        key_res = []
        dihedrals = []
        distances = []
        mean_dist = inf.features(args.pdb, args.chain, args.coord, ligand,
                                 numbering)
        my_kinase = kinase_model.Kinase(
            args.pdb, args.chain, kinase_id, name, struct_id, ligand,
            pocket_seq, numbering, key_res, dihedrals, distances, mean_dist)
    elif args.feature == "both":
        import protein_features as pf
        import interact_features as inf
        (kinase_id, name, struct_id, ligand, pocket_seq,
         numbering) = inf.basics(args.pdb, args.chain)
        (kinase_id, name, struct_id, pocket_seq, numbering,
         key_res) = pf.basics(args.pdb, args.chain)
        (dihedrals, distances) = pf.features(args.pdb, args.chain, args.coord,
                                             numbering)
        mean_dist = inf.features(args.pdb, args.chain, args.coord, ligand,
                                 numbering)
        my_kinase = kinase_model.Kinase(
            args.pdb, args.chain, kinase_id, name, struct_id, ligand,
            pocket_seq, numbering, key_res, dihedrals, distances, mean_dist)
    return my_kinase


if __name__ == "__main__":
    # parsing command-line arguments
    parser = argparse.ArgumentParser(
        prog='kinomodel',
        description='Featurize kinase structures',
    )
    parser.add_argument(
        '--pdb',
        required=True,
        action='store',
        nargs='?',
        type=str,
        help='the PDB code of the structure to featurize')
    parser.add_argument(
        '--chain',
        required=True,
        action='store',
        nargs='?',
        type=str,
        help='the chain index of the structure to featurize')
    parser.add_argument(
        '--feature',
        required=True,
        action='store',
        choices=['conf', 'interact', 'both'],
        nargs='?',
        help=
        'compute collective variables related to protein conformation, protein-ligand interaction, or both'
    )
    parser.add_argument(
        '--coord',
        required=True,
        action='store',
        choices=['pdb','dcd']
        nargs='?',
        type=str,
        help=
        'the coordinates (a pdb file) or trajectories (e.g. a dcd file with associated topology info) to featurize.' 
        'Default is the pdb coordinates under the given PDB code.'
    )
    args = parser.parse_args()
    # Do something if this file is invoked on its own
    main(args)

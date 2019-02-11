"""
kinomodel.py
This is a tool for kinase conformational/ligand binding featurization through the entire Kinome.
Handles the primary functions

"""

def featurize_cli(**kwargs):
    """Checks whether the user has provided any arguments directly to
    the main function through the kwargs keyword, and if not, attempt
    to find them from the command line arguments using argparse.

    Notes
    -----
    This method requires the pdb, chain, feature and coord args to be either
    all present in the kwargs dictionary, or given on the command line.

    Parameters
    ----------
    kwargs: dict of str, values, optional
        A dictionary of passed arguments (which may be empty), where the key is
        the name of the argument, and the value is the argument value.

    Returns
    -------
    argparse.Namespace

        A namespace containing the arguments which either were directly passed
        in as kwargs, or this which were taken as input from the comman line.

        e.g. Namespace(chain='A', coord='pdb', feature='conf', pdb='3PP0')
    """

    arguments = None

    if not kwargs:
        import argparse

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
            help='compute collective variables related to protein conformation, protein-ligand interaction, or both'
        )
        parser.add_argument(
            '--coord',
            required=True,
            action='store',
            choices=['pdb', 'dcd'],
            nargs='?',
            type=str,
            help='the coordinates (a pdb file) or trajectories (e.g. a dcd file with associated topology info) to '
                 'featurize. Default is the pdb coordinates under the given PDB code.'
        )

        arguments = parser.parse_args()

    else:

        # Check to make sure that all of the required args are present
        # in the kwargs dictionary.

        # You might want to add your own validation logic here.
        assert 'pdb' in kwargs
        assert 'chain' in kwargs
        assert 'feature' in kwargs
        assert 'coord' in kwargs

        arguments = argparse.Namespace(**kwargs)

    return arguments

# main function
def dispatch_featurize(**kwargs):
    """
    Compute structural and/or interaction features for a given PDB structure.

    Parameters
    ----------
    args: a Namespace object from argparse
        Information from parsing the command line
        e.g. Namespace(chain='A', coord='pdb', feature='conf', pdb='3PP0')

    Returns
    -------
    my_kinase: kinomodel.models.Kinase
        Contains basic information of a kinase including PDB code, chain id, kinase id, name, structure id,
        ligand name, the 85 pocket residues and their numbering, key dihedrals and distances that define
        the conformational state of a kinase, the mean distance between ligand heavy atoms and the pocket residues etc.

    """
    from kinomodel.models import Kinase
    from kinomodel.features import protein as pf
    from kinomodel.features import interaction as inf

    args = _parse_arguments(**kwargs)

    my_kinase = None

    if args.feature == "conf":
        (kinase_id, name, struct_id, pocket_seq, numbering, key_res) = pf.basics(args.pdb, args.chain)
        ligand = 'N/A'
        mean_dist = 0
        (dihedrals, distances) = pf.features(args.pdb, args.chain, args.coord, numbering)

    elif args.feature == "interact":
        (kinase_id, name, struct_id, ligand, pocket_seq, numbering) = inf.basics(args.pdb, args.chain)
        key_res = []
        dihedrals = []
        distances = []
        mean_dist = inf.features(args.pdb, args.chain, args.coord, ligand, numbering)

    elif args.feature == "both":
        (kinase_id, name, struct_id, ligand, pocket_seq, numbering) = inf.basics(args.pdb, args.chain)
        (kinase_id, name, struct_id, pocket_seq, numbering, key_res) = pf.basics(args.pdb, args.chain)
        (dihedrals, distances) = pf.features(args.pdb, args.chain, args.coord, numbering)
        mean_dist = inf.features(args.pdb, args.chain, args.coord, ligand, numbering)

    else:
        raise Exception("Unknown feature '{}'".format(args.feature))


    my_kinase = Kinase(
        args.pdb, args.chain, kinase_id, name, struct_id, ligand,
        pocket_seq, numbering, key_res, dihedrals, distances, mean_dist)

    return my_kinase

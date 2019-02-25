"""
Utilities for featurizing kinase-ligand interactions

"""

# Setup general logging (guarantee output/error message in case of interruption)
# TODO: Can we log to the terminal instead?
import logging
logger = logging.getLogger(__name__)
logging.root.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("urllib3").setLevel(logging.WARNING)

def compute_simple_interaction_features(pdbid, chainid, coordfile, ligand_name, resids):
    """
    This function takes the PDB code, chain id, certain coordinates, ligand name and the numbering of
    pocket residues of a kinase from a command line and returns its structural features.

    Parameters
    ----------
    pdbid: str
        The PDB code of the query kinase.
    chainid: str
        The chain index of the query kinase.
    coordfile: str
        Specifies the source of coordinates ('pdb' or 'dcd')
    ligand_name: str
        Specifies the ligand name of the complex.
    resids: list of int
        Protein residue indices to use in computing simple interaction features.

    Returns
    -------
    mean_dist: float
            A float (one frame) or a list of floats (multiple frames), which is the mean pairwise distance
            between ligand heavy atoms and the CAs of the 85 pocket residues.

    .. todo :: Instead of a PDB file or a trj/dcd/h5, accept an MDTraj.Trajectory---this will be much more flexible.

    .. todo :: Use kwargs with sensible defaults instead of relying only on positional arguments.

    """
    import tempfile
    import os
    import mdtraj as md
    import numpy as np

    pdb_file = None

    # A safer way to download files as wget may not exist on systems such MacOS
    # TODO: Since we retrieve the PDB file in multiple pieces of code, let's refactor this into one utility function
    # to avoid code duplication.
    import urllib
    with urllib.request.urlopen('http://www.pdb.org/pdb/files/{}.pdb'.format(pdbid)) as response:
        pdb_file = response.read()

    with tempfile.TemporaryDirectory() as pdb_directory:
        pdb = os.path.join(pdb_directory,'{}.pdb'.format(pdbid))
        with open(pdb, 'w') as file:
            file.write(pdb_file.decode())
            # load traj before the temp pdb file was removed
            if coordfile == 'pdb':
                traj = md.load(pdb)
            # get topology info from the structure
            topology = md.load(pdb).topology

    table, bonds = topology.to_dataframe()
    atoms = table.values
    # translate a letter chain id into a number index (A->0, B->1 etc)
    # TODO: This may not be robust, since chains aren't always in sequence from A to Z
    chain_index = ord(str(chainid).lower()) - 97

    #np.set_printoptions(threshold=np.nan)
    #print (atoms)
    # get the array of atom indices for the calculation of:
    #       * mean of pairwise distances between each ligand atom and CA of 85 binding pocket
    #residues (an (85*n*2) array where n = # of ligand heavy atoms (usually <= 100)
    #and each row contains indices of the two atoms for each distance)
    dis = np.zeros(shape=(8500, 2), dtype=int, order='C')

    # parse the topology info
    pocket_atm = []
    chain_num = 0
    atm_count = 0
    count = 0
    for line in atoms:
        # find the key binding pocket atoms in the protein
        # for CA of the pocket residues in the specified protein chain
        if (line[5] == chain_index) and (line[3] in resids) and (line[1] == 'CA'):
            pocket_atm.append(int(line[0]))
        # start of ligand
        if line[4] == ligand_name:
            # check whether the pocket_atm list is complete
            if (len(pocket_atm) < 85) and (0 in resids):
                pocket_atm.insert(resids.index(0), 0)
            # for the specified chain and ligand heavy atoms
            if line[5] == int(chain_num + 1) and 'H' not in line[1]:
                for i in range(85):
                    dis[count * 85 + i][0] = atm_count
                    dis[count * 85 + i][1] = pocket_atm[i]
                count += 1
        else:
            chain_num = line[5]
        atm_count += 1

    # clean array and remove empty lines
    import sys
    np.set_printoptions(threshold=sys.maxsize)
    dis = dis[~np.all(dis == 0, axis=1)]
    # check if there is any missing coordinates;
    # if so, skip distance calculation for those residues
    del_lst = []
    # find out lines with 0 at the protein residue position
    for i in range(len(dis)):
        if dis[i][1] == 0:
            dis[i][0] = 0
    for i in range(len(dis)):
        if dis[i][0] and dis[i][1]:
        # the atom indices fed to mdtraj should be 0-based
            dis[i] -= 1
    #if check_dis == 0:
        #logging.info(
        #    "There is no missing coordinates.  All distances will be computed."
        #)
    #else:
        #logging.info(
        #    "Some of the pairwise distances will not be calculated due to missing coordinates."
        #)

    # calculate the distances for the user-specifed structure (a static structure or an MD trajectory)
    if coordfile == 'dcd':
        traj = md.load(str(pdbid) + '.dcd',top = str(pdbid) + '_fixed_solvated.pdb')
    mean_dist = []
    for frame in md.compute_distances(traj, dis):
        mean_dist.append(np.mean(frame))

    #logging.info(
    #    "The mean distance between ligand heavy atoms and CAs of the 85 binding pocket residues for PDB# "
    #    + str(pdbid) + ", chain " + str(chainid) + " is: " +
    #    str(mean_dist))

    # clean up
    # TODO: This is dangerous! Instead, rely on using the "with tempfile.TemporaryDirectory" context manager idiom to create and clean up temporary directories
    #rm_file = 'rm ./' + str(pdb) + '.pdb*'
    #import subprocess
    #subprocess.call(rm_file, shell=True)
    del traj, dis

    return mean_dist

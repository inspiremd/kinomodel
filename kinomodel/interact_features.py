'''
interact_features.py
This is a tool to featurize kinase-ligand interaction through the entire Kinome.

'''
import logging
import sys
import requests
import ast
import urllib.request
import simtk.openmm as mm
import simtk.unit as unit
import numpy as np
import mdtraj as md
import subprocess

## Setup general logging (guarantee output/error message in case of interruption)
logger = logging.getLogger(__name__)
logging.root.setLevel(logging.DEBUG)
logging.basicConfig(level=logging.DEBUG, filename="kinomodel_interaction.log", filemode="a+", format="%(message)s")
logging.getLogger("urllib3").setLevel(logging.WARNING)

def basics(pdb, chain):
    """
    This function takes the PDB code and chain id of a kinase from a command line and returns its basic information.
    
    Parameters
    ----------
    pdb: str
        The PDB code of the inquiry kinase.
    chain: str
        The chain index of the inquiry kinase.

    Returns
    -------
    kinase_id: int
        The standard ID of a kinase enforced by the KLIFS database.
    name: str
        The standard name of the kinase used by the KLIFS database.
    pocket_seq: str
        The 85 discontinuous residues (from multisequence alignment) that define the binding pocket of a kinase.
    struct_id: int
        The ID associated with a specific chain in the pdb structure of a kinase.
    ligand: str
        The ligand name as it appears in the pdb file.
    numbering: list of int
        The residue indices of the 85 pocket residues specific to the structure.

    """

    # get information of the query kinase from the KLIFS database and gives values 
    # of kinase_id, name and pocket_seq (numbering)
    url = "http://klifs.vu-compmedchem.nl/api/structures_pdb_list?pdb-codes=" + str(pdb)

    # check to make to sure the search returns valid info
    # if return is empty
    if len(requests.get(url).text) == 0:
        raise ValueError("No matched pdb structure found in KLIFS. Please make sure to provide a valid PDB code.")
    else:
        # clean up the info from KLIFS
        clean = requests.get(url).text.replace('true', 'True').replace('false', 'False')

    # each pdb code corresponds to multiple structures
    found = 0
    for structure in ast.literal_eval(clean):
        # find the specific chain
        if structure['chain'] == str(chain):
            kinase_id = int(structure['kinase_ID'])
            name = str(structure['kinase'])
            pocket_seq = str(structure['pocket'])
            struct_id = int(structure['structure_ID'])
            # make sure the specified structure is not an apo structure
            if structure['ligand'] != 0:
                ligand = str(structure['ligand'])
            else:
                raise ValueError(
                    "The PDB code you provide corresponds to an apo protein (no ligand) so the receptor-ligand interaction cannot be computed. Please double check."
                )
            found = 1
    if not found:
        raise ValueError("No matched chain found. Please make sure you provide a capital letter (A, B, C, ...) as a chain ID.")
    # Get the numbering of the 85 pocket residues
    cmd = "http://klifs.vu-compmedchem.nl/details.php?structure_id=" + str(
        struct_id)
    preload = urllib.request.urlopen(cmd)
    info = urllib.request.urlopen(cmd)
    for line_number, line in enumerate(info):
        line = line.decode()
        if 'pocketResidues=[' in line:
            numbering = ast.literal_eval(
                (line[line.find('=') + 1:line.find(';')]))
    # check if there is gaps/missing residues among the pocket residues.
    # If so, enforce their indices as 0 and avoid using them to compute collective variables.
    for i in range(len(numbering)):
        if numbering[i] == -1:
            logging.info(
                "Warning: There is a gap/missing residue at position: " +
                str(i + 1) +
                ". Its index will be enforced as 0 and it will not be used to compute collective variables."
            )
            numbering[i] = 0

    # print out kinase information
    logging.info("Kinase ID: " + str(kinase_id))
    logging.info("Kinase name: " + str(name))
    logging.info("Pocket residues: " + str(pocket_seq))
    logging.info("Structure ID: " + str(struct_id))
    logging.info("Ligand name: " + str(ligand))
    logging.info("Numbering of the 85 pocket residues: " + str(numbering))

    return kinase_id, name, struct_id, ligand, pocket_seq, numbering


def features(pdb, chain, coord, ligand, numbering):
    """
    This function takes the PDB code, chain id, certain coordinates, ligand name and the numbering of
    pocket residues of a kinase from a command line and returns its structural features.
    
    Parameters
    ----------
    pdb: str
        The PDB code of the inquiry kinase.
    chain: str
        The chain index of the inquiry kinase.
    coord: str
        Specifies the file constaining the kinase coordinates (either a pdb file or a trajectory, i.e. trj, dcd, h5)
    ligand: str
        Specifies the ligand name of the complex.
    numbering: list of int
        The residue indices of the 85 pocket residues specific to the structure.    

    Returns
    -------
    mean_dist: float
            A float (one frame) or a list of floats (multiple frames), which is the mean pairwise distance 
            between ligand heavy atoms and the CAs of the 85 pocket residues.
    
    """
    # download the pdb structure
    cmd = 'wget -q http://www.pdb.org/pdb/files/' + str(
        pdb) + '.pdb'
    subprocess.call(cmd, shell=True)

    # get topology info from the structure
    topology = md.load(str(pdb) + '.pdb').topology
    table, bonds = topology.to_dataframe()
    atoms = table.values
    # translate a letter chain id into a number index (A->0, B->1 etc)
    chain_index = ord(str(chain).lower()) - 97 

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
        if line[5] == chain_index and line[3] in numbering and line[1] == 'CA':
            pocket_atm.append(int(line[0]))
        # start of ligand
        if line[4] == ligand:
            # check whether the pocket_atm list is complete
            if len(pocket_atm) < 85 and 0 in numbering:
                pocket_atm.insert(numbering.index(0), 0)
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
    np.set_printoptions(threshold=np.nan)
    dis = dis[~np.all(dis == 0, axis=1)]
    # check if there is any missing coordinates;
    # if so, skip distance calculation for those residues
    check_flag = 1
    del_lst = []
    # find out lines with 0 at the protein residue position
    for i in range(len(dis)):
        if dis[i][1] == 0:
            del_lst.append(i)

    if del_lst:
        # delete them all at once
        dis = np.delete(dis, (del_lst), axis=0)
        check_flag = 0
    for i in range(len(dis)):
        # the atom indices fed to mdtraj should be 0-based
        dis[i] -= 1
    if check_flag:
        logging.info(
            "There is no missing coordinates.  All distances will be computed."
        )
    else:
        logging.info(
            "Some of the pairwise distances will not be calculated due to missing coordinates."
        )

    # calculate the distances for the user-specifed structure (a static structure or an MD trajectory)
    if coord == 'pdb':
        traj = md.load(str(pdb) + '.pdb')
    elif coord == 'dcd':
        traj = md.load(str(pdb) + '.dcd',top = str(pdb) + '_fixed_solvated.pdb')
    mean_dist = []
    for frame in md.compute_distances(traj, dis):
        mean_dist.append(np.mean(frame))

    logging.info(
        "The mean distance between ligand heavy atoms and CAs of the 85 binding pocket residues for PDB# "
        + str(pdb) + ", chain " + str(chain) + " is: " +
        str(mean_dist))

    # clean up
    rm_file = 'rm ./' + str(pdb) + '.pdb*'
    subprocess.call(rm_file, shell=True)
    del traj, dis

    return mean_dist


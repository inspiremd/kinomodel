'''
interact_features.py
This is a tool for featurization of ligand interaction to kinases through the entire Kinome.

'''

import sys
import requests
import ast
import urllib.request
import simtk.openmm as mm
import simtk.unit as unit
import numpy as np
import mdtraj as md
import subprocess


def interact(args=None):
    '''
    When given a PDB code plus a chain index, this script computes distances as well as features such as intermolecular H-bonding that together define protein-ligand interaction.

    '''
    # get the the relevant PDB code and chain index from user input
    input_info = input(
        'Please input the PDB code and chain index of the structure you initialized or will initialize your simulation with) e.g. (3pp0, A): '
    ).replace(' ', '').split(',')
    pdb_chainid = tuple(input_info)
    #pdb_chainid = ('3rcd','A')

    # get information of the query kinase from the KLIFS database and gives values of kinase_id, name and pocket_seq (numbering)
    url = "http://klifs.vu-compmedchem.nl/api/structures_pdb_list?pdb-codes=" + str(
        pdb_chainid[0])  # form the query command
    clean = requests.get(url).text.replace('true', 'True').replace(
        'false', 'False')  # clean up the info from KLIFS
    for structure in ast.literal_eval(
            clean):  # each pdb code corresponds to multiple structures
        if structure['chain'] == str(
                pdb_chainid[1]):  # find the specific chain
            kinase_id = int(structure['kinase_ID'])
            name = str(structure['kinase'])
            pocket_seq = str(structure['pocket'])
            struct_id = int(structure['structure_ID'])
            if structure[
                    'ligand'] != 0:  # make sure the specified structure is not an apo structure
                ligand = str(structure['ligand'])
            else:
                raise ValueError(
                    "The PDB code you provide corresponds to an apo protein (no ligand) so the receptor-ligand interaction cannot be computed. Please double check."
                )
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
    # check if there is gaps/missing residues among the pocket residues. If so, enforce their indices as 0 and avoid using them to compute collective variables.
    for i in range(len(numbering)):
        if numbering[i] == -1:
            print(
                "Warning: There is a gap/missing residue at position: " +
                str(i + 1) +
                ". Its index will be enforced as 0 and it will not be used to compute collective variables."
            )
            numbering[i] = 0

    # print out kinase information
    print("---------------------Results----------------------")
    print("Kinase ID: " + str(kinase_id))
    print("Kinase name: " + str(name))
    print("Pocket residues: " + str(pocket_seq))
    print("Structure ID: " + str(struct_id))
    print("Ligand name: " + str(ligand))
    print("Numbering of the 85 pocket residues: " + str(numbering))
    '''    
    Download the pdb structure corresponding to the given PDB code and chain index, where the atom indices of the relevant atoms will be inferred and used to calculate dihedrals and distances as collective variables. 
    '''
    # download the pdb structure
    cmd = 'wget http://www.pdb.org/pdb/files/' + str(pdb_chainid[0]) + '.pdb'
    subprocess.call(cmd, shell=True)

    # get topology info from the structure
    topology = md.load(str(pdb_chainid[0]) + '.pdb').topology
    table, bonds = topology.to_dataframe()
    atoms = table.values
    chain_index = ord(str(pdb_chainid[1]).lower(
    )) - 97  # translate a letter chain id into a number index (A->0, B->1 etc)

    #np.set_printoptions(threshold=np.nan)
    #print (atoms)
    # get the array of atom indices for the calculation of:
    #       * mean of pairwise distances between each ligand atom and CA of 85 binding pocket residues (an (85*n*2) array where n = # of ligand heavy atoms (usually <= 100) and each row contains indices of the two atoms for each distance)
    dis = np.zeros(shape=(8500, 2), dtype=int, order='C')

    # parse the topology info
    pocket_atm = []
    chain_num = 0
    count = 0
    for line in atoms:
        # find the key binding pocket atoms in the protein
        if line[5] == chain_index and line[3] in numbering and line[
                1] == 'CA':  # for CA of the pocket residues in the specified protein chain
            pocket_atm.append(int(line[0]))
        if line[4] == ligand:  # start of ligand
            # check whether the pocket_atm list is complete
            if len(pocket_atm) < 85 and 0 in numbering:
                pocket_atm.insert(numbering.index(0), 0)
            # for the specified chain and ligand heavy atoms
            if line[5] == int(chain_num + 1) and 'H' not in line[1]:
                for i in range(85):
                    dis[count * 85 + i][0] = line[0]
                    dis[count * 85 + i][1] = pocket_atm[i]
                count += 1
        else:
            chain_num = line[5]

    # clean array and remove empty lines
    np.set_printoptions(threshold=np.nan)
    dis = dis[~np.all(dis == 0, axis=1)]
    # check if there is any missing coordinates; if so, skip distance calculation for those residues
    check_flag = 1
    del_lst = []
    for i in range(
            len(dis)):  # find out lines with 0 at the protein residue position
        if dis[i][1] == 0:
            del_lst.append(i)

    if del_lst:  # delete them all at once
        dis = np.delete(dis, (del_lst), axis=0)
        check_flag = 0
    if check_flag:
        print(
            "There is no missing coordinates.  All distances will be computed."
        )
    else:
        print(
            "Some of the pairwise distances will not be calculated due to missing coordinates."
        )

    # calculate the distances for the user-specifed structure (a static structure or an MD trajectory)
    user_input = str(
        input(
            'Please specify the (full path to) the trajectory (and if necessary, also the topology file) to analyze in the following format (trajectory) or (trajectory,topology); otherwise type "pdb" if you want the current pdb structure to be analyzed): '
        )).replace(' ', '')
    if ',' in user_input:  # if both trajectory and topology files are input
        user_input = user_input.split(',')
        traj_top = tuple(user_input)
        traj = md.load(str(traj_top[0]), top=str(traj_top[1]))
    elif user_input == 'pdb':
        traj = md.load(str(pdb_chainid[0]) + '.pdb')
    mean_dist = np.mean(md.compute_distances(traj, dis)[0])

    print(
        "The mean distance between ligand heavy atoms and CAs of the 85 binding pocket residues for PDB# "
        + str(pdb_chainid[0]) + ", chain " + str(pdb_chainid[1]) + " is: " +
        str(mean_dist))

    # clean up
    rm_file = 'rm ./' + str(pdb_chainid[0]) + '.*'
    subprocess.call(rm_file, shell=True)
    del traj, dis

    return pdb_chainid, kinase_id, name, struct_id, ligand, pocket_seq, numbering, mean_dist

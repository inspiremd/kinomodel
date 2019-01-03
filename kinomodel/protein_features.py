"""
protein_features.py
This is a tool for featurization of kinase conformational changes through the entire Kinome.

"""
import sys
import requests
import ast
import urllib.request
import simtk.openmm as mm
import simtk.unit as unit
import numpy as np
import mdtraj as md
import subprocess

# clean traceback msgs
sys.tracebacklimit = 0


def basics(args=None):
    '''
    When given a PDB code plus a chain index, this script collects information including (1) the kinase ID, (2) the standard name, (3) the KLIFS-defined structure ID, (4) the KLIFS-defined sequence of 85 binding pocket residues (http://klifs.vu-compmedchem.nl/index.php), (5) the indices of the 85 residues in the corresponding structure, (6) residue indices involved in collective variables for kinase conformational changes and ligand interactions, (7) key dihedral angles and (8) intramolecular distances that are relevant to kinase conformation.

    '''
    # get the the relevant PDB code and chain index from user input
    input_info = str(
        input(
            'Please input the PDB code and chain index of the structure you initialized or will initialize your simulation with) e.g. (3pp0, A): '
        )).replace(' ', '').split(',')

    # make sure the input format is expect
    if len(input_info) != 2:
        raise ValueError(
            "The input must only be PDB_code and chain_id, separated by a comma."
        )
    elif type(input_info[0]) == str and type(input_info[1]) == str:
        pdb_chainid = tuple(input_info)
    else:
        raise ValueError("The input must be a string (PDB_code,chain_id).")

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

    # define indices of the residues relevant to a list of 12 collective variables relevant to kinase conformational changes. These variables include: angle between aC and aE helices, the key K-E salt bridge, DFG-Phe conformation (two distances), X-DFG-Phi, X-DFG-Psi, DFG-Asp-Phi, DFG-Asp-Psi, DFG-Phe-Phi, DFG-Phe-Psi, DFG-Phe-Chi1, and the FRET L-S distance. All features are under the current numbering of the structure provided.
    key_res = []
    # angle between aC and aE helices
    key_res.append(numbering[20])  # residue 21 (res1 in aC)
    key_res.append(numbering[28])  # res29 (res2 in aC)
    key_res.append(numbering[60])  # res61 (res1 in aE)
    key_res.append(numbering[62])  # res63 (res2 in aE)

    # key salt bridge
    key_res.append(numbering[16])  # res17 (K in beta3)
    key_res.append(numbering[23])  # res24 (E in aC)

    # DFG conformation and Phe conformation
    key_res.append(numbering[27])  # res28 (ExxxX)
    key_res.append(numbering[81])  # res82 (DFG-Phe)

    # X-DFG Phi/Psi
    key_res.append(numbering[79])  # res80 (X-DFG)

    # DFG-Asp Phi/Psi
    key_res.append(numbering[80])  # res81 (DFG-Asp)

    # FRET distance
    key_res.append(
        numbering[16] + 120
    )  # not in the list of 85 (equivalent to Aura"S284"), use the 100% conserved beta III K as a reference

    key_res.append(
        numbering[16] + 61
    )  # not in the list of 85 (equivalent to Aura"L225"), use the 100% conserved beta III K as a reference

    # print out kinase information
    print("---------------------Results----------------------")
    print("Kinase ID: " + str(kinase_id))
    print("Kinase name: " + str(name))
    print("Pocket residues: " + str(pocket_seq))
    print("Structure ID: " + str(struct_id))
    print("Numbering of the 85 pocket residues: " + str(numbering))
    print("Residues involved in collective variables: " + str(key_res))

    return pdb_chainid, kinase_id, name, struct_id, pocket_seq, numbering, key_res


def features(pdb_chainid, numbering):
    '''    
    Download the pdb structure corresponding to the given PDB code and chain index, where the atom indices of the relevant atoms will be inferred and used to calculate dihedrals and distances as collective variables. 
    '''
    # download the pdb structure
    cmd = 'wget -q http://www.pdb.org/pdb/files/' + str(
        pdb_chainid[0]) + '.pdb'
    subprocess.call(cmd, shell=True)

    # get topology info from the structure
    topology = md.load(str(pdb_chainid[0]) + '.pdb').topology
    table, bonds = topology.to_dataframe()
    atoms = table.values
    chain_index = ord(str(pdb_chainid[1]).lower(
    )) - 97  # translate a letter chain id into a number index (A->0, B->1 etc)

    # get the array of atom indices for the calculation of:
    #       * eight dihedrals (a 8*4 array where each row contains indices of the four atoms for each dihedral)
    #       * five ditances (a 5*2 array where each row contains indices of the two atoms for each dihedral)
    dih = np.zeros(shape=(8, 4), dtype=int, order='C')
    dis = np.zeros(shape=(5, 2), dtype=int, order='C')

    # name list of the dihedrals and distances
    dih_names = [
        'aC_rot', 'xDFG_phi', 'xDFG_psi', 'dFG_phi', 'dFG_psi', 'DfG_phi',
        'DfG_psi', 'DfG_chi'
    ]
    dis_names = ['K_E1', 'K_E2', 'DFG_conf1', 'DFG_conf2', 'fret']

    # parse the topology info
    for line in atoms:
        # for the specified chain
        if line[5] == chain_index:
            # dihedral 1: between aC and aE helices
            dih[0][0] = line[0] if line[3] == numbering[20] and line[
                1] == 'CA' else dih[0][0]
            dih[0][1] = line[0] if line[3] == numbering[28] and line[
                1] == 'CA' else dih[0][1]
            dih[0][2] = line[0] if line[3] == numbering[60] and line[
                1] == 'CA' else dih[0][2]
            dih[0][3] = line[0] if line[3] == numbering[62] and line[
                1] == 'CA' else dih[0][3]

            # dihedral 2 & 3: X-DFG Phi & Psi
            dih[1][0] = line[0] if line[3] == numbering[78] and line[
                1] == 'C' else dih[1][0]
            dih[1][1] = line[0] if line[3] == numbering[79] and line[
                1] == 'N' else dih[1][1]
            dih[1][2] = line[0] if line[3] == numbering[79] and line[
                1] == 'CA' else dih[1][2]
            dih[1][3] = line[0] if line[3] == numbering[79] and line[
                1] == 'C' else dih[1][3]
            dih[2][0] = dih[1][1]
            dih[2][1] = dih[1][2]
            dih[2][2] = dih[1][3]
            dih[2][3] = line[0] if line[3] == numbering[80] and line[
                1] == 'N' else dih[2][3]

            # dihedral 4 & 5: DFG-Asp Phi & Psi
            dih[3][0] = dih[1][3]
            dih[3][1] = dih[2][3]
            dih[3][2] = line[0] if line[3] == numbering[80] and line[
                1] == 'CA' else dih[3][2]
            dih[3][3] = line[0] if line[3] == numbering[80] and line[
                1] == 'C' else dih[3][3]
            dih[4][0] = dih[3][1]
            dih[4][1] = dih[3][2]
            dih[4][2] = dih[3][3]
            dih[4][3] = line[0] if line[3] == numbering[81] and line[
                1] == 'N' else dih[4][3]

            # dihedral 6 & 7: DFG-Phe Phi & Psi
            dih[5][0] = dih[3][3]
            dih[5][1] = dih[4][3]
            dih[5][2] = line[0] if line[3] == numbering[81] and line[
                1] == 'CA' else dih[5][2]
            dih[5][3] = line[0] if line[3] == numbering[81] and line[
                1] == 'C' else dih[5][3]
            dih[6][0] = dih[5][1]
            dih[6][1] = dih[5][2]
            dih[6][2] = dih[5][3]
            dih[6][3] = line[0] if line[3] == numbering[82] and line[
                1] == 'N' else dih[6][3]

            # dihedral 8: DFG-Phe Chi
            dih[7][0] = dih[5][1]
            dih[7][1] = dih[5][2]
            dih[7][2] = line[0] if line[3] == numbering[81] and line[
                1] == 'CB' else dih[7][2]
            dih[7][3] = line[0] if line[3] == numbering[81] and line[
                1] == 'CG' else dih[7][3]

            # distance 1 & 2: K-E salt bridge
            dis[0][0] = line[0] if line[3] == numbering[16] and line[
                1] == 'NZ' else dis[0][0]
            dis[0][1] = line[0] if line[3] == numbering[23] and line[
                1] == 'OE1' else dis[0][1]
            dis[1][0] = line[0] if line[3] == numbering[16] and line[
                1] == 'NZ' else dis[1][0]
            dis[1][1] = line[0] if line[3] == numbering[23] and line[
                1] == 'OE2' else dis[1][1]

            # distance 3 & 4: DFG conformation-related distances
            dis[2][0] = line[0] if line[3] == numbering[27] and line[
                1] == 'CA' else dis[2][0]
            dis[2][1] = line[0] if line[3] == numbering[81] and line[
                1] == 'CZ' else dis[2][1]
            dis[3][0] = line[0] if line[3] == numbering[16] and line[
                1] == 'CA' else dis[3][0]
            dis[3][1] = dis[2][1]

            # distance 5: FRET distance
            dis[4][0] = line[0] if line[3] == int(
                numbering[80] + 10) and line[1] == 'CA' else dis[4][0]
            dis[4][1] = line[0] if line[3] == int(
                numbering[80] - 20) and line[1] == 'CA' else dis[4][1]

        if line[5] > chain_index:
            break

    # check if there is any missing coordinates; if so, skip dihedral/distance calculation for those residues
    check_flag = 1
    for i in range(len(dih)):
        if 0 in dih[i]:
            dih = np.delete(dih, (i), axis=0)
            print(
                'The "' + str(dih_name[i]) +
                '" dihedral will not be computed due to missing coordinates.')
            dih_names.remove(dih_names[i])
            check_flag = 0
        else: # the atom indices given to mdtraj should be 0-based
            dih[i] = dih[i]-1
    for i in range(len(dis)):
        if 0 in dis[i]:
            dis = np.delete(dis, (i), axis=0)
            print(
                'The "' + str(dis_names[i]) +
                '" distance will not be calculated due to missing coordinates.'
            )
            dis_names.remove(dis_names[i])
            check_flag = 0
        else: # the atom indices given to mdtraj should be 0-based
            dis[i] = dis[i]-1
    if check_flag:
        print(
            "There is no missing coordinates.  All dihedrals and distances will be computed."
        )

    # calculate the dihedrals and distances for the user-specifed structure (a static structure or an MD trajectory)
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
    dihedrals = md.compute_dihedrals(traj, dih)
    distances = md.compute_distances(traj, dis)
    print("Key dihedrals relevant to kinase conformation are as follows:")
    print(dih_names)
    #print(dihedrals/np.pi*180) # dihedrals in degrees
    print(dihedrals) # dihedrals in radians
    print("Key distances relevant to kinase conformation are as follows:")
    print(dis_names)
    print(distances)

    # clean up
    rm_file = 'rm ./' + str(pdb_chainid[0]) + '.pdb'
    subprocess.call(rm_file, shell=True)
    del traj, dih, dis

    return dihedrals, distances


if __name__ == "pf":
    (pdb_chainid, kinase_id, name, struct_id, pocket_seq, numbering,
     key_res) = basics()
#if __name__ == "__features__":
    features(pdb_chainid, numbering)

"""
protein_features.py
This is a tool to featurize kinase conformational changes through the entire Kinome.

"""

# Setup general logging (guarantee output/error message in case of interruption)
# TODO: Can we log to the terminal instead?
import logging
logger = logging.getLogger(__name__)
logging.root.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("urllib3").setLevel(logging.WARNING)


def key_klifs_residues(numbering):
    """
    Retrieve a list of PDB residue indices relevant to key kinase conformations mapped via KLIFS.

    Define indices of the residues relevant to a list of 12 collective variables relevant to
    kinase conformational changes. These variables include: angle between aC and aE helices,
    the key K-E salt bridge, DFG-Phe conformation (two distances), X-DFG-Phi, X-DFG-Psi,
    DFG-Asp-Phi, DFG-Asp-Psi, DFG-Phe-Phi, DFG-Phe-Psi, DFG-Phe-Chi1, and the FRET L-S distance.
    All features are under the current numbering of the structure provided.

    Parameters
    ----------
    numbering : list of int
        numbering[klifs_index] is the residue number for the given PDB file corresponding to KLIFS residue index 'klifs_index'

    Returns
    -------
    key_res : list of int
        Key residue indices

    """

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
    # not in the list of 85 (equivalent to Aura"S284"), use the 100% conserved beta III K as a reference
    key_res.append(numbering[16] + 120)

    # not in the list of 85 (equivalent to Aura"L225"), use the 100% conserved beta III K as a reference
    key_res.append(numbering[16] + 61)

    return key_res

def compute_simple_protein_features(pdbid, chainid, coordfile, numbering):
    """
    This function takes the PDB code, chain id and certain coordinates of a kinase from
    a command line and returns its structural features.

    Parameters
    ----------
    pdbid : str
        The PDB code of the inquiry kinase.
    chainid : str
        The chain index of the inquiry kinase.
    coordfile : str
        Specifies the file constaining the kinase coordinates (either a pdb file or a trajectory, i.e. trj, dcd, h5)
    numbering : list of int
        The residue indices of the 85 pocket residues specific to the structure.

    Returns
    -------
    dihedrals: list of floats
        A list (one frame) or lists (multiple frames) of dihedrals relevant to kinase conformation.
    distances: list of floats
        A list (one frame) or lists (multiple frames) of intramolecular distances relevant to kinase conformation.

    .. todo ::

       Instead of featurizing on dihedrals (which are discontinuous), it's often better to use sin() and cos()
       of the dihedrals or some other non-discontinuous representation.

    .. todo :: Instead of a PDB file or a trj/dcd/h5, accept an MDTraj.Trajectory---this will be much more flexible.

    .. todo :: Use kwargs with sensible defaults instead of relying only on positional arguments.


    """
    import tempfile
    import os
    import mdtraj as md
    import numpy as np
    import sys
    np.set_printoptions(threshold=sys.maxsize)

    pdb_file = None

    # A safer way to download files as wget may not exist on systems such MacOS
    # TODO: Since we retrieve the PDB file in multiple pieces of code, let's refactor this into one utility function
    # to avoid code duplication.
    import urllib

    # get toppology info either from fixed pdb or original pdb file (based on input) 
    # if analyzing a trajectory 
    if coordfile == 'dcd':
        traj = md.load(subprocess.check_output('ls *dcd', shell=True).strip(b'\n').decode("utf-8"),top = str(pdbid) + '_minimized.pdb')
        #traj = md.load(str(pdbid) + '.dcd',top = str(pdbid) + '_fixed_solvated.pdb')
        topology = md.load(str(pdbid)+'_fixed.pdb').topology

    chain_lst = [] # get a list of chains to identify the chain of interest
    import string
    with urllib.request.urlopen('http://www.pdb.org/pdb/files/{}.pdb'.format(pdbid)) as response:
        pdb_file = response.read()
        # check whether chains are in reverse order (e.g. B, A)
        check = pdb_file.decode().split()
        for i in range(len(check)):
            if check[i] == 'ATOM':
                if check[i+4] in list(string.ascii_uppercase):
                    if check[i+4] not in chain_lst:
                        chain_lst.append(check[i+4])
            elif check[i] == 'HETATM':
                if check[i+4] in list(string.ascii_uppercase):
                    if 'HET{}'.format(check[i+4]) not in chain_lst:
                        chain_lst.append('HET{}'.format(check[i+4]))
            elif check[i] == 'CONECT':
                break
    with tempfile.TemporaryDirectory() as pdb_directory:
        pdb = os.path.join(pdb_directory,'{}.pdb'.format(pdbid))
        with open(pdb, 'w') as file:
            file.write(pdb_file.decode())
            # load traj before the temp pdb file was removed
            if coordfile == 'pdb':
                print("loading top from pdb")
                traj = md.load(pdb)
                topology = md.load(pdb).topology
    coord = traj.xyz
    table, bonds = topology.to_dataframe()
    atoms = table.values
    #print(atoms)

    # translate a letter chain id into a number index based on chain order
    # identify the index of the chain of interest based on it's position in chain list
    chain_index = chain_lst.index(chainid) 
    
    # get the array of atom indices for the calculation of:
    #       * eight dihedrals (a 8*4 array where each row contains indices of the four atoms for each dihedral)
    #       * five ditances (a 5*2 array where each row contains indices of the two atoms for each dihedral)
    dih = np.zeros(shape=(10, 4), dtype=int, order='C')
    dis = np.zeros(shape=(5, 2), dtype=int, order='C')

    # name list of the dihedrals and distances
    dih_names = [
        'aC_rot', 'xDFG_phi', 'xDFG_psi', 'dFG_phi', 'dFG_psi', 'dFG_chi1', 'dFG_chi2', 'DfG_phi',
        'DfG_psi', 'DfG_chi1'
    ]
    dis_names = ['K_E1', 'K_E2', 'DFG_conf1', 'DFG_conf2', 'fret']

    # parse the topology info
    '''
    The coordinates are located by row number (usually is atom index minus one, which is also why it's zero-based)
    by mdtraj but when the atom indices are not continuous there is a problem so a safer way to locate the coordinates
    is through row number (as a fake atom index) in case the atom indices are not continuous.
    '''

    # dihedral 0: between aC and aE helices
    dih[0][0] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[20]))
    dih[0][1] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[28]))
    dih[0][2] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[60]))
    dih[0][3] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[62]))

    # dihedral 1 & 2: X-DFG Phi & Psi
    dih[1][0] = topology.select("chainid {} and residue {} and name C".format(chain_index,numbering[78]))
    dih[1][1] = topology.select("chainid {} and residue {} and name N".format(chain_index,numbering[79]))
    dih[1][2] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[79]))
    dih[1][3] = topology.select("chainid {} and residue {} and name C".format(chain_index,numbering[79]))

    dih[2][0] = dih[1][1]
    dih[2][1] = dih[1][2]
    dih[2][2] = dih[1][3]
    dih[2][3] = topology.select("chainid {} and residue {} and name N".format(chain_index,numbering[80]))

    # dihedral 3 & 4: DFG-Asp Phi & Psi
    dih[3][0] = dih[1][3]
    dih[3][1] = dih[2][3]
    dih[3][2] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[80]))
    dih[3][3] = topology.select("chainid {} and residue {} and name C".format(chain_index,numbering[80])) 

    dih[4][0] = dih[3][1]
    dih[4][1] = dih[3][2]
    dih[4][2] = dih[3][3]
    dih[4][3] = topology.select("chainid {} and residue {} and name N".format(chain_index,numbering[81]))

    # dihedral 5 & 6: DFG-Asp Chi1 & Chi2
    dih[5][0] = dih[2][3] # DFG-Asp N
    dih[5][1] = dih[3][2] # DFG-Asp CA
    dih[5][2] = topology.select("chainid {} and residue {} and name CB".format(chain_index,numbering[80])) # DFG-Asp CB
    dih[5][3] = topology.select("chainid {} and residue {} and name CG".format(chain_index,numbering[80])) # DFG-Asp CG

    dih[6][0] = dih[5][1]
    dih[6][1] = dih[5][2]
    dih[6][2] = dih[5][3]
    dih[6][3] = topology.select("chainid {} and residue {} and name OD1".format(chain_index,numbering[80])) #DFG-Asp OD1

    # dihedral 7 & 8: DFG-Phe Phi & Psi
    dih[7][0] = dih[3][3]
    dih[7][1] = dih[4][3]
    dih[7][2] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[81]))              
    dih[7][3] = topology.select("chainid {} and residue {} and name C".format(chain_index,numbering[81]))
    
    dih[8][0] = dih[7][1]
    dih[8][1] = dih[7][2]
    dih[8][2] = dih[7][3]
    dih[8][3] = topology.select("chainid {} and residue {} and name N".format(chain_index,numbering[82]))

    # dihedral 9: DFG-Phe Chi1
    dih[9][0] = dih[4][3] #DFG-Phe N
    dih[9][1] = dih[7][2] #DFG-Phe CA
    dih[9][2] = topology.select("chainid {} and residue {} and name CB".format(chain_index,numbering[81])) # DFG-Phe CB
    dih[9][3] = topology.select("chainid {} and residue {} and name CG".format(chain_index,numbering[81])) # DFG-Phe CG

    # distance 0 & 1: K-E salt bridge
    dis[0][0] = topology.select("chainid {} and residue {} and name NZ".format(chain_index,numbering[16]))
    dis[0][1] = topology.select("chainid {} and residue {} and name OE1".format(chain_index,numbering[23]))
    dis[1][0] = dis[0][0]
    dis[1][1] = topology.select("chainid {} and residue {} and name OE2".format(chain_index,numbering[23]))

    # distance 2 & 3: DFG conformation-related distances
    dis[2][0] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[27]))
    dis[2][1] = topology.select("chainid {} and residue {} and name CZ".format(chain_index,numbering[81]))
    dis[3][0] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[16]))
    dis[3][1] = dis[2][1]

    # distance 4: FRET distance
    dis[4][0] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[16]+120))
    dis[4][1] = topology.select("chainid {} and residue {} and name CA".format(chain_index,numbering[16]+61))

    # check if there is any missing coordinates; if so, skip dihedral/distance calculation for those residues
    check_flag = 1
    for i in range(len(dih)):
        if 0 in dih[i]:
            dih[i] = [0,0,0,0]
            #logging.info(
            #    'The "' + str(dih_names[i]) +
            #    '" dihedral will not be computed due to missing coordinates.')
            check_flag = 0
    for i in range(len(dis)):
        if 0 in dis[i]:
            dis[i] = [0,0]
            #logging.info(
            #    'The "' + str(dis_names[i]) +
            #    '" distance will not be calculated due to missing coordinates.'
            #)
            check_flag = 0
    #if check_flag:
        #logging.info(
        #    "There is no missing coordinates.  All dihedrals and distances will be computed."
        #)
    # calculate the dihedrals and distances for the user-specifed structure (a static structure or an MD trajectory)
    dihedrals = md.compute_dihedrals(traj, dih)/np.pi*180
    distances = md.compute_distances(traj, dis)

    # option to log the results
    '''
    logging.info("Key dihedrals relevant to kinase conformation are as follows:")
    logging.info(dih_names)
    #logging.info(dihedrals/np.pi*180) # dihedrals in degrees
    logging.info(dihedrals)  # dihedrals in radians
    logging.info("Key distances relevant to kinase conformation are as follows:")
    logging.info(dis_names)
    logging.info(distances)
    '''

    # clean up
    del traj, dih, dis

    return dihedrals, distances

def hybrid_docking(receptor_path, molecules_path, docked_molecules_path, n_poses=10):
    """Automated hybrid docking of small molecules to a receptor.

    Parameters
    ----------
    receptor_path : str
        Path to PDB file containing receptor and reference ligand, or pre-prepared receptor for docking
    molecules_path : str
        Path to file containing one or more molecules (in OpenEye readable format) to be docked.
        (For example, list of SMILES)
    docked_molecules_path : str
        Path to output file to be created to contain docked molecules
        Uses OpenEye recognized file extension, such as .mol2 or .sdf
    n_poses : int, optional, default=1
        Number of docked poses to generate
    receptor_filename : str, optional, default=None
        If not None, the pre-prepared receptor is loaded

    TODO: How can this API be improved?

    """
    from .docking import create_receptor, load_receptor, pose_molecule
    from openeye import oedocking, oechem
    #import openmoltools as moltools # TODO: Bring these methods into this module

    # Try to load pre-prepared receptor from specified file
    receptor = oechem.OEGraphMol()
    print('Attempting to load receptor from {}...'.format(receptor_path))
    if not oedocking.OEReadReceptorFile(receptor, receptor_path):
        # Load complex of receptor and reference ligand
        complex_istream = oechem.oemolistream(receptor_path)
        complex = oechem.OEGraphMol()
        oechem.OEReadMolecule(complex_istream, complex)

        # Attempt to split into components and build receptor based on reference ligand
        print('Attempting to split complex into components...')
        ligand = oechem.OEGraphMol()
        protein = oechem.OEGraphMol()
        water = oechem.OEGraphMol()
        other = oechem.OEGraphMol()
        if oechem.OESplitMolComplex(ligand, protein, water, other, complex):
            # Create receptor using bound ligand reference
            print('Creating receptor using reference ligand...')
            oedocking.OEMakeReceptor(receptor, protein, ligand)
            # TODO: We can store prepared receptor file if desired
            oedocking.OEWriteReceptorFile(receptor, '/home/guoj1/projects/INSPIRE/kinomodel/kinomodel/data/docking/prepared_receptor.oeb')

        else:
            raise Exception('Could not split specified PDB file {} into receptor and reference ligand'.format(receptor_path))

    # Open file for writing docked molecules
    docked_molecules_ostream = oechem.oemolostream(docked_molecules_path)

    # Configure omega
    # From canonical recipe: https://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
    from openeye import oeomega
    omega = oeomega.OEOmega()
    omega.SetIncludeInput(False)
    omega.SetCanonOrder(False)
    omega.SetSampleHydrogens(True)
    eWindow = 15.0
    omega.SetEnergyWindow(eWindow)
    omega.SetMaxConfs(800)
    omega.SetRMSThreshold(1.0)

    # Dock all molecules requested
    dock_method = oedocking.OEDockMethod_Hybrid2
    dock_resolution = oedocking.OESearchResolution_Standard
    dock = oedocking.OEDock(dock_method, dock_resolution)
    dock.Initialize(receptor)
    molecules_istream = oechem.oemolistream(molecules_path)
    molecule = oechem.OEGraphMol()
    for molecule in molecules_istream.GetOEMols():
        print("docking", molecule.GetTitle())
        #docked_molecules = pose_molecule(receptor, molecule, n_poses=n_poses)
        #molecule = moltools.openeye.get_charges(molecule, keep_confs=10)

        # Generate conformers
        if not omega(molecule):
            continue

        # Apply charges
        from openeye import oequacpac
        oequacpac.OEAssignCharges(molecule, oequacpac.OEAM1BCCELF10Charges())

        # Dock
        docked_molecule = oechem.OEGraphMol()
        dock.DockMultiConformerMolecule(docked_molecule, molecule)
        sdtag = oedocking.OEDockMethodGetName(dock_method)
        oedocking.OESetSDScore(docked_molecule, dock, sdtag)
        dock.AnnotatePose(docked_molecule)
        oechem.OEWriteMolecule(docked_molecules_ostream, docked_molecule)

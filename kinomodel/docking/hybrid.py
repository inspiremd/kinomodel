def hybrid_docking(receptor_path, molecules_path, docked_molecules_path, n_poses=1):
    """Automated hybrid docking of small molecules to a receptor.

    Parameters
    ----------
    receptor_path : str
        Path to file containing receptor and reference ligand.
    molecules_path : str
        Path to file containing one or more molecules (in OpenEye readable format) to be docked.
        (For example, list of SMILES)
    docked_molecules_path : str
        Path to output file to be created to contain docked molecules
        Uses OpenEye recognized file extension, such as .mol2 or .sdf
    n_poses : int, optional, default=1
        Number of docked poses to generate

    TODO: How can this API be improved?

    """
    from .docking import create_receptor, load_receptor, pose_molecule
    from openeye import oedocking, oechem
    import openmoltools as moltools

    # Load receptor
    complex_istream = oechem.oemolistream(receptor_path)
    complex = oechem.OEGraphMol()
    oechem.OEReadMolecule(complex_istream, complex)

    # for input molecule inmol ...
    ligand = oechem.OEGraphMol()
    protein = oechem.OEGraphMol()
    water = oechem.OEGraphMol()
    other = oechem.OEGraphMol()

    if oechem.OESplitMolComplex(ligand, protein, water, other, complex):
        # Create receptor using bound ligand reference
        receptor = oechem.OEGraphMol()
        oedocking.OEMakeReceptor(receptor, protein, ligand)

        # Open file for writing docked molecules
        docked_molecules_ostream = oechem.oemolostream(docked_molecules_path)

        # Dock all molecules requestd
        molecules_istream = oechem.oemolistream(molecules_path)
        molecule = oechem.OEGraphMol()
        while oechem.OEReadMolecule(molecules_istream, molecule):
            docked_molecules = pose_molecule(receptor, molecule, n_poses=n_poses)
            oechem.OEWriteMolecule(docked_molecules_ostream, docked_molecules)

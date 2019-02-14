def mmpbsa(inpcrd_filenames, prmtop_filenames, niterations=1000):
    """
    Run simple MM-PBSA like workflow to estimate free energy of binding from enthalpies of interaction.

    This method can accept either explicit or implicit solvent.

    .. todo :: Allow user to leave out 'receptor' if desired.

    Parameters
    ----------
    inpcrd_filename : dict of str
        inpcrd_filename[phase] is the Amber inpcrd filename for phase in 'complex', 'ligand', and (optionally) 'receptor'
    prmtop_filenames : dict of str
        prmtop_filename[phase] is the Amber prmtop filename for phase in 'complex', 'ligand', and (optionally) 'receptor'
    niterations : int, optional, default=100
        Maximum number of iterations to run.

    Returns
    -------
    DeltaG, dDeltaG : float
        Estimated free energy of binding (in kT)

    """
    from simtk import unit, openmm
    from simtk.openmm import app

    # Thermodynamic and simulation parameters
    # TODO: Make these options
    temperature = 300 * unit.kelvin
    collision_rate = 1.0 / unit.picoseconds
    hydrogen_mass = 3.5 * unit.amu
    timestep = 4.0 * unit.femtoseconds
    pressure = 1.0 * unit.atmospheres

    # Store thermal energy
    kT = unit.kB_BOLTZMANN_CONSTANT * temperature

    # TODO: Autodetect nonbondedMethod

    # Termination criteria
    nsteps_per_iteration = 250 # 1 picosecond
    # TODO: Use statistical error threshold to drive this to convergence

    enthalpies = dict()
    for phase in phases:
        print('Simulating phase {}'.format(phase))

        # Load Amber system
        prmtop = app.AmberPrmtopFile(prmtop_filenames[phase])
        inpcrd = app.AmberInpcrdFile(inpcrd_filenames[phase])
        # Create system
        system = prmtop.createSystem(nonbondedMethod=app.PME, constraints=app.Hbonds, hydrogenMass=hydrogen_mass)
        # Create integrator and context
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator)
        context.setPositions(inpcrd.getPositions())
        context.setPeriodicBoxVectors(inpcrd.getPeriodicBoxVectors())
        # Sample and store enthalpies
        enthalpies[phase] = np.zeros([niterations])
        for iteration in range(niterations):
            integrator.step(nsteps)
            state = context.getState(getEnergy=True)
            potential_energy = state.getPotentialEnergy()
            volume = state.getPeriodicBoxVolume()
            enthalpies[phase][iteration] = (potential_energy + pressure*volume) / kT

    # Use automatic equilibration detection and pymbar.timeseries to subsample
    from pymbar import timeseries
    for phase in phases:
        [t0, g, Neff_max] = timeseries.detectEquilibration(enthalpies[phase]) # compute indices of uncorrelated timeseries
        enthalpies[phase] = enthalpies[phase][t0:]
        indices = timeseries.subsampleCorrelatedData(enthalpies[phase], g=g)
        enthalpies[phase] = enthalpies[phase][indices]

    # Compute binding free energy from entalpies estimate
    DeltaH = enthalpies['complex'].mean() - enthalpies['ligand'].mean()
    dDeltaH = np.sqrt(enthalpies['complex'].std()**2/len(enthalpies['complex']) + enthalpies['ligand'].std()**2/len(enthalpies['complex']))
    if 'receptor' in enthalpies:
        DeltaH -= enthalpies['receptor'].mean()
        dDeltaH = np.sqrt(dDeltaH**2 + enthalpies['receptor'].std()/len(enthalpies['receptor']))

    return DeltaH, dDeltaH

def mmpbsa_cli():
    """
    MM-PBSA command-line driver.
    """
    # TODO
    pass

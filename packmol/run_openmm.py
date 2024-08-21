import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *

# Load the merged PDB file
#pdb = PDBFile('merged_structure.pdb')
pdb = PDBFile('1l2y.pdb')

# Load the force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Create a Modeller object
modeller = Modeller(pdb.topology, pdb.positions)


modeller.addMembrane(forcefield, lipidType='POPC', minimumPadding=1*nanometer)

# Solvate the system with water
#modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometers)


# Save the solvated system to a PDB file for verification
with open('solvated_system.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, file=f)

# Create the system
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)

# Define the integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Create the simulation object
simulation = Simulation(modeller.topology, system, integrator)

# Set the initial positions
simulation.context.setPositions(modeller.positions)

# Minimize energy
print('Minimizing...')
simulation.minimizeEnergy()

# Equilibrate the system
simulation.context.setVelocitiesToTemperature(300*kelvin)
print('Running Equilibration...')
simulation.step(10000)

# Production run
print('Running Production...')
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter('data.log', 1000, step=True, potentialEnergy=True, temperature=True))

simulation.step(100000)
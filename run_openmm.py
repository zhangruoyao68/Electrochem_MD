import os
from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np

# Parameters
pdb_file = 'tethered_polymer_system.pdb'
output_pdb_file = 'output.pdb'
output_log_file = 'output.log'
temperature = 300 * kelvin
friction = 1 / picosecond
timestep = 2 * femtoseconds
simulation_steps = 10000  # Number of steps to run the simulation
box_padding = 1.0 * nanometers  # Padding around the system for the water box
ionic_strength = 0.1 * molar  # Ionic strength for adding ions

# Load the cleaned PDB file
pdb = PDBFile(pdb_file)

'''
# Create a force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Create a modeller
modeller = Modeller(pdb.topology, pdb.positions)

# Add solvent (water)
modeller.addSolvent(forcefield, model='tip3p', padding=box_padding, ionicStrength=ionic_strength)

# Create a system
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, constraints=HBonds)

# Define an integrator
integrator = LangevinIntegrator(temperature, friction, timestep)

# Create a simulation
simulation = Simulation(modeller.topology, system, integrator)

# Set the initial positions
simulation.context.setPositions(modeller.positions)

# Minimize energy
print('Minimizing energy...')
simulation.minimizeEnergy()

# Set up reporters to output data
simulation.reporters.append(PDBReporter(output_pdb_file, 1000))  # Save PDB every 1000 steps
simulation.reporters.append(StateDataReporter(output_log_file, 1000, time=True, potentialEnergy=True, temperature=True))

# Run the simulation
print('Running simulation...')
simulation.step(simulation_steps)

print('Simulation complete')'''
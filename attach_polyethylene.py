from ase import Atoms
from ase.build import fcc111
from ase.constraints import FixAtoms
import numpy as np

# Parameters
a = 4.05  # Lattice constant for aluminum
size = 10  # Number of unit cells in both x and y directions for a square base
size_z = 1  # Number of unit cells in the z direction (for the slab)
polymer_length = 20  # Number of monomers in the polymer chain
bond_length = 1.54  # C-C bond length in Angstroms
bond_distance_to_surface = 2.0  # Bond distance between C and Al in Angstroms

# Create a bulk aluminum slab with a square base
aluminum_slab = fcc111('Al', size=(size, size, size_z), a=a, vacuum=10.0, orthogonal=True)

# Identify the bottom layer atoms (assuming z=0 is the bottom)
z_positions = aluminum_slab.positions[:, 2]
bottom_layer_indices = np.where(z_positions == z_positions.min())[0]

# Apply a constraint to fix the bottom layer atoms
constraint = FixAtoms(indices=bottom_layer_indices)
aluminum_slab.set_constraint(constraint)

# Create a polyethylene chain (simplified model)
polymer_positions = np.array([[0, 0, i * bond_length] for i in range(polymer_length)])
polymer = Atoms('C' * polymer_length, positions=polymer_positions)

# Select a surface atom to tether the polymer to
tether_point = aluminum_slab.positions[bottom_layer_indices[0]]

# Position the first carbon atom of the polymer at the bonding distance from the tether point
polymer.positions[0] = tether_point + np.array([0, 0, bond_distance_to_surface])

# Adjust the rest of the polymer chain to account for the new position of the first atom
for i in range(1, polymer_length):
    polymer.positions[i] = polymer.positions[0] + np.array([0, 0, i * bond_length])

# Merge the polymer with the slab
combined_system = aluminum_slab + polymer

# Assign unique residue numbers to each monomer in the polymer chain
# Using ASE's info dictionary to store residue information
residue_numbers = list(range(1, polymer_length + 1))
for i in range(polymer_length):
    combined_system.info[f'residue_{i + 1}'] = residue_numbers[i]

# Print the combined system for verification
print("Combined system with covalently bonded polymer and unique residue numbers:")
print(combined_system)

# Visualize the system (optional, requires an appropriate viewer like ASE's default viewer)
# view(combined_system)

# Save the structure
# combined_system.write('tethered_polymer_system.xyz') # To an XYZ file
combined_system.write('tethered_polymer_system.pdb') # To a PDB file
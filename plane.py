from ase import Atoms
from ase.build import bulk, fcc111
from ase.constraints import FixAtoms
import numpy as np

# Parameters
a = 4.05  # Lattice constant for aluminum
size = 10  # Number of unit cells in both x and y directions for a square base
size_z = 1  # Number of unit cells in the z direction (for the slab)

# Create a bulk aluminum slab with a square base
aluminum_slab = fcc111('Al', size=(size, size, size_z), a=a, vacuum=10.0, orthogonal=True)

# Identify the bottom layer atoms (assuming z=0 is the bottom)
z_positions = aluminum_slab.positions[:, 2]
bottom_layer_indices = np.where(z_positions == z_positions.min())[0]

# Apply a constraint to fix the bottom layer atoms
constraint = FixAtoms(indices=bottom_layer_indices)
aluminum_slab.set_constraint(constraint)

# Print the atoms and constraints for verification
print("Aluminum slab with square base and inert bottom surface:")
print(aluminum_slab)
print("Fixed atom indices:", bottom_layer_indices)

# Save the structure to a file (e.g., XYZ format)
aluminum_slab.write('aluminum_slab_square_base.xyz')
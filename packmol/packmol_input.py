import numpy as np
import os
from ase import Atoms
from ase.build import fcc111
from ase.io import write
import requests

def download_pdb(pdb_id, output_file):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(output_file, 'w') as file:
            file.write(response.text)
        print(f"PDB file {pdb_id} downloaded successfully.")
    else:
        print(f"Failed to download PDB file {pdb_id}.")

# Example: Download Bovine Pancreatic Trypsin Inhibitor (BPTI)
download_pdb("1l2y", "1l2y.pdb")

# Parameters
a = 4.05  # Lattice constant for aluminum
size = 10  # Number of unit cells in both x and y directions for a square base
size_z = 1  # Number of unit cells in the z direction (for the slab)
vacuum = 10.0  # Vacuum layer thickness
bond_distance_to_surface = 2.0 # Distance from the surface to the first atom of the protein

# Create a bulk aluminum slab with a square base
aluminum_slab = fcc111('Al', size=(size, size, size_z), a=a, vacuum=vacuum, orthogonal=True)
# Save the slab structure to a PDB file
write('aluminum_slab.pdb', aluminum_slab)

# Calculate the dimensions of the aluminum slab
slab_x_min = np.min(aluminum_slab.positions[:, 0])
slab_x_max = np.max(aluminum_slab.positions[:, 0])
slab_y_min = np.min(aluminum_slab.positions[:, 1])
slab_y_max = np.max(aluminum_slab.positions[:, 1])
slab_z_max = np.max(aluminum_slab.positions[:, 2])

# Create Packmol input script
packmol_input = f"""
tolerance 2.0

filetype pdb

output merged_structure.pdb

structure aluminum_slab.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure

structure 1l2y.pdb
  number 1
  inside box {slab_x_min} {slab_y_min} {slab_z_max + bond_distance_to_surface} {slab_x_max} {slab_y_max} {slab_z_max + bond_distance_to_surface + 30.0}
end structure
"""

# Write the Packmol input script to a file
with open('packmol_input.inp', 'w') as f:
    f.write(packmol_input)

# Run Packmol
os.system('packmol < packmol_input.inp')
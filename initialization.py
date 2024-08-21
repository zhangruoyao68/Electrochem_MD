import numpy as np
from openmm.app import PDBFile, Modeller, ForceField, Element
from openmm import unit
from pdbfixer import PDBFixer
from openmm.vec3 import Vec3

# Define the peptide sequence
sequence = "SKGPG"

# Create a simple PDB file with alpha carbons
def create_initial_pdb(sequence, filename):
    with open(filename, 'w') as f:
        for i, residue in enumerate(sequence):
            x, y, z = i * 3.8, 0, 0  # Simple linear structure
            f.write(f"ATOM  {i+1:5d}  CA  {residue} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n")
        f.write("END\n")

initial_pdb = "initial_ca.pdb"
create_initial_pdb(sequence, initial_pdb)

# Use PDBFixer to add missing atoms
fixer = PDBFixer(filename=initial_pdb)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)

# Get the positions and create a new PDB file
positions = fixer.positions
topology = fixer.topology

# Translate the peptide so that the C-terminal glycine is at z=0
min_z = min(pos[2] for pos in positions)
translation = unit.Quantity(np.array([0, 0, -min_z]), unit.nanometers)
new_positions = [pos + translation for pos in positions]

# Create a surface below the peptide
surface_z = -0.3 * unit.nanometers  # 3 Angstroms below the peptide
surface_size = 5 * unit.nanometers
surface_spacing = 0.1 * unit.nanometers
x_range = np.arange(-surface_size/2, surface_size/2, surface_spacing)
y_range = np.arange(-surface_size/2, surface_size/2, surface_spacing)
surface_positions = []

# Add surface atoms to the topology and positions
modeller = Modeller(topology, new_positions)
chain = modeller.topology.addChain()
for x in x_range:
    for y in y_range:
        residue = modeller.topology.addResidue('SRF', chain)
        atom = modeller.topology.addAtom('X', Element.getBySymbol('C'), residue)
        surface_positions.append(Vec3(x, y, surface_z))

modeller.add(modeller.topology, surface_positions)

# Write the new structure to a PDB file
with open('initial_peptide_on_surface_all_atom.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

print("All-atom initial configuration saved as 'initial_peptide_on_surface_all_atom.pdb'")
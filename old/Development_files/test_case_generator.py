from pymatgen.core import Structure, Molecule, Lattice
from pymatgen.io.cif import CifWriter

# Define test case 1: Water dimer with hydrogen bonding
water_dimer = Molecule(
    species=["O", "H", "H", "O", "H", "H"],
    coords=[
        [0.0, 0.0, 0.0],      # Oxygen of the first water
        [0.9572, 0.0, 0.0],   # Hydrogen 1 of the first water
        [-0.239988, 0.927297, 0.0],  # Hydrogen 2 of the first water
        [2.8, 0.0, 0.0],      # Oxygen of the second water (within H-bonding distance)
        [3.7572, 0.0, 0.0],   # Hydrogen 1 of the second water
        [2.560012, 0.927297, 0.0],  # Hydrogen 2 of the second water
    ],
)

# Define test case 2: Two non-interacting molecules (water and methane far apart)
non_interacting = Molecule(
    species=["O", "H", "H", "C", "H", "H", "H", "H"],
    coords=[
        [0.0, 0.0, 0.0],      # Oxygen of the water molecule
        [0.9572, 0.0, 0.0],   # Hydrogen 1 of water
        [-0.239988, 0.927297, 0.0],  # Hydrogen 2 of water
        [10.0, 0.0, 0.0],     # Carbon of methane (far away from water)
        [10.6, 0.6, 0.0],     # Hydrogen 1 of methane
        [10.6, -0.6, 0.0],    # Hydrogen 2 of methane
        [9.4, 0.6, 0.0],      # Hydrogen 3 of methane
        [9.4, -0.6, 0.0],     # Hydrogen 4 of methane
    ],
)

# Define test case 3: Hydrogen bonding between a carboxyl group of an amino acid and water
amino_acid_water = Molecule(
    species=["C", "C", "O", "O", "H", "N", "H", "H", "H", "O", "H", "H"],
    coords=[
        [0.0, 0.0, 0.0],       # Alpha carbon
        [1.5, 0.0, 0.0],       # Carbon in carboxyl group
        [2.3, 0.8, 0.0],       # Oxygen in carboxyl group (single bond)
        [1.5, -1.2, 0.0],      # Oxygen in carboxyl group (double bond)
        [0.0, 0.9, 0.9],       # Hydrogen on alpha carbon
        [-1.1, 0.0, 0.0],      # Nitrogen in amine group
        [-1.1, 0.9, 0.9],      # Hydrogen on nitrogen
        [-1.1, -0.9, 0.9],     # Another hydrogen on nitrogen
        [0.0, 0.0, 1.2],       # Hydrogen on alpha carbon
        [3.2, 0.8, 0.0],       # Oxygen from water (within hydrogen bonding distance)
        [3.9, 0.8, 0.0],       # Hydrogen 1 of water
        [3.2, 1.6, 0.0],       # Hydrogen 2 of water
    ],
)

# Combine all test cases into a structure
lattice = Lattice.cubic(30.0)  # A large cubic lattice to house all molecules
species = []
coords = []

# Add species and coordinates from all test cases
for molecule in [water_dimer, non_interacting, amino_acid_water]:
    species.extend([site.specie for site in molecule.sites])
    coords.extend([site.coords for site in molecule.sites])

# Create the structure
structure = Structure(lattice, species, coords)

# Write to CIF format
cif_writer = CifWriter(structure)
cif_writer.write_file("./hbond_test_cases_fixed.cif")

print("CIF file saved as 'hbond_test_cases_fixed.cif'")

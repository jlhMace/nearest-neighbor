#### Reads trajectory file and performs nearest neighbor algorithm (test version)
#### Input: .traj file with one Atoms object of n atoms,
####    index of particular atom in Atoms object,
####    neighbor radius in Angstroms
#### Output: Atoms object of nearby atoms,
####    number of neighbors (printed to terminal)

from ase.io.trajectory import Trajectory
from ase.io import read
from ase import Atoms, Atom


# inputs
atom_index = 1 # chosen specific atom in Atoms object
radius = 2 # in Angstroms
# read .traj file
traj = Trajectory('100000atoms-100A.traj')
atoms = traj[-1]

# parameters
box_size = atoms.get_cell_lengths_and_angles()

print(box_size)
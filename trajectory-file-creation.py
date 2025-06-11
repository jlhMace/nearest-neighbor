#### Creates a trajectory file containing 100000 Atom objects
#### Input: none
#### Output: .traj file at ~/Documents/code-projects/pyamff/{file-name}.traj
####    containing 100,000 Atom objects at randomized 3D coordinates
#### in an orthorhombic unit cell of side 100 Angstroms

from ase.io.trajectory import Trajectory
from ase.io import write
from ase import Atoms, Atom
from random import uniform

# set parameters
cutoff = 10
element = 'H'
n = 7           # bins per dimension
a = n*cutoff    # dimensions of unit cell
b = n*cutoff
c = n*cutoff
num = pow(n,3)*50    # number of total atoms
print(num,pow(n,3),num/pow(n,3))

# creates list of Atoms objects
atom_list = Atoms(cell=[a,b,c])
for e in range(num):
    x = uniform(0,a)
    y = uniform(0,b)
    z = uniform(0,c)
    atom_list.append(Atom(element,position=(x,y,z)))

# Configures traj file
traj = Trajectory('trajprop-10Acutoff-70A.traj', 'w', atoms=atom_list)
traj.write(atom_list)
traj.close()

'''
# Test print out atom_list
traj = Trajectory('10000atoms-100A.traj')
for atoms in traj:
    for atom in atoms:
        print(atom.position)
traj.close()
'''

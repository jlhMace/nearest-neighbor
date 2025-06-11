from ase.io.trajectory import Trajectory
from ase.io import read
from ase import Atoms, Atom
from math import fmod
from ase.geometry import distance,get_distances
from timer import Timer


def neighbor_list(atoms,atom_index,cutoff):
    '''Reads trajectory file and performs nearest neighbor algorithm
        Input: Atoms object of n atoms,
            index of particular atom in Atoms object,
            list of 
            neighbor cutoff in Angstroms
        Output: Atoms object of nearby atoms,
            number of neighbors (printed to terminal)'''

    #print(len(atoms),len(bin_list))
    cutoff=cutoff/2

    ### Create Nearest Neighbor list ###
    neighbor_atoms = Atoms()
    for index in range(len(atoms)):
        if atoms.get_distances(index,atom_index)<=cutoff:
            neighbor_atoms.append(atoms[index])

    #print(len(neighbor_atoms))
    return neighbor_atoms


def main():
    ### Timing classes for testing ###
    time_total = Timer("total")
    time_nn = Timer("nn")
    print("Timing start...")
    time_total.start()

    traj_file = 'trajprop-10Acutoff-70A.traj'
    #traj_file = '10000atoms-100A.traj'
    traj = Trajectory(traj_file)  # read traj file
    atoms = traj[-1]
    cutoff = 10 # Angstroms

    print('Atoms:',len(atoms))

    print("Nearest neighbor search:")
    time_nn.start()
    #for index in range(len(atoms)):
    for index in range(len(atoms)):
        neighbor_list(atoms,index,10)  # 10 Angstom cutoff radius
    time_nn.stop()

    print("Total:")
    time_total.stop()


if __name__=='__main__':
    main()
from ase.io.trajectory import Trajectory
from ase.io import read
from ase import Atoms, Atom
from math import fmod
from ase.geometry import distance,get_distances
from timer import Timer
import numpy as np

         

def bin_sort(atoms,cutoff):
    '''Sorts atoms into bins based on cutoff
        Input: Atoms object of n atoms,
            neighbor cutoff in Angstroms
        Output: h*k*l nested lists of atomic indexes in each orthagonal bin'''
    
    ### Initialization ###
    cell_size = atoms.cell.cellpar()  # unit cell size
    
    ### Create Bins ###
    bin_num = [int(cell_size[0]/(cutoff)),int(cell_size[1]/(cutoff)),int(cell_size[2]/(cutoff))]  # number of bins in each dimension
    bin_dim = (cell_size[0]/bin_num[0],cell_size[1]/bin_num[1],cell_size[2]/bin_num[2])  # size of each bin
    #print('Bin num =', bin_num[0]*bin_num[1]*bin_num[2], 'Bin dimension =', bin_dim[0],bin_dim[1],bin_dim[2])

    # make matrix of bins
    bins = []
    atom_list = []
    for i in range(bin_num[0]):
        bins.append([]) # Append an empty sublist inside the list
        for j in range(bin_num[1]):
            bins[i].append([]) # Second dimension
            for k in range(bin_num[2]):
                bins[i][j].append([]) # Third dimension

    # Sort atom indexes into bins 
    for i in range(len(atoms)):
        [x,y,z] = np.floor_divide(atoms[i].position,bin_dim).astype(int)
        bins[x][y][z].append(i)
        atom_list.append((x,y,z))

    return bins, atom_list, bin_num


def bin_cull(index,atom_list,bins,bin_num):
    '''Return list of atom indexes of surrounding grid of index bin
        Input: 3-D nested list of bins with cooresponding atomic indexes,
            index of particular atom in Atoms object,
            list of bin coordinates by atom index
        Output: list of atomic indexes in neighboring bins'''
    
    [x,y,z] = atom_list[index]
    [a,b,c] = bin_num
    [x0,x1]=[(x-1)%(a), (x+1)%(a)]
    [y0,y1]=[(y-1)%(b), (y+1)%(b)]
    [z0,z1]=[(z-1)%(c), (z+1)%(c)]
    #print([x,y,z], [x0,x1],[y0,y1],[z0,z1])
    bin_list = []
    for i in [x0,x,x1]:
        for j in [y0,y,y1]:
            for k in [z0,z,z1]:
                #print(i,j,k)
                for e in bins[i][j][k]:
                    bin_list.append(e)

    return bin_list


def neighbor_list(atoms,atom_index,bin_list,cutoff):
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
    for index in bin_list:
        if atoms.get_distances(index,atom_index)<=cutoff:
            neighbor_atoms.append(atoms[index])

    #print(len(neighbor_atoms))
    return neighbor_atoms



def main():
    ### Timing classes for testing ###
    time_total = Timer("total")
    time_bin = Timer("bin")
    time_nn = Timer("nn")
    print("Timing start...")
    time_total.start()

    traj_file = 'trajprop-10Acutoff-5.traj'
    #traj_file = '10000atoms-100A.traj'
    traj = Trajectory(traj_file)  # read traj file
    atoms = traj[-1]
    cutoff = 10 # Angstroms

    print("Bin sorting:")
    time_bin.start()
    bins, atoms_list, bin_num = bin_sort(atoms,cutoff)
    time_bin.stop()

    print('Atoms, # of bins, atoms/bin')
    n=len(atoms)
    number=bin_num[0]*bin_num[1]*bin_num[2]
    print(n,number,n/number)

    print("Nearest neighbor search:")
    time_nn.start()
    #for index in range(len(atoms)):
    for index in range(len(atoms)):
        bin_list = bin_cull(index,atoms_list,bins,bin_num)
        neighbor_list(atoms,index,bin_list,10)  # 10 Angstom cutoff radius
    time_nn.stop()

    print("Total:")
    time_total.stop()


if __name__=='__main__':
    main()
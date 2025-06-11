from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.geometry import distance,get_distances
from timer import Timer
import numpy as np

         

def bin_sort(atoms,cutoff):
    '''Sorts atoms into bins based on cutoff
        Input:  atoms:      Atoms object of n atoms,
                cutoff:     neighbor cutoff in Angstroms
        Output: list_index_by_bin:       h*k*l nested list of atomic indexes in each orthagonal bin,
                atom_list:  list of 3D tuple indicating bin location by atomic index,
                bin_num:    3x1 list of number of bins in unit cell'''
    
    ### Initialization ###
    cell_size = atoms.cell.cellpar()  # unit cell size
    
    ### Create Bins ###
    bin_num = [int(cell_size[0]/(cutoff)),int(cell_size[1]/(cutoff)),int(cell_size[2]/(cutoff))]  # number of bins in each dimension
    bin_dim = (cell_size[0]/bin_num[0],cell_size[1]/bin_num[1],cell_size[2]/bin_num[2])  # size of each bin
    #print('Bin num =', bin_num[0]*bin_num[1]*bin_num[2], 'Bin dimension =', bin_dim[0],bin_dim[1],bin_dim[2])

    # make matrix of bins
    list_index_by_bin = []
    atom_list = []
    for i in range(bin_num[0]):
        list_index_by_bin.append([]) # Append an empty sublist inside the list
        for j in range(bin_num[1]):
            list_index_by_bin[i].append([]) # Second dimension
            for k in range(bin_num[2]):
                list_index_by_bin[i][j].append([]) # Third dimension

    # Sort atom indexes into bins 
    for i in range(len(atoms)):
        [x,y,z] = np.floor_divide(atoms[i].position,bin_dim).astype(int)
        list_index_by_bin[x][y][z].append(i)
        atom_list.append((x,y,z))

    return list_index_by_bin, atom_list, bin_num


def bin_cull(index,atom_list,list_index_by_bin,bin_num):
    '''Return list of atom indexes of surrounding grid of index bin
        Input:  index:      3-D nested list of bins with cooresponding atomic indexes,
                atom_list:  index of particular atom in Atoms object,
                list_index_by_bin:       h*k*l nested list of atom indicies by bin coordinate
                bin_num:    3x1 list of number of bins in unit cell
        Output: bin_list:   list of atomic indexes in neighboring bins'''
    
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
                for e in list_index_by_bin[i][j][k]:
                    bin_list.append(e)

    if a or b or c < 3:
        bin_list = np.unique(bin_list)

    return bin_list


def neighbor_list(atoms,atom_index,bin_list,cutoff):
    '''Reads trajectory file and performs nearest neighbor algorithm
        Input:  atoms:          Atoms object of n atoms,
                atom_index:     index of particular atom in Atoms object,
                bin_list:       list of 
                cutoff:         neighbor cutoff in Angstroms
        Output: neighbor_atoms: Atoms object of nearby atoms'''

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


    ### Trajectory file import and setup ###
    traj_file = 'trajprop-10Acutoff-10A.traj'
    traj = Trajectory(traj_file)  # read traj file
    atoms = traj[-1]
    cutoff = 10 # Angstroms


    ### Sorting and neighbor list functions ###
    print("Bin sorting:")
    time_bin.start()
    list_index_by_bin, atoms_list, bin_num = bin_sort(atoms,cutoff)
    print(list_index_by_bin)
    time_bin.stop()

    print('Atoms, # of bins, atoms/bin')
    n=len(atoms)
    number=bin_num[0]*bin_num[1]*bin_num[2]
    print(n,'   ',number,'     ',n/number)

    print("Nearest neighbor search:")
    time_nn.start()
    #for index in range(len(atoms)):
    for index in range(len(atoms)):
        bin_list = bin_cull(index,atoms_list,list_index_by_bin,bin_num)
        neighbor_list(atoms,index,bin_list,10)  # 10 Angstom cutoff radius
    time_nn.stop()

    print("Total:")
    time_total.stop()


if __name__=='__main__':
    main()
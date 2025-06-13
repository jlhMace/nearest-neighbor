from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.neighborlist import NeighborList
from timer import Timer
import numpy as np
import csv

         

def bin_sort(atoms,cutoff):
    '''Sorts atoms into bins based on cutoff
        Input:  atoms:              Atoms object of n atoms,
                cutoff:             neighbor cutoff in Angstroms
        Output: list_index_by_bin:  h*k*l nested list of atomic indexes in each orthagonal bin,
                atom_list:          list of 3D tuple indicating bin location by atomic index,
                bin_num:            3x1 list of number of bins in unit cell'''
    
    ### Initialization ###
    cell_size = atoms.cell.cellpar()  # unit cell size and angles

    ##### 

    ######

    
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
        Input:  index:              index of particular atom in Atoms object,
                atom_list:          3-D nested list of bins with cooresponding atomic indexes,
                list_index_by_bin:  h*k*l nested list of atom indicies by bin coordinate
                bin_num:            3x1 list of dimensional number  of bins in unit cell
        Output: bin_list:           list of atomic indexes in neighboring bins'''
    
    [x,y,z] = atom_list[index]
    [a,b,c] = bin_num
    [x0,x1]=[(x-1)%(a), (x+1)%(a)]
    [y0,y1]=[(y-1)%(b), (y+1)%(b)]
    [z0,z1]=[(z-1)%(c), (z+1)%(c)]
    bin_list = []
    for i in np.unique([x0,x,x1]):
        for j in np.unique([y0,y,y1]):
            for k in np.unique([z0,z,z1]):
                for e in list_index_by_bin[i][j][k]:
                    bin_list.append(e)

    return bin_list


def neighbor_list(atoms,atom_index,bin_list,cutoff):
    '''Reads trajectory file and performs nearest neighbor algorithm
        Input:  atoms:          Atoms object of n atoms,
                atom_index:     index of particular atom in Atoms object,
                bin_list:       list of atomic indexes in neighboring bins
                cutoff:         neighbor cutoff in Angstroms
        Output: neighbor_atoms: Atoms object of nearby atoms'''

    cutoff=cutoff/2

    ### Create Nearest Neighbor list ###
    neighbor_atoms = Atoms()
    for index in bin_list:
        if atoms.get_distances(index,atom_index)<=cutoff:
            neighbor_atoms.append(atoms[index])

    #print(len(neighbor_atoms))
    return neighbor_atoms
    
    # From PyAMFF
    #nl = NeighborList(cutoffs=[cutoff]*len(atoms), sorted=False, self_interaction=False, bothways=True, skin=0.)
    #nl.update(atoms)
    #return nl

def neighbor_list_binless(atoms,atom_index,cutoff):
    '''Reads trajectory file and performs nearest neighbor algorithm without binning
        Input:  atoms:          Atoms object of n atoms,
                atom_index:     index of particular atom in Atoms object,
                cutoff:         neighbor cutoff in Angstroms
        Output: neighbor_atoms: Atoms object of nearby atoms'''

    cutoff=cutoff/2

    ### Create Nearest Neighbor list ###
    neighbor_atoms = Atoms()
    for index in range(len(atoms)):
        if atoms.get_distances(index,atom_index)<=cutoff:
            neighbor_atoms.append(atoms[index])

    #print(len(neighbor_atoms))
    return neighbor_atoms


def run_nn(traj_file,outfile=None):
    '''Runs nn with bins'''
    ### Timing classes for testing ###
    time_total = Timer("total")
    time_bin = Timer("bin")
    time_nn = Timer("nn")
    print("Timing start...")
    time_total.start()

    ### Trajectory file import and setup ###
    traj = Trajectory(traj_file)  # read traj file
    atoms = traj[-1]
    cutoff = 10 # Angstroms

    ### Sorting and neighbor list functions ###
    print("Bin sorting:")
    time_bin.start()
    list_index_by_bin, atoms_list, bin_num = bin_sort(atoms,cutoff)
    time_bin_dat = time_bin.stop()

    print('Atoms, # of bins, atoms/bin')
    n=len(atoms)
    number=bin_num[0]*bin_num[1]*bin_num[2]
    print(n,'   ',number,'     ',n/number)

    print("Nearest neighbor search:")
    time_nn.start()
    for index in range(len(atoms)):
        bin_list = bin_cull(index,atoms_list,list_index_by_bin,bin_num)
        neighbor_list(atoms,index,bin_list,10)  # 10 Angstom cutoff radius
    time_nn_dat = time_nn.stop()

    print("Total:")
    time_total_dat = time_total.stop()

    if outfile:
        # write to csv
        fields = [2,len(atoms),bin_num[0]*bin_num[1]*bin_num[2],time_bin_dat,time_nn_dat,time_total_dat]
        with open(outfile,'a') as f:
            writer = csv.writer(f)
            writer.writerow(fields)


def run_nn_binless(traj_file,outfile=None):
    '''Runs nn without bins'''
    ### Timing classes for testing ###
    time_total = Timer("total")
    time_bin = Timer("bin")
    time_nn = Timer("nn")
    print("Timing start...")
    time_total.start()

    ### Trajectory file import and setup ###
    traj = Trajectory(traj_file)  # read traj file
    atoms = traj[-1]
    cutoff = 10 # Angstroms

    ### Neighbor list functions ###
    print("Nearest neighbor search:")
    time_nn.start()
    for index in range(len(atoms)):
        neighbor_list_binless(atoms,index,10)  # 10 Angstom cutoff radius
    time_nn_dat = time_nn.stop()

    print("Total:")
    time_total_dat = time_total.stop()

    if outfile:
        # write to csv
        fields = [1,len(atoms),time_nn_dat,time_total_dat]
        with open(outfile,'a') as f:
            writer = csv.writer(f)
            writer.writerow(fields)



def main():
    #run_nn('','')
    run_nn_binless('data-trajectory-files/uniform_orthorhombic/10Acutoff-flat-40A.traj','data-graphing/orthorhombic_flat_nobin.csv')




if __name__=='__main__':
    main()
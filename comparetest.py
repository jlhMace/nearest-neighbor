import nn_main as nn
import pickle
from ase.neighborlist import NeighborList, build_neighbor_list, get_connectivity_matrix
from ase.io import Trajectory
import filecmp


def create_reference(trajfile,outfile):
    '''Creates neighborlist file for nearest neighbors of each atom of a trajectory using binless nearest neighbor.'''

    with open(outfile,'wb') as f:
        cutoff = 10/2
        atoms = Trajectory(trajfile)[-1]
        nl = build_neighbor_list(atoms, cutoffs=[cutoff]*len(atoms), sorted=False, self_interaction=False, bothways=True, skin=0.)
        matrix = nl.get_connectivity_matrix()
        print(type(matrix))
        pickle.dump(matrix,f)

def create_neighborlist(trajfile,outfile):
    '''Creates neighborlist file for nearest neighbors of each atom of a trajectory using binless nearest neighbor.'''

    with open(outfile,'w') as f:
        cutoff = 10/2
        atoms = Trajectory(trajfile)[-1]
        nl = nn.run_nn(trajfile)
        matrix = nl.get_connectivity_matrix()
        print(type(matrix))
        pickle.dump(matrix,f)

def compare_files(f1,f2):
    '''Compares two files using filecmp. True if contents are the same, False if different.'''
    return filecmp(f1,f2,shallow=False)


def main():
    trajfile = 'data-testing/10Acutoff-thin-40A.traj'
    outfile = 'data-testing/10Acutoff-thin-40A.pkl'
    outtest = 'data-testing/10Acutoff-thin-40A_test.pkl'
    create_reference(trajfile,outfile)
    #create_neighborlist(trajfile,outtest)


if __name__=='__main__':
    main()
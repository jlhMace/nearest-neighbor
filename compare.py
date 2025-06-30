import nn_main as nn
import pickle
from ase.neighborlist import build_neighbor_list
from ase.io import Trajectory
import filecmp
import numpy as np
import pandas as pd


def create_reference(trajfile,outfile):
    '''Creates neighborlist file for nearest neighbors of each atom of a trajectory using binless nearest neighbor.'''

    with open(outfile,'wb') as f:
        cutoff = 10/2
        atoms = Trajectory(trajfile)[-1]
        nlist = []
        offlist = []
        for i in range(len(atoms)):
            nlist.append([])
            offlist.append([])
            
        nl = build_neighbor_list(atoms, cutoffs=[cutoff]*len(atoms), sorted=False, self_interaction=False, bothways=True, skin=0.)
        for i in range(0,len(atoms)):
            indices, offsets = nl.get_neighbors(i)
            nlist[i].append(np.sort(indices))
            offlist[i].append(offsets)

        print(nlist[0])
        pickle.dump(nlist,f)

def create_neighborlist(trajfile,outfile):
    '''Creates neighborlist file for nearest neighbors of each atom of a trajectory using binless nearest neighbor.'''

    with open(outfile,'wb') as f:
        cutoff = 10
        atoms = Trajectory(trajfile)[-1]
        nlist = []
        offlist = []
        for i in range(len(atoms)):
            nlist.append([])
            offlist.append([])
            
        list_index_by_bin, atoms_list, bin_num = nn.bin_sort(atoms,cutoff)
        for i in range(0,len(atoms)):
            bin_nlist,pointers = nn.bin_cull(i,atoms,atoms_list,list_index_by_bin,bin_num)
            indices, offsets = nn.neighbor_list(bin_nlist,i,cutoff,bin_nlist,pointers)
            nlist[i].append(np.sort(indices))
            offlist[i].append(offsets)

        print(nlist[0])
        pickle.dump(nlist,f)

def compare_files(f1,f2):
    '''Compares two files using filecmp. True if contents are the same, False if different.'''
    return filecmp.cmp(f1,f2,shallow=False)

def view_pickle(file1,file2):
    f1 = pd.read_pickle(file1)
    f2 = pd.read_pickle(file2)
    print(f1[100])
    print(f2[100])


def compare_nlist(trajfile):
    cutoff = 10
    atoms = Trajectory(trajfile)[-1]

    nlist_b, offlist_b = nn.run_nn(trajfile)
    nlist_nb, offlist_nb = nn.run_nn_binless(trajfile)

    for i in range(len(nlist_b)):
        differ = list(set(nlist_b[i]).difference(nlist_nb[i]))
        if np.size(differ)!=0:
            print(f'At index {i}: {differ}')






def main():
    trajfile = 'data-trajectory-files/uniform_cubic/trajprop-10Acutoff-20A.traj'
    outfile = 'data-testing/10Acutoff-thin-20A.pkl'
    outtest = 'data-testing/10Acutoff-thin-20A_test.pkl'
    #create_reference(trajfile,outfile)
    #print('reference ^^')
    #create_neighborlist(trajfile,outtest)
    #print('nlist ^^')
    #print(compare_files(outfile,outtest))
    #print('Done')
    #view_pickle('data-testing/10Acutoff-cubic-10A.pkl','data-testing/10Acutoff-cubic-10A_test.pkl')
    #view_pickle('data-testing/10Acutoff-thin-30A.pkl','data-testing/10Acutoff-thin-30A_test.pkl')
    compare_nlist(trajfile)


if __name__=='__main__':
    main()
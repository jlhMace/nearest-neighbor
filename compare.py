import nn_main as nn
import pickle
from ase.neighborlist import build_neighbor_list
from ase.io import Trajectory
import filecmp
import numpy as np
import pandas as pd


def create_reference(trajfile,outfile):
    '''Creates neighborlist file for nearest neighbors of each atom of a trajectory using binless nearest neighbor.'''

    print('Creating reference Nlist...')

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

        #print(nlist[0])
        pickle.dump(nlist,f)

    print('Reference complete.')


def create_neighborlist(trajfile,outfile):
    '''Creates neighborlist file for nearest neighbors of each atom of a trajectory using binless nearest neighbor.'''
    
    print('Creating Test Nlist...')

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

        #print(nlist[0])
        pickle.dump(nlist,f)
    
    print('Test complete.')


def compare_files(f1,f2):
    '''Compares two files using filecmp. True if contents are the same, False if different.'''
    compare = filecmp.cmp(f1,f2,shallow=False)
    print(f'Are these files the same? {compare}')
    return compare


def view_pickle(file1,file2,index):
    f1 = pd.read_pickle(file1)
    f2 = pd.read_pickle(file2)
    print(f1[index])
    print(f2[index])


def compare_nlist(trajfile,reffile,trajformat='pkl',refformat='pkl'):
    '''Compare the output of a run with the reference nlist,
        from either a .pkl input from the create_neighborlist function,
        or from a fresh run_nn run'''

    if trajformat==('pkl' or 'pickle'):
        nlist_b = pd.read_pickle(trajfile)
    elif trajformat==('traj' or 'Trajectory' or 'trajectory'):
        nlist_b, offlist_b = nn.run_nn(trajfile)
    else:
        print("Trajectory file format not supported, please indicate \'pkl\' or \'traj\'.")

    if refformat==('pkl' or 'pickle'):
        nlist_nb = pd.read_pickle(reffile)
    elif refformat==('traj' or 'Trajectory' or 'trajectory'):
        nlist_nb, offlist_nb = nn.run_nn_binless(reffile)
    else:
        print("Reference file format not supported, please indicate \'pkl\' or \'traj\'.")

    for i in range(len(nlist_b)):
        differ = np.setdiff1d(nlist_b[i],nlist_nb[i])
        if np.size(differ)!=0:
            print(f'At index {i}: {differ}')



def main():
    trajfile = 'data-testing/10Acutoff-cubic-40A.traj'
    outfile = 'data-testing/10Acutoff-cubic-40A.pkl'
    outtest = 'data-testing/10Acutoff-cubic-40A_test.pkl'
    create_reference(trajfile,outfile)
    create_neighborlist(trajfile,outtest)
    compare_files(outfile,outtest)
    #view_pickle('data-testing/10Acutoff-cubic-10A.pkl','data-testing/10Acutoff-cubic-10A_test.pkl',100)
    #view_pickle('data-testing/10Acutoff-thin-30A.pkl','data-testing/10Acutoff-thin-30A_test.pkl',100)
    #compare_nlist(outtest,outfile)  # Compare nonbins to bins
    compare_nlist(outfile,outtest) # Compare bins to nonbins
    print('Done')


if __name__=='__main__':
    main()
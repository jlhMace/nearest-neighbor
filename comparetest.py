import nn_main as nn
import json


def create_reference(trajfile,outfile):
    '''Creates neighborlist file for nearest neighbors of each atom of a trajectory using binless nearest neighbor.'''

    with open(outfile,'w') as f:
        neighborlist = nn.run_nn_binless(trajfile)
        json.dump(neighborlist,f,indent=2)

def create_neighborlist(trajfile,outfile):
    '''Creates neighborlist file for nearest neighbors of each atom of a trajectory using binless nearest neighbor.'''

    with open(outfile,'w') as f:
        neighborlist = nn.run_nn(trajfile)
        json.dump(neighborlist,f,indent=2)


def main():
    trajfile = 'data-trajectory-files/uniform_cubic/trajprop-10Acutoff-10A.traj'
    outfile = 'data-testing/uniform_cubic_10A.json'
    create_reference(trajfile,str(outfile+'.json'))
    create_neighborlist(trajfile,str(outfile+'_bintest.json'))


if __name__=='__main__':
    main()
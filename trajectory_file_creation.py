#### Function to create a trajectory file of uniformly distributed atoms
####    using input unit cell dimensions

from ase.io.trajectory import Trajectory
from ase import Atoms, Atom
from random import uniform
import numpy as np
from math import cos,sin,sqrt,pi,radians


def create_trajectory(filename,cutoff,a,b,c,alpha=90,beta=90,gamma=90):
    '''Create 1-step trajectory test file of uniformly distributed hydrogen atoms
    Input:  a,b,c               dimensions of unit cell
            alpha,beta,gamma    angles of unit cell
    Output: {filename}.traj     Trajectory file created at {filename}.traj'''

    # parameters and density check
    element = 'H'
    num = int((a*b*c)*50/pow(cutoff,3))    # number of total atoms
    bins = int((a*b*c)/pow(cutoff,3))
    print(num,bins,num/bins)

    # Creates Atom object
    atom_list = Atoms(cell=[a,b,c,alpha,beta,gamma])

    # switch to degrees because ASE hates me
    alpha = radians(alpha)
    beta = radians(beta)
    gamma = radians(gamma)

    # Transformation matrix
    cosalp = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma))
    A = np.array([[a,b*cos(gamma),c*cos(beta)],
                [0,b*sin(gamma),-c*sin(beta)*cosalp],
                [0,0,c*sin(beta)*sqrt(1-pow(cosalp,2))]])

    # populates Atoms object
    for e in range(num):
        # generates uniform fractional coordinates
        x = uniform(0,1)
        y = uniform(0,1)
        z = uniform(0,1)
        # Each point is transformed from fractional to orthogonal coordinates
        ortho_coord = np.matmul(A,np.array([x,y,z]))
        atom_list.append(Atom(element,position=ortho_coord))

    # Configures traj file
    traj = Trajectory(f'{filename}', 'w', atoms=atom_list)
    traj.write(atom_list)
    traj.close()

def main():
    #create_trajectory('test2',10,100,100,100,90,90,70)
    pass


if __name__=='__main__':
    main()
#### Function to create a trajectory file of uniformly distributed atoms
####    using input unit cell dimensions

from ase.io.trajectory import Trajectory
from ase import Atoms, Atom
from random import uniform
import numpy as np
from math import cos,sin,sqrt,pi,radians,degrees
from ase.visualize import view


def create_trajectory_uniform(filename,cutoff,a,b,c,alpha=90,beta=90,gamma=90):
    '''Create 1-step trajectory test file of uniformly distributed hydrogen atoms
            of volume based on the orthorhombic shape
        Input:  a,b,c               dimensions of unit cell
                alpha,beta,gamma    angles of unit cell
        Output: {filename}.traj     Trajectory file created at {filename}.traj'''

    # parameters and density check
    element = 'H'
    density = 50    # atoms per bin
    num = int((a*b*c)*density/pow(cutoff,3))    # number of total atoms
    bins = int((a*b*c)/pow(cutoff,3))
    print(num,bins,num/bins)

    # Creates Atom object
    atom_list = Atoms(cell=[a,b,c,alpha,beta,gamma],pbc=True)

    # switch to degrees because ASE hates me
    alpha = radians(alpha)
    beta = radians(beta)
    gamma = radians(gamma)
    
    if a!=b and a!=c:
        print(f'I haven\'t coded for uneven cell sizes yet, defaulting to a=b=c.')

    # define new cell lengths for parallelepiped of same volume in relation to orthorhombic
    volume = a*b*c
    abc = sqrt((pow(volume,2)/(1 + 2*cos(alpha)*cos(beta)*cos(gamma) - cos(alpha)**2 - cos(beta)**2 -cos(gamma)**2))**(1./3.))
    a = abc
    b = abc
    c = abc
    atom_list.set_cell([a,b,c,degrees(alpha),degrees(beta),degrees(gamma)])

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

    print(atom_list.cell.cellpar(),atom_list.cell.volume)
        
    # Configures traj file
    traj = Trajectory(f'{filename}', 'w', atoms=atom_list)
    traj.write(atom_list)
    traj.close()
    #view(atom_list)


def create_custom(filename):
    '''Creates a hardcoded defined system'''
    n=40-39.46785238
    atomsys = Atoms('LiHH', positions=[(13.10486985, 15.95092176, 9.01085498), (15.23491011, 16.41032859, -n), (15.23491011, 16.41032859, 39.46785238)], pbc=True, cell=[40.0, 40.0, 40.0])

    traj = Trajectory(f'{filename}', 'w', atoms=atomsys)
    traj.write(atomsys)
    traj.close()

def main():
    create_trajectory_uniform('data-testing/10Acutoff-cubic-40A.traj',10,40,40,40,90,90,90)
    #create_custom('test2')


if __name__=='__main__':
    main()
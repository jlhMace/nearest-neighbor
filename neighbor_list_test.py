import numpy as np
from ase.neighborlist import NeighborList
from ase import Atoms
import ase

'''
Isolated Python loop_images for grid-search testing
Input: atoms list (single .traj index), atom index a_index
Output: 
'''

class NeighborListsTest:

    def __init__(self, cutoff):
        self.cutoff = cutoff / 2.0

        # Big dictionary containing all data
        self.data = {}

        self.elements = None
        self.index = 1
        self.unitVects = []
        self.vects = []
        self.dists = []
        self.nlist = []
        self.offsets = []
        self.symbols = []
        self.neighborSymbols = []
        self.angles = []
        self.coords = []

        self.a_index = 0
    
    def loop_images(self, atoms, a_index):
        r_cut_2 = self.cutoff # radial cutoff / 2

        for i in range(len(atoms)):
            self.unitVects.append([])
            self.dists.append([])
            self.vects.append([])
            self.nlist.append([])
            self.offsets.append([])
            self.symbols.append([])
            self.neighborSymbols.append([])
            self.angles.append([])

        nl = NeighborList(cutoffs=[r_cut_2]*len(atoms), sorted=False, self_interaction=False, bothways=True, skin=0.)
        nl.update(atoms)
        self.coords.append(atoms.get_positions())
        for i in range(0,len(atoms)):
            indices,offsets = nl.get_neighbors(i)
            self.nlist[a_index].append(np.sort(indices))
            #self.nlist[a_index].append(indices)
            self.offsets[a_index].append(offsets)

        for i in range(0,len(self.nlist[a_index])):
            dists = []
            vects = []
            unit_vects = []
            symbols = []
            neighborSymbols = []
            coords = []
            for j in range(0,len(self.nlist[a_index][i])):
                distance = atoms.get_distance(i, self.nlist[a_index][i][j], mic=True)
                dists.append((distance))
                vector = atoms.get_distance(i, self.nlist[a_index][i][j], mic=True, vector=True)
                vects.append(vector)
                unit_vector = vector / distance
                unit_vects.append(unit_vector)
                neighborSymbols.append(atoms.symbols[self.nlist[a_index][i][j]])
            self.symbols[a_index].append(atoms.symbols[i])
            self.dists[a_index].append(dists)
            self.vects[a_index].append(vects)
            self.unitVects[a_index].append(unit_vects)
            self.neighborSymbols[a_index].append(neighborSymbols)

        # Generate angle data
        for aa in range(0,len(self.nlist[a_index])):
            angles = []
            for aaa in range(0,len(self.nlist[a_index][aa])):
                for aaaa in range(aaa+1,len(self.nlist[a_index][aa])):
                    angs = atoms.get_angle(self.nlist[a_index][aa][aaa], aa, self.nlist[a_index][aa][aaaa], mic=True)
                    angles.append(angs)

            self.angles[a_index].append(angles)
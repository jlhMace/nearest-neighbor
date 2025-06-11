from ase.io.trajectory import Trajectory
from neighbor_list import NeighborLists
from neighbor_list_test import NeighborListsTest
from timer import Timer

traj = Trajectory('100000atoms-100A.traj')
atoms = traj[-1]
nl = NeighborLists(1.9)

testTime = Timer("test")
testTime.start()
nl.calculate(traj)
#nl.loop_images(atoms,0)
testTime.stop()
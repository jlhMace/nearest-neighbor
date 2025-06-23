import io
import os
import unittest
#import nn_main as nn
from ase.io import Trajectory
import re
from comparetest import create_neighborlist, compare_files


class TestReferenceFiles():

    def __init__(self):
        self.directory = 'data-testing'
    
    def file_list(self):
        filelist = []
        for file in os.listdir(self.directory):
            filename = os.fsdecode(file)
            if filename.endswith('.traj'):
                filelist.append(filename)
        return filelist
    
    def check_file_exists(self,filename):
        filename = re.sub('.traj','.pkl',filename)
        file = str(self.directory + '/' + filename)
        with io.open(file) as f:
            assert os.path.isfile(file)

    def test_gen_file_exists(self):
        filelist = self.file_list()
        for e in filelist:
            yield self.check_file_exists, e


class TestBinAccuracy():

    def __init__(self):
        self.directory = 'data-testing'

    def file_list(self):
        filelist = []
        for file in os.listdir(self.directory):
            filename = os.fsdecode(file)
            if filename.endswith('.traj'):
                filelist.append(self.directory + '/' + filename)
        return filelist
    
    def setUp(self):
        filelist = self.file_list()
        for trajfile in filelist:
            trajfile = str(self.directory + '/' + trajfile)
            outfile = re.sub('.traj','_test.pkl',trajfile)
            create_neighborlist(trajfile,outfile)

    def test_gen_file_compare(self):
        filelist = self.file_list()
        for e in filelist:
            f1 = re.sub('.traj','.pkl',e)
            f2 = re.sub('.pkl','_test.pkl',f1)
            print(f1,f2,self.compare_output(f1,f2)==True)
            yield self.compare_output, f1, f2

    def compare_output(self,f1,f2):
        compare_files(f1,f2)

    
            


if __name__=='__main__':
    unittest.main()

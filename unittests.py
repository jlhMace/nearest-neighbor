import io
import os
import unittest
#import nn_main as nn
from ase.io import Trajectory
import re
from compare import create_neighborlist, compare_files
import pytest

class TestReferenceFiles():
    
    @pytest.fixture
    def filelist(self):
        directory = 'data-testing'
        filelist = []
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith('.traj'):
                filelist.append(filename)
        return filelist
    
    def check_file_exists(self,filelist):
        for filename in filelist:
            filename = re.sub('.traj','.pkl',filename)
            file = str('data-testing/' + filename)
            with io.open(file) as f:
                assert os.path.isfile(file)


class TestBinAccuracy():

    @classmethod
    def setup(self):
        self.directory = 'data-testing'

    @pytest.fixture
    def file_list(self):
        filelist = []
        for file in os.listdir(self.directory):
            filename = os.fsdecode(file)
            if filename.endswith('.traj'):
                filelist.append(self.directory + '/' + filename)
        return filelist
    
    @pytest.mark.skip(reason="WIP")
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
            assert self.compare_output(f1, f2)

    @pytest.fixture
    def compare_output(self,f1,f2):
        compare_files(f1,f2)

    
            


if __name__=='__main__':
    unittest.main()

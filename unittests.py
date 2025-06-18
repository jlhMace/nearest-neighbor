import io
import os
import unittest
#import nn_main as nn
from ase.io import Trajectory
import re


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
        filename = re.sub('.traj','.txt',filename)
        file = str(self.directory + '/' + filename)
        with io.open(file) as f:
            assert os.path.isfile(file)

    def test_gen_file_exists(self):
        filelist = self.file_list()
        for e in filelist:
            yield self.check_file_exists, e


class TestBinAccuracy():

    def setUp(self):
        pass

    @unittest.skip("Base test case")
    def test_nn_base(self,test_path,ref_path):
        with io.open(test_path) as tp:
            with io.open(ref_path) as rp:
                self.assertEqual(list(tp),list(rp))

    def test_nn_cubic(self):
            #nn.run_nn()
            self.test_nn_base()
            


if __name__=='__main__':
    unittest.main()

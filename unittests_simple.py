import io
import os
import unittest
#import nn_main as nn
from ase.io import Trajectory
import re
from comparetest import create_neighborlist, compare_files
import pytest
import filecmp

def test_check_file_same_10Acubic():
    file = str('data-testing/10Acutoff-cubic-10A.pkl')
    with io.open(file) as f:
        f2 = str(re.sub('.pkl','_test.pkl',file))
        assert filecmp.cmp(file,f2,shallow=False)

def test_check_file_same_30Athin():
    file = str('data-testing/10Acutoff-thin-30A.pkl')
    with io.open(file) as f:
        f2 = str(re.sub('.pkl','_test.pkl',file))
        assert filecmp.cmp(file,f2,shallow=False)
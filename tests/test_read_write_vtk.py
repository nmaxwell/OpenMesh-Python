import unittest
import openmesh
import os

import numpy as np

class ReadWriteVTK(unittest.TestCase):

    def setUp(self):
        if not os.path.exists('OutFiles'):
            os.makedirs('OutFiles')

    def test_read_write_read(self):
        mesh1 = openmesh.read_trimesh("TestFiles/cube-minimal.obj")
        openmesh.write_mesh('OutFiles/test_read_write_read.vtk', mesh1)
        # since there is no vtk reader, all we can do is check that
        # write_mesh() doesn't throw

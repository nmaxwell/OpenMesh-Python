import unittest
import openmesh

import numpy as np

class Python(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        openmesh.read_mesh(self.mesh, 'TestFiles/cube-minimal.obj')

    def test_vertex_property_array(self):
        prop = openmesh.VPropHandle()
        self.mesh.add_property(prop)
        arr1 = np.random.rand(self.mesh.n_vertices(), 10)
        self.mesh.set_property_array(prop, arr1)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1, arr2))

    def test_halfedge_property_array(self):
        prop = openmesh.HPropHandle()
        self.mesh.add_property(prop)
        arr1 = np.random.rand(self.mesh.n_halfedges(), 10)
        self.mesh.set_property_array(prop, arr1)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1, arr2))

    def test_edge_property_array(self):
        prop = openmesh.EPropHandle()
        self.mesh.add_property(prop)
        arr1 = np.random.rand(self.mesh.n_edges(), 10)
        self.mesh.set_property_array(prop, arr1)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1, arr2))

    def test_face_property_array(self):
        prop = openmesh.FPropHandle()
        self.mesh.add_property(prop)
        arr1 = np.random.rand(self.mesh.n_faces(), 10)
        self.mesh.set_property_array(prop, arr1)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1, arr2))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

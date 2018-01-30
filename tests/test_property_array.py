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
        # c_contiguous
        arr1 = np.random.rand(self.mesh.n_vertices(), 10)
        self.mesh.set_property_array(prop, arr1)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1, arr2))
        for vh in self.mesh.vertices():
            arr3 = self.mesh.property(prop, vh)
            self.assertTrue(np.allclose(arr1[vh.idx()], arr3))
        # f_contiguous
        arr1 = np.random.rand(10, self.mesh.n_vertices())
        self.mesh.set_property_array(prop, arr1.T)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1.T, arr2))
        for vh in self.mesh.vertices():
            arr3 = self.mesh.property(prop, vh)
            self.assertTrue(np.allclose(arr1.T[vh.idx()], arr3))

    def test_halfedge_property_array(self):
        prop = openmesh.HPropHandle()
        self.mesh.add_property(prop)
        # c_contiguous
        arr1 = np.random.rand(self.mesh.n_halfedges(), 10)
        self.mesh.set_property_array(prop, arr1)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1, arr2))
        for hh in self.mesh.halfedges():
            arr3 = self.mesh.property(prop, hh)
            self.assertTrue(np.allclose(arr1[hh.idx()], arr3))
        # f_contiguous
        arr1 = np.random.rand(10, self.mesh.n_halfedges())
        self.mesh.set_property_array(prop, arr1.T)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1.T, arr2))
        for hh in self.mesh.halfedges():
            arr3 = self.mesh.property(prop, hh)
            self.assertTrue(np.allclose(arr1.T[hh.idx()], arr3))

    def test_edge_property_array(self):
        prop = openmesh.EPropHandle()
        self.mesh.add_property(prop)
        # c_contiguous
        arr1 = np.random.rand(self.mesh.n_edges(), 10)
        self.mesh.set_property_array(prop, arr1)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1, arr2))
        for eh in self.mesh.edges():
            arr3 = self.mesh.property(prop, eh)
            self.assertTrue(np.allclose(arr1[eh.idx()], arr3))
        # f_contiguous
        arr1 = np.random.rand(10, self.mesh.n_edges())
        self.mesh.set_property_array(prop, arr1.T)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1.T, arr2))
        for eh in self.mesh.edges():
            arr3 = self.mesh.property(prop, eh)
            self.assertTrue(np.allclose(arr1.T[eh.idx()], arr3))

    def test_face_property_array(self):
        prop = openmesh.FPropHandle()
        self.mesh.add_property(prop)
        # c_contiguous
        arr1 = np.random.rand(self.mesh.n_faces(), 10)
        self.mesh.set_property_array(prop, arr1)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1, arr2))
        for fh in self.mesh.faces():
            arr3 = self.mesh.property(prop, fh)
            self.assertTrue(np.allclose(arr1[fh.idx()], arr3))
        # f_contiguous
        arr1 = np.random.rand(10, self.mesh.n_faces())
        self.mesh.set_property_array(prop, arr1.T)
        arr2 = self.mesh.property_array(prop)
        self.assertTrue(np.allclose(arr1.T, arr2))
        for fh in self.mesh.faces():
            arr3 = self.mesh.property(prop, fh)
            self.assertTrue(np.allclose(arr1.T[fh.idx()], arr3))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

import unittest
import openmesh

import numpy as np

class Python(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        openmesh.read_mesh(self.mesh, 'TestFiles/cube-minimal.obj')

    def test_vertex_property_array(self):
        # c_contiguous
        arr1 = np.random.rand(self.mesh.n_vertices(), 10)
        self.mesh.set_vertex_property_array('random', arr1)
        arr2 = self.mesh.vertex_property_array('random')
        self.assertTrue(np.allclose(arr1, arr2))
        for vh in self.mesh.vertices():
            arr3 = self.mesh.vertex_property('random', vh)
            self.assertTrue(np.allclose(arr1[vh.idx()], arr3))
        # f_contiguous
        arr1 = np.random.rand(10, self.mesh.n_vertices())
        self.mesh.set_vertex_property_array('random', arr1.T)
        arr2 = self.mesh.vertex_property_array('random')
        self.assertTrue(np.allclose(arr1.T, arr2))
        for vh in self.mesh.vertices():
            arr3 = self.mesh.vertex_property('random', vh)
            self.assertTrue(np.allclose(arr1.T[vh.idx()], arr3))

    def test_halfedge_property_array(self):
        # c_contiguous
        arr1 = np.random.rand(self.mesh.n_halfedges(), 10)
        self.mesh.set_halfedge_property_array('random', arr1)
        arr2 = self.mesh.halfedge_property_array('random')
        self.assertTrue(np.allclose(arr1, arr2))
        for hh in self.mesh.halfedges():
            arr3 = self.mesh.halfedge_property('random', hh)
            self.assertTrue(np.allclose(arr1[hh.idx()], arr3))
        # f_contiguous
        arr1 = np.random.rand(10, self.mesh.n_halfedges())
        self.mesh.set_halfedge_property_array('random', arr1.T)
        arr2 = self.mesh.halfedge_property_array('random')
        self.assertTrue(np.allclose(arr1.T, arr2))
        for hh in self.mesh.halfedges():
            arr3 = self.mesh.halfedge_property('random', hh)
            self.assertTrue(np.allclose(arr1.T[hh.idx()], arr3))

    def test_edge_property_array(self):
        # c_contiguous
        arr1 = np.random.rand(self.mesh.n_edges(), 10)
        self.mesh.set_edge_property_array('random', arr1)
        arr2 = self.mesh.edge_property_array('random')
        self.assertTrue(np.allclose(arr1, arr2))
        for eh in self.mesh.edges():
            arr3 = self.mesh.edge_property('random', eh)
            self.assertTrue(np.allclose(arr1[eh.idx()], arr3))
        # f_contiguous
        arr1 = np.random.rand(10, self.mesh.n_edges())
        self.mesh.set_edge_property_array('random', arr1.T)
        arr2 = self.mesh.edge_property_array('random')
        self.assertTrue(np.allclose(arr1.T, arr2))
        for eh in self.mesh.edges():
            arr3 = self.mesh.edge_property('random', eh)
            self.assertTrue(np.allclose(arr1.T[eh.idx()], arr3))

    def test_face_property_array(self):
        # c_contiguous
        arr1 = np.random.rand(self.mesh.n_faces(), 10)
        self.mesh.set_face_property_array('random', arr1)
        arr2 = self.mesh.face_property_array('random')
        self.assertTrue(np.allclose(arr1, arr2))
        for fh in self.mesh.faces():
            arr3 = self.mesh.face_property('random', fh)
            self.assertTrue(np.allclose(arr1[fh.idx()], arr3))
        # f_contiguous
        arr1 = np.random.rand(10, self.mesh.n_faces())
        self.mesh.set_face_property_array('random', arr1.T)
        arr2 = self.mesh.face_property_array('random')
        self.assertTrue(np.allclose(arr1.T, arr2))
        for fh in self.mesh.faces():
            arr3 = self.mesh.face_property('random', fh)
            self.assertTrue(np.allclose(arr1.T[fh.idx()], arr3))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

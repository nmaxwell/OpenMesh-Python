import unittest
import openmesh

import numpy as np

class Python(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        openmesh.read_mesh(self.mesh, 'TestFiles/cube-minimal.obj')

    def test_face_vertex_indices_trimesh(self):
        indices1 = self.mesh.face_vertex_indices()
        indices2 = self.mesh.fv_indices()
        for fh in self.mesh.faces():
            fv_it = self.mesh.fv(fh)
            vh1 = next(fv_it)
            vh2 = next(fv_it)
            vh3 = next(fv_it)
            self.assertEqual(indices1[fh.idx(), 0], vh1.idx())
            self.assertEqual(indices1[fh.idx(), 1], vh2.idx())
            self.assertEqual(indices1[fh.idx(), 2], vh3.idx())
            self.assertEqual(indices2[fh.idx(), 0], vh1.idx())
            self.assertEqual(indices2[fh.idx(), 1], vh2.idx())
            self.assertEqual(indices2[fh.idx(), 2], vh3.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_faces(), 3))
        self.assertEqual(indices2.shape, (self.mesh.n_faces(), 3))

    def test_face_vertex_indices_polymesh(self):
        self.mesh = openmesh.PolyMesh()
        vh1 = self.mesh.add_vertex(np.array([0, 0, 0]))
        vh2 = self.mesh.add_vertex(np.array([1, 0, 0]))
        vh3 = self.mesh.add_vertex(np.array([1, 1, 0]))
        vh4 = self.mesh.add_vertex(np.array([0, 1, 0]))
        vh5 = self.mesh.add_vertex(np.array([2, 0, 0]))
        vh6 = self.mesh.add_vertex(np.array([3, 0, 0]))
        vh7 = self.mesh.add_vertex(np.array([3, 1, 0]))
        self.mesh.add_face(vh1, vh2, vh3, vh4)
        self.mesh.add_face(vh5, vh6, vh7)

        # Test setup:
        #  3 === 2         6
        #  |     |       / |
        #  |     |      /  |
        #  |     |     /   |
        #  0 === 1   4 === 5

        indices1 = self.mesh.face_vertex_indices()
        indices2 = self.mesh.fv_indices()
        correct = np.array([[0, 1, 2, 3], [4, 5, 6, -1]])
        self.assertTrue(np.array_equal(indices1, correct))
        self.assertTrue(np.array_equal(indices2, correct))

    def test_edge_vertex_indices(self):
        indices1 = self.mesh.edge_vertex_indices()
        indices2 = self.mesh.ev_indices()
        for eh in self.mesh.edges():
            heh = self.mesh.halfedge_handle(eh, 0)
            vh1 = self.mesh.from_vertex_handle(heh)
            vh2 = self.mesh.to_vertex_handle(heh)
            self.assertEqual(indices1[eh.idx(), 0], vh1.idx())
            self.assertEqual(indices1[eh.idx(), 1], vh2.idx())
            self.assertEqual(indices2[eh.idx(), 0], vh1.idx())
            self.assertEqual(indices2[eh.idx(), 1], vh2.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_edges(), 2))
        self.assertEqual(indices2.shape, (self.mesh.n_edges(), 2))

    def test_halfedge_vertex_indices(self):
        indices1 = self.mesh.halfedge_vertex_indices()
        indices2 = self.mesh.hv_indices()
        for heh in self.mesh.halfedges():
            vh1 = self.mesh.from_vertex_handle(heh)
            vh2 = self.mesh.to_vertex_handle(heh)
            self.assertEqual(indices1[heh.idx(), 0], vh1.idx())
            self.assertEqual(indices1[heh.idx(), 1], vh2.idx())
            self.assertEqual(indices2[heh.idx(), 0], vh1.idx())
            self.assertEqual(indices2[heh.idx(), 1], vh2.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_halfedges(), 2))
        self.assertEqual(indices2.shape, (self.mesh.n_halfedges(), 2))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

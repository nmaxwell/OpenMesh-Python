import unittest
import openmesh

import numpy as np

class Python(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.read_trimesh('TestFiles/cube-minimal.obj')

    def test_vertex_property_list(self):
        self.assertFalse(self.mesh.has_vertex_property('random'))
        res = self.mesh.vertex_property('random')
        self.assertTrue(all(val is None for val in res))
        rnd = [np.random.rand(10) for h in self.mesh.vertices()]
        self.mesh.set_vertex_property('random', rnd)
        res = self.mesh.vertex_property('random')
        self.assertTrue(np.allclose(np.stack(rnd), np.stack(res)))
        for i in range(self.mesh.n_vertices()):
            rnd[i][:] = np.random.rand(10)
        res = self.mesh.vertex_property('random')
        self.assertTrue(np.allclose(np.stack(rnd), np.stack(res)))

    def test_halfedge_property_list(self):
        self.assertFalse(self.mesh.has_halfedge_property('random'))
        res = self.mesh.halfedge_property('random')
        self.assertTrue(all(val is None for val in res))
        rnd = [np.random.rand(10) for h in self.mesh.halfedges()]
        self.mesh.set_halfedge_property('random', rnd)
        res = self.mesh.halfedge_property('random')
        self.assertTrue(np.allclose(np.stack(rnd), np.stack(res)))
        for i in range(self.mesh.n_halfedges()):
            rnd[i][:] = np.random.rand(10)
        res = self.mesh.halfedge_property('random')
        self.assertTrue(np.allclose(np.stack(rnd), np.stack(res)))

    def test_edge_property_list(self):
        self.assertFalse(self.mesh.has_edge_property('random'))
        res = self.mesh.edge_property('random')
        self.assertTrue(all(val is None for val in res))
        rnd = [np.random.rand(10) for h in self.mesh.edges()]
        self.mesh.set_edge_property('random', rnd)
        res = self.mesh.edge_property('random')
        self.assertTrue(np.allclose(np.stack(rnd), np.stack(res)))
        for i in range(self.mesh.n_edges()):
            rnd[i][:] = np.random.rand(10)
        res = self.mesh.edge_property('random')
        self.assertTrue(np.allclose(np.stack(rnd), np.stack(res)))

    def test_face_property_list(self):
        self.assertFalse(self.mesh.has_face_property('random'))
        res = self.mesh.face_property('random')
        self.assertTrue(all(val is None for val in res))
        rnd = [np.random.rand(10) for h in self.mesh.faces()]
        self.mesh.set_face_property('random', rnd)
        res = self.mesh.face_property('random')
        self.assertTrue(np.allclose(np.stack(rnd), np.stack(res)))
        for i in range(self.mesh.n_faces()):
            rnd[i][:] = np.random.rand(10)
        res = self.mesh.face_property('random')
        self.assertTrue(np.allclose(np.stack(rnd), np.stack(res)))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

import unittest
import openmesh

import numpy as np

class TriMeshCirculatorVertexEdge(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

        # Add some vertices
        self.vhandle = []

        self.vhandle.append(self.mesh.add_vertex(np.array([0, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([1, 0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([0,-1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2,-1, 0])))
        
        # Add four faces
        self.mesh.add_face(self.vhandle[0], self.vhandle[1], self.vhandle[2])
        self.mesh.add_face(self.vhandle[1], self.vhandle[3], self.vhandle[4])
        self.mesh.add_face(self.vhandle[0], self.vhandle[3], self.vhandle[1])
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[4])

        '''
        Test setup:
            0 ==== 2
            |\  0 /|
            | \  / |
            |2  1 3|
            | /  \ |
            |/  1 \|
            3 ==== 4
        '''

    def test_vertex_edge_iter_without_holes_increment(self):
        # Iterate around vertex 1 at the middle
        ve_it = openmesh.VertexEdgeIter(self.mesh, self.vhandle[1])
        self.assertEqual(next(ve_it).idx(), 5)
        self.assertEqual(next(ve_it).idx(), 3)
        self.assertEqual(next(ve_it).idx(), 0)
        self.assertEqual(next(ve_it).idx(), 1)
        with self.assertRaises(StopIteration):
            next(ve_it)
        
    def test_vertex_edge_iter_boundary_increment(self):
        # Iterate around vertex 2 at the boundary
        ve_it = openmesh.VertexEdgeIter(self.mesh, self.vhandle[2])
        self.assertEqual(next(ve_it).idx(), 7)
        self.assertEqual(next(ve_it).idx(), 1)
        self.assertEqual(next(ve_it).idx(), 2)
        with self.assertRaises(StopIteration):
            next(ve_it)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshCirculatorVertexEdge)
    unittest.TextTestRunner(verbosity=2).run(suite)

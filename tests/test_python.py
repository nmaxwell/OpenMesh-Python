import unittest
import openmesh

import numpy as np

class Python(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

        # Add some vertices
        self.vhandle = []

        self.vhandle.append(self.mesh.add_vertex(np.array([0, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([1, 0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([0,-1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2,-1, 0])))

        # Add four faces using Python lists
        vertex_list = [self.vhandle[0], self.vhandle[1], self.vhandle[2]]
        self.mesh.add_face(vertex_list)
        vertex_list = [self.vhandle[1], self.vhandle[3], self.vhandle[4]]
        self.mesh.add_face(vertex_list)
        vertex_list = [self.vhandle[0], self.vhandle[3], self.vhandle[1]]
        self.mesh.add_face(vertex_list)
        vertex_list = [self.vhandle[2], self.vhandle[1], self.vhandle[4]]
        self.mesh.add_face(vertex_list)

        # Test setup:
        #  0 ==== 2
        #  |\  0 /|
        #  | \  / |
        #  |2  1 3|
        #  | /  \ |
        #  |/  1 \|
        #  3 ==== 4

        self.assertEqual(self.mesh.n_faces(), 4)

    def test_python_iterator(self):
        # Iterate over all vertices
        indices = [0, 1, 2, 3, 4]
        for v, idx in zip(self.mesh.vertices(), indices):
            self.assertEqual(v.idx(), idx)

    def test_python_circulator(self):
        # Iterate around vertex 1 at the middle
        indices = [4, 3, 0, 2]
        for v, idx in zip(self.mesh.vv(self.vhandle[1]), indices):
             self.assertEqual(v.idx(), idx)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

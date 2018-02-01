import unittest
import openmesh

import numpy as np

class TriMeshCirculatorFaceEdge(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()

        # Add some vertices
        self.vhandle = []

        self.vhandle.append(self.mesh.add_vertex(np.array([0, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([1, 0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([3, 0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([4, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2,-1, 0])))
        
        # Add four faces
        self.mesh.add_face(self.vhandle[0], self.vhandle[1], self.vhandle[2])
        self.mesh.add_face(self.vhandle[2], self.vhandle[1], self.vhandle[3])
        self.mesh.add_face(self.vhandle[2], self.vhandle[3], self.vhandle[4])
        self.mesh.add_face(self.vhandle[1], self.vhandle[5], self.vhandle[3])

        '''
        Test setup:
            0 ------ 2 ------ 4
             \      / \      /
              \  0 /   \  2 /
               \  /  1  \  /
                1 ------- 3
                 \       /
                  \  3  /
                   \   /
                    \ /
                     5
        '''

    def test_face_edge_iter_without_holes_increment(self):
        # Iterate around face 1 at the middle
        fe_it = openmesh.FaceEdgeIter(self.mesh, self.mesh.face_handle(1))
        self.assertEqual(next(fe_it).idx(), 4)
        self.assertEqual(next(fe_it).idx(), 1)
        self.assertEqual(next(fe_it).idx(), 3)
        with self.assertRaises(StopIteration):
            next(fe_it)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshCirculatorFaceEdge)
    unittest.TextTestRunner(verbosity=2).run(suite)

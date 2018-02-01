import unittest
import openmesh

import numpy as np

class TriMeshCirculatorFaceVertex(unittest.TestCase):

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
        self.fh0 = self.mesh.add_face(self.vhandle[0], self.vhandle[1], self.vhandle[2])
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

    def test_face_vertex_iter_without_increment(self):
        self.assertEqual(self.fh0.idx(), 0)
    
        # Iterate around face 0 at the top
        fv_it = openmesh.FaceVertexIter(self.mesh, self.fh0)
        self.assertEqual(next(fv_it).idx(), 0)
        self.assertEqual(next(fv_it).idx(), 1)
        self.assertEqual(next(fv_it).idx(), 2)
        with self.assertRaises(StopIteration):
            next(fv_it)
        

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TriMeshCirculatorFaceVertex)
    unittest.TextTestRunner(verbosity=2).run(suite)

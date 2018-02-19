import unittest
import openmesh

import numpy as np

class TrimeshCirculatorFaceFace(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        self.vhandle = []
    
    def test_face_face_iter_with_holes(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(np.array([0, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([1, 0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([3, 0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([4, 1, 0])))

        # Add three faces
        face_vhandles = []

        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)

        face_vhandles = []

        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)

        # Test setup:
        #
        # 0 ------ 2 ------ 4
        #  \      / \      /
        #   \  0 /   \  2 /
        #    \  /  1  \  /
        #     1 ------- 3

        ff_it = self.mesh.ff(self.mesh.face_handle(1))

        self.assertEqual(next(ff_it).idx(), 2)
        self.assertEqual(next(ff_it).idx(), 0)
        with self.assertRaises(StopIteration):
            next(ff_it)

    def test_face_face_iter_without_holes(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(np.array([0,  1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([1,  0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2,  1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([3,  0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([4,  1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([2, -1, 0])))
        
        # Add four faces
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[2])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[4])
        self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[1])
        face_vhandles.append(self.vhandle[5])
        face_vhandles.append(self.vhandle[3])
        self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #
        # 0 ------ 2 ------ 4
        #  \      / \      /
        #   \  0 /   \  2 /
        #    \  /  1  \  /
        #     1 ------- 3
        #      \       /
        #       \  3  /
        #        \   /
        #         \ /
        #          5

        ff_it = self.mesh.ff(self.mesh.face_handle(1))
        
        self.assertEqual(next(ff_it).idx(), 2)
        self.assertEqual(next(ff_it).idx(), 0)
        self.assertEqual(next(ff_it).idx(), 3)
        with self.assertRaises(StopIteration):
            next(ff_it)

    def test_face_face_iter_without_holes(self):
        # Add some vertices
        self.vhandle.append(self.mesh.add_vertex(np.array([0, 0, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([1, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([1, 1, 0])))
        self.vhandle.append(self.mesh.add_vertex(np.array([1, 0, 0])))
        
        # Add two faces
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[2])
        face_vhandles.append(self.vhandle[1])
        fh1 = self.mesh.add_face(face_vhandles)
        
        face_vhandles = []
        
        face_vhandles.append(self.vhandle[0])
        face_vhandles.append(self.vhandle[3])
        face_vhandles.append(self.vhandle[2])
        fh2 = self.mesh.add_face(face_vhandles)
        
        # Test setup:
        #
        #  1 -------- 2
        #  | f0  /    |
        #  |    / f1  |
        #  0 -------- 3
        
        # Check setup
        self.assertEqual(self.mesh.n_vertices(), 4)
        self.assertEqual(self.mesh.n_faces(), 2)

        face_iter = self.mesh.ff(fh1)

        # Get the face via the handle
        faceHandle1 = next(face_iter)
        self.assertEqual(faceHandle1.idx(), 1)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TrimeshCirculatorFaceFace)
    unittest.TextTestRunner(verbosity=2).run(suite)

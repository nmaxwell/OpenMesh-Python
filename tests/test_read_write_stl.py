import unittest
import openmesh
import os

import numpy as np

class ReadWriteSTL(unittest.TestCase):

    def setUp(self):
        if not os.path.exists('OutFiles'):
            os.makedirs('OutFiles')

    def test_read_write_read(self):
        mesh1 = openmesh.read_trimesh("TestFiles/cube-minimal.stl")
        openmesh.write_mesh('OutFiles/test_read_write_read.stl', mesh1)
        mesh2 = openmesh.read_trimesh('OutFiles/test_read_write_read.stl')
        self.assertTrue(np.allclose(mesh1.points(), mesh2.points()))
        self.assertTrue(np.array_equal(mesh1.face_vertex_indices(), mesh2.face_vertex_indices()))

    def test_load_simple_stl_file(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube1.stl")
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)
    
    def test_load_simple_stl_file_with_normals(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube1.stl", face_normal=True)
        
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[0], -0.038545)
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[1], -0.004330)
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[2],  0.999247)
        
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)
        
        self.mesh.release_face_normals()

    def test_load_simple_stl_binary_file(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube1Binary.stl")
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)

    def test_load_simple_stl_binary_file_with_normals(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube1Binary.stl", face_normal=True, binary=True)
        
        self.assertTrue(self.mesh.has_face_normals())
        self.assertFalse(self.mesh.has_vertex_normals())
        
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[0], -0.038545, 5)
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[1], -0.004330, 5)
        self.assertAlmostEqual(self.mesh.normal(self.mesh.face_handle(0))[2],  0.999247, 5)
        
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)

        self.mesh.release_face_normals()

    def test_read_nonexistent_stl(self):
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_trimesh("TestFiles/nonexistent.stl")
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_polymesh("TestFiles/nonexistent.stl")


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteSTL)
    unittest.TextTestRunner(verbosity=2).run(suite)

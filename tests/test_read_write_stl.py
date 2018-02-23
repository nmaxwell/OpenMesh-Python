import unittest
import openmesh

class ReadWriteSTL(unittest.TestCase):

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


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteSTL)
    unittest.TextTestRunner(verbosity=2).run(suite)

import unittest
import openmesh
import os

import numpy as np

class ReadWriteOM(unittest.TestCase):

    def setUp(self):
        if not os.path.exists('OutFiles'):
            os.makedirs('OutFiles')

    def test_read_write_read(self):
        mesh1 = openmesh.read_trimesh("TestFiles/cube-minimal.obj")
        openmesh.write_mesh('OutFiles/test_read_write_read.om', mesh1)
        mesh2 = openmesh.read_trimesh('OutFiles/test_read_write_read.om')
        self.assertTrue(np.allclose(mesh1.points(), mesh2.points()))
        self.assertTrue(np.array_equal(mesh1.face_vertex_indices(), mesh2.face_vertex_indices()))

    def test_load_simple_om_force_vertex_colors_although_not_available(self):
        with self.assertRaises(RuntimeError):
            openmesh.read_trimesh("TestFiles/cube-minimal.om", vertex_color=True)

    def test_load_simple_om_with_texcoords(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal-texCoords.om", vertex_tex_coord=True)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(0))[0], 10.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(0))[1], 10.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(2))[0], 6.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(2))[1], 6.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(4))[0], 9.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(4))[1], 9.0)
        
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(7))[0], 12.0)
        self.assertEqual(self.mesh.texcoord2D(self.mesh.vertex_handle(7))[1], 12.0)
        
        self.assertFalse(self.mesh.has_vertex_normals())
        self.assertTrue(self.mesh.has_vertex_texcoords1D())
        self.assertTrue(self.mesh.has_vertex_texcoords2D())
        self.assertTrue(self.mesh.has_vertex_texcoords3D())
        self.assertFalse(self.mesh.has_vertex_colors())
        
        self.mesh.release_vertex_texcoords2D()

    def test_load_simple_om_with_vertex_colors(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal-vertexColors.om", vertex_color=True)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[2], 1.0)
        
        self.assertFalse(self.mesh.has_vertex_normals())
        self.assertFalse(self.mesh.has_vertex_texcoords1D())
        self.assertFalse(self.mesh.has_vertex_texcoords2D())
        self.assertFalse(self.mesh.has_vertex_texcoords3D())
        self.assertTrue(self.mesh.has_vertex_colors())
        
        self.mesh.release_vertex_colors()

    def test_write_triangle(self):
        self.mesh = openmesh.TriMesh()
        
        # Generate data
        v1 = self.mesh.add_vertex(np.array([1.0, 0.0, 0.0]))
        v2 = self.mesh.add_vertex(np.array([0.0, 1.0, 0.0]))
        v3 = self.mesh.add_vertex(np.array([0.0, 0.0, 1.0]))
        self.mesh.add_face(v1, v2, v3)
        
        # Save
        filename = "OutFiles/triangle-minimal.om"
        openmesh.write_mesh(filename, self.mesh)
        
        # Load
        self.mesh = openmesh.read_trimesh(filename)
        
        # Compare
        self.assertEqual(self.mesh.n_vertices(), 3)
        self.assertEqual(self.mesh.n_edges(), 3)
        self.assertEqual(self.mesh.n_faces(), 1)
        
        self.assertTrue(np.allclose(self.mesh.point(v1), np.array([1.0, 0.0, 0.0])))
        self.assertTrue(np.allclose(self.mesh.point(v2), np.array([0.0, 1.0, 0.0])))
        self.assertTrue(np.allclose(self.mesh.point(v3), np.array([0.0, 0.0, 1.0])))
        
        # Cleanup
        os.remove(filename)

    def test_write_triangle_vertex_integer_color(self):
        self.mesh = openmesh.TriMesh()
        
        # Generate data
        v1 = self.mesh.add_vertex(np.array([1.0, 0.0, 0.0]))
        v2 = self.mesh.add_vertex(np.array([0.0, 1.0, 0.0]))
        v3 = self.mesh.add_vertex(np.array([0.0, 0.0, 1.0]))
        self.mesh.add_face(v1, v2, v3)
        
        c1 = np.array([0.00, 0.00, 0.50, 1.00])
        c2 = np.array([0.25, 0.00, 0.00, 1.00])
        c3 = np.array([0.00, 0.75, 0.00, 1.00])
    
        self.mesh.set_color(v1, c1)
        self.mesh.set_color(v2, c2)
        self.mesh.set_color(v3, c3)
            
        # Save
        filename = "OutFiles/triangle-minimal-ColorsPerVertex.om"
        openmesh.write_mesh(filename, self.mesh, vertex_color=True, color_float=True)
            
        self.mesh.release_vertex_colors()
            
        # Load
        cmpMesh = openmesh.read_trimesh(filename, vertex_color=True, color_float=True)
            
        self.assertTrue(cmpMesh.has_vertex_colors())
            
        # Compare
        self.assertEqual(self.mesh.n_vertices(), 3)
        self.assertEqual(self.mesh.n_edges(), 3)
        self.assertEqual(self.mesh.n_faces(), 1)
        
        self.assertTrue(np.allclose(cmpMesh.point(v1), np.array([1.0, 0.0, 0.0])))
        self.assertTrue(np.allclose(cmpMesh.point(v2), np.array([0.0, 1.0, 0.0])))
        self.assertTrue(np.allclose(cmpMesh.point(v3), np.array([0.0, 0.0, 1.0])))
        
        self.assertAlmostEqual(cmpMesh.color(v1)[0], c1[0], 2)
        self.assertAlmostEqual(cmpMesh.color(v1)[1], c1[1], 2)
        self.assertAlmostEqual(cmpMesh.color(v1)[2], c1[2], 2)
        self.assertAlmostEqual(cmpMesh.color(v1)[3], c1[3], 2)
        
        self.assertAlmostEqual(cmpMesh.color(v2)[0], c2[0], 2)
        self.assertAlmostEqual(cmpMesh.color(v2)[1], c2[1], 2)
        self.assertAlmostEqual(cmpMesh.color(v2)[2], c2[2], 2)
        self.assertAlmostEqual(cmpMesh.color(v2)[3], c2[3], 2)
        
        self.assertAlmostEqual(cmpMesh.color(v3)[0], c3[0], 2)
        self.assertAlmostEqual(cmpMesh.color(v3)[1], c3[1], 2)
        self.assertAlmostEqual(cmpMesh.color(v3)[2], c3[2], 2)
        self.assertAlmostEqual(cmpMesh.color(v3)[3], c3[3], 2)
        
        # Clean up
        cmpMesh.release_vertex_colors()
        os.remove(filename)

    # TODO property tests

    def test_read_nonexistent_om(self):
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_trimesh("TestFiles/nonexistent.om")
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_polymesh("TestFiles/nonexistent.om")


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteOM)
    unittest.TextTestRunner(verbosity=2).run(suite)

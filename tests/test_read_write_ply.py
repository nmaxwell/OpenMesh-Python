import unittest
import openmesh
import os

import numpy as np

class ReadWritePLY(unittest.TestCase):

    def setUp(self):
        if not os.path.exists('OutFiles'):
            os.makedirs('OutFiles')

    def test_read_write_read(self):
        mesh1 = openmesh.read_trimesh("TestFiles/cube-minimal.obj")
        openmesh.write_mesh('OutFiles/test_read_write_read.ply', mesh1)
        mesh2 = openmesh.read_trimesh('OutFiles/test_read_write_read.ply')
        self.assertTrue(np.allclose(mesh1.points(), mesh2.points()))
        self.assertTrue(np.array_equal(mesh1.face_vertex_indices(), mesh2.face_vertex_indices()))

    def test_load_simple_point_ply_file_with_bad_encoding(self):
        self.mesh = openmesh.read_trimesh("TestFiles/pointCloudBadEncoding.ply")
        self.assertEqual(self.mesh.n_vertices(), 10)
        self.assertEqual(self.mesh.n_edges(), 0)
        self.assertEqual(self.mesh.n_faces(), 0)

    def test_load_simple_point_ply_file_with_good_encoding(self):
        self.mesh = openmesh.read_trimesh("TestFiles/pointCloudGoodEncoding.ply")
        self.assertEqual(self.mesh.n_vertices(), 10)
        self.assertEqual(self.mesh.n_edges(), 0)
        self.assertEqual(self.mesh.n_faces(), 0)

    def test_load_simple_ply(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal.ply")
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)

    def test_load_simple_ply_force_vertex_colors_although_not_available(self):
        with self.assertRaises(RuntimeError):
            openmesh.read_trimesh("TestFiles/cube-minimal.ply", vertex_color=True)

    def test_load_simple_ply_with_vertex_colors(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal-vertexColors.ply", vertex_color=True)
        
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

    def test_load_ply_from_mesh_lab_with_vertex_colors(self):
        self.mesh = openmesh.read_trimesh("TestFiles/meshlab.ply", vertex_color=True)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
    
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
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

    def test_write_and_read_binary_ply_with_vertex_colors(self):
        self.mesh = openmesh.read_trimesh("TestFiles/meshlab.ply", vertex_color=True)
        openmesh.write_mesh("OutFiles/meshlab_binary.ply", self.mesh, vertex_color=True, binary=True)
        self.mesh = openmesh.read_trimesh("OutFiles/meshlab_binary.ply", vertex_color=True, binary=True)

        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)

        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 1.0)

        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)

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

    def test_write_and_read_ply_with_float_vertex_colors(self):
        self.mesh = openmesh.read_trimesh("TestFiles/meshlab.ply", vertex_color=True)
        openmesh.write_mesh("OutFiles/meshlab_float.ply", self.mesh, vertex_color=True, color_float=True)
        self.mesh = openmesh.read_trimesh("OutFiles/meshlab_float.ply", vertex_color=True, color_float=True)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
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

    def test_write_and_read_binary_ply_with_float_vertex_colors(self):
        self.mesh = openmesh.read_trimesh("TestFiles/meshlab.ply", vertex_color=True)
        openmesh.write_mesh("OutFiles/meshlab_binary_float.ply", self.mesh, vertex_color=True, color_float=True, binary=True)
        self.mesh = openmesh.read_trimesh("OutFiles/meshlab_binary_float.ply", vertex_color=True, color_float=True, binary=True)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
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

    def test_load_simple_ply_with_texcoords(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal-texCoords.ply", vertex_tex_coord=True)
        
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

    def test_load_simple_ply_with_normals(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal-normals.ply", vertex_normal=True)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertTrue(self.mesh.has_vertex_normals())
        self.assertFalse(self.mesh.has_vertex_texcoords1D())
        self.assertFalse(self.mesh.has_vertex_texcoords2D())
        self.assertFalse(self.mesh.has_vertex_texcoords3D())
        self.assertFalse(self.mesh.has_vertex_colors())
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[2], 1.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[0], 1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[1], 0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[2], 0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[0], 1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[2], 1.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[0], 1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[1], 1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[2], 2.0)
        
        self.mesh.release_vertex_normals()

    def test_read_nonexistent_ply(self):
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_trimesh("TestFiles/nonexistent.ply")
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_polymesh("TestFiles/nonexistent.ply")


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWritePLY)
    unittest.TextTestRunner(verbosity=2).run(suite)

import unittest
import openmesh
import os

import numpy as np

class ReadWriteOFF(unittest.TestCase):

    def setUp(self):
        if not os.path.exists('OutFiles'):
            os.makedirs('OutFiles')

    def test_read_write_read(self):
        mesh1 = openmesh.read_trimesh("TestFiles/cube-minimal.obj")
        openmesh.write_mesh('OutFiles/test_read_write_read.off', mesh1)
        mesh2 = openmesh.read_trimesh('OutFiles/test_read_write_read.off')
        self.assertTrue(np.allclose(mesh1.points(), mesh2.points()))
        self.assertTrue(np.array_equal(mesh1.face_vertex_indices(), mesh2.face_vertex_indices()))

    def test_load_simple_off_file(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube1.off")
        self.assertEqual(self.mesh.n_vertices(), 7526)
        self.assertEqual(self.mesh.n_edges(), 22572)
        self.assertEqual(self.mesh.n_faces(), 15048)

    def test_write_and_read_vertex_colors_to_and_from_off_file(self):
        self.mesh = openmesh.TriMesh()
        self.mesh.request_vertex_colors()

        self.mesh.add_vertex(np.array([0, 0, 1]))
        self.mesh.add_vertex(np.array([0, 1, 0]))
        self.mesh.add_vertex(np.array([0, 1, 1]))
        self.mesh.add_vertex(np.array([1, 0, 1]))
        
        # Using the default openmesh Python color type
        testColor = np.array([1.0, 0.5, 0.25, 1.0])
        
        # Setting colors (different from black)
        for v in self.mesh.vertices():
            self.mesh.set_color(v, testColor)
        
        # Check if the colors are correctly set
        count = 0
        for v in self.mesh.vertices():
            color = self.mesh.color(v)
            if color[0] != testColor[0] or color[1] != testColor[1] or color[2] != testColor[2]:
                count += 1

        self.assertEqual(count, 0)

        openmesh.write_mesh("OutFiles/temp.off", self.mesh, vertex_color=True, color_float=True)
        self.mesh = openmesh.read_trimesh("OutFiles/temp.off", vertex_color=True, color_float=True)

        # Check if vertices still have the same color
        count = 0
        for v in self.mesh.vertices():
            color = self.mesh.color(v)
            if color[0] != testColor[0] or color[1] != testColor[1] or color[2] != testColor[2]:
                count += 1

        self.assertEqual(count, 0)

        self.mesh.release_vertex_colors()

    def test_write_and_read_float_vertex_colors_to_and_from_off_file(self):
        # Read the mesh
        self.mesh = openmesh.read_trimesh("TestFiles/meshlab.ply", vertex_color=True)

        # Write the mesh
        openmesh.write_mesh("OutFiles/cube_floating.off", self.mesh, vertex_color=True, color_float=True)
        
        # Read the mesh again
        self.mesh = openmesh.read_trimesh("OutFiles/cube_floating.off", vertex_color=True, color_float=True)

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

    def test_write_and_read_binary_float_vertex_colors_to_and_from_off_file(self):
        # Read the mesh
        self.mesh = openmesh.read_trimesh("TestFiles/meshlab.ply", vertex_color=True)

        # Write the mesh
        openmesh.write_mesh("OutFiles/cube_floating_binary.off", self.mesh, vertex_color=True, binary=True, color_float=True)

        # Read the mesh again
        self.mesh = openmesh.read_trimesh("OutFiles/cube_floating_binary.off", vertex_color=True, binary=True, color_float=True)
        
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
        self.assertFalse(self.mesh.has_face_colors())
        self.assertTrue(self.mesh.has_vertex_colors())
        
        self.mesh.release_vertex_colors()

    def test_read_nonexistent_off(self):
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_trimesh("TestFiles/nonexistent.off")
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_polymesh("TestFiles/nonexistent.off")


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteOFF)
    unittest.TextTestRunner(verbosity=2).run(suite)

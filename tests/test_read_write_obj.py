import unittest
import openmesh
import os

import numpy as np

class ReadWriteOBJ(unittest.TestCase):

    def setUp(self):
        if not os.path.exists('OutFiles'):
            os.makedirs('OutFiles')

    def test_read_write_read(self):
        mesh1 = openmesh.read_trimesh("TestFiles/cube-minimal.obj")
        openmesh.write_mesh('OutFiles/test_read_write_read.obj', mesh1)
        mesh2 = openmesh.read_trimesh('OutFiles/test_read_write_read.obj')
        self.assertTrue(np.allclose(mesh1.points(), mesh2.points()))
        self.assertTrue(np.array_equal(mesh1.face_vertex_indices(), mesh2.face_vertex_indices()))

    def test_load_simple_obj(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal.obj")
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)

    def test_load_simple_obj_check_halfedge_and_vertex_normals(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal.obj", vertex_normal=True)
        self.mesh.update_halfedge_normals()
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        self.assertEqual(self.mesh.n_halfedges(), 36)
        
        # =====================================================
        # Check vertex normals
        # =====================================================
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[1], -1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(0))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(3))[2],  1.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[1], -1.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(4))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.vertex_handle(7))[2],  1.0)
        
        # =====================================================
        # Check halfedge normals
        # =====================================================
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle( 0))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle( 0))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle( 0))[2], -1.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(10))[0], -1.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(10))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(10))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(19))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(19))[1],  1.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(19))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(24))[0],  1.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(24))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(24))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(30))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(30))[1], -1.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(30))[2],  0.0)
        
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(35))[0],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(35))[1],  0.0)
        self.assertEqual(self.mesh.normal(self.mesh.halfedge_handle(35))[2],  1.0)
        
        self.mesh.release_vertex_normals()
        self.mesh.release_halfedge_normals()

    def test_load_simple_obj_force_vertex_colors_although_not_available(self):
        with self.assertRaises(RuntimeError):
            openmesh.read_trimesh("TestFiles/cube-minimal.obj", vertex_color=True)

    def test_load_obj_with_material(self):
        self.mesh = openmesh.read_trimesh("TestFiles/square_material.obj", face_color=True)
        fh = self.mesh.face_handle(self.mesh.halfedge_handle(0))

        self.assertTrue(fh.is_valid())

        self.assertAlmostEqual(self.mesh.color(fh)[0], 0.5, 2)
        self.assertAlmostEqual(self.mesh.color(fh)[1], 0.5, 2)
        self.assertAlmostEqual(self.mesh.color(fh)[2], 0.5, 2)

        self.mesh.release_face_colors()

    def test_load_simple_obj_with_vertex_colors_after_vertices(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal-vertex-colors-after-vertex-definition.obj", vertex_color=True)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[1], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[2], 1.0)
        
        self.mesh.release_vertex_colors()

    def test_load_simple_obj_with_vertex_colors_as_vc_lines(self):
        self.mesh = openmesh.read_trimesh("TestFiles/cube-minimal-vertex-colors-as-vc-lines.obj", vertex_color=True)
        
        self.assertEqual(self.mesh.n_vertices(), 8)
        self.assertEqual(self.mesh.n_edges(), 18)
        self.assertEqual(self.mesh.n_faces(), 12)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(0))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[0], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[1], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(3))[2], 1.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[1], 0.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(4))[2], 0.0)
        
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[0], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[1], 1.0)
        self.assertEqual(self.mesh.color(self.mesh.vertex_handle(7))[2], 1.0)
        
        self.mesh.release_vertex_colors()

    def test_read_nonexistent_obj(self):
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_trimesh("TestFiles/nonexistent.obj")
        with self.assertRaises(RuntimeError):
            self.mesh = openmesh.read_polymesh("TestFiles/nonexistent.obj")


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ReadWriteOBJ)
    unittest.TextTestRunner(verbosity=2).run(suite)

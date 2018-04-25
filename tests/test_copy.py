import unittest
import openmesh
import copy

class Python(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.TriMesh()
        self.vh1 = self.mesh.add_vertex([0, 0, 0])
        self.vh2 = self.mesh.add_vertex([1, 0, 0])
        self.vh3 = self.mesh.add_vertex([1, 1, 0])
        self.fh = self.mesh.add_face(self.vh1, self.vh2, self.vh3)
        self.mesh.set_vertex_property('test', self.vh1, self.mesh)
        self.mesh.set_vertex_property('test', self.vh2, self.mesh)
        self.one_two_three = [1, 2, 3]
        self.mesh.set_vertex_property('test', self.vh3, self.one_two_three)
        self.mesh.set_face_property('test', self.fh, self.one_two_three)

    def test_shallowcopy(self):
        mesh_copy = copy.copy(self.mesh)
        mesh_copy.vertex_property('test', self.vh3)[:] = [4, 5, 6]
        self.assertIsNot(self.mesh, mesh_copy)
        self.assertIs(mesh_copy.vertex_property('test', self.vh1), self.mesh)
        self.assertIs(mesh_copy.vertex_property('test', self.vh2), self.mesh)
        self.assertIs(mesh_copy.vertex_property('test', self.vh3), self.one_two_three)
        self.assertIs(mesh_copy.face_property('test', self.fh), self.one_two_three)
        self.assertIs(mesh_copy.vertex_property('test', self.vh3), mesh_copy.face_property('test', self.fh))
        self.assertEqual(self.mesh.vertex_property('test', self.vh3), [4, 5, 6])
        self.assertEqual(mesh_copy.vertex_property('test', self.vh3), [4, 5, 6])

    def test_deepcopy(self):
        mesh_copy = copy.deepcopy(self.mesh)
        mesh_copy.vertex_property('test', self.vh3)[:] = [4, 5, 6]
        self.assertIsNot(self.mesh, mesh_copy)
        self.assertIs(mesh_copy.vertex_property('test', self.vh1), mesh_copy)
        self.assertIs(mesh_copy.vertex_property('test', self.vh2), mesh_copy)
        self.assertIsNot(mesh_copy.vertex_property('test', self.vh3), self.one_two_three)
        self.assertIsNot(mesh_copy.face_property('test', self.fh), self.one_two_three)
        self.assertIs(mesh_copy.vertex_property('test', self.vh3), mesh_copy.face_property('test', self.fh))
        self.assertListEqual(self.mesh.vertex_property('test', self.vh3), [1, 2, 3])
        self.assertListEqual(mesh_copy.vertex_property('test', self.vh3), [4, 5, 6])


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

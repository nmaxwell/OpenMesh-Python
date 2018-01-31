import unittest
import openmesh

import numpy as np

class Python(unittest.TestCase):

    def setUp(self):
        pass

    def one_triangle(self):
        mesh = openmesh.TriMesh()
        vh0 = mesh.add_vertex([0, 0, 0])
        vh1 = mesh.add_vertex([1, 0, 0])
        vh2 = mesh.add_vertex([1, 1, 0])
        mesh.add_face(vh0, vh1, vh2)
        return mesh

    def test_on_demand_getter(self):
        # check if standard properties are requested on demand,
        # i.e. if it doesn't crash it works
        mesh = self.one_triangle()
        mesh.normal(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.normal(openmesh.HalfedgeHandle(0))

        mesh = self.one_triangle()
        mesh.normal(openmesh.FaceHandle(0))

        mesh = self.one_triangle()
        mesh.color(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.color(openmesh.HalfedgeHandle(0))

        mesh = self.one_triangle()
        mesh.color(openmesh.EdgeHandle(0))

        mesh = self.one_triangle()
        mesh.color(openmesh.FaceHandle(0))

        mesh = self.one_triangle()
        mesh.color(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.texcoord1D(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.texcoord1D(openmesh.HalfedgeHandle(0))

        mesh = self.one_triangle()
        mesh.texcoord2D(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.texcoord2D(openmesh.HalfedgeHandle(0))

        mesh = self.one_triangle()
        mesh.texcoord3D(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.texcoord3D(openmesh.HalfedgeHandle(0))

    def test_on_demand_setter(self):
        # check if standard properties are requested on demand,
        # i.e. if it doesn't crash it works
        mesh = self.one_triangle()
        mesh.set_normal(openmesh.VertexHandle(0), np.array([0, 0, 0]))

        mesh = self.one_triangle()
        mesh.set_normal(openmesh.HalfedgeHandle(0), np.array([0, 0, 0]))

        mesh = self.one_triangle()
        mesh.set_normal(openmesh.FaceHandle(0), np.array([0, 0, 0]))

        mesh = self.one_triangle()
        mesh.set_color(openmesh.VertexHandle(0), np.array([0, 0, 0, 0]))

        mesh = self.one_triangle()
        mesh.set_color(openmesh.HalfedgeHandle(0), np.array([0, 0, 0, 0]))

        mesh = self.one_triangle()
        mesh.set_color(openmesh.EdgeHandle(0), np.array([0, 0, 0, 0]))

        mesh = self.one_triangle()
        mesh.set_color(openmesh.FaceHandle(0), np.array([0, 0, 0, 0]))

        mesh = self.one_triangle()
        mesh.set_texcoord1D(openmesh.VertexHandle(0), np.array([0]))

        mesh = self.one_triangle()
        mesh.set_texcoord1D(openmesh.HalfedgeHandle(0), np.array([0]))

        mesh = self.one_triangle()
        mesh.set_texcoord2D(openmesh.VertexHandle(0), np.array([0, 0]))

        mesh = self.one_triangle()
        mesh.set_texcoord2D(openmesh.HalfedgeHandle(0), np.array([0, 0]))

        mesh = self.one_triangle()
        mesh.set_texcoord3D(openmesh.VertexHandle(0), np.array([0, 0, 0]))

        mesh = self.one_triangle()
        mesh.set_texcoord3D(openmesh.HalfedgeHandle(0), np.array([0, 0, 0]))

    def test_on_demand_arrays(self):
        # check if standard properties are requested on demand,
        # i.e. if it doesn't crash it works
        mesh = self.one_triangle()
        mesh.vertex_normals()

        mesh = self.one_triangle()
        mesh.vertex_colors()

        mesh = self.one_triangle()
        mesh.vertex_texcoords1D()

        mesh = self.one_triangle()
        mesh.vertex_texcoords2D()

        mesh = self.one_triangle()
        mesh.vertex_texcoords3D()

        mesh = self.one_triangle()
        mesh.halfedge_normals()

        mesh = self.one_triangle()
        mesh.halfedge_colors()

        mesh = self.one_triangle()
        mesh.halfedge_texcoords1D()

        mesh = self.one_triangle()
        mesh.halfedge_texcoords2D()

        mesh = self.one_triangle()
        mesh.halfedge_texcoords3D()

        mesh = self.one_triangle()
        mesh.edge_colors()

        mesh = self.one_triangle()
        mesh.face_normals()

        mesh = self.one_triangle()
        mesh.face_colors()

    def test_on_demand_delete(self):
        # check if standard properties are requested on demand,
        # i.e. if it doesn't crash it works
        mesh = self.one_triangle()
        mesh.delete_vertex(openmesh.VertexHandle(0))
        mesh.garbage_collection()

        mesh = self.one_triangle()
        mesh.delete_edge(openmesh.EdgeHandle(0))
        mesh.garbage_collection()

        mesh = self.one_triangle()
        mesh.delete_face(openmesh.FaceHandle(0))
        mesh.garbage_collection()

    def test_on_demand_update(self):
        # check if standard properties are requested on demand,
        # i.e. if it doesn't crash it works
        mesh = self.one_triangle()
        mesh.update_normal(openmesh.FaceHandle(0))

        mesh = self.one_triangle()
        mesh.update_face_normals()

        mesh = self.one_triangle()
        mesh.update_normal(openmesh.HalfedgeHandle(0))

        mesh = self.one_triangle()
        mesh.update_halfedge_normals()

        mesh = self.one_triangle()
        mesh.update_normal(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.update_vertex_normals()

        mesh = self.one_triangle()
        mesh.update_normals()

    def test_on_demand_calc(self):
        # check if standard properties are requested on demand,
        # i.e. if it doesn't crash it works
        mesh = self.one_triangle()
        mesh.calc_halfedge_normal(openmesh.HalfedgeHandle(0))

        mesh = self.one_triangle()
        mesh.calc_vertex_normal(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.calc_vertex_normal_fast(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.calc_vertex_normal_correct(openmesh.VertexHandle(0))

        mesh = self.one_triangle()
        mesh.calc_vertex_normal_loop(openmesh.VertexHandle(0))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

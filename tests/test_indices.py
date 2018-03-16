import unittest
import openmesh

import numpy as np

class Python(unittest.TestCase):

    def setUp(self):
        self.mesh = openmesh.PolyMesh()
        vh0 = self.mesh.add_vertex(np.array([0, 0, 0]))
        vh1 = self.mesh.add_vertex(np.array([1, 0, 0]))
        vh2 = self.mesh.add_vertex(np.array([1, 1, 0]))
        vh3 = self.mesh.add_vertex(np.array([0, 1, 0]))
        vh4 = self.mesh.add_vertex(np.array([2, 0, 0]))
        vh5 = self.mesh.add_vertex(np.array([3, 0, 0]))
        vh6 = self.mesh.add_vertex(np.array([3, 1, 0]))
        vh7 = self.mesh.add_vertex(np.array([2, 1, 0]))
        vh8 = self.mesh.add_vertex(np.array([4, 0, 0]))
        vh9 = self.mesh.add_vertex(np.array([5, 0, 0]))
        vh10 = self.mesh.add_vertex(np.array([5, 1, 0]))
        self.mesh.add_face(vh0, vh1, vh2, vh3)
        self.mesh.add_face(vh4, vh5, vh6)
        self.mesh.add_face(vh6, vh7, vh4)
        self.mesh.add_face(vh8, vh9, vh10)

        # Test setup:
        #  3 === 2   7 === 6        10
        #  |     |   |   / |       / |
        #  |     |   |  /  |      /  |
        #  |     |   | /   |     /   |
        #  0 === 1   4 === 5   8 === 9

    def delete_vertices(self):
        for vh in self.mesh.vertices():
            self.mesh.delete_vertex(vh)

    def test_vertex_vertex_indices(self):
        indices1 = self.mesh.vertex_vertex_indices()
        indices2 = self.mesh.vv_indices()
        for h1 in self.mesh.vertices():
            correct = [h2.idx() for h2 in self.mesh.vv(h1)]
            while len(correct) < 3:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_vertices(), 3))
        self.assertEqual(indices2.shape, (self.mesh.n_vertices(), 3))
        
    def test_vertex_vertex_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.vertex_vertex_indices)
        self.assertRaises(RuntimeError, self.mesh.vv_indices)
        
    def test_vertex_face_indices(self):
        indices1 = self.mesh.vertex_face_indices()
        indices2 = self.mesh.vf_indices()
        for h1 in self.mesh.vertices():
            correct = [h2.idx() for h2 in self.mesh.vf(h1)]
            while len(correct) < 2:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_vertices(), 2))
        self.assertEqual(indices2.shape, (self.mesh.n_vertices(), 2))
        
    def test_vertex_face_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.vertex_face_indices)
        self.assertRaises(RuntimeError, self.mesh.vf_indices)
        
    def test_vertex_edge_indices(self):
        indices1 = self.mesh.vertex_edge_indices()
        indices2 = self.mesh.ve_indices()
        for h1 in self.mesh.vertices():
            correct = [h2.idx() for h2 in self.mesh.ve(h1)]
            while len(correct) < 3:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_vertices(), 3))
        self.assertEqual(indices2.shape, (self.mesh.n_vertices(), 3))
        
    def test_vertex_edge_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.vertex_edge_indices)
        self.assertRaises(RuntimeError, self.mesh.ve_indices)
        
    def test_vertex_outgoing_halfedge_indices(self):
        indices1 = self.mesh.vertex_outgoing_halfedge_indices()
        indices2 = self.mesh.voh_indices()
        for h1 in self.mesh.vertices():
            correct = [h2.idx() for h2 in self.mesh.voh(h1)]
            while len(correct) < 3:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_vertices(), 3))
        self.assertEqual(indices2.shape, (self.mesh.n_vertices(), 3))
        
    def test_vertex_outgoing_halfedge_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.vertex_outgoing_halfedge_indices)
        self.assertRaises(RuntimeError, self.mesh.voh_indices)
        
    def test_vertex_incoming_halfedge_indices(self):
        indices1 = self.mesh.vertex_incoming_halfedge_indices()
        indices2 = self.mesh.vih_indices()
        for h1 in self.mesh.vertices():
            correct = [h2.idx() for h2 in self.mesh.vih(h1)]
            while len(correct) < 3:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_vertices(), 3))
        self.assertEqual(indices2.shape, (self.mesh.n_vertices(), 3))
        
    def test_vertex_incoming_halfedge_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.vertex_incoming_halfedge_indices)
        self.assertRaises(RuntimeError, self.mesh.vih_indices)
        
    def test_face_vertex_indices_trimesh(self):
        self.mesh = openmesh.read_trimesh('TestFiles/cube-minimal.obj')
        indices1 = self.mesh.face_vertex_indices()
        indices2 = self.mesh.fv_indices()
        for fh in self.mesh.faces():
            fv_it = self.mesh.fv(fh)
            vh1 = next(fv_it)
            vh2 = next(fv_it)
            vh3 = next(fv_it)
            self.assertEqual(indices1[fh.idx(), 0], vh1.idx())
            self.assertEqual(indices1[fh.idx(), 1], vh2.idx())
            self.assertEqual(indices1[fh.idx(), 2], vh3.idx())
            self.assertEqual(indices2[fh.idx(), 0], vh1.idx())
            self.assertEqual(indices2[fh.idx(), 1], vh2.idx())
            self.assertEqual(indices2[fh.idx(), 2], vh3.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_faces(), 3))
        self.assertEqual(indices2.shape, (self.mesh.n_faces(), 3))
        
    def test_face_vertex_indices_trimesh_delete(self):
        self.mesh = openmesh.read_trimesh('TestFiles/cube-minimal.obj')
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.face_vertex_indices)
        self.assertRaises(RuntimeError, self.mesh.fv_indices)
        
    def test_face_vertex_indices_polymesh(self):
        indices1 = self.mesh.face_vertex_indices()
        indices2 = self.mesh.fv_indices()
        for h1 in self.mesh.faces():
            correct = [h2.idx() for h2 in self.mesh.fv(h1)]
            while len(correct) < 4:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_faces(), 4))
        self.assertEqual(indices2.shape, (self.mesh.n_faces(), 4))
        
    def test_face_vertex_indices_polymesh_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.face_vertex_indices)
        self.assertRaises(RuntimeError, self.mesh.fv_indices)
        
    def test_face_face_indices(self):
        indices1 = self.mesh.face_face_indices()
        indices2 = self.mesh.ff_indices()
        for h1 in self.mesh.faces():
            correct = [h2.idx() for h2 in self.mesh.ff(h1)]
            while len(correct) < 1:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_faces(), 1))
        self.assertEqual(indices2.shape, (self.mesh.n_faces(), 1))
        
    def test_face_face_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.face_face_indices)
        self.assertRaises(RuntimeError, self.mesh.ff_indices)
        
    def test_face_edge_indices(self):
        indices1 = self.mesh.face_edge_indices()
        indices2 = self.mesh.fe_indices()
        for h1 in self.mesh.faces():
            correct = [h2.idx() for h2 in self.mesh.fe(h1)]
            while len(correct) < 4:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_faces(), 4))
        self.assertEqual(indices2.shape, (self.mesh.n_faces(), 4))
        
    def test_face_edge_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.face_edge_indices)
        self.assertRaises(RuntimeError, self.mesh.fe_indices)
        
    def test_face_halfedge_indices(self):
        indices1 = self.mesh.face_halfedge_indices()
        indices2 = self.mesh.fh_indices()
        for h1 in self.mesh.faces():
            correct = [h2.idx() for h2 in self.mesh.fh(h1)]
            while len(correct) < 4:
                correct.append(-1)
            correct = np.array(correct)
            self.assertTrue(np.array_equal(indices1[h1.idx()], correct))
            self.assertTrue(np.array_equal(indices2[h1.idx()], correct))
        self.assertEqual(indices1.shape, (self.mesh.n_faces(), 4))
        self.assertEqual(indices2.shape, (self.mesh.n_faces(), 4))
        
    def test_face_halfedge_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.face_halfedge_indices)
        self.assertRaises(RuntimeError, self.mesh.fh_indices)
        
    def test_edge_vertex_indices(self):
        indices1 = self.mesh.edge_vertex_indices()
        indices2 = self.mesh.ev_indices()
        for eh in self.mesh.edges():
            heh = self.mesh.halfedge_handle(eh, 0)
            vh1 = self.mesh.from_vertex_handle(heh)
            vh2 = self.mesh.to_vertex_handle(heh)
            self.assertEqual(indices1[eh.idx(), 0], vh1.idx())
            self.assertEqual(indices1[eh.idx(), 1], vh2.idx())
            self.assertEqual(indices2[eh.idx(), 0], vh1.idx())
            self.assertEqual(indices2[eh.idx(), 1], vh2.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_edges(), 2))
        self.assertEqual(indices2.shape, (self.mesh.n_edges(), 2))
        
    def test_edge_vertex_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.edge_vertex_indices)
        self.assertRaises(RuntimeError, self.mesh.ev_indices)
        
    def test_edge_face_indices(self):
        indices1 = self.mesh.edge_face_indices()
        indices2 = self.mesh.ef_indices()
        for eh in self.mesh.edges():
            heh1 = self.mesh.halfedge_handle(eh, 0)
            heh2 = self.mesh.halfedge_handle(eh, 1)
            fh1 = self.mesh.face_handle(heh1)
            fh2 = self.mesh.face_handle(heh2)
            self.assertEqual(indices1[eh.idx(), 0], fh1.idx())
            self.assertEqual(indices1[eh.idx(), 1], fh2.idx())
            self.assertEqual(indices2[eh.idx(), 0], fh1.idx())
            self.assertEqual(indices2[eh.idx(), 1], fh2.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_edges(), 2))
        self.assertEqual(indices2.shape, (self.mesh.n_edges(), 2))
        
    def test_edge_face_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.edge_face_indices)
        self.assertRaises(RuntimeError, self.mesh.ef_indices)
        
    def test_edge_halfedge_indices(self):
        indices1 = self.mesh.edge_halfedge_indices()
        indices2 = self.mesh.eh_indices()
        for eh in self.mesh.edges():
            heh1 = self.mesh.halfedge_handle(eh, 0)
            heh2 = self.mesh.halfedge_handle(eh, 1)
            self.assertEqual(indices1[eh.idx(), 0], heh1.idx())
            self.assertEqual(indices1[eh.idx(), 1], heh2.idx())
            self.assertEqual(indices2[eh.idx(), 0], heh1.idx())
            self.assertEqual(indices2[eh.idx(), 1], heh2.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_edges(), 2))
        self.assertEqual(indices2.shape, (self.mesh.n_edges(), 2))
        
    def test_edge_halfedge_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.edge_halfedge_indices)
        self.assertRaises(RuntimeError, self.mesh.eh_indices)

    def test_halfedge_vertex_indices(self):
        indices1 = self.mesh.halfedge_vertex_indices()
        indices2 = self.mesh.hv_indices()
        for heh in self.mesh.halfedges():
            vh1 = self.mesh.from_vertex_handle(heh)
            vh2 = self.mesh.to_vertex_handle(heh)
            self.assertEqual(indices1[heh.idx(), 0], vh1.idx())
            self.assertEqual(indices1[heh.idx(), 1], vh2.idx())
            self.assertEqual(indices2[heh.idx(), 0], vh1.idx())
            self.assertEqual(indices2[heh.idx(), 1], vh2.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_halfedges(), 2))
        self.assertEqual(indices2.shape, (self.mesh.n_halfedges(), 2))
        
    def test_halfedge_vertex_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.halfedge_vertex_indices)
        self.assertRaises(RuntimeError, self.mesh.hv_indices)
        
    def test_halfedge_to_vertex_indices(self):
        indices1 = self.mesh.halfedge_to_vertex_indices()
        indices2 = self.mesh.htv_indices()
        for heh in self.mesh.halfedges():
            vh = self.mesh.to_vertex_handle(heh)
            self.assertEqual(indices1[heh.idx()], vh.idx())
            self.assertEqual(indices2[heh.idx()], vh.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_halfedges(),))
        self.assertEqual(indices2.shape, (self.mesh.n_halfedges(),))
        
    def test_halfedge_to_vertex_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.halfedge_to_vertex_indices)
        self.assertRaises(RuntimeError, self.mesh.htv_indices)
        
    def test_halfedge_from_vertex_indices(self):
        indices1 = self.mesh.halfedge_from_vertex_indices()
        indices2 = self.mesh.hfv_indices()
        for heh in self.mesh.halfedges():
            vh = self.mesh.from_vertex_handle(heh)
            self.assertEqual(indices1[heh.idx()], vh.idx())
            self.assertEqual(indices2[heh.idx()], vh.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_halfedges(),))
        self.assertEqual(indices2.shape, (self.mesh.n_halfedges(),))
        
    def test_halfedge_from_vertex_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.halfedge_from_vertex_indices)
        self.assertRaises(RuntimeError, self.mesh.hfv_indices)
        
    def test_halfedge_face_indices(self):
        indices1 = self.mesh.halfedge_face_indices()
        indices2 = self.mesh.hf_indices()
        for heh in self.mesh.halfedges():
            fh = self.mesh.face_handle(heh)
            self.assertEqual(indices1[heh.idx()], fh.idx())
            self.assertEqual(indices2[heh.idx()], fh.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_halfedges(),))
        self.assertEqual(indices2.shape, (self.mesh.n_halfedges(),))
        
    def test_halfedge_face_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.halfedge_face_indices)
        self.assertRaises(RuntimeError, self.mesh.hf_indices)
        
    def test_halfedge_edge_indices(self):
        indices1 = self.mesh.halfedge_edge_indices()
        indices2 = self.mesh.he_indices()
        for heh in self.mesh.halfedges():
            eh = self.mesh.edge_handle(heh)
            self.assertEqual(indices1[heh.idx()], eh.idx())
            self.assertEqual(indices2[heh.idx()], eh.idx())
        self.assertEqual(indices1.shape, (self.mesh.n_halfedges(),))
        self.assertEqual(indices2.shape, (self.mesh.n_halfedges(),))
        
    def test_halfedge_edge_indices_delete(self):
        self.delete_vertices()
        self.assertRaises(RuntimeError, self.mesh.halfedge_edge_indices)
        self.assertRaises(RuntimeError, self.mesh.he_indices)
        
    def test_init_with_arrays(self):
        points = self.mesh.points()
        face_vertex_indices = self.mesh.face_vertex_indices()
        # init polymesh
        polymesh = openmesh.PolyMesh(points, face_vertex_indices)
        self.assertEqual(polymesh.n_vertices(), self.mesh.n_vertices())
        self.assertEqual(polymesh.n_faces(), self.mesh.n_faces())
        # init trimesh (one face will be triangulated)
        trimesh = openmesh.TriMesh(points, face_vertex_indices)
        self.assertEqual(trimesh.n_vertices(), self.mesh.n_vertices())
        self.assertEqual(trimesh.n_faces(), self.mesh.n_faces() + 1)
        # init with empty points and faces
        trimesh = openmesh.TriMesh(np.empty((0, 3)), np.empty((0, 4)))
        self.assertEqual(trimesh.n_vertices(), 0)
        self.assertEqual(trimesh.n_faces(), 0)
        # init with empty points
        trimesh = openmesh.TriMesh(np.empty((0, 3)), face_vertex_indices)
        self.assertEqual(trimesh.n_vertices(), 0)
        self.assertEqual(trimesh.n_faces(), 0)
        # init with empty faces
        trimesh = openmesh.TriMesh(points, np.empty((0, 4)))
        self.assertEqual(trimesh.n_vertices(), self.mesh.n_vertices())
        self.assertEqual(trimesh.n_faces(), 0)
        # init with points only
        trimesh = openmesh.TriMesh(points)
        self.assertEqual(trimesh.n_vertices(), self.mesh.n_vertices())
        self.assertEqual(trimesh.n_faces(), 0)
        # init with wrong points shape
        with self.assertRaises(RuntimeError):
            openmesh.TriMesh(points[:, :2])
        # init with wrong faces shape
        with self.assertRaises(RuntimeError):
            openmesh.TriMesh(points, face_vertex_indices[:, :2])
        # init with points and invalid faces
        face_vertex_indices[1] = [-1, -1, -1, -1]
        face_vertex_indices[3] = [-1, -1, -1, -1]
        polymesh = openmesh.PolyMesh(points, face_vertex_indices)
        self.assertEqual(polymesh.n_vertices(), self.mesh.n_vertices())
        self.assertEqual(polymesh.n_faces(), self.mesh.n_faces() - 2)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Python)
    unittest.TextTestRunner(verbosity=2).run(suite)

import unittest
import openmesh

import numpy as np

class DeleteItems(unittest.TestCase):
    
    def setUp(self):
        self.mesh = openmesh.TriMesh()
        
        # Test setup:
        #
        # 5 ======== 4 ======== 3
        # |  0    /  |  2    /  |
        # |     /    |     /    |
        # |   /      |   /      |
        # | /     1  | /     3  |
        # 0 ======== 1 ======== 2
        
        vh0 = self.mesh.add_vertex([0, 0, 0])
        vh1 = self.mesh.add_vertex([1, 0, 0])
        vh2 = self.mesh.add_vertex([2, 0, 0])
        vh3 = self.mesh.add_vertex([2, 1, 0])
        vh4 = self.mesh.add_vertex([1, 1, 0])
        vh5 = self.mesh.add_vertex([0, 1, 0])
        
        fh0 = self.mesh.add_face(vh0, vh4, vh5)
        fh1 = self.mesh.add_face(vh0, vh1, vh4)
        fh2 = self.mesh.add_face(vh1, vh3, vh4)
        fh3 = self.mesh.add_face(vh1, vh2, vh3)

    def test_delete_vertex_with_idx_1(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete
        self.mesh.delete_vertex(self.mesh.vertex_handle(1))
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 3)
        self.assertEqual(self.mesh.n_halfedges(), 6)
        self.assertEqual(self.mesh.n_edges(), 3)
        self.assertEqual(self.mesh.n_faces(), 1)
        
    def test_delete_vertex_with_idx_1_without_isolated(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete
        self.mesh.delete_vertex(self.mesh.vertex_handle(1), False)
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 5)
        self.assertEqual(self.mesh.n_halfedges(), 6)
        self.assertEqual(self.mesh.n_edges(), 3)
        self.assertEqual(self.mesh.n_faces(), 1)
        
    def test_delete_vertex_with_idx_2(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete
        self.mesh.delete_vertex(self.mesh.vertex_handle(2))
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 5)
        self.assertEqual(self.mesh.n_halfedges(), 14)
        self.assertEqual(self.mesh.n_edges(), 7)
        self.assertEqual(self.mesh.n_faces(), 3)
        
    def test_delete_vertex_with_idx_2_without_isolated(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete
        self.mesh.delete_vertex(self.mesh.vertex_handle(2), False)
        self.mesh.garbage_collection()

        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 5)
        self.assertEqual(self.mesh.n_halfedges(), 14)
        self.assertEqual(self.mesh.n_edges(), 7)
        self.assertEqual(self.mesh.n_faces(), 3)
        
    def test_delete_edge_with_idx_1(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_edge(self.mesh.edge_handle(1))
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 5)
        self.assertEqual(self.mesh.n_halfedges(), 14)
        self.assertEqual(self.mesh.n_edges(), 7)
        self.assertEqual(self.mesh.n_faces(), 3)
        
    def test_delete_edge_with_idx_1_without_isolated(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_edge(self.mesh.edge_handle(1), False)
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 14)
        self.assertEqual(self.mesh.n_edges(), 7)
        self.assertEqual(self.mesh.n_faces(), 3)
        
    def test_delete_edge_with_idx_4(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_edge(self.mesh.edge_handle(4))
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 12)
        self.assertEqual(self.mesh.n_edges(), 6)
        self.assertEqual(self.mesh.n_faces(), 2)
        
    def test_delete_edge_with_idx_4_without_isolated(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_edge(self.mesh.edge_handle(4), False)
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 12)
        self.assertEqual(self.mesh.n_edges(), 6)
        self.assertEqual(self.mesh.n_faces(), 2)
        
    def test_delete_edge_with_idx_5(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_edge(self.mesh.edge_handle(5))
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 4)
        self.assertEqual(self.mesh.n_halfedges(), 10)
        self.assertEqual(self.mesh.n_edges(), 5)
        self.assertEqual(self.mesh.n_faces(), 2)
        
    def test_delete_edge_with_idx_5_without_isolated(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_edge(self.mesh.edge_handle(5), False)
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 10)
        self.assertEqual(self.mesh.n_edges(), 5)
        self.assertEqual(self.mesh.n_faces(), 2)
        
    def test_delete_face_with_idx_1(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_face(self.mesh.face_handle(1))
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 16)
        self.assertEqual(self.mesh.n_edges(), 8)
        self.assertEqual(self.mesh.n_faces(), 3)
        
    def test_delete_face_with_idx_1_without_isolated(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_face(self.mesh.face_handle(1), False)
        self.mesh.garbage_collection()

        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 16)
        self.assertEqual(self.mesh.n_edges(), 8)
        self.assertEqual(self.mesh.n_faces(), 3)
        
    def test_delete_face_with_idx_3(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_face(self.mesh.face_handle(3))
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 5)
        self.assertEqual(self.mesh.n_halfedges(), 14)
        self.assertEqual(self.mesh.n_edges(), 7)
        self.assertEqual(self.mesh.n_faces(), 3)
        
    def test_delete_face_with_idx_3_without_isolated(self):
        # check setup
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 18)
        self.assertEqual(self.mesh.n_edges(), 9)
        self.assertEqual(self.mesh.n_faces(), 4)
        
        # delete center edge
        self.mesh.delete_face(self.mesh.face_handle(3), False)
        self.mesh.garbage_collection()
        
        # check mesh after deletion
        self.assertEqual(self.mesh.n_vertices(), 6)
        self.assertEqual(self.mesh.n_halfedges(), 14)
        self.assertEqual(self.mesh.n_edges(), 7)
        self.assertEqual(self.mesh.n_faces(), 3)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(DeleteItems)
    unittest.TextTestRunner(verbosity=2).run(suite)

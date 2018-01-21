#include "Bindings.hh"
#include "Miscellaneous.hh"
#include "Vector.hh"

#include <pybind11/pybind11.h>
namespace py = pybind11;


PYBIND11_MODULE(openmesh, m) {
	expose_items(m);
	expose_handles(m);
	expose_status_bits_and_info(m);

	expose_vec<float,  2>(m, "Vec2f");
	expose_vec<float,  3>(m, "Vec3f");
	expose_vec<float,  4>(m, "Vec4f");
	expose_vec<double, 2>(m, "Vec2d");
	expose_vec<double, 3>(m, "Vec3d");
	expose_vec<double, 4>(m, "Vec4d");

//	expose_mesh<PolyMesh>("PolyMesh");
//	expose_mesh<TriMesh>("TriMesh");

//	expose_iterator<OpenMesh::PolyConnectivity::VertexIter, &OpenMesh::ArrayKernel::n_vertices>("VertexIter");
//	expose_iterator<OpenMesh::PolyConnectivity::HalfedgeIter, &OpenMesh::ArrayKernel::n_halfedges>("HalfedgeIter");
//	expose_iterator<OpenMesh::PolyConnectivity::EdgeIter, &OpenMesh::ArrayKernel::n_edges>("EdgeIter");
//	expose_iterator<OpenMesh::PolyConnectivity::FaceIter, &OpenMesh::ArrayKernel::n_faces>("FaceIter");

//	expose_circulator<OpenMesh::PolyConnectivity::VertexVertexIter, VertexHandle>("VertexVertexIter");
//	expose_circulator<OpenMesh::PolyConnectivity::VertexIHalfedgeIter, VertexHandle>("VertexIHalfedgeIter");
//	expose_circulator<OpenMesh::PolyConnectivity::VertexOHalfedgeIter, VertexHandle>("VertexOHalfedgeIter");
//	expose_circulator<OpenMesh::PolyConnectivity::VertexEdgeIter, VertexHandle>("VertexEdgeIter");
//	expose_circulator<OpenMesh::PolyConnectivity::VertexFaceIter, VertexHandle>("VertexFaceIter");

//	expose_circulator<OpenMesh::PolyConnectivity::FaceVertexIter, FaceHandle>("FaceVertexIter");
//	expose_circulator<OpenMesh::PolyConnectivity::FaceHalfedgeIter, FaceHandle>("FaceHalfedgeIter");
//	expose_circulator<OpenMesh::PolyConnectivity::FaceEdgeIter, FaceHandle>("FaceEdgeIter");
//	expose_circulator<OpenMesh::PolyConnectivity::FaceFaceIter, FaceHandle>("FaceFaceIter");

//	expose_circulator<OpenMesh::PolyConnectivity::HalfedgeLoopIter, HalfedgeHandle>("HalfedgeLoopIter");

//	typedef IteratorWrapperT<PolyConnectivity::VertexIter, &ArrayKernel::n_vertices> VertexIterWrapper;
//	typedef IteratorWrapperT<PolyConnectivity::HalfedgeIter, &ArrayKernel::n_halfedges> HalfedgeIterWrapper;
//	typedef IteratorWrapperT<PolyConnectivity::EdgeIter, &ArrayKernel::n_edges> EdgeIterWrapper;
//	typedef IteratorWrapperT<PolyConnectivity::FaceIter, &ArrayKernel::n_faces> FaceIterWrapper;

//	expose_property_manager<VPropHandleT<object>, VertexHandle, VertexIterWrapper>("VPropertyManager");
//	expose_property_manager<HPropHandleT<object>, HalfedgeHandle, HalfedgeIterWrapper>("HPropertyManager");
//	expose_property_manager<EPropHandleT<object>, EdgeHandle, EdgeIterWrapper>("EPropertyManager");
//	expose_property_manager<FPropHandleT<object>, FaceHandle, FaceIterWrapper>("FPropertyManager");

//	expose_io();

//	expose_decimater<PolyMesh>("PolyMesh");
//	expose_decimater<TriMesh>("TriMesh");
}

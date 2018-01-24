#include "MeshTypes.hh"
#include "Miscellaneous.hh"
#include "Vector.hh"
#include "Mesh.hh"
#include "Iterator.hh"
#include "Circulator.hh"
#include "PropertyManager.hh"
#include "InputOutput.hh"

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace OM = OpenMesh;


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

	expose_mesh<PolyMesh>(m, "PolyMesh");
	expose_mesh<TriMesh>(m, "TriMesh");

	expose_iterator<OM::PolyConnectivity::VertexIter, &OM::ArrayKernel::n_vertices>(m, "VertexIter");
	expose_iterator<OM::PolyConnectivity::HalfedgeIter, &OM::ArrayKernel::n_halfedges>(m, "HalfedgeIter");
	expose_iterator<OM::PolyConnectivity::EdgeIter, &OM::ArrayKernel::n_edges>(m, "EdgeIter");
	expose_iterator<OM::PolyConnectivity::FaceIter, &OM::ArrayKernel::n_faces>(m, "FaceIter");

	expose_circulator<OM::PolyConnectivity::VertexVertexIter, OM::VertexHandle>(m, "VertexVertexIter");
	expose_circulator<OM::PolyConnectivity::VertexIHalfedgeIter, OM::VertexHandle>(m, "VertexIHalfedgeIter");
	expose_circulator<OM::PolyConnectivity::VertexOHalfedgeIter, OM::VertexHandle>(m, "VertexOHalfedgeIter");
	expose_circulator<OM::PolyConnectivity::VertexEdgeIter, OM::VertexHandle>(m, "VertexEdgeIter");
	expose_circulator<OM::PolyConnectivity::VertexFaceIter, OM::VertexHandle>(m, "VertexFaceIter");

	expose_circulator<OM::PolyConnectivity::FaceVertexIter, OM::FaceHandle>(m, "FaceVertexIter");
	expose_circulator<OM::PolyConnectivity::FaceHalfedgeIter, OM::FaceHandle>(m, "FaceHalfedgeIter");
	expose_circulator<OM::PolyConnectivity::FaceEdgeIter, OM::FaceHandle>(m, "FaceEdgeIter");
	expose_circulator<OM::PolyConnectivity::FaceFaceIter, OM::FaceHandle>(m, "FaceFaceIter");

	expose_circulator<OM::PolyConnectivity::HalfedgeLoopIter, OM::HalfedgeHandle>(m, "HalfedgeLoopIter");

	typedef IteratorWrapperT<OM::PolyConnectivity::VertexIter, &OM::ArrayKernel::n_vertices> VertexIterWrapper;
	typedef IteratorWrapperT<OM::PolyConnectivity::HalfedgeIter, &OM::ArrayKernel::n_halfedges> HalfedgeIterWrapper;
	typedef IteratorWrapperT<OM::PolyConnectivity::EdgeIter, &OM::ArrayKernel::n_edges> EdgeIterWrapper;
	typedef IteratorWrapperT<OM::PolyConnectivity::FaceIter, &OM::ArrayKernel::n_faces> FaceIterWrapper;

	expose_property_manager<OM::VPropHandleT<py::object>, OM::VertexHandle, VertexIterWrapper>(m, "VPropertyManager");
	expose_property_manager<OM::HPropHandleT<py::object>, OM::HalfedgeHandle, HalfedgeIterWrapper>(m, "HPropertyManager");
	expose_property_manager<OM::EPropHandleT<py::object>, OM::EdgeHandle, EdgeIterWrapper>(m, "EPropertyManager");
	expose_property_manager<OM::FPropHandleT<py::object>, OM::FaceHandle, FaceIterWrapper>(m, "FPropertyManager");

	expose_io(m);

//	expose_decimater<PolyMesh>("PolyMesh");
//	expose_decimater<TriMesh>("TriMesh");
}

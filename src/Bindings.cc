#include "MeshTypes.hh"
#include "Miscellaneous.hh"
#include "Mesh.hh"
#include "Iterator.hh"
#include "Circulator.hh"
#include "InputOutput.hh"
#include "Decimater.hh"

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace OM = OpenMesh;


PYBIND11_MODULE(openmesh, m) {
	expose_handles(m);

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

	expose_io(m);

	expose_decimater<PolyMesh>(m, "PolyMesh");
	expose_decimater<TriMesh>(m, "TriMesh");
}

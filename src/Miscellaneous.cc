#include "Miscellaneous.hh"

#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh/Core/Utils/Property.hh>

namespace OM = OpenMesh;


/**
 * Expose item and property handles to %Python.
 */
void expose_handles(py::module& m) {
	py::class_<OM::BaseHandle>(m, "BaseHandle")
		.def(py::init<>())
		.def(py::init<int>())
		.def("idx", &OM::BaseHandle::idx)
		.def("is_valid", &OM::BaseHandle::is_valid)
		.def("reset", &OM::BaseHandle::reset)
		.def("invalidate", &OM::BaseHandle::invalidate)
		.def("__eq__", &OM::BaseHandle::operator ==)
		.def("__ne__", &OM::BaseHandle::operator !=)
		.def("__lt__", &OM::BaseHandle::operator <)
		;

	py::class_<OM::VertexHandle, OM::BaseHandle>(m, "VertexHandle")
		.def(py::init<>())
		.def(py::init<int>())
		;
	py::class_<OM::HalfedgeHandle, OM::BaseHandle>(m, "HalfedgeHandle")
		.def(py::init<>())
		.def(py::init<int>())
		;
	py::class_<OM::EdgeHandle, OM::BaseHandle>(m, "EdgeHandle")
		.def(py::init<>())
		.def(py::init<int>())
		;
	py::class_<OM::FaceHandle, OM::BaseHandle>(m, "FaceHandle")
		.def(py::init<>())
		.def(py::init<int>())
		;
}

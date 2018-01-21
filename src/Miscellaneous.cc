#include "Miscellaneous.hh"

#include <OpenMesh/Core/Mesh/ArrayItems.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh/Core/Utils/Property.hh>

namespace OM = OpenMesh;


/**
 * Expose mesh items to %Python.
 */
void expose_items(py::module& m) {
	py::class_<OM::ArrayItems::Vertex>(m, "Vertex");
	py::class_<OM::ArrayItems::Halfedge>(m, "Halfedge");
	py::class_<OM::ArrayItems::Edge>(m, "Edge");
	py::class_<OM::ArrayItems::Face>(m, "Face");
}

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

	py::class_<OM::BasePropHandleT<py::object>, OM::BaseHandle>(m, "BasePropHandle")
		.def(py::init<>())
		.def(py::init<int>())
		;

	py::class_<OM::VPropHandleT<py::object>, OM::BasePropHandleT<py::object> >(m, "VPropHandle")
		.def(py::init<>())
		.def(py::init<int>())
		.def(py::init<const OM::BasePropHandleT<py::object>&>())
		;
	py::class_<OM::HPropHandleT<py::object>, OM::BasePropHandleT<py::object> >(m, "HPropHandle")
		.def(py::init<>())
		.def(py::init<int>())
		.def(py::init<const OM::BasePropHandleT<py::object>&>())
		;
	py::class_<OM::EPropHandleT<py::object>, OM::BasePropHandleT<py::object> >(m, "EPropHandle")
		.def(py::init<>())
		.def(py::init<int>())
		.def(py::init<const OM::BasePropHandleT<py::object>&>())
		;
	py::class_<OM::FPropHandleT<py::object>, OM::BasePropHandleT<py::object> >(m, "FPropHandle")
		.def(py::init<>())
		.def(py::init<int>())
		.def(py::init<const OM::BasePropHandleT<py::object>&>())
		;
	py::class_<OM::MPropHandleT<py::object>, OM::BasePropHandleT<py::object> >(m, "MPropHandle")
		.def(py::init<>())
		.def(py::init<int>())
		.def(py::init<const OM::BasePropHandleT<py::object>&>())
		;
}


/**
 * Expose the StatusBits enum and StatusInfo class to %Python.
 */
void expose_status_bits_and_info(py::module& m) {
	using OM::Attributes::StatusBits;
	using OM::Attributes::StatusInfo;

	py::enum_<StatusBits>(m, "StatusBits", py::arithmetic())
		.value("DELETED", OM::Attributes::DELETED)
		.value("LOCKED", OM::Attributes::LOCKED)
		.value("SELECTED", OM::Attributes::SELECTED)
		.value("HIDDEN", OM::Attributes::HIDDEN)
		.value("FEATURE", OM::Attributes::FEATURE)
		.value("TAGGED", OM::Attributes::TAGGED)
		.value("TAGGED2", OM::Attributes::TAGGED2)
		.value("FIXEDNONMANIFOLD", OM::Attributes::FIXEDNONMANIFOLD)
		.value("UNUSED", OM::Attributes::UNUSED)
		;

	py::class_<StatusInfo>(m, "StatusInfo")
		.def("deleted", &StatusInfo::deleted)
		.def("set_deleted", &StatusInfo::set_deleted)
		.def("locked", &StatusInfo::locked)
		.def("set_locked", &StatusInfo::set_locked)
		.def("selected", &StatusInfo::selected)
		.def("set_selected", &StatusInfo::set_selected)
		.def("hidden", &StatusInfo::hidden)
		.def("set_hidden", &StatusInfo::set_hidden)
		.def("feature", &StatusInfo::feature)
		.def("set_feature", &StatusInfo::set_feature)
		.def("tagged", &StatusInfo::tagged)
		.def("set_tagged", &StatusInfo::set_tagged)
		.def("tagged2", &StatusInfo::tagged2)
		.def("set_tagged2", &StatusInfo::set_tagged2)
		.def("fixed_nonmanifold", &StatusInfo::fixed_nonmanifold)
		.def("set_fixed_nonmanifold", &StatusInfo::set_fixed_nonmanifold)
		.def("bits", &StatusInfo::bits)
		.def("set_bits", &StatusInfo::set_bits)
		.def("is_bit_set", &StatusInfo::is_bit_set)
		.def("set_bit", &StatusInfo::set_bit)
		.def("unset_bit", &StatusInfo::unset_bit)
		.def("change_bit", &StatusInfo::change_bit)
		;
}

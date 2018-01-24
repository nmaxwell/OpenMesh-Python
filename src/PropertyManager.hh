#ifndef OPENMESH_PYTHON_PROPERTYMANAGER_HH
#define OPENMESH_PYTHON_PROPERTYMANAGER_HH

#include "MeshTypes.hh"
#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace OM = OpenMesh;


/**
 * Implementation of %Python's \_\_getitem\_\_ magic method.
 *
 * @tparam PropertyManager A property manager type.
 * @tparam IndexHandle The appropriate handle type.
 *
 * @param _self The property manager instance that is to be used.
 * @param _handle The index of the property value to be returned.
 *
 * @return The requested property value.
 */
template <class PropertyManager, class IndexHandle>
py::object propman_get_item(PropertyManager& _self, IndexHandle _handle) {
	return _self[_handle];
}

/**
 * Implementation of %Python's \_\_setitem\_\_ magic method.
 *
 * @tparam PropertyManager A property manager type.
 * @tparam IndexHandle The appropriate handle type.
 *
 * @param _self The property manager instance that is to be used.
 * @param _handle The index of the property value to be set.
 * @param _value The property value to be set.
 */
template <class PropertyManager, class IndexHandle>
void propman_set_item(PropertyManager& _self, IndexHandle _handle, py::object _value) {
	_self[_handle] = _value;
}

/**
 * Conveniently set the property value for an entire range of mesh items
 * using a %Python iterator.
 *
 * @tparam PropertyManager A property manager type.
 * @tparam Iterator A %Python iterator type.
 *
 * @param _self The property manager instance that is to be used.
 * @param _it An iterator that iterates over the items in the range.
 * @param _value The value the range will be set to.
 */
template <class PropertyManager, class Iterator>
void propman_set_range(PropertyManager& _self, Iterator _it, py::object _value) {
	try {
		while (true) {
			_self[_it.next()] = _value;
		}
	}
	catch (py::stop_iteration exception) {
		// This is expected behavior
		PyErr_Clear();
	}
}

/**
 * Thin wrapper for propertyExists.
 *
 * @tparam PropertyManager A property manager type.
 * @tparam Mesh A mesh type.
 *
 * @param _mesh The mesh that is used to check if the property exists.
 * @param _propname The name of the property.
 */
template <class PropertyManager, class Mesh>
bool property_exists(Mesh& _mesh, const char *_propname) {
	return PropertyManager::propertyExists(_mesh, _propname);
}

/**
 * Expose a property manager type to %Python.
 *
 * This function template is used to expose property managers to %Python. The
 * template parameters are used to instantiate the appropriate property manager
 * type.
 *
 * @tparam PropHandle A property handle type (e.g. %VPropHandle\<object\>).
 * @tparam IndexHandle The appropriate handle type (e.g. %VertexHandle for
 * %VPropHandle\<object\>).
 * @tparam Iterator A %Python iterator type. This type is used to instantiate
 * the propman_set_range function.
 *
 * @param _name The name of the property manager type to be exposed.
 */
template <class PropHandle, class IndexHandle, class Iterator>
void expose_property_manager(py::module& m, const char *_name) {
	// Convenience typedef
	typedef OM::PropertyManager<PropHandle, OM::PolyConnectivity> PropertyManager;

	// Function pointers
	py::object (*getitem)(PropertyManager&, IndexHandle            ) = &propman_get_item;
	void       (*setitem)(PropertyManager&, IndexHandle, py::object) = &propman_set_item;

	void (*set_range)(PropertyManager&, Iterator, py::object) = &propman_set_range;

	bool (*property_exists_poly)(PolyMesh&, const char *) = &property_exists<PropertyManager, PolyMesh>;
	bool (*property_exists_tri )(TriMesh&,  const char *) = &property_exists<PropertyManager, TriMesh >;

	// Expose property manager
	py::class_<PropertyManager>(m, _name)
		.def(py::init<PolyMesh&, const char *>(), py::keep_alive<1,2>())
		.def(py::init<PolyMesh&, const char *, bool>(), py::keep_alive<1,2>())
		.def(py::init<TriMesh&,  const char *>(), py::keep_alive<1,2>())
		.def(py::init<TriMesh&,  const char *, bool>(), py::keep_alive<1,2>())

		.def("swap", &PropertyManager::swap)
		.def("is_valid", &PropertyManager::isValid)

		.def("__bool__", &PropertyManager::operator bool)
		.def("__nonzero__", &PropertyManager::operator bool)

		.def("get_raw_property", &PropertyManager::getRawProperty, py::return_value_policy::copy)
		.def("get_name", &PropertyManager::getName, py::return_value_policy::copy)
		.def("get_mesh", &PropertyManager::getMesh, py::return_value_policy::reference)

		.def("retain", &PropertyManager::retain, py::arg("do_retain")=false)

		.def("__getitem__", getitem)
		.def("__setitem__", setitem)

		.def("set_range", set_range)

		.def_static("property_exists", property_exists_poly)
		.def_static("property_exists", property_exists_tri)
		;
}

#endif

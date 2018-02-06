/** @file */

#ifndef OPENMESH_PYTHON_MESHTYPES_HH
#define OPENMESH_PYTHON_MESHTYPES_HH

#define OM_STATIC_BUILD

#include "Utilities.hh"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


struct MeshTraits : public OpenMesh::DefaultTraits {
	/** Use double precision points */
	typedef OpenMesh::Vec3d Point;

	/** Use double precision normals */
	typedef OpenMesh::Vec3d Normal;

	/** Use RGBA colors */
	typedef OpenMesh::Vec4f Color;

	/** Use double precision texcoords */
	typedef double TexCoord1D;
	typedef OpenMesh::Vec2d TexCoord2D;
	typedef OpenMesh::Vec3d TexCoord3D;
};


template <class Mesh>
class MeshWrapperT : public Mesh {
public:

	typedef OpenMesh::VPropHandleT<py::none> VPropHandle;
	typedef OpenMesh::HPropHandleT<py::none> HPropHandle;
	typedef OpenMesh::EPropHandleT<py::none> EPropHandle;
	typedef OpenMesh::FPropHandleT<py::none> FPropHandle;

	template <class Handle, class PropHandle>
	py::none py_property(const std::string& _name, Handle _h) {
		const auto prop = py_prop_on_demand<Handle, PropHandle>(_name);
		return Mesh::property(prop, _h);
	}

	template <class Handle, class PropHandle>
	void py_set_property(const std::string& _name, Handle _h, py::object _val) {
		const auto prop = py_prop_on_demand<Handle, PropHandle>(_name);
		Mesh::property(prop, _h) = _val;
	}

	template <class Handle>
	bool py_has_property(const std::string& _name) {
		auto& prop_map = py_prop_map(Handle());
		return prop_map.count(_name);
	}

	template <class Handle>
	void py_remove_property(const std::string& _name) {
		auto& prop_map = py_prop_map(Handle());
		if (prop_map.count(_name) != 0) {
			Mesh::remove_property(prop_map.at(_name));
			prop_map.erase(_name);
		}
	}

	template <class Handle, class PropHandle>
	py::list py_property_generic(const std::string& _name) {
		const size_t n = py_n_items(Handle());
		const auto prop = py_prop_on_demand<Handle, PropHandle>(_name);

		py::list res;
		for (size_t i = 0; i < n; ++i) {
			res.append(Mesh::property(prop, Handle(i)));
		}
		return res;
	}

	template <class Handle, class PropHandle>
	void py_set_property_generic(const std::string& _name, py::list _list) {
		const size_t n = py_n_items(Handle());
		const auto prop = py_prop_on_demand<Handle, PropHandle>(_name);

		if (_list.size() != n) {
			return;
		}
		for (size_t i = 0; i < n; ++i) {
			Mesh::property(prop, Handle(i)) = py::object(_list[i]);
		}
	}

	template <class Handle, class PropHandle>
	py::array_t<double> py_property_array(const std::string& _name) {
		const size_t n = py_n_items(Handle());
		const auto prop = py_prop_on_demand<Handle, PropHandle>(_name);

		// assume that all arrays have the same size and
		// retrieve the size of the first array
		const py::object tmp_obj = Mesh::property(prop, Handle(0));
		py::array_t<double> tmp_arr;
		try {
			tmp_arr = tmp_obj.cast<py::array_t<double> >();
		}
		catch (py::error_already_set& e) {
			return py::array_t<double>();
		}
		const size_t size = tmp_arr.size();

		// better check this now
		if (size == 0) {
			return py::array_t<double>();
		}

		// allocate memory
		double *data = new double[size * n];

		// copy one array at a time
		for (size_t i = 0; i < n; ++i) {
			const Handle hnd(i);
			const py::object obj = Mesh::property(prop, hnd);
			try {
				const auto arr = make_c_style(obj.cast<py::array_t<double> >());
				if (arr.size() != size) {
					throw py::error_already_set();
				}
				std::copy(arr.data(0), arr.data(0) + size, &data[size * i]);
			}
			catch (py::error_already_set& e) {
				delete[] data;
				return py::array_t<double>();
			}
		}

		// make numpy array
		const auto shape = {n, size};
		const auto strides = {size * sizeof(double), sizeof(double)};
		py::capsule base = free_when_done(data);
		return py::array_t<double>(shape, strides, data, base);
	}

	template <class Handle, class PropHandle>
	void py_set_property_array(const std::string& _name, py::array_t<double, py::array::c_style | py::array::forcecast> _arr) {
		const size_t n = py_n_items(Handle());
		const auto prop = py_prop_on_demand<Handle, PropHandle>(_name);

		// array cannot be empty and its shape has to be (_n, m,...)
		if (_arr.size() == 0 || _arr.ndim() < 2 || _arr.shape(0) != n) {
			return;
		}

		// copy one array at a time
		const size_t size = _arr.strides(0) / sizeof(double);
		for (size_t i = 0; i < n; ++i) {
			double *data = new double[size];
			std::copy(_arr.data(i), _arr.data(i) + size, data);
			const auto shape = {size};
			const auto strides = {sizeof(double)};
			py::capsule base = free_when_done(data);
			py::array_t<double> tmp(shape, strides, data, base);
			Mesh::property(prop, Handle(i)) = tmp;
		}
	}

	template <class Handle, class PropHandle>
	void py_copy_property(const std::string& _name, Handle _from, Handle _to) {
		auto prop = py_prop_on_demand<Handle, PropHandle>(_name);
		Mesh::copy_property(prop, _from, _to);
	}

	size_t py_n_items(OpenMesh::VertexHandle) const { return Mesh::n_vertices(); }
	size_t py_n_items(OpenMesh::HalfedgeHandle) const { return Mesh::n_halfedges(); }
	size_t py_n_items(OpenMesh::EdgeHandle) const { return Mesh::n_edges(); }
	size_t py_n_items(OpenMesh::FaceHandle) const { return Mesh::n_faces(); }

private:

	template <class Handle, class PropHandle>
	PropHandle py_prop_on_demand(const std::string& _name) {
		auto& prop_map = py_prop_map(Handle());
		if (prop_map.count(_name) == 0) {
			PropHandle prop;
			Mesh::add_property(prop, _name);
			prop_map[_name] = prop;
		}
		return prop_map.at(_name);
	}

	std::map<std::string, VPropHandle>& py_prop_map(OpenMesh::VertexHandle) { return vprop_map; }
	std::map<std::string, HPropHandle>& py_prop_map(OpenMesh::HalfedgeHandle) { return hprop_map; }
	std::map<std::string, EPropHandle>& py_prop_map(OpenMesh::EdgeHandle) { return eprop_map; }
	std::map<std::string, FPropHandle>& py_prop_map(OpenMesh::FaceHandle) { return fprop_map; }

	std::map<std::string, VPropHandle> vprop_map;
	std::map<std::string, HPropHandle> hprop_map;
	std::map<std::string, EPropHandle> eprop_map;
	std::map<std::string, FPropHandle> fprop_map;
};


typedef MeshWrapperT<OpenMesh::TriMesh_ArrayKernelT<MeshTraits> > TriMesh;
typedef MeshWrapperT<OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> > PolyMesh;

#endif

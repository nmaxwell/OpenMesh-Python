#ifndef OPENMESH_PYTHON_MESH_HH
#define OPENMESH_PYTHON_MESH_HH

#include "MeshTypes.hh"
#include "Iterator.hh"
#include "Circulator.hh"

#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace OM = OpenMesh;


/**
 * Return value policy for functions that return references to objects that are
 * managed by %OpenMesh.
 */
#define OPENMESH_PYTHON_DEFAULT_POLICY py::return_value_policy::copy

/**
 * Set the status of an item.
 *
 * @tparam Mesh A mesh type.
 * @tparam PropHandle A handle type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _h The handle of the item whose status is to be set.
 * @param _info The status to be set.
 *
 * Depending on @ref OPENMESH_PYTHON_DEFAULT_POLICY, Mesh::status may
 * return by value instead of reference. This function ensures that the
 * status of an item can be changed nonetheless.
 */
template <class Mesh, class IndexHandle>
void set_status(Mesh& _self, IndexHandle _h, const OpenMesh::Attributes::StatusInfo& _info) {
	_self.status(_h) = _info;
}

/**
 * Set the value of a property of an item.
 *
 * @tparam Mesh A mesh type.
 * @tparam PropHandle A property handle type.
 * @tparam IndexHandle The appropriate handle type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _ph The property that is to be set.
 * @param _h The handle of the item whose property is to be set.
 * @param _value The value to be set.
 *
 * Depending on @ref OPENMESH_PYTHON_DEFAULT_POLICY, Mesh::property may
 * return by value instead of reference. This function ensures that the
 * property value of an item can be changed nonetheless.
 */
template <class Mesh, class PropHandle, class IndexHandle>
void set_property(Mesh& _self, PropHandle _ph, IndexHandle _h, const py::object& _value) {
	_self.property(_ph, _h) = _value;
}

/**
 * Set the value of a mesh property.
 *
 * @tparam Mesh A mesh type.
 * @tparam PropHandle A property handle type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _ph The property that is to be set.
 * @param _value The value to be set.
 *
 * Depending on @ref OPENMESH_PYTHON_DEFAULT_POLICY, Mesh::property may
 * return by value instead of reference. This function ensures that the
 * property value of an item can be changed nonetheless.
 */
template <class Mesh, class PropHandle>
void set_property(Mesh& _self, PropHandle _ph, const py::object& _value) {
	_self.property(_ph) = _value;
}

/**
 * Thin wrapper for assign_connectivity.
 *
 * @tparam Mesh A mesh type.
 * @tparam OtherMesh A mesh type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _other The mesh from which the connectivity is to be copied.
 */
template <class Mesh, class OtherMesh>
void assign_connectivity(Mesh& _self, const OtherMesh& _other) {
	_self.assign_connectivity(_other);
}

/**
 * Get an iterator.
 */
template <class Mesh, class Iterator, size_t (OM::ArrayKernel::*n_items)() const>
IteratorWrapperT<Iterator, n_items> get_iterator(Mesh& _self) {
	return IteratorWrapperT<Iterator, n_items>(_self, typename Iterator::value_type(0));
}

/**
 * Get a skipping iterator.
 */
template <class Mesh, class Iterator, size_t (OM::ArrayKernel::*n_items)() const>
IteratorWrapperT<Iterator, n_items> get_skipping_iterator(Mesh& _self) {
	return IteratorWrapperT<Iterator, n_items>(_self, typename Iterator::value_type(0), true);
}

/**
 * Get a circulator.
 *
 * @tparam Mesh A Mesh type.
 * @tparam Circulator A circulator type.
 * @tparam CenterEntityHandle The appropriate handle type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _handle The handle of the item to circulate around.
 */
template <class Mesh, class Circulator, class CenterEntityHandle>
CirculatorWrapperT<Circulator, CenterEntityHandle> get_circulator(Mesh& _self, CenterEntityHandle _handle) {
	return CirculatorWrapperT<Circulator, CenterEntityHandle>(_self, _handle);
}

/**
 * Garbage collection using lists instead of vectors to keep track of a set of
 * handles.
 *
 * @tparam Mesh A Mesh type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _vh_to_update The list of vertex handles to be updated.
 * @param _hh_to_update The list of halfedge handles to be updated.
 * @param _fh_to_update The list of face handles to be updated.
 * @param _v Remove deleted vertices?
 * @param _e Remove deleted edges?
 * @param _f Remove deleted faces?
 */
template <class Mesh>
void garbage_collection(Mesh& _self, py::list& _vh_to_update, py::list& _hh_to_update, py::list& _fh_to_update, bool _v = true, bool _e = true, bool _f = true) {
	// Convert list of handles to vector of pointers
	std::vector<OM::VertexHandle*> vh_vector;
	for (auto item : _vh_to_update) {
		if (py::isinstance<OM::VertexHandle>(item)) {
			vh_vector.push_back(item.cast<OM::VertexHandle*>());
		}
	}

	// Convert list of handles to vector of pointers
	std::vector<OM::HalfedgeHandle*> hh_vector;
	for (auto item : _hh_to_update) {
		if (py::isinstance<OM::HalfedgeHandle>(item)) {
			hh_vector.push_back(item.cast<OM::HalfedgeHandle*>());
		}
	}

	// Convert list of handles to vector of pointers
	std::vector<OM::FaceHandle*> fh_vector;
	for (auto item : _fh_to_update) {
		if (py::isinstance<OM::FaceHandle>(item)) {
			fh_vector.push_back(item.cast<OM::FaceHandle*>());
		}
	}

	// Call garbage collection
	_self.garbage_collection(vh_vector, hh_vector, fh_vector, _v, _e, _f);
}

/**
 * Add a new face from a %Python list of vertex handles.
 *
 * @tparam Mesh A Mesh type.
 *
 * @param _self The mesh instance that is to be used.
 * @param _vhandles The list of vertex handles.
 */
template<class Mesh>
OM::FaceHandle add_face(Mesh& _self, const py::list& _vhandles) {
	// TODO Use automatic conversion instead?
	std::vector<OM::VertexHandle> vector;
	for (auto item : _vhandles) {
		if (py::isinstance<OM::VertexHandle>(item)) {
			vector.push_back(item.cast<OM::VertexHandle>());
		}
	}
	return _self.add_face(vector);
}

/**
 * Converts OpenMesh vectors to numpy arrays.
 *
 * The returned array references the vector's underlying data, i.e. changes
 * made to the returned array affect the original mesh.
 *
 * @tparam Mesh A Mesh type.
 * @tparam vector A Vector type.
 *
 * @param _mesh The mesh that owns the vector's underlying memory. In order
 * to avaoid dangling pointers, the lifetime of this mesh is tied to the
 * lifetime of the returned numpy array.
 * @param _vec The vector to be converted.
 */
template<class Mesh, class Vector>
py::array_t<typename Vector::value_type> vec2numpy(Mesh& _mesh, Vector& _vec, size_t _n = 1) {
	typedef typename Vector::value_type dtype;
	std::vector<size_t> shape;
	std::vector<size_t> strides;
	if (_n == 1) {
		shape = {_vec.size()};
		strides = {sizeof(dtype)};
	}
	else {
		shape = {_n, _vec.size()};
		strides = {_vec.size() * sizeof(dtype), sizeof(dtype)};
	}
	return py::array_t<dtype>(shape, strides, _vec.data(), py::cast(_mesh));
}

template<class Mesh>
py::array_t<float> flt2numpy(Mesh& _mesh, const float& _flt, size_t _n = 1) {
	return py::array_t<float>({_n}, {sizeof(float)}, &_flt, py::cast(_mesh));
}

template<class Mesh>
py::array_t<double> flt2numpy(Mesh& _mesh, const double& _flt, size_t _n = 1) {
	return py::array_t<double>({_n}, {sizeof(double)}, &_flt, py::cast(_mesh));
}

/**
 * Returns the face indices of this mesh as a numpy array.
 *
 * Note that the array is constructed on the fly and does not
 * reference the underlying mesh.
 */
py::array_t<int> face_indices_tri(TriMesh& _self) {
	if (_self.n_faces() == 0) {
		return py::array_t<int>();
	}
	int *indices = new int[_self.n_faces() * 3];
	for (auto fh : _self.faces()) {
		auto fv_it = _self.fv_iter(fh);
		indices[fh.idx() * 3 + 0] = fv_it->idx(); ++fv_it;
		indices[fh.idx() * 3 + 1] = fv_it->idx(); ++fv_it;
		indices[fh.idx() * 3 + 2] = fv_it->idx();
	}
	const auto shape = {_self.n_faces(), size_t(3)};
	const auto strides = {3 * sizeof(int), sizeof(int)};
	py::capsule free_when_done(indices, [](void *f) {
		int *ptr = reinterpret_cast<int *>(f);
		delete[] ptr;
	});
	return py::array_t<int>(shape, strides, indices, free_when_done);
}

/**
 * Returns the face indices of this mesh as a numpy array.
 *
 * If the faces of the mesh have different valences, the array
 * is padded with -1 entries.
 *
 * Note that the array is constructed on the fly and does not
 * reference the underlying mesh.
 */
py::array_t<int> face_indices_poly(PolyMesh& _self) {
	if (_self.n_faces() == 0) {
		return py::array_t<int>();
	}
	int max_valence = 0;
	for (auto fh : _self.faces()) {
		max_valence = std::max(max_valence, int(_self.valence(fh)));
	}
	int *indices = new int[_self.n_faces() * max_valence];
	for (auto fh : _self.faces()) {
		auto fv_it = _self.fv_iter(fh);
		const int valence = _self.valence(fh);
		for (int i = 0; i < valence; ++i) {
			indices[fh.idx() * max_valence + i] = fv_it->idx();
			++fv_it;
		}
		for (int i = valence; i < max_valence; ++i) {
			indices[fh.idx() * max_valence + i] = -1;
		}
	}
	const auto shape = {_self.n_faces(), size_t(max_valence)};
	const auto strides = {max_valence * sizeof(int), sizeof(int)};
	py::capsule free_when_done(indices, [](void *f) {
		int *ptr = reinterpret_cast<int *>(f);
		delete[] ptr;
	});
	return py::array_t<int>(shape, strides, indices, free_when_done);
}

/**
 * Attempts to return a custom property for all mesh items at once using a
 * numpy array. Returns an empty array if the property contains elements that
 * are not numpy arrays or at least one of the arrays has a non-uniform shape.
 *
 * @note The returned array is constructed on the fly and contains copies
 * of the actual property values.
 */
template <class Mesh, class PropHandle, class IndexHandle>
py::array_t<double> property_array(Mesh& _self, PropHandle _ph, size_t _n) {
	// assume that all arrays have the same size and
	// retrieve the size of the first array
	const py::object tmp_obj = _self.property(_ph, IndexHandle(0));
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
	double *data = new double[size * _n];

	// copy one array at a time
	for (size_t i = 0; i < _n; ++i) {
		const IndexHandle ih(i);
		const py::object obj = _self.property(_ph, ih);
		try {
			const auto arr = obj.cast<py::array_t<double> >();
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
	const auto shape = {_n, size};
	const auto strides = {size * sizeof(double), sizeof(double)};
	py::capsule free_when_done(data, [](void *f) {
		double *ptr = reinterpret_cast<double *>(f);
		delete[] ptr;
	});
	return py::array_t<double>(shape, strides, data, free_when_done);
}

/**
 * Attempts to set a custom property for all mesh items at once using a
 * numpy array.
 *
 * @note The property is set to copies of slices of _arr.
 */
template <class Mesh, class PropHandle, class IndexHandle>
void set_property_array(Mesh& _self, PropHandle _ph, py::array_t<double> _arr, size_t _n) {
	// array cannot be empty and its shape has to be (_n, m,...)
	if (_arr.size() == 0 || _arr.ndim() < 2 || _arr.shape(0) != _n) {
		return;
	}

	// copy one array at a time
	const size_t size = _arr.strides(0) / sizeof(double);
	for (size_t i = 0; i < _n; ++i) {
		double *data = new double[size];
		std::copy(_arr.data(i), _arr.data(i) + size, data);
		const auto shape = {size};
		const auto strides = {sizeof(double)};
		py::capsule free_when_done(data, [](void *f) {
			double *ptr = reinterpret_cast<double *>(f);
			delete[] ptr;
		});
		py::array_t<double> tmp(shape, strides, data, free_when_done);
		_self.property(_ph, IndexHandle(i)) = tmp;
	}
}

/**
 * This function template is used to expose mesh member functions that are only
 * available for a specific type of mesh (i.e. they are available for polygon
 * meshes or triangle meshes, but not both).
 *
 * @tparam Class A pybind11::class type.
 *
 * @param _class The pybind11::class instance for which the member
 * functions are to be defined.
 */
template <class Class>
void expose_type_specific_functions(Class& _class) {
	// See the template specializations below
}

/**
 * Function template specialization for polygon meshes.
 */
template <>
void expose_type_specific_functions(py::class_<PolyMesh>& _class) {
	typedef PolyMesh::Scalar Scalar;
	typedef PolyMesh::Point  Point;
	typedef PolyMesh::Normal Normal;
	typedef PolyMesh::Color  Color;

	OM::FaceHandle (PolyMesh::*add_face_3_vh)(OM::VertexHandle, OM::VertexHandle, OM::VertexHandle                  ) = &PolyMesh::add_face;
	OM::FaceHandle (PolyMesh::*add_face_4_vh)(OM::VertexHandle, OM::VertexHandle, OM::VertexHandle, OM::VertexHandle) = &PolyMesh::add_face;
	OM::FaceHandle (*add_face_list)(PolyMesh&, const py::list&) = &add_face;

	void (PolyMesh::*split_eh_pt)(OM::EdgeHandle, const Point&    ) = &PolyMesh::split;
	void (PolyMesh::*split_eh_vh)(OM::EdgeHandle, OM::VertexHandle) = &PolyMesh::split;
	void (PolyMesh::*split_fh_pt)(OM::FaceHandle, const Point&    ) = &PolyMesh::split;
	void (PolyMesh::*split_fh_vh)(OM::FaceHandle, OM::VertexHandle) = &PolyMesh::split;

	Normal (PolyMesh::*calc_face_normal_pt)(const Point&, const Point&, const Point&) const = &PolyMesh::calc_face_normal;

	_class
		.def("add_face", add_face_3_vh)
		.def("add_face", add_face_4_vh)
		.def("add_face", add_face_list)

		.def("split", split_eh_pt)
		.def("split", split_eh_vh)
		.def("split", split_fh_pt)
		.def("split", split_fh_vh)

		.def("split_copy", &PolyMesh::split_copy)
		.def("calc_face_normal", calc_face_normal_pt)
		.def("insert_edge", &PolyMesh::insert_edge)

		.def("face_indices", &face_indices_poly)
		;
}

/**
 * Function template specialization for triangle meshes.
 */
template <>
void expose_type_specific_functions(py::class_<TriMesh>& _class) {
	typedef TriMesh::Scalar Scalar;
	typedef TriMesh::Point  Point;
	typedef TriMesh::Normal Normal;
	typedef TriMesh::Color  Color;

	OM::FaceHandle (TriMesh::*add_face_3_vh)(OM::VertexHandle, OM::VertexHandle, OM::VertexHandle) = &TriMesh::add_face;
	OM::FaceHandle (*add_face_list)(TriMesh&, const py::list&) = &add_face;

	OM::VertexHandle (TriMesh::*split_eh_pt)(OM::EdgeHandle, const Point&    ) = &TriMesh::split;
	void             (TriMesh::*split_eh_vh)(OM::EdgeHandle, OM::VertexHandle) = &TriMesh::split;
	OM::VertexHandle (TriMesh::*split_fh_pt)(OM::FaceHandle, const Point&    ) = &TriMesh::split;
	void             (TriMesh::*split_fh_vh)(OM::FaceHandle, OM::VertexHandle) = &TriMesh::split;

	OM::VertexHandle (TriMesh::*split_copy_eh_pt)(OM::EdgeHandle, const Point&    ) = &TriMesh::split_copy;
	void             (TriMesh::*split_copy_eh_vh)(OM::EdgeHandle, OM::VertexHandle) = &TriMesh::split_copy;
	OM::VertexHandle (TriMesh::*split_copy_fh_pt)(OM::FaceHandle, const Point&    ) = &TriMesh::split_copy;
	void             (TriMesh::*split_copy_fh_vh)(OM::FaceHandle, OM::VertexHandle) = &TriMesh::split_copy;

	OM::HalfedgeHandle (TriMesh::*vertex_split_pt)(Point,            OM::VertexHandle, OM::VertexHandle, OM::VertexHandle) = &TriMesh::vertex_split;
	OM::HalfedgeHandle (TriMesh::*vertex_split_vh)(OM::VertexHandle, OM::VertexHandle, OM::VertexHandle, OM::VertexHandle) = &TriMesh::vertex_split;

	_class
		.def("add_face", add_face_3_vh)
		.def("add_face", add_face_list)

		.def("split", split_eh_pt)
		.def("split", split_eh_vh)
		.def("split", split_fh_pt)
		.def("split", split_fh_vh)

		.def("split_copy", split_copy_eh_pt)
		.def("split_copy", split_copy_eh_vh)
		.def("split_copy", split_copy_fh_pt)
		.def("split_copy", split_copy_fh_vh)

		.def("opposite_vh", &TriMesh::opposite_vh)
		.def("opposite_he_opposite_vh", &TriMesh::opposite_he_opposite_vh)

		.def("vertex_split", vertex_split_pt)
		.def("vertex_split", vertex_split_vh)

		.def("is_flip_ok", &TriMesh::is_flip_ok)
		.def("flip", &TriMesh::flip)

		.def("face_indices", &face_indices_tri)
		;
}


/**
 * Expose a mesh type to %Python.
 *
 * @tparam Mesh A mesh type.
 *
 * @param _name The name of the mesh type to be exposed.
 */
template <class Mesh>
void expose_mesh(py::module& m, const char *_name) {
	using OpenMesh::Attributes::StatusInfo;

	typedef typename Mesh::Scalar Scalar;
	typedef typename Mesh::Point  Point;
	typedef typename Mesh::Normal Normal;
	typedef typename Mesh::Color  Color;

	typedef typename Mesh::TexCoord1D TexCoord1D;
	typedef typename Mesh::TexCoord2D TexCoord2D;
	typedef typename Mesh::TexCoord3D TexCoord3D;

	//======================================================================
	//  KernelT Function Pointers
	//======================================================================

	// Get the i'th item
	OM::VertexHandle   (Mesh::*vertex_handle_uint  )(unsigned int) const = &Mesh::vertex_handle;
	OM::HalfedgeHandle (Mesh::*halfedge_handle_uint)(unsigned int) const = &Mesh::halfedge_handle;
	OM::EdgeHandle     (Mesh::*edge_handle_uint    )(unsigned int) const = &Mesh::edge_handle;
	OM::FaceHandle     (Mesh::*face_handle_uint    )(unsigned int) const = &Mesh::face_handle;

	// Delete items
	void (Mesh::*garbage_collection_bools)(bool, bool, bool) = &Mesh::garbage_collection;
	void (*garbage_collection_lists_bools)(Mesh&, py::list&, py::list&, py::list&, bool, bool, bool) = &garbage_collection;

	// Vertex connectivity
	OM::HalfedgeHandle (Mesh::*halfedge_handle_vh)(OM::VertexHandle) const = &Mesh::halfedge_handle;
	OM::HalfedgeHandle (Mesh::*halfedge_handle_fh)(OM::FaceHandle  ) const = &Mesh::halfedge_handle;

	// Halfedge connectivity
	OM::FaceHandle     (Mesh::*face_handle_hh         )(OM::HalfedgeHandle) const = &Mesh::face_handle;
	OM::HalfedgeHandle (Mesh::*prev_halfedge_handle_hh)(OM::HalfedgeHandle) const = &Mesh::prev_halfedge_handle;
	OM::EdgeHandle     (Mesh::*edge_handle_hh         )(OM::HalfedgeHandle) const = &Mesh::edge_handle;

	// Edge connectivity
	OM::HalfedgeHandle (Mesh::*halfedge_handle_eh_uint)(OM::EdgeHandle, unsigned int) const = &Mesh::halfedge_handle;

	// Set halfedge
	void (Mesh::*set_halfedge_handle_vh_hh)(OM::VertexHandle, OM::HalfedgeHandle) = &Mesh::set_halfedge_handle;
	void (Mesh::*set_halfedge_handle_fh_hh)(OM::FaceHandle,   OM::HalfedgeHandle) = &Mesh::set_halfedge_handle;

	// Handle -> Item
	const typename Mesh::Vertex&   (Mesh::*vertex  )(OM::VertexHandle  ) const = &Mesh::vertex;
	const typename Mesh::Halfedge& (Mesh::*halfedge)(OM::HalfedgeHandle) const = &Mesh::halfedge;
	const typename Mesh::Edge&     (Mesh::*edge    )(OM::EdgeHandle    ) const = &Mesh::edge;
	const typename Mesh::Face&     (Mesh::*face    )(OM::FaceHandle    ) const = &Mesh::face;

	// Item -> Handle
	OM::VertexHandle   (Mesh::*handle_v)(const typename Mesh::Vertex&  ) const = &Mesh::handle;
	OM::HalfedgeHandle (Mesh::*handle_h)(const typename Mesh::Halfedge&) const = &Mesh::handle;
	OM::EdgeHandle     (Mesh::*handle_e)(const typename Mesh::Edge&    ) const = &Mesh::handle;
	OM::FaceHandle     (Mesh::*handle_f)(const typename Mesh::Face&    ) const = &Mesh::handle;

	// Get value of a standard property (point, normal, color)
	const typename Mesh::Point&  (Mesh::*point_vh )(OM::VertexHandle  ) const = &Mesh::point;
	const typename Mesh::Normal& (Mesh::*normal_vh)(OM::VertexHandle  ) const = &Mesh::normal;
	const typename Mesh::Normal& (Mesh::*normal_hh)(OM::HalfedgeHandle) const = &Mesh::normal;
	const typename Mesh::Normal& (Mesh::*normal_fh)(OM::FaceHandle    ) const = &Mesh::normal;
	const typename Mesh::Color&  (Mesh::*color_vh )(OM::VertexHandle  ) const = &Mesh::color;
	const typename Mesh::Color&  (Mesh::*color_hh )(OM::HalfedgeHandle) const = &Mesh::color;
	const typename Mesh::Color&  (Mesh::*color_eh )(OM::EdgeHandle    ) const = &Mesh::color;
	const typename Mesh::Color&  (Mesh::*color_fh )(OM::FaceHandle    ) const = &Mesh::color;

	// Get value of a standard property (texture coordinate)
	const typename Mesh::TexCoord1D& (Mesh::*texcoord1D_vh)(OM::VertexHandle  ) const = &Mesh::texcoord1D;
	const typename Mesh::TexCoord1D& (Mesh::*texcoord1D_hh)(OM::HalfedgeHandle) const = &Mesh::texcoord1D;
	const typename Mesh::TexCoord2D& (Mesh::*texcoord2D_vh)(OM::VertexHandle  ) const = &Mesh::texcoord2D;
	const typename Mesh::TexCoord2D& (Mesh::*texcoord2D_hh)(OM::HalfedgeHandle) const = &Mesh::texcoord2D;
	const typename Mesh::TexCoord3D& (Mesh::*texcoord3D_vh)(OM::VertexHandle  ) const = &Mesh::texcoord3D;
	const typename Mesh::TexCoord3D& (Mesh::*texcoord3D_hh)(OM::HalfedgeHandle) const = &Mesh::texcoord3D;

	// Get value of a standard property (status)
	const StatusInfo& (Mesh::*status_vh)(OM::VertexHandle  ) const = &Mesh::status;
	const StatusInfo& (Mesh::*status_hh)(OM::HalfedgeHandle) const = &Mesh::status;
	const StatusInfo& (Mesh::*status_eh)(OM::EdgeHandle    ) const = &Mesh::status;
	const StatusInfo& (Mesh::*status_fh)(OM::FaceHandle    ) const = &Mesh::status;

	// Set value of a standard property (point, normal, color)
	void (Mesh::*set_normal_vh)(OM::VertexHandle,   const typename Mesh::Normal&) = &Mesh::set_normal;
	void (Mesh::*set_normal_hh)(OM::HalfedgeHandle, const typename Mesh::Normal&) = &Mesh::set_normal;
	void (Mesh::*set_normal_fh)(OM::FaceHandle,     const typename Mesh::Normal&) = &Mesh::set_normal;
	void (Mesh::*set_color_vh )(OM::VertexHandle,   const typename Mesh::Color& ) = &Mesh::set_color;
	void (Mesh::*set_color_hh )(OM::HalfedgeHandle, const typename Mesh::Color& ) = &Mesh::set_color;
	void (Mesh::*set_color_eh )(OM::EdgeHandle,     const typename Mesh::Color& ) = &Mesh::set_color;
	void (Mesh::*set_color_fh )(OM::FaceHandle,     const typename Mesh::Color& ) = &Mesh::set_color;

	// Set value of a standard property (texture coordinate)
	void (Mesh::*set_texcoord1D_vh)(OM::VertexHandle,   const typename Mesh::TexCoord1D&) = &Mesh::set_texcoord1D;
	void (Mesh::*set_texcoord1D_hh)(OM::HalfedgeHandle, const typename Mesh::TexCoord1D&) = &Mesh::set_texcoord1D;
	void (Mesh::*set_texcoord2D_vh)(OM::VertexHandle,   const typename Mesh::TexCoord2D&) = &Mesh::set_texcoord2D;
	void (Mesh::*set_texcoord2D_hh)(OM::HalfedgeHandle, const typename Mesh::TexCoord2D&) = &Mesh::set_texcoord2D;
	void (Mesh::*set_texcoord3D_vh)(OM::VertexHandle,   const typename Mesh::TexCoord3D&) = &Mesh::set_texcoord3D;
	void (Mesh::*set_texcoord3D_hh)(OM::HalfedgeHandle, const typename Mesh::TexCoord3D&) = &Mesh::set_texcoord3D;

	// Set value of a standard property (status)
	void (*set_status_vh)(Mesh&, OM::VertexHandle,   const StatusInfo&) = &set_status;
	void (*set_status_hh)(Mesh&, OM::HalfedgeHandle, const StatusInfo&) = &set_status;
	void (*set_status_eh)(Mesh&, OM::EdgeHandle,     const StatusInfo&) = &set_status;
	void (*set_status_fh)(Mesh&, OM::FaceHandle,     const StatusInfo&) = &set_status;

	// Property management - add property
	void (Mesh::*add_property_vph)(OM::VPropHandleT<py::none>&, const std::string&) = &Mesh::add_property;
	void (Mesh::*add_property_eph)(OM::EPropHandleT<py::none>&, const std::string&) = &Mesh::add_property;
	void (Mesh::*add_property_hph)(OM::HPropHandleT<py::none>&, const std::string&) = &Mesh::add_property;
	void (Mesh::*add_property_fph)(OM::FPropHandleT<py::none>&, const std::string&) = &Mesh::add_property;
	void (Mesh::*add_property_mph)(OM::MPropHandleT<py::none>&, const std::string&) = &Mesh::add_property;

	// Property management - remove property
	void (Mesh::*remove_property_vph)(OM::VPropHandleT<py::none>&) = &Mesh::remove_property;
	void (Mesh::*remove_property_eph)(OM::EPropHandleT<py::none>&) = &Mesh::remove_property;
	void (Mesh::*remove_property_hph)(OM::HPropHandleT<py::none>&) = &Mesh::remove_property;
	void (Mesh::*remove_property_fph)(OM::FPropHandleT<py::none>&) = &Mesh::remove_property;
	void (Mesh::*remove_property_mph)(OM::MPropHandleT<py::none>&) = &Mesh::remove_property;

	// Property management - get property by name
	bool (Mesh::*get_property_handle_vph)(OM::VPropHandleT<py::none>&, const std::string&) const = &Mesh::get_property_handle;
	bool (Mesh::*get_property_handle_eph)(OM::EPropHandleT<py::none>&, const std::string&) const = &Mesh::get_property_handle;
	bool (Mesh::*get_property_handle_hph)(OM::HPropHandleT<py::none>&, const std::string&) const = &Mesh::get_property_handle;
	bool (Mesh::*get_property_handle_fph)(OM::FPropHandleT<py::none>&, const std::string&) const = &Mesh::get_property_handle;
	bool (Mesh::*get_property_handle_mph)(OM::MPropHandleT<py::none>&, const std::string&) const = &Mesh::get_property_handle;

	// Property management - get property value for an item
	const py::none& (Mesh::*property_vertex  )(OM::VPropHandleT<py::none>, OM::VertexHandle  ) const = &Mesh::property;
	const py::none& (Mesh::*property_edge    )(OM::EPropHandleT<py::none>, OM::EdgeHandle    ) const = &Mesh::property;
	const py::none& (Mesh::*property_halfedge)(OM::HPropHandleT<py::none>, OM::HalfedgeHandle) const = &Mesh::property;
	const py::none& (Mesh::*property_face    )(OM::FPropHandleT<py::none>, OM::FaceHandle    ) const = &Mesh::property;
	const py::none& (Mesh::*property_mesh    )(OM::MPropHandleT<py::none>                    ) const = &Mesh::property;

	// Property management - set property value for an item
	void (*set_property_vertex  )(Mesh&, OM::VPropHandleT<py::none>, OM::VertexHandle,   const py::object&) = &set_property;
	void (*set_property_edge    )(Mesh&, OM::EPropHandleT<py::none>, OM::EdgeHandle,     const py::object&) = &set_property;
	void (*set_property_halfedge)(Mesh&, OM::HPropHandleT<py::none>, OM::HalfedgeHandle, const py::object&) = &set_property;
	void (*set_property_face    )(Mesh&, OM::FPropHandleT<py::none>, OM::FaceHandle,     const py::object&) = &set_property;
	void (*set_property_mesh    )(Mesh&, OM::MPropHandleT<py::none>,                     const py::object&) = &set_property;

	// Low-level adding new items
	OM::VertexHandle (Mesh::*new_vertex_void )(void                        ) = &Mesh::new_vertex;
	OM::VertexHandle (Mesh::*new_vertex_point)(const typename Mesh::Point& ) = &Mesh::new_vertex;
	OM::FaceHandle   (Mesh::*new_face_void   )(void                        ) = &Mesh::new_face;
	OM::FaceHandle   (Mesh::*new_face_face   )(const typename Mesh::Face&  ) = &Mesh::new_face;

	// Kernel item iterators
	IteratorWrapperT<typename Mesh::VertexIter,   &Mesh::n_vertices > (*vertices )(Mesh&) = &get_iterator;
	IteratorWrapperT<typename Mesh::HalfedgeIter, &Mesh::n_halfedges> (*halfedges)(Mesh&) = &get_iterator;
	IteratorWrapperT<typename Mesh::EdgeIter,     &Mesh::n_edges    > (*edges    )(Mesh&) = &get_iterator;
	IteratorWrapperT<typename Mesh::FaceIter,     &Mesh::n_faces    > (*faces    )(Mesh&) = &get_iterator;

	IteratorWrapperT<typename Mesh::VertexIter,   &Mesh::n_vertices > (*svertices )(Mesh&) = &get_skipping_iterator;
	IteratorWrapperT<typename Mesh::HalfedgeIter, &Mesh::n_halfedges> (*shalfedges)(Mesh&) = &get_skipping_iterator;
	IteratorWrapperT<typename Mesh::EdgeIter,     &Mesh::n_edges    > (*sedges    )(Mesh&) = &get_skipping_iterator;
	IteratorWrapperT<typename Mesh::FaceIter,     &Mesh::n_faces    > (*sfaces    )(Mesh&) = &get_skipping_iterator;

	//======================================================================
	//  BaseKernel Function Pointers
	//======================================================================

	// Copy property
	void (Mesh::*copy_property_vprop)(OM::VPropHandleT<py::none>&, OM::VertexHandle,   OM::VertexHandle  ) = &Mesh::copy_property;
	void (Mesh::*copy_property_hprop)(OM::HPropHandleT<py::none>,  OM::HalfedgeHandle, OM::HalfedgeHandle) = &Mesh::copy_property;
	void (Mesh::*copy_property_eprop)(OM::EPropHandleT<py::none>,  OM::EdgeHandle,     OM::EdgeHandle    ) = &Mesh::copy_property;
	void (Mesh::*copy_property_fprop)(OM::FPropHandleT<py::none>,  OM::FaceHandle,     OM::FaceHandle    ) = &Mesh::copy_property;

	// Copy all properties
	void (Mesh::*copy_all_properties_vh_vh_bool)(OM::VertexHandle,   OM::VertexHandle,   bool) = &Mesh::copy_all_properties;
	void (Mesh::*copy_all_properties_hh_hh_bool)(OM::HalfedgeHandle, OM::HalfedgeHandle, bool) = &Mesh::copy_all_properties;
	void (Mesh::*copy_all_properties_eh_eh_bool)(OM::EdgeHandle,     OM::EdgeHandle,     bool) = &Mesh::copy_all_properties;
	void (Mesh::*copy_all_properties_fh_fh_bool)(OM::FaceHandle,     OM::FaceHandle,     bool) = &Mesh::copy_all_properties;

	//======================================================================
	//  PolyConnectivity Function Pointers
	//======================================================================

	// Assign connectivity
	void (*assign_connectivity_poly)(Mesh&, const PolyMesh&) = &assign_connectivity;
	void (*assign_connectivity_tri )(Mesh&, const TriMesh& ) = &assign_connectivity;

	// Vertex and face valence
	unsigned int (Mesh::*valence_vh)(OM::VertexHandle) const = &Mesh::valence;
	unsigned int (Mesh::*valence_fh)(OM::FaceHandle  ) const = &Mesh::valence;

	// Triangulate face or mesh
	void (Mesh::*triangulate_fh  )(OM::FaceHandle) = &Mesh::triangulate;
	void (Mesh::*triangulate_void)(              ) = &Mesh::triangulate;

	// Vertex and Face circulators
	CirculatorWrapperT<typename Mesh::VertexVertexIter,    OM::VertexHandle  > (*vv )(Mesh&, OM::VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::VertexIHalfedgeIter, OM::VertexHandle  > (*vih)(Mesh&, OM::VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::VertexOHalfedgeIter, OM::VertexHandle  > (*voh)(Mesh&, OM::VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::VertexEdgeIter,      OM::VertexHandle  > (*ve )(Mesh&, OM::VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::VertexFaceIter,      OM::VertexHandle  > (*vf )(Mesh&, OM::VertexHandle  ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::FaceVertexIter,      OM::FaceHandle    > (*fv )(Mesh&, OM::FaceHandle    ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::FaceHalfedgeIter,    OM::FaceHandle    > (*fh )(Mesh&, OM::FaceHandle    ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::FaceEdgeIter,        OM::FaceHandle    > (*fe )(Mesh&, OM::FaceHandle    ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::FaceFaceIter,        OM::FaceHandle    > (*ff )(Mesh&, OM::FaceHandle    ) = &get_circulator;
	CirculatorWrapperT<typename Mesh::HalfedgeLoopIter,    OM::HalfedgeHandle> (*hl )(Mesh&, OM::HalfedgeHandle) = &get_circulator;

	// Boundary and manifold tests
	bool (Mesh::*is_boundary_hh)(OM::HalfedgeHandle  ) const = &Mesh::is_boundary;
	bool (Mesh::*is_boundary_eh)(OM::EdgeHandle      ) const = &Mesh::is_boundary;
	bool (Mesh::*is_boundary_vh)(OM::VertexHandle    ) const = &Mesh::is_boundary;
	bool (Mesh::*is_boundary_fh)(OM::FaceHandle, bool) const = &Mesh::is_boundary;

	// Generic handle derefertiation
	const typename Mesh::Vertex&   (Mesh::*deref_vh)(OM::VertexHandle  ) const = &Mesh::deref;
	const typename Mesh::Halfedge& (Mesh::*deref_hh)(OM::HalfedgeHandle) const = &Mesh::deref;
	const typename Mesh::Edge&     (Mesh::*deref_eh)(OM::EdgeHandle    ) const = &Mesh::deref;
	const typename Mesh::Face&     (Mesh::*deref_fh)(OM::FaceHandle    ) const = &Mesh::deref;

	//======================================================================
	//  PolyMeshT Function Pointers
	//======================================================================

	void (Mesh::*calc_edge_vector_eh_normal)(OM::EdgeHandle,     Normal&) const = &Mesh::calc_edge_vector;
	void (Mesh::*calc_edge_vector_hh_normal)(OM::HalfedgeHandle, Normal&) const = &Mesh::calc_edge_vector;

	Normal (Mesh::*calc_edge_vector_eh)(OM::EdgeHandle    ) const = &Mesh::calc_edge_vector;
	Normal (Mesh::*calc_edge_vector_hh)(OM::HalfedgeHandle) const = &Mesh::calc_edge_vector;

	Scalar (Mesh::*calc_edge_length_eh)(OM::EdgeHandle    ) const = &Mesh::calc_edge_length;
	Scalar (Mesh::*calc_edge_length_hh)(OM::HalfedgeHandle) const = &Mesh::calc_edge_length;

	Scalar (Mesh::*calc_edge_sqr_length_eh)(OM::EdgeHandle    ) const = &Mesh::calc_edge_sqr_length;
	Scalar (Mesh::*calc_edge_sqr_length_hh)(OM::HalfedgeHandle) const = &Mesh::calc_edge_sqr_length;

	Scalar (Mesh::*calc_dihedral_angle_fast_hh)(OM::HalfedgeHandle) const = &Mesh::calc_dihedral_angle_fast;
	Scalar (Mesh::*calc_dihedral_angle_fast_eh)(OM::EdgeHandle    ) const = &Mesh::calc_dihedral_angle_fast;

	Scalar (Mesh::*calc_dihedral_angle_hh)(OM::HalfedgeHandle) const = &Mesh::calc_dihedral_angle;
	Scalar (Mesh::*calc_dihedral_angle_eh)(OM::EdgeHandle    ) const = &Mesh::calc_dihedral_angle;

	unsigned int (Mesh::*find_feature_edges)(Scalar) = &Mesh::find_feature_edges;

	void (Mesh::*split_fh_vh)(OM::FaceHandle, OM::VertexHandle) = &Mesh::split;
	void (Mesh::*split_eh_vh)(OM::EdgeHandle, OM::VertexHandle) = &Mesh::split;

	Normal (Mesh::*calc_face_normal    )(OM::FaceHandle            ) const = &Mesh::calc_face_normal;

	void  (Mesh::*calc_face_centroid_fh_point)(OM::FaceHandle, Point&) const = &Mesh::calc_face_centroid;
	Point (Mesh::*calc_face_centroid_fh      )(OM::FaceHandle        ) const = &Mesh::calc_face_centroid;

	//======================================================================
	//  Mesh Type
	//======================================================================

	py::class_<Mesh> class_mesh(m, _name);

	class_mesh
		.def(py::init<>())

		//======================================================================
		//  KernelT
		//======================================================================

		.def("reserve", &Mesh::reserve)

		.def("vertex", vertex, py::return_value_policy::reference)
		.def("halfedge", halfedge, py::return_value_policy::reference)
		.def("edge", edge, py::return_value_policy::reference)
		.def("face", face, py::return_value_policy::reference)

		.def("handle", handle_v)
		.def("handle", handle_h)
		.def("handle", handle_e)
		.def("handle", handle_f)

		.def("vertex_handle", vertex_handle_uint)
		.def("halfedge_handle", halfedge_handle_uint)
		.def("edge_handle", edge_handle_uint)
		.def("face_handle", face_handle_uint)

		.def("clear", &Mesh::clear)
		.def("clean", &Mesh::clean)
		.def("garbage_collection", garbage_collection_bools,
			py::arg("v")=true, py::arg("e")=true, py::arg("f")=true)
		.def("garbage_collection", garbage_collection_lists_bools,
			py::arg("vh_to_update"), py::arg("hh_to_update"), py::arg("fh_to_update"),
			py::arg("v")=true, py::arg("e")=true, py::arg("f")=true)

		.def("n_vertices", &Mesh::n_vertices)
		.def("n_halfedges", &Mesh::n_halfedges)
		.def("n_edges", &Mesh::n_edges)
		.def("n_faces", &Mesh::n_faces)
		.def("vertices_empty", &Mesh::vertices_empty)
		.def("halfedges_empty", &Mesh::halfedges_empty)
		.def("edges_empty", &Mesh::edges_empty)
		.def("faces_empty", &Mesh::faces_empty)

		.def("halfedge_handle", halfedge_handle_vh)
		.def("set_halfedge_handle", set_halfedge_handle_vh_hh)

		.def("to_vertex_handle", &Mesh::to_vertex_handle)
		.def("from_vertex_handle", &Mesh::from_vertex_handle)
		.def("set_vertex_handle", &Mesh::set_vertex_handle)
		.def("face_handle", face_handle_hh)
		.def("set_face_handle", &Mesh::set_face_handle)
		.def("next_halfedge_handle", &Mesh::next_halfedge_handle)
		.def("set_next_halfedge_handle", &Mesh::set_next_halfedge_handle)
		.def("prev_halfedge_handle", prev_halfedge_handle_hh)
		.def("opposite_halfedge_handle", &Mesh::opposite_halfedge_handle)
		.def("ccw_rotated_halfedge_handle", &Mesh::ccw_rotated_halfedge_handle)
		.def("cw_rotated_halfedge_handle", &Mesh::cw_rotated_halfedge_handle)
		.def("edge_handle", edge_handle_hh)

		.def("halfedge_handle", halfedge_handle_eh_uint)

		.def("halfedge_handle", halfedge_handle_fh)
		.def("set_halfedge_handle", set_halfedge_handle_fh_hh)

		.def("point_vec", point_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_point_vec", &Mesh::set_point)
		.def("normal_vec", normal_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_normal_vec", set_normal_vh)
		.def("normal_vec", normal_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_normal_vec", set_normal_hh)
		.def("color_vec", color_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_color_vec", set_color_vh)
		.def("texcoord1D_vec", texcoord1D_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord1D_vec", set_texcoord1D_vh)
		.def("texcoord2D_vec", texcoord2D_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord2D_vec", set_texcoord2D_vh)
		.def("texcoord3D_vec", texcoord3D_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord3D_vec", set_texcoord3D_vh)
		.def("texcoord1D_vec", texcoord1D_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord1D_vec", set_texcoord1D_hh)
		.def("texcoord2D_vec", texcoord2D_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord2D_vec", set_texcoord2D_hh)
		.def("texcoord3D_vec", texcoord3D_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_texcoord3D_vec", set_texcoord3D_hh)
		.def("status", status_vh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_status", set_status_vh)
		.def("status", status_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_status", set_status_hh)
		.def("color_vec", color_hh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_color_vec", set_color_hh)
		.def("color_vec", color_eh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_color_vec", set_color_eh)
		.def("status", status_eh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_status", set_status_eh)
		.def("normal_vec", normal_fh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_normal_vec", set_normal_fh)
		.def("color_vec", color_fh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_color_vec", set_color_fh)
		.def("status", status_fh, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("set_status", set_status_fh)

		.def("request_vertex_normals", &Mesh::request_vertex_normals)
		.def("request_vertex_colors", &Mesh::request_vertex_colors)
		.def("request_vertex_texcoords1D", &Mesh::request_vertex_texcoords1D)
		.def("request_vertex_texcoords2D", &Mesh::request_vertex_texcoords2D)
		.def("request_vertex_texcoords3D", &Mesh::request_vertex_texcoords3D)
		.def("request_vertex_status", &Mesh::request_vertex_status)
		.def("request_halfedge_status", &Mesh::request_halfedge_status)
		.def("request_halfedge_normals", &Mesh::request_halfedge_normals)
		.def("request_halfedge_colors", &Mesh::request_halfedge_colors)
		.def("request_halfedge_texcoords1D", &Mesh::request_halfedge_texcoords1D)
		.def("request_halfedge_texcoords2D", &Mesh::request_halfedge_texcoords2D)
		.def("request_halfedge_texcoords3D", &Mesh::request_halfedge_texcoords3D)
		.def("request_edge_status", &Mesh::request_edge_status)
		.def("request_edge_colors", &Mesh::request_edge_colors)
		.def("request_face_normals", &Mesh::request_face_normals)
		.def("request_face_colors", &Mesh::request_face_colors)
		.def("request_face_status", &Mesh::request_face_status)
		.def("request_face_texture_index", &Mesh::request_face_texture_index)

		.def("release_vertex_normals", &Mesh::release_vertex_normals)
		.def("release_vertex_colors", &Mesh::release_vertex_colors)
		.def("release_vertex_texcoords1D", &Mesh::release_vertex_texcoords1D)
		.def("release_vertex_texcoords2D", &Mesh::release_vertex_texcoords2D)
		.def("release_vertex_texcoords3D", &Mesh::release_vertex_texcoords3D)
		.def("release_vertex_status", &Mesh::release_vertex_status)
		.def("release_halfedge_status", &Mesh::release_halfedge_status)
		.def("release_halfedge_normals", &Mesh::release_halfedge_normals)
		.def("release_halfedge_colors", &Mesh::release_halfedge_colors)
		.def("release_halfedge_texcoords1D", &Mesh::release_halfedge_texcoords1D)
		.def("release_halfedge_texcoords2D", &Mesh::release_halfedge_texcoords2D)
		.def("release_halfedge_texcoords3D", &Mesh::release_halfedge_texcoords3D)
		.def("release_edge_status", &Mesh::release_edge_status)
		.def("release_edge_colors", &Mesh::release_edge_colors)
		.def("release_face_normals", &Mesh::release_face_normals)
		.def("release_face_colors", &Mesh::release_face_colors)
		.def("release_face_status", &Mesh::release_face_status)
		.def("release_face_texture_index", &Mesh::release_face_texture_index)

		.def("has_vertex_normals", &Mesh::has_vertex_normals)
		.def("has_vertex_colors", &Mesh::has_vertex_colors)
		.def("has_vertex_texcoords1D", &Mesh::has_vertex_texcoords1D)
		.def("has_vertex_texcoords2D", &Mesh::has_vertex_texcoords2D)
		.def("has_vertex_texcoords3D", &Mesh::has_vertex_texcoords3D)
		.def("has_vertex_status", &Mesh::has_vertex_status)
		.def("has_halfedge_status", &Mesh::has_halfedge_status)
		.def("has_halfedge_normals", &Mesh::has_halfedge_normals)
		.def("has_halfedge_colors", &Mesh::has_halfedge_colors)
		.def("has_halfedge_texcoords1D", &Mesh::has_halfedge_texcoords1D)
		.def("has_halfedge_texcoords2D", &Mesh::has_halfedge_texcoords2D)
		.def("has_halfedge_texcoords3D", &Mesh::has_halfedge_texcoords3D)
		.def("has_edge_status", &Mesh::has_edge_status)
		.def("has_edge_colors", &Mesh::has_edge_colors)
		.def("has_face_normals", &Mesh::has_face_normals)
		.def("has_face_colors", &Mesh::has_face_colors)
		.def("has_face_status", &Mesh::has_face_status)
		.def("has_face_texture_index", &Mesh::has_face_texture_index)

		.def("add_property", add_property_vph, py::arg("ph"), py::arg("name")="<vprop>")
		.def("add_property", add_property_eph, py::arg("ph"), py::arg("name")="<eprop>")
		.def("add_property", add_property_hph, py::arg("ph"), py::arg("name")="<hprop>")
		.def("add_property", add_property_fph, py::arg("ph"), py::arg("name")="<fprop>")
		.def("add_property", add_property_mph, py::arg("ph"), py::arg("name")="<mprop>")

		.def("remove_property", remove_property_vph)
		.def("remove_property", remove_property_eph)
		.def("remove_property", remove_property_hph)
		.def("remove_property", remove_property_fph)
		.def("remove_property", remove_property_mph)

		.def("get_property_handle", get_property_handle_vph)
		.def("get_property_handle", get_property_handle_eph)
		.def("get_property_handle", get_property_handle_hph)
		.def("get_property_handle", get_property_handle_fph)
		.def("get_property_handle", get_property_handle_mph)

		.def("property", property_vertex, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("property", property_edge, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("property", property_halfedge, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("property", property_face, OPENMESH_PYTHON_DEFAULT_POLICY)
		.def("property", property_mesh, OPENMESH_PYTHON_DEFAULT_POLICY)

		.def("set_property", set_property_vertex)
		.def("set_property", set_property_edge)
		.def("set_property", set_property_halfedge)
		.def("set_property", set_property_face)
		.def("set_property", set_property_mesh)

		.def("new_vertex", new_vertex_void)
		.def("new_vertex", new_vertex_point)
		.def("new_edge", &Mesh::new_edge)
		.def("new_face", new_face_void)
		.def("new_face", new_face_face)

		.def("vertices", vertices)
		.def("halfedges", halfedges)
		.def("edges", edges)
		.def("faces", faces)

		.def("svertices", svertices)
		.def("shalfedges", shalfedges)
		.def("sedges", sedges)
		.def("sfaces", sfaces)

		//======================================================================
		//  BaseKernel
		//======================================================================

		.def("copy_property", copy_property_vprop)
		.def("copy_property", copy_property_hprop)
		.def("copy_property", copy_property_eprop)
		.def("copy_property", copy_property_fprop)

		.def("copy_all_properties", copy_all_properties_vh_vh_bool,
			py::arg("vh_from"), py::arg("vh_to"), py::arg("copy_build_in")=false)
		.def("copy_all_properties", copy_all_properties_hh_hh_bool,
			py::arg("hh_from"), py::arg("hh_to"), py::arg("copy_build_in")=false)
		.def("copy_all_properties", copy_all_properties_eh_eh_bool,
			py::arg("eh_from"), py::arg("eh_to"), py::arg("copy_build_in")=false)
		.def("copy_all_properties", copy_all_properties_fh_fh_bool,
			py::arg("fh_from"), py::arg("fh_to"), py::arg("copy_build_in")=false)

		//======================================================================
		//  PolyConnectivity
		//======================================================================

		.def("assign_connectivity", assign_connectivity_poly)
		.def("assign_connectivity", assign_connectivity_tri)

		.def("opposite_face_handle", &Mesh::opposite_face_handle)
		.def("adjust_outgoing_halfedge", &Mesh::adjust_outgoing_halfedge)
		.def("find_halfedge", &Mesh::find_halfedge)
		.def("valence", valence_vh)
		.def("valence", valence_fh)
		.def("collapse", &Mesh::collapse)
		.def("is_simple_link", &Mesh::is_simple_link)
		.def("is_simply_connected", &Mesh::is_simply_connected)
		.def("remove_edge", &Mesh::remove_edge)
		.def("reinsert_edge", &Mesh::reinsert_edge)
		.def("triangulate", triangulate_fh)
		.def("triangulate", triangulate_void)
		.def("split_edge", &Mesh::split_edge)
		.def("split_edge_copy", &Mesh::split_edge_copy)

		.def("add_vertex", &Mesh::add_vertex)
		.def("add_vertex", [](Mesh& _self, py::array_t<typename Point::value_type> _arr)
			{ return _self.add_vertex(Point(_arr.at(0), _arr.at(1), _arr.at(2))); })

		.def("is_collapse_ok",  &Mesh::is_collapse_ok)

		.def("delete_vertex", [](Mesh& _self, OM::VertexHandle _vh, bool _delete_isolated) {
				if (!_self.has_vertex_status()) _self.request_vertex_status();
				if (!_self.has_face_status()) _self.request_face_status();
				_self.delete_vertex(_vh, _delete_isolated);
			}, py::arg("vh"), py::arg("delete_isolated_vertices")=true)
		.def("delete_edge", [](Mesh& _self, OM::EdgeHandle _eh, bool _delete_isolated) {
				if (!_self.has_vertex_status() && _delete_isolated) _self.request_vertex_status();
				if (!_self.has_edge_status()) _self.request_edge_status();
				if (!_self.has_face_status()) _self.request_face_status();
				_self.delete_edge(_eh, _delete_isolated);
			}, py::arg("eh"), py::arg("delete_isolated_vertices")=true)
		.def("delete_face", [](Mesh& _self, OM::FaceHandle _fh, bool _delete_isolated) {
				if (!_self.has_vertex_status() && _delete_isolated) _self.request_vertex_status();
				if (!_self.has_face_status()) _self.request_face_status();
				_self.delete_face(_fh, _delete_isolated);
			}, py::arg("fh"), py::arg("delete_isolated_vertices")=true)

		.def("vv", vv)
		.def("vih", vih)
		.def("voh", voh)
		.def("ve", ve)
		.def("vf", vf)

		.def("fv", fv)
		.def("fh", fh)
		.def("fe", fe)
		.def("ff", ff)

		.def("hl", hl)

		.def("is_boundary", is_boundary_hh)
		.def("is_boundary", is_boundary_eh)
		.def("is_boundary", is_boundary_vh)
		.def("is_boundary", is_boundary_fh, py::arg("fh"), py::arg("check_vertex")=false)
		.def("is_manifold", &Mesh::is_manifold)

		.def("deref", deref_vh, py::return_value_policy::reference)
		.def("deref", deref_hh, py::return_value_policy::reference)
		.def("deref", deref_eh, py::return_value_policy::reference)
		.def("deref", deref_fh, py::return_value_policy::reference)

		.def_static("is_triangles", &Mesh::is_triangles)

		.def_readonly_static("InvalidVertexHandle", &Mesh::InvalidVertexHandle)
		.def_readonly_static("InvalidHalfedgeHandle", &Mesh::InvalidHalfedgeHandle)
		.def_readonly_static("InvalidEdgeHandle", &Mesh::InvalidEdgeHandle)
		.def_readonly_static("InvalidFaceHandle", &Mesh::InvalidFaceHandle)

		//======================================================================
		//  PolyMeshT
		//======================================================================

		.def("add_vertex", &Mesh::add_vertex)

		.def("calc_edge_vector", calc_edge_vector_eh_normal)
		.def("calc_edge_vector", calc_edge_vector_eh)
		.def("calc_edge_vector", calc_edge_vector_hh_normal)
		.def("calc_edge_vector", calc_edge_vector_hh)

		.def("calc_edge_length", calc_edge_length_eh)
		.def("calc_edge_length", calc_edge_length_hh)
		.def("calc_edge_sqr_length", calc_edge_sqr_length_eh)
		.def("calc_edge_sqr_length", calc_edge_sqr_length_hh)

		.def("calc_sector_vectors", &Mesh::calc_sector_vectors)
		.def("calc_sector_angle", &Mesh::calc_sector_angle)
		.def("calc_sector_normal", &Mesh::calc_sector_normal)
		.def("calc_sector_area", &Mesh::calc_sector_area)

		.def("calc_dihedral_angle_fast", calc_dihedral_angle_fast_hh)
		.def("calc_dihedral_angle_fast", calc_dihedral_angle_fast_eh)
		.def("calc_dihedral_angle", calc_dihedral_angle_hh)
		.def("calc_dihedral_angle", calc_dihedral_angle_eh)

		.def("find_feature_edges", find_feature_edges, py::arg("angle_tresh")=OM::deg_to_rad(44.0))

		.def("split", split_fh_vh)
		.def("split", split_eh_vh)

		.def("update_normals", [](Mesh& _self) {
				if (!_self.has_face_normals()) {
					_self.request_face_normals();
				}
				if (!_self.has_halfedge_normals()) {
					_self.request_halfedge_normals();
				}
				if (!_self.has_vertex_normals()) {
					_self.request_vertex_normals();
				}
				_self.update_normals();
			})

		.def("update_normal", [](Mesh& _self, OM::FaceHandle _fh) {
				if (!_self.has_face_normals()) _self.request_face_normals();
				_self.update_normal(_fh);
			})
		.def("update_face_normals", [](Mesh& _self) {
				if (!_self.has_face_normals()) _self.request_face_normals();
				_self.update_face_normals();
			})

		.def("update_normal", [](Mesh& _self, OM::HalfedgeHandle _hh, double _feature_angle) {
				if (!_self.has_face_normals()) {
					_self.request_face_normals();
					_self.update_face_normals();
				}
				if (!_self.has_halfedge_normals()) {
					_self.request_halfedge_normals();
				}
				_self.update_normal(_hh, _feature_angle);
			}, py::arg("heh"), py::arg("feature_angle")=0.8)
		.def("update_halfedge_normals", [](Mesh& _self, double _feature_angle) {
				if (!_self.has_face_normals()) {
					_self.request_face_normals();
					_self.update_face_normals();
				}
				if (!_self.has_halfedge_normals()) {
					_self.request_halfedge_normals();
				}
				_self.update_halfedge_normals(_feature_angle);
			}, py::arg("feature_angle")=0.8)

		.def("update_normal", [](Mesh& _self, OM::VertexHandle _vh) {
				if (!_self.has_face_normals()) {
					_self.request_face_normals();
					_self.update_face_normals();
				}
				if (!_self.has_vertex_normals()) {
					_self.request_vertex_normals();
				}
				_self.update_normal(_vh);
			})
		.def("update_vertex_normals", [](Mesh& _self) {
				if (!_self.has_face_normals()) {
					_self.request_face_normals();
					_self.update_face_normals();
				}
				if (!_self.has_vertex_normals()) {
					_self.request_vertex_normals();
				}
				_self.update_vertex_normals();
			})

		.def("calc_face_normal", calc_face_normal)
		.def("calc_halfedge_normal", &Mesh::calc_halfedge_normal, py::arg("heh"), py::arg("feature_angle")=0.8)

		.def("calc_face_centroid", calc_face_centroid_fh_point)
		.def("calc_face_centroid", calc_face_centroid_fh)

		.def("is_estimated_feature_edge", &Mesh::is_estimated_feature_edge)

		.def("calc_vertex_normal", &Mesh::calc_vertex_normal)
		.def("calc_vertex_normal_fast", &Mesh::calc_vertex_normal_fast)
		.def("calc_vertex_normal_correct", &Mesh::calc_vertex_normal_correct)
		.def("calc_vertex_normal_loop", &Mesh::calc_vertex_normal_loop)

		.def_static("is_polymesh", &Mesh::is_polymesh)
		.def("is_trimesh", &Mesh::is_trimesh)

		//======================================================================
		//  numpy vector getter
		//======================================================================

		.def("point", [](Mesh& _self, OM::VertexHandle _h) {
				return vec2numpy(_self, _self.point(_h));
			})

		.def("normal", [](Mesh& _self, OM::VertexHandle _h) {
				if (!_self.has_vertex_normals()) _self.request_vertex_normals();
				return vec2numpy(_self, _self.normal(_h));
			})
		.def("normal", [](Mesh& _self, OM::HalfedgeHandle _h) {
				if (!_self.has_halfedge_normals()) _self.request_halfedge_normals();
				return vec2numpy(_self, _self.normal(_h));
			})
		.def("normal", [](Mesh& _self, OM::FaceHandle _h) {
				if (!_self.has_face_normals()) _self.request_face_normals();
				return vec2numpy(_self, _self.normal(_h));
			})

		.def("color", [](Mesh& _self, OM::VertexHandle _h) {
				if (!_self.has_vertex_colors()) _self.request_vertex_colors();
				return vec2numpy(_self, _self.color(_h));
			})
		.def("color", [](Mesh& _self, OM::HalfedgeHandle _h) {
				if (!_self.has_halfedge_colors()) _self.request_halfedge_colors();
				return vec2numpy(_self, _self.color(_h));
			})
		.def("color", [](Mesh& _self, OM::EdgeHandle _h) {
				if (!_self.has_edge_colors()) _self.request_edge_colors();
				return vec2numpy(_self, _self.color(_h));
			})
		.def("color", [](Mesh& _self, OM::FaceHandle _h) {
				if (!_self.has_face_colors()) _self.request_face_colors();
				return vec2numpy(_self, _self.color(_h));
			})

		.def("texcoord1D", [](Mesh& _self, OM::VertexHandle _h) {
				if (!_self.has_vertex_texcoords1D()) _self.request_vertex_texcoords1D();
				return flt2numpy(_self, _self.texcoord1D(_h));
			})
		.def("texcoord1D", [](Mesh& _self, OM::HalfedgeHandle _h) {
				if (!_self.has_halfedge_texcoords1D()) _self.request_halfedge_texcoords1D();
				return flt2numpy(_self, _self.texcoord1D(_h));
			})
		.def("texcoord2D", [](Mesh& _self, OM::VertexHandle _h) {
				if (!_self.has_vertex_texcoords2D()) _self.request_vertex_texcoords2D();
				return vec2numpy(_self, _self.texcoord2D(_h));
			})
		.def("texcoord2D", [](Mesh& _self, OM::HalfedgeHandle _h) {
				if (!_self.has_halfedge_texcoords2D()) _self.request_halfedge_texcoords2D();
				return vec2numpy(_self, _self.texcoord2D(_h));
			})
		.def("texcoord3D", [](Mesh& _self, OM::VertexHandle _h) {
				if (!_self.has_vertex_texcoords3D()) _self.request_vertex_texcoords3D();
				return vec2numpy(_self, _self.texcoord3D(_h));
			})
		.def("texcoord3D", [](Mesh& _self, OM::HalfedgeHandle _h) {
				if (!_self.has_halfedge_texcoords3D()) _self.request_halfedge_texcoords3D();
				return vec2numpy(_self, _self.texcoord3D(_h));
			})

		//======================================================================
		//  numpy vector setter
		//======================================================================

		.def("set_point", [](Mesh& _self, OM::VertexHandle _h, py::array_t<typename Point::value_type> _arr) {
				_self.point(_h) = Point(_arr.at(0), _arr.at(1), _arr.at(2));
			})

		.def("set_normal", [](Mesh& _self, OM::VertexHandle _h, py::array_t<typename Normal::value_type> _arr) {
				if (!_self.has_vertex_normals()) _self.request_vertex_normals();
				_self.set_normal(_h, typename Mesh::Normal(_arr.at(0), _arr.at(1), _arr.at(2)));
			})
		.def("set_normal", [](Mesh& _self, OM::HalfedgeHandle _h, py::array_t<typename Normal::value_type> _arr) {
				if (!_self.has_halfedge_normals()) _self.request_halfedge_normals();
				_self.set_normal(_h, typename Mesh::Normal(_arr.at(0), _arr.at(1), _arr.at(2)));
			})
		.def("set_normal", [](Mesh& _self, OM::FaceHandle _h, py::array_t<typename Normal::value_type> _arr) {
				if (!_self.has_face_normals()) _self.request_face_normals();
				_self.set_normal(_h, typename Mesh::Normal(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		.def("set_color", [](Mesh& _self, OM::VertexHandle _h, py::array_t<typename Color::value_type> _arr) {
				if(!_self.has_vertex_colors()) _self.request_vertex_colors();
				_self.set_color(_h, typename Mesh::Color(_arr.at(0), _arr.at(1), _arr.at(2), _arr.at(3)));
			})
		.def("set_color", [](Mesh& _self, OM::HalfedgeHandle _h, py::array_t<typename Color::value_type> _arr) {
				if(!_self.has_halfedge_colors()) _self.request_halfedge_colors();
				_self.set_color(_h, typename Mesh::Color(_arr.at(0), _arr.at(1), _arr.at(2), _arr.at(3)));
			})
		.def("set_color", [](Mesh& _self, OM::EdgeHandle _h, py::array_t<typename Color::value_type> _arr) {
				if(!_self.has_edge_colors()) _self.request_edge_colors();
				_self.set_color(_h, typename Mesh::Color(_arr.at(0), _arr.at(1), _arr.at(2), _arr.at(3)));
			})
		.def("set_color", [](Mesh& _self, OM::FaceHandle _h, py::array_t<typename Color::value_type> _arr) {
				if(!_self.has_face_colors()) _self.request_face_colors();
				_self.set_color(_h, typename Mesh::Color(_arr.at(0), _arr.at(1), _arr.at(2), _arr.at(3)));
			})

		.def("set_texcoord1D", [](Mesh& _self, OM::VertexHandle _h, py::array_t<TexCoord1D> _arr) {
				if (!_self.has_vertex_texcoords1D()) _self.request_vertex_texcoords1D();
				_self.set_texcoord1D(_h, _arr.at(0));
			})
		.def("set_texcoord1D", [](Mesh& _self, OM::HalfedgeHandle _h, py::array_t<TexCoord1D> _arr) {
				if (!_self.has_halfedge_texcoords1D()) _self.request_halfedge_texcoords1D();
				_self.set_texcoord1D(_h, _arr.at(0));
			})

		.def("set_texcoord2D", [](Mesh& _self, OM::VertexHandle _h, py::array_t<typename TexCoord2D::value_type> _arr) {
				if (!_self.has_vertex_texcoords2D()) _self.request_vertex_texcoords2D();
				_self.set_texcoord2D(_h, typename Mesh::TexCoord2D(_arr.at(0), _arr.at(1)));
			})
		.def("set_texcoord2D", [](Mesh& _self, OM::HalfedgeHandle _h, py::array_t<typename TexCoord2D::value_type> _arr) {
				if (!_self.has_halfedge_texcoords2D()) _self.request_halfedge_texcoords2D();
				_self.set_texcoord2D(_h, typename Mesh::TexCoord2D(_arr.at(0), _arr.at(1)));
			})

		.def("set_texcoord3D", [](Mesh& _self, OM::VertexHandle _h, py::array_t<typename TexCoord3D::value_type> _arr) {
				if (!_self.has_vertex_texcoords3D()) _self.request_vertex_texcoords3D();
				_self.set_texcoord3D(_h, typename Mesh::TexCoord3D(_arr.at(0), _arr.at(1), _arr.at(2)));
			})
		.def("set_texcoord3D", [](Mesh& _self, OM::HalfedgeHandle _h, py::array_t<typename TexCoord3D::value_type> _arr) {
				if (!_self.has_halfedge_texcoords3D()) _self.request_halfedge_texcoords3D();
				_self.set_texcoord3D(_h, typename Mesh::TexCoord3D(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		//======================================================================
		//  numpy matrix getter
		//======================================================================

		.def("points", [](Mesh& _self) {
				return vec2numpy(_self, _self.point(OM::VertexHandle(0)), _self.n_vertices());
			})

		.def("vertex_normals", [](Mesh& _self) {
				if (!_self.has_vertex_normals()) _self.request_vertex_normals();
				return vec2numpy(_self, _self.normal(OM::VertexHandle(0)), _self.n_vertices());
			})
		.def("vertex_colors", [](Mesh& _self) {
				if (!_self.has_vertex_colors()) _self.request_vertex_colors();
				return vec2numpy(_self, _self.color(OM::VertexHandle(0)), _self.n_vertices());
			})
		.def("vertex_texcoords1D", [](Mesh& _self) {
				if (!_self.has_vertex_texcoords1D()) _self.request_vertex_texcoords1D();
				return flt2numpy(_self, _self.texcoord1D(OM::VertexHandle(0)), _self.n_vertices());
			})
		.def("vertex_texcoords2D", [](Mesh& _self) {
				if (!_self.has_vertex_texcoords2D()) _self.request_vertex_texcoords2D();
				return vec2numpy(_self, _self.texcoord2D(OM::VertexHandle(0)), _self.n_vertices());
			})
		.def("vertex_texcoords3D", [](Mesh& _self) {
				if (!_self.has_vertex_texcoords3D()) _self.request_vertex_texcoords3D();
				return vec2numpy(_self, _self.texcoord3D(OM::VertexHandle(0)), _self.n_vertices());
			})

		.def("halfedge_normals", [](Mesh& _self) {
				if (!_self.has_halfedge_normals()) _self.request_halfedge_normals();
				return vec2numpy(_self, _self.normal(OM::HalfedgeHandle(0)), _self.n_halfedges());
			})
		.def("halfedge_colors", [](Mesh& _self) {
				if (!_self.has_halfedge_colors()) _self.request_halfedge_colors();
				return vec2numpy(_self, _self.color(OM::HalfedgeHandle(0)), _self.n_halfedges());
			})
		.def("halfedge_texcoords1D", [](Mesh& _self) {
				if (!_self.has_halfedge_texcoords1D()) _self.request_halfedge_texcoords1D();
				return flt2numpy(_self, _self.texcoord1D(OM::HalfedgeHandle(0)), _self.n_halfedges());
			})
		.def("halfedge_texcoords2D", [](Mesh& _self) {
				if (!_self.has_halfedge_texcoords2D()) _self.request_halfedge_texcoords2D();
				return vec2numpy(_self, _self.texcoord2D(OM::HalfedgeHandle(0)), _self.n_halfedges());
			})
		.def("halfedge_texcoords3D", [](Mesh& _self) {
				if (!_self.has_halfedge_texcoords3D()) _self.request_halfedge_texcoords3D();
				return vec2numpy(_self, _self.texcoord3D(OM::HalfedgeHandle(0)), _self.n_halfedges());
			})

		.def("edge_colors", [](Mesh& _self) {
				if (!_self.has_edge_colors()) _self.request_edge_colors();
				return vec2numpy(_self, _self.color(OM::EdgeHandle(0)), _self.n_edges());
			})

		.def("face_normals", [](Mesh& _self) {
				if (!_self.has_face_normals()) _self.request_face_normals();
				return vec2numpy(_self, _self.normal(OM::FaceHandle(0)), _self.n_faces());
			})
		.def("face_colors",  [](Mesh& _self) {
				if (!_self.has_face_colors()) _self.request_face_colors();
				return vec2numpy(_self, _self.color (OM::FaceHandle(0)), _self.n_faces());
			})

		//======================================================================
		//  property_array
		//======================================================================

		.def("property_array", [] (Mesh& _self, OM::VPropHandleT<py::none> _ph) {
				return property_array<Mesh, OM::VPropHandleT<py::none>, OM::VertexHandle>(_self, _ph, _self.n_vertices());
			})
		.def("property_array", [] (Mesh& _self, OM::HPropHandleT<py::none> _ph) {
				return property_array<Mesh, OM::HPropHandleT<py::none>, OM::HalfedgeHandle>(_self, _ph, _self.n_halfedges());
			})
		.def("property_array", [] (Mesh& _self, OM::EPropHandleT<py::none> _ph) {
				return property_array<Mesh, OM::EPropHandleT<py::none>, OM::EdgeHandle>(_self, _ph, _self.n_edges());
			})
		.def("property_array", [] (Mesh& _self, OM::FPropHandleT<py::none> _ph) {
				return property_array<Mesh, OM::FPropHandleT<py::none>, OM::FaceHandle>(_self, _ph, _self.n_faces());
			})

		//======================================================================
		//  set_property_array
		//======================================================================

		.def("set_property_array", [] (Mesh& _self, OM::VPropHandleT<py::none> _ph, py::array_t<double, py::array::c_style | py::array::forcecast> _arr) {
				return set_property_array<Mesh, OM::VPropHandleT<py::none>, OM::VertexHandle>(_self, _ph, _arr, _self.n_vertices());
			})
		.def("set_property_array", [] (Mesh& _self, OM::HPropHandleT<py::none> _ph, py::array_t<double, py::array::c_style | py::array::forcecast> _arr) {
				return set_property_array<Mesh, OM::HPropHandleT<py::none>, OM::HalfedgeHandle>(_self, _ph, _arr, _self.n_halfedges());
			})
		.def("set_property_array", [] (Mesh& _self, OM::EPropHandleT<py::none> _ph, py::array_t<double, py::array::c_style | py::array::forcecast> _arr) {
				return set_property_array<Mesh, OM::EPropHandleT<py::none>, OM::EdgeHandle>(_self, _ph, _arr, _self.n_edges());
			})
		.def("set_property_array", [] (Mesh& _self, OM::FPropHandleT<py::none> _ph, py::array_t<double, py::array::c_style | py::array::forcecast> _arr) {
				return set_property_array<Mesh, OM::FPropHandleT<py::none>, OM::FaceHandle>(_self, _ph, _arr, _self.n_faces());
			})
		;

	expose_type_specific_functions(class_mesh);

	//======================================================================
	//  Nested Types
	//======================================================================

	class_mesh.attr("Point") = m.attr("Vec3d");
	class_mesh.attr("Normal") = m.attr("Vec3d");
	class_mesh.attr("Color") = m.attr("Vec4f");
	class_mesh.attr("TexCoord2D") = m.attr("Vec2f");
	class_mesh.attr("TexCoord3D") = m.attr("Vec3f");

}

#endif

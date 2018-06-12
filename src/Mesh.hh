#ifndef OPENMESH_PYTHON_MESH_HH
#define OPENMESH_PYTHON_MESH_HH

#include "Utilities.hh"
#include "MeshTypes.hh"
#include "Iterator.hh"
#include "Circulator.hh"

#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace OM = OpenMesh;


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
 * Converts OpenMesh vectors to numpy arrays.
 *
 * @tparam vector A Vector type.
 * @param _vec The vector to be converted.
 */
template<class Vector>
py::array_t<typename Vector::value_type> vec2numpy(const Vector& _vec) {
	typedef typename Vector::value_type dtype;
	dtype *data = new dtype[_vec.size()];
	std::copy_n(_vec.data(), _vec.size(), data);
	py::capsule base = free_when_done(data);
	return py::array_t<dtype>({_vec.size()}, {sizeof(dtype)}, data, base);
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

py::array_t<int> face_vertex_indices_trimesh(TriMesh& _self) {
	if (_self.n_faces() == 0) {
		return py::array_t<int>();
	}

	const bool has_status = _self.has_face_status();

	int *indices = new int[_self.n_faces() * 3];
	py::capsule base = free_when_done(indices);

	for (auto fh : _self.all_faces()) {
		if (has_status && _self.status(fh).deleted()) {
			PyErr_SetString(PyExc_RuntimeError, "Mesh has deleted items. Please call garbage_collection() first.");
			throw py::error_already_set();
		}
		auto fv_it = _self.fv_iter(fh);
		indices[fh.idx() * 3 + 0] = fv_it->idx(); ++fv_it;
		indices[fh.idx() * 3 + 1] = fv_it->idx(); ++fv_it;
		indices[fh.idx() * 3 + 2] = fv_it->idx();
	}
	const auto shape = {_self.n_faces(), size_t(3)};
	const auto strides = {3 * sizeof(int), sizeof(int)};
	return py::array_t<int>(shape, strides, indices, base);
}

struct FuncEdgeVertex {
	static void call(const OM::ArrayKernel& _mesh, OM::EdgeHandle _eh, int *_ptr) {
		const auto heh = _mesh.halfedge_handle(_eh, 0);
		_ptr[0] = _mesh.from_vertex_handle(heh).idx();
		_ptr[1] = _mesh.to_vertex_handle(heh).idx();
	}
};

struct FuncEdgeFace {
	static void call(const OM::ArrayKernel& _mesh, OM::EdgeHandle _eh, int *_ptr) {
		const auto heh1 = _mesh.halfedge_handle(_eh, 0);
		const auto heh2 = _mesh.halfedge_handle(_eh, 1);
		_ptr[0] = _mesh.face_handle(heh1).idx();
		_ptr[1] = _mesh.face_handle(heh2).idx();
	}
};

struct FuncEdgeHalfedge {
	static void call(const OM::ArrayKernel& _mesh, OM::EdgeHandle _eh, int *_ptr) {
		_ptr[0] = _mesh.halfedge_handle(_eh, 0).idx();
		_ptr[1] = _mesh.halfedge_handle(_eh, 1).idx();
	}
};

struct FuncHalfedgeToVertex {
	static void call(const OM::ArrayKernel& _mesh, OM::HalfedgeHandle _heh, int *_ptr) {
		*_ptr = _mesh.to_vertex_handle(_heh).idx();
	}
	static size_t dim() { return 1; }
};

struct FuncHalfedgeFromVertex {
	static void call(const OM::ArrayKernel& _mesh, OM::HalfedgeHandle _heh, int *_ptr) {
		*_ptr = _mesh.from_vertex_handle(_heh).idx();
	}
	static size_t dim() { return 1; }
};

struct FuncHalfedgeFace {
	static void call(const OM::ArrayKernel& _mesh, OM::HalfedgeHandle _heh, int *_ptr) {
		*_ptr = _mesh.face_handle(_heh).idx();
	}
	static size_t dim() { return 1; }
};

struct FuncHalfedgeEdge {
	static void call(const OM::ArrayKernel& _mesh, OM::HalfedgeHandle _heh, int *_ptr) {
		*_ptr = _mesh.edge_handle(_heh).idx();
	}
	static size_t dim() { return 1; }
};

struct FuncHalfedgeVertex {
	static void call(const OM::ArrayKernel& _mesh, OM::HalfedgeHandle _heh, int *_ptr) {
		_ptr[0] = _mesh.from_vertex_handle(_heh).idx();
		_ptr[1] = _mesh.to_vertex_handle(_heh).idx();
	}
	static size_t dim() { return 2; }
};

template <class Mesh, class CopyFunc>
py::array_t<int> edge_other_indices(Mesh& _self) {
	if (_self.n_edges() == 0) {
		return py::array_t<int>();
	}

	const bool has_status = _self.has_edge_status();

	int *indices = new int[_self.n_edges() * 2];
	py::capsule base = free_when_done(indices);

	for (auto eh : _self.all_edges()) {
		if (has_status && _self.status(eh).deleted()) {
			PyErr_SetString(PyExc_RuntimeError, "Mesh has deleted items. Please call garbage_collection() first.");
			throw py::error_already_set();
		}
		CopyFunc::call(_self, eh, &indices[eh.idx() * 2]);
	}
	const auto shape = {_self.n_edges(), size_t(2)};
	const auto strides = {2 * sizeof(int), sizeof(int)};
	return py::array_t<int>(shape, strides, indices, base);
}

template <class Mesh, class CopyFunc>
py::array_t<int> halfedge_other_indices(Mesh& _self) {
	if (_self.n_halfedges() == 0) {
		return py::array_t<int>();
	}

	const bool has_status = _self.has_halfedge_status();
	const size_t dim = CopyFunc::dim();

	int *indices = new int[_self.n_halfedges() * dim];
	py::capsule base = free_when_done(indices);

	for (auto heh : _self.all_halfedges()) {
		if (has_status && _self.status(heh).deleted()) {
			PyErr_SetString(PyExc_RuntimeError, "Mesh has deleted items. Please call garbage_collection() first.");
			throw py::error_already_set();
		}
		CopyFunc::call(_self, heh, &indices[heh.idx() * dim]);
	}

	std::vector<size_t> shape;
	std::vector<size_t> strides;
	if (dim == 1) {
		shape = {_self.n_halfedges()};
		strides = {sizeof(int)};
	}
	else {
		shape = {_self.n_halfedges(), dim};
		strides = {dim * sizeof(int), sizeof(int)};
	}

	return py::array_t<int>(shape, strides, indices, base);
}

template <class Mesh, class Handle, class Circulator>
py::array_t<int> indices(Mesh& _self) {
	const size_t n = _self.py_n_items(Handle());
	if (n == 0) return py::array_t<int>();

	const bool has_status = _self.py_has_status(Handle());

	// find max valence and check status
	int max_valence = 0;
	for (size_t i = 0; i < n; ++i) {
		Handle hnd(i);
		if (has_status && _self.status(hnd).deleted()) {
			PyErr_SetString(PyExc_RuntimeError, "Mesh has deleted items. Please call garbage_collection() first.");
			throw py::error_already_set();
		}
		int valence = 0;
		for (auto it = Circulator(_self, hnd); it.is_valid(); ++it) {
			valence++;
		}
		max_valence = std::max(max_valence, valence);
	}

	// allocate memory
	int *indices = new int[n * max_valence];

	// copy indices
	for (size_t i = 0; i < n; ++i) {
		int valence = 0;
		for (auto it = Circulator(_self, Handle(i)); it.is_valid(); ++it) {
			indices[i * max_valence + valence] = it->idx();
			valence++;
		}
		for (size_t j = valence; j < max_valence; ++j) {
			indices[i * max_valence + j] = -1;
		}
	}

	// make numpy array
	const auto shape = {n, size_t(max_valence)};
	const auto strides = {max_valence * sizeof(int), sizeof(int)};
	py::capsule base = free_when_done(indices);
	return py::array_t<int>(shape, strides, indices, base);
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

	typedef py::array_t<typename Point::value_type> np_point_t;

	OM::FaceHandle (PolyMesh::*add_face_4_vh)(OM::VertexHandle, OM::VertexHandle, OM::VertexHandle, OM::VertexHandle) = &PolyMesh::add_face;

	_class
		.def("add_face", add_face_4_vh)

		.def("split", [](PolyMesh& _self, OM::EdgeHandle _eh, np_point_t _arr) {
				_self.split(_eh, Point(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		.def("split", [](PolyMesh& _self, OM::FaceHandle _fh, np_point_t _arr) {
				_self.split(_fh, Point(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		.def("insert_edge", &PolyMesh::insert_edge)

		.def("face_vertex_indices", &indices<PolyMesh, OM::FaceHandle, PolyMesh::FaceVertexIter>)
		.def("fv_indices", &indices<PolyMesh, OM::FaceHandle, PolyMesh::FaceVertexIter>)

		.def("calc_face_normal", [](PolyMesh& _self, np_point_t _p0, np_point_t _p1, np_point_t _p2) {
				const Point p0(_p0.at(0), _p0.at(1), _p0.at(2));
				const Point p1(_p1.at(0), _p1.at(1), _p1.at(2));
				const Point p2(_p2.at(0), _p2.at(1), _p2.at(2));
				return vec2numpy(_self.calc_face_normal(p0, p1, p2));
			})
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

	typedef py::array_t<typename Point::value_type> np_point_t;

	void (TriMesh::*split_copy_eh_vh)(OM::EdgeHandle, OM::VertexHandle) = &TriMesh::split_copy;
	OM::HalfedgeHandle (TriMesh::*vertex_split_vh)(OM::VertexHandle, OM::VertexHandle, OM::VertexHandle, OM::VertexHandle) = &TriMesh::vertex_split;

	_class
		.def("split", [](TriMesh& _self, OM::EdgeHandle _eh, np_point_t _arr) {
				return _self.split(_eh, Point(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		.def("split", [](TriMesh& _self, OM::FaceHandle _fh, np_point_t _arr) {
				return _self.split(_fh, Point(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		.def("split_copy", split_copy_eh_vh)

		.def("split_copy", [](TriMesh& _self, OM::EdgeHandle _eh, np_point_t _arr) {
				return _self.split_copy(_eh, Point(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		.def("split_copy", [](TriMesh& _self, OM::FaceHandle _fh, np_point_t _arr) {
				return _self.split_copy(_fh, Point(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		.def("opposite_vh", &TriMesh::opposite_vh)
		.def("opposite_he_opposite_vh", &TriMesh::opposite_he_opposite_vh)

		.def("vertex_split", vertex_split_vh)

		.def("vertex_split", [](TriMesh& _self, np_point_t _arr, OM::VertexHandle _v1, OM::VertexHandle _vl, OM::VertexHandle _vr) {
				return _self.vertex_split(Point(_arr.at(0), _arr.at(1), _arr.at(2)), _v1, _vl, _vr);
			})

		.def("is_flip_ok", &TriMesh::is_flip_ok)
		.def("flip", &TriMesh::flip)

		.def("face_vertex_indices", &face_vertex_indices_trimesh)
		.def("fv_indices", &face_vertex_indices_trimesh)
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
	typedef typename Mesh::Scalar Scalar;
	typedef typename Mesh::Point  Point;
	typedef typename Mesh::Normal Normal;
	typedef typename Mesh::Color  Color;

	typedef typename Mesh::TexCoord1D TexCoord1D;
	typedef typename Mesh::TexCoord2D TexCoord2D;
	typedef typename Mesh::TexCoord3D TexCoord3D;
	typedef typename Mesh::TextureIndex TextureIndex;

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

	// Low-level adding new items
	OM::VertexHandle (Mesh::*new_vertex_void )(void                        ) = &Mesh::new_vertex;
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

	// Adding items to a mesh
	OM::FaceHandle (Mesh::*add_face_3_vh)(OM::VertexHandle, OM::VertexHandle, OM::VertexHandle) = &Mesh::add_face;
	OM::FaceHandle (Mesh::*add_face_list)(const std::vector<OM::VertexHandle>&) = &Mesh::add_face;

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

	//======================================================================
	//  PolyMeshT Function Pointers
	//======================================================================

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

	void (Mesh::*split_copy_fh_vh)(OM::FaceHandle, OM::VertexHandle) = &Mesh::split_copy;

	//======================================================================
	//  Mesh Type
	//======================================================================

	py::class_<Mesh> class_mesh(m, _name);

	class_mesh
		.def(py::init<>())

		.def(py::init([](py::array_t<typename Point::value_type> _points, py::array_t<int> _faces) {
				Mesh mesh;

				// return if _points is empty
				if (_points.size() == 0) {
					return mesh;
				}

				// _points is not empty, throw if _points has wrong shape
				if (_points.ndim() != 2 || _points.shape(1) != 3) {
					PyErr_SetString(PyExc_RuntimeError, "Array 'points' must have shape (n, 3)");
					throw py::error_already_set();
				}

				for (ssize_t i = 0; i < _points.shape(0); ++i) {
					mesh.add_vertex(Point(_points.at(i, 0), _points.at(i, 1), _points.at(i, 2)));
				}

				// return if _faces is empty
				if (_faces.size() == 0) {
					return mesh;
				}

				// _faces is not empty, throw if _faces has wrong shape
				if (_faces.ndim() != 2 || _faces.shape(1) < 3) {
					PyErr_SetString(PyExc_RuntimeError, "Array 'face_vertex_indices' must have shape (n, m) with m > 2");
					throw py::error_already_set();
				}

				for (ssize_t i = 0; i < _faces.shape(0); ++i) {
					std::vector<OM::VertexHandle> vhandles;
					for (ssize_t j = 0; j < _faces.shape(1); ++j) {
						if (_faces.at(i, j) >= 0 && _faces.at(i, j) < _points.shape(0)) {
							vhandles.push_back(OM::VertexHandle(_faces.at(i, j)));
						}
					}
					if (vhandles.size() >= 3) {
						mesh.add_face(vhandles);
					}
				}

				return mesh;
			}), py::arg("points"), py::arg("face_vertex_indices")=py::array_t<int>())

		//======================================================================
		//  Copy interface
		//======================================================================

		.def("__copy__", &Mesh::py_copy)
		.def("__deepcopy__", &Mesh::py_deepcopy)

		//======================================================================
		//  KernelT
		//======================================================================

		.def("reserve", &Mesh::reserve)

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

		.def("is_deleted", [](Mesh& _self, OM::VertexHandle _h) {
				if (!_self.has_vertex_status()) return false;
				return _self.status(_h).deleted();
			})

		.def("set_deleted", [](Mesh& _self, OM::VertexHandle _h, bool _val) {
				if (!_self.has_vertex_status()) _self.request_vertex_status();
				_self.status(_h).set_deleted(_val);
			})

		.def("is_deleted", [](Mesh& _self, OM::HalfedgeHandle _h) {
				if (!_self.has_halfedge_status()) return false;
				return _self.status(_h).deleted();
			})

		.def("set_deleted", [](Mesh& _self, OM::HalfedgeHandle _h, bool _val) {
				if (!_self.has_halfedge_status()) _self.request_halfedge_status();
				_self.status(_h).set_deleted(_val);
			})

		.def("is_deleted", [](Mesh& _self, OM::EdgeHandle _h) {
				if (!_self.has_edge_status()) return false;
				return _self.status(_h).deleted();
			})

		.def("set_deleted", [](Mesh& _self, OM::EdgeHandle _h, bool _val) {
				if (!_self.has_edge_status()) _self.request_edge_status();
				_self.status(_h).set_deleted(_val);
			})

		.def("is_deleted", [](Mesh& _self, OM::FaceHandle _h) {
				if (!_self.has_face_status()) return false;
				return _self.status(_h).deleted();
			})

		.def("set_deleted", [](Mesh& _self, OM::FaceHandle _h, bool _val) {
				if (!_self.has_face_status()) _self.request_face_status();
				_self.status(_h).set_deleted(_val);
			})

		.def("request_vertex_normals", &Mesh::request_vertex_normals)
		.def("request_vertex_colors", &Mesh::request_vertex_colors)
		.def("request_vertex_texcoords1D", &Mesh::request_vertex_texcoords1D)
		.def("request_vertex_texcoords2D", &Mesh::request_vertex_texcoords2D)
		.def("request_vertex_texcoords3D", &Mesh::request_vertex_texcoords3D)
		.def("request_halfedge_normals", &Mesh::request_halfedge_normals)
		.def("request_halfedge_colors", &Mesh::request_halfedge_colors)
		.def("request_halfedge_texcoords1D", &Mesh::request_halfedge_texcoords1D)
		.def("request_halfedge_texcoords2D", &Mesh::request_halfedge_texcoords2D)
		.def("request_halfedge_texcoords3D", &Mesh::request_halfedge_texcoords3D)
		.def("request_edge_colors", &Mesh::request_edge_colors)
		.def("request_face_normals", &Mesh::request_face_normals)
		.def("request_face_colors", &Mesh::request_face_colors)
		.def("request_face_texture_index", &Mesh::request_face_texture_index)

		.def("release_vertex_normals", &Mesh::release_vertex_normals)
		.def("release_vertex_colors", &Mesh::release_vertex_colors)
		.def("release_vertex_texcoords1D", &Mesh::release_vertex_texcoords1D)
		.def("release_vertex_texcoords2D", &Mesh::release_vertex_texcoords2D)
		.def("release_vertex_texcoords3D", &Mesh::release_vertex_texcoords3D)
		.def("release_halfedge_normals", &Mesh::release_halfedge_normals)
		.def("release_halfedge_colors", &Mesh::release_halfedge_colors)
		.def("release_halfedge_texcoords1D", &Mesh::release_halfedge_texcoords1D)
		.def("release_halfedge_texcoords2D", &Mesh::release_halfedge_texcoords2D)
		.def("release_halfedge_texcoords3D", &Mesh::release_halfedge_texcoords3D)
		.def("release_edge_colors", &Mesh::release_edge_colors)
		.def("release_face_normals", &Mesh::release_face_normals)
		.def("release_face_colors", &Mesh::release_face_colors)
		.def("release_face_texture_index", &Mesh::release_face_texture_index)

		.def("has_vertex_normals", &Mesh::has_vertex_normals)
		.def("has_vertex_colors", &Mesh::has_vertex_colors)
		.def("has_vertex_texcoords1D", &Mesh::has_vertex_texcoords1D)
		.def("has_vertex_texcoords2D", &Mesh::has_vertex_texcoords2D)
		.def("has_vertex_texcoords3D", &Mesh::has_vertex_texcoords3D)
		.def("has_halfedge_normals", &Mesh::has_halfedge_normals)
		.def("has_halfedge_colors", &Mesh::has_halfedge_colors)
		.def("has_halfedge_texcoords1D", &Mesh::has_halfedge_texcoords1D)
		.def("has_halfedge_texcoords2D", &Mesh::has_halfedge_texcoords2D)
		.def("has_halfedge_texcoords3D", &Mesh::has_halfedge_texcoords3D)
		.def("has_edge_colors", &Mesh::has_edge_colors)
		.def("has_face_normals", &Mesh::has_face_normals)
		.def("has_face_colors", &Mesh::has_face_colors)
		.def("has_face_texture_index", &Mesh::has_face_texture_index)

		.def("new_vertex", new_vertex_void)
		.def("new_vertex", [](Mesh& _self, py::array_t<typename Point::value_type> _arr) {
				return _self.new_vertex(Point(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

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

		.def("texture_index", [](Mesh& _self, OM::FaceHandle _h) {
				if (!_self.has_face_texture_index()) _self.request_face_texture_index();
				return _self.texture_index(_h);
			})

		.def("set_texture_index", [](Mesh& _self, OM::FaceHandle _h, TextureIndex _idx) {
				if (!_self.has_face_texture_index()) _self.request_face_texture_index();
				_self.set_texture_index(_h, _idx);
			})

		.def("texture_name", [](Mesh& _self, TextureIndex _idx) {
				OM::MPropHandleT<std::map<TextureIndex, std::string> > prop;
				if (_self.get_property_handle(prop, "TextureMapping")) {
					const auto map =  _self.property(prop);
					if (map.count(_idx) == 0) {
						throw py::index_error();
					}
					else {
						return map.at(_idx);
					}
				}
				else {
					PyErr_SetString(PyExc_RuntimeError, "Mesh has no textures.");
					throw py::error_already_set();
				}
			})

		//======================================================================
		//  BaseKernel
		//======================================================================

		.def("copy_all_properties", copy_all_properties_vh_vh_bool,
			py::arg("vh_from"), py::arg("vh_to"), py::arg("copy_build_in")=false)
		.def("copy_all_properties", copy_all_properties_hh_hh_bool,
			py::arg("hh_from"), py::arg("hh_to"), py::arg("copy_build_in")=false)
		.def("copy_all_properties", copy_all_properties_eh_eh_bool,
			py::arg("eh_from"), py::arg("eh_to"), py::arg("copy_build_in")=false)
		.def("copy_all_properties", copy_all_properties_fh_fh_bool,
			py::arg("fh_from"), py::arg("fh_to"), py::arg("copy_build_in")=false)

		//======================================================================
		//  ArrayKernel
		//======================================================================

		.def("is_valid_handle", (bool (Mesh::*)(OM::VertexHandle) const) &Mesh::is_valid_handle)
		.def("is_valid_handle", (bool (Mesh::*)(OM::HalfedgeHandle) const) &Mesh::is_valid_handle)
		.def("is_valid_handle", (bool (Mesh::*)(OM::EdgeHandle) const) &Mesh::is_valid_handle)
		.def("is_valid_handle", (bool (Mesh::*)(OM::FaceHandle) const) &Mesh::is_valid_handle)

		.def("delete_isolated_vertices", [](Mesh& _self) {
				if (!_self.has_vertex_status()) _self.request_vertex_status();
				_self.delete_isolated_vertices();
			})

		//======================================================================
		//  PolyConnectivity
		//======================================================================

		.def("assign_connectivity", assign_connectivity_poly)
		.def("assign_connectivity", assign_connectivity_tri)

		.def("add_face", add_face_3_vh)
		.def("add_face", add_face_list)

		.def("opposite_face_handle", &Mesh::opposite_face_handle)
		.def("adjust_outgoing_halfedge", &Mesh::adjust_outgoing_halfedge)
		.def("find_halfedge", &Mesh::find_halfedge)
		.def("valence", valence_vh)
		.def("valence", valence_fh)
		.def("is_simple_link", &Mesh::is_simple_link)
		.def("is_simply_connected", &Mesh::is_simply_connected)
		.def("remove_edge", &Mesh::remove_edge)
		.def("reinsert_edge", &Mesh::reinsert_edge)
		.def("triangulate", triangulate_fh)
		.def("triangulate", triangulate_void)
		.def("split_edge", &Mesh::split_edge)
		.def("split_edge_copy", &Mesh::split_edge_copy)

		.def("is_collapse_ok", [](Mesh& _self, OM::HalfedgeHandle _heh) {
				if (!_self.has_vertex_status()) _self.request_vertex_status();
				if (!_self.has_halfedge_status()) _self.request_halfedge_status();
				if (!_self.has_edge_status()) _self.request_edge_status();
				if (!_self.has_face_status()) _self.request_face_status();
				return _self.is_collapse_ok(_heh);
			})

		.def("collapse", [](Mesh& _self, OM::HalfedgeHandle _heh) {
				if (!_self.has_vertex_status()) _self.request_vertex_status();
				if (!_self.has_halfedge_status()) _self.request_halfedge_status();
				if (!_self.has_edge_status()) _self.request_edge_status();
				if (!_self.has_face_status()) _self.request_face_status();
				_self.collapse(_heh);
			})

		.def("add_vertex", [](Mesh& _self, py::array_t<typename Point::value_type> _arr) {
				return _self.add_vertex(Point(_arr.at(0), _arr.at(1), _arr.at(2)));
			})

		.def("delete_vertex", [](Mesh& _self, OM::VertexHandle _vh, bool _delete_isolated) {
				if (!_self.has_vertex_status()) _self.request_vertex_status();
				if (!_self.has_halfedge_status()) _self.request_halfedge_status();
				if (!_self.has_edge_status()) _self.request_edge_status();
				if (!_self.has_face_status()) _self.request_face_status();
				_self.delete_vertex(_vh, _delete_isolated);
			}, py::arg("vh"), py::arg("delete_isolated_vertices")=true)

		.def("delete_edge", [](Mesh& _self, OM::EdgeHandle _eh, bool _delete_isolated) {
				if (!_self.has_vertex_status() && _delete_isolated) _self.request_vertex_status();
				if (!_self.has_halfedge_status()) _self.request_halfedge_status();
				if (!_self.has_edge_status()) _self.request_edge_status();
				if (!_self.has_face_status()) _self.request_face_status();
				_self.delete_edge(_eh, _delete_isolated);
			}, py::arg("eh"), py::arg("delete_isolated_vertices")=true)

		.def("delete_face", [](Mesh& _self, OM::FaceHandle _fh, bool _delete_isolated) {
				if (!_self.has_vertex_status() && _delete_isolated) _self.request_vertex_status();
				if (!_self.has_halfedge_status()) _self.request_halfedge_status();
				if (!_self.has_edge_status()) _self.request_edge_status();
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

		.def_static("is_triangles", &Mesh::is_triangles)

		.def_readonly_static("InvalidVertexHandle", &Mesh::InvalidVertexHandle)
		.def_readonly_static("InvalidHalfedgeHandle", &Mesh::InvalidHalfedgeHandle)
		.def_readonly_static("InvalidEdgeHandle", &Mesh::InvalidEdgeHandle)
		.def_readonly_static("InvalidFaceHandle", &Mesh::InvalidFaceHandle)

		//======================================================================
		//  PolyMeshT
		//======================================================================

		.def("calc_edge_length", calc_edge_length_eh)
		.def("calc_edge_length", calc_edge_length_hh)
		.def("calc_edge_sqr_length", calc_edge_sqr_length_eh)
		.def("calc_edge_sqr_length", calc_edge_sqr_length_hh)

		.def("calc_sector_angle", &Mesh::calc_sector_angle)
		.def("calc_sector_area", &Mesh::calc_sector_area)

		.def("calc_dihedral_angle_fast", calc_dihedral_angle_fast_hh)
		.def("calc_dihedral_angle_fast", calc_dihedral_angle_fast_eh)
		.def("calc_dihedral_angle", calc_dihedral_angle_hh)
		.def("calc_dihedral_angle", calc_dihedral_angle_eh)

		.def("find_feature_edges", find_feature_edges, py::arg("angle_tresh")=OM::deg_to_rad(44.0))

		.def("split", split_fh_vh)
		.def("split", split_eh_vh)

		.def("split_copy", split_copy_fh_vh)

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

		.def("is_estimated_feature_edge", &Mesh::is_estimated_feature_edge)

		.def_static("is_polymesh", &Mesh::is_polymesh)
		.def("is_trimesh", &Mesh::is_trimesh)

		//======================================================================
		//  numpy calc_*
		//======================================================================

		.def("calc_face_normal", [](Mesh& _self, OM::FaceHandle _fh) {
				return vec2numpy(_self.calc_face_normal(_fh));
			})

		.def("calc_halfedge_normal", [](Mesh& _self, OM::HalfedgeHandle _heh, double _feature_angle) {
				if (!_self.has_face_normals()) {
					_self.request_face_normals();
					_self.update_face_normals();
				}
				return vec2numpy(_self.calc_halfedge_normal(_heh, _feature_angle));
			}, py::arg("heh"), py::arg("feature_angle")=0.8)

		.def("calc_vertex_normal", [](Mesh& _self, OM::VertexHandle _vh) {
				if (!_self.has_face_normals()) {
					_self.request_face_normals();
					_self.update_face_normals();
				}
				return vec2numpy(_self.calc_vertex_normal(_vh));
			})

		.def("calc_vertex_normal_fast", [](Mesh& _self, OM::VertexHandle _vh) {
				if (!_self.has_face_normals()) {
					_self.request_face_normals();
					_self.update_face_normals();
				}
				typename Mesh::Normal n;
				_self.calc_vertex_normal_fast(_vh, n);
				return vec2numpy(n);
			})

		.def("calc_vertex_normal_correct", [](Mesh& _self, OM::VertexHandle _vh) {
				typename Mesh::Normal n;
				_self.calc_vertex_normal_correct(_vh, n);
				return vec2numpy(n);
			})

		.def("calc_vertex_normal_loop", [](Mesh& _self, OM::VertexHandle _vh) {
				typename Mesh::Normal n;
				_self.calc_vertex_normal_loop(_vh, n);
				return vec2numpy(n);
			})

		.def("calc_face_centroid", [](Mesh& _self, OM::FaceHandle _fh) {
				return vec2numpy(_self.calc_face_centroid(_fh));
			})

		.def("calc_edge_vector", [](Mesh& _self, OM::EdgeHandle _eh) {
				return vec2numpy(_self.calc_edge_vector(_eh));
			})

		.def("calc_edge_vector", [](Mesh& _self, OM::HalfedgeHandle _heh) {
				return vec2numpy(_self.calc_edge_vector(_heh));
			})

		.def("calc_sector_vectors", [](Mesh& _self, OM::HalfedgeHandle _heh) {
				typename Mesh::Normal vec0;
				typename Mesh::Normal vec1;
				_self.calc_sector_vectors(_heh, vec0, vec1);
				return std::make_tuple(vec2numpy(vec0), vec2numpy(vec1));
			})

		.def("calc_sector_normal", [](Mesh& _self, OM::HalfedgeHandle _heh) {
				typename Mesh::Normal n;
				_self.calc_sector_normal(_heh, n);
				return vec2numpy(n);
			})

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
		//  numpy indices
		//======================================================================

		.def("vertex_vertex_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexVertexIter>)
		.def("vv_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexVertexIter>)
		.def("vertex_face_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexFaceIter>)
		.def("vf_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexFaceIter>)
		.def("vertex_edge_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexEdgeIter>)
		.def("ve_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexEdgeIter>)
		.def("vertex_outgoing_halfedge_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexOHalfedgeIter>)
		.def("voh_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexOHalfedgeIter>)
		.def("vertex_incoming_halfedge_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexIHalfedgeIter>)
		.def("vih_indices", &indices<Mesh, OM::VertexHandle, typename Mesh::VertexIHalfedgeIter>)

		.def("face_face_indices", &indices<Mesh, OM::FaceHandle, typename Mesh::FaceFaceIter>)
		.def("ff_indices", &indices<Mesh, OM::FaceHandle, typename Mesh::FaceFaceIter>)
		.def("face_edge_indices", &indices<Mesh, OM::FaceHandle, typename Mesh::FaceEdgeIter>)
		.def("fe_indices", &indices<Mesh, OM::FaceHandle, typename Mesh::FaceEdgeIter>)
		.def("face_halfedge_indices", &indices<Mesh, OM::FaceHandle, typename Mesh::FaceHalfedgeIter>)
		.def("fh_indices", &indices<Mesh, OM::FaceHandle, typename Mesh::FaceHalfedgeIter>)

		.def("edge_vertex_indices", &edge_other_indices<Mesh, FuncEdgeVertex>)
		.def("ev_indices", &edge_other_indices<Mesh, FuncEdgeVertex>)
		.def("edge_face_indices", &edge_other_indices<Mesh, FuncEdgeFace>)
		.def("ef_indices", &edge_other_indices<Mesh, FuncEdgeFace>)
		.def("edge_halfedge_indices", &edge_other_indices<Mesh, FuncEdgeHalfedge>)
		.def("eh_indices", &edge_other_indices<Mesh, FuncEdgeHalfedge>)

		.def("halfedge_vertex_indices", &halfedge_other_indices<Mesh, FuncHalfedgeVertex>)
		.def("hv_indices", &halfedge_other_indices<Mesh, FuncHalfedgeVertex>)

		.def("halfedge_to_vertex_indices", &halfedge_other_indices<Mesh, FuncHalfedgeToVertex>)
		.def("htv_indices", &halfedge_other_indices<Mesh, FuncHalfedgeToVertex>)
		.def("halfedge_from_vertex_indices", &halfedge_other_indices<Mesh, FuncHalfedgeFromVertex>)
		.def("hfv_indices", &halfedge_other_indices<Mesh, FuncHalfedgeFromVertex>)
		.def("halfedge_face_indices", &halfedge_other_indices<Mesh, FuncHalfedgeFace>)
		.def("hf_indices", &halfedge_other_indices<Mesh, FuncHalfedgeFace>)
		.def("halfedge_edge_indices", &halfedge_other_indices<Mesh, FuncHalfedgeEdge>)
		.def("he_indices", &halfedge_other_indices<Mesh, FuncHalfedgeEdge>)

		//======================================================================
		//  new property interface: single item
		//======================================================================

		.def("vertex_property", &Mesh::template py_property<OM::VertexHandle, typename Mesh::VPropHandle>)
		.def("halfedge_property", &Mesh::template py_property<OM::HalfedgeHandle, typename Mesh::HPropHandle>)
		.def("edge_property", &Mesh::template py_property<OM::EdgeHandle, typename Mesh::EPropHandle>)
		.def("face_property", &Mesh::template py_property<OM::FaceHandle, typename Mesh::FPropHandle>)

		.def("set_vertex_property", &Mesh::template py_set_property<OM::VertexHandle, typename Mesh::VPropHandle>)
		.def("set_halfedge_property", &Mesh::template py_set_property<OM::HalfedgeHandle, typename Mesh::HPropHandle>)
		.def("set_edge_property", &Mesh::template py_set_property<OM::EdgeHandle, typename Mesh::EPropHandle>)
		.def("set_face_property", &Mesh::template py_set_property<OM::FaceHandle, typename Mesh::FPropHandle>)

		.def("has_vertex_property", &Mesh::template py_has_property<OM::VertexHandle>)
		.def("has_halfedge_property", &Mesh::template py_has_property<OM::HalfedgeHandle>)
		.def("has_edge_property", &Mesh::template py_has_property<OM::EdgeHandle>)
		.def("has_face_property", &Mesh::template py_has_property<OM::FaceHandle>)

		.def("remove_vertex_property", &Mesh::template py_remove_property<OM::VertexHandle>)
		.def("remove_halfedge_property", &Mesh::template py_remove_property<OM::HalfedgeHandle>)
		.def("remove_edge_property", &Mesh::template py_remove_property<OM::EdgeHandle>)
		.def("remove_face_property", &Mesh::template py_remove_property<OM::FaceHandle>)

		//======================================================================
		//  new property interface: generic
		//======================================================================

		.def("vertex_property", &Mesh::template py_property_generic<OM::VertexHandle, typename Mesh::VPropHandle>)
		.def("halfedge_property", &Mesh::template py_property_generic<OM::HalfedgeHandle, typename Mesh::HPropHandle>)
		.def("edge_property", &Mesh::template py_property_generic<OM::EdgeHandle, typename Mesh::EPropHandle>)
		.def("face_property", &Mesh::template py_property_generic<OM::FaceHandle, typename Mesh::FPropHandle>)

		.def("set_vertex_property", &Mesh::template py_set_property_generic<OM::VertexHandle, typename Mesh::VPropHandle>)
		.def("set_halfedge_property", &Mesh::template py_set_property_generic<OM::HalfedgeHandle, typename Mesh::HPropHandle>)
		.def("set_edge_property", &Mesh::template py_set_property_generic<OM::EdgeHandle, typename Mesh::EPropHandle>)
		.def("set_face_property", &Mesh::template py_set_property_generic<OM::FaceHandle, typename Mesh::FPropHandle>)

		//======================================================================
		//  new property interface: array
		//======================================================================

		.def("vertex_property_array", &Mesh::template py_property_array<OM::VertexHandle, typename Mesh::VPropHandle>)
		.def("halfedge_property_array", &Mesh::template py_property_array<OM::HalfedgeHandle, typename Mesh::HPropHandle>)
		.def("edge_property_array", &Mesh::template py_property_array<OM::EdgeHandle, typename Mesh::EPropHandle>)
		.def("face_property_array", &Mesh::template py_property_array<OM::FaceHandle, typename Mesh::FPropHandle>)

		.def("set_vertex_property_array", &Mesh::template py_set_property_array<OM::VertexHandle, typename Mesh::VPropHandle>)
		.def("set_halfedge_property_array", &Mesh::template py_set_property_array<OM::HalfedgeHandle, typename Mesh::HPropHandle>)
		.def("set_edge_property_array", &Mesh::template py_set_property_array<OM::EdgeHandle, typename Mesh::EPropHandle>)
		.def("set_face_property_array", &Mesh::template py_set_property_array<OM::FaceHandle, typename Mesh::FPropHandle>)

		//======================================================================
		//  new property interface: copy
		//======================================================================

		.def("copy_property", &Mesh::template py_copy_property<OM::VertexHandle, typename Mesh::VPropHandle>)
		.def("copy_property", &Mesh::template py_copy_property<OM::HalfedgeHandle, typename Mesh::HPropHandle>)
		.def("copy_property", &Mesh::template py_copy_property<OM::EdgeHandle, typename Mesh::EPropHandle>)
		.def("copy_property", &Mesh::template py_copy_property<OM::FaceHandle, typename Mesh::FPropHandle>)
		;

	expose_type_specific_functions(class_mesh);

}

#endif

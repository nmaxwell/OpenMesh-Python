#ifndef OPENMESH_PYTHON_ITERATOR_HH
#define OPENMESH_PYTHON_ITERATOR_HH

#include "MeshTypes.hh"

#include <pybind11/pybind11.h>
namespace py = pybind11;


/**
 * Wrapper for mesh item iterators.
 *
 * This class template is used to wrap mesh item iterators for %Python. It
 * implements %Python's iterator protocol (the magic methods \_\_iter\_\_ and
 * \_\_next\_\_).
 *
 * @tparam Iterator An iterator type.
 * @tparam n_items A member function pointer that points to the mesh function
 * that returns the number of items to iterate over (e.g. n_vertices).
 */
template<class Iterator, size_t (OpenMesh::ArrayKernel::*n_items)() const>
class IteratorWrapperT {
	public:

		/**
		 * Constructor
		 *
		 * @param _mesh The mesh that contains the items to iterate over.
		 * @param _hnd The handle of the first item to iterate over.
		 * @param _skip Specifies if deleted/hidden elements are skipped.
		 */
		IteratorWrapperT(const PolyMesh& _mesh, typename Iterator::value_type _hnd, bool _skip = false) :
			mesh_(_mesh), n_items_(n_items),
			iterator_(_mesh, _hnd, _skip),
			iterator_end_(_mesh, typename Iterator::value_type(int((_mesh.*n_items)()))) {
		}

		/**
		 * Constructor
		 *
		 * @param _mesh The mesh that contains the items to iterate over.
		 * @param _hnd The handle of the first item to iterate over.
		 * @param _skip Specifies if deleted/hidden elements are skipped.
		 */
		IteratorWrapperT(const TriMesh& _mesh, typename Iterator::value_type _hnd, bool _skip = false) :
			mesh_(_mesh), n_items_(n_items),
			iterator_(_mesh, _hnd, _skip),
			iterator_end_(_mesh, typename Iterator::value_type(int((_mesh.*n_items)()))) {
		}

		/**
		 * Implementation of %Python's \_\_iter\_\_ magic method.
		 *
		 * @return This iterator.
		 */
		IteratorWrapperT iter() const {
			return *this;
		}

		/**
		 * Implementation of %Python's \_\_next\_\_ magic method.
		 *
		 * @return The next item. Raises a %Python StopIteration exception if
		 * there are no more items.
		 */
		typename Iterator::value_type next() {
			if (iterator_ != iterator_end_) {
				typename Iterator::value_type res = *iterator_;
				++iterator_;
				return res;
			}
			else {
				throw py::stop_iteration();
			}
			return typename Iterator::value_type();
		}

		/**
		 * Implementation of %Python's \_\_len\_\_ magic method.
		 *
		 * @return The number of items in the mesh.
		 */
		unsigned int len() const {
			return (mesh_.*n_items_)();
		}

	private:
		const OpenMesh::PolyConnectivity& mesh_;
		size_t (OpenMesh::ArrayKernel::*n_items_)() const;
		Iterator iterator_;
		Iterator iterator_end_;
};

/**
 * Expose an iterator type to %Python.
 *
 * @tparam Iterator An iterator type.
 * @tparam n_items A member function pointer that points to the mesh function
 * that returns the number of items to iterate over (e.g. n_vertices).
 *
 * @param _name The name of the iterator type to be exposed.
 *
 * @note %Iterators are wrapped by IteratorWrapperT before they are exposed to
 * %Python, i.e. they are not exposed directly. This means that iterators
 * that are passed from %Python to C++ are instances of IteratorWrapperT.
 */
template<class Iterator, size_t (OpenMesh::ArrayKernel::*n_items)() const>
void expose_iterator(py::module& m, const char *_name) {
	py::class_<IteratorWrapperT<Iterator, n_items> >(m, _name)
		.def(py::init<PolyMesh&, typename Iterator::value_type>())
		.def(py::init<PolyMesh&, typename Iterator::value_type, bool>())
		.def(py::init<TriMesh&, typename Iterator::value_type>())
		.def(py::init<TriMesh&, typename Iterator::value_type, bool>())
		.def("__iter__", &IteratorWrapperT<Iterator, n_items>::iter)
		.def("__next__", &IteratorWrapperT<Iterator, n_items>::next)
		.def("__len__", &IteratorWrapperT<Iterator, n_items>::len)
		;
}

#endif

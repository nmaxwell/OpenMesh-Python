#ifndef OPENMESH_PYTHON_UTILITIES_HH
#define OPENMESH_PYTHON_UTILITIES_HH

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


template<class dtype>
py::capsule free_when_done(dtype *data) {
	return py::capsule(data, [](void *f) {
		dtype *ptr = reinterpret_cast<dtype *>(f);
		delete[] ptr;
	});
}

template<class dtype>
py::array_t<dtype, py::array::c_style | py::array::forcecast> make_c_style(py::array_t<dtype> _arr) {
	return _arr;
}

#endif

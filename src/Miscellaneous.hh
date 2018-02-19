#ifndef OPENMESH_PYTHON_MISCELLANEOUS_HH
#define OPENMESH_PYTHON_MISCELLANEOUS_HH

#include <pybind11/pybind11.h>
namespace py = pybind11;


void expose_handles(py::module& m);
void expose_status_bits_and_info(py::module& m);

#endif

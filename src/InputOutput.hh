#ifndef OPENMESH_PYTHON_INPUTOUTPUT_HH
#define OPENMESH_PYTHON_INPUTOUTPUT_HH

#include <pybind11/pybind11.h>
namespace py = pybind11;


void expose_io(py::module& m);

#endif

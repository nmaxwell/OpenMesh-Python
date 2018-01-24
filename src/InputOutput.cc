#include "InputOutput.hh"
#include "MeshTypes.hh"

#include <pybind11/operators.h>

namespace OM = OpenMesh;


/**
 * Expose the input/output functions and options to Python.
 */
void expose_io(py::module& m) {

	//======================================================================
	//  Options
	//======================================================================

	py::class_<OM::IO::Options> class_options(m, "Options");

	class_options
		.def(py::init<>())
		.def(py::init<OM::IO::Options::Flag>())
		.def("cleanup", &OM::IO::Options::cleanup)
		.def("clear", &OM::IO::Options::clear)
		.def("is_empty", &OM::IO::Options::is_empty)
		.def("check", &OM::IO::Options::check)
		.def("is_binary", &OM::IO::Options::is_binary)
		.def("vertex_has_normal", &OM::IO::Options::vertex_has_normal)
		.def("vertex_has_color", &OM::IO::Options::vertex_has_color)
		.def("vertex_has_texcoord", &OM::IO::Options::vertex_has_texcoord)
		.def("edge_has_color", &OM::IO::Options::edge_has_color)
		.def("face_has_normal", &OM::IO::Options::face_has_normal)
		.def("face_has_color", &OM::IO::Options::face_has_color)
		.def("face_has_texcoord", &OM::IO::Options::face_has_texcoord)
		.def("color_has_alpha", &OM::IO::Options::color_has_alpha)
		.def("color_is_float", &OM::IO::Options::color_is_float)

		.def(py::self == py::self)
		.def(py::self != py::self)
		.def(py::self -= OM::IO::Options::Flag())
		.def(py::self += OM::IO::Options::Flag())
		;

	py::enum_<OM::IO::Options::Flag>(class_options, "Flag", py::arithmetic())
		.value("Default", OM::IO::Options::Default)
		.value("Binary", OM::IO::Options::Binary)
		.value("MSB", OM::IO::Options::MSB)
		.value("LSB", OM::IO::Options::LSB)
		.value("Swap", OM::IO::Options::Swap)
		.value("VertexNormal", OM::IO::Options::VertexNormal)
		.value("VertexColor", OM::IO::Options::VertexColor)
		.value("VertexTexCoord", OM::IO::Options::VertexTexCoord)
		.value("EdgeColor", OM::IO::Options::EdgeColor)
		.value("FaceNormal", OM::IO::Options::FaceNormal)
		.value("FaceColor", OM::IO::Options::FaceColor)
		.value("FaceTexCoord", OM::IO::Options::FaceTexCoord)
		.value("ColorAlpha", OM::IO::Options::ColorAlpha)
		.value("ColorFloat", OM::IO::Options::ColorFloat)
		.export_values()
		;

	//======================================================================
	//  Functions
	//======================================================================

	bool (*read_mesh_poly        )(PolyMesh&, const std::string&                        ) = &OM::IO::read_mesh;
	bool (*read_mesh_poly_options)(PolyMesh&, const std::string&, OM::IO::Options&, bool) = &OM::IO::read_mesh;
	bool (*read_mesh_tri         )(TriMesh&,  const std::string&                        ) = &OM::IO::read_mesh;
	bool (*read_mesh_tri_options )(TriMesh&,  const std::string&, OM::IO::Options&, bool) = &OM::IO::read_mesh;

	bool (*write_mesh_poly)(const PolyMesh&, const std::string&, OM::IO::Options, std::streamsize) = &OM::IO::write_mesh;
	bool (*write_mesh_tri )(const TriMesh&,  const std::string&, OM::IO::Options, std::streamsize) = &OM::IO::write_mesh;

	m.def("read_mesh", read_mesh_poly);
	m.def("read_mesh", read_mesh_poly_options,
		py::arg("mesh"), py::arg("filename"), py::arg("opt"), py::arg("clear")=true);
	m.def("read_mesh", read_mesh_tri);
	m.def("read_mesh", read_mesh_tri_options,
		py::arg("mesh"), py::arg("filename"), py::arg("opt"), py::arg("clear")=true);

	m.def("write_mesh", write_mesh_poly,
		py::arg("mesh"), py::arg("filename"),
		py::arg("opt")=OM::IO::Options(), py::arg("precision")=6);
	m.def("write_mesh", write_mesh_tri,
		py::arg("mesh"), py::arg("filename"),
		py::arg("opt")=OM::IO::Options(), py::arg("precision")=6);
}


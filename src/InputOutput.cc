#include "InputOutput.hh"
#include "MeshTypes.hh"

#include <pybind11/operators.h>

namespace OM = OpenMesh;


template <class Mesh>
void def_read_mesh(py::module& m, const char *_name) {
	m.def(_name,
		[](
			const std::string& _filename,
			bool _binary,
			bool _msb,
			bool _lsb,
			bool _swap,
			bool _vertex_normal,
			bool _vertex_color,
			bool _vertex_tex_coord,
			bool _halfedge_tex_coord,
			bool _edge_color,
			bool _face_normal,
			bool _face_color,
			bool _face_texture_index,
			bool _color_alpha,
			bool _color_float
		)
		{
			Mesh mesh;
			OM::IO::Options options;

			if (_binary) options += OM::IO::Options::Binary;
			if (_msb) options += OM::IO::Options::MSB;
			if (_lsb) options += OM::IO::Options::LSB;
			if (_swap) options += OM::IO::Options::Swap;

			if (_vertex_normal) {
				options += OM::IO::Options::VertexNormal;
				mesh.request_vertex_normals();
			}
			if (_vertex_color) {
				options += OM::IO::Options::VertexColor;
				mesh.request_vertex_colors();
			}
			if (_vertex_tex_coord) {
				options += OM::IO::Options::VertexTexCoord;
				mesh.request_vertex_texcoords1D();
				mesh.request_vertex_texcoords2D();
				mesh.request_vertex_texcoords3D();
			}
			if (_halfedge_tex_coord) {
				options += OM::IO::Options::FaceTexCoord;
				mesh.request_halfedge_texcoords1D();
				mesh.request_halfedge_texcoords2D();
				mesh.request_halfedge_texcoords3D();
			}
			if (_edge_color) {
				options += OM::IO::Options::EdgeColor;
				mesh.request_edge_colors();
			}
			if (_face_normal) {
				options += OM::IO::Options::FaceNormal;
				mesh.request_face_normals();
			}
			if (_face_color) {
				options += OM::IO::Options::FaceColor;
				mesh.request_face_colors();
			}
			if (_face_texture_index) {
				mesh.request_face_texture_index();
			}

			if (_color_alpha) options += OM::IO::Options::ColorAlpha;
			if (_color_float) options += OM::IO::Options::ColorFloat;

			const bool ok = OM::IO::read_mesh(mesh, _filename, options);

			if (!ok) {
				const std::string msg = "File could not be read: " + _filename;
				PyErr_SetString(PyExc_RuntimeError, msg.c_str());
				throw py::error_already_set();
			}
			if (_vertex_normal && !options.vertex_has_normal()) {
				PyErr_SetString(PyExc_RuntimeError, "Vertex normals could not be read.");
				throw py::error_already_set();
			}
			if (_vertex_color && !options.vertex_has_color()) {
				PyErr_SetString(PyExc_RuntimeError, "Vertex colors could not be read.");
				throw py::error_already_set();
			}
			if (_vertex_tex_coord && !options.vertex_has_texcoord()) {
				PyErr_SetString(PyExc_RuntimeError, "Vertex texcoords could not be read.");
				throw py::error_already_set();
			}
			if (_edge_color && !options.edge_has_color()) {
				PyErr_SetString(PyExc_RuntimeError, "Edge colors could not be read.");
				throw py::error_already_set();
			}
			if (_face_normal && !options.face_has_normal()) {
				PyErr_SetString(PyExc_RuntimeError, "Face normals could not be read.");
				throw py::error_already_set();
			}
			if (_face_color && !options.face_has_color()) {
				PyErr_SetString(PyExc_RuntimeError, "Face colors could not be read.");
				throw py::error_already_set();
			}
			if (_halfedge_tex_coord && !options.face_has_texcoord()) {
				PyErr_SetString(PyExc_RuntimeError, "Halfedge texcoords could not be read.");
				throw py::error_already_set();
			}

			return mesh;
		},
		py::arg("filename"),
		py::arg("binary")=false,
		py::arg("msb")=false,
		py::arg("lsb")=false,
		py::arg("swap")=false,
		py::arg("vertex_normal")=false,
		py::arg("vertex_color")=false,
		py::arg("vertex_tex_coord")=false,
		py::arg("halfedge_tex_coord")=false,
		py::arg("edge_color")=false,
		py::arg("face_normal")=false,
		py::arg("face_color")=false,
		py::arg("face_texture_index")=false,
		py::arg("color_alpha")=false,
		py::arg("color_float")=false
	);
}

template <class Mesh>
void def_write_mesh(py::module& m) {
	m.def("write_mesh",
		[](
			const std::string& _filename,
			const Mesh& _mesh,
			bool _binary,
			bool _msb,
			bool _lsb,
			bool _swap,
			bool _vertex_normal,
			bool _vertex_color,
			bool _vertex_tex_coord,
			bool _halfedge_tex_coord,
			bool _edge_color,
			bool _face_normal,
			bool _face_color,
			bool _color_alpha,
			bool _color_float
		)
		{
			OM::IO::Options options;

			if (_binary) options += OM::IO::Options::Binary;
			if (_msb) options += OM::IO::Options::MSB;
			if (_lsb) options += OM::IO::Options::LSB;
			if (_swap) options += OM::IO::Options::Swap;

			if (_vertex_normal) options += OM::IO::Options::VertexNormal;
			if (_vertex_color) options += OM::IO::Options::VertexColor;
			if (_vertex_tex_coord) options += OM::IO::Options::VertexTexCoord;
			if (_halfedge_tex_coord) options += OM::IO::Options::FaceTexCoord;
			if (_edge_color) options += OM::IO::Options::EdgeColor;
			if (_face_normal) options += OM::IO::Options::FaceNormal;
			if (_face_color) options += OM::IO::Options::FaceColor;

			if (_color_alpha) options += OM::IO::Options::ColorAlpha;
			if (_color_float) options += OM::IO::Options::ColorFloat;

			const bool ok = OM::IO::write_mesh(_mesh, _filename, options);

			if (!ok) {
				const std::string msg = "File could not be written: " + _filename;
				PyErr_SetString(PyExc_RuntimeError, msg.c_str());
				throw py::error_already_set();
			}
		},
		py::arg("filename"),
		py::arg("mesh"),
		py::arg("binary")=false,
		py::arg("msb")=false,
		py::arg("lsb")=false,
		py::arg("swap")=false,
		py::arg("vertex_normal")=false,
		py::arg("vertex_color")=false,
		py::arg("vertex_tex_coord")=false,
		py::arg("halfedge_tex_coord")=false,
		py::arg("edge_color")=false,
		py::arg("face_normal")=false,
		py::arg("face_color")=false,
		py::arg("color_alpha")=false,
		py::arg("color_float")=false
	);
}

void expose_io(py::module& m) {
	def_read_mesh<TriMesh>(m, "read_trimesh");
	def_read_mesh<PolyMesh>(m, "read_polymesh");
	def_write_mesh<TriMesh>(m);
	def_write_mesh<PolyMesh>(m);
}


/** @file */

#ifndef OPENMESH_PYTHON_BINDINGS_HH
#define OPENMESH_PYTHON_BINDINGS_HH

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>


/**
 * Return value policy for functions that return references to objects that are
 * managed by %OpenMesh.
 */
#define OPENMESH_PYTHON_DEFAULT_POLICY py::return_value_policy::copy


struct MeshTraits : public OpenMesh::DefaultTraits {
	/** Use double precision points */
	typedef OpenMesh::Vec3d Point;

	/** Use double precision normals */
	typedef OpenMesh::Vec3d Normal;

	/** Use RGBA colors */
	typedef OpenMesh::Vec4f Color;
};

typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> TriMesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> PolyMesh;

#endif

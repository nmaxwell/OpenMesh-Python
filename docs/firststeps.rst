
***********************************
First Steps: Creating a simple mesh
***********************************

This section demonstrates how to create a new mesh, add some vertices and faces
to it and then modify the newly inserted points.

First, we will import the openmesh and numpy modules:

.. code:: python

	import openmesh as om
	import numpy as np

Next, we can create an empty mesh:

.. code:: python

	mesh = om.TriMesh()

OpenMesh provides two mesh types: One for polygonal meshes (PolyMesh) and one
for triangle meshes (TriMesh). You should use triangle meshes whenever
possible, since they are usually more efficient. In addition, some algorithms
are only implemented for triangle meshes while triangle meshes inherit the full
functionality of polygonal meshes.

Now that we have our empty mesh object we can add a couple of vertices:

.. code:: python

	vh0 = mesh.add_vertex([0, 1, 0])
	vh1 = mesh.add_vertex([1, 0, 0])
	vh2 = mesh.add_vertex([2, 1, 0])
	vh3 = mesh.add_vertex([0,-1, 0])
	vh4 = mesh.add_vertex([2,-1, 0])

The :func:`~openmesh.TriMesh.add_vertex` member function takes numpy arrays with
shape (3,) as point coordinates and returns a handle to the newly inserted
vertex. As shown in the code above we can also pass lists with 3 elements as
point coordinates. The lists are automatically converted to numpy arrays.

In order to add a new face to our mesh we have to call
:func:`~openmesh.TriMesh.add_face`. This function takes the handles of the
vertices that make up the new face and returns a handle to the newly inserted
face:

.. code:: python

	fh0 = mesh.add_face(vh0, vh1, vh2)
	fh1 = mesh.add_face(vh1, vh3, vh4)
	fh2 = mesh.add_face(vh0, vh3, vh1)

We can also pass a list of vertex handles to :func:`~openmesh.TriMesh.add_face`:

.. code:: python

	vh_list = [vh2, vh1, vh4]
	fh3 = mesh.add_face(vh_list)

Our mesh should now look like this:

.. code:: python

	#  0 ==== 2
	#  |\  0 /|
	#  | \  / |
	#  |2  1 3|
	#  | /  \ |
	#  |/  1 \|
	#  3 ==== 4

We can access the point coordinates of each vertex by calling
:func:`~openmesh.TriMesh.point`. This member function takes a vertex handle and
returns a numpy array with shape (3,):

.. code:: python

	point = mesh.point(vh0)

We can also get an array containing all points of a mesh by calling
:func:`~openmesh.TriMesh.points`. The returned array has shape (n, 3), where n
is the number of vertices:

.. code:: python

	point_array = mesh.points()

The latter is useful if we want to update all points of a mesh at once. For
example, we can translate our mesh along the x-axis like this:

.. code:: python

	point_array += np.array([1, 0, 0])

The arrays returned by :func:`~openmesh.TriMesh.point` and
:func:`~openmesh.TriMesh.points` both reference the underlying mesh data. This
means that changes made to either one of these arrays affect the original mesh.

The complete source for this section looks like this:

.. code:: python

	import openmesh as om
	import numpy as np

	mesh = om.TriMesh()

	# add a a couple of vertices to the mesh
	vh0 = mesh.add_vertex([0, 1, 0])
	vh1 = mesh.add_vertex([1, 0, 0])
	vh2 = mesh.add_vertex([2, 1, 0])
	vh3 = mesh.add_vertex([0,-1, 0])
	vh4 = mesh.add_vertex([2,-1, 0])

	# add a couple of faces to the mesh
	fh0 = mesh.add_face(vh0, vh1, vh2)
	fh1 = mesh.add_face(vh1, vh3, vh4)
	fh2 = mesh.add_face(vh0, vh3, vh1)

	# add another face to the mesh, this time using a list
	vh_list = [vh2, vh1, vh4]
	fh3 = mesh.add_face(vh_list)

	#  0 ==== 2
	#  |\  0 /|
	#  | \  / |
	#  |2  1 3|
	#  | /  \ |
	#  |/  1 \|
	#  3 ==== 4

	# get the point with vertex handle vh0
	point = mesh.point(vh0)

	# get all points of the mesh
	point_array = mesh.points()

	# translate the mesh along the x-axis
	point_array += np.array([1, 0, 0])

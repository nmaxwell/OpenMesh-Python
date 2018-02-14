.. _iterators:

*************************
Iterators and Circulators
*************************

This section demonstrates how to use mesh iterators and circulators. The
example outputs on this page are based on the mesh from the previous section,
which looks like this:

.. code:: python

    #  0 ==== 2
    #  |\  0 /|
    #  | \  / |
    #  |2  1 3|
    #  | /  \ |
    #  |/  1 \|
    #  3 ==== 4

Iterators
#########

Iterators make it possible to enumerate the items of a mesh. For example, the
code below iterates over all vertices of a mesh:

.. code:: python

    for vh in mesh.vertices():
        print(vh.idx())

Using the mesh from the previous section, this will produce the following
output:

.. code:: python

    0
    1
    2
    3
    4

.. note:: Iterators and circulators return handles to mesh items instead of
    the items themself. For example, the vertex iterator returns vertex handles
    instead of actual vertices/points. You can access the vertex coordinates
    by calling :func:`~openmesh.TriMesh.point` with the appropriate vertex
    handle.

We can also iterate over all halfedges, edges and faces by calling
:func:`~openmesh.TriMesh.halfedges`, :func:`~openmesh.TriMesh.edges` and
:func:`~openmesh.TriMesh.faces` respectively:

.. code:: python

    # iterate over all halfedges
    for heh in mesh.halfedges():
        print(heh.idx())

    # iterate over all edges
    for eh in mesh.edges():
        print(eh.idx())

    # iterate over all faces
    for fh in mesh.faces():
        print(fh.idx())

Circulators
###########

Circulators provide the means to iterate over items adjacent to another item.
For example, to iterate over the 1-ring of a vertex we can call
:func:`~openmesh.TriMesh.vv`, which is short for vertex-vertex circulator, and
pass the handle of the center vertex:

.. code:: python

    for vh in mesh.vv(vh1):
        print(vh.idx())

Using the mesh from the previous section, this will produce the following
output:

.. code:: python

    4
    3
    0
    2

We can also iterate over the adjacent halfedges, edges and faces of a vertex:

.. code:: python

    # iterate over all incoming halfedges
    for heh in mesh.vih(vh1):
        print(heh.idx())

    # iterate over all outgoing halfedges
    for heh in mesh.voh(vh1):
        print(heh.idx())

    # iterate over all adjacent edges
    for eh in mesh.ve(vh1):
        print(eh.idx())

    # iterate over all adjacent faces
    for fh in mesh.vf(vh1):
        print(fh.idx())

To iterate over the items adjacent to a face we can use the following functions:

.. code:: python

    # iterate over the face's vertices
    for vh in mesh.fv(fh0):
        print(vh.idx())

    # iterate over the face's halfedges
    for heh in mesh.fh(fh0):
        print(heh.idx())

    # iterate over the face's edges
    for eh in mesh.fe(fh0):
        print(eh.idx())

    # iterate over all edge-neighboring faces
    for fh in mesh.ff(fh0):
        print(fh.idx())

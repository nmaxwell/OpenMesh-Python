
***************************
The Halfedge Data Structure
***************************

This section describes the underlying data structure that is used to store the
mesh vertices, edges and faces, as well as their connectivity information.

There are many popular data structures used to represent polygonal meshes. For
a detailed comparison of them refer to the papers at the end of this section.

The data structure used in this project is the so called halfedge data
structure. While face-based structures store their connectivity in faces
referencing their vertices and neighbors, edge-based structures put the
connectivity information into the edges. Each edge references its two vertices,
the faces it belongs to and the two next edges in these faces. If one now splits
the edges (i.e. an edge connecting vertex A and vertex B becomes two directed
halfeges from A to B and vice versa) one gets a halfedge-based data structure.
The following figure illustrates the way connectivity is stored in this
structure:

.. image:: img/halfedge_structure.gif
    :width: 225
    :align: center

- Each vertex references one outgoing halfedge, i.e. a halfedge that starts at
  this vertex (1).
- Each face references one of the halfedges bounding it (2).
- Each halfedge provides a handle to

    - the vertex it points to (3),
    - the face it belongs to (4),
    - the next halfedge inside the face (ordered counter-clockwise) (5),
    - the opposite halfedge (6),
    - the previous halfedge in the face (7).

Having these links between the items, it is now possible to circulate around a
face in order to enumerate all its vertices, halfedges, or neighboring faces.
When starting at a vertex' halfedge and iterating to the opposite of its
previous one, one can easily circulate around this vertex and get all its
one-ring neighbors, the incoming/outgoing halfedges, or the adjacent faces.
All this functionality is encapsulated into the so-called circulators,
described in :ref:`iterators`.

.. note:: In order to efficiently classify a boundary vertex, the outgoing
    halfedge of these vertices must be a boundary halfedge (see
    :func:`~openmesh.TriMesh.is_boundary`). Whenever you modify the topology
    using low-level topology changing functions, be sure to guarantee this
    behaviour (see :func:`~openmesh.TriMesh.adjust_outgoing_halfedge`).

While the halfedge-based structures usually consume more memory than their
face-based counter-parts they have the following important advantages:

- It is easy to mix faces of arbitrary vertex count in one mesh.
- We now have an explicit representation of vertices, faces, and
  edges/halfedges. This becomes extremely useful if one has to store data per
  edge/halfedge since this can easily be modelled by member variables of these
  types.
- Circulating around a vertex in order to get its one-ring neighbors is an
  important operation for many kinds of algorithms on polygonal meshes. For
  face-based structures this leads to many if-then branchings, the halfedge
  structure provides this funcionality without conditional branching in
  constant time.

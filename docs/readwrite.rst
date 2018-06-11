
*************
I/O Functions
*************

OpenMesh provides functions that read and write meshes from and to files:
:func:`~openmesh.read_trimesh`, :func:`~openmesh.read_polymesh` and :func:`~openmesh.write_mesh`

.. code:: python

    import openmesh as om

    trimesh = om.read_trimesh("bunny.ply")
    polymesh = om.read_polymesh("bunny.ply")
    # modify mesh ...
    om.write_mesh(trimesh, "bunny.ply")

OpenMesh currently supports five file types: .obj, .off, .ply, .stl and .om

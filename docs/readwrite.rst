
*************
I/O Functions
*************

OpenMesh provides two functions that read and write meshes from and to files:
:func:`~openmesh.read_mesh` and :func:`~openmesh.write_mesh`

.. code:: python

    import openmesh as om

    mesh = om.TriMesh()

    om.read_mesh(mesh, "bunny.ply")
    # modify mesh ...
    om.write_mesh(mesh, "bunny.ply")

The file type is automatically deduced from the file extension. OpenMesh
currently supports five file types: .obj, .off, .ply, .stl and .om

The behaviour of the I/O functions can be fine-tuned by passing an instance of
the :class:`~openmesh.Options` class to either :func:`~openmesh.read_mesh` or
:func:`~openmesh.write_mesh`. When reading a file the options are used as hints,
i.e. depending on the format we can help the reader to interpret the data
correctly. When writing a file the options determine whether or not to use the
binary variant of the respective file format and the desired byte-ordering.

.. code:: python

    import openmesh as om

    mesh = om.TriMesh()

    # hint: read vertex normals
    options = om.Options()
    options += om.Options.VertexNormal
    om.read_mesh(mesh, "bunny.ply", options)

    # write binary file
    options = om.Options()
    options += om.Options.Binary
    om.write_mesh(mesh, "bunny_binary.ply", options)

The :class:`~openmesh.Options` class controls the behaviour of the I/O functions
by means of enabled/disabled bits in a bitset. The following list contains all
available option bits:

- mode bits - control binary reading/writing

    - Options.Binary
    - Options.MSB
    - Options.LSB
    - Options.Swap (MSB|LSB)

- property bits - controls which standard properties to read/write

    - Options.VertexNormal
    - Options.VertexTexCoord
    - Options.VertexColor
    - Options.FaceNormal
    - Options.FaceColor
    - Options.ColorAlpha
    - Options.ColorFloat

Multiple options can be combined using simple arithmetic:

.. code:: python

    options = om.Options()
    options += om.Options.VertexNormal
    options += om.Options.VertexColor

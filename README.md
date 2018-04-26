# OpenMesh Python Bindings

[![](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/badges/master/pipeline.svg)](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/commits/master)

OpenMesh is a versatile halfedge-based data structure for representing and manipulating polygon meshes.
The OpenMesh Python bindings are are tightly integrated with [numpy](http://www.numpy.org/) and are implemented using [pybind11](https://github.com/pybind/pybind11). 

## Example
```python
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

# write and read meshes
om.write_mesh('test.off', mesh)
mesh_2 = om.read_trimesh('test_off')
```
For further examples see the documentation or refer to the [unit tests](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/tree/master/tests).

## Installation

### Using `pip`

    pip install openmesh

### Prebuilt Binaries

We provide prebuilt wheels for manual installation with `pip` for the following configurations:
#### Linux
* [Python 2.7](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-2.7-linux)
* [Python 3.5](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-3.5-linux)

#### macOS 10.13
* [Python 2.7](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-2.7-macos)
* [Python 3.5](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-3.5-macos)

#### Windows
* [Python 3.6](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-3.6-VS2017)

### Building from source
1. recursively clone the repo
2. `cd` to repo dir
3. ensure the correct virtualenv is activated
4. `pip install -e .`

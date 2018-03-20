# OpenMesh Python Bindings
[![pipeline status](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/badges/master/pipeline.svg)](https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/commits/master)

OpenMesh python bindings implemented with
[pybind11](https://github.com/pybind/pybind11) that are tightly integrated with
[numpy](http://www.numpy.org/).

## Installing

### Prebuild Binaries

We provide prebuild wheels for installation with pip for the following configurations:
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
2. cd to repo dir
3. `pip install -e .` (or `pip install -e . --user` if you are not root or in a
   virtualenv)

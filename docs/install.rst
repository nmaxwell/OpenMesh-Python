************
Installation
************

Using `pip`
----------

    pip install openmesh

Prebuilt Binaries
-----------------

We provide prebuilt wheels for manual installation with `pip` for the following configurations:

Linux
^^^^^
* `Python 2.7 <https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-2.7-linux>`_
* `Python 3.5 <https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-3.5-linux>`_

macOS 10.13
^^^^^^^^^^^
* `Python 2.7 <https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-2.7-macos>`_
* `Python 3.5 <https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-3.5-macos>`_

Windows
^^^^^^^
* `Python 3.6 <https://www.graphics.rwth-aachen.de:9000/OpenMesh/openmesh-python/-/jobs/artifacts/master/browse/release?job=deploy-3.6-VS2017>`_

Building from source
^^^^^^^^^^^^^^^^^^^^
1. recursively clone the repo
2. `cd` to repo dir
3. ensure the correct virtualenv is activated
4. `pip install -e .`

..
    Running the tests
    #################
    
    In your cmake build directory (e.g. build/):
    
    .. code:: python
    
        ctest --verbose

cmake_minimum_required(VERSION 2.8.12)
project(openmesh)

set(OPENMESH_BUILD_PYTHON_BINDINGS OFF)

add_subdirectory(OpenMesh)
add_subdirectory(pybind11)

acg_append_files(HEADERS "src/*.hh" .)
acg_append_files(SOURCES "src/*.cc" .)

include_directories(${OPENMESH_INCLUDE_DIRS})
pybind11_add_module(openmesh ${HEADERS} ${SOURCES})
target_link_libraries(openmesh PRIVATE ${OPENMESH_LIBRARIES})

enable_testing()
add_test(
	NAME discover_tests
	COMMAND ${PYTHON_EXECUTABLE} -m unittest
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
)
set_tests_properties(
	discover_tests PROPERTIES
	ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}
)

cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

project(openmesh)

find_package(OpenMesh)
include_directories(${OPENMESH_INCLUDE_DIRS})

add_subdirectory(pybind11)
pybind11_add_module(openmesh src/Bindings.cc src/InputOutput.cc src/Miscellaneous.cc)

target_link_libraries(openmesh PRIVATE ${OPENMESH_LIBRARIES})

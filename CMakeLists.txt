cmake_minimum_required(VERSION 3.10.0)
project(FVM VERSION 1.0.0)

enable_language(Fortran)
enable_language(CXX)

find_package(MPI REQUIRED)
include_directories(/opt/OpenMpi/include)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_subdirectory(src)
add_subdirectory(inc)
# add_subdirectory(test)


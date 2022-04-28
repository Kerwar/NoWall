cmake_minimum_required(VERSION 3.2.0)

find_package(MPI REQUIRED)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_subdirectory(src)
add_subdirectory(inc)


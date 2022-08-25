cmake_minimum_required(VERSION 3.14.0)
project(NoWall VERSION 1.0.0)

enable_language(Fortran)

find_package(MPI REQUIRED)
include_directories(/opt/OpenMpi/include)

set(CMAKE_CXX_STANDARD 17)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(src)
add_subdirectory(inc)
add_subdirectory(test)
add_subdirectory(dataprocess)
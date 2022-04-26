cmake_minimum_required(VERSION 3.2.0)
project(FVM VERSION 0.1.0)

find_package(MPI REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=c++17")

include_directories(${MPI_INCLUDE_PATH})

SET(SOURCES main.cpp Grid.cpp Paralel.cpp Field.cpp FileWriter.cpp FiniteMatrix.cpp Equation.cpp FileReader.cpp Problem.cpp Variable.cpp)
SET(HEADERS Grid.hpp Paralel.hpp ForAllOperators.hpp Field.hpp FileWriter.hpp FiniteMatrix.hpp FiniteVolumeOperations.hpp Equation.hpp FileReader.hpp Problem.hpp Variable.hpp)

add_executable(FVM ${SOURCES} ${HEADERS})
target_link_libraries(FVM ${MPI_LIBRARIES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
enable_testing()


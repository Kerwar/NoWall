cmake_minimum_required(VERSION 3.2.0)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wunused-parameter -std=c++17")

include_directories(${MPI_INCLUDE_PATH})
include_directories(../inc)

SET(SOURCES main.cpp Grid.cpp Paralel.cpp Field.cpp FileWriter.cpp FiniteMatrix.cpp Equation.cpp FileReader.cpp Problem.cpp Variable.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
enable_testing()

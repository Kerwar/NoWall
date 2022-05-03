
find_package(Catch2 REQUIRED)
find_package(MPI REQUIRED)

include_directories(${CATCH2_INCLUDE_PATH})
include_directories(${MPI_INCLUDE_PATH})
include_directories("../src")
include_directories("../inc")

SET(SOURCES ../src/Grid.cpp ../src/Paralel.cpp ../src/Field.cpp test-main.cpp mpiHelpers.hpp testGrid.cpp )#test.cpp testField.cpp testGrid.cpp)
add_executable(tests ${SOURCES} ) 

target_link_libraries(tests Catch2::Catch2WithMain Catch2::Catch2 ${MPI_LIBRARIES})
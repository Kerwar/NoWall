set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wunused-parameter -Wshadow -pedantic -Wsign-compare -Wtype-limits -Wignored-qualifiers -Wempty-body -Wclobbered -Wuninitialized -std=c++17 ")

include_directories("../src")
include_directories("../inc")

SET(SOURCES datamain.cpp ../src/FileReader.cpp ../src/Field.cpp)
add_executable(datas ${SOURCES} ) 
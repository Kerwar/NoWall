set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wunused-parameter -Wshadow -pedantic -Wsign-compare -Wtype-limits -Wignored-qualifiers -Wempty-body -Wclobbered -Wuninitialized -std=c++17 ")

include_directories("../src")
include_directories("../inc")

SET(SOURCESQVSXF qvsxf.cpp datamanager.cpp ../src/filereader.cpp ../src/field.cpp)
SET(SOURCESAVSXF avsxf.cpp datamanager.cpp ../src/filereader.cpp ../src/field.cpp)
SET(SOURCESLEWIS lewis.cpp datamanager.cpp ../src/filereader.cpp ../src/field.cpp)
SET(SOURCESLEVSXF xfvsle.cpp datamanager.cpp ../src/filereader.cpp ../src/field.cpp)

add_executable(q_vs_xf ${SOURCESQVSXF} ) 
add_executable(a_vs_xf ${SOURCESAVSXF} ) 
add_executable(lewis ${SOURCESLEWIS} ) 
add_executable(le_vs_xf ${SOURCESLEVSXF} ) 

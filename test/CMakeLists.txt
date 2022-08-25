
find_package(MPI REQUIRED)

enable_testing()
include(CTest)
include(GoogleTest)

include_directories(../inc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wunused-parameter -Wshadow -pedantic -Wsign-compare -Wtype-limits -Wignored-qualifiers -Wempty-body -Wclobbered -Wuninitialized -std=c++17 -larmadillo -Werror")

macro(package_add_test_with_libraries TESTEXECUTABLE FILES LIBRARIES TEST_WORKING_DIRECTORY)
    add_executable(${TESTEXECUTABLE} ${FILES})
    target_link_libraries(${TESTEXECUTABLE} gtest gmock gtest_main ${LIBRARIES})
    gtest_discover_tests(${TESTEXECUTABLE}
        WORKING_DIRECTORY ${TEST_WORKING_DIRECTORY}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${TEST_WORKING_DIRECTORY}"
    )
    set_target_properties(${TESTEXECUTABLE} PROPERTIES FOLDER tests)
endmacro()

package_add_test_with_libraries(testField testField.cpp ${PROJECT_NAME}_lib "${PROJECT_DIR}")
package_add_test_with_libraries(testGrid testGrid.cpp ${PROJECT_NAME}_lib "${PROJECT_DIR}")
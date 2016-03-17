SET(TEST_SOURCES
    test.cpp
    etl/std.cpp
    pdb/Model.cpp
    nuc2d/SSTree.cpp
    nuc3d/Env.cpp
    nuc3d/C2A.cpp
)

INCLUDE_DIRECTORIES(.)

ADD_EXECUTABLE(nsp_test ${TEST_SOURCES})

TARGET_LINK_LIBRARIES(nsp_test ${EXTRA_LIBS})

ADD_CUSTOM_TARGET(run_test COMMAND ${CMAKE_BINARY_DIR}/test/nsp_test DEPENDS nsp_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

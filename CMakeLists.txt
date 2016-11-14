cmake_minimum_required(VERSION 2.8.7)
project(nsp)

#if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
#    message(FATAL_ERROR "g++ version must be greater than 4.8, but now is ${CMAKE_CXX_COMPILER_VERSION}")
#endif()

set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread")

include_directories(include)
include_directories(lib)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(lib)

add_custom_target(para-debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=mpicxx ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to para-debug"
)

add_custom_target(serial-debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to serial-debug"
)

add_custom_target(para-release
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to para-release"
)

add_custom_target(serial-release
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to serial-release"
)


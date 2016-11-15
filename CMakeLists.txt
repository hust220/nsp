cmake_minimum_required(VERSION 2.8.7)
project(nsp)

message("g++ version is ${CMAKE_CXX_COMPILER_VERSION}")

#if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
#    message(FATAL_ERROR "g++ version must be greater than 4.8, but now is ${CMAKE_CXX_COMPILER_VERSION}")
#endif()

IF(DEFINED P)
    SET(CMAKE_INSTALL_PREFIX ${P})
ENDIF()

SET(CMAKE_BUILD_TYPE "Release")

IF(V MATCHES "debug")
    SET(CMAKE_BUILD_TYPE "Debug")
ENDIF()

IF(V MATCHES "para")
    ADD_DEFINITIONS(-DJN_PARA)
    SET(CMAKE_CXX_COMPILER "$ENV{MPI_HOME}/bin/mpicxx")
ENDIF()

set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread")

include_directories(include)
include_directories(lib)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(lib)
cmake_minimum_required(VERSION 2.8.7)
project(nsp)

if(DEFINED P)
    set(CMAKE_INSTALL_PREFIX ${P})
endif()

set(CMAKE_BUILD_TYPE "Release")

if(V MATCHES "debug")
    SET(CMAKE_BUILD_TYPE "Debug")
endif()

if(V MATCHES "para")
    add_definitions(-DJN_PARA)
    set(CMAKE_CXX_COMPILER "$ENV{MPI_HOME}/bin/mpicxx")
endif()

message("C++ compiler is ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
        message(FATAL_ERROR "g++ version must be greater than 4.8, but now is ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
endif()

set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pthread")

include_directories(include)
include_directories(lib)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(lib)
add_subdirectory(test)

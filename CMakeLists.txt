cmake_minimum_required(VERSION 2.8.7)
project(nsp)

find_package(Boost 1.55)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    set(EXTRA_LIBS ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif()

include_directories(include)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++14")
if (NOT MSVC)
    ADD_COMPILE_OPTIONS("-std=c++14")
    ADD_COMPILE_OPTIONS("-g")
endif()

add_subdirectory(src)
add_subdirectory(test)

add_custom_target(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
)


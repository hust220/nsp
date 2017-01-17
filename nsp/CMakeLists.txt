file(GLOB_RECURSE NSP_SOURCES *.cpp)

#link_directories(${CMAKE_BINARY_DIR}/jnpacks) 
add_executable(nsp ${NSP_SOURCES} ${jnpacks_sources})
target_link_libraries(nsp ${EXTRA_LIBS})
#target_link_libraries(nsp libjnpacks.a)

install(TARGETS nsp DESTINATION bin)


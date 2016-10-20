file(GLOB_RECURSE jian_sources jian/*.cpp)
add_library(jian SHARED ${jian_sources})
install(TARGETS jian DESTINATION lib)
install(DIRECTORY jian DESTINATION include FILES_MATCHING PATTERN *.hpp)

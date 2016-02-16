#ifndef JIAN_ETL_LIB
#define JIAN_ETL_LIB

#include <iostream>
#include <string>
#include <cstdlib>

namespace jian {

inline int die(const std::string &str, int ret_val = 1) {
    std::cerr << str << std::endl;
    std::exit(ret_val);
}

inline std::string env(const std::string &str) {
    char *path = std::getenv(str.c_str());
    return (path ? path : throw "Please set the '"+str+"' environment variable!");
}

} // namespace jian

#endif






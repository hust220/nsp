#ifndef JIAN_UTIL_FILE
#define JIAN_UTIL_FILE

#include "std.h"

namespace jian {
namespace file {

inline std::string name(const std::string &file_path) {
    std::smatch result;
    if (std::regex_match(file_path, result, std::regex("(.+)(\\.[^.]+)"))) { return result[1];
    } else return file_path;
}

inline std::string type(const std::string &file_path) {
    std::smatch result;
    if (std::regex_match(file_path, result, std::regex("(.+)(\\.[^.]+)"))) { return result[2];
    } else return "";
}

} // namespace file
} // namespace jian

#endif


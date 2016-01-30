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

template<typename Fn>
inline void each_line(const std::string &file_name, Fn &&f) {
    std::ifstream ifile(file_name.c_str());
    std::string line; int num_line = 0; 
    while (ifile) {
        num_line++;
        std::getline(ifile, line);
        if (!f(line, num_line)) break;
    }
    ifile.close();
}

} // namespace file
} // namespace jian

#endif


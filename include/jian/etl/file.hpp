#ifndef JIAN_UTIL_FILE
#define JIAN_UTIL_FILE

#include <regex>
#include "string.hpp"

namespace jian {
namespace file {

#define EACH_LINE(f, c) \
    {\
        std::ifstream ifile(f);\
        int N = 0;\
        std::string L;\
        while (std::getline(ifile, L)) {\
            ifile >> L;\
            c;\
            N++;\
        }\
        ifile.close();\
    }\

#define EACH_SPLIT_LINE(f, t, c) \
    {\
        std::ifstream ifile(f);\
        int N = 0;\
        std::string L;\
        std::vector<std::string> F;\
        while (std::getline(ifile, L)) {\
            jian::tokenize(L, F, t);\
            c;\
            N++;\
        }\
        ifile.close();\
    }\

//inline std::string name(const std::string &file_path) {
//    std::smatch result;
//    if (std::regex_match(file_path, result, std::regex("(.+)(\\.[^.]+)"))) { return result[1];
//    } else return file_path;
//}

inline std::string name(const std::string &file_path) {
    std::smatch result;
    if (std::regex_match(file_path, result, std::regex("(.*/)*([^/]+)(\\.[^.]+)"))) { return result[2];
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


#pragma once

#include "string.hpp"
#include <regex>
#include <fstream>

namespace jian {

struct file {
    static std::string name(const std::string &file_path);
    static std::string type(const std::string &file_path);
    static void clear_file(const std::string &file_name);

//    template<typename Fn>
//    void each_line(const std::string &file_name, Fn &&f) {
//        std::ifstream ifile(file_name.c_str());
//        std::string line; int num_line = 0; 
//        while (ifile) {
//            num_line++;
//            std::getline(ifile, line);
//            if (!f(line, num_line)) break;
//        }
//        ifile.close();
//    }
};

#define EACH_LINE(f, c) {\
    std::ifstream ifile(f);\
    int N = 0;\
    std::string L;\
    while (std::getline(ifile, L)) {\
        c;\
        N++;\
    }\
    ifile.close();\
}\

#define EACH_SPLIT_LINE(f, t, c) {\
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

} // namespace jian


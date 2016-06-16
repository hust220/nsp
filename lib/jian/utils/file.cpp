#include <string>
#include <regex>
#include <iostream>
#include <fstream>
#include "file.hpp"

namespace jian {

std::string file::name(const std::string &file_path) {
    std::smatch result;
    if (std::regex_match(file_path, result, std::regex("(.*/)*([^/]+)(\\.[^.]+)"))) { 
        return result[2];
    } else {
        return file_path;
    }
}

std::string file::type(const std::string &file_path) {
    std::smatch result;
    if (std::regex_match(file_path, result, std::regex("(.+)(\\.[^.]+)"))) {
        return result[2];
    } else {
        return "";
    }
}

void file::clear_file(const std::string &file_name) {
    std::ofstream ofile(file_name.c_str());
    ofile.close();
}

} // namespace jian


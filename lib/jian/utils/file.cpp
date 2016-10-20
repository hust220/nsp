#include <mutex>
#include <string>
#include <regex>
#include <iostream>
#include <fstream>
#include "file.hpp"

namespace jian {

std::string file::name(const std::string &file_path) {
    int pos1 = file_path.find_last_of('/');
    int pos2 = file_path.find_last_of('.');
    if (pos1 == std::string::npos) {
        pos1 = -1;
    }
    if (pos2 == std::string::npos || pos2 < pos1) {
        pos2 = file_path.size();
    }
    return file_path.substr(pos1+1, pos2-pos1-1);
}

std::string file::type(const std::string &file_path) {
    int pos = file_path.find_last_of('.');
    if (pos == std::string::npos) {
        return "";
    } else {
        return file_path.substr(pos+1);
    }
}

void file::clean(const std::string &file_name) {
    std::ofstream ofile(file_name.c_str());
    ofile.close();
}

} // namespace jian


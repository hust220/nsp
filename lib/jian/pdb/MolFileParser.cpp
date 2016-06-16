#include <map>
#include <string>
#include "MolFileParser.hpp"

namespace jian {

std::map<std::string, std::function<MolFileParser *(const std::string &)>> MolFileParser::s_parsers;

MolFileParser::MolFileParser(const std::string &f) {
    ifile.open(f);
    if (!ifile) {
        throw std::string("Open file '") + f + "' failed!";
    }
}

MolFileParser::~MolFileParser() {
    ifile.close();
}

} // namespace jian


#include <map>
#include <string>
#include "MolFileParser.hpp"

namespace jian {

MolFileParser::factory_t MolFileParser::s_parsers;

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


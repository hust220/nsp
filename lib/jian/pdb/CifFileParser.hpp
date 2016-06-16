#pragma once

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "MolFileParser.hpp"

namespace jian {

class Cif;

class CifFileParser : public MolFileParser {
public:
    CifFileParser(const std::string &);
    ~CifFileParser();
    MolParsedLine *line();
private:
    Cif *_cif;
};

} // namespace jian


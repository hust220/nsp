#pragma once

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "MolParser.hpp"

namespace jian {

class Cif;

class CifFileParser : public MolParser {
public:
    CifFileParser(const std::string &);
    ~CifFileParser();
    MolParsedLine *getline();
private:
    Cif *_cif;
};

} // namespace jian


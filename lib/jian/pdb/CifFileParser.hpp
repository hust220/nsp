#pragma once

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "molstream.hpp"

namespace jian {

class Cif;

class CifFileParser : public molstream {
public:
    CifFileParser(const std::string &);
    ~CifFileParser();
    MolParsedLine *getline();
private:
    Cif *_cif;
};

} // namespace jian


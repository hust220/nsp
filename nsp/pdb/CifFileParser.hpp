#pragma once

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "MolParser.hpp"

BEGIN_JN

class Cif;

class CifFileParser : public MolParser {
public:
    CifFileParser(const S &);
    ~CifFileParser();
    MolParsedLine *getline();
private:
    Cif *_cif;
};

END_JN


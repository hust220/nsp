#include <string>
#include "pdb_parser.hpp"

namespace jian {

class PdbFileParser : public MolParser {
public:
    PdbFileParser(const S &f);
    MolParsedLine *getline();
};

}


#include <string>
#include "MolParser.hpp"

namespace jian {

class PdbFileParser : public MolParser {
public:
    PdbFileParser(const std::string &f);
    MolParsedLine *getline();
};

} // namespace jian


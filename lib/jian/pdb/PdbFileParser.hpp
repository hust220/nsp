#include <string>
#include "MolFileParser.hpp"

namespace jian {

class PdbFileParser : public MolFileParser {
public:
    PdbFileParser(const std::string &f);
    MolParsedLine *line();
};

} // namespace jian


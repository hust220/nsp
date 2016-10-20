#include <string>
#include "molstream.hpp"

namespace jian {

class PdbFileParser : public molstream {
public:
    PdbFileParser(const std::string &f);
    MolParsedLine *getline();
};

} // namespace jian


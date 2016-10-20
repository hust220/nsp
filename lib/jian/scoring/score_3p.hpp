#pragma

#include "../pdb/Model.hpp"

namespace jian {
namespace scoring {

double score_stacking_3p(const Residue &r1, const Residue &r2);

double score_pairing_3p(const Residue &r1, const Residue &r2);

void train_3p(const std::string &s);

} // namespace scoring
} // namespace jian


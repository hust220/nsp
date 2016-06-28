#pragma once

#include "../pdb/Model.hpp"

namespace jian {
namespace scoring {

double score_stacking_psb(const Residue &r1, const Residue &r2);

double score_pairing_psb(const Residue &r1, const Residue &r2);

void train_psb(const std::string &s);

void print_freqs_psb();

} // namespace scoring
} // namespace jian


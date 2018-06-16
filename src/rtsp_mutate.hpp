#pragma once

#include <string>
#include "pdb.hpp"

namespace jian {

Residue mutate(const Residue &r, const S &name);

Chain mutate(const Chain &m, const S &seq, const S &type);

}


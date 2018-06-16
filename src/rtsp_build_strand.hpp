#pragma once

#include "pdb.hpp"

namespace jian {

Chain build_strand(Str seq, Str ss);

void sample_strand(Chain &strand, Str ss);

}


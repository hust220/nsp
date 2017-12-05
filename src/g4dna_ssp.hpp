#pragma once

#include "jian.hpp"

BEGIN_JN

using g4dna_ss_info_t = Deque<Deque<Int>>;

Deque<g4dna_ss_info_t> g4dna_ssp(Str seq);

END_JN


#pragma once

#include "pdb.hpp"

BEGIN_JN

Chain build_loop(Str seq, Str ss);

Chain init_loop(Str seq, Str ss);

void sample_loop(Chain &l, Str ss);

END_JN

#pragma once

#include "../pdb.hpp"

BEGIN_JN

Chain build_strand(Str seq, Str ss);

void sample_strand(Chain &strand, Str ss);

END_JN


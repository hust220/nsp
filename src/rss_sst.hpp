#pragma once

#include "yield_tree.hpp"

namespace jian {

class SSE;

struct SST {
    SSE *head;
};

SST sst_new(Str seq, Str ss, Int h = 1);

TreeNodes<SSE> sst_nodes(SST sst);

void sst_print(SST sst, std::ostream &stream);

}

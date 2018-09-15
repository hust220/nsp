#pragma once

#include "yield_tree.hpp"

namespace jian {

class SSE;

template<typename Node_, typename Cb_>
void tree_each_node(Node_ *node, Cb_ &&cb) {
    if (node != NULL) {
        cb(node);
        tree_each_node(node->son, cb);
        tree_each_node(node->bro, cb);
    }
}

struct SST {
    SSE *head;

    template<typename Cb_>
    void each_sse(Cb_ &&cb) const {
        tree_each_node(head, cb);
    }

};

SST sst_new(Str seq, Str ss, Int h = 1);

TreeNodes<SSE> sst_nodes(SST sst);

void sst_print(SST sst, std::ostream &stream);

}

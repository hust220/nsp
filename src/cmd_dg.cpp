#include "nsp.hpp"
#include "pdb.hpp"
#include "dg.hpp"
#include "topology.hpp"

namespace jian {

static Mat make_bound(Topology *top, const std::string &seq, const std::string &ss) {
    // Make the secondary structure tree
    SST sst(seq, ss, 2);

    // Init the bound matrix
    int n = top->get_num_atoms(seq);
    Mat bound(n, n);

    // Set the bound matrix by helices and loops, respectively
    sst.each_sse([&bound, &top](SSE *sse){
        if (sse->is_helix()) {
            set_bound_helix(bound, sse, top);
        }
        else if (sse->is_loop()) {
            set_bound_loop(bound, sse, top);
        }
        else {
        }
    });

    return bound;
}

static void set_bound_helix(Mat &bound, SSE *sse, Topology *top) {
}

static void set_bound_loop(Mat &bound, SSE *sse, Topology *top) {
}

static Chain make_chain(Topology *top, const Mat &c) {
}

template<typename Chain_>
static void pdb_write_chain(const Chain_ &chain, const std::string &ofname) {
}

REGISTER_NSP_COMPONENT(dg) {
    std::string seq = par.get("seq");
    std::string ss = par.get("ss");
    std::string ofname = par.get("o");
    
    // Read topology type
    Topology *top = make_topology(par.get("top"));

    // Make bound matrix from sequence and secondary structure
    auto bound = make_bound(top, seq, ss);

    // Construct DG
    auto dg = std::make_unique<Dg>(bound, top->min_inter_distance);

    // Sample a coordinate
    auto c = dg->sample();

    // Make a chain from the coordinates
    Chain chain = make_chain(top, c);

    // Write Chain
    pdb_write_chain(chain, ofname);

    // Free topology
    delete top;
}

}


#include "nsp.hpp"
#include <jian/nuc3d/Ass.hpp>
#include <jian/nuc3d/MCpsb.hpp>
#include <jian/nuc3d/triple/TriPred.hpp>
#include <thread>
#include <jian/pdb/utils/cluster_chains.hpp>

namespace jian {

static void refine(const Par &par, const Chain &chain, int i) {
    std::ostringstream stream;
    nuc3d::MCpsb mc(par);
    mc._pred_chain = chain;
    mc.predict();
    stream << mc._name << ".3ddna." << i+1 << ".pdb";
    residues_to_file(mc._pred_chain, stream.str());
}

static void mc(const Par &par, int i) {
    nuc3d::triple::TriPred<nuc3d::mc::MCpsb> tri(par);
    tri.predict();
    std::ostringstream stream;
    stream << tri._name << ".3ddna." << i+1 << ".pdb";
    residues_to_file(tri._pred_chain, stream.str());
}

REGISTER_NSP_COMPONENT(3ddna) {
    par._pars["type"] = Par::val_t{"DNA"};
    nuc3d::Assemble ass(par);
    ass.select_templates();
    auto &ss = ass._ss;
    int n = ass._num;
    assert(n <= 8);
    auto &pred_type = par["pred_type"][0];
    if (pred_type == "duplex") {
        if (std::find_if(ss.begin(), ss.end(), [](auto &&c){return c != '.' && c != '(' && c != ')';}) != ss.end() || ass.lack_templates()) {
            std::thread t[n];
            for (int i = 0; i < n; i++) {
                ass.assemble();
                t[i] = std::thread(refine, par, ass._pred_chain, i);
                ass.sample_all_templates();
            }
            for (int i = 0; i < n; i++) {
                t[i].join();
            }
        } else {
            std::deque<Chain> chains;
            for (int i = 0; i < ass._num_sampling; i++) {
                ass.assemble();
                chains.push_back(std::move(ass._pred_chain));
                ass.sample_one_template();
            }
            auto result = pdb::cluster_chains(chains, n);
            for (int i = 0; i < n; i++) {
                std::ostringstream stream;
                stream << ass._name << ".3ddna." << i+1 << ".pdb";
                residues_to_file(chains[result[i][0]], stream.str());
            }
        }
    } else if (pred_type == "triplex") {
        std::thread t[n];
        for (int i = 0; i < n; i++) {
            t[i] = std::thread(mc, par, i);
        }
        for (int i = 0; i < n; i++) {
            t[i].join();
        }
    }
}

} // namespace jian


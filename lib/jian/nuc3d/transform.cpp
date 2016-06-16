#include "../pdb.hpp"
#include "Convert.hpp"
#include "transform.hpp"

namespace jian {

static Convert cvt;

Model to_rna(const Model &model, const std::string &seq) {
    Model m = model;
    int i = 0;
    for (auto &&chain: m) for (auto &&residue: chain) {
        cvt(residue, std::string() + seq[i]);
        i++;
    }
    return m;
}

Model to_dna(const Model &model, const std::string &seq) {
    Model m = model;
    int i = 0;
    for (auto &&chain: m) for (auto &&residue: chain) {
        cvt(residue, std::string("D") + seq[i]);
        i++;
    }
    return m;
}

Model transform(const Model &m, const std::string &seq, const std::string &type) {
    if (type == "RNA") {
        return to_rna(m, seq);
    } else if (type == "DNA") {
        return to_dna(m, seq);
    }
}

} // namespace jian


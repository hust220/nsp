#include "../pdb.hpp"
#include "Convert.hpp"
#include "transform.hpp"

BEGIN_JN

Model to_rna(const Model &model, const S &seq) {
    Model m;
    int i = 0;
    for (auto && chain : model) {
        Chain c;
        for (auto && res : chain) {
            if (seq[i] != 'X') {
                c.push_back(convert_res(res, std::string() + seq[i]));
            }
            i++;
        }
        m.push_back(c);
    }
    return m;
}

Model to_dna(const Model &model, const S &seq) {
    Model m;
    int i = 0;
    for (auto && chain : model) {
        Chain c;
        for (auto && res : chain) {
            if (seq[i] != 'X') {
                c.push_back(convert_res(res, std::string("D") + seq[i]));
            }
            i++;
        }
        m.push_back(c);
    }
    return m;
}

Model transform(const Model &m, const S &seq, const S &type) {
    if (type == "RNA") {
        return to_rna(m, seq);
    }
	else if (type == "DNA") {
        return to_dna(m, seq);
    }
	else {
		throw "jian::transform error!";
	}
}

END_JN


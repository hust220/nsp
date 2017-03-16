#include <set>
#include <vector>
#include <map>
#include "mutate.hpp"
#include <jian/geom.hpp>
#include <jian/matrix.hpp>
#include <jian/utils/Env.hpp>

BEGIN_JN

class ConvertImpl {
    public:
        using names_t = std::vector<std::string>;

        Residue sugar_rna;
        Mat mat_sugar_rna;
        Residue sugar_dna;
        Mat mat_sugar_dna;
        std::map<std::string, Residue> bases;
        std::map<std::string, Mat> mat_bases;

        Mat mat_res(const Residue &r, const names_t &names) {
            Mat m(3, 3);
            for (auto && atom : r) {
                auto t = std::find(names.begin(), names.end(), atom.name);
                if (t != names.end()) {
                    int n = std::distance(names.begin(), t);
                    for (int i = 0; i < 3; i++) m(n, i) = atom[i];
                }
            }
            return m;
        }

        Mat mat_sugar(const Residue &r) {
            return mat_res(r, { "C5*", "O3*", "C1*" });
        }

        Mat mat_base(const Residue &r) {
            auto &name = r.name;
            if (name == "A" || name == "G" || name == "DA" || name == "DG") {
                return mat_res(r, { "N9", "C4", "C8" });
            }
            else if (name == "U" || name == "C" || name == "DT" || name == "DC") {
                return mat_res(r, { "N1", "C2", "C6" });
            }
            else {
                throw "Convert error!";
            }
        }

        Residue read_res(const S &name) {
            S file_name = Env::lib() + "/RNA/pars/nuc3d/Convert/" + name + ".pdb";
            Chain chain;
            chain_read_model(chain, file_name);
            return chain[0];
        }

        ConvertImpl() {
            sugar_rna = read_res("sugar.rna");
            mat_sugar_rna = mat_sugar(sugar_rna);
            sugar_dna = read_res("sugar.dna");
            mat_sugar_dna = mat_sugar(sugar_dna);
            for (auto && name : { "A", "U", "G", "C", "DA", "DT", "DG", "DC" }) {
                bases[name] = read_res(name);
                mat_bases[name] = mat_base(bases[name]);
            }
        }

        Residue convert_res(const Residue &res, const S &name) {
            if (res.name == name) return res;
            Residue r;
            r.name = name;
            add_phos(r, res);
            add_sugar(r, res, name);
            add_base(r, res, name);
            return r;
        }

        void add_phos(Residue &r, const Residue &res) {
            for (auto && atom : res) {
                auto &name = atom.name;
                if (std::find(name.begin(), name.end(), 'P') != name.end()) {
                    r.push_back(atom);
                }
            }
        }

        void add_sugar(Residue &r, const Residue &res, const S &name) {
            Mat a, b = mat_sugar(res);
            Residue sugar;
            if (name[0] == 'D') {
                sugar = sugar_dna;
                a = mat_sugar_dna;
            }
            else {
                sugar = sugar_rna;
                a = mat_sugar_rna;
            }
            auto sp = geom::suppos(a, b);
            for (auto && atom : sugar) {
                sp.apply(atom);
                r.push_back(atom);
            }
        }

        void add_base(Residue &r, const Residue &res, const S &name) {
            Mat a = mat_bases[name], b = mat_base(res);
            Residue base = bases[name];
            auto sp = geom::suppos(a, b);
            for (auto && atom : base) {
                sp.apply(atom);
                if (atom.name != "C1*") {
                    r.push_back(atom);
                }
            }
        }

};

Residue mutate(const Residue &res, const S &name) {
    static ConvertImpl convert;
    return convert.convert_res(res, name);
}

static Chain to_rna(const Chain &chain, const S &seq) {
    int i = 0;
    Chain c;
    for (auto && res : chain) {
        if (seq[i] != 'X') {
            c.push_back(mutate(res, std::string() + seq[i]));
        }
        i++;
    }
    return c;
}

static Chain to_dna(const Chain &chain, const S &seq) {
    int i = 0;
    Chain c;
    for (auto && res : chain) {
        if (seq[i] != 'X') {
            c.push_back(mutate(res, std::string("D") + seq[i]));
        }
        i++;
    }
    return c;
}

Chain mutate(const Chain &m, const S &seq, const S &type) {
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


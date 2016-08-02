#include <set>
#include <vector>
#include <map>
#include "Convert.hpp"
#include "../geom.hpp"
#include "../matrix.hpp"
#include "../utils/Env.hpp"

namespace jian {

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
        return mat_res(r, {"C5*", "O3*", "C1*"});
    }

    Mat mat_base(const Residue &r) {
        auto &name = r.name;
        if (name == "A" || name == "G" || name == "DA" || name == "DG") {
            return mat_res(r, {"N9", "C4", "C8"});
        } else if (name == "U" || name == "C" || name == "DT" || name == "DC") {
            return mat_res(r, {"N1", "C2", "C6"});
        }
    }

    Residue read_res(const std::string &name) {
        std::string file_name = Env::lib() + "/RNA/pars/nuc3d/Convert/" + name + ".pdb";
        return residues_from_file(file_name)[0];
    }

    ConvertImpl() {
        sugar_rna = read_res("sugar.rna");
        mat_sugar_rna = mat_sugar(sugar_rna);
        sugar_dna = read_res("sugar.dna");
        mat_sugar_dna = mat_sugar(sugar_dna);
        for (auto && name : {"A", "U", "G", "C", "DA", "DT", "DG", "DC"}) {
            bases[name] = read_res(name);
            mat_bases[name] = mat_base(bases[name]);
        }
    }

    Residue convert_res(const Residue &res, const std::string &name) {
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

    void add_sugar(Residue &r, const Residue &res, const std::string &name) {
        Mat a, b = mat_sugar(res);
        Residue sugar;
        if (name[0] == 'D') {
            sugar = sugar_dna;
            a = mat_sugar_dna;
        } else {
            sugar = sugar_rna;
            a = mat_sugar_rna;
        }
        auto sp = geom::suppos(a, b);
        INIT_SUPPOS(sp);
        for (auto && atom : sugar) {
            APPLY_SUPPOS(atom, sp);
            r.push_back(atom);
        }
    }

    void add_base(Residue &r, const Residue &res, const std::string &name) {
        Mat a = mat_bases[name], b = mat_base(res);
        Residue base = bases[name];
        auto sp = geom::suppos(a, b);
        INIT_SUPPOS(sp);
        for (auto && atom : base) {
            APPLY_SUPPOS(atom, sp);
            if (atom.name != "C1*") {
                r.push_back(atom);
            }
        }
    }

};

ConvertImpl l_convert_impl;

Residue convert_res(const Residue &res, const std::string &name) {
    return l_convert_impl.convert_res(res, name);
}

} // namespace jian


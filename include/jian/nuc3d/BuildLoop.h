#ifndef JIAN_NUC3D_BUILDLOOP
#define JIAN_NUC3D_BUILDLOOP

#include "../etl.h"
#include "../pdb.h"
#include "../geom.h"
#include "../nuc2d/util.h"
#include "ParseHelix.h"
#include "LM.h"

namespace jian {
namespace nuc3d {

class BuildLoop : public virtual Rand {
public:
    using val_t = double;
    using Mat = Eigen::Matrix<val_t, -1, -1>;
    using Hinge = std::array<int, 4>;
    using Hinges = std::deque<Hinge>;
    using Helix = struct {val_t theta, phi;};
    using Helices = std::vector<Helix>;

    std::string _lib = env("NSP");
    std::map<std::string, Hinges> _cache;
    std::string type = "RNA";

    ParseHelix parse_helix;
    LM lm;
    Log log;

    Model operator ()(const std::string &seq, const std::string &ss) {
        auto hinges = ss_to_hinges(nuc2d::hinge_ss(ss));
        //print_hinges(hinges);
        if (hinges.size() == 2) return lm(seq, ss)[0];
        else return make_junction(hinges);
    }

    Model operator ()(const std::string &ss) {
        std::string seq(ss.size(), 'A');
        return (*this)(seq, ss);
    }

    void print_hinges(const Hinges &hinges) {
        for (auto &loop : hinges) {
            std::cout << loop[0];
            for (int i = 1; i < 4; i++) std::cout << '-' << loop[i];
            std::cout << ' ';
        }
        std::cout << std::endl;
    }

    Hinges ss_to_hinges(const std::string &ss) {
        if (_cache.count(ss)) return _cache[ss];
        Hinges hinges; std::deque<char> stack; std::deque<int> stack2;
        for (int i = 0; i < ss.size(); i++) {
            stack.push_back(ss[i]); stack2.push_back(i);
            int l = stack.size();
            if (l>=4 && stack[l-4]=='(' && stack[l-3]=='(' && stack[l-2]==')' && stack[l-1]==')') {
                hinges.push_back({stack2[l-4], stack2[l-3], stack2[l-2], stack2[l-1]});
                for (int j = 0; j < 4; j++) {stack.pop_back(); stack2.pop_back();}
            }
        }
        if (!stack.empty()) throw ("jian::nuc3d::BuildLoop::ss_to_hinges error! ss: "s + ss).c_str();
        return hinges;
    }

    Model make_internal_loop(const Hinges &hinges) {
        return Model();
    }

    Model make_junction(const Hinges &hinges) {
        int len = hinges.size(); Helices helices(len);
        std::vector<int> v(len); std::iota(v.begin(), v.end(), 0);
        std::shuffle(v.begin(), v.end(), _rand_engine);
        helices[v[0]] = {0, 0}; helices[v[1]] = {PI, 0};
        for (int i = 2; i < len; i++) {
            helices[v[i]] = {int(rand()*6)*PI/6, int(rand()*12)*PI/6};
        }
        for (auto &helix : helices) log(helix.theta, ':', helix.phi, ' '); log('\n');
        return helices_to_model(helices, hinges);
    }

    Model helices_to_model(const Helices &helices, const Hinges &hinges) {
        Model model; model.resize(1); model[0].resize(hinges.size()*4);
        auto pairs = read_standard_pairs();
        auto result = parse_helix(pairs);
        val_t radius = 15;
        for (int i = 0; i < helices.size(); i++) {
            auto helix = pairs;
            adjust_helix(helix, result.origin, std::vector<val_t>{
                radius*std::sin(helices[i].theta)*std::cos(helices[i].phi), 
                radius*std::sin(helices[i].theta)*std::sin(helices[i].phi), 
                radius*std::cos(helices[i].theta)},
                result.theta, result.phi, helices[i].theta, helices[i].phi);
            append_helix(model, helix, hinges[i]);
        }
        return model;
    }

    template<typename K, typename U, typename F>
    void append_helix(K &model, const U &helix, const F &hinge) {
        each_residue(helix, [&](const auto &r, int i) { model[0][hinge[i]] = r; });
    }

    Model read_standard_pairs() {
        if (type == "RNA") {
            return Model(_lib + "/RNA/pars/nuc3d/BuildLoop/pairs.pdb");
        } else if (type == "DNA") {
            return Model(_lib + "/DNA/pars/nuc3d/BuildLoop/pairs.pdb");
        }
    }

    template<typename O, typename N>
    void adjust_helix(Model &model, const O &o, const N &n, val_t theta_o, val_t phi_o, val_t theta_n, val_t phi_n) {
        auto rot = geom::suppos_axis_polar<Mat>(theta_o, phi_o, theta_n, phi_n);
        std::vector<val_t> v1 {-o[0], -o[1], -o[2]}, v2 {n[0], n[1], n[2]};
        log(o[0], ':', o[1], ':', o[2], ' '); log(n[0], ':', n[1], ':', n[2], '\n');
        for (auto &chain : model) for (auto &residue : chain) for (auto &atom : residue) {
            geom::translate(atom, v1); geom::rotate(atom, rot); geom::translate(atom, v2);
        }
    }

};

} // namespace nuc3d
} // namespace jian

#endif



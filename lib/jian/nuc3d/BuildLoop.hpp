#pragma once

#include "../utils/Debug.hpp"
#include "../utils/rand.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../nuc2d.hpp"
#include "../utils/Env.hpp"
#include "ParseHelix.hpp"

namespace jian {

class BuildLoop {
public:
    using val_t = double;
    using Mat = Eigen::Matrix<val_t, -1, -1>;
    using Hinge = std::array<int, 4>;
    using Hinges = std::deque<Hinge>;
    using Helix = struct {val_t theta, phi;};
    using Helices = std::vector<Helix>;

    std::string _lib = Env::lib();
    std::map<std::string, Hinges> _cache;
    std::string type = "RNA";

    Model operator ()(const std::string &seq, const std::string &ss) {
        return make_junction(ss_to_hinges(NucSS::hinge_ss(ss)));
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
        if (!stack.empty()) throw (std::string("jian::ss_to_hinges error! ss: ") + ss).c_str();
        return hinges;
    }

    Model make_internal_loop(const Hinges &hinges) {
        return Model();
    }

    Model make_junction(const Hinges &hinges) {
        DEBUG_IN;
        int len = hinges.size(); Helices helices(len);
        std::vector<int> v(len); std::iota(v.begin(), v.end(), 0);
        std::random_shuffle(v.begin(), v.end());
        helices[v[0]] = {0, 0}; helices[v[1]] = {PI, 0};
        for (int i = 2; i < len; i++) {
            helices[v[i]] = {int(rand()*6)*PI/6, int(rand()*12)*PI/6};
        }
        for (auto &helix : helices) Debug::print(helix.theta, ':', helix.phi, ' '); Debug::print('\n');
        DEBUG_OUT;
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
                    radius*std::cos(helices[i].theta)
                },
                result.theta, result.phi, helices[i].theta, helices[i].phi
            );
            append_helix(model, helix, hinges[i]);
        }
        return model;
    }

    void append_helix(Model &model, const Model &helix, const Hinge &hinge) {
        EACH_RES(helix, model[0][hinge[N_RES]] = RES);
    }

    Model read_standard_pairs() {
        return Model(_lib + "/" + type + "/pars/nuc3d/BuildLoop/pairs.pdb");
    }

    template<typename O, typename N>
    void adjust_helix(Model &model, const O &o, const N &n, val_t theta_o, val_t phi_o, val_t theta_n, val_t phi_n) {
        DEBUG_IN;
        auto rot = geom::suppos_axis_polar(theta_o, phi_o, theta_n, phi_n);
        std::vector<val_t> v1 {-o[0], -o[1], -o[2]}, v2 {n[0], n[1], n[2]};
        Debug::print(o[0], ':', o[1], ':', o[2], ' '); Debug::print(n[0], ':', n[1], ':', n[2], '\n');
        for (auto &chain : model) for (auto &residue : chain) for (auto &atom : residue) {
            geom::translate(atom, v1); geom::rotate(atom, rot); geom::translate(atom, v2);
        }
        DEBUG_OUT;
    }

};

} // namespace jian


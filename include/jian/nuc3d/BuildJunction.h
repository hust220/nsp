#ifndef JIAN_NUC3D_BUILDJUNCTION
#define JIAN_NUC3D_BUILDJUNCTION

#include "../pdb/Molecule.h"
#include "HelixParser.h"
#include "../util/Log.h"
#include "../util/rand.h"
#include "../geom.h"

namespace jian {
namespace nuc3d {

template<typename ModelType = pdb::RNA, 
         std::enable_if_t<std::is_same<ModelType, pdb::RNA>::value || 
                          std::is_same<ModelType, pdb::DNA>::value, int> = 42>
class BuildJunction : public virtual Rand {
public:
    using val_t = double;
    using Mat = Eigen::Matrix<val_t, -1, -1>;
    using Hinge = std::array<int, 4>;
    using Hinges = std::deque<Hinge>;
    using Helix = struct {val_t theta, phi;};
    using Helices = std::vector<Helix>;

    std::string _lib = env("NSP");
    std::map<std::string, Hinges> _cache;

    HelixParser parse_helix;
    Log log;

    auto operator ()(const std::string &ss) {
        auto hinges = ss_to_hinges(ss);
        print_hinges(hinges);
        if (hinges.size() == 2) return make_internal_loop(hinges);
        else return make_junction(hinges);
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
        if (!stack.empty()) throw "jian::nuc3d::BuildJunction::ss_to_hinges error!";
        return hinges;
    }

    auto make_internal_loop(const Hinges &hinges) {
        return ModelType();
    }

    auto make_junction(const Hinges &hinges) {
        Helices helices(hinges.size());
        helices[0].theta = 0; helices[0].phi = 0;
        std::vector<int> v(helices.size()); std::iota(v.begin(), v.end(), 0);
        for (int i = 0; i < helices.size(); i++) {
            else if (i == 1) {helices[i].theta = PI; helices[i].phi = 0;}
            else {helices[i].theta = int(rand()*6)*PI/6; helices[i].phi = int(rand()*12)*PI/6;}
        }
        for (auto &helix : helices) log(helix.theta, ':', helix.phi, ' '); log('\n');
        return helices_to_model(helices, hinges);
    }

    auto helices_to_model(const Helices &helices, const Hinges &hinges) {
        ModelType model; model.resize(1); model[0].resize(hinges.size()*4);
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

    template<typename T, typename U, typename F>
    void append_helix(T &model, const U &helix, const F &hinge) {
        pdb::each_residue(helix, [&](const auto &r, int i) { model[0][hinge[i]] = r; });
    }

    auto read_standard_pairs() {
        if (std::is_same<pdb::RNA, ModelType>::value) {
            return ModelType(_lib + "/RNA/pars/nuc3d/BuildJunction/pairs.pdb");
        } else if (std::is_same<pdb::DNA, ModelType>::value) {
            return ModelType(_lib + "/DNA/pars/nuc3d/BuildJunction/pairs.pdb");
        }
    }

    template<typename O, typename N>
    void adjust_helix(ModelType &model, const O &o, const N &n, val_t theta_o, val_t phi_o, val_t theta_n, val_t phi_n) {
        auto rot = geom::suppos_axis_polar<Mat>(theta_o, phi_o, theta_n, phi_n);
        std::vector<val_t> v1 {-o[0], -o[1], -o[2]}, v2 {n[0], n[1], n[2]};
        log(o[0], ':', o[1], ':', o[2], ' '); log(n[0], ':', n[1], ':', n[2], '\n');
        for (auto &chain : model) for (auto &residue : chain) for (auto &atom : residue) {
            geom::translate(atom, v1);
            geom::rotate(atom, rot);
            geom::translate(atom, v2);
        }
    }

};

} // namespace nuc3d
} // namespace jian

#endif



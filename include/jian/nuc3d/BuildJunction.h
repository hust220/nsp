#ifndef JIAN_NUC3D_BUILDJUNCTION
#define JIAN_NUC3D_BUILDJUNCTION

#include "../pdb/Molecule.h"
#include "HelixParser.h"
#include "../util/Log.h"

namespace jian {
namespace nuc3d {

template<typename ModelType = pdb::RNA, 
         enable_if_t<std::is_same<ModelType, pdb::RNA>::value || 
                     std::is_same<ModelType, pdb::DNA>::value, int> = 42>
class BuildJunction {
public:
    using Hinge = std::array<int, 4>;
    using Hinges = std::list<Hinge>;
    using Helix = struct {double theta, phi;};
    using Helices = std::vector<Helix>;

    std::mt19937 _rand_engine{11};
    std::uniform_real_distribution<double> _unif_distr{0, 1};
    std::string _lib = env("NSP");

    HelixParser parse_helix;
    Log log;

    auto operator ()(const std::string &ss) {
        auto hinges = ss_to_hinges(ss);
        print_hinges(hinges);
        if (hinges.size() == 2) return make_internal_loop(hinges);
        else return make_junction(hinges);
    }

    double rand() {
        return _unif_distr(_rand_engine);
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
        std::vector<Helix> helices(hinges.size());
        helices[0].theta = 0; helices[0].phi = 0;
        for (int i = 1; i < helices.size(); i++) {
            helices[i].theta = rand()*PI; helices[i].phi = rand()*2*PI;
        }
        return helices_to_model(helices, hinges);
    }

    auto helices_to_model(const Helices &helices, const Hinges &hinges) {
        ModelType model; model.resize(1); model[0].resize(helices.size()*4);
        auto pairs = read_standard_pairs();
        auto result = parse_helix(pairs);
        for (int i = 0; i < helices.size(); i++) {
            auto helix = pairs;
            translate(helix, center, std::vector<double>{0, 0, 0});
            rotate(helix, direction, helices[i]);
        }
        return model;
    }

    auto read_standard_pairs() {
        if (std::is_same<pdb::RNA, ModelType>::value) {
            return ModelType(_lib + "/RNA/pars/nuc3d/BuildJunction/pairs.pdb");
        } else if (std::is_same<pdb::DNA, ModelType>::value) {
            return ModelType(_lib + "/DNA/pars/nuc3d/BuildJunction/pairs.pdb");
        }
    }

    template<typename T1, typename T2>
    void translate(const ModelType &model, const T1 &old_center, const T2 &new_center) {}

    void rotate(const ModelType &model, const Helix &old_direct, const Helix &new_direct) {}

};

} // namespace nuc3d
} // namespace jian

#endif



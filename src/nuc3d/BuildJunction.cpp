#include "BuildJunction.h"

namespace jian {

namespace nuc3d {

Model JunctBuild::operator ()(std::string seq, std::string ss) {
//    extend_helix(_helix_len);
    Model model;
    return model;
} 

void JunctBuild::move() {
    
}

void JunctBuild::rollback() {
    
}

double JunctBuild::energy() {
    
}

double JunctBuild::min_energy() {
    
}

void JunctBuild::set_min_state() {
    
}

void JunctBuild::train(const Model &model, std::string ss) {
    std::vector<Model> helices;
    std::vector<int> strand_lens;
    std::tie(helices, strand_lens) = get_helices(model, ss);
    for (auto &&helix: helices) {
        _helices.push_back(helix_par(helix));
    }
    for (int i = 0; i < _helices.size(); i++) {
        Vector3f orig1, orig2, axis1, axis2, direct1, direct2;
        std::tie(orig1, axis1, direct1) = _helices[i];
        std::tie(orig2, axis2, direct2) = (i != _helices.size() - 1 ? _helices[i + 1] : _helices[0]);
        std::cout << strand_lens[i] << ' ' << geometry::distance(orig1, orig2) << ' ' << geometry::angle(axis1, Point(0, 0, 0), axis2) / 3.14159 * 180 << ' ' << geometry::angle(direct1, Point(0, 0, 0), direct2) / 3.14159 * 180 << std::endl;
    }
    for (auto &&helix: _helices) {
    }
}

std::pair<std::vector<Model>, std::vector<int>> JunctBuild::get_helices(const Model &model, std::string ss) {
    std::vector<Residue> residues;
    for (auto &&chain: model.chains) {
        for (auto &&residue: chain.residues) {
            residues.push_back(residue);
        }
    }

    std::vector<Model> helices;
    std::vector<int> strand_lens;
    int i;
    for (i = 0; i < ss.size(); i++) {
        if (ss[i] == '(' && ss[ss.size() - 1 - i] == ')') {
            continue;
        } else {
            if (i == 0) break;
            Model helix;
            helix.chains.resize(2);
            std::copy(residues.begin(), std::next(residues.begin(), i), std::back_inserter(helix.chains[0].residues));
            std::copy(std::next(residues.begin(), ss.size() - i), residues.end(), std::back_inserter(helix.chains[1].residues));
            helices.push_back(helix);
            break;
        }
    }

    int helix_len = i;
    int flag = i;
    for (; i < ss.size() - 1; i++) {
        if (ss[i] == '(' && ss[i + 1] == ')') {
            for (int j = i - 1; j >= 0; j--) {
                if (ss[j] == '(' && ss[2 * i - j] == ')') {
                    continue;
                } else {
                    strand_lens.push_back(j + 1 - flag);
                    flag = 2 * i - j + 1;
                    Model helix;
                    helix.chains.resize(2);
                    std::copy(std::next(residues.begin(), j + 1), std::next(residues.begin(), i + 1), std::back_inserter(helix.chains[0].residues));
                    std::copy(std::next(residues.begin(), i + 1), std::next(residues.begin(), 2 * i - j + 1), std::back_inserter(helix.chains[1].residues));
                    helices.push_back(helix);
                    break;
                }
            }
        }
    }
    strand_lens.push_back(ss.size() - helix_len - flag);

    return std::make_pair(helices, strand_lens);
}

std::tuple<Vector3f, Vector3f, Vector3f> JunctBuild::helix_par(const Model &model) {
    auto axis1 = pdb::normal_vector(model[0][0]);
    auto axis2 = pdb::normal_vector(model[1].residues.back());
    Vector3f axis = axis1 - axis2;

    auto pivot1 = pdb::pivot(model[0][0]);
    auto pivot2 = pdb::pivot(model[1].residues.back());
    Vector3f orig = 0.5 * (pivot1 + pivot2);
    Vector3f diret = pivot2 - pivot1;
    return std::make_tuple(orig, axis, diret);
}

} /// namespace nuc3d

} /// namespace jian


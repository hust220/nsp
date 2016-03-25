#pragma once

#include "../pdb.hpp"
#include "Env.hpp"

namespace jian {

class CG2AA {
public:
    std::deque<std::shared_ptr<MatrixXd>> _frags;
    std::deque<std::string> _names;
    std::string _path;
    int _num_residues_frag = 4;
    std::vector<std::string> _coarse_grained_atoms {"C5*", "O3*", "C1*"};

    CG2AA() {
        int num_atoms = _num_residues_frag * _coarse_grained_atoms.size();
        _path = Env::lib() + "/RNA/pars/nuc3d/CG2AA/";
        EACH_SPLIT_LINE((_path + "inf.txt").c_str(), " ",
            auto m = std::make_shared<MatrixXd>(num_atoms, 3);
            _names.push_back(F[0]);
            FOR((i, num_atoms), FOR((j, 3), (*m)(i, j) = JN_DBL(F[i*3+j+1])));
            _frags.push_back(m);
        );
    }

    static CG2AA &instance() {
        static CG2AA c2a;
        return c2a;
    }

    template<typename T, typename U>
    static auto cg2aa(T &&coord, U &&frag) {
        return instance()(coord, frag);
    }

    template<typename T>
    auto get_residues(int i, T &&c) {
        Trace::log(_path + _names[i] + ".pdb\n");
        auto residues = residues_from_file(_path + _names[i] + ".pdb");
        auto s = geom::suppos(*(_frags[i]), c);
        INIT_SUPPOS(s);
        EACH(res, residues, EACH(atom, res, APPLY_SUPPOS(atom, s)));
        return residues;
    }

    template<typename T, typename U>
    auto operator ()(T &&coord, U &&frag) {
        Trace::log("## CG2AA\n");
        int num_atoms = _num_residues_frag * _coarse_grained_atoms.size();
        int len = frag[1] - frag[0] + 1; 
        Chain chain;
        for (int i = 0; i < len - num_atoms + 1; i += _coarse_grained_atoms.size()) {
            MatrixXd c(num_atoms, 3); 
            FOR((j, num_atoms), FOR((k, 3), c(j, k) = coord(frag[0]+i+j, k)));
            std::deque<double> scores; 
            EACH((j, n), _frags, 
                auto r = geom::suppos(*j, c); 
                scores.push_back(r.rmsd);
            );
            auto min = std::min_element(scores.begin(), scores.end()); 
            int index = std::distance(scores.begin(), min);
            auto residues = get_residues(index, c);
            if (i == 0) FOR((j, _num_residues_frag-1), chain.push_back(residues[j]));
            chain.push_back(residues[_num_residues_frag-1]);
        }
        return chain;
    }

    template<typename T>
    static void extract_frags(T &&model) {
        std::deque<Residue> dq; 
        std::ofstream ofile(instance()._path + "inf.txt");
        auto &names = instance()._coarse_grained_atoms;
        EACH_RES(model,
            dq.push_back(RES);
            if (dq.size() >= 2 && geom::distance(dq.back()["C4*"], dq[dq.size()-2]["C4*"]) > 10) {
                dq.clear();
            } else if (dq.size() == instance()._num_residues_frag) {
                std::string name = model.name + "-" + JN_STR(N_RES); 
                ofile << name;
                for (int i = 0; i < dq.size(); i++) {
                    for (auto && atom : dq[i]) {
                        if (std::find(names.begin(), names.end(), atom.name) != names.end()) {
                            for (int j = 0; j < 3; j++) {
                                ofile << ' ' << atom[j];
                            }
                        }
                    }
                }
                ofile << '\n';
                residues_to_file(dq, instance()._path + name + ".pdb");
                dq.pop_front();
            }
        );
        ofile.close();
    }

};

} // namespace jian


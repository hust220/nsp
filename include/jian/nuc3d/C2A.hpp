#pragma once

#include "../pdb.hpp"
#include "Env.hpp"

namespace jian {

class C2A {
public:
    std::deque<std::shared_ptr<MatrixXd>> _frags;
    std::deque<std::string> _names;
    std::string _path;
    int _len_frag = 4;

    C2A() {
        _path = Env::lib() + "/RNA/pars/nuc3d/C2A/";
        EACH_SPLIT_LINE((_path + "inf.txt").c_str(), " ",
            auto m = std::make_shared<MatrixXd>(_len_frag, 3);
            _names.push_back(F[0]);
            FOR((i, _len_frag), FOR((j, 3), (*m)(i, j) = JN_DBL(F[i*3+j+1])));
            _frags.push_back(m);
        );
    }

    static C2A &instance() {
        static C2A c2a;
        return c2a;
    }

    template<typename T, typename U>
    static auto c2a(T &&coord, U &&frag) {
        return instance()(coord, frag);
    }

    template<typename T>
    auto get_residues(int i, T &&c) {
        Trace::log(_path + _names[i] + ".pdb\n");
        auto residues = residues_from_file(_path + _names[i] + ".pdb");
        auto s = geom::suppos(*(_frags[i]), c);
        EACH(res, residues, EACH(atom, res, geom::apply_suppos(atom, s)));
        return residues;
    }

    template<typename T, typename U>
    auto operator ()(T &&coord, U &&frag) {
        Trace::log("## C2A\n");
        int len = frag[1] - frag[0] + 1; 
        Chain chain;
        for (int i = 0; i < len - _len_frag + 1; i++) {
            MatrixXd c(_len_frag, 3); 
            FOR((j, _len_frag), FOR((k, 3), c(j, k) = coord(frag[0]+i+j, k)));
            std::deque<double> scores; 
            EACH((j, n), _frags, auto r = geom::suppos(*j, c); scores.push_back(r.rmsd));
            auto min = std::min_element(scores.begin(), scores.end()); 
            int index = std::distance(scores.begin(), min);
            auto residues = get_residues(index, c);
            if (i == 0) FOR((j, _len_frag-1), chain.push_back(residues[j]));
            chain.push_back(residues[_len_frag-1]);
        }
        return chain;
    }

    template<typename T>
    static void extract_frags(T &&model) {
        std::deque<Residue> dq; std::ofstream ofile(instance()._path + "inf.txt");
        EACH_RES(model,
            dq.push_back(RES);
            if (dq.size() >= 2 && geom::distance(dq.back()["C4*"], dq[dq.size()-2]["C4*"]) > 10) {
                dq.clear();
            } else if (dq.size() == instance()._len_frag) {
                std::string name = model.name + "-" + JN_STR(N_RES); ofile << name;
                FOR((i, instance()._len_frag), FOR((j, 3), auto &&atom = dq[i]["C4*"]; ofile << ' ' << atom[j]));
                ofile << '\n';
                residues_to_file(dq, instance()._path + name + ".pdb");
                dq.pop_front();
            }
        );
        ofile.close();
    }

};

} // namespace jian


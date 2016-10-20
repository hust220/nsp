#include "CG2AA.hpp"
#include "../utils/Debug.hpp"
#include "../pp.hpp"
#include "../utils/file.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../utils/Env.hpp"
#include <deque>
#include <vector>
#include <memory>

namespace jian {

class CG2AA {
public:
    std::deque<std::shared_ptr<Eigen::MatrixXd>> _frags;
    std::deque<std::string> _names;
    std::string _path;
    int _num_residues_frag = 4;
    std::vector<std::string> _coarse_grained_atoms {"C5*", "O3*", "C1*"};

    CG2AA() {
        int num_atoms = _num_residues_frag * _coarse_grained_atoms.size();
        _path = Env::lib() + "/RNA/pars/nuc3d/CG2AA/";
        EACH_SPLIT_LINE((_path + "inf.txt").c_str(), " ",
            auto m = std::make_shared<Eigen::MatrixXd>(num_atoms, 3);
            _names.push_back(F[0]);
			for (int i = 0; i < num_atoms; i++) for (int j = 0; j < 3; j++) (*m)(i, j) = JN_DBL(F[i * 3 + j + 1]);
            //FOR((i, num_atoms), FOR((j, 3), (*m)(i, j) = JN_DBL(F[i*3+j+1])));
            _frags.push_back(m);
        );
    }

    template<typename T>
    auto get_residues(int i, T &&c) {
        Debug::print(_path + _names[i] + ".pdb\n");
        Chain residues;
        chain_read_model(residues, _path + _names[i] + ".pdb");
        auto s = geom::suppos(*(_frags[i]), c);
        INIT_SUPPOS(s);
		for (auto && res : residues) for (auto && atom : res) {
			APPLY_SUPPOS(atom, s);
		}
        //EACH(res, residues, EACH(atom, res, APPLY_SUPPOS(atom, s)));
        return residues;
    }

    template<typename T, typename U>
    auto run(T &&coord, U &&frag) {
        Debug::print("## CG2AA\n");
        int num_atoms = _num_residues_frag * _coarse_grained_atoms.size();
        int len = frag[1] - frag[0] + 1; 
        Chain chain;
        for (int i = 0; i < len - num_atoms + 1; i += _coarse_grained_atoms.size()) {
            Eigen::MatrixXd c(num_atoms, 3); 
			for (int j = 0; j < num_atoms; j++) for (int k = 0; k < 3; k++) c(j, k) = coord(frag[0] + i + j, k);
            //FOR((j, num_atoms), FOR((k, 3), c(j, k) = coord(frag[0]+i+j, k)));
            std::deque<double> scores; 
			for (auto && j : _frags) {
				//EACH((j, n), _frags,
				auto r = geom::suppos(*j, c);
				scores.push_back(r.rmsd);
			}
            auto min = std::min_element(scores.begin(), scores.end()); 
            int index = std::distance(scores.begin(), min);
            auto residues = get_residues(index, c);
			if (i == 0) for (int j = 0; j < _num_residues_frag - 1; j++) chain.push_back(residues[j]);
            //if (i == 0) FOR((j, _num_residues_frag-1), chain.push_back(residues[j]));
            chain.push_back(residues[_num_residues_frag-1]);
        }
        return chain;
    }

    template<typename T>
    void extract_frags(T &&model) {
        Chain dq; 
        std::ofstream ofile(_path + "inf.txt");
        auto &names = _coarse_grained_atoms;
		int n_res = 0;
		for (auto && chain : model) for (auto && res : chain) {
			dq.push_back(res);
			if (dq.size() >= 2 && geom::distance(dq.back()["C4*"], dq[dq.size() - 2]["C4*"]) > 10) {
				dq.clear();
			}
			else if (dq.size() == _num_residues_frag) {
				std::string name = model.name + "-" + JN_STR(n_res);
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
				mol_write(dq, _path + name + ".pdb");
				dq.pop_front();
			}
			n_res++;
		}
        ofile.close();
    }

};

Chain cg2aa(const Eigen::MatrixXd &c, int beg, int end) {
    CG2AA cg2aa;
    return cg2aa.run(c, std::vector<int>{beg, end});
}

} // namespace jian


#include <deque>
#include <vector>
#include "CGpsb.hpp"
#include "../pp.hpp"
#include "../pdb.hpp"
#include "../geom.hpp"
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "../utils/Debug.hpp"

namespace jian {

class PSB2AA {
public:
    std::deque<std::shared_ptr<Eigen::MatrixXd>> _frags;
    std::deque<std::string> _names;
    std::string _path;
    int _num_residues_frag = 4;

    PSB2AA() {
        int num_atoms = _num_residues_frag * 3;
        _path = Env::lib() + "/RNA/pars/nuc3d/psb2aa/";
        EACH_SPLIT_LINE((_path + "inf.txt").c_str(), " ",
            auto m = std::make_shared<Eigen::MatrixXd>(num_atoms, 3);
            _names.push_back(F[0]);
            FOR((i, num_atoms), FOR((j, 3), (*m)(i, j) = std::stod(F[i*3+j+1])));
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
        EACH(res, residues, EACH(atom, res, APPLY_SUPPOS(atom, s)));
        return residues;
    }

    template<typename T, typename U>
    auto run(T &&coord, U &&frag) {
        std::cout << "## PSB2AA\n" << std::endl;
        int num_atoms = _num_residues_frag * 3;
        int len = frag[1] - frag[0] + 1; 
        Chain chain;
        for (int i = 0; i < len - num_atoms + 1; i += 3) {
            Eigen::MatrixXd c(num_atoms, 3); 
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

    void extract_frags(const std::string &pdb) {
        std::deque<int> dq; 
        Chain full_chain;
        chain_read_model(full_chain, pdb);
        Chain chain = CGpsb::chain(full_chain);
        std::ofstream ofile(_path + "inf.txt");
        for (int n = 0; n < chain.size(); n++) {
            dq.push_back(n);
            if (dq.size() >= 2 && geom::distance(chain[dq.back()][0], chain[dq[dq.size()-2]][0]) > 10) {
                dq.clear();
            } else if (dq.size() == _num_residues_frag) {
                std::string name = chain.model_name + "-" + std::to_string(n); 
                ofile << name;
                for (int i = 0; i < dq.size(); i++) {
                    for (auto && atom : chain[dq[i]]) {
                        for (int j = 0; j < 3; j++) {
                            ofile << ' ' << atom[j];
                        }
                    }
                }
                ofile << '\n';
                Chain temp_chain;
                for (auto && i : dq) temp_chain.push_back(full_chain[i]);
                mol_write(temp_chain, _path + name + ".pdb");
                dq.pop_front();
            }
        }
        ofile.close();
    }

};

thread_local static PSB2AA psb;

Chain CGpsb::aa(const Eigen::MatrixXd &c, int beg, int end) {
    return psb.run(c, std::vector<int>{beg, end});
}

void CGpsb::extract_frags(const std::string &pdb) {
    psb.extract_frags(pdb);
}

} // namespace jian


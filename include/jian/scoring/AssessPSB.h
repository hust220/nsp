#ifndef JIAN_SCORING_ASSESSPSB
#define JIAN_SCORING_ASSESSPSB

#include "../pdb/Molecule.h"
#include "../geom/geometry.h"
#include "../fpl.h"

namespace jian {
namespace scoring {

class AssessPSB {
public:
    using val_t = double;

    pdb::PSB _model;
    val_t _bond_length_score, _bond_angle_score, _bond_dihedral_score;
    val_t _base_pairing_score, _base_stacking_score;

    std::array<val_t, 3> _std_bond_length;
    std::array<val_t, 4> _std_bond_angle;
    std::array<val_t, 3> _std_bond_dihedral;
    std::array<val_t, 3> _weight_bond_length;
    std::array<val_t, 4> _weight_bond_angle;
    std::array<val_t, 3> _weight_bond_dihedral;
    std::array<val_t, 3> _score_bond_length;
    std::array<val_t, 4> _score_bond_angle;
    std::array<val_t, 3> _score_bond_dihedral;
    MatrixXd _pars_base_stacking;
    MatrixXd _pars_base_pairing;
    val_t _weight_base_stacking;
    val_t _weight_base_pairing;
    val_t _score_base_stacking;
    val_t _score_base_pairing;

    AssessPSB() {
        auto init_array = [](auto &arr){for (auto &&i : arr){i = 0;}};
        init_array(_std_bond_length);
        init_array(_std_bond_angle);
        init_array(_std_bond_dihedral);
        init_array(_weight_bond_length);
        init_array(_weight_bond_angle);
        init_array(_weight_bond_dihedral);
        init_array(_score_bond_length);
        init_array(_score_bond_angle);
        init_array(_score_bond_dihedral);
        _score_base_stacking = 0;
        _score_base_pairing = 0;
    }

    template<typename ModelType>
    val_t operator ()(ModelType model) { return assess(model); }

    template<typename ModelType>
    val_t assess(const ModelType &model) {
        init(model);
        assess();
        return score_bond_length() + score_bond_angle() + score_bond_dihedral() + score_base_pairing() + score_base_stacking();
    }

    template<typename ModelType>
    void analyze(const ModelType &model) {
        init(model);
        analyze();
    }

    template<typename ModelType>
    void init(const ModelType &model) {
        _model = model;
        _bond_length_score = _bond_angle_score = _bond_dihedral_score = 0;
        _base_pairing_score = _base_stacking_score = 0;
    }

    void assess() {
        assess_residue();
        assess_two_consecutive_residues();
        assess_three_consecutive_residues();
        assess_two_separate_residues();
    }

    void analyze() {
        analyze_residue();
        analyze_two_consecutive_residues();
        analyze_three_consecutive_residues();
        analyze_two_separate_residues();
    }

    template<typename Res1, typename Res2>
    val_t dist(Res1 &&res1, Res2 &&res2) {
        return geometry::distance(res1.atom("S"), res2.atom("P"));
    }

    template<typename Fn>
    void each_residue(Fn &&f) {
        for (auto &&chain : _model) for (auto &&res : chain) {
            if (res.size() == 2) continue;
            val_t len1 = geometry::distance(res[0], res[1]);
            val_t len2 = geometry::distance(res[1], res[2]);
            val_t ang1 = geometry::angle(res[0], res[1], res[2]);
            f(len1, len2, ang1);
        }
    }

    void assess_residue() {
        each_residue([&](val_t len1, val_t len2, val_t ang1){
            _score_bond_length[0] += en_bond_length(_std_bond_length[0], len1);
            _score_bond_length[1] += en_bond_length(_std_bond_length[1], len2);
            _score_bond_angle[0] += en_bond_angle(_std_bond_angle[0], ang1);
        });
    }

    void analyze_residue() {
        each_residue([](val_t len1, val_t len2, val_t ang1){
            std::cout << len1 << ' ' << len2 << ' ' << ang1 << std::endl;
        });
    }

    template<typename Fn>
    void each_two_consecutive_residues(Fn &&f) {
        for (auto &&chain : _model) {
            for (int i = 1; i < chain.size(); i++) {
                if (chain[i-1].size() == 2 or chain[i].size() == 2 or dist(chain[i-1], chain[i]) > 4.2) continue;
                auto len3 = geometry::distance(chain[i-1][1], chain[i][0]);
                auto ang2 = geometry::angle(chain[i-1][0], chain[i-1][1], chain[i][0]);
                auto ang3 = geometry::angle(chain[i-1][2], chain[i-1][1], chain[i][0]);
                auto ang4 = geometry::angle(chain[i-1][1], chain[i][0],   chain[i][1]);
                auto dih1 = geometry::dihedral(chain[i-1][0], chain[i-1][1], chain[i-1][2], chain[i][0]);
                auto dih2 = geometry::dihedral(chain[i-1][0], chain[i-1][1], chain[i][0],   chain[i][1]);
                MatrixXd dists_stacking = MatrixXd::Zero(2, 2);
                for (auto m : {0, 1}) for (auto n : {0, 1}) {
                    dists_stacking(m, n) = geometry::distance(chain[i-1][m+1], chain[i][n+1]);
                }
                val_t dih_stacking = geometry::dihedral(chain[i-1][2], chain[i-1][1], chain[i][1], chain[i][2]);
                f(len3, ang2, ang3, ang4, dih1, dih2, dists_stacking, dih_stacking);
            }
        }
    }

    void assess_two_consecutive_residues() {
        each_two_consecutive_residues([&](val_t len3, val_t ang2, val_t ang3, val_t ang4, 
                                         val_t dih1, val_t dih2, const auto &dists_stacking, val_t dih_stacking){
            _score_bond_length[2] += en_bond_length(_std_bond_length[2], len3);
            _score_bond_angle[1] += en_bond_angle(_std_bond_angle[1], ang2);
            _score_bond_angle[2] += en_bond_angle(_std_bond_angle[2], ang3);
            _score_bond_angle[3] += en_bond_angle(_std_bond_angle[3], ang4);
            _score_bond_dihedral[0] += en_bond_dihedral(_std_bond_dihedral[0], dih1);
            _score_bond_dihedral[1] += en_bond_dihedral(_std_bond_dihedral[1], dih2);
            _score_base_stacking += en_base_stacking(dists_stacking, dih_stacking);
        });
    }

    void analyze_two_consecutive_residues() {
        each_two_consecutive_residues([](val_t len3, val_t ang2, val_t ang3, val_t ang4, 
                                      val_t dih1, val_t dih2, const auto &dists_stacking, val_t dih_stacking){
            std::cout << len3 << ' ' << ang2 << ' ' << ang3 << ' ' << ang4 << ' ' << dih1 << ' ' << dih2 << std::endl;
            std::cout << dists_stacking << std::endl;
            std::cout << dih_stacking << std::endl;
        });
    }

    template<typename Fn>
    void each_three_consecutive_residues(Fn &&f) {
        for (auto &&chain : _model) for (int i = 2; i < chain.size(); i++) {
            if (chain[i-2].size() == 2 or chain[i-1].size() == 2 or chain[i].size() == 2 or 
                dist(chain[i-2], chain[i-1]) > 4.2 or dist(chain[i-1], chain[i]) > 4.2) continue;
            auto dih3 = geometry::dihedral(chain[i-2][1], chain[i-1][0], chain[i-1][1], chain[i][0]);
            f(dih3);
        }
    }

    void assess_three_consecutive_residues() {
        each_three_consecutive_residues([&](val_t dih3){
            _score_bond_dihedral[2] += en_bond_dihedral(_std_bond_dihedral[2], dih3);
        });
    }

    void analyze_three_consecutive_residues() {
        each_three_consecutive_residues([](val_t dih3){
            std::cout << dih3 << std::endl;
        });
    }

    template<typename Fn>
    void each_two_separate_residues(Fn &&f) {
        int num_res1 = 0;
        for (auto &&chain1 : _model) for (auto &&res1 : chain1) {
            if (res1.size() != 2) {
                int num_res2 = 0;
                for (auto &&chain2 : _model) for (auto &&res2 : chain2) {
                    if (num_res1 < num_res2 and res2.size() != 2) {
                        MatrixXd dists_pairing = MatrixXd::Zero(2, 2);
                        for (auto m : {0, 1}) for (auto n : {0, 1}) {
                            dists_pairing(m, n) = geometry::distance(res1[m+1], res2[n+1]);
                        }
                        val_t dih_pairing = geometry::dihedral(res1[1], res1[2], res2[2], res2[1]);
                        f(dists_pairing, dih_pairing);
                    }
                    num_res2++;
                }
            }
            num_res1++;
        }
    }

    void assess_two_separate_residues() {
        each_two_separate_residues([&](const auto &dists_pairing, val_t dih_pairing){
            _score_base_pairing += en_base_pairing(dists_pairing, dih_pairing);
        });
    }

    void analyze_two_separate_residues() {
        each_two_separate_residues([](const auto &dists_pairing, val_t dih_pairing){
            std::cout << dists_pairing << std::endl;
            std::cout << dih_pairing << std::endl;
        });
    }

    val_t score_bond_length() {
        return fpl::fold([&](val_t sum, val_t i, val_t j){return sum+i*j;}, 0, _score_bond_length, _weight_bond_length);
    }

    val_t score_bond_angle() {
        return fpl::fold([&](val_t sum, val_t i, val_t j){return sum+i*j;}, 0, _score_bond_angle, _weight_bond_angle);
    }

    val_t score_bond_dihedral() {
        return fpl::fold([&](val_t sum, val_t i, val_t j){return sum+i*j;}, 0, _score_bond_dihedral, _weight_bond_dihedral);
    }

    val_t score_base_pairing() {
        return _score_base_pairing * _weight_base_pairing;
    }

    val_t score_base_stacking() {
        return _score_base_stacking * _weight_base_stacking;
    }

    val_t en_bond_length(val_t standard, val_t actual) {
        return square(actual - standard);
    }

    val_t en_bond_angle(val_t standard, val_t actual) {
        return square(actual - standard);
    }

    val_t en_bond_dihedral(val_t standard, val_t actual) {
        return square(actual - standard);
    }

    template<typename Mat>
    val_t en_base_pairing(Mat &&dists, val_t dih) {
        return 0;
    }

    template<typename Mat>
    val_t en_base_stacking(Mat &&dists, val_t dih) {
        return 0;
    }

};

} // namespace scoring
} // namespace jian

#endif


#ifndef JIAN_SCORING_ASSESSPSB
#define JIAN_SCORING_ASSESSPSB

#include "../pdb/Molecule.h"
#include "../geom/core.h"
#include "../etl.h"

namespace jian {
namespace scoring {

class AssessPSB {
public:
    using val_t = double;
    using diff_t = struct {val_t dist, ang1, ang2, dih;};

    Log log;
    pdb::PSB _model;

    std::array<val_t, 3> _weight_bond_length;
    std::array<val_t, 4> _weight_bond_angle;
    std::array<val_t, 3> _weight_bond_dihedral;
    val_t _weight_base_stacking;
    val_t _weight_base_pairing;

    std::array<val_t, 3> _score_bond_length;
    std::array<val_t, 4> _score_bond_angle;
    std::array<val_t, 3> _score_bond_dihedral;
    val_t _score_base_stacking;
    val_t _score_base_pairing;

    val_t _cutoff = 20;

    AssessPSB() {
        auto init_array = [](auto &arr){for (auto &&i : arr){i = 0;}};

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
        init_score(model);
        assess();
        return score_bond_length() + score_bond_angle() + score_bond_dihedral() + score_base_pairing() + score_base_stacking();
    }

    template<typename ModelType>
    void analyze(const ModelType &model) {
        _model = model;
        analyze();
    }

    template<typename ModelType>
    void init_score(const ModelType &model) {
        _score_bond_length = _score_bond_angle = _score_bond_dihedral = 0;
        _score_base_stacking = _score_base_pairing = 0;
    }

    void assess() {
        assess_residue();
        assess_two_consecutive_residues();
        assess_three_consecutive_residues();
        assess_two_separate_residues();
    }

    void analyze() {
        log(">Residue", '\n');
        analyze_residue();
        log(">Two consecutive residues", '\n');
        analyze_two_consecutive_residues();
        log(">Three consecutive residues", '\n');
        analyze_three_consecutive_residues();
        log(">Two separate residues", '\n');
        analyze_two_separate_residues();
    }

    template<typename Res1, typename Res2>
    val_t dist(Res1 &&res1, Res2 &&res2) {
        return geom::distance(pdb::atom(res1, "S"), pdb::atom(res2, "P"));
    }

    template<typename Fn>
    void each_residue(Fn &&f) {
        for (auto &&chain : _model) for (auto &&res : chain) {
            if (res.size() == 2) continue;
            val_t len1 = geom::distance(res[0], res[1]);
            val_t len2 = geom::distance(res[1], res[2]);
            val_t ang1 = geom::angle(res[0], res[1], res[2]);
            f(res, len1, len2, ang1);
        }
    }

    void assess_residue() {
        each_residue([&](const auto &res, val_t len1, val_t len2, val_t ang1){
            _score_bond_length[0] += en_bond_length(len1);
            _score_bond_length[1] += en_bond_length(len2);
            _score_bond_angle[0] += en_bond_angle(ang1);
        });
    }

    void analyze_residue() {
        each_residue([&](const auto &res, val_t len1, val_t len2, val_t ang1){
            log(res._name, ' ', len1, ' ', len2, ' ', ang1, '\n');
        });
    }

    template<typename Fn>
    void each_two_consecutive_residues(Fn &&f) {
        for (auto &&chain : _model) {
            for (int i = 1; i < chain.size(); i++) {
                if (chain[i-1].size() == 2 or chain[i].size() == 2 or dist(chain[i-1], chain[i]) > 4.2) continue;
                auto len3 = geom::distance(chain[i-1][1], chain[i][0]);
                auto ang2 = geom::angle(chain[i-1][0], chain[i-1][1], chain[i][0]);
                auto ang3 = geom::angle(chain[i-1][2], chain[i-1][1], chain[i][0]);
                auto ang4 = geom::angle(chain[i-1][1], chain[i][0],   chain[i][1]);
                auto dih1 = geom::dihedral(chain[i-1][0], chain[i-1][1], chain[i-1][2], chain[i][0]);
                auto dih2 = geom::dihedral(chain[i-1][0], chain[i-1][1], chain[i][0],   chain[i][1]);
                auto diff_v = diff(chain[i-1], chain[i]);
                f(chain[i-1], chain[i], len3, ang2, ang3, ang4, dih1, dih2, diff_v);
            }
        }
    }

    void assess_two_consecutive_residues() {
        each_two_consecutive_residues([&](const auto &res1, const auto &res2, 
                                          val_t len3, val_t ang2, val_t ang3, val_t ang4, val_t dih1, val_t dih2, 
                                          const auto &diff_v) {
            _score_bond_length[2] += en_bond_length(len3);
            _score_bond_angle[1] += en_bond_angle(ang2);
            _score_bond_angle[2] += en_bond_angle(ang3);
            _score_bond_angle[3] += en_bond_angle(ang4);
            _score_bond_dihedral[0] += en_bond_dihedral(dih1);
            _score_bond_dihedral[1] += en_bond_dihedral(dih2);
        });
    }

    void analyze_two_consecutive_residues() {
        each_two_consecutive_residues([&](const auto &res1, const auto &res2, 
                                          val_t len3, val_t ang2, val_t ang3, val_t ang4, val_t dih1, val_t dih2, 
                                          const auto &diff_v) {
            log(res1._name, ' ', res2._name, ' ', len3, ' ', ang2, ' ', ang3, ' ', ang4, ' ', dih1, ' ', dih2, ' ',
                diff_v.dist, ' ', diff_v.ang1, ' ', diff_v.ang2, ' ', diff_v.dih, '\n');
        });
    }

    template<typename Fn>
    void each_three_consecutive_residues(Fn &&f) {
        for (auto &&chain : _model) for (int i = 2; i < chain.size(); i++) {
            if (chain[i-2].size() == 2 || chain[i-1].size() == 2 || chain[i].size() == 2 ||
                dist(chain[i-2], chain[i-1]) > 4.2 || dist(chain[i-1], chain[i]) > 4.2) continue;
            auto dih3 = geom::dihedral(chain[i-2][1], chain[i-1][0], chain[i-1][1], chain[i][0]);
            f(chain[i-2], chain[i-1], chain[i], dih3);
        }
    }

    void assess_three_consecutive_residues() {
        each_three_consecutive_residues([&](const auto &res1, const auto &res2, const auto &res3, val_t dih3){
            _score_bond_dihedral[2] += en_bond_dihedral(dih3);
        });
    }

    void analyze_three_consecutive_residues() {
        each_three_consecutive_residues([&](const auto &res1, const auto &res2, const auto &res3, val_t dih3){
            log(res1._name, ' ', res2._name, ' ', res3._name, ' ', dih3, '\n');
        });
    }

    template<typename Fn>
    void each_two_separate_residues(Fn &&f) {
        int num_res1 = 0; for (auto &&chain1 : _model) for (auto &&res1 : chain1) {
            int num_res2 = 0; for (auto &&chain2 : _model) for (auto &&res2 : chain2) {
                if (num_res1 < num_res2) {
                    auto diff_v = diff(res1, res2);
                    if (diff_v.dist < _cutoff) f(res1, res2, diff_v);
                }
                num_res2++;
            }
            num_res1++;
        }
    }

    void assess_two_separate_residues() {
        each_two_separate_residues([&](const auto &res1, const auto &res2, const auto &diff_v){
        });
    }

    void analyze_two_separate_residues() {
        each_two_separate_residues([&](const auto &res1, const auto &res2, const auto &diff_v){
            log(res1._name, ' ', res2._name, ' ',
                diff_v.dist, ' ', diff_v.ang1, ' ', diff_v.ang2, ' ', diff_v.dih, '\n');
        });
    }

    val_t score_bond_length() {
        return fold([&](val_t sum, val_t i, val_t j){return sum+i*j;}, 0, _score_bond_length, _weight_bond_length);
    }

    val_t score_bond_angle() {
        return fold([&](val_t sum, val_t i, val_t j){return sum+i*j;}, 0, _score_bond_angle, _weight_bond_angle);
    }

    val_t score_bond_dihedral() {
        return fold([&](val_t sum, val_t i, val_t j){return sum+i*j;}, 0, _score_bond_dihedral, _weight_bond_dihedral);
    }

    val_t score_base_pairing() {
        return _score_base_pairing * _weight_base_pairing;
    }

    val_t score_base_stacking() {
        return _score_base_stacking * _weight_base_stacking;
    }

    val_t en_bond_length(val_t actual) {
    }

    val_t en_bond_angle(val_t actual) {
    }

    val_t en_bond_dihedral(val_t actual) {
    }

    template<typename Mat>
    val_t en_base_pairing(Mat &&dists, val_t dih) {
        return 0;
    }

    template<typename Mat>
    val_t en_base_stacking(Mat &&dists, val_t dih) {
        return 0;
    }

    template<typename T, typename U>
    diff_t diff(T &&t, U &&u) {
        diff_t d; 
        int t1 = t.size() - 2, t2 = t1 + 1, u1 = u.size() - 2, u2 = u1 + 1;
        d.dist = geom::distance(t[t1], u[u1]);
        d.ang1 = geom::angle(t[t2], t[t1], u[u1]);
        d.ang2 = geom::angle(u[u2], u[u1], t[t1]);
        d.dih = geom::dihedral(t[t2], t[t1], u[u1], u[u2]);
        return d;
    }

};

} // namespace scoring
} // namespace jian

#endif


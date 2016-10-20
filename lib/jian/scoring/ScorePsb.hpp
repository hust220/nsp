#pragma once

#include <iostream>
#include <memory>
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "../geom.hpp"
#include "../cg/CGpsb.hpp"
#include "../pdb.hpp"
#include "Score.hpp"

namespace jian {

template<>
class Score<CGpsb> {
public:
    Mat m_freqs_stacking;
    Mat m_freqs_pairing;
    Mati m_counts_stacking;
    Mati m_counts_pairing;
    std::map<std::string, int> m_map {{"A", 0}, {"U", 1}, {"G", 2}, {"C", 3}};
    double m_cutoff = 20;
    double m_bin = 0.2;
    double m_bins;
    std::vector<int> m_indices;
    const Chain *m_chain;

    Score() {
        m_bins = m_cutoff / m_bin;
        read_counts();
        set_freqs();
    }

    void train(const std::string &s) {
        m_counts_stacking = Mati::Zero(144, m_bins);
        m_counts_pairing = Mati::Zero(144, m_bins);
        EACH_SPLIT_LINE(s.c_str(), " ",
            std::cout << F[0] << std::endl;
            Chain chain;
            chain_read_model(chain, F[0]);
            update_counts(CGpsb::chain(chain));
        );
        std::cout << m_counts_stacking << std::endl;
        std::cout << m_counts_pairing << std::endl;
    }

    void update_counts(const Chain &chain) {
        m_chain = &chain;
        set_indices();
        int len = chain.size();
        for (int i = 0; i < len - 1; i++) {
            update_counts_stacking(i, i+1);
            for (int j = i + 2; j < len; j++) {
                update_counts_pairing(i, j);
            }
        }
    }

    void set_indices() {
        int len = m_chain->size();
        m_indices.resize(len);
        for (int i = 0; i < len; i++) {
            m_indices[i] = m_map[m_chain->at(i).name];
        }
    }

    void update_counts_stacking(int n1, int n2) {
        double d;
        int n, a = m_indices[n1], b = m_indices[n2];
        const Residue &r1 = m_chain->at(n1);
        const Residue &r2 = m_chain->at(n2);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                d = geom::distance(r1[i], r2[j]);
                if (d < 20) {
                    n = (a * 3 + i) * 12 + (b * 3 + j);
                    m_counts_stacking(n, int(d / m_bin))++;
                }
            }
        }
    }

    void update_counts_pairing(int n1, int n2) {
        double d;
        int n, a = m_indices[n1], b = m_indices[n2];
        const Residue &r1 = m_chain->at(n1);
        const Residue &r2 = m_chain->at(n2);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                d = geom::distance(r1[i], r2[j]);
                if (d < 20) {
                    n = (a * 3 + i) * 12 + (b * 3 + j);
                    m_counts_pairing(n, int(d / m_bin))++;
                }
            }
        }
    }

    void read_counts() {
        m_counts_stacking = Mati::Zero(144, m_bins);
        m_counts_pairing = Mati::Zero(144, m_bins);
        std::string name = Env::lib() + "/RNA/pars/scoring/score_psb/counts";
        std::ifstream ifile(name.c_str());
        for (int i = 0; i < 144; i++) {
            for (int j = 0; j < m_bins; j++) {
                ifile >> m_counts_stacking(i, j);
            }
        }
        for (int i = 0; i < 144; i++) {
            for (int j = 0; j < m_bins; j++) {
                ifile >> m_counts_pairing(i, j);
            }
        }
        ifile.close();
    }

    void set_freqs() {
        m_freqs_stacking = Mat::Zero(144, m_bins);
        m_freqs_pairing = Mat::Zero(144, m_bins);
        Vec vs = Vec::Zero(m_bins);
        Vec vp = Vec::Zero(m_bins);
        int n, sum_s, sum_p, total_sum_s = 0, total_sum_p = 0;
        for (int i = 0; i < 144; i++) {
            sum_s = 0;
            sum_p = 0;
            for (int j = 0; j < m_bins; j++) {
                n = m_counts_stacking(i, j);
                sum_s += n;
                total_sum_s += n;
                vs[j] += n;
                n = m_counts_pairing(i, j);
                sum_p += n;
                total_sum_p += n;
                vp[j] += n;
            }
            for (int j = 0; j < m_bins; j++) {
                m_freqs_stacking(i, j) = m_counts_stacking(i, j) / double(sum_s);
                m_freqs_pairing(i, j) = m_counts_pairing(i, j) / double(sum_p);
            }
        }
        for (int i = 0; i < m_bins; i++) {
            vs[i] /= double(total_sum_s);
            vp[i] /= double(total_sum_p);
        }
        for (int i = 0; i < 144; i++) {
            for (int j = 0; j < m_bins; j++) {
                if (vs[j] > 0.0003) {
                    m_freqs_stacking(i, j) /= vs[j];
                } else {
                    m_freqs_stacking(i, j) = 0;
                }
                if (vp[j] > 0.0003) {
                    m_freqs_pairing(i, j) /= vp[j];
                } else {
                    m_freqs_pairing(i, j) = 0;
                }
            }
        }
//        std::cout << m_freqs_stacking << std::endl;
    }

    double en_stacking(const Residue &r1, const Residue &r2) {
        double e = 0, d, f;
        int a, b, t1 = m_map[r1.name], t2 = m_map[r2.name];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                d = geom::distance(r1[i], r2[j]);
                if (d < 20) {
                    a = t1 * 3 + i;
                    b = t2 * 3 + j;
                    f = m_freqs_stacking(a * 12 + b, int(d/ m_bin));
                    if (f != 0) {
                        e += -std::log(f);
                    } else {
                        e += 0;
                    }
                }
            }
        }
        return e;
    }

    double en_pairing(const Residue &r1, const Residue &r2) {
        double e = 0, d, f;
        int a, b, t1 = m_map[r1.name], t2 = m_map[r2.name];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                d = geom::distance(r1[i], r2[j]);
                if (d < 20) {
                    a = t1 * 3 + i;
                    b = t2 * 3 + j;
                    f = m_freqs_pairing(a * 12 + b, int(d/ m_bin));
                    if (f != 0) {
                        e += -std::log(f);
                    } else {
                        e += 0;
                    }
                }
            }
        }
        return e;
    }

    void print_freqs() {
        std::cout << m_freqs_stacking << std::endl;
        std::cout << m_freqs_pairing << std::endl;
    }

};

} // namespace jian


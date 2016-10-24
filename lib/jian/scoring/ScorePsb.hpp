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
class Score<CGpsb> : public ScoreBase {
public:
	using cg_t = CGpsb;

    Mat m_freqs_stacking;
    Mat m_freqs_pairing;
    Mati m_counts_stacking;
    Mati m_counts_pairing;
    std::map<std::string, int> m_map {{"A", 0}, {"U", 1}, {"G", 2}, {"C", 3}};
    double m_cutoff = 20;
    double m_bin = 0.2;
    int m_bins;
    std::vector<int> m_indices;
    const Chain *m_chain;

	virtual void init();

	virtual void run(const Chain & chain);

	virtual void train(const Chain & chain);

	virtual void print_counts(std::ostream &) const;

	//void update_counts(const Chain &chain);

	void set_indices();

	void update_counts_stacking(int n1, int n2);

	void update_counts_pairing(int n1, int n2);

	void read_counts();

	void set_freqs();

	double en_stacking(const Residue &r1, const Residue &r2);

	double en_pairing(const Residue &r1, const Residue &r2);

	void print_freqs_stacking(std::ostream &) const;

	void print_freqs_pairing(std::ostream &) const;

	void print_counts_stacking(std::ostream &) const;

	void print_counts_pairing(std::ostream &) const;

};

} // namespace jian


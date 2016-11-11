#pragma once

#include <iostream>
#include <memory>
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "../geom.hpp"
#include "../cg/CG6p.hpp"
#include "../pdb.hpp"
#include "ScoreBase.hpp"
#include "ParBp.hpp"

namespace jian {

	class Score : public ScoreBase {
	public:
		std::shared_ptr<ParBp> parbp;
		Mat m_freqs_stacking;
		Mat m_freqs_pairing;
		Mat m_freqs_wc;
		Mat m_freqs_nwc;
		Mati m_counts_stacking;
		Mati m_counts_pairing;
		Mati m_counts_wc;
		Mati m_counts_nwc;
		std::array<std::array<Mat3i, 4>, 4> m_counts_bp;
		std::array<std::array<Mat3, 4>, 4> m_freqs_bp;
		std::array<std::array<Mat3i, 4>, 4> m_counts_st53;
		std::array<std::array<Mat3, 4>, 4> m_freqs_st53;
		std::array<std::array<Mat3i, 4>, 4> m_counts_st35;
		std::array<std::array<Mat3, 4>, 4> m_freqs_st35;
		std::array<std::array<double, 4>, 4> m_weights_bp;
		std::array<std::array<double, 4>, 4> m_weights_st53;
		std::array<std::array<double, 4>, 4> m_weights_st35;
		std::vector<double> m_bond_len_std;
		std::vector<double> m_bond_angle_std;
		std::vector<double> m_bond_dihedral_std;
		double m_cutoff_stacking;
		double m_cutoff_pairing;
		std::map<std::string, int> m_map{ { "A", 0 },{ "U", 1 },{ "G", 2 },{ "C", 3 } };
		double m_cutoff;
		double m_bin;
		int m_bins;
		int m_res_size;
		int m_num_types;
		int m_rows;
		std::vector<int> m_indices;
		const Chain *m_chain;

		virtual void init();

		virtual void run(const Chain & chain);

		virtual void train(const Chain & chain);

		virtual void print_counts(std::ostream &) const;

		virtual void print_freqs(std::ostream &) const;

		virtual double en_stacking(const Residue &r1, const Residue &r2);

		virtual double en_pairing(const Residue &r1, const Residue &r2);

		virtual bool in_base(int);

		virtual ScoreBase &en_bp(const Residue &r1, const Residue &r2);

		virtual double en_len(const Chain &c, int beg);

		virtual double en_ang(const Chain &c, int beg);

		virtual double en_dih(const Chain &c, int beg);

		virtual double en_crash(const Residue &r1, const Residue &r2);

		void set_indices();

		void update_counts_stacking(int n1, int n2);

		void update_counts_pairing(int n1, int n2);

		void read_counts();

		void read_pars();

		void set_freqs();

	};

} // namespace jian


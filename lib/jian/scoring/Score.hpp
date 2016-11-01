#pragma once

#include <iostream>
#include <memory>
#include "../utils/Env.hpp"
#include "../utils/file.hpp"
#include "../geom.hpp"
#include "../cg/CG6p.hpp"
#include "../pdb.hpp"
#include "ScoreBase.hpp"

namespace jian {

	class Score : public ScoreBase {
	public:
		Mat m_freqs_stacking;
		Mat m_freqs_pairing;
		Mati m_counts_stacking;
		Mati m_counts_pairing;
		double m_bond_len_std;
		std::vector<double> m_bond_angle_std;
		double m_bond_dihedral_std;
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

		virtual double en_len(const Residue &r1, const Residue &r2);

		virtual double en_ang(const Residue &r1, const Residue &r2, const Residue &r3);

		virtual double en_dih(const Residue &r1, const Residue &r2, const Residue &r3, const Residue &r4);

		virtual double en_crash(const Residue &r1, const Residue &r2);

		void set_indices();

		void update_counts_stacking(int n1, int n2);

		void update_counts_pairing(int n1, int n2);

		void read_counts(Mati & mat, std::string path, int cutoff);

		void read_counts();

		void read_pars();

		void set_freqs(Mat & m, const Mati & c);

		void set_freqs();

		void print_counts(std::ostream &, const Mati &) const;

		void print_freqs(std::ostream &, const Mat &) const;

	};

} // namespace jian


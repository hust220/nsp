#pragma once

#include <array>
#include <string>

namespace jian {

	class Chain;

	class DistAnal {
	public:
		using Point = std::array<double, 4>;

		double score = 0;
		double interval;
		double cutoff;
		int bins;
		double penalty = 0;

		int len = 0;
		std::vector<int> num, type, size_nt;
		std::vector<Point *> m_coords;

		std::vector<double> m_freqs;
		std::vector<int> m_counts;

		DistAnal & init(double = 0.5, double = 20);
		~DistAnal();

		void read_mol(const Chain &);
		DistAnal & train(const Chain &);
		DistAnal & run(const Chain &);

		void read_freqs(std::string);
		void print_freqs(std::ostream &) const;
		void print_counts(std::ostream &) const;
		bool in_base(int type);
		int res_type(std::string name);
		int atom_type(const Residue &r, int k);
		double en_stacking(const Residue &r1, const Residue &r2);
		double en_pairing(const Residue &r1, const Residue &r2);
	};

}



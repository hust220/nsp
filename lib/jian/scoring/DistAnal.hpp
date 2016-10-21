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
		unsigned int bins;
		double penalty = 0;

		unsigned int len = 0;
		std::vector<unsigned int> num, type, size_nt;
		std::vector<Point *> m_coords;

		std::vector<double> m_freqs;
		std::vector<unsigned int> m_counts;

		DistAnal & init(double = 0.5, int = 20);
		~DistAnal();

		void read_mol(const Chain &);
		DistAnal & train(const Chain &);
		DistAnal & run(const Chain &);

		void read_freqs(std::string);
		void print_freqs() const;
		void print_counts() const;
		bool in_base(unsigned int type);
		int res_type(std::string name);

	};

}


#pragma once

#include "../cg.hpp"
#include "../geom.hpp"
#include "../pp.hpp"
#include "../utils/log.hpp"
#include "../utils/Env.hpp"

#define MEM_EN_EN len, ang, dih, crash, cons, vdw, stacking, pairing
#define DEF_MEM_EN_EN(a) double a = 0;
#define SUM_MEM_EN_EN(a) + a
#define PRINT_MEM_EN(a) << a << PP_STRING3((a)) << ' '

#define MEM_EN \
	bond_length_weight,\
    bond_angle_weight, bond_angle_std,\
    bond_dihedral_weight, bond_dihedral_std, \
    pairing_weight, stacking_weight, \
    constraints_weight, crash_weight, vdw_weight

#define DEF_MEM_EN(a) double PP_CAT3(m_, a);

namespace jian {

	template<typename T>
	class En;

	template<>
	class En<CGpsb> {
	public:
		struct en_t {
			JN_MAP(DEF_MEM_EN_EN, MEM_EN_EN)
			double sum() const { return 0 JN_MAP(SUM_MEM_EN_EN, MEM_EN_EN); }
			void print() const { LOG << sum() << "(total) " JN_MAP(PRINT_MEM_EN, MEM_EN_EN) << std::endl; }
		};

		// Define members
		JN_MAP(DEF_MEM_EN, MEM_EN)

		Mat m_stacking_pars;
		Mat m_pairing_pars;

		void init();

		void read_stacking_pairing_parameters();

		void read_parameters();

		void print_parameters() const;

		double en_stacking(const Mat &m, int iTypeRes1, int iTypeRes2) const;

		double en_pairing(const Mat & m, int iTypeRes1, int iTypeRes2) const;

		double en_crash(const Residue & a, const Residue & b, Mat & arr) const;

	};

	using EnPsb = En<CGpsb>;

} // namespace jian


#pragma once

#include "McsmBase.hpp"
#include "../cg.hpp"
#include "../scoring/Score.hpp"

#define MEM_EN_MCPSB len, ang, dih, crash, cons, vdw, stacking, pairing, wc, nwc
#define DEF_MEM_EN_MCPSB(a) double a = 0;
#define SUM_MEM_EN_MCPSB(a) + a
#define PRINT_MEM_EN_MCPSB(a) << a << PP_STRING3((a)) << ' '

namespace jian {
	namespace nuc3d {
		namespace mc {

			class MCSM : public MCBase {
			public:
				struct en_t {
					JN_MAP(DEF_MEM_EN_MCPSB, MEM_EN_MCPSB);
					double sum() const { return 0 JN_MAP(SUM_MEM_EN_MCPSB, MEM_EN_MCPSB); }
					void print() const { LOG << sum() << "(total) " JN_MAP(PRINT_MEM_EN_MCPSB, MEM_EN_MCPSB) << std::endl; }
				};

				std::vector<int> m_indices;
				ScoreBase *m_scorer;

				MCSM() = default;

				void init(const Par &par);

				void set_indices();

				void print_pairing();

				virtual double mc_partial_energy();

				void mc_total_energy(en_t &e);

				void mc_partial_energy_crash(en_t &e);
				void mc_total_energy_crash(en_t &e);

				void mc_partial_energy_bond(en_t &e);
				void mc_total_energy_bond(en_t &e);

				void mc_partial_energy_angle(en_t &e);
				void mc_total_energy_angle(en_t &e);

				void mc_partial_energy_dihedral(en_t &e);
				void mc_total_energy_dihedral(en_t &e);

				void mc_partial_energy_constraints(en_t &e);
				void mc_total_energy_constraints(en_t &e);

				virtual double dist_two_res(const Residue &r1, const Residue &r2) const;

				double total_energy();

				virtual void write_en();

				virtual std::string file_parameters() const;

				virtual void finish_run();

				virtual bool is_selected(const int &i) const = 0;

				virtual Vec rotating_center() const = 0;

			};

		} // namespace mc
	} // namespace nuc3d
} // namespace jian


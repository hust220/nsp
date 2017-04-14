#pragma once

#include "McsmBase.hpp"
#include "../cg.hpp"
#include "../scoring/Score.hpp"

#define MEM_EN_MCPSB len, ang, dih, crash, cons, vdw, stacking, pairing, rg
#define DEF_MEM_EN_MCPSB(a) double a = 0;
#define INIT_MEM_EN_MCPSB(a) a = 0;
#define SUM_MEM_EN_MCPSB(a) + a
#define PRINT_MEM_EN_MCPSB(a) << a << PP_STRING3((a)) << ' '

BEGIN_JN

class MCSM : public MCBase {
    public:
        struct en_t {
            JN_MAP(DEF_MEM_EN_MCPSB, MEM_EN_MCPSB);
            en_t() {
                JN_MAP(INIT_MEM_EN_MCPSB, MEM_EN_MCPSB);
            }
            double sum() const { return 0 JN_MAP(SUM_MEM_EN_MCPSB, MEM_EN_MCPSB); }
            void print() const { LOG << sum() << "(total) " JN_MAP(PRINT_MEM_EN_MCPSB, MEM_EN_MCPSB) << std::endl; }
        };

        std::vector<int> m_indices;
        Score *m_scorer;
        std::vector<int> m_bps;
        std::array<std::array<Mat, 4>, 4> m_bp_min_distances_table;
        std::array<std::array<Mat, 4>, 4> m_bp_max_distances_table;

        MCSM() = default;

        void init(const Par &par);

        void set_indices();

        void print_pairing();

        virtual double mc_partial_energy();

        virtual double mc_total_energy();

        void mc_energy_crash(en_t &e, bool is_total);
        void mc_energy_bond(en_t &e, bool is_total);
        void mc_energy_angle(en_t &e, bool is_total);
        void mc_energy_dihedral(en_t &e, bool is_total);
        void mc_energy_constraints(en_t &e, bool is_total);
        void mc_energy_rg(en_t &e, bool is_total);

//        virtual double dist_two_res(const Residue &r1, const Residue &r2) const;

        void set_total_energy(en_t &e);

        virtual void write_en();

        virtual void finish_run();

        virtual void mc_select() = 0;

        virtual void mc_rollback() = 0;

        virtual void mc_backup() = 0;

        virtual void mc_sample() = 0;

        virtual bool is_selected(const int &i) const = 0;

//        virtual Vec rotating_center() const = 0;

};

END_JN


#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "rtsp.hpp"
#include "rtsp_templ_rec.hpp"
#include "cg.hpp"
#include "cg_res_frags.hpp"
#include "env.hpp"
#include "file.hpp"
#include "mc.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "rss.hpp"
#include "pp.hpp"
#include "rtsp_build_chain.hpp"
#include "rtsp_build_line.hpp"
#include "score_res_conf.hpp"
#include "score_frag_conf.hpp"
#include "bmmc_mvel.hpp"

#define JN_MCXP_PARS1 \
    heat_steps, cool_steps, cycle_steps, write_steps, heat_rate, dec_rate,\
init_tempr, queue, lowest_rate, lowest_en, lowest_tempr, highest_tempr

#define JN_MCXP_PARS2 \
    bond_length_weight, bond_angle_weight, bond_angle_std, bond_dihedral_weight, bond_dihedral_std, \
pairing_weight, wc_weight, nwc_weight, stacking_weight,  constraints_weight, crash_weight, rg_weight, \
contacts_weight, vdw_weight, max_shift

#define JN_MCXP_DEF_PAR(a) Num PP_CAT3(_mc_, a);

BEGIN_JN

// MCBase: Base of Monte Carlo Molecule Simulation
class MCBase : public TSP, public MC {
    public:
        using SpaceItem = Atom;
        using SpaceVal = List<Int>;
        using Space = Map<Int, Map<Int, Map<Int, SpaceVal>>>;
        using SpaceVals = Vector<SpaceVal *>;

        enum SampleMode {
            SAMPLE_SSE,
            SAMPLE_TREE
        };

        struct atom_acc_t {
            int n_res;
            int n_atom;
        };

        struct distance_constraints_t {
            atom_acc_t atom1, atom2;
            double min, max;
        };

        // grow
        Bool m_sample_tree;
        Bool m_del_pk;
        Bool m_pk_ahead;
        Bool m_save_ss;
        Bool m_grow_mode;
        Int m_grow_steps;
        Int m_grow_length;

        Deque<distance_constraints_t> m_distance_constraints;

        ResConf::MapConfs m_res_confs;
        FragConf<3>::MapConfs m_frag_confs;
        SampleMode m_sample_mode;
        Bool m_cal_en_constraints;
        Bool m_will_write_traj;
        Space m_space;
        SpaceVals m_item_space;
        Int m_box;
        Num m_box_size;
        Deque<Atom> _moved_atoms;
        Str m_traj;
        Deque<Num> m_scores;
        Int mc_num_writing = 1;
        Vector<Int> m_continuous_pts;
        Vector<Int> m_ang_pts;
        Vector<Int> m_dih_pts;
        Vector<Int> m_brk_pts;
        Str m_par_file;
        MolWriter m_writer;
        Num m_max_angle;
        Str m_init_sfile;
        Str m_alignfile;
        Deque<Array<Int, 2>> m_align;

        //Deque<MvEl *> m_mvels;
        Vector<MvEl *> m_base_mvels;
        //MvEl *m_selected_mvel = NULL;
        List<Vector<Bool>> m_fixed_areas;

        JN_MAP(JN_MCXP_DEF_PAR, JN_MCXP_PARS2);

        MCBase() = default;

        void init(const Par &par);

        void set_align();

        void align_to_fixed_areas();

        void restore_raw();

        void set_traj_name();

        void validate_constraints();

        void read_parameters();

        void set_parameters(const Par &par);

        void print_parameters();

        void set_continuous_pts();

        void mc_write();

        void write_traj();

        void init_space();

        int space_index(double n) const;

        SpaceItem &item(int i);

        SpaceVal &space_val(int i);

        void space_update_item(int i);

        void run();

        void score();

        void transform();

        void print_final_constraints();

        virtual void mc_next_step();

        virtual double mc_partial_energy() = 0;

        virtual double mc_total_energy() = 0;

        virtual void write_en() = 0;

        virtual void before_run();

        virtual void finish_run();

        virtual void save_fixed_ranges();

        virtual void restore_fixed_ranges();

        virtual void mc_select() = 0;

        virtual void mc_rollback() = 0;

        virtual void mc_backup() = 0;

        virtual void mc_sample() = 0;

        virtual bool is_selected(const int &i) const = 0;

};

END_JN


#pragma once

#include <Eigen/Dense>
#include <map>
#include <deque>
#include <vector>
#include <numeric>
#include "g4dna_module.hpp"
#include "pdb.hpp"
#include "geom.hpp"
#include "bmmc.hpp"
#include "factory.hpp"
#include "par.hpp"
#include "env.hpp"

namespace jian {
namespace qhmc {

    using fac_t = Factory<Module::cons_t>;

    class QHMC : public MCSM {
        public:
            using res_t = struct {char seq; char ss; int num;};
            using res_list_t = std::deque<res_t>;
            using indices_t = std::deque<int>;
            using related_residues_t = std::deque<indices_t *>;

            Tree m_tree;
            std::deque<Module *> m_modules;
            int m_selected_index;
            related_residues_t m_related_residues,
                               m_unrelated_residues,
                               m_all_indices;
            std::array<bool, 4> m_directions; // true: parallel, false: antiparallel
            std::array<int, 4> m_arrangement;
            std::map<indices_t *, Chain> m_fixed_ranges;

            QHMC(): MCSM() {};

            void init(const Par &par);

            ~QHMC();

            void print_modules();

            void set_modules();

            void set_res_list(res_list_t &res_list);

            void ss_to_tree();

            void build_initial_scaffold();

            Chain build_helix(int len);

            Chain load_quadruple_helix(int n);

            void set_coords_residue(Mat &c1, int m, const Residue &r);

            Chain connect_quadruple_helix(Chain &c1, Chain &c2);

            void shrink_to_fit(const Chain &c);

            void tuples_from_ss(Tuples &tuples);

            template<typename T, typename U>
                bool adjacent(T &&t1, U &&t2) {
                    return abs(t1[0] - t2[0]) == 1 &&
                        abs(t1[1] - t2[1]) == 1 &&
                        abs(t1[2] - t2[2]) == 1 &&
                        abs(t1[3] - t2[3]) == 1;
                }

            void tuples_to_tree(Tuples tuples);

            void print_tuple(const Tuple &tuple);

            void print_helix(const Tuples &helix);

            void print_tree();

            // mc-related functions

            void set_unrelated_residues();

            void print_related_residues(const related_residues_t &r);

            void set_related_residues();

            void set_related_and_unrelated_residues();

            virtual void before_run();

            virtual void mc_sample();

            virtual void mc_rollback();

            virtual void mc_backup();

            virtual void mc_select();

            virtual bool is_selected(const int &i) const;

            virtual Vec rotating_center() const;

            void save_helix();

            virtual void save_fixed_ranges();

            void restore_helix(indices_t * indices);

            virtual void restore_fixed_ranges();

            virtual void read_ss();

    };

} // namespace qhmc

void tsp_g4(const Par &par);

}



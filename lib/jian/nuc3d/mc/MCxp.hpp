#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../../mc/MC.hpp"
#include "../JobPredict3D.hpp"
#include "../transform.hpp"
#include "../../utils/Env.hpp"
#include "../../pp.hpp"
#include "../../utils/file.hpp"
#include "../../pdb.hpp"
#include "../../geom.hpp"

namespace jian {
namespace nuc3d {
namespace mc {

template<typename CG_T>
class MCxp : public JobPredict3D, public MC {
public:    
    using cg_t = CG_T;

    using item_t = Atom;
    using space_val_t = std::list<int>;
    using space_t = std::map<int, std::map<int, std::map<int, space_val_t>>>;
    using item_space_t = std::vector<space_val_t *>;

    space_t m_space;
    item_space_t m_item_space;
    int m_box = 3;
    double m_box_size = 6;
    std::deque<Atom> _moved_atoms;
    Chain _pred_chain;
    std::string m_traj;
    std::string m_out_pdb;
    int mc_num_writing = 1;
    std::vector<int> m_continuous_pts;
    std::vector<int> m_ang_pts;
    std::vector<int> m_dih_pts;
    std::vector<int> m_brk_pts;
    std::string m_par_file;

    #define JN_MCXP_PARS1 heat_steps, cool_steps, cycle_steps, write_steps, heat_rate, dec_rate, init_tempr
    #define JN_MCXP_PARS2 bond_length_weight, bond_angle_weight, bond_angle_std, bond_dihedral_weight, bond_dihedral_std, \
                        constraints_weight, crash_weight, pairing_weight, stacking_weight, vdw_weight, max_shift

    #define JN_MCXP_DEF_PAR(a) double PP_CAT(_mc_, a);
    JN_MAP(JN_MCXP_DEF_PAR, JN_MCXP_PARS2)

    MCxp(const Par &par) : MC(), JobPredict3D(par) {
        LOG << "# Set the file of trajectory..." << std::endl;
        par.set(m_traj, "traj");

        LOG << "# Set continuous points..." << std::endl;
        set_continuous_pts();

        LOG << "# Print continuous, angel, dihedral points..." << std::endl;
        for (auto && i : m_continuous_pts) LOG << i << ' '; LOG << std::endl;
        for (auto && i : m_ang_pts) LOG << i << ' '; LOG << std::endl;
        for (auto && i : m_dih_pts) LOG << i << ' '; LOG << std::endl;

        LOG << "# Read initial structure" << std::endl;
        if (par.has("pdb")) _pred_chain = residues_from_file(par["pdb"][0]);

    }

    void read_parameters() {
        LOG << "Reading parameters of " << m_par_file << "..." << std::endl;
        Par temp_par(Env::lib() + "/RNA/pars/nuc3d/mc/" + m_par_file + ".par");
        #define JN_MCXP_TEMPPAR_SET(a) temp_par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
        JN_MAP(JN_MCXP_TEMPPAR_SET, JN_MCXP_PARS1, JN_MCXP_PARS2)
    }

    void set_parameters(const Par &par) {
        #define JN_MCXP_PAR_SET(a) par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
        JN_MAP(JN_MCXP_PAR_SET, JN_MCXP_PARS1, JN_MCXP_PARS2)
    }

    void print_parameters() {
        #define JN_MCXP_TEMP(a) LOG << PP_STRING3(PP_CAT(mc_, a)) << ' ' << PP_CAT(_mc_, a) << std::endl;
        JN_MAP(JN_MCXP_TEMP, JN_MCXP_PARS1, JN_MCXP_PARS2)
    }

    void set_continuous_pts() {
        int i = 0;
        for (auto && c : _ss) {
            if (c != '&') {
                i++;
            } else {
                m_brk_pts.push_back(i-1);
            }
        }
        for (int i = 0; i < _seq.size() - 1; i++) {
            if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end()) {
                m_continuous_pts.push_back(i);
            }
        }
        for (int i = 0; i < _seq.size() - 2; i++) {
            if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() && 
                std::find(m_brk_pts.begin(), m_brk_pts.end(), i+1) == m_brk_pts.end()) {
                m_ang_pts.push_back(i);
            }
        }
        for (int i = 0; i < _seq.size() - 3; i++) {
            if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
                std::find(m_brk_pts.begin(), m_brk_pts.end(), i+1) == m_brk_pts.end() &&
                std::find(m_brk_pts.begin(), m_brk_pts.end(), i+2) == m_brk_pts.end()) {
                m_dih_pts.push_back(i);
            }
        }
    }

    void mc_write() {
        write_traj();
        write_en();
        mc_num_writing++;
    }

    void write_traj() {
        if (!(m_traj.empty())) {
            if (mc_num_writing == 1) file::clear_file(m_traj);
            append_chain_to_file(_pred_chain, m_traj, mc_num_writing);
        }
    }

    void mc_sample() {
        backup();
        if (rand() < 0.5) {
            // translate
            int index = int(rand() * 3);
            double dist = (rand() - 0.5) * 2 * _mc_max_shift;
            for (int i = 0; i < _seq.size(); i++) {
                if (is_selected(i)) {
                    for (auto && atom : _pred_chain[i]) {
                        atom[index] += dist;
                        space_update_item(i);
                    }
                }
            }
        } else {
            // rotate
            int index = int(rand() * 3);
            double dih = (rand() - 0.5) * PI / 6;
            auto &&rot = geom::rot_mat(index, dih);
            auto &&origin = rotating_center();
            for (int i = 0; i < _seq.size(); i++) {
                if (is_selected(i)) {
                    for (auto && atom : _pred_chain[i]) {
                        geom::rotate(atom, origin, rot);
                    }
                    space_update_item(i);
                }
            }
        }
    }

    void mc_back() {
        for (int i = 0; i < _seq.size(); i++) {
            if (is_selected(i)) {
                for (auto && atom : _pred_chain[i]) {
                    atom = _moved_atoms.front();
                    _moved_atoms.pop_front();
                }
                space_update_item(i);
            }
        }
    }

    void backup() {
        _moved_atoms.clear();
        for (int i = 0; i < _seq.size(); i++) {
            if (is_selected(i)) {
                for (auto && atom : _pred_chain[i]) {
                    _moved_atoms.push_back(atom);
                }
            }
        }
    }

    void init_space() {
        m_item_space.resize(_seq.size());
        for (int i = 0; i < _seq.size(); i++) {
            space_val_t &s = space_val(i);
            s.push_back(i);
            m_item_space[i] = &s;
        }
    }

    int space_index(double n) const {
        return (n+1000)/m_box_size;
    }

    item_t &item(int i) {
        return _pred_chain[i][0];
    }

    space_val_t &space_val(int i) {
        item_t &a = item(i);
        return m_space[space_index(a[0])][space_index(a[1])][space_index(a[2])];
    }

    void space_update_item(int i) {
        space_val_t &n = space_val(i);
        space_val_t &o = *(m_item_space[i]);
        if (&o != &n) {
            o.erase(std::find(o.begin(), o.end(), i));
            m_item_space[i] = &n;
            n.push_back(i);
        }
    }

    void run_job() {
        LOG << "# Displaying starting information..." << std::endl;
        display_start_information();

        run();

        LOG << "# Writing final model..." << std::endl;
        write_final_model();

        LOG << "# Displaying ending information..." << std::endl;
        display_end_information();
    }

    void write_final_model() {
        std::ostringstream stream;
        stream << _name << ".mc.pdb";
        m_out_pdb = stream.str();
        _par->set(m_out_pdb, "out", "out_pdb");
        residues_to_file(_pred_chain, m_out_pdb);
    }

    void run() {
        LOG << "# Read parameters..." << std::endl;
        m_par_file = file_parameters();
        _par->set(m_par_file, "par_file");
        read_parameters();

        LOG << "# Set parameters..." << std::endl;
        set_parameters(*_par);

        LOG << "# Print parameters..." << std::endl;
        print_parameters();

        LOG << "# Initializing running..." << std::endl;
        init_run();

        LOG << "# Carrying on CG processing with the Chain..." << std::endl;
        _pred_chain = cg_t::chain(_pred_chain);

        LOG << "# Init space..." << std::endl;
        init_space();

        LOG << "# MC..." << std::endl;
        mc_heat();
        mc_cool();

        LOG << "# Print Constraints..." << std::endl;
        print_constraints();

        LOG << "# Finishing running..." << std::endl;
        finish_run();

        LOG << "# Coarsed Grained To All Atom..." << std::endl;
        cg_to_aa(chain_to_coords());

        LOG << "# Transform..." << std::endl;
        this->transform();

    }

    void transform() {
        Model m;
        m.push_back(_pred_chain);
        _pred_chain = std::move(jian::transform(m, _seq, _type)[0]);
    }

    Mat chain_to_coords() {
        int n_atoms = num_atoms(_pred_chain);
        Mat c(n_atoms, 3);
        int n_atom = 0;
        for (int i = 0; i < _pred_chain.size(); i++) {
            for (auto && atom : _pred_chain[i]) {
                for (int j = 0; j < 3; j++) {
                    c(n_atom, j) = atom[j];
                }
                n_atom++;
            }
        }
        return c;
    }

    void print_constraints() {
        double d;
        int i, j;
        for (auto && ct : _constraints.distances) {
            i = ct.key[0];
            j = ct.key[1];
            d = dist_two_res(_pred_chain[i], _pred_chain[j]);
            LOG << i << ' ' << j << " value:" << ct.value << " weight:" << ct.weight << " dist:" << d << std::endl;
        }
    }

    void cg_to_aa(const Mat &c) {
        _pred_chain = cg_t::aa(c, 0, c.rows()-1);
    }

    virtual double mc_partial_energy() = 0;
    virtual double dist_two_res(const Residue &, const Residue &) const = 0;
    virtual void write_en() = 0;

    virtual void init_run() {}
    virtual void finish_run() {}

    virtual std::string file_parameters() const {
        return "3drna";
    }

    virtual void mc_select() = 0;
    virtual bool is_selected(const int &i) const = 0;
    virtual Vec rotating_center() const = 0;

};

} // namespace mc
} // namespace nuc3d
} // namespace jian


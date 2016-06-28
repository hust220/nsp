#pragma once

#include <iostream>
#include <set>
#include <memory>
#include <sstream>
#include "../../mc/MC.hpp"
#include "../JobPredict3D.hpp"
#include "../C2A.hpp"
#include "../transform.hpp"
#include "../../utils/Env.hpp"
#include "../../pp.hpp"
#include "../../utils/file.hpp"
#include "../../pdb.hpp"
#include "../../geom.hpp"

namespace jian {
namespace nuc3d {
namespace mc {

class MCxp : public JobPredict3D, public MC {
public:    
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

    #define JN_MC_PARS1 heat_steps, cool_steps, cycle_steps, write_steps, heat_rate, dec_rate, init_tempr
    #define JN_MC_PARS2 bond_length_weight, bond_angle_weight, bond_angle_std, bond_dihedral_weight, bond_dihedral_std, \
                        constraints_weight, crash_weight, pairing_weight, stacking_weight, vdw_weight, max_shift

    #define JN_MCPSB_DEF_PAR(a) double PP_CAT(_mc_, a);
    JN_MAP(JN_MCPSB_DEF_PAR, JN_MC_PARS2)

    MCxp(const Par &par) : MC(), JobPredict3D(par) {
        std::cout << "# Read parameters..." << std::endl;
        read_parameters();

        std::cout << "# Set parameters..." << std::endl;
        #define JN_MCPSB_PAR_SET(a) par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
        JN_MAP(JN_MCPSB_PAR_SET, JN_MC_PARS1, JN_MC_PARS2)

        std::cout << "# Print parameters..." << std::endl;
        print_parameters();

    }

    void read_parameters() {
        Par temp_par(Env::lib() + "/RNA/pars/mc/" + m_cmd + ".par");
        #define JN_MCPSB_TEMPPAR_SET(a) temp_par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
        JN_MAP(JN_MCPSB_TEMPPAR_SET, JN_MC_PARS1, JN_MC_PARS2)
    }

    void print_parameters() {
        #define JN_MCPSB_TEMP(a) std::cout << PP_STRING3(PP_CAT(mc_, a)) << ' ' << PP_CAT(_mc_, a) << std::endl;
        JN_MAP(JN_MCPSB_TEMP, JN_MC_PARS1, JN_MC_PARS2)
    }

    void mc_write() {
        static int n = 1;
        std::ostringstream stream;
        stream << _name << "." << m_cmd << "." << m_seed << ".traj.pdb";
        std::string name = stream.str();

        if (n == 1) file::clear_file(name);
        append_chain_to_file(_pred_chain, name, n);
        write_en();
        n++;
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

    void run() {
        std::cout << "# Init chain..." << std::endl;
        init_chain();

        std::cout << "# Init space..." << std::endl;
        init_space();

        std::cout << "# MC..." << std::endl;
        mc_heat();
        mc_cool();

        std::cout << "# Print Constraints..." << std::endl;
        print_constraints();

        std::cout << "# Coarsed Grained To All Atom..." << std::endl;
        cg_to_aa(chain_to_coords());

        std::cout << "# Transform..." << std::endl;
        this->transform();

        std::cout << "# Writing to file..." << std::endl;
        std::ostringstream stream;
        stream << _name << "." << m_cmd << "." << m_seed << ".pdb";
        residues_to_file(_pred_chain, stream.str());

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
            std::cout << i << ' ' << j << " value:" << ct.value << " weight:" << ct.weight << " dist:" << d << std::endl;
        }
    }

    virtual void cg_to_aa(const Mat &) = 0;
    virtual double mc_partial_energy() = 0;
    virtual double dist_two_res(const Residue &, const Residue &) const = 0;
    virtual void write_en() = 0;

    virtual void init_chain() = 0;
    virtual void mc_select() = 0;
    virtual bool is_selected(const int &i) const = 0;
    virtual Vec rotating_center() const = 0;

};

} // namespace mc
} // namespace nuc3d
} // namespace jian


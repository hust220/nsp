#include <functional>
#include "../nuc3d/Assemble.hpp"
#include "../rtsp/mutate.hpp"
#include "McsmBase.hpp"

#define JN_MCXP_TEMPPAR_SET(a) temp_par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
#define JN_MCXP_PAR_SET(a) par.set(PP_CAT(_mc_, a), PP_STRING3(PP_CAT(mc_, a)));
#define JN_MCXP_TEMP(a) log << PP_STRING3(PP_CAT(mc_, a)) << ' ' << PP_CAT(_mc_, a) << std::endl;

BEGIN_JN

void MCBase::set_traj_name() {
#ifdef JN_PARA
    //m_traj = (mpi_size() == 1 ? to_str(_name, ".traj.pdb") : to_str(_name, ".", mpi_rank() + 1, ".traj.pdb"));
    m_traj = to_str(m_out_dir, '/', _name, ".traj.p", mpi_rank() + 1, ".pdb");
#else
    m_traj = to_str(m_out_dir, '/', _name, ".traj.p1.pdb");
#endif
    log << "# Trajectory file: " << m_traj << std::endl;

}

void MCBase::init(const Par &par) {
    TSP::init(par);

    // grow
    m_grow_mode = par.has("grow");
    m_grow_steps = 20000;
    m_grow_length = 3;

    // save ss
    m_del_pk = par.has("del_pk");
    m_pk_ahead = par.has("pk_ahead");
    m_save_ss = par.has("save_ss");
    m_sample_tree = par.has("sample_tree");

    if (m_del_pk) {
        for (auto && c : _ss) {
            if (c != '(' && c != ')' && c != '&') {
                c = '.';
            }
        }
    }

    m_selected_mvel = NULL;
    m_sample_mode = SAMPLE_SSE;
    m_cal_en_constraints = true;
    m_max_angle = PI * 0.5;
    m_box = 2;
    m_box_size = 12;
    m_will_write_traj = !_name.empty();

    log << "# Extract residue conformations..." << std::endl;
    for_each_model(to_str(Env::lib(), "/RNA/pars/cg/CG2AA/templates.pdb"), [this](const Model &m, int n) {
            ResConf::extract(m_res_confs, m.residues());
            });
    for_each_model(to_str(Env::lib(), "/RNA/pars/cg/CG2AA/templates.pdb"), [this](const Model &m, int n) {
            FragConf<3>::extract(m_frag_confs, m.residues());
            });
    log << "confs: ";
    for (auto && p : m_frag_confs) log << p.first << '(' << size(p.second) << ") ";
    log << std::endl;
    //MolReader reader(to_str(Env::lib(), "/RNA/pars/cg/CG2AA/templates.pdb"));
    //for (auto it = reader.model_begin(); it != reader.model_end(); it++) {
    //	STD_ cerr << *it << STD_ endl;
    //	ResConf::extract(m_res_confs, it->residues());
    //}
    for (auto && pair : m_res_confs) {
        for (auto && conf : pair.second) {
            conf.res = m_cg->to_cg(conf.res);
        }
    }

    for (auto && pair : m_frag_confs) {
        for (auto && conf : pair.second) {
            for (auto && r : conf.frag) {
                r = m_cg->to_cg(r);
            }
        }
    }

    if (m_will_write_traj) {
        set_traj_name();
    }

    log << "# Read parameters..." << std::endl;
    m_par_file = m_cg_type;
    _par->set(m_par_file, "par_file");
    read_parameters();

    log << "# Set parameters..." << std::endl;
    set_parameters(*_par);

    log << "# Print parameters..." << std::endl;
    print_parameters();

    log << "# Set continuous points..." << std::endl;
    set_continuous_pts();

    log << "# Print continuous, angel, dihedral points..." << std::endl;
    for (auto && i : m_continuous_pts) log << i << ' '; log << std::endl;
    for (auto && i : m_ang_pts) log << i << ' '; log << std::endl;
    for (auto && i : m_dih_pts) log << i << ' '; log << std::endl;

    log << "# Read initial structure" << std::endl;
    par.set(m_init_sfile, "init");
    if (m_init_sfile.empty()) {
        if (par.has("init:chain")) {
            _pred_chain = BuildChain()(_seq.size()).m_chain;
        }
        else if (par.has("init:line")) {
            _pred_chain = BuildLine()(_seq.size()).m_chain;
        }
        else if (par.has("init:raw")) {
            nuc3d::Assemble assemble(Par(par)("seq", _seq)("ss", _ss)("loop_building", "partial_raw")("name", ""));
            assemble.predict_one();
            _pred_chain = assemble._pred_chain;
        }
    }
    else {
        chain_read_model(_pred_chain, m_init_sfile);
    }

    log << "# Read alignment file" << STD_ endl;
    par.set(m_alignfile, "align");
    if (!m_alignfile.empty()) {
        set_align();
        align_to_fixed_areas();
    }

    log << "# Check constraints" << std::endl;
    validate_constraints();

    //_mc_queue = "heat:30000:20+cool:1000000";
    par.set(_mc_queue, "queue");
}

void MCBase::set_align() {
    Str seq1, seq2;
    Int i, n1, n2, l;
    if (!m_alignfile.empty()) {
        for (auto &&it : FileLines(m_alignfile)) {
            if (it.n == 0) seq1 = it.line;
            else if (it.n == 1) seq2 = it.line;
        }
        if (size(seq1) != size(seq2)) die("Bad align file! The lengths of the two aligned sequences should be equal!");
        l = size(seq1);
        for (i = 0, n1 = 0, n2 = 0; i < l; i++) {
            if (seq1[i] != '-' && seq2[i] != '-') m_align.push_back({ n1, n2 });
            if (seq1[i] != '-') n1++;
            if (seq2[i] != '-') n2++;
        }
    }
#ifndef NDEBUG
    log << "# Print align" << STD_ endl;
    for (auto && pair : m_align) {
        log << pair[0] << ' ' << pair[1] << STD_ endl;
    }
#endif
}

void MCBase::align_to_fixed_areas() {
    Vector<Bool> v(_seq.size(), false);
    for (auto && pair : m_align) {
        v[pair[1]] = true;
    }
    m_fixed_areas.push_back(STD_ move(v));
}

void MCBase::mc_next_step() {
    _mc_step++;
    if (m_grow_length < _seq.size() && _mc_step % m_grow_steps == 0) {
        //for (int i = 0; i < _seq.size(); i++) std::cout << is_selected(i) << std::endl;
        m_grow_length += 3;
    }
    if (_mc_step >= 100000) {
        //m_sample_mode = SAMPLE_TREE;
        m_cal_en_constraints = true;
    }
}

void MCBase::validate_constraints() {
    int l = _seq.size();
    int i, j;
    for (auto && ct : _constraints.distances) {
        i = ct.key[0];
        j = ct.key[1];
        if (i < 0 || i >= l || j < 0 || j >= l) {
            std::cerr << "Illegal constraints pair (" << i << ", " << j << ") in your constraints file!" << std::endl;
            std::exit(1);
        }
    }
}

void MCBase::read_parameters() {
    Str filename = to_str(Env::lib(), "/RNA/pars/nuc3d/mc/", m_par_file, ".par");
    log << "Reading parameters file: " << filename << std::endl;
    Par temp_par(filename);
    JN_MAP(JN_MCXP_TEMPPAR_SET, JN_MCXP_PARS1, JN_MCXP_PARS2)
}

void MCBase::set_parameters(const Par &par) {
    JN_MAP(JN_MCXP_PAR_SET, JN_MCXP_PARS1, JN_MCXP_PARS2)
}

void MCBase::print_parameters() {
    JN_MAP(JN_MCXP_TEMP, JN_MCXP_PARS1, JN_MCXP_PARS2)
}

void MCBase::set_continuous_pts() {
    int i, n;

    const auto &bks = NASS::instance().break_keys;

    // set broken points
    // from chain lens
    for (n = 0, i = 0; i < m_chain_lens.size(); i++) {
        n += m_chain_lens[i];
        m_brk_pts.push_back(n - 1);
    }
    // from secondary structure
    for (n = 0, i = 0; i < _ss.size(); i++) {
        if (std::find(bks.begin(), bks.end(), _ss[i]) == bks.end()) n++;
        else if (std::find(m_brk_pts.begin(), m_brk_pts.end(), n - 1) == m_brk_pts.end()) m_brk_pts.push_back(n - 1);
    }

    // set continuous points
    for (i = 0; i < _seq.size() - 1; i++) {
        if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end()) {
            m_continuous_pts.push_back(i);
        }
    }
    for (i = 0; i < _seq.size() - 2; i++) {
        if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
                std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 1) == m_brk_pts.end()) {
            m_ang_pts.push_back(i);
        }
    }
    for (i = 0; i < _seq.size() - 3; i++) {
        if (std::find(m_brk_pts.begin(), m_brk_pts.end(), i) == m_brk_pts.end() &&
                std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 1) == m_brk_pts.end() &&
                std::find(m_brk_pts.begin(), m_brk_pts.end(), i + 2) == m_brk_pts.end()) {
            m_dih_pts.push_back(i);
        }
    }
}

void MCBase::mc_write() {
    if (m_will_write_traj) {
        write_traj();
    }
    write_en();
    mc_num_writing++;
}

void MCBase::write_traj() {
    if (mc_num_writing == 1) {
        file::clean(m_traj);
    }
    std::ofstream output(m_traj.c_str(), std::ios::app);
    m_writer.bind_stream(output);
    m_writer.write_model([&]() {
            this->m_writer.write(this->_pred_chain);
            });
    output.close();
}

void MCBase::mc_back() {
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

void MCBase::backup() {
    _moved_atoms.clear();

    for (int i = 0; i < _seq.size(); i++) {
        if (is_selected(i)) {
            for (auto && atom : _pred_chain[i]) {
                _moved_atoms.push_back(atom);
            }
        }
    }
}

void MCBase::init_space() {
    m_item_space.resize(_seq.size());
    for (int i = 0; i < _seq.size(); i++) {
        SpaceVal &s = space_val(i);
        s.push_back(i);
        m_item_space[i] = &s;
    }
}

int MCBase::space_index(double n) const {
    return int((n + 1000) / m_box_size);
}

MCBase::SpaceItem &MCBase::item(int i) {
    return _pred_chain[i][2];
}

MCBase::SpaceVal &MCBase::space_val(int i) {
    SpaceItem &a = item(i);
    return m_space[space_index(a[0])][space_index(a[1])][space_index(a[2])];
}

void MCBase::space_update_item(int i) {
    SpaceVal &n = space_val(i);
    SpaceVal &o = *(m_item_space[i]);
    if (&o != &n) {
        o.erase(std::find(o.begin(), o.end(), i));
        m_item_space[i] = &n;
        n.push_back(i);
    }
}

void MCBase::restore_raw() {
    Mat dist;
    Int l = size(_seq);
    dg_dist_init(dist, l);
    //dg_dist_read_ss(dist, _seq, _ss);
    Vector<Int> v(l);
    for (Int i = 0; i < l; i++) {
        auto it = STD_ find_if(m_align.begin(), m_align.end(), [i](auto &&p) {return p[0] == i; });
        if (it == m_align.end()) v[i] = -1;
        else v[i] = it->at(1);
    }
    dg_dist_read_chain(dist, _pred_chain, v);
    Mat && c = DG(dist)();
    auto cg = CG::fac_t::make("1p");
    _pred_chain = cg->to_aa(c, 0, c.rows() - 1);
}

void MCBase::run() {
    log << "# Check initial structure..." << std::endl;
    if (_pred_chain.empty()) die("Please give me an initial structure for optimizing!");
    if (!m_alignfile.empty()) restore_raw();

    log << "# Initializing running..." << std::endl;
    before_run();

    if (size(m_mvels) == 0) {
        log << "# Nothing to move, exit..." << std::endl;
        return;
    }

    log << "# Save fixed ranges..." << std::endl;
    save_fixed_ranges();

    log << "# Carrying out CG processing with the Chain..." << std::endl;
    _pred_chain = m_cg->to_cg(_pred_chain);

    log << "# Init space..." << std::endl;
    init_space();

    log << "# MC..." << std::endl;
    mc_run();

    log << "# Print Constraints and Distances..." << std::endl;
    print_final_constraints();

    log << "# Finishing running..." << std::endl;
    finish_run();

    log << "# Coarsed Grained To All Atom..." << std::endl;
    _pred_chain = m_cg->to_aa(_pred_chain);

    log << "# Restore fixed ranges..." << std::endl;
    restore_fixed_ranges();

    log << "# Transform sequences..." << std::endl;
    transform();

}

void MCBase::transform() {
    int i = 0;
    for (Residue &r : _pred_chain) {
        r = mutate(r, to_str(_seq[i]));
        i++;
    }
}

void MCBase::print_final_constraints() {
    double d;
    int i, j;

    log << "Print residue contacts:" << std::endl;
    for (auto && c : _constraints.contacts) {
        i = c.key[0];
        j = c.key[1];
        d = geom::distance(item(i), item(j));
        //d = dist_two_res(_pred_chain[i], _pred_chain[j]);
        log << i << ' ' << j << " weight:" << c.weight << " dist:" << d << std::endl;
    }

    log << "Print residue distances:" << std::endl;
    for (auto && c : _constraints.distances) {
        i = c.key[0];
        j = c.key[1];
        d = geom::distance(item(i), item(j));
        //d = dist_two_res(_pred_chain[i], _pred_chain[j]);
        log << i << ' ' << j << " value:" << c.value << " weight:" << c.weight << " dist:" << d << std::endl;
    }

    // base pairing
    log << "Print atom distances:" << std::endl;
    for (auto && c : m_distance_constraints) {
        d = geom::distance(_pred_chain[c.atom1.n_res][c.atom1.n_atom], _pred_chain[c.atom2.n_res][c.atom2.n_atom]);
        log << "Residue " << c.atom1.n_res << " Atom " << c.atom1.n_atom
            << ", Residue " << c.atom2.n_res << " Atom " << c.atom2.n_atom
            << ", min: " << c.min << ", max: " << c.max
            << ", value: " << d << std::endl;
    }
}

void MCBase::before_run() {}
void MCBase::finish_run() {}

//Str MCBase::file_parameters() const {
//    return "3drna";
//}
//
void MCBase::save_fixed_ranges() {}
void MCBase::restore_fixed_ranges() {}

END_JN

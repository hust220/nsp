#include <numeric>
#include "../nuc2d/NASS.hpp"
#include "DHMC.hpp"

BEGIN_JN

void DHMC::init(const Par &par) {
    MCSM::init(par);

    m_set_mvel_pk = false;
    m_pk_ahead = par.has("pk_ahead");
    m_sample_frag = par.has("frag");
    m_sample_all_res = par.has("sample_all_res");
    m_not_sample_hp = par.has("not_sample_hp");
    m_not_sample_il = par.has("not_sample_il");
    m_all_free = par.has("all_free");

    log << "# Load bp distances..." << std::endl;
    for (auto &&it : FileLines(to_str(Env::lib(), "/RNA/pars/nuc3d/bp_distances.txt"))) {
        if (size(it.arr) == 7) {
            for (int i = 0; i < 6; i++) {
                m_bp_distances[it.arr[0]][i] = lexical_cast<Num>(it.arr[i+1]);
            }
        }
    }

    log << "# Set bps" << std::endl;
    set_bps();

    log << "# Print bps" << std::endl;
    for (int i = 0; i < size(_seq); i++) if (m_bps[i] != -1) log << i << '-' << m_bps[i] << ' ';
    log << std::endl;

    log << "# Transform bps to constraints..." << std::endl;
    bps_to_constraints();

    log << "# Set 2D trees" << std::endl;
    set_trees();

    if (m_sample_frag) {
        log << "# Set fragments" << std::endl;
        m_frag_size = 3;
        par.set(m_frag_size, "frag_size");
        m_frags = &(ResFrags::instance(m_cg->m_cg, m_frag_size));
    }
    else {
        m_frag_size = 1;
    }

    log << "# Set moving elements" << std::endl;
    set_mvels();

    //log << "# Set moving elements of each base" << std::endl;
    //set_base_mvels();

    log << "# Print moving elements" << std::endl;
    print_mvels();

    //log << "# Remove useless constraints" << std::endl;
    //remove_useless_constraints();

    log << "# Print constraints" << std::endl;
    print_constraints();

}

DHMC::~DHMC() {
    for (auto && el : m_mvels) {
        delete el;
    }
    for (auto && h : m_saved_helices) {
        delete h.second;
    }
}

void DHMC::set_base_mvels() {
    m_base_mvels.resize(_seq.size());
    for (auto && el : m_mvels) {
        for (auto && frag : el->range) {
            for (int i = frag[0]; i <= frag[1]; i++) {
                m_base_mvels[i] = el;
            }
        }
    }
}

void DHMC::remove_useless_constraints() {
    Constraints cs;
    for (auto && c : _constraints.contacts) {
        if (m_base_mvels[c.key[0]] != m_base_mvels[c.key[1]]) {
            cs.contacts.push_back(c);
        }
    }
    for (auto && c : _constraints.distances) {
        if (m_base_mvels[c.key[0]] != m_base_mvels[c.key[1]]) {
            cs.distances.push_back(c);
        }
    }
    _constraints = cs;
}

void DHMC::print_constraints() {
    List<Array<Int, 2>> ls;
    for (auto && c : _constraints.distances) {
        ls.push_back({c.key[0], c.key[1]});
    }
    for (auto && c : _constraints.contacts) {
        ls.push_back({c.key[0], c.key[1]});
    }
    ls.sort([](auto &&a, auto &&b){return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);});
    for (auto && arr : ls) log << arr[0] << '-' << arr[1] << ' ';
    log << std::endl;
}

void DHMC::print_mvels() {
    for (auto && el : m_mvels) {
        log << *el << ' ';
    }
    log << std::endl;
}

void DHMC::set_trees() {
    const auto &pks = NASS::instance().paired_keys;
    const auto &bks = NASS::instance().break_keys;

    auto is_pks = [&pks](Char c) {
        return std::find_if(pks.begin(), pks.end(), [&c](auto &&p) {
                return p.first == c || p.second == c;
                }) != pks.end();
    };

    auto partial_ss = [&is_pks](S ss, const Pair<Char, Char> &pair) {
        for (auto && c : ss) {
            if (is_pks(c)) {
                c = (c == pair.first ? '(' : (c == pair.second ? ')' : '.'));
            }
        }
        return ss;
    };

    m_trees.push_back(std::make_shared<SSTree>());
    m_trees.back()->make_b(_seq, _ss, 2);
    for (auto it = pks.begin() + 1; it != pks.end(); it++) {
        auto ss = partial_ss(_ss, *it);
        if (std::any_of(ss.begin(), ss.end(), [](auto &&c) {return c != '.' && c != '&'; })) {
            m_trees.push_back(std::make_shared<SSTree>());
            m_trees.back()->make_b(_seq, ss, 1);
        }
        else {
            break;
        }
    }
}

void DHMC::set_bps() {
    m_bps = NASS::get_bps(_ss);
}

void DHMC::bps_to_constraints() {
    int i, j;

    auto foo = [this](int n1, int n2) {
        auto &dists = m_bp_distances[to_str(_seq[n1], _seq[n2])];
        //for (int i = 0; i < 1; i++) {
        for (int i = 0; i < 6; i++) {
            m_distance_constraints.push_back({ { n1, i },{ n2, i }, dists[i], dists[i] });
        }
    };

    for (i = 0; i < _seq.size(); i++) {
        j = m_bps[i];
        if (i < j) {
            foo(i, j);
            //_constraints.add_contact(i, j);
        }
    }
}

void DHMC::set_pseudo_knots() {
    auto translate_pseudo_knots_helix = [this](Model &m, const List<Int> &nums) {
        int n = m_cg->res_size() * nums.size() / 2;
        Mat x(n, 3), y(n, 3);
        int i = 0;
        int l = 0;
        for (auto && j : nums) {
            if (i < nums.size() / 2) {
                Residue && res1 = m_cg->to_cg(m[0][i]);
                Residue && res2 = m_cg->to_cg(_pred_chain[j]);
                for (int k = 0; k < m_cg->res_size(); k++) {
                    for (int t = 0; t < 3; t++) {
                        x(l, t) = res1[k][t];
                        y(l, t) = res2[k][t];
                    }
                    l++;
                }
            }
            i++;
        }
        geom::Superposition<Num> sp(x, y);
        for (auto && res : m[0]) {
            for (auto && atom : res) {
                sp.apply(atom);
            }
        }
    };

    auto foo = [this, &translate_pseudo_knots_helix](const Helix &h) {
        auto && seq = h.seq();
        auto && m = build_helix(seq);
        auto && nums = h.nums();

        assert(nums.size() >= 2 && nums.size() % 2 == 0);
        translate_pseudo_knots_helix(m, nums);

        int i = 0;
        for (auto && n : nums) {
            _pred_chain[n] = std::move(m[0][i]);
            i++;
        }
    };

    if (m_pk_ahead) {
        log << "# Set pseudo-knots" << std::endl;

        auto it = m_trees.begin();

        // 设置第一个二级结构树中长度为1的helix
        for (auto &&sse : **it) {
            if (sse.has_helix() && size(sse.helix) == 1) {
                foo(sse.helix);
            }
        }

        // 设置除了第一个二级结构树以外的其它的树中的helix
        for (it = m_trees.begin() + 1; it != m_trees.end(); it++) {
            for (auto &&sse : **it) {
                if (sse.has_helix()) {
                    foo(sse.helix);
                }
            }
        }
    }
}

void DHMC::transform_saved_helix(Chain *h, const Li & nums) {
    int n = m_cg->res_size() * nums.size();
    Mat x(n, 3), y(n, 3);
    int i = 0;
    int l = 0;
    for (auto && j : nums) {
        Residue && res1 = m_cg->to_cg(h->at(i));
        Residue && res2 = m_cg->to_cg(_pred_chain[j]);
        for (int k = 0; k < m_cg->res_size(); k++) {
            mat_set_rows(x, l, res1[k]);
            mat_set_rows(y, l, res2[k]);
            l++;
        }
        i++;
    }
    geom::Superposition<Num> sp(x, y);
    for (auto && res : *h) {
        for (auto && atom : res) {
            sp.apply(atom);
        }
    }
}

void DHMC::restore_helix(SSE *sse) {
    auto && seq = sse->helix.seq();
    auto && nums = sse->helix.nums();
    Chain *saved = m_saved_helices[sse];

    transform_saved_helix(saved, nums);

    int i = 0;
    for (auto && n : nums) {
        _pred_chain[n] = std::move(saved->at(i));
        i++;
    }
}

void DHMC::restore_fixed_ranges() {
    if (!m_all_free) {
        log << "# " << __FUNCTION__ << std::endl;

        auto end = (m_pk_ahead ? m_trees.end() : std::next(m_trees.begin()));
        for (auto it = m_trees.begin(); it != end; it++) {
            for (auto &&sse : **it) {
                if (sse.has_helix()) {
                    restore_helix(&sse);
                }
            }
        }
    }
}

void DHMC::set_mvels() {
    const auto &keys = NASS::instance().paired_keys;
    const auto &bks = NASS::instance().break_keys;
    Deque<Char> ss;
    for (auto && c : _ss) if (std::find(bks.begin(), bks.end(), c) == bks.end()) ss.push_back(c);

    auto is_paired = [this, &ss](int i, int j) {
        return m_bps[i] == j && (ss[i] == '(' || m_pk_ahead);
    };

    auto not_in_helix_ranges = [this, &ss, &is_paired, &keys](int i, int j) {
        if (is_paired(i, j) && i > 0 && j < size(_seq) - 1 && is_paired(i-1, j+1)) return false;
        auto end = std::next(keys.begin(), 1);
        std::vector<int> v(std::distance(keys.begin(), end), 0);
        int m;
        for (int n = i; n <= j; n++) {
            char c = ss[n];
            auto it1 = std::find_if(keys.begin(), end, [&c](const NASS::PairedKey &pair) {
                    return pair.first == c;
                    });
            auto it2 = std::find_if(keys.begin(), end, [&c](const NASS::PairedKey &pair) {
                    return pair.second == c;
                    });
            if (it1 != end) {
                m = std::distance(keys.begin(), it1);
                v[m]++;
            }
            else if (it2 != end) {
                m = std::distance(keys.begin(), it2);
                v[m]--;
                if (v[m] < 0) return false;
            }
        }
        return std::all_of(v.begin(), v.end(), [](int s) {return s == 0; });
    };

    auto area_min = [](const auto &area) {
        Int i = 0;
        for (; i < size(area); i++) if (area[i]) break;
        return i;
    };

    auto area_max = [](const auto &area) {
        Int i = size(area) - 1;
        for (; i >= 0; i--) if (area[i]) break;
        return i;
    };

    auto not_in_fixed_areas = [this, &area_min, &area_max](Int i, Int j) {
        for (auto && area : m_fixed_areas) {
            Int min = area_min(area);
            Int max = area_max(area);
            if (i <= min && j >= max) return true;
            for (Int n = i; n <= j; n++) if (area[n]) return false;
        }
        return true;
    };

    Int i, j, l = size(_seq);
    for (i = 0; i + 2 < l; i++) {
        m_mvels.push_back(new MvEl(i, i + 2, MvEl::MVEL_FG));
    }
//    for (i = 0; i < l; i++) {
//        for (j = i; j < l; j++) {
//            if (not_in_fixed_areas(i, j) && (m_all_free || not_in_helix_ranges(i, j))) {
//                m_mvels.push_back(new MvEl(i, j, MvEl::MVEL_FG));
//            }
//        }
//    }
}

// MC related methods

void DHMC::mc_sample() {
    if (m_selected_mvel->type == MvEl::MVEL_FG && m_sample_frag) {
        mc_sample_frag();
    }
    else {
        mc_sample_res();
    }
}

void DHMC::mc_sample_frag() {
    int i, j, k, id;
    Mat m, m1, m2;
    Mat *mat;

    if (m_selected_mvel->type != MvEl::MVEL_FG)
        throw "jian::DHMC::mc_sample_frag error!";

    backup();

    std::ostringstream stream;
    Frag &f = m_selected_mvel->range[0];
    for (i = f[0]; i <= f[1]; i++) {
        stream << _seq[i];
    }

    auto &ids = m_frags->m_ids[stream.str()];
    id = ids[int(ids.size() * rand())];
    mat = m_frags->m_mats[id];

    m = *mat;
    m1.resize(2 * m_cg->res_size(), 3);
    m2.resize(2 * m_cg->res_size(), 3);
    for (j = 0; j < m_cg->res_size(); j++) {
        for (k = 0; k < 3; k++) {
            m1(j, k) = m(j, k);
            m1(m_cg->res_size() + j, k) = m((m_frag_size - 1) * m_cg->res_size() + j, k);
            m2(j, k) = _pred_chain[f[0]][j][k];
            m2(m_cg->res_size() + j, k) = _pred_chain[f[0] + (m_frag_size - 1)][j][k];
        }
    }
    geom::Superposition<Num> sp(m1, m2);
    sp.apply_m(m);
    for (i = 0; i < m_frag_size; i++) {
        if (m_is_free[f[0] + i]) {
            for (j = 0; j < m_cg->res_size(); j++) {
                for (k = 0; k < 3; k++) {
                    _pred_chain[f[0] + i][j][k] = m(i * m_cg->res_size() + j, k);
                }
            }
            space_update_item(f[0] + i);
        }
    }
}

void DHMC::save_helix(SSE *sse) {
    auto && nums = sse->helix.nums();

    Chain *c = new Chain;

    for (auto && n : nums) {
        c->push_back(_pred_chain[n]);
    }

    m_saved_helices[sse] = c;
}

void DHMC::save_fixed_ranges() {
    if (!m_all_free) {
        log << "# " << __FUNCTION__ << std::endl;

        auto end = (m_pk_ahead ? m_trees.end() : std::next(m_trees.begin()));
        for (auto it = m_trees.begin(); it != end; it++) {
            for (auto &&sse : **it) {
                if (sse.has_helix()) {
                    save_helix(&sse);
                }
            }
        }
    }
}

void DHMC::before_run() {
    set_pseudo_knots();
}

void DHMC::mc_select() {
    m_selected_mvel = m_mvels[int(rand() * m_mvels.size())];
}

bool DHMC::is_selected(const I &i) const {
    if (m_selected_mvel == NULL) {
        return false;
    }
    else if (m_sample_mode == SAMPLE_SSE) {
        return m_selected_mvel->has(i);
    }
    else if (m_sample_mode == SAMPLE_TREE) {
        return m_selected_mvel->minmax_has(i);
    }
    else {
        throw "jian::DHMC::is_selected error! Illegal sample mode!";
    }
}

Vec DHMC::rotating_center() const {
    Vec origin = Vec::Zero(3);

    if (m_selected_mvel->type == MvEl::MVEL_FG) {
        for (I i = 0; i < 3; i++) origin[i] = _pred_chain[m_selected_mvel->range[0][0]][0][i];
    }
    else {
        int beg = m_selected_mvel->min();
        int end = m_selected_mvel->max();
        auto &r1 = _pred_chain[beg];
        auto &r2 = _pred_chain[end];
        int n_atom = 0;
        for (auto && atom : r1) {
            for (I i = 0; i < 3; i++) origin[i] += atom[i];
            n_atom++;
        }
        for (auto && atom : r2) {
            for (I i = 0; i < 3; i++) origin[i] += atom[i];
            n_atom++;
        }
        for (I i = 0; i < 3; i++) origin[i] /= n_atom;
    }
    return origin;
}

END_JN

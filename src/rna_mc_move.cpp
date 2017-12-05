#include "rna_mc_move.hpp"

BEGIN_JN

static auto rand_rot_mat(Num max_angle) {
    int index = int(rand() * 3);
    double dih = (rand() - 0.5) * max_angle;
    return geom::rot_mat(index, dih);
}

static auto helix_center(const Chain &chain, dhmc_mvel_t el) {
    Vec v(3);
    Num n;
    const Array<Int, 2> arr = *(Array<Int, 2> *) el.p;

    for (Int i = 0; i < 3; i++) {
        v[i] = 0;
        v[i] += chain[arr[0]][0][i];
        v[i] += chain[arr[1]][0][i];
        v[i] /= 2.0;
    }
    return v;
}

static void update_fragment(DHMC &m) {
    Array<Int, 2> arr = *(Array<Int, 2> *) m.m_selected_mvel.p;
    int min = arr[0];
    int max = arr[1];
//    int min = m.m_selected_mvel->min();
//    int max = m.m_selected_mvel->max();
//    auto type = m.m_selected_mvel->type;
    Str &_seq = m._seq;
    Str &_ss = m._ss;
    Chain &_pred_chain = m._pred_chain;

    if (m.m_grow_mode && max >= m.m_grow_length) return;

    std::stringstream stream;
    for (Int i = min; i <= max; i++) stream << m._pred_chain[i].name;
    Str frag_name = stream.str();

    FragConf<3>::Confs &confs = m.m_frag_confs[frag_name];
    Int l = size(confs);
    Int n = Int(rand()*l);
    //std::cout << frag_name << " size: " << l << ' ' << n << std::endl;
    FragConf<3> &conf = confs[n];
    auto frag = conf.frag;
    geom::Superposition<Num> sp;
    if (min == 0) {
        const Atom &p = m._pred_chain[max+1]["P"];
        for (int i = min; i <= max; i++) {
            if (m.is_dependent(i)) continue;
            int j = 0;
            for (auto && atom : frag[i-min]) {
                for (int k = 0; k < 3; k++) m._pred_chain[i][j][k] = atom[k] - conf.p2[k] + p[k];
                j++;
            }
        }
    }
    else if (max == _seq.size() - 1) {
        const Atom &p = m._pred_chain[min]["P"];
        for (int i = min; i <= max; i++) {
            if (m.is_dependent(i)) continue;
            int j = 0;
            for (auto && atom : frag[i-min]) {
                for (int k = 0; k < 3; k++) m._pred_chain[i][j][k] = atom[k] - conf.p1[k] + p[k];
                j++;
            }
        }
    }
    else {
        Mat x(2, 3);
        mat_set_rows(x, 0, m._pred_chain[min]["P"], m._pred_chain[max + 1]["P"]);
        sp.init(conf.pp, x);
        for (Int i = min; i <= max; i++) {
            if (m.is_dependent(i)) continue;
            for (auto && atom : frag[i-min]) sp.apply(atom);
            set_atoms(m._pred_chain[i], frag[i-min]);
        }
    }

    if (min == 0 || max == size(_seq) - 1) {
        int t = (min == 0 ? max : min);
        int index = int(rand() * 3);
        double dih = (rand() - 0.5) * m.m_max_angle;
        auto &&rot = geom::rot_mat(index, dih);
        Vec origin(3);
        vec_set(origin, _pred_chain[t][0]);
        for (int i = min; i <= max; i++) {
            if (m.is_dependent(i)) continue;
            for (auto && atom : _pred_chain[i]) {
                geom::rotate(atom, origin, rot);
            }
//            m.space_update_item(i);
        }
    }
    else {
        geom::RotateAlong<double> rotate_along(_pred_chain[min][0], _pred_chain[max + 1][0], m.m_max_angle * (rand() - 0.5));
        for (int i = min; i <= max; i++) {
            if (m.is_dependent(i)) continue;
            for (auto && atom : _pred_chain[i]) {
                rotate_along(atom);
            }
//            m.space_update_item(i);
        }
    }

}

// translate
static void translate_mvel(DHMC &m) {
    Chain &chain = m._pred_chain;
    Str &seq = m._seq;

    int index = int(rand() * 3);
    double dist = (rand() - 0.5) * 1 * m._mc_max_shift;
    for (int i = 0; i < size(seq); i++) {
        if (m.is_selected(i)) {
            for (auto && atom : chain[i]) {
                atom[index] += dist;
            }
//            m.space_update_item(i);
        }
    }
}

// rotate about center
static void rotate_about_center(DHMC &m) {
    Chain &chain = m._pred_chain;
    Str &seq = m._seq;

    auto &&rot = rand_rot_mat(m.m_max_angle);
    auto &&origin = helix_center(chain, m.m_selected_mvel);
    for (int i = 0; i < seq.size(); i++) {
        if (m.is_selected(i)) {
            for (auto && atom : chain[i]) {
                geom::rotate(atom, origin, rot);
            }
//            m.space_update_item(i);
        }
    }
}

static void load_template(DHMC &m, Str name) {
    auto cg = CG::fac_t::make("6p");
    auto & map = m.m_templates;
    Chain chain;
    map[name] = cg->to_cg(read_model_to_chain(to_str(Env::lib(), "/RNA/templates/", name, ".pdb")));
}

static Chain &get_template(DHMC &m, Str name) {
    auto & map = m.m_templates;
    if (map.find(name) == map.end()) {
        load_template(m, name);
    }
    return map[name];
}

static Chain &rand_helix_template(DHMC &m, SSE *sse) {
    auto & ls = m.m_helix_records[sse];
    Int n = size(ls);
    Int i = Int(JN_ rand() * n);
    Str name = ls[i].name;
    auto && helix = get_template(m, name);
    if (size(helix)/2 != sse->helix.size()) {
        std::cout << name << std::endl;
        throw "error!";
    }
    return helix;
//    return get_template(m, name);
}

static Chain &rand_loop_template(DHMC &m, SSE *sse) {
    auto & ls = m.m_loop_records[sse];
    Int n = size(ls);
    Int i = Int(JN_ rand() * n);
    Str name = ls[i].name;
    return get_template(m, name);
}

template<typename _Nums>
static void chain_set(Chain &chain, const Chain &c, _Nums && nums) {
    Int i = 0;
    for (auto && n : nums) {
        chain[n] = c[i];
        i++;
    }
}

template<typename _Nums>
static Mat chain_mat(const Chain &chain, _Nums && nums) {
    Mat m(4*6, 3);
    Int i = 0;
    for (auto && n : nums) {
        for (Int j = 0; j < 6; j++) {
            for (Int k = 0; k < 3; k++) {
                m(i*6+j, k) = chain[n][j][k];
            }
        }
        i++;
    }
    return m;
}

static Vector<Int> helix_head_inds(SSE *sse) {
    auto & h = sse->helix;
    Int l = size(h);
    return {0, 1, 2*l-2, 2*l-1};
}

static Vector<Int> helix_head_nums(SSE *sse) {
    auto & h = sse->helix;
    return {h[0].res1.num-1, h[1].res1.num-1, h[1].res2.num-1, h[0].res2.num-1};
}

static Vector<Int> helix_tail_inds(SSE *sse) {
    auto & h = sse->helix;
    Int l = size(h);
    return {l-2, l-1, l, l+1};
}

static Vector<Int> helix_tail_nums(SSE *sse) {
    auto & h = sse->helix;
    Int l = size(h);
    return {h[l-2].res1.num-1, h[l-1].res1.num-1, h[l-1].res2.num-1, h[l-2].res2.num-1};
}

static void dhmc_move_helix(DHMC &m) {
    auto el = m.m_selected_mvel;
    SSE *sse = (SSE *) el.p;
    auto & h = sse->helix;
    Int l = size(h);

    Chain & helix = rand_helix_template(m, sse);
    Chain & chain = m._pred_chain;

    auto sp = geom::suppos(chain_mat(helix, helix_head_inds(sse)), chain_mat(chain, helix_head_nums(sse)));
    for (auto && r : helix) for (auto && a : r) sp.apply(a);

    sp = geom::suppos(chain_mat(chain, helix_tail_nums(sse)), chain_mat(helix, helix_tail_inds(sse)));
    for (Int i = h[l-2].res1.num-1; i <= h[l-2].res2.num-1; i++) for (auto && a : chain[i]) sp.apply(a);

    chain_set(chain, helix, h.nums());
}

static Vector<Int> loop_head_inds(SSE *sse) {
    auto & l = sse->loop;
    Int sz = size(l);
    return {0, 1, sz-2, sz-1};
}

static Vector<Int> loop_head_nums(SSE *sse) {
    auto & l = sse->loop;
    Int sz = size(l);
    return {l[0].num-1, l[1].num-1, l[sz-2].num-1, l[sz-1].num-1};
}

static Vector<Int> loop_son_inds(SSE *sse, SSE *son) {
    Int n = son->helix[0].res1.num;
    auto & l = sse->loop;
    auto it = std::find_if(l.begin(), l.end(), [&n](auto && r){return r.num == n;});
    Int i = std::distance(l.begin(), it);
    return {i, i+1, i+2, i+3};
}

static void dhmc_move_loop(DHMC &m) {
    auto el = m.m_selected_mvel;
    SSE * sse = (SSE *) el.p;
    auto & l = sse->loop;

    Chain & loop = rand_loop_template(m, sse);
    Chain & chain = m._pred_chain;

    if (sse->has_helix()) {
        auto sp = geom::suppos(chain_mat(loop, loop_head_inds(sse)), chain_mat(chain, loop_head_nums(sse)));
        for (auto && r : loop) for (auto && a : r) sp.apply(a);
    }

    for_ (son, tree_sons(sse)) {
        auto nums = helix_head_nums(son);
        auto sp = geom::suppos(chain_mat(chain, nums), chain_mat(loop, loop_son_inds(sse, son)));
        for (Int i = nums[0]; i <= nums[3]; i++) for (auto && a : chain[i]) sp.apply(a);
    } end_;

    chain_set(chain, loop, l.nums());
}

static void dhmc_move_frag(DHMC &m) {
    //auto el = m.m_selected_mvel;
    if (JN_ rand() < 0.5) translate_mvel(m);
    else rotate_about_center(m);
}

static void dhmc_move_frag3(DHMC &m) {
    update_fragment(m);
}

static void dhmc_space_update(DHMC &m) {
    Int l = size(m._seq);

    for (int i = 0; i < l; i++) {
        if (m.is_selected(i)) {
            m.space_update_item(i);
        }
    }
}

void dhmc_move(DHMC &m) {
    auto el = m.m_selected_mvel;
    Int type = el.t;
    if (type == DHMC_MVEL_HELIX) dhmc_move_helix(m);
    else if (type == DHMC_MVEL_LOOP) dhmc_move_loop(m);
    else if (type == DHMC_MVEL_FRAG) dhmc_move_frag(m);
    else if (type == DHMC_MVEL_FRAG3) dhmc_move_frag3(m);
    else throw "error!";
    dhmc_space_update(m);
}

//void dhmc_sample_res(DHMC &m) {
//    m.backup();
//
//    MvEl *mvel = m.m_selected_mvel;
//    auto type = mvel->type;
//
//    if (type == MVEL_LOOP) mvel_loop_move();
//    else if (type == MVEL_FRAG) mvel_frag_move();
//    else if (type == MVEL_FRAG3) mvel_frag3_move();
//    else throw "Unknown mvel type!";
//
//}

//Bool is_dependent(DHMC &m, Int i) {
//    if (m.m_save_ss) {
//        if (m._ss[i] == '(' || m._ss[i] == ')') {
//            if ((i+1 < size(m._seq) && m.m_bps[i] - 1 == m.m_bps[i+1]) || (i > 0 && m.m_bps[i-1] - 1 == m.m_bps[i])) {
//                return true;
//            }
//            else {
//                return false;
//            }
//        }
//        else {
//            if (m.m_pk_ahead) {
//                return m._ss[i] != '.';
//            }
//            else {
//                return false;
//                //return m._ss[i] == '(' || m._ss[i] == ')';
//            }
//        }
//    }
//    else {
//        return false;
//    }
//}

//static auto get_base_axis(const Residue &residue) {
//    int t = pdb::res_type(residue.name);
//    std::array<Vec, 2> axis;
//    axis[0].resize(3);
//    axis[1].resize(3);
//    vec_set(axis[0], residue[1]);
//    if (t == 1 || t == 3)  vec_set(axis[1], std::plus<>{}, residue[3], residue[5]);
//    else if (t == 0 || t == 2) vec_set(axis[1], residue[5]);
//    else throw "jian::DHMC::mc_sample_res::get_base_axis error!";
//    return axis;
//}

END_JN


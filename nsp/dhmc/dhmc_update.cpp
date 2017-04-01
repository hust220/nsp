#include "dhmc_update.hpp"

BEGIN_JN

void dhmc_sample_res(DHMC &m) {
    m.backup();

    MvEl *mvel = m.m_selected_mvel;
    auto type = mvel->type;
    if (type == MvEl::MVEL_FG && mvel->range[0][1] - mvel->range[0][0] == 2) {
        update_fragment(m);
    } else {
        if (JN_ rand() < 0.5) {
            translate_mvel(m);
        }
        else {
            rotate_about_center(m);
        }
    }
}

auto rand_rot_mat(Num max_angle) {
    int index = int(rand() * 3);
    double dih = (rand() - 0.5) * max_angle;
    return geom::rot_mat(index, dih);
}

auto helix_center(const Chain &chain, MvEl *mvel) {
    Vec v(3);
    Num n;

    for (Int i = 0; i < 3; i++) {
        v[i] = 0;
        n = 0;
        for (auto && frag : mvel->range) {
            v[i] += chain[frag[0]][0][i];
            v[i] += chain[frag[1]][0][i];
            n += 2;
        }
        v[i] /= n;
    }
    return v;
}

Bool is_fixed(DHMC &m, Int i) {
    if (m.m_save_ss) {
        if (m._ss[i] == '(' || m._ss[i] == ')') {
            if ((i+1 < size(m._seq) && m.m_bps[i] - 1 == m.m_bps[i+1]) || (i > 0 && m.m_bps[i-1] - 1 == m.m_bps[i])) {
                return true;
            }
            else {
                return false;
            }
        }
        else {
            if (m.m_pk_ahead) {
                return m._ss[i] != '.';
            }
            else {
                return false;
                //return m._ss[i] == '(' || m._ss[i] == ')';
            }
        }
    }
    else {
        return false;
    }
}

auto get_base_axis(const Residue &residue) {
    int t = pdb::res_type(residue.name);
    std::array<Vec, 2> axis;
    axis[0].resize(3);
    axis[1].resize(3);
    vec_set(axis[0], residue[1]);
    if (t == 1 || t == 3)  vec_set(axis[1], std::plus<>{}, residue[3], residue[5]);
    else if (t == 0 || t == 2) vec_set(axis[1], residue[5]);
    else throw "jian::DHMC::mc_sample_res::get_base_axis error!";
    return axis;
}

void update_fragment(DHMC &m) {
    int min = m.m_selected_mvel->min();
    int max = m.m_selected_mvel->max();
    auto type = m.m_selected_mvel->type;
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
    FragConf<3> &conf = confs[n];
    auto frag = conf.frag;
    geom::Superposition<Num> sp;
    if (min == 0) {
        const Atom &p = m._pred_chain[max+1]["P"];
        for (int i = min; i <= max; i++) {
            if (is_fixed(m, i)) continue;
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
            if (is_fixed(m, i)) continue;
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
            if (is_fixed(m, i)) continue;
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
            if (is_fixed(m, i)) continue;
            for (auto && atom : _pred_chain[i]) {
                geom::rotate(atom, origin, rot);
            }
            m.space_update_item(i);
        }
    }
    else {
        geom::RotateAlong<double> rotate_along(_pred_chain[min][0], _pred_chain[max + 1][0], m.m_max_angle * (rand() - 0.5));
        for (int i = min; i <= max; i++) {
            if (is_fixed(m, i)) continue;
            for (auto && atom : _pred_chain[i]) {
                rotate_along(atom);
            }
            m.space_update_item(i);
        }
    }

}

// translate
void translate_mvel(DHMC &m) {
    Chain &chain = m._pred_chain;
    MvEl *mvel = m.m_selected_mvel;
    Str &seq = m._seq;

    int index = int(rand() * 3);
    double dist = (rand() - 0.5) * 1 * m._mc_max_shift;
    for (int i = 0; i < size(seq); i++) {
        if (m.is_selected(i)) {
            for (auto && atom : chain[i]) {
                atom[index] += dist;
            }
            m.space_update_item(i);
        }
    }
}

// rotate about center
void rotate_about_center(DHMC &m) {
    Chain &chain = m._pred_chain;
    MvEl *mvel = m.m_selected_mvel;
    Str &seq = m._seq;

    auto &&rot = rand_rot_mat(m.m_max_angle);
    auto &&origin = helix_center(chain, mvel);
    for (int i = 0; i < seq.size(); i++) {
        if (m.is_selected(i)) {
            for (auto && atom : chain[i]) {
                geom::rotate(atom, origin, rot);
            }
            m.space_update_item(i);
        }
    }
}

//void DHMC::mc_sample_frag() {
//    int i, j, k, id;
//    Mat m, m1, m2;
//    Mat *mat;
//
//    if (m_selected_mvel->type != MvEl::MVEL_FG)
//        throw "jian::DHMC::mc_sample_frag error!";
//
//    backup();
//
//    std::ostringstream stream;
//    Frag &f = m_selected_mvel->range[0];
//    for (i = f[0]; i <= f[1]; i++) {
//        stream << _seq[i];
//    }
//
//    auto &ids = m_frags->m_ids[stream.str()];
//    id = ids[int(ids.size() * rand())];
//    mat = m_frags->m_mats[id];
//
//    m = *mat;
//    m1.resize(2 * m_cg->res_size(), 3);
//    m2.resize(2 * m_cg->res_size(), 3);
//    for (j = 0; j < m_cg->res_size(); j++) {
//        for (k = 0; k < 3; k++) {
//            m1(j, k) = m(j, k);
//            m1(m_cg->res_size() + j, k) = m((m_frag_size - 1) * m_cg->res_size() + j, k);
//            m2(j, k) = _pred_chain[f[0]][j][k];
//            m2(m_cg->res_size() + j, k) = _pred_chain[f[0] + (m_frag_size - 1)][j][k];
//        }
//    }
//    geom::Superposition<Num> sp(m1, m2);
//    sp.apply_m(m);
//    for (i = 0; i < m_frag_size; i++) {
//        if (m_is_free[f[0] + i]) {
//            for (j = 0; j < m_cg->res_size(); j++) {
//                for (k = 0; k < 3; k++) {
//                    _pred_chain[f[0] + i][j][k] = m(i * m_cg->res_size() + j, k);
//                }
//            }
//            space_update_item(f[0] + i);
//        }
//    }
//}

END_JN

//int min = m_selected_mvel->min();
//int max = m_selected_mvel->max();
//if (m_sample_mode == SAMPLE_SSE) {
//	if (min == max) {
//		actions[std::vector<int>{0, 1, 2}[int(rand() * 3)]]();
//	}
//	else {
//		actions[std::vector<int>{0, 2}[int(rand()*2)]]();
//	}
//}
//else if (m_sample_mode == SAMPLE_TREE) {
//	if (min == max) {
//		actions[std::vector<int>{0, 1}[int(rand() * 2)]]();
//	}
//	else {
//		actions[0]();
//	}
//}
//else {
//	throw "jian::DHMC::sample_res error!";
//}

    /*
       if (min == 0 && max == _seq.size() - 1) {
       return;
       }
       else if (min == max) {
       ResConf::Confs &confs = m_res_confs[_pred_chain[min].name];
       int l = size(confs);
       int n = int(rand()*l);
       geom::Superposition<Num> sp;
       if (min == 0) {
       Mat x(1, 3);
       mat_set_rows(x, 0, _pred_chain[min + 1]["P"]);
       sp.init(confs[n].p2, x);
       }
       else if (max == _seq.size() - 1) {
       Mat x(1, 3);
       mat_set_rows(x, 0, _pred_chain[min]["P"]);
       sp.init(confs[n].p1, x);
       }
       else {
       Mat x(2, 3);
       mat_set_rows(x, 0, _pred_chain[min]["P"], _pred_chain[min + 1]["P"]);
       sp.init(confs[n].pp, x);
       }
       Residue r = confs[n].res;
       for (auto && atom : r) sp.apply(atom);
       set_atoms(_pred_chain[min], r);
       space_update_item(min);
       }
       else {
       Num d1 = -1;
       Num d2 = -1;
       if (min > 0) d1 = geom::distance(_pred_chain[min - 1][1], _pred_chain[min][0]);
       if (max < size(_seq) - 1) d2 = geom::distance(_pred_chain[max][1], _pred_chain[max + 1][0]);
       if (min == 0 || max == size(_seq) - 1 || d1 > 4.0 || d2 > 4.0) {
       int t = ((min == 0 || d1 > 4.0) ? max : min);
       int index = int(rand() * 3);
       double dih = (rand() - 0.5) * m_max_angle;
       auto &&rot = geom::rot_mat(index, dih);
       Vec origin(3);
       vec_set(origin, _pred_chain[t][0]);
       for (int i = min; i <= max; i++) {
       for (auto && atom : _pred_chain[i]) {
       geom::rotate(atom, origin, rot);
       }
       space_update_item(i);
       }
       }
       else {
       geom::RotateAlong<double> rotate_along(_pred_chain[min][0], _pred_chain[max + 1][0], m_max_angle * (rand() - 0.5));
       for (int i = min; i <= max; i++) {
       for (auto && atom : _pred_chain[i]) {
       rotate_along(atom);
       }
       space_update_item(i);
       }
       }
       }
       */

// rotate base
//void rotate_base(DHMC &m) {
//    for (int i = 0; i < _seq.size(); i++) {
//        if (is_selected(i)) {
//            auto axis = get_base_axis(_pred_chain[i]);
//            geom::RotateAlong<double> rotate_along(axis[0], axis[1], m_max_angle * (rand() - 0.5));
//            for (int j = 3; j < 6; j++) {
//                rotate_along(_pred_chain[i][j]);
//            }
//            space_update_item(i);
//        }
//    }
//}



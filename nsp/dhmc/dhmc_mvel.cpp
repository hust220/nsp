#include "DHMC.hpp"
#include "dhmc_mvel.hpp"

BEGIN_JN

std::ostream &operator <<(std::ostream &stream, const dhmc_mvel_t &el) {
    if (el.t == DHMC_MVEL_LOOP) {
        SSE *sse = (SSE *) el.p;
        stream << "DHMC_MVEL_LOOP(" << join(", ", sse->loop.nums()) << ")";
    }
    else if (el.t == DHMC_MVEL_HELIX) {
        SSE *sse = (SSE *) el.p;
        stream << "DHMC_MVEL_HELIX(" << join(", ", sse->helix.nums()) << ")";
    }
    else if (el.t == DHMC_MVEL_FRAG) {
        Array<Int, 2> *arr = (Array<Int, 2> *) el.p;
        stream << "DHMC_MVEL_FRAG(" << arr->at(0) << '-' << arr->at(1) << ")";
    }
    else if (el.t == DHMC_MVEL_FRAG3) {
        Array<Int, 2> *arr = (Array<Int, 2> *) el.p;
        stream << "DHMC_MVEL_FRAG3(" << arr->at(0) << '-' << arr->at(1) << ")";
    }
    else {
        throw "error!";
    }
    return stream;
}

void dhmc_mvel_free(dhmc_mvel_t el) {
    Int t = el.t;
    if (t == DHMC_MVEL_HELIX || t == DHMC_MVEL_LOOP) {
        SSE *sse = (SSE *) el.p;
        delete sse;
    }
    else if (t = DHMC_MVEL_FRAG || t == DHMC_MVEL_FRAG3) {
        Array<Int, 2> *p = (Array<Int, 2> *) el.p;
        delete p;
    }
    else {
        throw "dhmc_mvel_free error!";
    }
}

static Deque<Int> loop_nums(SSE *sse) {
    Deque<Int> dq;
    Int a = 0, b = 0;
    Int i = 0;
    for (auto && res : sse->loop) {
        if (res.type == '(') {
            a++;
            if (a % 2 != 0) {
                dq.push_back(res.num - 1);
            }
        } else if (res.type == ')') {
            b++;
            if (b % 2 != 1) {
                dq.push_back(res.num - 1);
            }
        } else {
            dq.push_back(res.num - 1);
        }
        i++;
    }
    if (sse->has_helix()) {
        dq.pop_front();
        dq.pop_back();
    }
    return dq;
}

static dhmc_mvel_t dhmc_mvel_helix_new(SSE * sse) {
    dhmc_mvel_t el;
    el.p = (void *) sse;
    el.t = DHMC_MVEL_HELIX;
    return el;
}

static dhmc_mvel_t dhmc_mvel_loop_new(SSE * sse) {
    dhmc_mvel_t el;
    el.p = (void *) sse;
    el.t = DHMC_MVEL_LOOP;
    return el;
}

static dhmc_mvel_t dhmc_mvel_frag_new(Int i, Int j) {
    dhmc_mvel_t el;
    Array<Int, 2> *p = new Array<Int, 2>{i, j};
    el.p = (void *) p;
    el.t = DHMC_MVEL_FRAG;
    return el;
}

static dhmc_mvel_t dhmc_mvel_frag3_new(Int i, Int j) {
    dhmc_mvel_t el;
    Array<Int, 2> *p = new Array<Int, 2>{i, j};
    el.p = (void *) p;
    el.t = DHMC_MVEL_FRAG3;
    return el;
}

static void set_mvels_helix(DHMC &m, SSE * sse) {
    m.m_mvels.push_back(dhmc_mvel_helix_new(sse));
    for (auto && i : sse->helix.nums()) m.m_dependent[i] = true;
}

static void set_mvels_loop(DHMC &m, SSE * sse) {
    m.m_mvels.push_back(dhmc_mvel_loop_new(sse));
    for (auto && i : sse->loop.nums()) m.m_dependent[i] = true;
}

static void set_mvels_frag(DHMC &m, SSE *sse) {
    auto nums = loop_nums(sse);
    Int l = size(nums);
    for (Int i = 0; i < l; i++) {
        for (Int j = i + 1; j < l; j++) {
            if (nums[j] - nums[i] > 2 && m._ss[nums[i]] != ')' && m._ss[nums[j]] != '(') {
                m.m_mvels.push_back(dhmc_mvel_frag_new(nums[i], nums[j]));
            }
        }
    }
}

static void set_mvels_frag3(DHMC &m) {
    Int i, j, l = size(m._seq);
    auto & v = m.m_dependent;
    for (i = 0; i + 2 < l; i++) {
        if (!v[i] || !v[i+1] || !v[i+2]) {
            m.m_mvels.push_back(dhmc_mvel_frag3_new(i, i+2));
        }
    }
}

void dhmc_set_mvels(DHMC &m) {
//    Map<Str, Chain> templates_cache;
//    Map<SSE *, Deque<Str>> templates;
//    Deque<Str> excludes, sources;
//    sst_find_templates(sst, templates, m.m_excludes, m.m_sources);
    Int sz = size(m._seq);
    m.m_dependent.resize(sz);
    for (auto && i : m.m_dependent) i = false;
    for_ (sse, tree_nodes(m.m_sst.head)) {
//        std::cout << *sse << std::endl;
        auto & dq = m.m_loop_records[sse];
        Int n = size(dq);
//        std::cout << "size: " << n << std::endl;
//        std::cout << "set helix mvels" << std::endl;
        if (sse->has_helix() && m.m_sample_helix_templates) set_mvels_helix(m, sse);
//        std::cout << "set loop mvels" << std::endl;
        if (sse->has_loop() && m.m_sample_loop_templates && n > 0) set_mvels_loop(m, sse);
        else set_mvels_frag(m, sse);
    } end_;
    set_mvels_frag3(m);
}

//static Int frag_has_frag(const Frag &f1, const Frag &f2) {
//    if (f2[1] < f1[0] || f2[0] > f1[1]) return -1;
//    else if (f2[0] >= f1[0] && f2[1] <= f1[1]) return 1;
//    else return 0;
//}
//
//static Int mvel_has_frag(MvEl *el, const Frag &frag) {
//    for (auto && f : el->range) {
//        Int n = frag_has_frag(f, frag);
//        if (n == 0) return 0;
//        else if (n == 1) return 1;
//    }
//    return -1;
//}
//
//static Bool mvel_compatible_el(MvEl *el, const Frags &frags) {
//    Int n = 0;
//    for (auto && f : frags) {
//        n += mvel_has_frag(el, f);
//    }
//    return std::abs(n) == size(frags);
//}
//
//static Bool mvel_compatible_els(MvEl *el, const Deque<Frags> &els) {
//    for (auto && frags : els) {
//        if (!mvel_compatible_el(el, frags)) return false;
//    }
//    return true;
//}
//
//void mvels_set_fixed_els(Deque<MvEl *> &mvels, const Deque<Frags> &els) {
//    Deque<MvEl *> ls;
//    for (auto && el : mvels) {
//        if (mvel_compatible_els(el, els)) {
//            ls.push_back(el);
//        }
//    }
//    mvels = ls;
//}
//
//void set_mvels_helices(DHMC &m) {
//    // set mvels of helices
//    auto end = (m.m_pk_ahead ? m.m_trees.end() : std::next(m.m_trees.begin()));
//    for (auto it = m.m_trees.begin(); it != end; it++) {
//        for (auto &&sse : **it) {
//            if (sse.has_helix()) {
//                auto &front = sse.helix.front();
//                auto &back = sse.helix.back();
//                m.m_mvels.push_back(new MvEl(front.res1.num-1, back.res1.num-1, back.res2.num-1, front.res2.num-1, MvEl::MVEL_HL));
//            }
//        }
//    }
//}
//
//void set_mvels_frag3(DHMC &m) {
//    Int i, j, l = size(m._seq);
//    for (i = 0; i + 2 < l; i++) {
//        if (!is_fixed(m, i) || !is_fixed(m, i+1) || !is_fixed(m, i+2)) {
//            m.m_mvels.push_back(new MvEl(i, i + 2, MvEl::MVEL_FG));
//        }
//    }
//}
//
//Bool is_small_bulge(const SSE &sse) {
//    if (sse.is_il()) {
//        Deque<Int> dq;
//        Int n = 0;
//        Char c = '(';
//        for (auto && res : sse.loop) {
//            if (res.type == '(' || res.type == ')') {
//                if (c == '(' || c == ')') {
//                    // pass
//                }
//                else {
//                    dq.push_back(n);
//                    n = 0;
//                }
//            }
//            else {
//                n++;
//            }
//            c = res.type;
//        }
//        if (c != '(' && c != ')') {
//            dq.push_back(n);
//            n = 0;
//        }
//        return size(dq) == 1 && dq[0] < 4;
//    }
//    else {
//        return false;
//    }
//}
//
//void set_mvels_tree(DHMC &m) {
//    for (auto && sse : *(m.m_trees.front())) {
//        if (sse.has_loop() && !(m.m_save_bg && is_small_bulge(sse))) {
//            auto nums = loop_nums(sse);
//            Int l = size(nums);
//            for (Int i = 0; i < l; i++) {
//                for (Int j = i + 1; j < l; j++) {
//                    if (nums[j] - nums[i] > 2 && m._ss[nums[i]] != ')' && m._ss[nums[j]] != '(') {
//                        m.m_mvels.push_back(new MvEl(nums[i], nums[j], MvEl::MVEL_FG));
//                    }
//                }
//            }
//        }
//    }
//}
//
//void set_mvels_all(DHMC &m) {
//    Int l = size(m._seq);
//    for (Int i = 0; i < l; i++) {
//        for (Int j = i; j < l; j++) {
//            //if (not_in_fixed_areas(i, j) && (m_all_free || not_in_helix_ranges(i, j))) {
//            m.m_mvels.push_back(new MvEl(i, j, MvEl::MVEL_FG));
//        }
//    }
//}
//
//static void loop_free_global(DHMC &m) {
//    for (auto && sse : *(m.m_trees.front())) {
//        if (sse.has_loop() && !(m.m_save_bg && is_small_bulge(sse))) {
//            auto nums = loop_nums(sse);
//            Int l = size(nums);
//            for (Int i = 0; i < l; i++) {
//                for (Int j = i + 1; j < l; j++) {
//                    if (nums[j] - nums[i] > 2 && m._ss[nums[i]] != ')' && m._ss[nums[j]] != '(') {
//                        m.m_mvels.push_back(new MvEl(nums[i], nums[j], MvEl::MVEL_FG));
//                    }
//                }
//            }
//        }
//    }
//}
//
//static void loop_free_local(DHMC &m) {
//    Int l = size(m._seq);
//    Vector<Bool> v(l);
//    for (auto && i : v) i = true;
//    auto end = (m.m_pk_ahead ? m.m_trees.end() : std::next(m.m_trees.begin()));
//    for (auto it = m.m_trees.begin(); it != end; it++) {
//        for (auto &&sse : **it) {
//            if (sse.has_helix()) {
//                auto &front = sse.helix.front();
//                auto &back = sse.helix.back();
//                Int a1 = front.res1.num-1, a2 = back.res1.num-1, b1 = front.res2.num-1, b2 = back.res2.num-1;
//                m.m_mvels.push_back(new MvEl(a1, a2, b2, b1, MvEl::MVEL_HL));
//                for (Int i = a1; i < a2; i++) v[i] = false;
//                for (Int i = b2; i < b1; i++) v[i] = false;
//            }
//        }
//    }
//    for (Int i = 0; i < l; i++) if (v[i]) m.m_mvels.push_back(new MvEl(i, i, MvEl::MVEL_FG));
//}
//
//static void all_free_global(DHMC &m) {
//    Int l = size(m._seq);
//    for (Int i = 0; i < l; i++) {
//        for (Int j = i + 3; j < l; j++) {
//            m.m_mvels.push_back(new MvEl(i, j, MvEl::MVEL_FG));
//        }
//    }
//}
//
//static void all_free_local(DHMC &m) {
//    Int l = size(m._seq);
//    for (Int i = 0; i < l; i++) m.m_mvels.push_back(new MvEl(i, i, MvEl::MVEL_FG));
//}

END_JN


#include "dca.hpp"
#include "rss_ss_tree.hpp"
#include "rss_sse.hpp"
#include "nsp.hpp"

namespace jian {

namespace dca {

TreePath<SSE> direct_path(const TreePath<SSE> &p1, const TreePath<SSE> &p2) {
    auto pos = std::adjacent_find(p1.begin(), p1.end(), [&p2](const SSE *a1, const SSE *a2) {
        return std::adjacent_find(p2.begin(), p2.end(), [&a1, &a2](const SSE *b1, const SSE *b2) {
            return a1 == b1 && a2 != b2;
        }) != p2.end() || size(p2) == 1;
    });
    if (pos == p1.end() && size(p1) == 1) pos = p1.begin();
    TreePath<SSE> path;
    for (auto it = std::next(p1.rbegin()); it != p1.rend(); it++) {
        if (*it == *pos) {
            path.push_back(*it);
            break;
        }
        else {
            path.push_back(*it);
        }
    }
    bool flag = false;
    for (auto it = p2.begin(); std::next(it) != p2.end(); it++) {
        if (*it == *pos) {
            flag = true;
            continue;
        }
        if (flag) path.push_back(*it);
    }
    return path;
}

void print_path(const TreePath<SSE> &path) {
    for (auto && l : path) {
        LOG << l << ' ';
    }
    LOG << std::endl;
}

bool connect_by_ml(int a, int b, const SSTree &sst) {
    TreePath<SSE> path_a, path_b;
    for (auto &&sse : sst) {
        if (sse.has(a + 1)) path_a.push_back(&sse);
        if (sse.has(b + 1)) path_b.push_back(&sse);
    }
    auto path = direct_path(path_a, path_b);
    LOG << a << ' ' << b << std::endl;
    print_path(path_a);
    print_path(path_b);
    print_path(path);
    return std::any_of(path.begin(), path.end(), [](const SSE *l) {return l->is_ml(); });
}

void trim_di(const Par &par) {
    Str ss = par.get("ss");
    int l = size(ss);
    Str seq;
    seq.resize(l, 'A');
    Str difile = par.get("di");
    tuples_t tuples = tuples_from_file(difile, ss.size());
    tuples_t ls;
    SSTree sst(seq, ss, 1);
    V<SSTree::El *> v;
    v.resize(l);
    for (Int i = 0; i < l; i++) {
        v[i] = std::find_if(sst.begin(), sst.end(), [&i](auto &&sse) {return sse.has(i + 1); }).el;
    }
    for (auto && l : v) LOG << l << ' '; LOG << std::endl;
    for (auto && tuple : tuples) {
        if (v[tuple.a] != v[tuple.b]) {
            ls.push_back(tuple);
        }
    }
    print_tuples(JN_OUT, ls);
}

void sort_di(const Par &par) {
    S seq = par.get("seq");
    Num k = 0.2;
    S di_file = par.get("di");

    par.set(k, "k");

    //LOG << "Sequence:" << std::endl;
    //LOG << seq << std::endl;

    int l = int(size(seq) * k);
    //LOG << "Read first " << l << " pairs:" << std::endl;
    dca::pairs_t &&pairs = dca::pairs_from_file(di_file, l);
    dca::print_pairs(JN_OUT, pairs);
}

void pred_di(const Par &par) {
    //        int n = 1;
    float step = 1;
    float pw = 0.5;
    S mol_type = "RNA";
    S method = "mf";
    S out_file = par.get("o", "out");
    S fa_file = par.get("i", "in");

    par.set(method, "m", "method");
    //        par.set(n, "n");
    par.set(step, "step");
    par.set(mol_type, "t", "type");
    par.set(pw, "w");

    //        dca_run(fa_file, out_file, mol_type, method);
    Dca *dca = FacDca::create(method, mol_type, pw);
    //        if (method == "mp") dca->set_step(step);
    //        dca->run(fa_file, out_file, n - 1);
    dca->run(fa_file, out_file);
    delete dca;
}

void rm_fp(const Par &par) {
    Str di_file;
    Str seq;
    Num cutoff = 1000;
    Num k = 0.5;
    Num frac_reserved;
    Int nr;
    int i, j, l;
    Vb v;
    tuples_t::iterator it, it1, it2;
    tuples_t pairs, new_pairs;

    auto neighbor = [](auto it1, auto it2) {
        return 
            (it1->a + 1 == it2->a && it1->b - 1 == it2->b) ||
            (it1->a - 1 == it2->a && it1->b + 1 == it2->b);
    };

    di_file = par.get("di");
    seq = par.get("seq");

    par.set(k, "k");
    l = int(size(seq) * k);

    frac_reserved = 0;
    par.set(frac_reserved, "fr", "frac_reserved");
    nr = Int(frac_reserved * size(seq));

    par.set(cutoff, "c", "cutoff");

    pairs = tuples_from_file(di_file, l);
    std::sort(pairs.begin(), pairs.end(), [](auto && t1, auto && t2){return t1.c > t2.c;});
    v.resize(l, false);
    for (it1 = pairs.begin(), i = 0; it1 != pairs.end(); it1++, i++) {
        for (it2 = std::next(it1), j = i + 1; it2 != pairs.end(); it2++, j++) {
            if (neighbor(it1, it2)) {
                v[i] = true;
                v[j] = true;
            }
        }
    }
    for (i = 0, it = pairs.begin(); it != pairs.end(); i++, it++) {
        if (v[i] || i <= nr || it->c > cutoff) new_pairs.push_back(*it);
    }
    dca::print_tuples(JN_OUT, new_pairs);
}

REGISTER_NSP_COMPONENT(dca) {
    Par::pars_t global = par.getv("global");
    if (size(global) == 2) {
        if (global[1] == "trim_di") {
            trim_di(par);
        }
        else if (global[1] == "sort_di") {
            sort_di(par);
        }
        else if (global[1] == "rm_fp") {
            rm_fp(par);
        }
    }
    else {
        pred_di(par);
    }

}

} // namespace dca

}


#include "build_strand.hpp"
#include <jian/utils/rand.hpp>
#include <jian/geom.hpp>
#include <jian/utils/Env.hpp>
#include "../scoring/FragConf.hpp"

BEGIN_JN

using Sp = geom::Superposition<Num>;

struct FragConf2 {
    FragConf<2>::MapConfs map_confs;

    FragConf2() {
        for_each_model(to_str(Env::lib(), "/RNA/pars/cg/CG2AA/templates.pdb"), [this](const Model &m, int n) {
                FragConf<2>::extract(map_confs, m.residues());
                });
    }

    void set_frag(Chain &frag, Str name) {
        auto & confs = map_confs[name];
        Int l = size(confs);
        Int n = Int(rand()*l);
        auto & conf = confs[n];
        frag.resize(2);
        frag[0] = conf.frag[0];
        frag[1] = conf.frag[1];
    }
};

static void set_frag(Chain &frag, Char seq1, Char seq2, Char ss1, Char ss2) {
    static FragConf2 conf;
    //static Int n = 0;
    //n++;
    conf.set_frag(frag, to_str(seq1, seq2));
    //mol_write(frag, to_str("aa-", n, ".pdb"));
}

static Sp res_sp(const Residue &r1, const Residue &r2) {
    Mat m1, m2;
    Int l1 = size(r1), l2 = size(r2);
    if (l1 != l2) throw to_str("build_strand.cpp line: ", __LINE__, " error!");
    m1.resize(l1, 3);
    m2.resize(l2, 3);
    for (Int i = 0; i < l1; i++) {
        for (Int j = 0; j < 3; j++) {
            m1(i, j) = r1[i][j];
            m2(i, j) = r2[i][j];
        }
    }
    return geom::suppos(m2, m1);
}

static void strand_append(Chain &strand, Chain &frag) {
    Sp sp = res_sp(strand.back(), frag.front());
    for (auto && res : frag) {
        for (auto && atom : res) {
            sp.apply(atom);
        }
    }
    for (auto it = std::next(frag.begin()); it != frag.end(); it++) {
        strand.push_back(std::move(*it));
    }
}

Chain build_strand(Str seq, Str ss = "..") {
    Int l = size(ss);
    if (l < 2) throw "Shortest strand: 2 nt!";

    Chain strand, frag;
    set_frag(strand, seq[0], seq[1], ss[0], ss[1]);
    for (Int i = 1; i < l-1; i++) {
        set_frag(frag, seq[i], seq[i+1], ss[i], ss[i+1]);
        strand_append(strand, frag);
    }
    return strand;
}

static Bool is_helix(Char a, Char b) {
    return a == b && a != '.';
}

static Int select_position(Str ss) {
    Deque<Int> ls;
    Int l = size(ss);
    for (Int i = 0; i < l-1; i++) {
        if (!is_helix(ss[i], ss[i+1])) ls.push_back(i);
    }

    Int n = size(ls);
    if (n == 0) return -1;
    return ls[Int(rand()*n)];
}

static Str chain_seq(const Chain &chain) {
    Str seq;
    for (auto && res : chain) seq += res.name.back();
    return seq;
}

void sample_strand(Chain &strand, Str ss) {
    Int n = select_position(ss);
    if (n == -1) return;

    Int l = size(ss);

    Chain f;
    Str seq = chain_seq(strand);
    set_frag(f, seq[n], seq[n+1], ss[n], ss[n+1]);

    Sp sp = res_sp(strand[n], f[0]);
    for (auto && res : f) for (auto && atom : res) sp.apply(atom);

    sp = res_sp(f[1], strand[n+1]);
    for (Int i = n+1; i < l; i++) for (auto && atom : strand[i]) sp.apply(atom);
}

END_JN


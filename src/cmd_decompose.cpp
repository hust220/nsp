#include "nsp.hpp"
#include "pdb.hpp"
#include "rss.hpp"
#include "rss_get_ss.hpp"

BEGIN_JN

static Vector<Int> get_loop_nums(SSE *l) {
    Vector<Int> nums;
    for (auto && res : l->loop) {
        nums.push_back(res.num - 1);
    }
    return nums;
}

static Vector<Int> get_helix_nums(SSE *l) {
    Int len = size(l->helix);
    Vector<Int> nums(len * 2);
    Int index = 0;
    for (auto && b : l->helix) {
        nums[index] = b.res1.num - 1;
        nums[2 * len - 1 - index] = b.res2.num - 1;
        index++;
    }
    return nums;
}

static void parse_helix(SSE *l, Int &n_helix, const Model &m, Str dirname, Str fam = "other") {
    if (l->has_helix()) {
        n_helix++;
        Str helix_name = to_str(m.name, '.', n_helix, ".helix");
        mol_write(sub(m, get_helix_nums(l)), to_str(dirname, '/', helix_name, ".pdb"));
        JN_OUT << helix_name << ' ' << size(l->helix) << ' ' << l->helix.seq() << ' ' << l->helix.ss() << ' ' << fam << STD_ endl;
    }
}

static void parse_loop(SSE *l, Int &n_loop, const Model &m, Str dirname, Str fam = "other") {
    if (l->has_loop()) {
        n_loop++;
        Str loop_name = to_str(m.name, '.', n_loop, ".loop");
        mol_write(sub(m, get_loop_nums(l)), to_str(dirname, '/', loop_name, ".pdb"));
        JN_OUT << loop_name << ' ' << l->num_sons()+1 << ' ' << l->loop.seq() << ' ' << l->loop.ss() << ' ' << fam << std::endl;
    }
}

static void decompose(Str filename, Str dirname) {
    Model m = mol_read_to<Model>(filename);

    Str seq = JN_ seq(m);
    Str ss = get_ss(m.residues());

    SSTree ss_tree;
    Int hinge = 1;
    ss_tree.make(seq, ss, hinge);

    Int n_helix = 0;
    Int n_loop = 0;
    for (auto && sse : ss_tree) {
        parse_helix(&sse, n_helix, m, dirname);
        parse_loop(&sse, n_loop, m, dirname);
    }
}

REGISTER_NSP_COMPONENT(decompose) {
    auto g = par.getv("global");

    Str filename = g[1];
    Str dir = g[2];

    decompose(filename, dir);

}

END_JN


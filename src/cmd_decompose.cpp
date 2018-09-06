#include "nsp.hpp"
#include "pdb.hpp"
#include "rss.hpp"
#include "rss_get_ss.hpp"
#include "rss_sst.hpp"

namespace jian {

struct Decompose {
    // pars
    Int n_bps;
    Str filename;
    Str dirname;
    Bool is_full;

    // temporary variables
    Str seq, ss;

    void operator ()(const Par &par) {
        auto g = par.getv("global");

        filename = g[1];
        dirname = g[2];

        // Read the model
        Model m = mol_read_to<Model>(filename);

//        seq = jian::seq(m);
//        ss = get_ss(m.residues());

        // n_bps is the number of base pairs that would be retained for the loop.
        n_bps = 1;
        par.set(n_bps, "n", "n_bps");

        // Full refers to retaining all the basepairs of the loop rather than only n base pairs.
        // Please don't use the full option when building the templates library.
        is_full = par.has("full");

        // Set sequence
        seq = JN_ seq(m);

        // Set secondary structure
        if (par.has("ss")) {
            par.set(ss, "ss");
        }
        else {
            ss = get_ss(m.residues());
        }

        if (is_full) {
            parse_full(m);
        }
        else {
            parse(m);
        }

    }

    void parse_full(const Model &m) {
        SST sst = sst_new(seq, ss, n_bps);
        Int n = 0;
        for_ (sse, tree_nodes(sst.head)) {
            parse_sse_full(sse, n, m);
        } end_;
        tree_free(sst.head);
    }

    template<typename T, typename U>
    Str str_sub(T && nums, U && str) {
        std::stringstream stream;
        for (auto && i : nums) {
            stream << str[i];
        }
        return stream.str();
    }

    template<typename T>
    Str nums_str(T && nums) {
        List<Str> ls;
        for (auto && i : nums) {
            ls.push_back(to_str(i));
            ls.push_back("-");
        }
        ls.pop_back();
        Str str;
        for (auto && c : ls) str += c;
        return str;
    }

    Int num_sons(SSE *sse) {
        return sse->num_sons();
    }

    Int num_branches(SSE *sse) {
        if (sse->has_helix()) return num_sons(sse)+1;
        else return num_sons(sse);
    }

    void parse_sse_full(SSE *sse, Int &n, const Model &m) {
        if (sse->has_loop()) {
            n++;
            auto nums = get_sse_full_nums(sse);
            Str sse_name = to_str(m.name, '.', n, ".loop.full");
            mol_write(sub(m, nums), to_str(dirname, '/', sse_name, ".pdb"));
            JN_OUT << sse_name << ' ' << num_sons(sse) << ' ' << num_branches(sse) << ' ' << str_sub(nums, seq) << ' ' << str_sub(nums, ss) << ' ' << nums_str(nums) << std::endl;
        }
    }

    Set<Int> get_sse_full_nums(SSE *sse) {
        Set<Int> set;
        for (auto && i : get_helix_nums(sse)) set.insert(i);
        for (auto && i : get_loop_nums(sse)) set.insert(i);
        for_ (son, tree_sons(sse)) {
            for (auto && i : get_helix_nums(son)) set.insert(i);
        } end_;
        return set;
    }

    void parse(const Model &m) {
        SSTree ss_tree;
        ss_tree.make(seq, ss, n_bps);

        Int n_helix = 0;
        Int n_loop = 0;
        for (auto && sse : ss_tree) {
            parse_helix(&sse, n_helix, m);
            parse_loop(&sse, n_loop, m);
        }
    }

    Vector<Int> get_loop_nums(SSE *l) {
        Vector<Int> nums;
        for (auto && res : l->loop) {
            nums.push_back(res.num - 1);
        }
        return nums;
    }

    Vector<Int> get_helix_nums(SSE *l) {
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

    void parse_helix(SSE *l, Int &n_helix, const Model &m) {
        if (l->has_helix()) {
            n_helix++;
            Str helix_name = to_str(m.name, '.', n_helix, ".helix");
            mol_write(sub(m, get_helix_nums(l)), to_str(dirname, '/', helix_name, ".pdb"));
            JN_OUT << helix_name << ' ' << size(l->helix) << ' ' << l->helix.seq() << ' ' << l->helix.ss() << ' ' << "other" << std::endl;
        }
    }

    void parse_loop(SSE *l, Int &n_loop, const Model &m) {
        if (l->has_loop()) {
            n_loop++;
            Str loop_name = to_str(m.name, '.', n_loop, ".loop");
            mol_write(sub(m, get_loop_nums(l)), to_str(dirname, '/', loop_name, ".pdb"));
            JN_OUT << loop_name << ' ' << l->num_sons()+1 << ' ' << l->loop.seq() << ' ' << l->loop.ss() << ' ' << "other" << std::endl;
        }
    }

};

REGISTER_NSP_COMPONENT(decompose) {
    Decompose decompose;
    decompose(par);
}

}


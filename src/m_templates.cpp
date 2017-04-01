#include "nsp.hpp"
#include <nsp/pdb.hpp>
#include <nsp/nuc2d.hpp>
#include <nsp/nuc3d/TemplRec.hpp>

BEGIN_JN

namespace {

    template<typename _List>
    static void print_ls(_List && ls) {
        for (auto && i : ls) JN_OUT << i << ' '; JN_OUT << std::endl;
    }

    void find_templates(Str seq, Str ss) {
        auto && sst = ss_tree(seq, ss, 2);
        JN_OUT << sst;
        for (auto && sse : sst) {
            JN_OUT << sse << std::endl;
            auto && ls1 = find_helix_templates(sse.helix.seq());
            auto && ls2 = find_loop_templates(sse.loop.seq(), sse.loop.ss());
            for (auto && i : ls1) JN_OUT << i.name << ' '; JN_OUT << std::endl;
            for (auto && i : ls2) JN_OUT << i.name << ' '; JN_OUT << std::endl;
        }
    }

    static void print_ab(Int a, Int b) {
        if (a == b) JN_OUT << a+1;
        else JN_OUT << a+1 << '-' << b+1;
    }

    template<typename _List>
    void write_els(_List && ls, Bool is_last) {
        if (size(ls) == 0) return;

        Int a, b;
        a = *ls.begin();
        b = a;
        for (auto it = std::next(ls.begin()); it != ls.end(); it++) {
            if (*it - b == 1) {
                b = *it;
            }
            else {
                print_ab(a, b);
                JN_OUT << '+';
                a = *it;
                b = *it;
            }
        }
        print_ab(a, b);
        if (!is_last) JN_OUT << ':';
    }

    void find_els(Str seq, Str ss, Int n) {
        auto && sst = ss_tree(seq, ss, 2);
        for (auto it = sst.begin(); it != sst.end(); it++) {
            auto & sse = *it;
            Set<Int> set;
            for (auto && i : sse.helix.nums()) set.insert(i);
            auto && ls2 = find_loop_templates(sse.loop.seq(), sse.loop.ss());
            if (size(ls2) > n) {
                for (auto && i : sse.loop.nums()) set.insert(i);
            }
            write_els(set, std::next(it) == sst.end());
        }
    }

    REGISTER_NSP_COMPONENT(templates) {
        auto global = par.getv("global");
        auto sub = global[1];
        Str seq, ss;
        if (sub == "find") {
            ss = par.get("ss");
            if (!par.has("seq")) seq = Str(size(ss), 'X');
            Int n = 0;
            par.set(n, "n");
            if (par.has("write_els")) {
                find_els(seq, ss, n);
            }
            else {
                find_templates(seq, ss);
            }
        }
    }

}

END_JN


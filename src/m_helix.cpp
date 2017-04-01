#include "nsp.hpp"
#include <nsp/rtsp/build_helix.hpp>

BEGIN_JN

namespace {

    using Bp = Array<Int, 2>;
    using Helix = Deque<Bp>;
    using Helices = Deque<Helix>;
    using Inds = Vector<Int>;

    void print_helices(const Helices &hs) {
        Int n = 0;
        JN_OUT << "=== Helices ===" << std::endl;
        for (auto && h : hs) {
            JN_OUT << "Helix " << n + 1 << ": ";
            JN_OUT << h.front()[0]+1 << '-' << h.back()[0]+1 << ' ';
            JN_OUT << h.back()[1]+1 << '-' << h.front()[1]+1 << ' ';
            JN_OUT << '(' << size(h) << ')' << ' ';
            JN_OUT << '(' << h.back()[1] - h.back()[0] - 1 << ')' << std::endl;
            /*
            for (auto && bp : h) {
                JN_OUT << bp[0]+1 << '-' << bp[1]+1 << ' ';
            }
            JN_OUT << std::endl;
            */
            n++;
        }
    }

    Int code(Char c) {
        static Map<Char, Int> m {{'A', 1}, {'U', 2}, {'T', 2}, {'G', 4}, {'C', 8}};
        return m[c];
    }

    Bool is_bp(Int a, Int b) {
        Int n = a + b;
        return n == 3 || n == 12 || n == 6;
    }

    Inds seq_to_inds(Str seq) {
        Int l = size(seq);
        Inds inds(l);
        for (Int i = 0; i < l; i++) inds[i] = code(seq[i]);
        return inds;
    }

    Bool is_adjacent(Bp bp1, Bp bp2) {
        return bp1[0] + 1 == bp2[0] && bp1[1] - 1 == bp2[1];
    }

    Bool is_same(const Helix &h1, const Helix &h2) {
        return is_adjacent(h1.back(), h2.front());
    }

    void merge_helix(const Helix &h1, Helix &h2) {
        for (auto it = h1.rbegin(); it != h1.rend(); it++) {
            h2.push_front(*it);
        }
    }

    Helices all_helices(Str seq, Int cutoff) {
        Helices hs, hs2, hs3;
        Int l = size(seq);
        Inds && inds = seq_to_inds(seq);
        for (Int i = 0; i < l; i++) {
            for (Int j = i + 1; j < l; j++) {
                Int a = inds[i], b = inds[j];
                if (is_bp(a, b)) {
                    Bp bp{i, j};
                    hs3.push_back({bp});
                }
            }
            for (auto && h2 : hs2) {
                Bool has_same = false;
                for (auto && h3 : hs3) {
                    if (is_same(h2, h3)) {
                        merge_helix(h2, h3);
                        has_same = true;
                        break;
                    }
                }
                if (!has_same) {
                    if (size(h2) >= cutoff) hs.push_back(h2);
                }
            }
            hs2 = std::move(hs3);
        }
        return hs;
    }

    REGISTER_NSP_COMPONENT(helix) {
        auto globals = par.getv("global");
        Str m = globals[1];
        if (m == "build") {
            Str seq = globals[2];
            JN_OUT << build_helix(seq);
        }
        else if (m == "all") {
            Str seq = globals[2];
            Int cutoff = 3;
            par.set(cutoff, "cutoff", "c");
            Helices && helices = all_helices(seq, cutoff);
            print_helices(helices);
        }
    }

}

END_JN


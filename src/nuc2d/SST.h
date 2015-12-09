#ifndef JIAN_NUC2D_SST_H
#define JIAN_NUC2D_SST_H

#include "Tree.h"

namespace jian {
namespace nuc2d {

struct Res {
    int num;
    char seq;
    char ss;
};

bool operator ==(const Res &res1, const Res &res2) {
    return res1.num == res2.num && res1.seq == res2.seq && res1.ss == res2.ss;
}

class Frag : public std::vector<Res> {
public:
    Frag() {}
    Frag(const std::string &seq, const std::string &ss) {
        if (seq.size() != ss.size()) throw "jian::nuc2d::Frag::Frag(const std::string &, const std::string &) error!";
        reserve(seq.size());
        int ss_index = 0;
        for (int i = 0; i < seq.size(); i++) {
            while (ss[ss_index] == '&') ss_index++;
            push_back({i, seq[i], ss[ss_index]});
            ss_index++;
        }
    }

    std::string seq() const {
        std::string seq;
        seq.reserve(size());
        for (auto &&res: *this) seq.push_back(res.seq);
        return seq;
    }

    std::string ss() const {
        std::string ss;
        ss.reserve(size());
        for (auto &&res: *this) ss.push_back(res.ss);
        return ss;
    }
};

std::ostream &operator <<(std::ostream &out, const Frag &frag) {
    out << "Fragment(" << frag.size() << "nt) ";
    for (auto &&res: frag) out << res.num << '-' << res.seq << '-' << res.ss << " ";
    return out;
}

class Hairpin {
public:
    std::array<Frag, 2> _helix;
    std::list<Frag> _loop;

    int size() const {
        int size = _helix[0].size() + _helix[1].size();
        for (auto &&frag: _loop) size += frag.size();
        return size;
    }
};

std::ostream &operator <<(std::ostream &out, const Hairpin &hairpin) {
    out << "Hairpin:\n";
    out << "helix:\n" << hairpin._helix[0] << "\n" << hairpin._helix[1] << "\n";
    out << "loop:\n";
    for (auto &&frag: hairpin._loop) out << frag << "\n";
    return out;
}

typedef Tree<Hairpin> SST;

} // namespace nuc2d
} // namespace jian

#endif


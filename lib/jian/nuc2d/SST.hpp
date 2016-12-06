#pragma once

#include "Tree.hpp"

BEGIN_JN
namespace nuc2d {

struct Res {
    int num; char seq; char ss; 

    bool operator ==(const Res &res) const {
        return num == res.num && seq == res.seq && ss == res.ss;
    }

};


class Frag : public std::vector<Res> {
public:
    Frag() {}

    Frag(const S &seq, const S &ss) {
        if (seq.size() != ss.size()) throw "jian::Frag::Frag(const S &, const S &) error!";
        reserve(seq.size());
        int ss_index = 0;
        for (int i = 0; i < seq.size(); i++) {
            while (ss[ss_index] == '&') ss_index++;
            push_back({i, seq[i], ss[ss_index]});
            ss_index++;
        }
    }

    S seq() const {
        S seq;
        seq.reserve(size());
        for (auto &&res: *this) seq.push_back(res.seq);
        return seq;
    }

    S ss() const {
        S ss;
        ss.reserve(size());
        for (auto &&res: *this) ss.push_back(res.ss);
        return ss;
    }

};

class Helix : public std::array<Frag, 2> {
public:
    S seq() const {
        return at(0).seq() + at(1).seq();
    }

    S ss() const {
        return at(0).ss() + at(1).ss();
    }

};

class Loop : public std::list<Frag> {
public:
    S seq() const {
        return std::accumulate(begin(), end(), std::string(), [](const S &s, const Frag &frag) {
            return s + frag.seq();
        });
    }

    S ss() const {
        return std::accumulate(begin(), end(), std::string(), [](const S &s, const Frag &frag) {
            return s + frag.ss();
        });
    }

};

class SSE {
public:
    Helix _helix;
    Loop _loop;

    Helix &helix() {
        return _helix;
    }

    const Helix &helix() const {
        return _helix;
    }

    Loop &loop() {
        return _loop;
    }

    const Loop &loop() const {
        return _loop;
    }

    int size() const {
        int size = _helix[0].size() + _helix[1].size();
        for (auto &&frag: _loop) size += frag.size();
        return size;
    }

};

typedef Tree<SSE> SST;

template<typename T, std::enable_if_t<std::is_same<std::decay_t<T>, Frag>::value, int> = 42>
std::ostream &operator <<(std::ostream &out, T &&frag) {
    out << "Fragment(" << frag.size() << "nt) ";
    for (auto &&res: frag) out << res.num << '-' << res.seq << '-' << res.ss << " ";
    return out;
}

template<typename T, std::enable_if_t<std::is_same<std::decay_t<T>, SSE>::value, int> = 42>
std::ostream &operator <<(std::ostream &out, T &&hairpin) {
    out << "SSE:\n";
    out << "helix:\n" << hairpin._helix[0] << "\n" << hairpin._helix[1] << "\n";
    out << "loop:\n";
    for (auto &&frag: hairpin._loop) out << frag << "\n";
    return out;
}

} // namespace nuc2d
END_JN


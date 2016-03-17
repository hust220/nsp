#pragma once

#include "../etl.hpp"

namespace jian {
namespace DihBoundSpace {

using key_t = std::vector<int>;
using val_t = double;

struct hash {
    unsigned int operator ()(const key_t &v) const {
        std::hash<int> h; return JN_FOLD(_1 ^ (h(_2) << 1), (unsigned int, 0), v); 
    }
};

struct equal_to {
    bool operator ()(const key_t &vec1, const key_t &vec2) const{
        return vec1.size() == vec2.size() && ! JN_EXISTS(_1 != _2, vec1, vec2);
    }
};

class DihBound : public std::unordered_map<key_t, val_t, hash, equal_to> {
public:
    val_t &operator ()(int a, int b, int c, int d) {
        return (*this)[{a, b, c, d}];
    }

    bool exists(int a, int b, int c, int d) const {
        return (*this).count({a, b, c, d});
    }

    friend std::ostream &operator <<(std::ostream &out, const DihBound &dih_bound) {
        EACH(i, dih_bound, EACH(n, i.first, SEE(n, ' ')); SEELN(':', i.second));
    }

};

} // DihBoundSpace

using DihBound = DihBoundSpace::DihBound;

} // namespace jian


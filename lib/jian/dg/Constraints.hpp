#pragma once

namespace jian {
namespace dg {

class hash {
public:
    int operator ()(const vector<int> &vec) const {
        return fold([](int result, int i){return result ^ (std::hash<int>{}(i) << 1);}, 0, vec);
    }
};

class equal_to {
public:
    bool operator ()(const vector<int> &vec1, const vector<int> &vec2) const {
        if (vec1.size() != vec2.size()) return false;
        for (int i = 0; i < vec1.size(); i++) if (vec1[i] != vec2[i]) return false;
        return true;
    }
};

class Term {
public:
    double avg;
    double weight;
};

class Constraints {
public:
    std::unordered_map<std::vector<int>, Term, hash, equal_to> _len_terms;
    std::unordered_map<std::vector<int>, Term, hash, equal_to> _ang_terms;
    std::unordered_map<std::vector<int>, Term, hash, equal_to> _dih_terms;

    Term &operator ()(int a, int b) {
        return _len_terms[{a, b}];
    }

//    const Term &operator ()(int a, int b) const {
//        return _len_terms[{a, b}];
//    }
//
    Term &operator ()(int a, int b, int c) {
        return _ang_terms[{a, b, c}];
    }

//    const Term &operator ()(int a, int b, int c) const {
//        return _ang_terms[{a, b, c}];
//    }
//
    Term &operator ()(int a, int b, int c, int d) {
        return _dih_terms[{a, b, c, d}];
    }

//    const Term &operator ()(int a, int b, int c, int d) const {
//        return _dih_terms[{a, b, c, d}];
//    }
//
    bool exists(int a, int b) const {
        return _len_terms.count({a, b});
    }

    bool exists(int a, int b, int c) const {
        return _ang_terms.count({a, b, c});
    }

    bool exists(int a, int b, int c, int d) const {
        return _dih_terms.count({a, b, c, d});
    }
};

} // namespace dg
} // namespace jian


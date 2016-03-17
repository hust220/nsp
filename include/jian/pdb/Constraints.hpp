#pragma once

namespace jian {

struct Constraint {
    std::vector<int> key;
    double value;
    double min, max;
};

inline auto make_contact(int i, int j) {
    Constraint c; append(c.key, i, j); return c;
}

inline auto make_distance(int i, int j, double value) {
    Constraint c; append(c.key, i, j); c.value = value; return c;
}

inline auto make_angle(int i, int j, int k, double value) {
    Constraint c; append(c.key, i, j, k); c.value = value; return c;
}

inline auto make_dihedral(int i, int j, int k, int l, double value) {
    Constraint c; append(c.key, i, j, k, l); c.value = value; return c;
}

struct Constraints {
    std::deque<Constraint> contacts;
    std::deque<Constraint> distances;
    std::deque<Constraint> angles;
    std::deque<Constraint> dihedrals;
};

template<typename T>
inline bool equal(T &&ls, int i,  int j) {
    return ls.size() == 2 && (ls[0] == i && ls[1] == j || ls[0] == j && ls[1] == i);
}

template<typename T>
inline bool equal(T &&ls, int i,  int j, int k) {
    return ls.size() == 3 && ls[1] == j && (ls[0] == i && ls[2] == k || ls[0] == k && ls[2] == i);
}

template<typename T>
inline bool equal(T &&ls, int i,  int j, int k,  int l) {
    return ls.size() == 4 && (ls[0] == i && ls[3] == l || ls[0] == l && ls[3] == i) && 
                             (ls[1] == j && ls[2] == k || ls[1] == k && ls[2] == j);
}

inline bool has_contact(const Constraints &cs, int i, int j) {
    for (auto &&c : cs.contacts) if (equal(c.key, i, j)) return true; return false;
}

inline bool has_distance(const Constraints &cs, int i, int j) {
    for (auto &&c : cs.distances) if (equal(c.key, i, j)) return true; return false;
}

inline bool has_angle(const Constraints &cs, int i, int j, int k) {
    for (auto &&c : cs.angles) if (equal(c.key, i, j, k)) return true; return false;
}

inline bool has_dihedral(const Constraints &cs, int i, int j, int k, int l) {
    for (auto &&c : cs.dihedrals) if (equal(c.key, i, j, k, l)) return true; return false;
}

inline void add_contact(Constraints &c, const Constraint &ct) {
    c.contacts.push_back(ct);
}

inline void add_distance(Constraints &c, const Constraint &ct) {
    c.distances.push_back(ct);
}

inline void add_angle(Constraints &c, const Constraint &ct) {
    c.angles.push_back(ct);
}

inline void add_dihedral(Constraints &c, const Constraint &ct) {
    c.dihedrals.push_back(ct);
}

inline Constraints read_constraints(const std::string &f) {
    Constraints c; if (f == "") return c;
    EACH_SPLIT_LINE(f.c_str(), " ", if (F.size() == 3) add_distance(c, make_distance(JN_INT(F[0]), JN_INT(F[1]), JN_DBL(F[2])));)
    return c;
}

} // namespace jian

